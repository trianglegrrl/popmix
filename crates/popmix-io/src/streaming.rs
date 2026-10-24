//! Streaming allele frequency computation.
//!
//! This module provides functions to compute allele frequencies directly from
//! genotype files without loading all data into memory. This matches the approach
//! used by admixtools2's cpp_packedancestrymap_to_afs.

use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

use ndarray::Array2;

use popmix_core::{AlleleFrequencies, Ploidy, PopmixError, PopulationId, Result, Sample, Snp};

/// Decode 4 genotypes from a single byte using PACKEDANCESTRYMAP encoding.
///
/// IMPORTANT: In PACKEDANCESTRYMAP format, the first individual (index 0) is stored
/// in the HIGH bits (6-7) of each byte, and the fourth individual (index 3) is in
/// the LOW bits (0-1). This matches the admixtools2 C++ implementation.
#[inline]
fn decode_packed_byte(byte: u8) -> [u8; 4] {
    [
        (byte >> 6) & 0b11, // Individual 0: bits 6-7 (high bits)
        (byte >> 4) & 0b11, // Individual 1: bits 4-5
        (byte >> 2) & 0b11, // Individual 2: bits 2-3
        byte & 0b11,        // Individual 3: bits 0-1 (low bits)
    ]
}

/// Bytes per record (the header and each SNP) of a PACKEDANCESTRYMAP file.
///
/// Records are padded to at least 48 bytes, so this is not ceil(n / 4) for
/// fewer than 189 samples. As admixtools2 (`cpp_packedancestrymap_to_afs`),
/// it is the file size over the number of records; a file whose size is not
/// a whole number of records long enough for every sample is an error.
fn record_length(geno_path: &str, n_samples: usize, n_snps: usize) -> Result<usize> {
    let len = std::fs::metadata(geno_path)
        .map_err(|e| PopmixError::io_with_path(e, geno_path))?
        .len() as usize;
    let records = n_snps + 1;
    let bytes = len / records;
    if !len.is_multiple_of(records) || bytes < n_samples.div_ceil(4) {
        return Err(PopmixError::InvalidFileFormat {
            file: geno_path.to_string(),
            details: format!(
                "{len} bytes is not {records} records (header and {n_snps} SNPs) of at least {} bytes",
                n_samples.div_ceil(4)
            ),
        });
    }
    Ok(bytes)
}

/// Compute allele frequencies by streaming through a PACKEDANCESTRYMAP file.
///
/// This reads the file SNP by SNP and computes allele frequencies without
/// ever loading the full genotype matrix. This matches admixtools2's approach.
///
/// # Arguments
/// * `geno_path` - Path to the .geno file
/// * `samples` - Sample metadata (with ploidy already set)
/// * `snps` - SNP metadata
/// * `populations` - List of population names to compute frequencies for
/// * `auto_only` - If true, only include autosomal SNPs (chr 1-22)
///
/// # Returns
/// AlleleFrequencies for the specified populations
pub fn compute_allele_frequencies_packedancestrymap(
    geno_path: &str,
    samples: &[Sample],
    snps: &[Snp],
    populations: &[&str],
    auto_only: bool,
) -> Result<AlleleFrequencies> {
    let n_samples = samples.len();
    let n_snps_total = snps.len();
    let n_pops = populations.len();

    // Filter SNPs if auto_only is set
    let (snp_indices, filtered_snps): (Vec<usize>, Vec<Snp>) = if auto_only {
        snps.iter()
            .enumerate()
            .filter(|(_, snp)| snp.is_autosomal())
            .map(|(idx, snp)| (idx, snp.clone()))
            .unzip()
    } else {
        snps.iter()
            .enumerate()
            .map(|(idx, snp)| (idx, snp.clone()))
            .unzip()
    };

    let n_snps_filtered = filtered_snps.len();

    // Build population -> (column index, sample indices with ploidy) mapping
    let mut pop_to_col: HashMap<&str, usize> = HashMap::new();
    let mut pop_sample_info: Vec<Vec<(usize, u32)>> = vec![Vec::new(); n_pops];

    for (col_idx, &pop_name) in populations.iter().enumerate() {
        pop_to_col.insert(pop_name, col_idx);
    }

    for (sample_idx, sample) in samples.iter().enumerate() {
        let pop_name = sample.population.as_str();
        if let Some(&col_idx) = pop_to_col.get(pop_name) {
            let ploidy = sample.ploidy.effective_ploidy() as u32;
            pop_sample_info[col_idx].push((sample_idx, ploidy));
        }
    }

    let bytes_per_snp = record_length(geno_path, n_samples, n_snps_total)?;

    // Open file
    let mut file = File::open(geno_path).map_err(|e| PopmixError::io_with_path(e, geno_path))?;

    // Verify magic number
    let mut magic = [0u8; 4];
    file.read_exact(&mut magic)
        .map_err(|e| PopmixError::io_with_path(e, geno_path))?;
    if &magic != b"GENO" {
        return Err(PopmixError::InvalidFileFormat {
            file: geno_path.to_string(),
            details: "Expected PACKEDANCESTRYMAP format".to_string(),
        });
    }

    // Seek past header (first row)
    file.seek(SeekFrom::Start(bytes_per_snp as u64))
        .map_err(|e| PopmixError::io_with_path(e, geno_path))?;

    // Initialize output arrays (only for filtered SNPs)
    let mut frequencies = Array2::from_elem((n_snps_filtered, n_pops), f64::NAN);
    let mut counts = Array2::zeros((n_snps_filtered, n_pops));

    // Buffer for one SNP's worth of data
    let mut snp_buffer = vec![0u8; bytes_per_snp];

    // Track which filtered SNP we're on
    let mut filtered_idx = 0;
    let mut next_snp_to_process = if snp_indices.is_empty() {
        usize::MAX
    } else {
        snp_indices[0]
    };

    // Process SNPs one at a time
    for file_snp_idx in 0..n_snps_total {
        // Read this SNP's data
        file.read_exact(&mut snp_buffer)
            .map_err(|e| PopmixError::io_with_path(e, geno_path))?;

        // Skip if this SNP is not in our filtered set
        if file_snp_idx != next_snp_to_process {
            continue;
        }

        // Compute allele frequencies for each population
        for (col_idx, sample_info) in pop_sample_info.iter().enumerate() {
            if sample_info.is_empty() {
                continue;
            }

            let mut alt_alleles: f64 = 0.0;
            let mut observed_alleles: f64 = 0.0;

            for &(sample_idx, ploidy) in sample_info {
                let byte_idx = sample_idx / 4;
                let pos_in_byte = sample_idx % 4;
                let geno_bits = decode_packed_byte(snp_buffer[byte_idx])[pos_in_byte];

                // geno_bits: 00=HomRef(0), 01=Het(1), 10=HomAlt(2), 11=Missing
                if geno_bits != 3 {
                    // Not missing
                    let geno_val = geno_bits as f64;
                    // admixtools2 formula: alt_alleles += val / (3.0 - ploidy)
                    // For diploid (ploidy=2): val / 1.0 = val
                    // For pseudohaploid (ploidy=1): val / 2.0 (so het=1 becomes 0.5)
                    alt_alleles += geno_val / (3.0 - ploidy as f64);
                    observed_alleles += ploidy as f64;
                }
            }

            if observed_alleles > 0.0 {
                frequencies[[filtered_idx, col_idx]] = alt_alleles / observed_alleles;
                counts[[filtered_idx, col_idx]] = observed_alleles as u32;
            }
        }

        // Move to next filtered SNP
        filtered_idx += 1;
        next_snp_to_process = if filtered_idx < snp_indices.len() {
            snp_indices[filtered_idx]
        } else {
            usize::MAX
        };
    }

    // Create population IDs
    let population_ids: Vec<PopulationId> =
        populations.iter().map(|&p| PopulationId::new(p)).collect();

    Ok(AlleleFrequencies {
        frequencies,
        counts,
        populations: population_ids,
        snps: filtered_snps,
    })
}

/// Detect ploidy for samples by scanning for heterozygous calls.
///
/// Reads the first n_snps_check SNPs and marks samples with Het calls as diploid.
pub fn detect_ploidy_packedancestrymap(
    geno_path: &str,
    samples: &mut [Sample],
    n_snps_check: usize,
    n_snps_total: usize,
) -> Result<()> {
    let n_samples = samples.len();
    let n_check = n_snps_check.min(n_snps_total);
    let bytes_per_snp = record_length(geno_path, n_samples, n_snps_total)?;

    let mut file = File::open(geno_path).map_err(|e| PopmixError::io_with_path(e, geno_path))?;

    // Skip header
    file.seek(SeekFrom::Start(bytes_per_snp as u64))
        .map_err(|e| PopmixError::io_with_path(e, geno_path))?;

    let mut has_het = vec![false; n_samples];
    let mut snp_buffer = vec![0u8; bytes_per_snp];

    for _ in 0..n_check {
        file.read_exact(&mut snp_buffer)
            .map_err(|e| PopmixError::io_with_path(e, geno_path))?;

        for (sample_idx, found_het) in has_het.iter_mut().enumerate() {
            if *found_het {
                continue; // Already found het for this sample
            }
            let byte_idx = sample_idx / 4;
            let pos_in_byte = sample_idx % 4;
            let geno_bits = decode_packed_byte(snp_buffer[byte_idx])[pos_in_byte];
            if geno_bits == 1 {
                // Het
                *found_het = true;
            }
        }
    }

    // Update ploidy
    for (sample_idx, sample) in samples.iter_mut().enumerate() {
        sample.ploidy = if has_het[sample_idx] {
            Ploidy::Diploid
        } else {
            Ploidy::Pseudohaploid
        };
    }

    Ok(())
}

/// Genotype rows streamed from a PACKEDANCESTRYMAP file: the SNPs marked in
/// `keep`, in file order, for the given samples (by index in the .ind file).
pub struct PackedRows {
    file: std::io::BufReader<File>,
    path: String,
    keep: Vec<bool>,
    samples: Vec<usize>,
    buffer: Vec<u8>,
    next_snp: usize,
}

impl PackedRows {
    pub fn open(
        geno_path: &str,
        n_samples: usize,
        keep: Vec<bool>,
        samples: Vec<usize>,
    ) -> Result<Self> {
        let bytes_per_snp = record_length(geno_path, n_samples, keep.len())?;
        if let Some(&bad) = samples.iter().find(|&&s| s >= n_samples) {
            return Err(PopmixError::InvalidFileFormat {
                file: geno_path.to_string(),
                details: format!("sample {bad} out of range for {n_samples} samples"),
            });
        }
        let mut file =
            File::open(geno_path).map_err(|e| PopmixError::io_with_path(e, geno_path))?;
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)
            .map_err(|e| PopmixError::io_with_path(e, geno_path))?;
        if &magic != b"GENO" {
            return Err(PopmixError::InvalidFileFormat {
                file: geno_path.to_string(),
                details: "Expected PACKEDANCESTRYMAP format".to_string(),
            });
        }
        file.seek(SeekFrom::Start(bytes_per_snp as u64))
            .map_err(|e| PopmixError::io_with_path(e, geno_path))?;
        Ok(PackedRows {
            file: std::io::BufReader::with_capacity(1 << 22, file),
            path: geno_path.to_string(),
            keep,
            samples,
            buffer: vec![0u8; bytes_per_snp],
            next_snp: 0,
        })
    }
}

impl popmix_core::GenotypeRows for PackedRows {
    fn next_row(&mut self, out: &mut [u8]) -> Result<()> {
        loop {
            if self.next_snp >= self.keep.len() {
                return Err(PopmixError::InvalidFileFormat {
                    file: self.path.clone(),
                    details: "read past the last SNP".to_string(),
                });
            }
            self.file
                .read_exact(&mut self.buffer)
                .map_err(|e| PopmixError::io_with_path(e, &self.path))?;
            self.next_snp += 1;
            if self.keep[self.next_snp - 1] {
                break;
            }
        }
        for (o, &s) in out.iter_mut().zip(&self.samples) {
            *o = decode_packed_byte(self.buffer[s / 4])[s % 4];
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eigenstrat::{read_ind_file, read_snp_file};
    use crate::PackedAncestryMapReader;

    const SMALL: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/golden/dates/family_chr20_22"
    );

    /// 54 samples: each record is padded to 48 bytes, not ceil(54 / 4) = 14.
    #[test]
    fn small_files_stream_with_48_byte_records() {
        let geno_path = format!("{SMALL}.geno");
        let mut samples = read_ind_file(&format!("{SMALL}.ind")).unwrap();
        let snps = read_snp_file(&format!("{SMALL}.snp")).unwrap();
        assert_eq!(samples.len(), 54);
        detect_ploidy_packedancestrymap(&geno_path, &mut samples, 1000, snps.len()).unwrap();
        let pops = ["French", "Yoruba"];
        let streamed =
            compute_allele_frequencies_packedancestrymap(&geno_path, &samples, &snps, &pops, false)
                .unwrap();

        let geno = PackedAncestryMapReader::new(SMALL)
            .unwrap()
            .read_all()
            .unwrap();
        let in_memory = AlleleFrequencies::from_genotypes_for_populations(&geno, Some(&pops));
        assert_eq!(streamed.n_snps(), in_memory.n_snps());
        for (col, pop) in pops.iter().enumerate() {
            let other = in_memory.pop_index(pop).unwrap();
            for snp in 0..streamed.n_snps() {
                let (a, b) = (
                    streamed.frequencies[[snp, col]],
                    in_memory.frequencies[[snp, other]],
                );
                assert!(
                    a == b || (a.is_nan() && b.is_nan()),
                    "{pop} SNP {snp}: {a} != {b}"
                );
            }
        }
    }

    #[test]
    fn a_truncated_file_is_an_error() {
        let dir = std::env::temp_dir().join(format!("popmix-streaming-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let geno_path = dir.join("short.geno");
        let full = std::fs::read(format!("{SMALL}.geno")).unwrap();
        std::fs::write(&geno_path, &full[..full.len() - 10]).unwrap();
        let samples = read_ind_file(&format!("{SMALL}.ind")).unwrap();
        let snps = read_snp_file(&format!("{SMALL}.snp")).unwrap();
        let err = compute_allele_frequencies_packedancestrymap(
            geno_path.to_str().unwrap(),
            &samples,
            &snps,
            &["French"],
            false,
        );
        assert!(err.is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn packed_rows_match_the_in_memory_genotypes() {
        use popmix_core::GenotypeRows;
        let geno = PackedAncestryMapReader::new(SMALL)
            .unwrap()
            .read_all()
            .unwrap();
        let keep: Vec<bool> = (0..geno.n_snps()).map(|i| i % 3 != 1).collect();
        let samples = vec![53, 0, 7, 20];
        let mut rows =
            PackedRows::open(&format!("{SMALL}.geno"), 54, keep.clone(), samples.clone()).unwrap();
        let mut row = vec![0u8; samples.len()];
        for snp in (0..geno.n_snps()).filter(|&i| keep[i]) {
            rows.next_row(&mut row).unwrap();
            let want: Vec<u8> = samples.iter().map(|&s| geno.snp_raw(snp)[s]).collect();
            assert_eq!(row, want, "SNP {snp}");
        }
        assert!(rows.next_row(&mut row).is_err());
    }
}
