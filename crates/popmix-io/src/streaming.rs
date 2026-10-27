//! Streaming allele frequency computation.
//!
//! This module provides functions to compute allele frequencies directly from
//! genotype files without loading all data into memory. This matches the approach
//! used by admixtools2's cpp_packedancestrymap_to_afs.

use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

use memmap2::Mmap;
use ndarray::{Array2, Axis};
use rayon::prelude::*;

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
    crate::packedancestrymap::check_header(geno_path, n_samples, n_snps)?;
    let len = std::fs::metadata(geno_path)
        .map_err(|e| PopmixError::io(e, geno_path))?
        .len() as usize;
    let records = n_snps + 1;
    let bytes = len / records;
    if !len.is_multiple_of(records) || bytes < n_samples.div_ceil(4) {
        return Err(PopmixError::InvalidFileFormat {
            path: geno_path.to_string().into(),
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
/// This maps the file and computes each SNP's allele frequencies, the SNPs in
/// parallel, without loading the full genotype matrix ([`SnpMajor`]). This matches admixtools2's approach.
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
    let file = SnpMajor::open_packed(geno_path, samples.len(), snps.len())?;
    Ok(file.allele_frequencies(samples, snps, populations, auto_only))
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

    let mut file = File::open(geno_path).map_err(|e| PopmixError::io(e, geno_path))?;

    // Skip header
    file.seek(SeekFrom::Start(bytes_per_snp as u64))
        .map_err(|e| PopmixError::io(e, geno_path))?;

    let mut has_het = vec![false; n_samples];
    let mut snp_buffer = vec![0u8; bytes_per_snp];

    for _ in 0..n_check {
        file.read_exact(&mut snp_buffer)
            .map_err(|e| PopmixError::io(e, geno_path))?;

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

/// Like [`detect_ploidy_packedancestrymap`], for the samples in `selected`
/// only, stopping once each has shown a heterozygous call: a sample with
/// none among the first `n_snps_check` SNPs is pseudo-haploid. Other samples'
/// ploidy is left as it is.
pub fn detect_ploidy_packedancestrymap_selected(
    geno_path: &str,
    samples: &mut [Sample],
    selected: &[usize],
    n_snps_check: usize,
    n_snps_total: usize,
) -> Result<()> {
    let file = SnpMajor::open_packed(geno_path, samples.len(), n_snps_total)?;
    file.detect_ploidy(samples, selected, n_snps_check);
    Ok(())
}

/// The bytes of a [`SnpMajor`] file: mapped, or held in memory (a VCF,
/// decompressed and packed).
enum Bytes {
    Mapped(Mmap),
    Owned(Vec<u8>),
}

impl std::ops::Deref for Bytes {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        match self {
            Bytes::Mapped(map) => map,
            Bytes::Owned(bytes) => bytes,
        }
    }
}

/// How a [`SnpMajor`] file stores a genotype.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Encoding {
    /// PACKEDANCESTRYMAP: 2 bits, the first sample in a byte's high bits
    Packed,
    /// PLINK .bed: 2 bits, the first sample in a byte's low bits
    Bed,
    /// Text EIGENSTRAT: one character, a line per SNP
    Text,
}

/// A memory-mapped file of genotypes stored SNP by SNP: a PACKEDANCESTRYMAP
/// .geno, a PLINK .bed (SNP-major) or a text EIGENSTRAT .geno whose lines all
/// hold one character per sample.
pub struct SnpMajor {
    path: String,
    map: Bytes,
    /// Offset of the first SNP's record.
    start: usize,
    bytes_per_snp: usize,
    n_samples: usize,
    n_snps: usize,
    encoding: Encoding,
}

/// popmix genotype codes (0, 1, 2; 3 missing) of PLINK's 2-bit values:
/// 00 homozygous A1, 01 missing, 10 heterozygous, 11 homozygous A2.
const BED_CODES: [u8; 4] = [2, 3, 1, 0];

/// popmix genotype codes of text EIGENSTRAT characters: 0, 1, 2, and
/// missing for 9 and anything else, as the text reader reads them.
const TEXT_CODES: [u8; 256] = {
    let mut codes = [3u8; 256];
    codes[b'0' as usize] = 0;
    codes[b'1' as usize] = 1;
    codes[b'2' as usize] = 2;
    codes
};

impl SnpMajor {
    fn map(path: &str) -> Result<Mmap> {
        let file = File::open(path).map_err(|e| PopmixError::io(e, path))?;
        // SAFETY: the map is only read; like every reader of these files,
        // popmix assumes no other process truncates the file while it is open.
        unsafe { Mmap::map(&file) }.map_err(|e| PopmixError::io(e, path))
    }

    /// A PACKEDANCESTRYMAP .geno of `n_samples` samples and `n_snps` SNPs.
    pub fn open_packed(geno_path: &str, n_samples: usize, n_snps: usize) -> Result<Self> {
        let bytes_per_snp = record_length(geno_path, n_samples, n_snps)?;
        let map = Self::map(geno_path)?;
        if map.len() < 4 || &map[..4] != b"GENO" {
            return Err(PopmixError::InvalidFileFormat {
                path: geno_path.to_string().into(),
                details: "Expected PACKEDANCESTRYMAP format".to_string(),
            });
        }
        // record_length checked that the file holds the header and every SNP
        Ok(SnpMajor {
            path: geno_path.to_string(),
            map: Bytes::Mapped(map),
            start: bytes_per_snp,
            bytes_per_snp,
            n_samples,
            n_snps,
            encoding: Encoding::Packed,
        })
    }

    /// A SNP-major PLINK .bed of `n_samples` samples and `n_snps` SNPs.
    pub fn open_bed(bed_path: &str, n_samples: usize, n_snps: usize) -> Result<Self> {
        let map = Self::map(bed_path)?;
        if map.len() < 3 || map[..3] != [0x6c, 0x1b, 0x01] {
            return Err(PopmixError::InvalidFileFormat {
                path: bed_path.to_string().into(),
                details: "not a SNP-major PLINK .bed (magic 6c 1b 01)".to_string(),
            });
        }
        let bytes_per_snp = n_samples.div_ceil(4);
        let expected = 3 + n_snps * bytes_per_snp;
        if map.len() != expected {
            return Err(PopmixError::UnexpectedFileSize {
                path: bed_path.to_string().into(),
                expected: expected as u64,
                actual: map.len() as u64,
            });
        }
        Ok(SnpMajor {
            path: bed_path.to_string(),
            map: Bytes::Mapped(map),
            start: 3,
            bytes_per_snp,
            n_samples,
            n_snps,
            encoding: Encoding::Bed,
        })
    }

    /// A text EIGENSTRAT .geno of `n_samples` samples and `n_snps` SNPs, if
    /// every line holds exactly `n_samples` characters and a `\n` (as popmix,
    /// convertf and admixtools2 write them); None otherwise, for the text
    /// reader, which also takes padded or CRLF lines.
    pub fn open_text(geno_path: &str, n_samples: usize, n_snps: usize) -> Result<Option<Self>> {
        let map = Self::map(geno_path)?;
        let line = n_samples + 1;
        if map.len() != n_snps * line || (0..n_snps).any(|i| map[i * line + n_samples] != b'\n') {
            return Ok(None);
        }
        Ok(Some(SnpMajor {
            path: geno_path.to_string(),
            map: Bytes::Mapped(map),
            start: 0,
            bytes_per_snp: line,
            n_samples,
            n_snps,
            encoding: Encoding::Text,
        }))
    }

    /// Genotypes held in memory: `n_snps` records of `n_samples.div_ceil(4)`
    /// bytes packed as PACKEDANCESTRYMAP packs them, with no header. `path`
    /// names the file they came from, for errors.
    pub fn from_packed_records(
        path: String,
        records: Vec<u8>,
        n_samples: usize,
        n_snps: usize,
    ) -> Self {
        let bytes_per_snp = n_samples.div_ceil(4);
        assert_eq!(
            records.len(),
            n_snps * bytes_per_snp,
            "packed records of {path}"
        );
        SnpMajor {
            path,
            map: Bytes::Owned(records),
            start: 0,
            bytes_per_snp,
            n_samples,
            n_snps,
            encoding: Encoding::Packed,
        }
    }

    pub fn n_snps(&self) -> usize {
        self.n_snps
    }

    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// The bytes of SNP `snp`.
    #[inline]
    pub fn record(&self, snp: usize) -> &[u8] {
        let at = self.start + snp * self.bytes_per_snp;
        &self.map[at..at + self.bytes_per_snp]
    }

    /// The genotype code (0, 1, 2; 3 missing) of sample `s` in a record.
    #[inline]
    pub fn code(&self, record: &[u8], s: usize) -> u8 {
        match self.encoding {
            Encoding::Bed => BED_CODES[((record[s / 4] >> (2 * (s % 4))) & 0b11) as usize],
            // the first sample of a byte is in its high bits
            Encoding::Packed => (record[s / 4] >> (6 - 2 * (s % 4))) & 0b11,
            Encoding::Text => TEXT_CODES[record[s] as usize],
        }
    }

    /// Set the ploidy of the samples in `selected`: pseudo-haploid if none
    /// of the first `n_snps_check` SNPs is heterozygous, reading only until
    /// each has shown one. Other samples' ploidy is left as it is.
    pub fn detect_ploidy(&self, samples: &mut [Sample], selected: &[usize], n_snps_check: usize) {
        // chunks of SNPs in parallel, a round at a time, until every sample
        // has shown a heterozygous call or the SNPs run out
        const CHUNK: usize = 1 << 14;
        let n_check = n_snps_check.min(self.n_snps);
        let mut pending: Vec<usize> = selected.to_vec();
        let round = CHUNK * rayon::current_num_threads().max(1);
        let mut start = 0;
        while start < n_check && !pending.is_empty() {
            let end = (start + round).min(n_check);
            let found: Vec<bool> = (start..end)
                .step_by(CHUNK)
                .collect::<Vec<_>>()
                .into_par_iter()
                .map(|lo| {
                    let mut found = vec![false; pending.len()];
                    for snp in lo..(lo + CHUNK).min(end) {
                        let record = self.record(snp);
                        for (f, &s) in found.iter_mut().zip(&pending) {
                            *f |= self.code(record, s) == 1;
                        }
                    }
                    found
                })
                .reduce(
                    || vec![false; pending.len()],
                    |a, b| a.iter().zip(&b).map(|(x, y)| *x || *y).collect(),
                );
            let mut k = 0;
            pending.retain(|_| {
                k += 1;
                !found[k - 1]
            });
            start = end;
        }
        for &s in selected {
            samples[s].ploidy = if pending.contains(&s) {
                Ploidy::Pseudohaploid
            } else {
                Ploidy::Diploid
            };
        }
    }

    /// Allele frequencies of `populations` (samples with their ploidy set),
    /// each SNP on its own, the SNPs in parallel; with `auto_only`, of the
    /// autosomal SNPs only.
    pub fn allele_frequencies(
        &self,
        samples: &[Sample],
        snps: &[Snp],
        populations: &[&str],
        auto_only: bool,
    ) -> AlleleFrequencies {
        self.allele_frequencies_at(samples, snps, 0, populations, auto_only)
    }

    /// [`allele_frequencies`](Self::allele_frequencies) of a run of SNPs,
    /// `snps` being the file's SNPs `first..first + snps.len()`.
    pub fn allele_frequencies_at(
        &self,
        samples: &[Sample],
        snps: &[Snp],
        first: usize,
        populations: &[&str],
        auto_only: bool,
    ) -> AlleleFrequencies {
        let n_pops = populations.len();
        let (snp_indices, filtered_snps): (Vec<usize>, Vec<Snp>) = snps
            .iter()
            .enumerate()
            .filter(|(_, snp)| !auto_only || snp.is_autosomal())
            .map(|(idx, snp)| (first + idx, snp.clone()))
            .unzip();
        let pop_to_col: HashMap<&str, usize> = populations
            .iter()
            .enumerate()
            .map(|(c, &p)| (p, c))
            .collect();
        let mut pop_sample_info: Vec<Vec<(usize, u32)>> = vec![Vec::new(); n_pops];
        for (sample_idx, sample) in samples.iter().enumerate() {
            if let Some(&col) = pop_to_col.get(sample.population.as_str()) {
                let ploidy = sample.ploidy.effective_ploidy() as u32;
                pop_sample_info[col].push((sample_idx, ploidy));
            }
        }

        let mut frequencies = Array2::from_elem((filtered_snps.len(), n_pops), f64::NAN);
        let mut counts = Array2::zeros((filtered_snps.len(), n_pops));
        frequencies
            .axis_iter_mut(Axis(0))
            .into_par_iter()
            .zip(counts.axis_iter_mut(Axis(0)))
            .zip(&snp_indices)
            .for_each(|((mut freq, mut count), &snp)| {
                let record = self.record(snp);
                for (col, sample_info) in pop_sample_info.iter().enumerate() {
                    // admixtools2 adds genotype / (3 - ploidy) alleles over
                    // ploidy observed: a pseudo-haploid het counts half.
                    // Twice that is an integer, and every sum is exact, so
                    // counting in integers gives the same frequencies.
                    let (mut twice_alt, mut observed) = (0u64, 0u64);
                    for &(sample_idx, ploidy) in sample_info {
                        let g = self.code(record, sample_idx);
                        // without a branch, as missing calls are common
                        let called = (g != 3) as u64 * ploidy as u64;
                        twice_alt += g as u64 * called;
                        observed += called;
                    }
                    if observed > 0 {
                        freq[col] = twice_alt as f64 / 2.0 / observed as f64;
                        count[col] = observed as u32;
                    }
                }
            });
        AlleleFrequencies {
            frequencies,
            counts,
            populations: populations.iter().map(|&p| PopulationId::new(p)).collect(),
            snps: filtered_snps,
        }
    }
}

impl SnpMajor {
    /// [`SampleGenotypes::fill_samples`] for `n_snps` SNPs, SNP `j` being
    /// this file's SNP `snp_of(j).0`, its alleles swapped (0 <-> 2) if
    /// `snp_of(j).1`.
    ///
    /// [`SampleGenotypes::fill_samples`]: crate::writers::SampleGenotypes::fill_samples
    pub fn fill_samples_of(
        &self,
        samples: std::ops::Range<usize>,
        n_snps: usize,
        out: &mut [u8],
        snp_of: &(dyn Fn(usize) -> (usize, bool) + Sync),
    ) {
        // tiles of SNPs in parallel, each record read once for the block
        const TILE: usize = 4096;
        let width = samples.len();
        let tiles: Vec<(usize, Vec<u8>)> = (0..n_snps)
            .step_by(TILE)
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|lo| {
                let hi = (lo + TILE).min(n_snps);
                let mut tile = vec![0u8; width * (hi - lo)];
                for j in lo..hi {
                    let (snp, flipped) = snp_of(j);
                    let record = self.record(snp);
                    for (k, s) in samples.clone().enumerate() {
                        let g = self.code(record, s);
                        tile[k * (hi - lo) + j - lo] = if flipped && g != 1 && g != 3 {
                            2 - g
                        } else {
                            g
                        };
                    }
                }
                (lo, tile)
            })
            .collect();
        for (lo, tile) in tiles {
            let len = tile.len() / width.max(1);
            for k in 0..width {
                out[k * n_snps + lo..k * n_snps + lo + len]
                    .copy_from_slice(&tile[k * len..(k + 1) * len]);
            }
        }
    }
}

impl crate::writers::SampleGenotypes for SnpMajor {
    fn fill_samples(&self, samples: std::ops::Range<usize>, n_snps: usize, out: &mut [u8]) {
        self.fill_samples_of(samples, n_snps, out, &|j| (j, false));
    }
}

/// Genotype rows streamed from a [`SnpMajor`] file: the SNPs marked in
/// `keep`, in file order, for the given samples (by index in the file).
pub struct SnpMajorRows {
    file: SnpMajor,
    keep: Vec<bool>,
    samples: Vec<usize>,
    next_snp: usize,
}

/// Rows of a PACKEDANCESTRYMAP file (the name before PLINK was streamed too).
pub type PackedRows = SnpMajorRows;

impl SnpMajorRows {
    /// Rows of `file`: `keep` has one entry per SNP.
    pub fn from_file(file: SnpMajor, keep: Vec<bool>, samples: Vec<usize>) -> Result<Self> {
        let n_samples = file.n_samples;
        if let Some(&bad) = samples.iter().find(|&&s| s >= n_samples) {
            return Err(PopmixError::InvalidFileFormat {
                path: file.path.clone().into(),
                details: format!("sample {bad} out of range for {n_samples} samples"),
            });
        }
        if keep.len() != file.n_snps {
            return Err(PopmixError::InvalidFileFormat {
                path: file.path.clone().into(),
                details: format!("{} SNPs selected of {}", keep.len(), file.n_snps),
            });
        }
        Ok(SnpMajorRows {
            file,
            keep,
            samples,
            next_snp: 0,
        })
    }

    /// Rows of a PACKEDANCESTRYMAP .geno with `keep.len()` SNPs.
    pub fn open(
        geno_path: &str,
        n_samples: usize,
        keep: Vec<bool>,
        samples: Vec<usize>,
    ) -> Result<Self> {
        let file = SnpMajor::open_packed(geno_path, n_samples, keep.len())?;
        Self::from_file(file, keep, samples)
    }

    /// Rows of a SNP-major PLINK .bed with `keep.len()` SNPs.
    pub fn open_bed(
        bed_path: &str,
        n_samples: usize,
        keep: Vec<bool>,
        samples: Vec<usize>,
    ) -> Result<Self> {
        let file = SnpMajor::open_bed(bed_path, n_samples, keep.len())?;
        Self::from_file(file, keep, samples)
    }
}

impl popmix_core::GenotypeRows for SnpMajorRows {
    fn next_rows(&mut self, out: &mut [u8], width: usize) -> Result<()> {
        // the SNPs of the rows, then the rows decoded in parallel
        let n = out.len() / width.max(1);
        let mut snps = Vec::with_capacity(n);
        while snps.len() < n {
            while self.next_snp < self.keep.len() && !self.keep[self.next_snp] {
                self.next_snp += 1;
            }
            if self.next_snp >= self.keep.len() {
                return Err(PopmixError::InvalidFileFormat {
                    path: self.file.path.clone().into(),
                    details: "read past the last SNP".to_string(),
                });
            }
            snps.push(self.next_snp);
            self.next_snp += 1;
        }
        out.par_chunks_mut(width.max(1))
            .zip(snps)
            .for_each(|(row, snp)| {
                let record = self.file.record(snp);
                for (o, &s) in row.iter_mut().zip(&self.samples) {
                    *o = self.file.code(record, s);
                }
            });
        Ok(())
    }

    fn next_row(&mut self, out: &mut [u8]) -> Result<()> {
        while self.next_snp < self.keep.len() && !self.keep[self.next_snp] {
            self.next_snp += 1;
        }
        if self.next_snp >= self.keep.len() {
            return Err(PopmixError::InvalidFileFormat {
                path: self.file.path.clone().into(),
                details: "read past the last SNP".to_string(),
            });
        }
        let record = self.file.record(self.next_snp);
        self.next_snp += 1;
        for (o, &s) in out.iter_mut().zip(&self.samples) {
            *o = self.file.code(record, s);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
