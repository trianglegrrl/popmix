//! PACKEDANCESTRYMAP binary format reader.
//!
//! PACKEDANCESTRYMAP uses the same .snp and .ind files as EIGENSTRAT,
//! but the .geno file is binary with 2-bit packed genotypes.
//!
//! Format:
//! - Record size: max(48, ceil(n_samples / 4)) bytes per record
//!   (minimum 48 bytes per the AdmixTools specification)
//! - Header: First record contains "GENO {n_samples} {n_snps} {hash}..." padded with zeros
//! - Data: SNP-major, 2-bit packed genotypes (4 genotypes per byte)
//!
//! Genotype encoding (2 bits per sample):
//! - 00 = Homozygous reference (0)
//! - 01 = Heterozygous (1)
//! - 10 = Homozygous alternate (2)
//! - 11 = Missing (9)

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;
use ndarray::Array2;

use popmix_core::{Genotype, GenotypeMatrix, Ploidy, PopmixError, Result, Sample, Snp};

/// Check if a .geno file is in packed binary format vs text EIGENSTRAT.
///
/// Packed format files start with "GENO" magic string.
pub fn is_packed_format(geno_path: &str) -> Result<bool> {
    use std::io::Read;

    let mut file = File::open(geno_path)
        .map_err(|e| PopmixError::io_with_path(e, geno_path))?;

    let mut magic = [0u8; 4];
    file.read_exact(&mut magic)
        .map_err(|e| PopmixError::io_with_path(e, geno_path))?;

    Ok(&magic == b"GENO")
}

/// Reader for PACKEDANCESTRYMAP binary format files.
pub struct PackedAncestryMapReader {
    prefix: String,
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    geno_mmap: Option<Mmap>,
    bytes_per_snp: usize,
}

impl PackedAncestryMapReader {
    /// Open a PACKEDANCESTRYMAP dataset by prefix.
    ///
    /// Uses the same .ind and .snp files as EIGENSTRAT.
    pub fn new(prefix: impl AsRef<Path>) -> Result<Self> {
        let prefix = prefix.as_ref().to_string_lossy().to_string();

        // Use EIGENSTRAT functions for .ind and .snp files (same format)
        let ind_path = format!("{}.ind", prefix);
        let samples = crate::eigenstrat::read_ind_file(&ind_path)?;

        let snp_path = format!("{}.snp", prefix);
        let snps = crate::eigenstrat::read_snp_file(&snp_path)?;

        // Verify the .geno file is in packed format
        let geno_path = format!("{}.geno", prefix);
        if !is_packed_format(&geno_path)? {
            return Err(PopmixError::InvalidFileFormat {
                file: geno_path,
                details: "Expected PACKEDANCESTRYMAP format, found EIGENSTRAT text format".to_string(),
            });
        }

        // Calculate bytes per SNP: max(48, ceil(n_samples / 4))
        // Per AdmixTools spec, minimum record size is 48 bytes
        let bytes_per_snp = samples.len().div_ceil(4).max(48);

        Ok(PackedAncestryMapReader {
            prefix,
            samples,
            snps,
            geno_mmap: None,
            bytes_per_snp,
        })
    }

    /// Get the number of samples.
    pub fn num_samples(&self) -> usize {
        self.samples.len()
    }

    /// Get the number of SNPs.
    pub fn num_snps(&self) -> usize {
        self.snps.len()
    }

    /// Get sample metadata.
    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }

    /// Get SNP metadata.
    pub fn snps(&self) -> &[Snp] {
        &self.snps
    }

    /// Get the file prefix.
    pub fn prefix(&self) -> &str {
        &self.prefix
    }

    /// Memory-map the .geno file for efficient reading.
    fn ensure_mmap(&mut self) -> Result<&Mmap> {
        if self.geno_mmap.is_none() {
            let geno_path = format!("{}.geno", self.prefix);
            let file = File::open(&geno_path)
                .map_err(|e| PopmixError::io_with_path(e, &geno_path))?;

            let mmap = unsafe {
                Mmap::map(&file).map_err(|e| PopmixError::io_with_path(e, &geno_path))?
            };

            // Verify magic number
            if mmap.len() < 4 || &mmap[0..4] != b"GENO" {
                return Err(PopmixError::InvalidFileFormat {
                    file: geno_path,
                    details: "Invalid PACKEDANCESTRYMAP magic number".to_string(),
                });
            }

            // Verify file size: header (1 row) + data (n_snps rows)
            let expected_size = self.bytes_per_snp + self.snps.len() * self.bytes_per_snp;
            if mmap.len() != expected_size {
                return Err(PopmixError::UnexpectedFileSize {
                    file: geno_path,
                    expected: expected_size as u64,
                    actual: mmap.len() as u64,
                });
            }

            self.geno_mmap = Some(mmap);
        }

        Ok(self.geno_mmap.as_ref().unwrap())
    }

    /// Read all genotypes into a GenotypeMatrix.
    pub fn read_all(&mut self) -> Result<GenotypeMatrix> {
        self.read_snp_range(0..self.snps.len())
    }

    /// Read genotypes for a range of SNPs.
    pub fn read_snp_range(&mut self, range: std::ops::Range<usize>) -> Result<GenotypeMatrix> {
        // Copy values before mmap borrow
        let n_samples = self.samples.len();
        let n_snps = range.len();
        let bytes_per_snp = self.bytes_per_snp;
        let snps: Vec<Snp> = range.clone().map(|i| self.snps[i].clone()).collect();
        let samples = self.samples.clone();

        let mmap = self.ensure_mmap()?;

        // Allocate output matrix
        let mut data = Array2::zeros((n_snps, n_samples));

        // Header is one row (bytes_per_snp bytes), then data starts
        let header_size = bytes_per_snp;

        // Read each SNP
        for (out_idx, snp_idx) in range.enumerate() {
            let offset = header_size + snp_idx * bytes_per_snp;
            let snp_bytes = &mmap[offset..offset + bytes_per_snp];

            // Decode genotypes
            for (sample_idx, geno) in decode_packed_genotypes(snp_bytes, n_samples).enumerate() {
                data[[out_idx, sample_idx]] = geno as u8;
            }
        }

        GenotypeMatrix::new(data, samples, snps)
    }

    /// Read genotypes for a single SNP.
    pub fn read_snp(&mut self, snp_idx: usize) -> Result<Vec<Genotype>> {
        // Copy values before mmap borrow
        let n_samples = self.samples.len();
        let bytes_per_snp = self.bytes_per_snp;

        let mmap = self.ensure_mmap()?;

        let header_size = bytes_per_snp;
        let offset = header_size + snp_idx * bytes_per_snp;
        let snp_bytes = &mmap[offset..offset + bytes_per_snp];

        Ok(decode_packed_genotypes(snp_bytes, n_samples).collect())
    }

    /// Detect ploidy for each sample based on heterozygosity.
    pub fn detect_ploidy(&mut self, n_snps_check: usize) -> Result<()> {
        // Copy values before mmap borrow
        let n_samples = self.samples.len();
        let n_check = n_snps_check.min(self.snps.len());
        let bytes_per_snp = self.bytes_per_snp;

        // Compute has_het in a scope where mmap is borrowed
        let has_het = {
            let mmap = self.ensure_mmap()?;
            let header_size = bytes_per_snp;
            let mut has_het = vec![false; n_samples];

            for snp_idx in 0..n_check {
                let offset = header_size + snp_idx * bytes_per_snp;
                let snp_bytes = &mmap[offset..offset + bytes_per_snp];

                for (sample_idx, geno) in decode_packed_genotypes(snp_bytes, n_samples).enumerate() {
                    if geno == Genotype::Het {
                        has_het[sample_idx] = true;
                    }
                }
            }
            has_het
        };

        // Update ploidy (mmap borrow is now dropped)
        for (sample_idx, sample) in self.samples.iter_mut().enumerate() {
            sample.ploidy = if has_het[sample_idx] {
                Ploidy::Diploid
            } else {
                Ploidy::Pseudohaploid
            };
        }

        Ok(())
    }
}

/// Decode 4 genotypes from a single byte using PACKEDANCESTRYMAP encoding.
///
/// PACKEDANCESTRYMAP stores genotypes as 2 bits each:
/// - 00 = Homozygous reference (0)
/// - 01 = Heterozygous (1)
/// - 10 = Homozygous alternate (2)
/// - 11 = Missing (9)
///
/// IMPORTANT: In PACKEDANCESTRYMAP format, the first individual (index 0) is stored
/// in the HIGH bits (6-7) of each byte, and the fourth individual (index 3) is in
/// the LOW bits (0-1). This matches the DATES/admixtools C implementation.
#[inline]
fn decode_packed_byte(byte: u8) -> [Genotype; 4] {
    [
        decode_2bits((byte >> 6) & 0b11), // Sample 0: bits 6-7 (HIGH)
        decode_2bits((byte >> 4) & 0b11), // Sample 1: bits 4-5
        decode_2bits((byte >> 2) & 0b11), // Sample 2: bits 2-3
        decode_2bits(byte & 0b11),        // Sample 3: bits 0-1 (LOW)
    ]
}

#[inline]
fn decode_2bits(bits: u8) -> Genotype {
    match bits {
        0b00 => Genotype::HomRef,
        0b01 => Genotype::Het,
        0b10 => Genotype::HomAlt,
        _ => Genotype::Missing, // 0b11
    }
}

/// Iterator over genotypes for a single SNP from raw bytes.
fn decode_packed_genotypes(bytes: &[u8], n_samples: usize) -> impl Iterator<Item = Genotype> + '_ {
    bytes
        .iter()
        .flat_map(|&byte| decode_packed_byte(byte))
        .take(n_samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_packed_byte() {
        // All HomRef (0b00000000)
        let genos = decode_packed_byte(0b00000000);
        assert!(genos.iter().all(|g| *g == Genotype::HomRef));

        // All Missing (0b11111111)
        let genos = decode_packed_byte(0b11111111);
        assert!(genos.iter().all(|g| *g == Genotype::Missing));

        // Mixed: Missing, HomAlt, Het, HomRef (reading left to right in bits)
        // Byte: 11_10_01_00 = 0xE4
        // Sample 0 (bits 6-7): 11 = Missing
        // Sample 1 (bits 4-5): 10 = HomAlt
        // Sample 2 (bits 2-3): 01 = Het
        // Sample 3 (bits 0-1): 00 = HomRef
        let genos = decode_packed_byte(0b11_10_01_00);
        assert_eq!(genos[0], Genotype::Missing); // bits 6-7 = 11
        assert_eq!(genos[1], Genotype::HomAlt);  // bits 4-5 = 10
        assert_eq!(genos[2], Genotype::Het);     // bits 2-3 = 01
        assert_eq!(genos[3], Genotype::HomRef);  // bits 0-1 = 00
    }
}
