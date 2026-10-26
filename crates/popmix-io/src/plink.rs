//! PLINK binary format (.bed, .bim, .fam) reader.
//!
//! PLINK stores genotypes in a binary format with 2 bits per genotype,
//! allowing 4 genotypes per byte. The format consists of three files:
//! - `.bed`: Binary genotype matrix with magic number header
//! - `.bim`: SNP metadata (one line per SNP)
//! - `.fam`: Sample metadata (one line per sample)
//!
//! The .bed file uses SNP-major mode (default) where each SNP's genotypes
//! are stored contiguously.

use std::fs::File;
use std::path::Path;

use std::sync::OnceLock;

use memmap2::Mmap;
use ndarray::Array2;

use popmix_core::{
    Chromosome, Genotype, GenotypeMatrix, Ploidy, PopmixError, PopulationId, Result, Sample, Sex,
    Snp,
};

/// PLINK .bed file magic number (first 3 bytes).
const PLINK_MAGIC: [u8; 3] = [0x6c, 0x1b, 0x01];

/// Reader for PLINK binary format files.
pub struct PlinkReader {
    prefix: String,
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    /// Mapped on first read; a `OnceLock` so reads take `&self`.
    bed_mmap: OnceLock<Mmap>,
    bytes_per_snp: usize,
}

impl PlinkReader {
    /// Open a PLINK dataset by prefix.
    ///
    /// The prefix should be the common part of the filename before the extension.
    /// For example, if files are `data.bed`, `data.bim`, `data.fam`, the prefix is `data`.
    pub fn new(prefix: impl AsRef<Path>) -> Result<Self> {
        let prefix = prefix.as_ref().to_string_lossy().to_string();

        // Read .fam file (samples)
        let fam_path = format!("{}.fam", prefix);
        let samples = read_fam_file(&fam_path)?;

        // Read .bim file (SNPs)
        let bim_path = format!("{}.bim", prefix);
        let snps = read_bim_file(&bim_path)?;

        // Calculate bytes per SNP (ceil(n_samples / 4))
        let bytes_per_snp = samples.len().div_ceil(4);

        Ok(PlinkReader {
            prefix,
            samples,
            snps,
            bed_mmap: OnceLock::new(),
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

    /// Memory-map the .bed file for efficient reading.
    fn ensure_mmap(&self) -> Result<&Mmap> {
        if let Some(map) = self.bed_mmap.get() {
            return Ok(map);
        }
        {
            let bed_path = format!("{}.bed", self.prefix);
            let file =
                File::open(&bed_path).map_err(|e| PopmixError::io_with_path(e, &bed_path))?;

            let mmap =
                // SAFETY: the map is only read; like every reader of these files,
            // popmix assumes no other process truncates the file while it is open.
            unsafe { Mmap::map(&file).map_err(|e| PopmixError::io_with_path(e, &bed_path))? };

            // Verify magic number and mode
            if mmap.len() < 3 {
                return Err(PopmixError::InvalidFileFormat {
                    file: bed_path,
                    details: "File too small to contain PLINK header".to_string(),
                });
            }

            if mmap[0..3] != PLINK_MAGIC {
                return Err(PopmixError::InvalidFileFormat {
                    file: bed_path,
                    details: format!(
                        "Invalid PLINK magic number: {:02x} {:02x} {:02x}",
                        mmap[0], mmap[1], mmap[2]
                    ),
                });
            }

            // Verify file size
            let expected_size = 3 + self.snps.len() * self.bytes_per_snp;
            if mmap.len() != expected_size {
                return Err(PopmixError::UnexpectedFileSize {
                    file: bed_path,
                    expected: expected_size as u64,
                    actual: mmap.len() as u64,
                });
            }

            // another thread may have mapped it meanwhile: either map serves
            let _ = self.bed_mmap.set(mmap);
        }

        Ok(self.bed_mmap.get().expect("mapped above"))
    }

    /// Read all genotypes into a GenotypeMatrix.
    ///
    /// This reads the entire .bed file and converts to a matrix.
    /// For very large files, consider using `read_snp_range` instead.
    pub fn read_all(&self) -> Result<GenotypeMatrix> {
        self.read_snp_range(0..self.snps.len())
    }

    /// Read genotypes for a range of SNPs.
    ///
    /// Returns a GenotypeMatrix with only the specified SNPs.
    pub fn read_snp_range(&self, range: std::ops::Range<usize>) -> Result<GenotypeMatrix> {
        crate::check_snp_range(&range, self.snps.len())?;
        // Copy values before mmap borrow
        let n_samples = self.samples.len();
        let n_snps = range.len();
        let bytes_per_snp = self.bytes_per_snp;
        let snps: Vec<Snp> = range.clone().map(|i| self.snps[i].clone()).collect();
        let samples = self.samples.clone();

        let mmap = self.ensure_mmap()?;

        // Allocate output matrix
        let mut data = Array2::zeros((n_snps, n_samples));

        // Read each SNP
        for (out_idx, snp_idx) in range.enumerate() {
            let offset = 3 + snp_idx * bytes_per_snp;
            let snp_bytes = &mmap[offset..offset + bytes_per_snp];

            // Decode genotypes
            for (sample_idx, geno) in decode_snp_genotypes(snp_bytes, n_samples).enumerate() {
                data[[out_idx, sample_idx]] = geno as u8;
            }
        }

        GenotypeMatrix::new(data, samples, snps)
    }

    /// Read genotypes for a single SNP.
    pub fn read_snp(&self, snp_idx: usize) -> Result<Vec<Genotype>> {
        crate::check_snp_range(&(snp_idx..snp_idx + 1), self.snps.len())?;
        // Copy values before mmap borrow
        let n_samples = self.samples.len();
        let bytes_per_snp = self.bytes_per_snp;

        let mmap = self.ensure_mmap()?;

        let offset = 3 + snp_idx * bytes_per_snp;
        let snp_bytes = &mmap[offset..offset + bytes_per_snp];

        Ok(decode_snp_genotypes(snp_bytes, n_samples).collect())
    }

    /// Detect ploidy for each sample based on heterozygosity.
    ///
    /// Reads the first `n_snps_check` SNPs to determine if each sample
    /// has any heterozygous calls.
    pub fn detect_ploidy(&mut self, n_snps_check: usize) -> Result<()> {
        // Copy values before mmap borrow
        let n_samples = self.samples.len();
        let n_check = n_snps_check.min(self.snps.len());
        let bytes_per_snp = self.bytes_per_snp;

        // Compute has_het in a scope where mmap is borrowed
        let has_het = {
            let mmap = self.ensure_mmap()?;
            let mut has_het = vec![false; n_samples];

            for snp_idx in 0..n_check {
                let offset = 3 + snp_idx * bytes_per_snp;
                let snp_bytes = &mmap[offset..offset + bytes_per_snp];

                for (sample_idx, geno) in decode_snp_genotypes(snp_bytes, n_samples).enumerate() {
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

/// Read a .fam file (PLINK sample metadata).
///
/// Format: FID IID PID MID SEX PHENO (space-separated)
/// We use FID as population and IID as sample ID.
fn read_fam_file(path: &str) -> Result<Vec<Sample>> {
    let content = std::fs::read_to_string(path).map_err(|e| PopmixError::io_with_path(e, path))?;

    let mut samples = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(PopmixError::parse(
                path,
                line_num + 1,
                format!("Expected at least 2 fields, got {}", fields.len()),
            ));
        }

        let population = PopulationId::new(fields[0]);
        let id = fields[1].to_string();
        let sex = fields.get(4).map_or(Sex::Unknown, |c| Sex::from_fam(c));

        samples.push(Sample::new(id, population).with_sex(sex));
    }

    Ok(samples)
}

/// Read a .bim file (PLINK SNP metadata).
///
/// Format: CHR SNP_ID GPOS BP A1 A2 (tab or space-separated). The genetic
/// position is in Morgans, or in centimorgans if any exceeds 10 (see
/// [`crate::to_morgans`]).
fn read_bim_file(path: &str) -> Result<Vec<Snp>> {
    let content = std::fs::read_to_string(path).map_err(|e| PopmixError::io_with_path(e, path))?;

    let mut snps = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(PopmixError::parse(
                path,
                line_num + 1,
                format!("Expected 6 fields, got {}", fields.len()),
            ));
        }

        let chromosome = Chromosome::parse(fields[0]);
        let id = fields[1].to_string();
        // Morgans, as EIGENSOFT and admixtools2 write and read it; rescaled
        // below if the file is in centimorgans
        let position_morgans = crate::parse_genetic_position(fields[2], path, line_num + 1)?;
        let position_bp = fields[3].parse::<u32>().map_err(|_| {
            PopmixError::parse(
                path,
                line_num + 1,
                format!("Invalid BP position: {}", fields[3]),
            )
        })?;

        let allele_ref = fields[4].chars().next().unwrap_or('N');
        let allele_alt = fields[5].chars().next().unwrap_or('N');

        snps.push(Snp::new(
            id,
            chromosome,
            position_bp,
            position_morgans,
            allele_ref,
            allele_alt,
        ));
    }

    crate::to_morgans(&mut snps);
    Ok(snps)
}

/// Decode 4 genotypes from a single byte using PLINK encoding.
///
/// PLINK stores genotypes as 2 bits each:
/// - 00 = Homozygous A1 (alternate) -> 2
/// - 01 = Missing -> 3
/// - 10 = Heterozygous -> 1
/// - 11 = Homozygous A2 (reference) -> 0
#[inline]
pub fn decode_plink_byte(byte: u8) -> [Genotype; 4] {
    [
        Genotype::from_plink_bits(byte & 0b11),
        Genotype::from_plink_bits((byte >> 2) & 0b11),
        Genotype::from_plink_bits((byte >> 4) & 0b11),
        Genotype::from_plink_bits((byte >> 6) & 0b11),
    ]
}

/// Iterator over genotypes for a single SNP from raw bytes.
fn decode_snp_genotypes(bytes: &[u8], n_samples: usize) -> impl Iterator<Item = Genotype> + '_ {
    bytes
        .iter()
        .flat_map(|&byte| decode_plink_byte(byte))
        .take(n_samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_genetic_position_that_is_not_a_number_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.bim");
        std::fs::write(&path, "1 rs1 0 100 A G\n1 rs2 x1 200 C T\n").unwrap();
        let err = read_bim_file(path.to_str().unwrap())
            .unwrap_err()
            .to_string();
        assert!(err.contains("line 2") && err.contains("x1"), "{err}");
    }

    #[test]
    fn test_decode_plink_byte() {
        // All HomRef (0b11111111)
        let genos = decode_plink_byte(0b11111111);
        assert!(genos.iter().all(|g| *g == Genotype::HomRef));

        // All HomAlt (0b00000000)
        let genos = decode_plink_byte(0b00000000);
        assert!(genos.iter().all(|g| *g == Genotype::HomAlt));

        // Mixed: HomAlt, Het, Missing, HomRef (reading right to left in bits)
        // Byte: 11_01_10_00 = 0xD8
        let genos = decode_plink_byte(0b11_01_10_00);
        assert_eq!(genos[0], Genotype::HomAlt); // 00
        assert_eq!(genos[1], Genotype::Het); // 10
        assert_eq!(genos[2], Genotype::Missing); // 01
        assert_eq!(genos[3], Genotype::HomRef); // 11
    }

    #[test]
    fn test_bytes_per_snp_calculation() {
        // 1-4 samples = 1 byte
        assert_eq!((1 + 3) / 4, 1);
        assert_eq!((4 + 3) / 4, 1);

        // 5-8 samples = 2 bytes
        assert_eq!((5 + 3) / 4, 2);
        assert_eq!((8 + 3) / 4, 2);

        // 100 samples = 25 bytes
        assert_eq!((100 + 3) / 4, 25);
    }
}
