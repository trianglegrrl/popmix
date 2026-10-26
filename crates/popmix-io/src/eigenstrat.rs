//! EIGENSTRAT text format (.geno, .snp, .ind) reader.
//!
//! EIGENSTRAT is a text-based format where:
//! - `.geno`: Genotypes as single digits (0, 1, 2, 9 for missing), one SNP per line
//! - `.snp`: SNP metadata (one line per SNP)
//! - `.ind`: Sample metadata (one line per sample)

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use ndarray::Array2;

use popmix_core::{
    Chromosome, Genotype, GenotypeMatrix, Ploidy, PopmixError, PopulationId, Result, Sample, Sex,
    Snp,
};

/// Reader for EIGENSTRAT text format files.
pub struct EigenstratReader {
    prefix: String,
    samples: Vec<Sample>,
    snps: Vec<Snp>,
}

impl EigenstratReader {
    /// Open an EIGENSTRAT dataset by prefix.
    pub fn new(prefix: impl AsRef<Path>) -> Result<Self> {
        let prefix = prefix.as_ref().to_string_lossy().to_string();

        // Read .ind file (samples)
        let ind_path = format!("{}.ind", prefix);
        let samples = read_ind_file(&ind_path)?;

        // Read .snp file (SNPs)
        let snp_path = format!("{}.snp", prefix);
        let snps = read_snp_file(&snp_path)?;

        // Verify .geno file exists
        let geno_path = format!("{}.geno", prefix);
        if !Path::new(&geno_path).exists() {
            return Err(PopmixError::FileNotFound { path: geno_path });
        }

        Ok(EigenstratReader {
            prefix,
            samples,
            snps,
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

    /// Read all genotypes into a GenotypeMatrix.
    ///
    /// This reads the entire .geno file and converts to a matrix.
    /// For very large files, consider using `read_snp_range` instead.
    pub fn read_all(&self) -> Result<GenotypeMatrix> {
        self.read_snp_range(0..self.snps.len())
    }

    /// Read genotypes for a range of SNPs.
    ///
    /// Returns a GenotypeMatrix with only the specified SNPs.
    pub fn read_snp_range(&self, range: std::ops::Range<usize>) -> Result<GenotypeMatrix> {
        crate::check_snp_range(&range, self.snps.len())?;
        let geno_path = format!("{}.geno", self.prefix);
        let file = File::open(&geno_path).map_err(|e| PopmixError::io_with_path(e, &geno_path))?;
        let reader = BufReader::new(file);

        let n_samples = self.samples.len();
        let n_snps = range.len();

        // Allocate output matrix
        let mut data = Array2::zeros((n_snps, n_samples));

        // Read through the file
        let mut out_idx = 0;
        for (line_num, line_result) in reader.lines().enumerate() {
            // Skip lines before range
            if line_num < range.start {
                continue;
            }
            // Stop after range
            if line_num >= range.end {
                break;
            }

            let line = line_result.map_err(|e| PopmixError::io_with_path(e, &geno_path))?;
            let line = line.trim();

            // Verify line length matches expected samples (using bytes since we're ASCII)
            let line_bytes = line.as_bytes();
            if line_bytes.len() != n_samples {
                return Err(PopmixError::parse(
                    &geno_path,
                    line_num + 1,
                    format!("Expected {} genotypes, got {}", n_samples, line_bytes.len()),
                ));
            }

            // Parse genotypes
            for (sample_idx, &byte) in line_bytes.iter().enumerate() {
                let geno = Genotype::from_eigenstrat_char(byte).unwrap_or(Genotype::Missing);
                data[[out_idx, sample_idx]] = geno as u8;
            }

            out_idx += 1;
        }

        // Verify we read the expected number of SNPs
        if out_idx != n_snps {
            return Err(PopmixError::parse(
                &geno_path,
                0,
                format!(
                    "Expected {} SNPs in range, but only read {}",
                    n_snps, out_idx
                ),
            ));
        }

        // Subset SNP metadata
        let snps: Vec<Snp> = range.map(|i| self.snps[i].clone()).collect();

        GenotypeMatrix::new(data, self.samples.clone(), snps)
    }

    /// Read genotypes for a single SNP.
    pub fn read_snp(&self, snp_idx: usize) -> Result<Vec<Genotype>> {
        let geno_path = format!("{}.geno", self.prefix);
        let file = File::open(&geno_path).map_err(|e| PopmixError::io_with_path(e, &geno_path))?;
        let reader = BufReader::new(file);

        let n_samples = self.samples.len();

        // Skip to the requested line
        for (line_num, line_result) in reader.lines().enumerate() {
            if line_num != snp_idx {
                continue;
            }

            let line = line_result.map_err(|e| PopmixError::io_with_path(e, &geno_path))?;
            let line = line.trim();

            let line_bytes = line.as_bytes();
            if line_bytes.len() != n_samples {
                return Err(PopmixError::parse(
                    &geno_path,
                    line_num + 1,
                    format!("Expected {} genotypes, got {}", n_samples, line_bytes.len()),
                ));
            }

            return Ok(line_bytes
                .iter()
                .map(|&b| Genotype::from_eigenstrat_char(b).unwrap_or(Genotype::Missing))
                .collect());
        }

        Err(PopmixError::parse(
            &geno_path,
            0,
            format!("SNP index {} out of range", snp_idx),
        ))
    }

    /// Detect ploidy for each sample based on heterozygosity.
    ///
    /// Reads the first `n_snps_check` SNPs to determine if each sample
    /// has any heterozygous calls.
    pub fn detect_ploidy(&mut self, n_snps_check: usize) -> Result<()> {
        let n_samples = self.samples.len();
        let n_check = n_snps_check.min(self.snps.len());

        let geno_path = format!("{}.geno", self.prefix);
        let file = File::open(&geno_path).map_err(|e| PopmixError::io_with_path(e, &geno_path))?;
        let reader = BufReader::new(file);

        let mut has_het = vec![false; n_samples];

        for (line_num, line_result) in reader.lines().enumerate() {
            if line_num >= n_check {
                break;
            }

            let line = line_result.map_err(|e| PopmixError::io_with_path(e, &geno_path))?;
            let line = line.trim();

            for (sample_idx, byte) in line.bytes().enumerate() {
                if sample_idx >= n_samples {
                    break;
                }
                if Genotype::from_eigenstrat_char(byte) == Some(Genotype::Het) {
                    has_het[sample_idx] = true;
                }
            }
        }

        // Update ploidy
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

/// Read a .ind file (EIGENSTRAT sample metadata).
///
/// Format: SAMPLE_ID SEX POPULATION (whitespace-separated)
pub fn read_ind_file(path: &str) -> Result<Vec<Sample>> {
    let content = std::fs::read_to_string(path).map_err(|e| PopmixError::io_with_path(e, path))?;

    let mut samples = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            return Err(PopmixError::parse(
                path,
                line_num + 1,
                format!("Expected at least 3 fields, got {}", fields.len()),
            ));
        }

        let id = fields[0].to_string();
        let population = PopulationId::new(fields[2]);

        samples.push(
            Sample::with_ploidy(id, population, Ploidy::Diploid).with_sex(Sex::from_ind(fields[1])),
        );
    }

    Ok(samples)
}

/// Read a .snp file (EIGENSTRAT SNP metadata).
///
/// Format: SNP_ID CHR CM BP A1 A2 (whitespace-separated)
/// Note: Column order differs from PLINK .bim (SNP_ID is first, CHR is second)
pub fn read_snp_file(path: &str) -> Result<Vec<Snp>> {
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

        let id = fields[0].to_string();
        let chromosome = Chromosome::parse(fields[1]);
        // SNP files store genetic position in Morgans - keep in Morgans to match DATES
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

    Ok(snps)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ind_line() {
        // Test that the module compiles and basic types work
        let sample = Sample::new("test", PopulationId::new("Pop1"));
        assert_eq!(sample.id, "test");
    }

    #[test]
    fn test_genotype_from_eigenstrat() {
        assert_eq!(Genotype::from_eigenstrat_char(b'0'), Some(Genotype::HomRef));
        assert_eq!(Genotype::from_eigenstrat_char(b'1'), Some(Genotype::Het));
        assert_eq!(Genotype::from_eigenstrat_char(b'2'), Some(Genotype::HomAlt));
        assert_eq!(
            Genotype::from_eigenstrat_char(b'9'),
            Some(Genotype::Missing)
        );
        assert_eq!(Genotype::from_eigenstrat_char(b'X'), None);
    }

    #[test]
    fn a_genetic_position_that_is_not_a_number_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        for bad in ["abc", "nan", "inf"] {
            let path = dir.path().join("t.snp");
            std::fs::write(&path, format!("rs1 1 0.01 100 A G\nrs2 1 {bad} 200 C T\n")).unwrap();
            let err = read_snp_file(path.to_str().unwrap())
                .unwrap_err()
                .to_string();
            assert!(err.contains("line 2") && err.contains(bad), "{bad}: {err}");
        }
    }

    #[test]
    fn a_snp_range_outside_the_file_is_an_error() {
        let fixture = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/fixtures/small_eigenstrat"
        );
        let reader = EigenstratReader::new(fixture).unwrap();
        let n = reader.snps().len();
        assert_eq!(reader.read_snp_range(1..n).unwrap().n_snps(), n - 1);
        for range in [
            0..n + 1,
            n + 5..n + 10,
            std::ops::Range { start: 3, end: 2 },
        ] {
            let err = reader
                .read_snp_range(range.clone())
                .unwrap_err()
                .to_string();
            assert!(err.contains("SNP range"), "{range:?}: {err}");
        }
    }
}
