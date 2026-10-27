//! ANCESTRYMAP format reader.
//!
//! ANCESTRYMAP is a sparse format where only non-missing genotypes are stored.
//! Files:
//! - `.ancestrymapgeno`: Sparse genotype data (SNP_ID, Sample_ID, Count)
//! - `.snp`: SNP metadata (same as EIGENSTRAT)
//! - `.ind`: Sample metadata (same as EIGENSTRAT)

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use popmix_core::error::{PopmixError, Result};
use popmix_core::genotype_matrix::GenotypeMatrix;
use popmix_core::genotypes::{Genotype, Sample, Snp};
use popmix_core::IoResultExt;

use crate::eigenstrat::{read_ind_file, read_snp_file};

/// Reader for ANCESTRYMAP format files.
pub struct AncestryMapReader {
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    geno_path: String,
    sample_index: HashMap<String, usize>,
    snp_index: HashMap<String, usize>,
}

impl AncestryMapReader {
    /// Create a new ANCESTRYMAP reader from a file prefix.
    ///
    /// Looks for files: `{prefix}.ancestrymapgeno`, `{prefix}.snp`, `{prefix}.ind`
    pub fn new(prefix: impl AsRef<Path>) -> Result<Self> {
        let prefix = prefix.as_ref();
        let prefix_str = prefix.to_string_lossy();

        // Check for .ancestrymapgeno file
        let geno_path = format!("{}.ancestrymapgeno", prefix_str);
        if !Path::new(&geno_path).exists() {
            return Err(PopmixError::FileNotFound {
                path: geno_path.into(),
            });
        }

        // Read SNP and sample metadata (same format as EIGENSTRAT)
        let snp_path = format!("{}.snp", prefix_str);
        let ind_path = format!("{}.ind", prefix_str);

        let snps = read_snp_file(&snp_path)?;
        let samples = read_ind_file(&ind_path)?;

        // Build index maps for fast lookup
        let sample_index: HashMap<String, usize> = samples
            .iter()
            .enumerate()
            .map(|(i, s)| (s.id.clone(), i))
            .collect();

        let snp_index: HashMap<String, usize> = snps
            .iter()
            .enumerate()
            .map(|(i, s)| (s.id.clone(), i))
            .collect();

        Ok(AncestryMapReader {
            samples,
            snps,
            geno_path,
            sample_index,
            snp_index,
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

    /// Read all genotypes into a GenotypeMatrix.
    pub fn read_all(&self) -> Result<GenotypeMatrix> {
        let n_samples = self.samples.len();
        let n_snps = self.snps.len();

        // Initialize with missing (3)
        let mut data = ndarray::Array2::<u8>::from_elem((n_snps, n_samples), 3);

        // Read sparse genotypes
        let file = File::open(&self.geno_path).with_path(&self.geno_path)?;
        let reader = BufReader::new(file);

        for (line_num, line) in reader.lines().enumerate() {
            let line = line.with_path(&self.geno_path)?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 3 {
                return Err(PopmixError::ParseError {
                    path: self.geno_path.clone().into(),
                    line: line_num + 1,
                    details: format!("Expected 3 columns, got {}", fields.len()),
                });
            }

            let snp_id = fields[0];
            let sample_id = fields[1];
            let count: u8 = fields[2].parse().map_err(|_| PopmixError::ParseError {
                path: self.geno_path.clone().into(),
                line: line_num + 1,
                details: format!("Invalid genotype count: {}", fields[2]),
            })?;

            // Look up indices
            let snp_idx = self
                .snp_index
                .get(snp_id)
                .ok_or_else(|| PopmixError::ParseError {
                    path: self.geno_path.clone().into(),
                    line: line_num + 1,
                    details: format!("Unknown SNP ID: {}", snp_id),
                })?;

            let sample_idx =
                self.sample_index
                    .get(sample_id)
                    .ok_or_else(|| PopmixError::ParseError {
                        path: self.geno_path.clone().into(),
                        line: line_num + 1,
                        details: format!("Unknown sample ID: {}", sample_id),
                    })?;

            // ANCESTRYMAP count is number of reference alleles (0, 1, or 2)
            // Our encoding: 0=HomRef, 1=Het, 2=HomAlt
            // So count=2 -> HomRef, count=1 -> Het, count=0 -> HomAlt
            let geno = match count {
                2 => Genotype::HomRef,
                1 => Genotype::Het,
                0 => Genotype::HomAlt,
                _ => Genotype::Missing,
            };

            data[[*snp_idx, *sample_idx]] = geno as u8;
        }

        GenotypeMatrix::new(data, self.samples.clone(), self.snps.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn write_test_files(dir: &TempDir) -> String {
        let prefix = dir.path().join("test");
        let prefix_str = prefix.to_string_lossy().to_string();

        // Write .ind file
        let mut ind = File::create(format!("{}.ind", prefix_str)).unwrap();
        writeln!(ind, "Sample1 M Pop1").unwrap();
        writeln!(ind, "Sample2 F Pop2").unwrap();

        // Write .snp file
        let mut snp = File::create(format!("{}.snp", prefix_str)).unwrap();
        writeln!(snp, "rs1 1 0.001 100 A G").unwrap();
        writeln!(snp, "rs2 1 0.002 200 C T").unwrap();

        // Write .ancestrymapgeno file
        let mut geno = File::create(format!("{}.ancestrymapgeno", prefix_str)).unwrap();
        writeln!(geno, "rs1 Sample1 2").unwrap(); // HomRef
        writeln!(geno, "rs1 Sample2 1").unwrap(); // Het
        writeln!(geno, "rs2 Sample1 0").unwrap(); // HomAlt
                                                  // rs2 Sample2 is missing (no entry)

        prefix_str
    }

    #[test]
    fn test_ancestrymap_reader() {
        let tmp = TempDir::new().unwrap();
        let prefix = write_test_files(&tmp);

        let reader = AncestryMapReader::new(&prefix).unwrap();
        assert_eq!(reader.num_samples(), 2);
        assert_eq!(reader.num_snps(), 2);

        let matrix = reader.read_all().unwrap();

        // rs1: Sample1=HomRef(0), Sample2=Het(1)
        assert_eq!(matrix.get(0, 0), Genotype::HomRef);
        assert_eq!(matrix.get(0, 1), Genotype::Het);

        // rs2: Sample1=HomAlt(2), Sample2=Missing(3)
        assert_eq!(matrix.get(1, 0), Genotype::HomAlt);
        assert_eq!(matrix.get(1, 1), Genotype::Missing);
    }
}
