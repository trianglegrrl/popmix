//! PED format reader.
//!
//! PED format (PLINK text) uses three files:
//! - `.ped`: One line per sample with genotype data
//! - `.map` or `.pedsnp`: SNP metadata
//! - `.pedind`: Sample metadata (optional, first columns of .ped)

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use popmix_core::error::{PopmixError, Result};
use popmix_core::genotype_matrix::GenotypeMatrix;
use popmix_core::genotypes::{Chromosome, Genotype, Sample, Sex, Snp};
use popmix_core::IoResultExt;

/// Reader for PED format files.
pub struct PedReader {
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    ped_path: String,
    /// Whether the .map lists each SNP's alleles; otherwise (PLINK's
    /// 4-column .map) they are taken from the genotypes.
    map_has_alleles: bool,
}

impl PedReader {
    /// Create a new PED reader from a file prefix.
    ///
    /// Looks for files: `{prefix}.ped`, `{prefix}.map` (or `.pedsnp`), and optionally `{prefix}.pedind`
    pub fn new(prefix: impl AsRef<Path>) -> Result<Self> {
        let prefix = prefix.as_ref();
        let prefix_str = prefix.to_string_lossy();

        // Find .ped file
        let ped_path = format!("{prefix_str}.ped");
        if !Path::new(&ped_path).exists() {
            return Err(PopmixError::FileNotFound {
                path: ped_path.into(),
            });
        }

        // Find .map or .pedsnp file
        let map_path = format!("{prefix_str}.map");
        let pedsnp_path = format!("{prefix_str}.pedsnp");
        let snp_path = if Path::new(&map_path).exists() {
            map_path
        } else if Path::new(&pedsnp_path).exists() {
            pedsnp_path
        } else {
            return Err(PopmixError::FileNotFound {
                path: format!("{prefix_str}.map or {prefix_str}.pedsnp").into(),
            });
        };

        // Read SNPs from map file
        let (snps, map_has_alleles) = read_map_file(&snp_path)?;

        // Read samples from .ped file (first 6 columns) or .pedind if present
        let pedind_path = format!("{prefix_str}.pedind");
        let samples = if Path::new(&pedind_path).exists() {
            read_pedind_file(&pedind_path)?
        } else {
            read_samples_from_ped(&ped_path)?
        };

        Ok(PedReader {
            samples,
            snps,
            ped_path,
            map_has_alleles,
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

        // Allele lookup from the .map, or (4-column .map) the first two
        // alleles each SNP's genotypes show, the first one counted
        let mut allele_lookup: Vec<(char, char)> = if self.map_has_alleles {
            self.snps
                .iter()
                .map(|s| (s.allele_ref, s.allele_alt))
                .collect()
        } else {
            vec![('\0', '\0'); n_snps]
        };

        // Initialize genotype matrix (SNP-major)
        let mut data = ndarray::Array2::<u8>::zeros((n_snps, n_samples));

        // Read .ped file
        let file = File::open(&self.ped_path).with_path(&self.ped_path)?;
        let reader = BufReader::new(file);

        for (sample_idx, line) in reader.lines().enumerate() {
            let line = line.with_path(&self.ped_path)?;
            let fields: Vec<&str> = line.split_whitespace().collect();

            // Skip first 6 columns (FID, IID, PID, MID, SEX, PHENO)
            // Or 7 columns if there's an extra column
            let geno_start = if fields.len() == 6 + 2 * n_snps {
                6
            } else if fields.len() == 7 + 2 * n_snps {
                7
            } else {
                return Err(PopmixError::ParseError {
                    path: self.ped_path.clone().into(),
                    line: sample_idx + 1,
                    details: format!(
                        "Expected {} or {} columns, got {}",
                        6 + 2 * n_snps,
                        7 + 2 * n_snps,
                        fields.len()
                    ),
                });
            };

            // Parse genotypes (2 alleles per SNP)
            for snp_idx in 0..n_snps {
                let a1 = fields[geno_start + snp_idx * 2];
                let a2 = fields[geno_start + snp_idx * 2 + 1];

                if !self.map_has_alleles {
                    let seen = &mut allele_lookup[snp_idx];
                    for c in [parse_allele(a1), parse_allele(a2)].into_iter().flatten() {
                        if seen.0 == '\0' {
                            seen.0 = c;
                        } else if seen.1 == '\0' && c != seen.0 {
                            seen.1 = c;
                        }
                    }
                }
                let geno = parse_ped_genotype(a1, a2, allele_lookup[snp_idx]);
                data[[snp_idx, sample_idx]] = geno as u8;
            }
        }

        let mut snps = self.snps.clone();
        if !self.map_has_alleles {
            for (snp, &(r, a)) in snps.iter_mut().zip(&allele_lookup) {
                snp.allele_ref = if r == '\0' { 'N' } else { r };
                snp.allele_alt = if a == '\0' { 'N' } else { a };
            }
        }
        GenotypeMatrix::new(data, self.samples.clone(), snps)
    }
}

/// Parse a PED genotype from two allele strings.
fn parse_ped_genotype(a1: &str, a2: &str, (ref_allele, alt_allele): (char, char)) -> Genotype {
    // Missing data is represented by 0
    if a1 == "0" || a2 == "0" {
        return Genotype::Missing;
    }

    // Parse alleles (can be ACGT or 1234)
    let c1 = parse_allele(a1);
    let c2 = parse_allele(a2);

    if c1.is_none() || c2.is_none() {
        return Genotype::Missing;
    }

    let c1 = c1.unwrap();
    let c2 = c2.unwrap();

    // a third allele (when the alleles are known) is missing data
    if alt_allele != '\0' && [c1, c2].iter().any(|&c| c != ref_allele && c != alt_allele) {
        return Genotype::Missing;
    }

    // Count reference alleles
    let ref_count = i32::from(c1 == ref_allele) + i32::from(c2 == ref_allele);

    match ref_count {
        0 => Genotype::HomAlt,
        1 => Genotype::Het,
        2 => Genotype::HomRef,
        _ => Genotype::Missing,
    }
}

/// Parse an allele string to a character.
fn parse_allele(s: &str) -> Option<char> {
    match s {
        "A" | "a" | "1" => Some('A'),
        "C" | "c" | "2" => Some('C'),
        "G" | "g" | "3" => Some('G'),
        "T" | "t" | "4" => Some('T'),
        "0" => None,
        _ => s.chars().next(),
    }
}

/// Read .map file (PED SNP format).
/// Format: CHR ID CM BP [A1 A2]
fn read_map_file(path: &str) -> Result<(Vec<Snp>, bool)> {
    let file = File::open(path).with_path(path)?;
    let reader = BufReader::new(file);
    let mut snps = Vec::new();
    let mut has_alleles = true;

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.with_path(path)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(PopmixError::ParseError {
                path: path.to_string().into(),
                line: line_num + 1,
                details: format!("Expected at least 4 columns, got {}", fields.len()),
            });
        }

        let chromosome = Chromosome::parse(fields[0]);
        let id = fields[1].to_string();
        let position_morgans: f64 = fields[2].parse().map_err(|_| PopmixError::ParseError {
            path: path.to_string().into(),
            line: line_num + 1,
            details: format!("Invalid genetic position: {}", fields[2]),
        })?;
        let position_bp: u32 = fields[3].parse().map_err(|_| PopmixError::ParseError {
            path: path.to_string().into(),
            line: line_num + 1,
            details: format!("Invalid physical position: {}", fields[3]),
        })?;

        // Alleles are optional in .map (PLINK writes 4 columns); without
        // them they come from the genotypes
        let (allele_ref, allele_alt) = if fields.len() >= 6 {
            (
                fields[4].chars().next().unwrap_or('A'),
                fields[5].chars().next().unwrap_or('G'),
            )
        } else {
            has_alleles = false;
            ('A', 'G')
        };

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
    Ok((snps, has_alleles))
}

/// Read .pedind file.
/// Format: FID IID PID MID SEX PHENO [1]
fn read_pedind_file(path: &str) -> Result<Vec<Sample>> {
    let file = File::open(path).with_path(path)?;
    let reader = BufReader::new(file);
    let mut samples = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.with_path(path)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(PopmixError::ParseError {
                path: path.to_string().into(),
                line: line_num + 1,
                details: format!("Expected at least 6 columns, got {}", fields.len()),
            });
        }

        // FID is used as population, IID as sample ID
        let population = fields[0].to_string();
        let id = fields[1].to_string();

        samples.push(Sample::new(id, population).with_sex(Sex::from_fam(fields[4])));
    }

    Ok(samples)
}

/// Read sample info from first columns of .ped file.
fn read_samples_from_ped(path: &str) -> Result<Vec<Sample>> {
    let file = File::open(path).with_path(path)?;
    let reader = BufReader::new(file);
    let mut samples = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.with_path(path)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(PopmixError::ParseError {
                path: path.to_string().into(),
                line: line_num + 1,
                details: format!("Expected at least 6 columns, got {}", fields.len()),
            });
        }

        // FID is used as population, IID as sample ID
        let population = fields[0].to_string();
        let id = fields[1].to_string();

        samples.push(Sample::new(id, population).with_sex(Sex::from_fam(fields[4])));
    }

    Ok(samples)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_alleles_come_from_the_genotypes_without_them_in_the_map() {
        // PLINK writes a 4-column .map; these SNPs are C/T and G/A
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x");
        std::fs::write(
            format!("{}.map", p.display()),
            "1 s1 0.01 100\n1 s2 0.02 200\n",
        )
        .unwrap();
        std::fs::write(
            format!("{}.ped", p.display()),
            "P a 0 0 0 -9 C C G A\nP b 0 0 0 -9 C T A A\nP c 0 0 0 -9 T T 0 0\nP d 0 0 0 -9 C G A A\n",
        )
        .unwrap();
        let g = PedReader::new(&p).unwrap().read_all().unwrap();
        // first allele seen is counted: C at s1, G at s2
        let s1: Vec<_> = (0..4).map(|i| g.get(0, i)).collect();
        assert_eq!(
            s1,
            vec![
                Genotype::HomRef,
                Genotype::Het,
                Genotype::HomAlt,
                Genotype::Missing
            ]
        );
        let s2: Vec<_> = (0..4).map(|i| g.get(1, i)).collect();
        assert_eq!(
            s2,
            vec![
                Genotype::Het,
                Genotype::HomAlt,
                Genotype::Missing,
                Genotype::HomAlt
            ]
        );
        assert_eq!((g.snps()[0].allele_ref, g.snps()[0].allele_alt), ('C', 'T'));
    }

    use super::*;

    #[test]
    fn test_parse_allele() {
        assert_eq!(parse_allele("A"), Some('A'));
        assert_eq!(parse_allele("1"), Some('A'));
        assert_eq!(parse_allele("C"), Some('C'));
        assert_eq!(parse_allele("2"), Some('C'));
        assert_eq!(parse_allele("G"), Some('G'));
        assert_eq!(parse_allele("3"), Some('G'));
        assert_eq!(parse_allele("T"), Some('T'));
        assert_eq!(parse_allele("4"), Some('T'));
        assert_eq!(parse_allele("0"), None);
    }

    #[test]
    fn test_parse_ped_genotype() {
        // A/G alleles, ref=A, alt=G
        assert_eq!(parse_ped_genotype("A", "A", ('A', 'G')), Genotype::HomRef);
        assert_eq!(parse_ped_genotype("A", "G", ('A', 'G')), Genotype::Het);
        assert_eq!(parse_ped_genotype("G", "A", ('A', 'G')), Genotype::Het);
        assert_eq!(parse_ped_genotype("G", "G", ('A', 'G')), Genotype::HomAlt);
        assert_eq!(parse_ped_genotype("0", "0", ('A', 'G')), Genotype::Missing);
        assert_eq!(parse_ped_genotype("A", "0", ('A', 'G')), Genotype::Missing);

        // Numeric alleles (1=A, 3=G)
        assert_eq!(parse_ped_genotype("1", "1", ('A', 'G')), Genotype::HomRef);
        assert_eq!(parse_ped_genotype("1", "3", ('A', 'G')), Genotype::Het);
        assert_eq!(parse_ped_genotype("3", "3", ('A', 'G')), Genotype::HomAlt);
    }
}
