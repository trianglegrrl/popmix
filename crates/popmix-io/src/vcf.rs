//! VCF (Variant Call Format) reader.
//!
//! VCF is a text format for storing genetic variants.
//! This module reads VCF files and converts them to genotype matrices.
//!
//! Reference: VCF v4.2 specification
//! Based on PLINK 2.0 VCF parser concepts by Shaun Purcell & Christopher Chang.
//! Reference: https://github.com/chrchang/plink-ng

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;
use ndarray::Array2;

use popmix_core::{
    Chromosome, Genotype, GenotypeMatrix, Ploidy, PopmixError, PopulationId, Result, Sample, Snp,
};

/// Reader for VCF (Variant Call Format) files.
///
/// Supports both plain text (.vcf) and gzip-compressed (.vcf.gz) files.
pub struct VcfReader {
    path: String,
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    /// Whether the file is gzip compressed
    is_gzipped: bool,
}

impl VcfReader {
    /// Open a VCF file for reading.
    ///
    /// The path can be either:
    /// - A plain text VCF file (.vcf)
    /// - A gzip-compressed VCF file (.vcf.gz)
    pub fn new(path: impl AsRef<Path>) -> Result<Self> {
        let path_str = path.as_ref().to_string_lossy().to_string();
        let is_gzipped = path_str.ends_with(".gz");

        // Parse header to get samples and SNP count
        let (samples, snps) = parse_vcf_header_and_variants(&path_str, is_gzipped)?;

        Ok(VcfReader {
            path: path_str,
            samples,
            snps,
            is_gzipped,
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
        read_vcf_genotypes(&self.path, self.is_gzipped, &self.samples, &self.snps)
    }
}

/// Parse a VCF GT field into a Genotype.
///
/// Supports both phased (|) and unphased (/) separators.
///
/// VCF GT encoding:
/// - 0/0 or 0|0 = Homozygous reference
/// - 0/1 or 0|1 or 1/0 or 1|0 = Heterozygous
/// - 1/1 or 1|1 = Homozygous alternate
/// - ./. or .|. = Missing
///
/// Returns None if the GT field cannot be parsed.
pub fn parse_vcf_gt(gt: &str) -> Option<Genotype> {
    // Handle missing
    if gt.starts_with('.') {
        return Some(Genotype::Missing);
    }

    // Split on / or |
    let parts: Vec<&str> = if gt.contains('|') {
        gt.split('|').collect()
    } else if gt.contains('/') {
        gt.split('/').collect()
    } else {
        return None;
    };

    if parts.len() != 2 {
        return None;
    }

    // Parse allele indices
    let a1 = parts[0].parse::<u8>().ok()?;
    let a2 = parts[1].parse::<u8>().ok()?;

    // Convert to genotype (count of ALT allele)
    let dosage = a1 + a2;
    match dosage {
        0 => Some(Genotype::HomRef), // 0/0
        1 => Some(Genotype::Het),    // 0/1 or 1/0
        2 => Some(Genotype::HomAlt), // 1/1
        _ => None,                   // Multi-allelic not supported
    }
}

/// Parse VCF header and count variants.
fn parse_vcf_header_and_variants(path: &str, is_gzipped: bool) -> Result<(Vec<Sample>, Vec<Snp>)> {
    let file = File::open(path).map_err(|e| PopmixError::io_with_path(e, path))?;

    let reader: Box<dyn BufRead> = if is_gzipped {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut samples = Vec::new();
    let mut snps = Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| PopmixError::io_with_path(e, path))?;

        // Skip meta-information lines (##)
        if line.starts_with("##") {
            continue;
        }

        // Parse header line (#CHROM)
        if line.starts_with("#CHROM") {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 10 {
                return Err(PopmixError::parse(
                    path,
                    line_num + 1,
                    "VCF header must have at least 10 columns (including sample columns)",
                ));
            }

            // Sample names start at column 9 (0-indexed)
            for sample_name in &fields[9..] {
                // VCF doesn't have population info - use sample ID as population
                samples.push(Sample::with_ploidy(
                    sample_name.to_string(),
                    PopulationId::new(sample_name.to_string()),
                    Ploidy::Diploid,
                ));
            }
            continue;
        }

        // Parse data lines
        if !line.starts_with('#') && !line.is_empty() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 8 {
                return Err(PopmixError::parse(
                    path,
                    line_num + 1,
                    format!(
                        "VCF data line must have at least 8 columns, got {}",
                        fields.len()
                    ),
                ));
            }

            let chrom = Chromosome::parse(fields[0]);
            let pos = fields[1].parse::<u32>().map_err(|_| {
                PopmixError::parse(
                    path,
                    line_num + 1,
                    format!("Invalid position: {}", fields[1]),
                )
            })?;
            let id = if fields[2] == "." {
                format!("{}:{}", fields[0], pos)
            } else {
                fields[2].to_string()
            };
            let ref_allele = fields[3].chars().next().unwrap_or('N');
            let alt_allele = fields[4].chars().next().unwrap_or('N');

            // VCF doesn't have genetic position - use 0.0
            snps.push(Snp::new(id, chrom, pos, 0.0, ref_allele, alt_allele));
        }
    }

    if samples.is_empty() {
        return Err(PopmixError::parse(
            path,
            0,
            "No samples found in VCF header",
        ));
    }

    Ok((samples, snps))
}

/// Read genotypes from VCF file.
fn read_vcf_genotypes(
    path: &str,
    is_gzipped: bool,
    samples: &[Sample],
    snps: &[Snp],
) -> Result<GenotypeMatrix> {
    let file = File::open(path).map_err(|e| PopmixError::io_with_path(e, path))?;

    let reader: Box<dyn BufRead> = if is_gzipped {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let n_samples = samples.len();
    let n_snps = snps.len();
    let mut data = Array2::zeros((n_snps, n_samples));

    let mut snp_idx = 0;

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| PopmixError::io_with_path(e, path))?;

        // Skip header lines
        if line.starts_with('#') {
            continue;
        }

        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 + n_samples {
            return Err(PopmixError::parse(
                path,
                line_num + 1,
                format!("Expected {} columns, got {}", 9 + n_samples, fields.len()),
            ));
        }

        // Find GT position in FORMAT field
        let format_field = fields[8];
        let gt_idx = format_field
            .split(':')
            .position(|f| f == "GT")
            .ok_or_else(|| {
                PopmixError::parse(path, line_num + 1, "No GT field in FORMAT column")
            })?;

        // Parse genotypes for each sample
        for (sample_idx, sample_field) in fields[9..].iter().enumerate() {
            if sample_idx >= n_samples {
                break;
            }

            let sample_parts: Vec<&str> = sample_field.split(':').collect();
            let gt_str = sample_parts.get(gt_idx).unwrap_or(&"./.");

            let geno = parse_vcf_gt(gt_str).unwrap_or(Genotype::Missing);
            data[[snp_idx, sample_idx]] = geno as u8;
        }

        snp_idx += 1;
    }

    if snp_idx != n_snps {
        return Err(PopmixError::parse(
            path,
            0,
            format!("Expected {} SNPs, read {}", n_snps, snp_idx),
        ));
    }

    GenotypeMatrix::new(data, samples.to_vec(), snps.to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Fixture paths relative to workspace root
    const VCF_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_vcf.vcf"
    );
    const VCF_MISSING_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_vcf_missing.vcf"
    );
    const VCF_PHASED_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_vcf_phased.vcf"
    );
    const EIGENSTRAT_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_eigenstrat"
    );

    // TDD RED: Test GT parsing
    #[test]
    fn test_parse_vcf_gt_homref() {
        assert_eq!(parse_vcf_gt("0/0"), Some(Genotype::HomRef));
        assert_eq!(parse_vcf_gt("0|0"), Some(Genotype::HomRef));
    }

    #[test]
    fn test_parse_vcf_gt_het() {
        assert_eq!(parse_vcf_gt("0/1"), Some(Genotype::Het));
        assert_eq!(parse_vcf_gt("1/0"), Some(Genotype::Het));
        assert_eq!(parse_vcf_gt("0|1"), Some(Genotype::Het));
        assert_eq!(parse_vcf_gt("1|0"), Some(Genotype::Het));
    }

    #[test]
    fn test_parse_vcf_gt_homalt() {
        assert_eq!(parse_vcf_gt("1/1"), Some(Genotype::HomAlt));
        assert_eq!(parse_vcf_gt("1|1"), Some(Genotype::HomAlt));
    }

    #[test]
    fn test_parse_vcf_gt_missing() {
        assert_eq!(parse_vcf_gt("./."), Some(Genotype::Missing));
        assert_eq!(parse_vcf_gt(".|."), Some(Genotype::Missing));
        assert_eq!(parse_vcf_gt("."), Some(Genotype::Missing));
    }

    #[test]
    fn test_parse_vcf_gt_invalid() {
        assert_eq!(parse_vcf_gt(""), None);
        assert_eq!(parse_vcf_gt("invalid"), None);
    }

    // TDD RED: Test VCF reader with fixture
    #[test]
    fn test_vcf_reader_samples() {
        let reader = VcfReader::new(VCF_FIXTURE).unwrap();
        assert_eq!(reader.num_samples(), 3);

        let sample_ids: Vec<&str> = reader.samples().iter().map(|s| s.id.as_str()).collect();
        assert_eq!(sample_ids, vec!["IND1", "IND2", "IND3"]);
    }

    #[test]
    fn test_vcf_reader_snps() {
        let reader = VcfReader::new(VCF_FIXTURE).unwrap();
        assert_eq!(reader.num_snps(), 5);

        // Check first SNP
        let snp = &reader.snps()[0];
        assert_eq!(snp.id, "rs1");
        assert_eq!(snp.chromosome, Chromosome::Autosome(1));
        assert_eq!(snp.position_bp, 100);
        assert_eq!(snp.allele_ref, 'A');
        assert_eq!(snp.allele_alt, 'G');
    }

    #[test]
    fn test_vcf_reader_genotypes() {
        let reader = VcfReader::new(VCF_FIXTURE).unwrap();
        let matrix = reader.read_all().unwrap();

        // Check dimensions
        assert_eq!(matrix.n_snps(), 5);
        assert_eq!(matrix.n_samples(), 3);

        // Verify genotypes match EIGENSTRAT fixture
        // EIGENSTRAT .geno row 1: 012 → IND1=HomRef, IND2=Het, IND3=HomAlt
        // VCF row 1: IND1=0/0 (HomRef), IND2=0/1 (Het), IND3=1/1 (HomAlt)
        assert_eq!(matrix.get(0, 0), Genotype::HomRef); // IND1, rs1
        assert_eq!(matrix.get(0, 1), Genotype::Het); // IND2, rs1
        assert_eq!(matrix.get(0, 2), Genotype::HomAlt); // IND3, rs1

        // EIGENSTRAT .geno row 2: 120 → IND1=Het, IND2=HomAlt, IND3=HomRef
        // VCF row 2: IND1=0/1 (Het), IND2=1/1 (HomAlt), IND3=0/0 (HomRef)
        assert_eq!(matrix.get(1, 0), Genotype::Het); // IND1, rs2
        assert_eq!(matrix.get(1, 1), Genotype::HomAlt); // IND2, rs2
        assert_eq!(matrix.get(1, 2), Genotype::HomRef); // IND3, rs2
    }

    #[test]
    fn test_vcf_reader_matches_eigenstrat() {
        // Read VCF
        let vcf_reader = VcfReader::new(VCF_FIXTURE).unwrap();
        let vcf_matrix = vcf_reader.read_all().unwrap();

        // Read EIGENSTRAT
        let eigenstrat_reader =
            crate::eigenstrat::EigenstratReader::new(EIGENSTRAT_FIXTURE).unwrap();
        let eigenstrat_matrix = eigenstrat_reader.read_all().unwrap();

        // Verify dimensions match
        assert_eq!(vcf_matrix.n_snps(), eigenstrat_matrix.n_snps());
        assert_eq!(vcf_matrix.n_samples(), eigenstrat_matrix.n_samples());

        // Both VCF and EIGENSTRAT use ALT allele count convention in our implementation:
        // VCF 0/0 (HomRef=0) ↔ EIGENSTRAT '0' (HomRef=0)
        // VCF 0/1 (Het=1) ↔ EIGENSTRAT '1' (Het=1)
        // VCF 1/1 (HomAlt=2) ↔ EIGENSTRAT '2' (HomAlt=2)
        //
        // However, the EIGENSTRAT genotype values in from_eigenstrat_char treat
        // the digit as "count of reference allele", which is INVERTED from VCF.
        //
        // EIGENSTRAT '0' → 0 copies of REF → HomAlt in VCF terms
        // EIGENSTRAT '2' → 2 copies of REF → HomRef in VCF terms
        //
        // BUT the existing from_eigenstrat_char maps:
        //   '0' → HomRef (treating 0 as "homozygous reference" directly)
        //   '2' → HomAlt (treating 2 as "homozygous alternate")
        //
        // So the existing code interprets EIGENSTRAT digits as ALT allele count,
        // which matches VCF convention. Our test data was created to match.
        for snp_idx in 0..vcf_matrix.n_snps() {
            for sample_idx in 0..vcf_matrix.n_samples() {
                let vcf_geno = vcf_matrix.get(snp_idx, sample_idx);
                let eig_geno = eigenstrat_matrix.get(snp_idx, sample_idx);
                assert_eq!(
                    vcf_geno, eig_geno,
                    "Mismatch at SNP {}, sample {}: VCF={:?}, EIGENSTRAT={:?}",
                    snp_idx, sample_idx, vcf_geno, eig_geno
                );
            }
        }
    }

    #[test]
    fn test_vcf_reader_missing_data() {
        let reader = VcfReader::new(VCF_MISSING_FIXTURE).unwrap();
        let matrix = reader.read_all().unwrap();

        // Check missing genotypes
        assert_eq!(matrix.get(0, 0), Genotype::Missing); // IND1, rs1
        assert_eq!(matrix.get(1, 1), Genotype::Missing); // IND2, rs2
        assert_eq!(matrix.get(2, 2), Genotype::Missing); // IND3, rs3
    }

    #[test]
    fn test_vcf_reader_phased() {
        let reader = VcfReader::new(VCF_PHASED_FIXTURE).unwrap();
        let matrix = reader.read_all().unwrap();

        // Phased genotypes should be parsed correctly (same data as unphased, with | instead of /)
        assert_eq!(matrix.get(0, 0), Genotype::HomRef); // 0|0
        assert_eq!(matrix.get(0, 1), Genotype::Het); // 0|1
        assert_eq!(matrix.get(0, 2), Genotype::HomAlt); // 1|1
    }
}
