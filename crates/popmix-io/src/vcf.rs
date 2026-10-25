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

use flate2::read::MultiGzDecoder;
use ndarray::Array2;

use popmix_core::{
    Chromosome, Genotype, GenotypeMatrix, Ploidy, PopmixError, PopulationId, Result, Sample, Snp,
};

/// Reader for VCF (Variant Call Format) files.
///
/// Supports both plain text (.vcf) and gzip-compressed (.vcf.gz) files.
/// Only SNPs with at most two alleles are read (a single-base REF, and a
/// single-base ALT or none, `.`); other records are skipped and counted
/// ([`VcfReader::skipped`]), as with PLINK 2's `--snps-only just-acgt
/// --max-alleles 2`. A site with no ALT is kept, monomorphic, with `.` as
/// its counted allele, as PLINK 2 and convertf write it.
pub struct VcfReader {
    path: String,
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    skipped: VcfSkipped,
    /// Whether the file is gzip compressed
    is_gzipped: bool,
}

/// Records a VCF reader skipped, by reason.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct VcfSkipped {
    /// More than one ALT allele.
    pub multiallelic: usize,
    /// An indel or a symbolic allele (not one base each).
    pub not_snp: usize,
}

impl VcfSkipped {
    /// Records skipped in all.
    pub fn total(&self) -> usize {
        self.multiallelic + self.not_snp
    }
}

/// Why a record is not a biallelic SNP, or None if it is one.
fn skip_reason(ref_allele: &str, alt: &str) -> Option<fn(&mut VcfSkipped)> {
    let base = |a: &str| {
        a.len() == 1
            && matches!(
                a.as_bytes()[0].to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T'
            )
    };
    if alt.contains(',') {
        Some(|k| k.multiallelic += 1)
    } else if !base(ref_allele) || !(base(alt) || alt == ".") {
        Some(|k| k.not_snp += 1)
    } else {
        None
    }
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
        let (samples, snps, skipped) = parse_vcf_header_and_variants(&path_str, is_gzipped)?;

        Ok(VcfReader {
            path: path_str,
            samples,
            snps,
            skipped,
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

    /// Records skipped because they are not biallelic SNPs.
    pub fn skipped(&self) -> VcfSkipped {
        self.skipped
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
        // A haploid call (male chrX, say) is homozygous, as plink2 reads it
        return match gt {
            "0" => Some(Genotype::HomRef),
            "1" => Some(Genotype::HomAlt),
            _ => None,
        };
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
fn parse_vcf_header_and_variants(
    path: &str,
    is_gzipped: bool,
) -> Result<(Vec<Sample>, Vec<Snp>, VcfSkipped)> {
    let file = File::open(path).map_err(|e| PopmixError::io_with_path(e, path))?;

    let reader: Box<dyn BufRead> = if is_gzipped {
        // bgzip (the usual .vcf.gz) is a series of gzip members
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut samples = Vec::new();
    let mut snps = Vec::new();
    let mut skipped = VcfSkipped::default();

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
            if let Some(count) = skip_reason(fields[3], fields[4]) {
                count(&mut skipped);
                continue;
            }
            let ref_allele = fields[3].chars().next().unwrap_or('N');
            let alt_allele = fields[4].chars().next().unwrap_or('N');

            // VCF doesn't have genetic position - use 0.0. Genotypes count ALT
            // copies, and the counted allele is listed first (EIGENSTRAT .snp
            // column 5, PLINK A1), as PLINK 2 and convertf list it.
            snps.push(Snp::new(id, chrom, pos, 0.0, alt_allele, ref_allele));
        }
    }

    if samples.is_empty() {
        return Err(PopmixError::parse(
            path,
            0,
            "No samples found in VCF header",
        ));
    }

    Ok((samples, snps, skipped))
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
        // bgzip (the usual .vcf.gz) is a series of gzip members
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
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

        // the records the header pass skipped
        if skip_reason(fields[3], fields[4]).is_some() {
            continue;
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

    #[test]
    fn a_haploid_call_is_homozygous() {
        // as plink2 reads one (male chrX, for instance)
        assert_eq!(parse_vcf_gt("0"), Some(Genotype::HomRef));
        assert_eq!(parse_vcf_gt("1"), Some(Genotype::HomAlt));
        assert_eq!(parse_vcf_gt("2"), None);
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
    fn test_bgzipped_vcf_reads_every_member() {
        // bgzip writes a series of gzip members; split the fixture across two
        use flate2::{write::GzEncoder, Compression};
        use std::io::Write;
        let text = std::fs::read_to_string(VCF_FIXTURE).unwrap();
        let cut = text
            .match_indices('\n')
            .nth(text.lines().count() - 3)
            .unwrap()
            .0
            + 1;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.vcf.gz");
        let mut bytes = Vec::new();
        for part in [&text[..cut], &text[cut..]] {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(part.as_bytes()).unwrap();
            bytes.extend(enc.finish().unwrap());
        }
        std::fs::write(&path, bytes).unwrap();
        let reader = VcfReader::new(path.to_str().unwrap()).unwrap();
        assert_eq!(reader.num_snps(), 5);
        let a = reader.read_all().unwrap();
        let b = VcfReader::new(VCF_FIXTURE).unwrap().read_all().unwrap();
        assert_eq!(a.n_snps(), b.n_snps());
        assert!((0..a.n_snps()).all(|j| a.snp_raw(j) == b.snp_raw(j)));
    }

    const VCF_MULTIALLELIC_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_vcf_multiallelic.vcf"
    );

    #[test]
    fn only_biallelic_snps_are_read() {
        // rs2 has two ALT alleles, rs3 and rs7 are indels, rs5 is a
        // symbolic allele: skipped and counted, not miscoded. rs4 has no ALT:
        // kept, monomorphic, as PLINK 2 keeps it
        let reader = VcfReader::new(VCF_MULTIALLELIC_FIXTURE).unwrap();
        let ids: Vec<&str> = reader.snps().iter().map(|s| s.id.as_str()).collect();
        assert_eq!(ids, ["rs1", "rs4", "rs6"]);
        assert_eq!(
            reader.skipped(),
            VcfSkipped {
                multiallelic: 1,
                not_snp: 3
            }
        );
        let m = reader.read_all().unwrap();
        assert_eq!(m.n_snps(), 3);
        assert_eq!(reader.snps()[1].allele_ref, '.');
        // rs6: C>T with 1/1, 0/1, 0/0, counted on ALT
        assert_eq!(m.get(2, 0), Genotype::HomAlt);
        assert_eq!(m.get(2, 1), Genotype::Het);
        assert_eq!(m.get(2, 2), Genotype::HomRef);
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
        // REF A, ALT G: the counted allele (listed first, as PLINK 2 and
        // convertf list it after importing a VCF) is ALT
        assert_eq!(snp.allele_ref, 'G');
        assert_eq!(snp.allele_alt, 'A');
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

        // A genotype value counts copies of the allele listed first (EIGENSTRAT
        // .snp column 5); compare copies of the same named allele, so a
        // reversed allele order cannot pass. The fixture's .snp lists each
        // VCF ALT allele first, as convertf does for an imported VCF.
        let copies = |snp: &popmix_core::genotypes::Snp, g: Genotype, allele: char| match g {
            Genotype::Missing => None,
            _ if snp.allele_ref == allele => Some(g as u8),
            _ if snp.allele_alt == allele => Some(2 - g as u8),
            _ => panic!("{} has no allele {allele}", snp.id),
        };
        for snp_idx in 0..vcf_matrix.n_snps() {
            for sample_idx in 0..vcf_matrix.n_samples() {
                let allele = eigenstrat_matrix.snps()[snp_idx].allele_ref;
                let vcf = copies(
                    &vcf_matrix.snps()[snp_idx],
                    vcf_matrix.get(snp_idx, sample_idx),
                    allele,
                );
                let eig = copies(
                    &eigenstrat_matrix.snps()[snp_idx],
                    eigenstrat_matrix.get(snp_idx, sample_idx),
                    allele,
                );
                assert_eq!(
                    vcf, eig,
                    "copies of {allele} at SNP {snp_idx}, sample {sample_idx}"
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
