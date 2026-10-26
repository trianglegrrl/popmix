//! Writers for genotype file formats.
//!
//! This module provides functions to write genotype data in various formats:
//! - ANCESTRYMAP (sparse format - one line per genotype)
//! - EIGENSTRAT (dense text - one line per SNP)
//! - PED (PLINK text - one line per sample)
//! - PACKEDPED (PLINK binary .bed/.bim/.fam)
//! - PACKEDANCESTRYMAP (packed binary with EIGENSTRAT snp/ind files)

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use popmix_core::error::Result;
use popmix_core::genotype_matrix::GenotypeMatrix;
use popmix_core::genotypes::{Genotype, Sample, Snp};

// ============================================================================
// ANCESTRYMAP Format
// ============================================================================

/// Write genotype data in ANCESTRYMAP format.
///
/// Creates three files:
/// - `{prefix}.ancestrymapgeno`: Sparse genotype data (one line per non-missing genotype)
/// - `{prefix}.snp`: SNP metadata (ID CHR GPOS PPOS)
/// - `{prefix}.ind`: Sample metadata (ID SEX POP)
///
/// ANCESTRYMAP format is sparse: each line contains SNP_ID, Sample_ID, and
/// genotype (0, 1, or 2). Missing genotypes are represented by absence of entry.
pub fn write_ancestrymap(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    let prefix = prefix.as_ref();

    // Write .ancestrymapgeno file
    let geno_path = format!("{}.ancestrymapgeno", prefix.display());
    write_ancestrymap_geno(&geno_path, matrix)?;

    // Write .snp file
    let snp_path = format!("{}.snp", prefix.display());
    write_snp_file(&snp_path, matrix.snps())?;

    // Write .ind file
    let ind_path = format!("{}.ind", prefix.display());
    write_ind_file(&ind_path, matrix.samples())?;

    Ok(())
}

// ============================================================================
// EIGENSTRAT Format
// ============================================================================

/// Write genotype data in EIGENSTRAT text format.
///
/// Creates three files:
/// - `{prefix}.eigenstratgeno` or `{prefix}.geno`: Text genotypes (one line per SNP)
/// - `{prefix}.snp`: SNP metadata (ID CHR CM BP A1 A2)
/// - `{prefix}.ind`: Sample metadata (ID SEX POP)
///
/// Each line in the genotype file contains one character per sample:
/// - '0' = zero copies of reference allele (HomAlt)
/// - '1' = one copy of reference allele (Het)
/// - '2' = two copies of reference allele (HomRef)
/// - '9' = missing data
pub fn write_eigenstrat(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    let prefix = prefix.as_ref();

    // Write .geno file (using .geno extension for compatibility)
    let geno_path = format!("{}.geno", prefix.display());
    write_eigenstrat_geno(&geno_path, matrix)?;

    // Write .snp file
    let snp_path = format!("{}.snp", prefix.display());
    write_snp_file(&snp_path, matrix.snps())?;

    // Write .ind file
    let ind_path = format!("{}.ind", prefix.display());
    write_ind_file(&ind_path, matrix.samples())?;

    Ok(())
}

// ============================================================================
// PED Format (PLINK text)
// ============================================================================

/// Write genotype data in PED format (PLINK text).
///
/// Creates three files:
/// - `{prefix}.ped`: One line per sample with 7 info columns + 2 allele columns per SNP
/// - `{prefix}.pedsnp` or `{prefix}.map`: SNP metadata (CHR ID GPOS PPOS)
/// - `{prefix}.pedind`: Sample metadata (first 7 columns of .ped file)
///
/// PED format uses allele codes (A/C/G/T or 1/2) with 0 for missing.
pub fn write_ped(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    let prefix = prefix.as_ref();

    // Write .ped file
    let ped_path = format!("{}.ped", prefix.display());
    write_ped_file(&ped_path, matrix)?;

    // Write .map file (PLINK uses .map extension)
    let map_path = format!("{}.map", prefix.display());
    write_map_file(&map_path, matrix.snps())?;

    // Write .pedind file (optional but useful)
    let ind_path = format!("{}.pedind", prefix.display());
    write_pedind_file(&ind_path, matrix.samples())?;

    Ok(())
}

// ============================================================================
// PACKEDPED Format (PLINK binary)
// ============================================================================

/// Write genotype data in PACKEDPED format (PLINK binary).
///
/// Creates three files:
/// - `{prefix}.bed`: Binary packed genotypes with magic number (2 bits per genotype)
/// - `{prefix}.bim`: SNP metadata (CHR ID CM BP A1 A2)
/// - `{prefix}.fam`: Sample metadata (FID IID PID MID SEX PHENO)
///
/// The .bed format is SNP-major with 4 samples per byte:
/// - Sample 0 in bits 0-1 (LOW)
/// - Sample 1 in bits 2-3
/// - Sample 2 in bits 4-5
/// - Sample 3 in bits 6-7 (HIGH)
///
/// Encoding: 00=HomAlt, 01=Missing, 10=Het, 11=HomRef
pub fn write_packedped(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    let prefix = prefix.as_ref();

    // Write .bed file
    let bed_path = format!("{}.bed", prefix.display());
    write_plink_bed(&bed_path, matrix)?;

    // Write .bim file
    let bim_path = format!("{}.bim", prefix.display());
    write_bim_file(&bim_path, matrix.snps())?;

    // Write .fam file
    let fam_path = format!("{}.fam", prefix.display());
    write_fam_file(&fam_path, matrix.samples())?;

    Ok(())
}

// Alias for backwards compatibility
pub use write_packedped as write_plink;

// ============================================================================
// PACKEDANCESTRYMAP Format
// ============================================================================

/// Write genotype data in PACKEDANCESTRYMAP binary format.
///
/// Creates three files:
/// - `{prefix}.geno`: Binary packed genotypes with "GENO" header
/// - `{prefix}.snp`: SNP metadata (same format as EIGENSTRAT)
/// - `{prefix}.ind`: Sample metadata (same format as EIGENSTRAT)
///
/// The binary format packs 4 genotypes per byte:
/// - Sample 0 in bits 6-7 (HIGH)
/// - Sample 1 in bits 4-5
/// - Sample 2 in bits 2-3
/// - Sample 3 in bits 0-1 (LOW)
///
/// Encoding: 00=HomRef, 01=Het, 10=HomAlt, 11=Missing
/// Note: This encoding differs from PACKEDPED/PLINK!
pub fn write_packedancestrymap(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    let prefix = prefix.as_ref();

    // Write .geno file (binary)
    let geno_path = format!("{}.geno", prefix.display());
    write_packedancestrymap_geno(&geno_path, matrix)?;

    // Write .snp file
    let snp_path = format!("{}.snp", prefix.display());
    write_snp_file(&snp_path, matrix.snps())?;

    // Write .ind file
    let ind_path = format!("{}.ind", prefix.display());
    write_ind_file(&ind_path, matrix.samples())?;

    Ok(())
}

/// Write a GenotypeMatrix as TGENO (.geno, .snp, .ind): the packed format
/// stored individual by individual, as AdmixTools 8 and AADR v66 write it
/// (see [`crate::tgeno`]).
pub fn write_tgeno(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    let prefix = prefix.as_ref();
    write_tgeno_geno(&format!("{}.geno", prefix.display()), matrix)?;
    write_snp_file(&format!("{}.snp", prefix.display()), matrix.snps())?;
    write_ind_file(&format!("{}.ind", prefix.display()), matrix.samples())?;
    Ok(())
}

// ============================================================================
// Internal helper functions - Genotype files
// ============================================================================

fn write_tgeno_geno(path: &str, matrix: &GenotypeMatrix) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    let (n_samples, n_snps) = (matrix.n_samples(), matrix.n_snps());
    let ind_ids: Vec<&str> = matrix.samples().iter().map(|s| s.id.as_str()).collect();
    let snp_ids: Vec<&str> = matrix.snps().iter().map(|s| s.id.as_str()).collect();
    let header_text = format!(
        "TGENO {} {} {:x} {:x}",
        n_samples,
        n_snps,
        crate::packedancestrymap::header_hash(&ind_ids),
        crate::packedancestrymap::header_hash(&snp_ids)
    );
    let mut header = vec![0u8; crate::tgeno::HEADER_LEN];
    header[..header_text.len()].copy_from_slice(header_text.as_bytes());
    writer.write_all(&header)?;

    // one record per individual: SNP 0 in the high bits of byte 0
    let mut record = vec![0u8; n_snps.div_ceil(4)];
    for sample_idx in 0..n_samples {
        record.fill(0);
        for snp_idx in 0..n_snps {
            let bits = matrix.get(snp_idx, sample_idx).to_packed_bits();
            record[snp_idx / 4] |= bits << (6 - (snp_idx % 4) * 2);
        }
        writer.write_all(&record)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_ancestrymap_geno(path: &str, matrix: &GenotypeMatrix) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for snp_idx in 0..matrix.n_snps() {
        let snp = &matrix.snps()[snp_idx];
        for sample_idx in 0..matrix.n_samples() {
            let sample = &matrix.samples()[sample_idx];
            let geno = matrix.get(snp_idx, sample_idx);

            // Only write non-missing genotypes (sparse format)
            if !geno.is_missing() {
                let dosage = match geno {
                    Genotype::HomRef => 2,
                    Genotype::Het => 1,
                    Genotype::HomAlt => 0,
                    Genotype::Missing => continue,
                };
                writeln!(writer, "{} {} {}", snp.id, sample.id, dosage)?;
            }
        }
    }

    Ok(())
}

fn write_eigenstrat_geno(path: &str, matrix: &GenotypeMatrix) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for snp_idx in 0..matrix.n_snps() {
        for sample_idx in 0..matrix.n_samples() {
            let geno = matrix.get(snp_idx, sample_idx);
            writer.write_all(&[geno.to_eigenstrat_char()])?;
        }
        writer.write_all(b"\n")?;
    }

    Ok(())
}

fn write_ped_file(path: &str, matrix: &GenotypeMatrix) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for sample_idx in 0..matrix.n_samples() {
        let sample = &matrix.samples()[sample_idx];

        // First 7 columns: FID IID PID MID SEX PHENO 1
        // Using population as FID, sample ID as IID
        // PID/MID = 0 (unknown parents)
        // PHENO = population label (for outputgroup=YES compatibility)
        write!(
            writer,
            "{}\t{}\t0\t0\t{}\t{}\t1",
            sample.population,
            sample.id,
            sample.sex.fam_code(),
            sample.population
        )?;

        // Two allele columns per SNP
        for snp_idx in 0..matrix.n_snps() {
            let snp = &matrix.snps()[snp_idx];
            let geno = matrix.get(snp_idx, sample_idx);

            let (a1, a2) = match geno {
                Genotype::HomRef => (snp.allele_ref, snp.allele_ref),
                Genotype::Het => (snp.allele_ref, snp.allele_alt),
                Genotype::HomAlt => (snp.allele_alt, snp.allele_alt),
                Genotype::Missing => ('0', '0'),
            };
            write!(writer, "\t{}\t{}", a1, a2)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn write_packedancestrymap_geno(path: &str, matrix: &GenotypeMatrix) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let n_samples = matrix.n_samples();
    let n_snps = matrix.n_snps();

    // Bytes per row (4 samples per byte, rounded up)
    // Minimum record size is 48 bytes per EIGENSTRAT spec
    let bytes_per_row = n_samples.div_ceil(4);
    let record_size = bytes_per_row.max(48);

    // Write header: "GENO" magic + metadata, padded to record_size
    // Format: GENO {n_samples} {n_snps} {ind_hash} {snp_hash}
    let ind_ids: Vec<&str> = matrix.samples().iter().map(|s| s.id.as_str()).collect();
    let snp_ids: Vec<&str> = matrix.snps().iter().map(|s| s.id.as_str()).collect();
    let header_text = format!(
        "GENO {} {} {:x} {:x}",
        n_samples,
        n_snps,
        crate::packedancestrymap::header_hash(&ind_ids),
        crate::packedancestrymap::header_hash(&snp_ids)
    );
    let mut header = vec![0u8; record_size];
    let header_bytes = header_text.as_bytes();
    header[..header_bytes.len()].copy_from_slice(header_bytes);
    writer.write_all(&header)?;

    // Write genotype rows
    for snp_idx in 0..n_snps {
        let mut row = vec![0u8; record_size];

        for sample_idx in 0..n_samples {
            let geno = matrix.get(snp_idx, sample_idx);
            let bits = geno.to_packed_bits();

            // Sample 0 in HIGH bits (6-7), Sample 3 in LOW bits (0-1)
            let byte_idx = sample_idx / 4;
            let bit_offset = 6 - (sample_idx % 4) * 2;

            row[byte_idx] |= bits << bit_offset;
        }

        writer.write_all(&row)?;
    }

    Ok(())
}

fn write_plink_bed(path: &str, matrix: &GenotypeMatrix) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Write magic number and mode byte
    // 0x6c 0x1b = PLINK magic
    // 0x01 = SNP-major mode
    writer.write_all(&[0x6c, 0x1b, 0x01])?;

    let n_samples = matrix.n_samples();
    let bytes_per_snp = n_samples.div_ceil(4);

    for snp_idx in 0..matrix.n_snps() {
        let mut row = vec![0u8; bytes_per_snp];

        for sample_idx in 0..n_samples {
            let geno = matrix.get(snp_idx, sample_idx);
            let bits = geno.to_plink_bits();

            // PLINK: Sample 0 in LOW bits (0-1), Sample 3 in HIGH bits (6-7)
            let byte_idx = sample_idx / 4;
            let bit_offset = (sample_idx % 4) * 2;

            row[byte_idx] |= bits << bit_offset;
        }

        writer.write_all(&row)?;
    }

    Ok(())
}

// ============================================================================
// Internal helper functions - Metadata files
// ============================================================================

/// Write .snp file (EIGENSTRAT/PACKEDANCESTRYMAP format)
/// Format: ID CHR GPOS BP A1 A2
/// Note: GPOS is in Morgans (not centiMorgans) per EIGENSTRAT spec
fn write_snp_file(path: &str, snps: &[Snp]) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for snp in snps {
        // .snp file uses Morgans (position_morgans is already in Morgans)
        writeln!(
            writer,
            "{} {} {:.6} {} {} {}",
            snp.id,
            snp.chromosome,
            snp.position_morgans,
            snp.position_bp,
            snp.allele_ref,
            snp.allele_alt
        )?;
    }

    Ok(())
}

/// Write .ind file (EIGENSTRAT/PACKEDANCESTRYMAP/ANCESTRYMAP format)
/// Format: ID SEX POP
fn write_ind_file(path: &str, samples: &[Sample]) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for sample in samples {
        writeln!(
            writer,
            "{} {} {}",
            sample.id,
            sample.sex.ind_code(),
            sample.population
        )?;
    }

    Ok(())
}

/// Write .map file (PED/PLINK format)
/// Format: CHR ID CM BP A1 A2
fn write_map_file(path: &str, snps: &[Snp]) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for snp in snps {
        // .map uses Morgans, not centiMorgans
        // Include alleles so PED reader can correctly interpret genotypes
        writeln!(
            writer,
            "{}\t{}\t{:.6}\t{}\t{}\t{}",
            snp.chromosome,
            snp.id,
            snp.position_morgans,
            snp.position_bp,
            snp.allele_ref,
            snp.allele_alt
        )?;
    }

    Ok(())
}

/// Write .bim file (PACKEDPED/PLINK format)
/// Format: CHR ID CM BP A1 A2
fn write_bim_file(path: &str, snps: &[Snp]) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for snp in snps {
        // Morgans, as EIGENSOFT convertf writes it and admixtools2 reads it,
        // in full precision
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            snp.chromosome,
            snp.id,
            snp.position_morgans,
            snp.position_bp,
            snp.allele_ref,
            snp.allele_alt
        )?;
    }

    Ok(())
}

/// Write .pedind file (PED format)
/// Format: FID IID PID MID SEX PHENO 1
fn write_pedind_file(path: &str, samples: &[Sample]) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for sample in samples {
        // Same as first 7 columns of .ped file
        writeln!(
            writer,
            "{}\t{}\t0\t0\t{}\t{}\t1",
            sample.population,
            sample.id,
            sample.sex.fam_code(),
            sample.population
        )?;
    }

    Ok(())
}

/// Write .fam file (PACKEDPED/PLINK format)
/// Format: FID IID PID MID SEX PHENO
fn write_fam_file(path: &str, samples: &[Sample]) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for sample in samples {
        // Using population as FID, sample ID as IID
        // PID/MID = 0 (unknown parents)
        // PHENO = -9 (missing)
        writeln!(
            writer,
            "{}\t{}\t0\t0\t{}\t-9",
            sample.population,
            sample.id,
            sample.sex.fam_code()
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use popmix_core::genotypes::{Chromosome, Sex};
    use tempfile::TempDir;

    fn make_test_samples() -> Vec<Sample> {
        vec![
            Sample::new("IND1", "PopA"),
            Sample::new("IND2", "PopA"),
            Sample::new("IND3", "PopB"),
        ]
    }

    fn make_test_snps() -> Vec<Snp> {
        vec![
            Snp::new("rs1", Chromosome::Autosome(1), 100, 0.0001, 'A', 'G'),
            Snp::new("rs2", Chromosome::Autosome(1), 200, 0.0002, 'C', 'T'),
            Snp::new("rs3", Chromosome::Autosome(2), 300, 0.0003, 'A', 'T'),
        ]
    }

    fn make_test_matrix() -> GenotypeMatrix {
        // 3 SNPs x 3 samples
        // SNP0: 0, 1, 2
        // SNP1: 1, 2, 0
        // SNP2: 2, 0, 1
        let data = array![[0, 1, 2], [1, 2, 0], [2, 0, 1]];
        GenotypeMatrix::new(data, make_test_samples(), make_test_snps()).unwrap()
    }

    #[test]
    fn test_eigenstrat_roundtrip() {
        use crate::eigenstrat::EigenstratReader;

        let tmp = TempDir::new().unwrap();
        let prefix = tmp.path().join("test");
        let matrix = make_test_matrix();

        // Write
        write_eigenstrat(&prefix, &matrix).unwrap();

        // Read back
        let reader = EigenstratReader::new(&prefix).unwrap();
        let reloaded = reader.read_all().unwrap();

        // Compare
        assert_eq!(reloaded.n_samples(), matrix.n_samples());
        assert_eq!(reloaded.n_snps(), matrix.n_snps());

        for snp_idx in 0..matrix.n_snps() {
            for sample_idx in 0..matrix.n_samples() {
                assert_eq!(
                    reloaded.get(snp_idx, sample_idx),
                    matrix.get(snp_idx, sample_idx),
                    "Mismatch at SNP {} sample {}",
                    snp_idx,
                    sample_idx
                );
            }
        }
    }

    #[test]
    fn test_packedancestrymap_roundtrip() {
        use crate::packedancestrymap::PackedAncestryMapReader;

        let tmp = TempDir::new().unwrap();
        let prefix = tmp.path().join("test");
        let matrix = make_test_matrix();

        // Write
        write_packedancestrymap(&prefix, &matrix).unwrap();

        // Read back
        let reader = PackedAncestryMapReader::new(&prefix).unwrap();
        let reloaded = reader.read_all().unwrap();

        // Compare
        assert_eq!(reloaded.n_samples(), matrix.n_samples());
        assert_eq!(reloaded.n_snps(), matrix.n_snps());

        for snp_idx in 0..matrix.n_snps() {
            for sample_idx in 0..matrix.n_samples() {
                assert_eq!(
                    reloaded.get(snp_idx, sample_idx),
                    matrix.get(snp_idx, sample_idx),
                    "Mismatch at SNP {} sample {}",
                    snp_idx,
                    sample_idx
                );
            }
        }
    }

    #[test]
    fn test_packedped_roundtrip() {
        use crate::plink::PlinkReader;

        let tmp = TempDir::new().unwrap();
        let prefix = tmp.path().join("test");
        let matrix = make_test_matrix();

        // Write
        write_packedped(&prefix, &matrix).unwrap();

        // Read back
        let reader = PlinkReader::new(&prefix).unwrap();
        let reloaded = reader.read_all().unwrap();

        // Compare
        assert_eq!(reloaded.n_samples(), matrix.n_samples());
        assert_eq!(reloaded.n_snps(), matrix.n_snps());

        for snp_idx in 0..matrix.n_snps() {
            for sample_idx in 0..matrix.n_samples() {
                assert_eq!(
                    reloaded.get(snp_idx, sample_idx),
                    matrix.get(snp_idx, sample_idx),
                    "Mismatch at SNP {} sample {}",
                    snp_idx,
                    sample_idx
                );
            }
        }
    }

    fn sexed_matrix() -> GenotypeMatrix {
        let samples = vec![
            Sample::new("IND1", "PopA").with_sex(Sex::Male),
            Sample::new("IND2", "PopA").with_sex(Sex::Female),
            Sample::new("IND3", "PopB"),
        ];
        let data = array![[0, 1, 2], [1, 2, 0], [2, 0, 1]];
        GenotypeMatrix::new(data, samples, make_test_snps()).unwrap()
    }

    fn sexes(matrix: &GenotypeMatrix) -> Vec<Sex> {
        matrix.samples().iter().map(|s| s.sex).collect()
    }

    #[test]
    fn test_sex_survives_ind_fam_and_ped() {
        use crate::eigenstrat::EigenstratReader;
        use crate::plink::PlinkReader;

        let tmp = TempDir::new().unwrap();
        let want = vec![Sex::Male, Sex::Female, Sex::Unknown];

        let eigen = tmp.path().join("eigen");
        write_eigenstrat(&eigen, &sexed_matrix()).unwrap();
        let ind = std::fs::read_to_string(eigen.with_extension("ind")).unwrap();
        assert_eq!(ind, "IND1 M PopA\nIND2 F PopA\nIND3 U PopB\n");
        let from_ind = EigenstratReader::new(&eigen).unwrap().read_all().unwrap();
        assert_eq!(sexes(&from_ind), want);

        let plink = tmp.path().join("plink");
        write_packedped(&plink, &from_ind).unwrap();
        let fam = std::fs::read_to_string(plink.with_extension("fam")).unwrap();
        let fam_sex: Vec<&str> = fam.lines().map(|l| l.split('\t').nth(4).unwrap()).collect();
        assert_eq!(fam_sex, ["1", "2", "0"]);
        let from_fam = PlinkReader::new(&plink).unwrap().read_all().unwrap();
        assert_eq!(sexes(&from_fam), want);

        let ped = tmp.path().join("ped");
        write_ped(&ped, &from_fam).unwrap();
        let from_ped = crate::ped::PedReader::new(&ped)
            .unwrap()
            .read_all()
            .unwrap();
        assert_eq!(sexes(&from_ped), want);
    }

    #[test]
    fn test_eigenstrat_to_packedped_conversion() {
        use crate::eigenstrat::EigenstratReader;
        use crate::plink::PlinkReader;

        let tmp = TempDir::new().unwrap();
        let eigen_prefix = tmp.path().join("eigen");
        let plink_prefix = tmp.path().join("plink");
        let matrix = make_test_matrix();

        // Write EIGENSTRAT
        write_eigenstrat(&eigen_prefix, &matrix).unwrap();

        // Read EIGENSTRAT
        let reader = EigenstratReader::new(&eigen_prefix).unwrap();
        let eigen_matrix = reader.read_all().unwrap();

        // Write PACKEDPED
        write_packedped(&plink_prefix, &eigen_matrix).unwrap();

        // Read PACKEDPED
        let reader = PlinkReader::new(&plink_prefix).unwrap();
        let plink_matrix = reader.read_all().unwrap();

        // Compare genotypes
        for snp_idx in 0..matrix.n_snps() {
            for sample_idx in 0..matrix.n_samples() {
                assert_eq!(
                    plink_matrix.get(snp_idx, sample_idx),
                    matrix.get(snp_idx, sample_idx),
                    "Mismatch at SNP {} sample {}",
                    snp_idx,
                    sample_idx
                );
            }
        }
    }

    #[test]
    fn test_packedped_to_packedancestrymap_conversion() {
        use crate::packedancestrymap::PackedAncestryMapReader;
        use crate::plink::PlinkReader;

        let tmp = TempDir::new().unwrap();
        let plink_prefix = tmp.path().join("plink");
        let packed_prefix = tmp.path().join("packed");
        let matrix = make_test_matrix();

        // Write PACKEDPED
        write_packedped(&plink_prefix, &matrix).unwrap();

        // Read PACKEDPED
        let reader = PlinkReader::new(&plink_prefix).unwrap();
        let plink_matrix = reader.read_all().unwrap();

        // Write PACKEDANCESTRYMAP
        write_packedancestrymap(&packed_prefix, &plink_matrix).unwrap();

        // Read PACKEDANCESTRYMAP
        let reader = PackedAncestryMapReader::new(&packed_prefix).unwrap();
        let packed_matrix = reader.read_all().unwrap();

        // Compare genotypes
        for snp_idx in 0..matrix.n_snps() {
            for sample_idx in 0..matrix.n_samples() {
                assert_eq!(
                    packed_matrix.get(snp_idx, sample_idx),
                    matrix.get(snp_idx, sample_idx),
                    "Mismatch at SNP {} sample {}",
                    snp_idx,
                    sample_idx
                );
            }
        }
    }

    #[test]
    fn test_with_missing_genotypes() {
        use crate::eigenstrat::EigenstratReader;

        let tmp = TempDir::new().unwrap();
        let prefix = tmp.path().join("test");

        // Include missing genotypes (value 3)
        let data = array![[0, 1, 3], [3, 2, 0], [2, 3, 1]];
        let matrix = GenotypeMatrix::new(data, make_test_samples(), make_test_snps()).unwrap();

        // Write and read back
        write_eigenstrat(&prefix, &matrix).unwrap();
        let reader = EigenstratReader::new(&prefix).unwrap();
        let reloaded = reader.read_all().unwrap();

        // Verify missing genotypes are preserved
        assert_eq!(reloaded.get(0, 2), Genotype::Missing);
        assert_eq!(reloaded.get(1, 0), Genotype::Missing);
        assert_eq!(reloaded.get(2, 1), Genotype::Missing);
    }

    #[test]
    fn test_ancestrymap_write() {
        let tmp = TempDir::new().unwrap();
        let prefix = tmp.path().join("test");

        // Include missing genotypes
        let data = array![[0, 1, 3], [3, 2, 0], [2, 3, 1]];
        let matrix = GenotypeMatrix::new(data, make_test_samples(), make_test_snps()).unwrap();

        // Write ANCESTRYMAP
        write_ancestrymap(&prefix, &matrix).unwrap();

        // Read the geno file and verify sparse format
        let geno_content =
            std::fs::read_to_string(format!("{}.ancestrymapgeno", prefix.display())).unwrap();
        let lines: Vec<&str> = geno_content.lines().collect();

        // Should have 6 lines (9 genotypes - 3 missing = 6)
        assert_eq!(lines.len(), 6);

        // Verify format: SNP_ID Sample_ID Dosage
        for line in &lines {
            let parts: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(parts.len(), 3);
            // Dosage should be 0, 1, or 2
            let dosage: u8 = parts[2].parse().unwrap();
            assert!(dosage <= 2);
        }
    }

    #[test]
    fn test_ped_write() {
        let tmp = TempDir::new().unwrap();
        let prefix = tmp.path().join("test");
        let matrix = make_test_matrix();

        // Write PED
        write_ped(&prefix, &matrix).unwrap();

        // Verify .ped file exists and has correct structure
        let ped_content = std::fs::read_to_string(format!("{}.ped", prefix.display())).unwrap();
        let lines: Vec<&str> = ped_content.lines().collect();

        // Should have 3 lines (one per sample)
        assert_eq!(lines.len(), 3);

        for line in &lines {
            let parts: Vec<&str> = line.split('\t').collect();
            // 7 info columns + 2 allele columns per SNP (3 SNPs = 6 allele columns)
            assert_eq!(parts.len(), 7 + 3 * 2);
        }

        // Verify .map file
        let map_content = std::fs::read_to_string(format!("{}.map", prefix.display())).unwrap();
        let map_lines: Vec<&str> = map_content.lines().collect();
        assert_eq!(map_lines.len(), 3);
    }
}
