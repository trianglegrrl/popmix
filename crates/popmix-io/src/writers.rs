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
use popmix_core::fstats::{GenotypeRows, MatrixRows};
use popmix_core::genotype_matrix::GenotypeMatrix;
use popmix_core::genotypes::{Genotype, Sample, Snp};
use popmix_core::IoResultExt;

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
/// Each line holds SNP_ID, Sample_ID and the genotype as EIGENSTRAT codes it
/// (copies of the first allele), -1 for a missing call, as AdmixTools
/// convertf writes it.
pub fn write_ancestrymap(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    write_ancestrymap_rows(
        prefix,
        matrix.samples(),
        matrix.snps(),
        &mut all_rows(matrix),
    )
}

/// [`write_ancestrymap`] from genotype rows (one per SNP, every sample).
pub fn write_ancestrymap_rows(
    prefix: impl AsRef<Path>,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let prefix = prefix.as_ref();
    write_ancestrymap_geno(
        &format!("{}.ancestrymapgeno", prefix.display()),
        samples,
        snps,
        rows,
    )?;
    write_snp_file(&format!("{}.snp", prefix.display()), snps)?;
    write_ind_file(&format!("{}.ind", prefix.display()), samples)?;
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
    write_eigenstrat_rows(
        prefix,
        matrix.samples(),
        matrix.snps(),
        &mut all_rows(matrix),
    )
}

/// [`write_eigenstrat`] from genotype rows (one per SNP, every sample).
pub fn write_eigenstrat_rows(
    prefix: impl AsRef<Path>,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let prefix = prefix.as_ref();
    // .geno, for compatibility
    write_eigenstrat_geno(
        &format!("{}.geno", prefix.display()),
        samples.len(),
        snps.len(),
        rows,
    )?;
    write_snp_file(&format!("{}.snp", prefix.display()), snps)?;
    write_ind_file(&format!("{}.ind", prefix.display()), samples)?;
    Ok(())
}

// ============================================================================
// PED Format (PLINK text)
// ============================================================================

/// Write genotype data in PED format (PLINK text).
///
/// Creates three files:
/// - `{prefix}.ped`: One line per sample with 6 info columns + 2 allele columns per SNP
/// - `{prefix}.pedsnp` or `{prefix}.map`: SNP metadata (CHR ID GPOS PPOS)
/// - `{prefix}.pedind`: Sample metadata (first 6 columns of .ped file)
///
/// PED format uses allele codes (A/C/G/T or 1/2) with 0 for missing.
pub fn write_ped(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    write_ped_samples(prefix, matrix.samples(), matrix.snps(), matrix)
}

/// [`write_ped`] from genotypes read a block of samples at a time, so that
/// at most [`BLOCK_BYTES`] of them are held at once.
pub fn write_ped_samples(
    prefix: impl AsRef<Path>,
    samples: &[Sample],
    snps: &[Snp],
    genotypes: &dyn SampleGenotypes,
) -> Result<()> {
    write_ped_blocks(prefix, samples, snps, genotypes, block_size(snps.len()))
}

fn write_ped_blocks(
    prefix: impl AsRef<Path>,
    samples: &[Sample],
    snps: &[Snp],
    genotypes: &dyn SampleGenotypes,
    block: usize,
) -> Result<()> {
    use rayon::prelude::*;
    let prefix = prefix.as_ref();
    let path = format!("{}.ped", prefix.display());
    let mut writer = BufWriter::new(File::create(&path).with_path(&path)?);
    let n_snps = snps.len();
    // a line is about 4 bytes per SNP: build a few at a time, in parallel
    let at_once = 2 * rayon::current_num_threads();
    for_sample_blocks(samples.len(), n_snps, block, genotypes, |lo, codes| {
        let rows: Vec<&[u8]> = codes.chunks(n_snps.max(1)).collect();
        for (c, chunk) in rows.chunks(at_once).enumerate() {
            let lines: Vec<Vec<u8>> = chunk
                .par_iter()
                .enumerate()
                .map(|(k, codes)| ped_line(&samples[lo + c * at_once + k], snps, codes))
                .collect();
            for line in lines {
                writer.write_all(&line).with_path(&path)?;
            }
        }
        Ok(())
    })?;
    writer.flush().with_path(&path)?;
    // .map, as PLINK names it, and .pedind
    write_map_file(&format!("{}.map", prefix.display()), snps)?;
    write_pedind_file(&format!("{}.pedind", prefix.display()), samples)?;
    Ok(())
}

/// A sample's line of a .ped: FID (the population), IID, no parents, sex,
/// the population again as phenotype (for convertf's outputgroup), then two
/// alleles per SNP, `0 0` when missing. The code counts copies of the
/// SNP's first allele.
fn ped_line(sample: &Sample, snps: &[Snp], codes: &[u8]) -> Vec<u8> {
    let mut line = format!(
        "{}\t{}\t0\t0\t{}\t{}",
        sample.population,
        sample.id,
        sample.sex.fam_code(),
        sample.population
    )
    .into_bytes();
    line.reserve(4 * snps.len() + 1);
    let mut buf = [0u8; 4];
    for (snp, &code) in snps.iter().zip(codes) {
        let (a1, a2) = match code {
            0 => (snp.allele_alt, snp.allele_alt),
            1 => (snp.allele_ref, snp.allele_alt),
            2 => (snp.allele_ref, snp.allele_ref),
            _ => ('0', '0'),
        };
        for allele in [a1, a2] {
            line.push(b'\t');
            line.extend_from_slice(allele.encode_utf8(&mut buf).as_bytes());
        }
    }
    line.push(b'\n');
    line
}

// ============================================================================
// Formats stored sample by sample (PED, TGENO)
// ============================================================================

/// Genotypes read a block of samples at a time, for the formats stored
/// sample by sample.
pub trait SampleGenotypes: Sync {
    /// The genotype codes (0, 1, 2; 3 missing) of samples `samples` at each
    /// of `n_snps` SNPs, sample-major: `out[(s - samples.start) * n_snps + snp]`.
    fn fill_samples(&self, samples: std::ops::Range<usize>, n_snps: usize, out: &mut [u8]);
}

impl SampleGenotypes for GenotypeMatrix {
    fn fill_samples(&self, samples: std::ops::Range<usize>, n_snps: usize, out: &mut [u8]) {
        use rayon::prelude::*;
        let data = self.data();
        out.par_chunks_mut(n_snps.max(1))
            .zip(samples)
            .for_each(|(row, s)| {
                for (snp, g) in row.iter_mut().enumerate() {
                    *g = data[[snp, s]];
                }
            });
    }
}

/// Genotypes held at once by the sample-major writers: 256 MB.
pub const BLOCK_BYTES: usize = 256 << 20;

/// Samples per block for `n_snps` SNPs.
fn block_size(n_snps: usize) -> usize {
    (BLOCK_BYTES / n_snps.max(1)).max(1)
}

/// Call `write(first sample, codes)` for each block of up to `block` samples,
/// in order, with their codes sample-major.
fn for_sample_blocks(
    n_samples: usize,
    n_snps: usize,
    block: usize,
    genotypes: &dyn SampleGenotypes,
    mut write: impl FnMut(usize, &[u8]) -> Result<()>,
) -> Result<()> {
    let mut codes = vec![0u8; block.min(n_samples) * n_snps];
    let mut lo = 0;
    while lo < n_samples {
        let hi = (lo + block).min(n_samples);
        let out = &mut codes[..(hi - lo) * n_snps];
        genotypes.fill_samples(lo..hi, n_snps, out);
        write(lo, out)?;
        lo = hi;
    }
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
    write_packedped_rows(
        prefix,
        matrix.samples(),
        matrix.snps(),
        &mut all_rows(matrix),
    )
}

/// [`write_packedped`] from genotype rows (one per SNP, every sample).
pub fn write_packedped_rows(
    prefix: impl AsRef<Path>,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let prefix = prefix.as_ref();
    write_plink_bed(
        &format!("{}.bed", prefix.display()),
        samples.len(),
        snps.len(),
        rows,
    )?;
    write_bim_file(&format!("{}.bim", prefix.display()), snps)?;
    write_fam_file(&format!("{}.fam", prefix.display()), samples)?;
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
    write_packedancestrymap_rows(
        prefix,
        matrix.samples(),
        matrix.snps(),
        &mut all_rows(matrix),
    )
}

/// [`write_packedancestrymap`] from genotype rows (one per SNP, every sample).
pub fn write_packedancestrymap_rows(
    prefix: impl AsRef<Path>,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let prefix = prefix.as_ref();
    write_packedancestrymap_geno(&format!("{}.geno", prefix.display()), samples, snps, rows)?;
    write_snp_file(&format!("{}.snp", prefix.display()), snps)?;
    write_ind_file(&format!("{}.ind", prefix.display()), samples)?;
    Ok(())
}

/// Every SNP of every sample of `matrix`, as rows.
fn all_rows(matrix: &GenotypeMatrix) -> MatrixRows<'_> {
    MatrixRows::new(
        matrix,
        (0..matrix.n_snps()).collect(),
        (0..matrix.n_samples()).collect(),
    )
}

/// A genotype code of a row (0, 1, 2; 3 missing) as a genotype.
fn genotype(code: u8) -> Genotype {
    match code {
        0 => Genotype::HomRef,
        1 => Genotype::Het,
        2 => Genotype::HomAlt,
        _ => Genotype::Missing,
    }
}

/// Write a GenotypeMatrix as TGENO (.geno, .snp, .ind): the packed format
/// stored individual by individual, as AdmixTools 8 and AADR v66 write it
/// (see [`crate::tgeno`]).
pub fn write_tgeno(prefix: impl AsRef<Path>, matrix: &GenotypeMatrix) -> Result<()> {
    write_tgeno_samples(prefix, matrix.samples(), matrix.snps(), matrix)
}

/// [`write_tgeno`] from genotypes read a block of samples at a time.
pub fn write_tgeno_samples(
    prefix: impl AsRef<Path>,
    samples: &[Sample],
    snps: &[Snp],
    genotypes: &dyn SampleGenotypes,
) -> Result<()> {
    write_tgeno_blocks(prefix, samples, snps, genotypes, block_size(snps.len()))
}

fn write_tgeno_blocks(
    prefix: impl AsRef<Path>,
    samples: &[Sample],
    snps: &[Snp],
    genotypes: &dyn SampleGenotypes,
    block: usize,
) -> Result<()> {
    use rayon::prelude::*;
    let prefix = prefix.as_ref();
    let path = format!("{}.geno", prefix.display());
    let mut writer = BufWriter::new(File::create(&path).with_path(&path)?);

    let (n_samples, n_snps) = (samples.len(), snps.len());
    let ind_ids: Vec<&str> = samples.iter().map(|s| s.id.as_str()).collect();
    let snp_ids: Vec<&str> = snps.iter().map(|s| s.id.as_str()).collect();
    let header_text = format!(
        "TGENO {} {} {:x} {:x}",
        n_samples,
        n_snps,
        crate::packedancestrymap::header_hash(&ind_ids),
        crate::packedancestrymap::header_hash(&snp_ids)
    );
    let mut header = vec![0u8; crate::tgeno::HEADER_LEN];
    header[..header_text.len()].copy_from_slice(header_text.as_bytes());
    writer.write_all(&header).with_path(&path)?;

    // one record per individual: SNP 0 in the high bits of byte 0
    let record_len = n_snps.div_ceil(4);
    for_sample_blocks(n_samples, n_snps, block, genotypes, |_, codes| {
        let n = codes.len() / n_snps.max(1);
        let mut records = vec![0u8; n * record_len];
        records
            .par_chunks_mut(record_len.max(1))
            .zip(codes.par_chunks(n_snps.max(1)))
            .for_each(|(record, codes)| {
                for (snp, &code) in codes.iter().enumerate() {
                    record[snp / 4] |= code << (6 - (snp % 4) * 2);
                }
            });
        writer.write_all(&records).with_path(&path)?;
        Ok(())
    })?;
    writer.flush().with_path(&path)?;
    write_snp_file(&format!("{}.snp", prefix.display()), snps)?;
    write_ind_file(&format!("{}.ind", prefix.display()), samples)?;
    Ok(())
}

fn write_ancestrymap_geno(
    path: &str,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let file = File::create(path).with_path(path)?;
    let mut writer = BufWriter::new(file);
    let mut row = vec![0u8; samples.len()];
    for snp in snps {
        rows.next_row(&mut row)?;
        for (sample, &code) in samples.iter().zip(&row) {
            // the EIGENSTRAT code, every genotype, -1 missing, as convertf
            let dosage: i8 = match genotype(code) {
                Genotype::HomRef => 0,
                Genotype::Het => 1,
                Genotype::HomAlt => 2,
                Genotype::Missing => -1,
            };
            writeln!(writer, "{} {} {}", snp.id, sample.id, dosage).with_path(path)?;
        }
    }
    writer.flush().with_path(path)?;
    Ok(())
}

// ============================================================================
// Internal helper functions - Metadata files
// ============================================================================

/// Write .snp file (EIGENSTRAT/PACKEDANCESTRYMAP format)
/// Format: ID CHR GPOS BP A1 A2
/// Note: GPOS is in Morgans (not centiMorgans) per EIGENSTRAT spec
fn write_snp_file(path: &str, snps: &[Snp]) -> Result<()> {
    let file = File::create(path).with_path(path)?;
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
        )
        .with_path(path)?;
    }

    Ok(())
}

/// Write .ind file (EIGENSTRAT/PACKEDANCESTRYMAP/ANCESTRYMAP format)
/// Format: ID SEX POP
fn write_ind_file(path: &str, samples: &[Sample]) -> Result<()> {
    let file = File::create(path).with_path(path)?;
    let mut writer = BufWriter::new(file);

    for sample in samples {
        writeln!(
            writer,
            "{} {} {}",
            sample.id,
            sample.sex.ind_code(),
            sample.population
        )
        .with_path(path)?;
    }

    Ok(())
}

/// Write .map file (PED/PLINK format)
/// Format: CHR ID CM BP A1 A2
fn write_map_file(path: &str, snps: &[Snp]) -> Result<()> {
    let file = File::create(path).with_path(path)?;
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
        )
        .with_path(path)?;
    }

    Ok(())
}

/// Write .bim file (PACKEDPED/PLINK format)
/// Format: CHR ID CM BP A1 A2
fn write_bim_file(path: &str, snps: &[Snp]) -> Result<()> {
    let file = File::create(path).with_path(path)?;
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
        )
        .with_path(path)?;
    }

    Ok(())
}

/// Write .pedind file (PED format)
/// Format: FID IID PID MID SEX PHENO
fn write_pedind_file(path: &str, samples: &[Sample]) -> Result<()> {
    let file = File::create(path).with_path(path)?;
    let mut writer = BufWriter::new(file);

    for sample in samples {
        // Same as the first 6 columns of the .ped file
        writeln!(
            writer,
            "{}\t{}\t0\t0\t{}\t{}",
            sample.population,
            sample.id,
            sample.sex.fam_code(),
            sample.population
        )
        .with_path(path)?;
    }

    Ok(())
}

/// Write .fam file (PACKEDPED/PLINK format)
/// Format: FID IID PID MID SEX PHENO
fn write_fam_file(path: &str, samples: &[Sample]) -> Result<()> {
    let file = File::create(path).with_path(path)?;
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
        )
        .with_path(path)?;
    }

    Ok(())
}

mod rows;
use rows::{write_eigenstrat_geno, write_packedancestrymap_geno, write_plink_bed};

#[cfg(test)]
mod tests;
