//! The SNP-major genotype files (EIGENSTRAT text .geno, PACKEDANCESTRYMAP
//! .geno, PLINK .bed) written from rows a batch at a time (TASKS.md P4):
//! the source decodes a batch (in parallel where it can), the records are
//! encoded in parallel and written in order, so a large dataset is
//! converted on every thread and with at most [`BLOCK_BYTES`] held.

use std::fs::File;
use std::io::{BufWriter, Write};

use popmix_core::error::Result;
use popmix_core::fstats::GenotypeRows;
use popmix_core::genotypes::Sample;
use popmix_core::genotypes::Snp;
use popmix_core::IoResultExt;
use rayon::prelude::*;

use super::{genotype, BLOCK_BYTES};

/// Write `n_snps` records of `record_len` bytes to `writer`, each encoded
/// by `encode(codes, record)` from a row of `width` codes.
fn write_records(
    writer: &mut impl Write,
    path: &str,
    n_snps: usize,
    width: usize,
    record_len: usize,
    rows: &mut dyn GenotypeRows,
    encode: impl Fn(&[u8], &mut [u8]) + Sync,
) -> Result<()> {
    if width == 0 {
        let mut record = vec![0u8; record_len];
        for _ in 0..n_snps {
            rows.next_row(&mut [])?;
            encode(&[], &mut record);
            writer.write_all(&record).with_path(path)?;
        }
        return Ok(());
    }
    let batch = (BLOCK_BYTES / 2 / width.max(record_len)).clamp(1, n_snps.max(1));
    let mut codes = vec![0u8; batch * width];
    let mut out = vec![0u8; batch * record_len];
    let mut done = 0;
    while done < n_snps {
        let n = batch.min(n_snps - done);
        let (codes, out) = (&mut codes[..n * width], &mut out[..n * record_len]);
        rows.next_rows(codes, width)?;
        out.par_chunks_mut(record_len)
            .zip(codes.par_chunks(width))
            .for_each(|(record, row)| encode(row, record));
        writer.write_all(out).with_path(path)?;
        done += n;
    }
    Ok(())
}

pub(super) fn write_eigenstrat_geno(
    path: &str,
    n_samples: usize,
    n_snps: usize,
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path).with_path(path)?);
    write_records(
        &mut writer,
        path,
        n_snps,
        n_samples,
        n_samples + 1,
        rows,
        |row, line| {
            for (c, &code) in line.iter_mut().zip(row) {
                *c = genotype(code).to_eigenstrat_char();
            }
            line[n_samples] = b'\n';
        },
    )?;
    writer.flush().with_path(path)?;
    Ok(())
}

pub(super) fn write_packedancestrymap_geno(
    path: &str,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path).with_path(path)?);
    let n_samples = samples.len();
    let n_snps = snps.len();

    // Bytes per row (4 samples per byte, rounded up)
    // Minimum record size is 48 bytes per EIGENSTRAT spec
    let record_size = n_samples.div_ceil(4).max(48);

    // Write header: "GENO" magic + metadata, padded to record_size
    // Format: GENO {n_samples} {n_snps} {ind_hash} {snp_hash}
    let ind_ids: Vec<&str> = samples.iter().map(|s| s.id.as_str()).collect();
    let snp_ids: Vec<&str> = snps.iter().map(|s| s.id.as_str()).collect();
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
    writer.write_all(&header).with_path(path)?;

    write_records(
        &mut writer,
        path,
        n_snps,
        n_samples,
        record_size,
        rows,
        |row, record| {
            record.fill(0);
            for (sample_idx, &code) in row.iter().enumerate() {
                // Sample 0 in HIGH bits (6-7), Sample 3 in LOW bits (0-1)
                let bit_offset = 6 - (sample_idx % 4) * 2;
                record[sample_idx / 4] |= genotype(code).to_packed_bits() << bit_offset;
            }
        },
    )?;
    writer.flush().with_path(path)?;
    Ok(())
}

pub(super) fn write_plink_bed(
    path: &str,
    n_samples: usize,
    n_snps: usize,
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path).with_path(path)?);

    // Write magic number and mode byte
    // 0x6c 0x1b = PLINK magic
    // 0x01 = SNP-major mode
    writer.write_all(&[0x6c, 0x1b, 0x01]).with_path(path)?;

    write_records(
        &mut writer,
        path,
        n_snps,
        n_samples,
        n_samples.div_ceil(4),
        rows,
        |row, record| {
            record.fill(0);
            for (sample_idx, &code) in row.iter().enumerate() {
                // PLINK: Sample 0 in LOW bits (0-1), Sample 3 in HIGH bits (6-7)
                record[sample_idx / 4] |= genotype(code).to_plink_bits() << ((sample_idx % 4) * 2);
            }
        },
    )?;
    writer.flush().with_path(path)?;
    Ok(())
}
