//! TGENO read in place a part at a time (TASKS.md P4), so that converting
//! it holds a part of the genotypes rather than all of them (AADR v66.p1,
//! 23,089 samples x 1.23M SNPs, is 28 GB as one code per genotype):
//! - [`TgenoRows`]: rows of SNPs for the writers of SNP-major formats, a
//!   chunk of SNPs at a time, read from every selected sample's record and
//!   transposed;
//! - [`SampleGenotypes`] for [`TgenoReader`]: blocks of samples for the
//!   sample-major writers (PED, TGENO), decoded from their records.

use std::fs::File;

use popmix_core::{GenotypeRows, PopmixError, Result};
use rayon::prelude::*;

use crate::tgeno::{TgenoReader, HEADER_LEN};
use crate::writers::{SampleGenotypes, BLOCK_BYTES};

/// The genotype code of `snp` in a TGENO record.
#[inline]
fn code(record: &[u8], snp: usize) -> u8 {
    (record[snp / 4] >> (6 - 2 * (snp % 4))) & 0b11
}

/// The four codes of a record byte, first SNP in the high bits.
const BYTE_CODES: [[u8; 4]; 256] = {
    let mut t = [[0u8; 4]; 256];
    let mut b = 0;
    while b < 256 {
        t[b] = [
            (b >> 6) as u8 & 3,
            (b >> 4) as u8 & 3,
            (b >> 2) as u8 & 3,
            b as u8 & 3,
        ];
        b += 1;
    }
    t
};

/// Record bytes per tile of the transpose: 64 SNPs.
const TILE: usize = 16;

/// Read `buf.len()` bytes of `file` at `offset`.
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    #[cfg(unix)]
    {
        std::os::unix::fs::FileExt::read_exact_at(file, buf, offset)
    }
    #[cfg(windows)]
    {
        let (mut done, mut offset) = (0, offset);
        while done < buf.len() {
            let n = std::os::windows::fs::FileExt::seek_read(file, &mut buf[done..], offset)?;
            if n == 0 {
                return Err(std::io::ErrorKind::UnexpectedEof.into());
            }
            done += n;
            offset += n as u64;
        }
        Ok(())
    }
}

/// Rows of a TGENO dataset: the SNPs `keep` marks, for `samples`. A chunk
/// of SNPs is read from each sample's record with a positioned read (not
/// through the map, whose pages would stay resident) and transposed.
pub struct TgenoRows {
    file: File,
    path: std::path::PathBuf,
    record_len: usize,
    keep: Vec<bool>,
    samples: Vec<usize>,
    /// A multiple of 4, so that a chunk starts at a record byte.
    chunk_snps: usize,
    /// The decoded chunk: SNPs `chunk_start..chunk_start + chunk_len`,
    /// SNP-major, one code per selected sample.
    chunk: Vec<u8>,
    chunk_start: usize,
    chunk_len: usize,
    /// The chunk's record bytes, sample by sample.
    raw: Vec<u8>,
    next_snp: usize,
}

impl TgenoRows {
    /// Rows of `reader`: `keep` has one entry per SNP. A chunk holds at
    /// most [`BLOCK_BYTES`] of codes and their record bytes.
    pub fn new(reader: TgenoReader, keep: Vec<bool>, samples: Vec<usize>) -> Result<Self> {
        let chunk_snps = BLOCK_BYTES / 2 / samples.len().max(1);
        Self::with_chunk(reader, keep, samples, chunk_snps)
    }

    fn with_chunk(
        reader: TgenoReader,
        keep: Vec<bool>,
        samples: Vec<usize>,
        chunk_snps: usize,
    ) -> Result<Self> {
        let n_samples = reader.samples().len();
        if let Some(&bad) = samples.iter().find(|&&s| s >= n_samples) {
            return Err(PopmixError::out_of_range("sample", bad, n_samples));
        }
        if keep.len() != reader.snps().len() {
            return Err(PopmixError::Internal(format!(
                "{} SNPs selected of {}",
                keep.len(),
                reader.snps().len()
            )));
        }
        let path = reader.geno_path().to_path_buf();
        let file = File::open(&path).map_err(|e| PopmixError::io(e, &path))?;
        Ok(TgenoRows {
            file,
            path,
            record_len: reader.record_len(),
            keep,
            samples,
            chunk_snps: chunk_snps.div_ceil(4).max(1) * 4,
            chunk: Vec::new(),
            chunk_start: 0,
            chunk_len: 0,
            raw: Vec::new(),
            next_snp: 0,
        })
    }

    /// Decode the chunk that starts at `start`, a multiple of 4.
    fn decode(&mut self, start: usize) -> Result<()> {
        let len = self.chunk_snps.min(self.keep.len() - start);
        let (first, n_bytes) = (start / 4, len.div_ceil(4));
        let width = self.samples.len();
        let (file, record_len) = (&self.file, self.record_len);
        self.raw.resize(width * n_bytes, 0);
        self.raw
            .par_chunks_mut(n_bytes.max(1))
            .zip(&self.samples)
            .try_for_each(|(buf, &s)| {
                read_at(file, buf, (HEADER_LEN + s * record_len + first) as u64)
            })
            .map_err(|e| PopmixError::io(e, &self.path))?;
        // tiles of TILE record bytes (4 * TILE SNPs): each sample's bytes
        // decoded into the tile's rows
        let raw = &self.raw;
        self.chunk.resize(width * len, 0);
        self.chunk
            .par_chunks_mut((4 * TILE * width).max(1))
            .enumerate()
            .for_each(|(t, rows)| {
                let lo = t * TILE;
                let hi = (lo + TILE).min(n_bytes);
                let n_rows = rows.len() / width.max(1);
                for j in 0..width {
                    for (k, &byte) in raw[j * n_bytes + lo..j * n_bytes + hi].iter().enumerate() {
                        for (i, &c) in BYTE_CODES[byte as usize].iter().enumerate() {
                            let r = 4 * k + i;
                            if r < n_rows {
                                rows[r * width + j] = c;
                            }
                        }
                    }
                }
            });
        self.chunk_start = start;
        self.chunk_len = len;
        Ok(())
    }
}

impl GenotypeRows for TgenoRows {
    fn next_row(&mut self, out: &mut [u8]) -> Result<()> {
        while self.next_snp < self.keep.len() && !self.keep[self.next_snp] {
            self.next_snp += 1;
        }
        let snp = self.next_snp;
        if snp >= self.keep.len() {
            return Err(PopmixError::Internal(
                "read past the last SNP of a TGENO".to_string(),
            ));
        }
        if snp >= self.chunk_start + self.chunk_len || snp < self.chunk_start {
            self.decode(snp - snp % 4)?;
        }
        let width = self.samples.len();
        let at = (snp - self.chunk_start) * width;
        out[..width].copy_from_slice(&self.chunk[at..at + width]);
        self.next_snp += 1;
        Ok(())
    }
}

/// Blocks of samples, read in order by the sample-major writers: the pages
/// of the records before a block are released once it is decoded.
impl SampleGenotypes for TgenoReader {
    fn fill_samples(&self, samples: std::ops::Range<usize>, n_snps: usize, out: &mut [u8]) {
        let end = samples.end;
        out.par_chunks_mut(n_snps.max(1))
            .zip(samples)
            .for_each(|(row, s)| {
                let record = self.record(s);
                for (snp, o) in row.iter_mut().enumerate() {
                    *o = code(record, snp);
                }
            });
        self.release_before(end);
    }
}

#[cfg(test)]
mod tests;
