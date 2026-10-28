//! An f2 directory written a part (a run of chromosomes) at a time, so that
//! extracting many populations holds one part's blocks rather than every
//! block of every pair (TASKS.md F4); and the f2 blocks of some of a
//! directory's populations, read without loading the others' pairs.
//!
//! The directory is the same bytes as `write_f2_blocks` of the whole.

use std::fs;
use std::path::Path;

use ndarray::Array2;
use popmix_core::fstats::{F2Arrays, FstBlocks};
use popmix_core::{F2Blocks, IoResultExt, PopmixError, PopulationId, Result};

use crate::f2_storage::{read_metadata, write_metadata, F2Metadata};
use crate::npy_stream::{read_npy_columns, NpyRowWriter};

/// The f2 directory being written.
pub struct F2DirWriter {
    dir: std::path::PathBuf,
    populations: Vec<PopulationId>,
    /// block_sums, block_raw_sums, block_weights, block_het_sums,
    /// block_het_weights, block_fst_num, block_fst_den
    arrays: Vec<NpyRowWriter>,
    block_lengths: Vec<usize>,
    fst_block_lengths: Vec<usize>,
}

const PAIR_ARRAYS: [&str; 3] = ["block_sums", "block_raw_sums", "block_weights"];
const POP_ARRAYS: [&str; 2] = ["block_het_sums", "block_het_weights"];
const FST_ARRAYS: [&str; 2] = ["block_fst_num", "block_fst_den"];

impl F2DirWriter {
    /// Start the directory `dir` (created if needed) for `populations`.
    pub fn create(dir: &Path, populations: &[PopulationId]) -> Result<Self> {
        fs::create_dir_all(dir).with_path(dir)?;
        let n = populations.len();
        let pairs = n * (n + 1) / 2;
        let mut arrays = Vec::new();
        for (names, cols) in [
            (&PAIR_ARRAYS[..], pairs),
            (&POP_ARRAYS[..], n),
            (&FST_ARRAYS[..], pairs),
        ] {
            for name in names {
                arrays.push(NpyRowWriter::create(
                    &dir.join(format!("{name}.npy")),
                    cols,
                )?);
            }
        }
        Ok(F2DirWriter {
            dir: dir.to_path_buf(),
            populations: populations.to_vec(),
            arrays,
            block_lengths: Vec::new(),
            fst_block_lengths: Vec::new(),
        })
    }

    /// Append the blocks of one part.
    pub fn push(&mut self, part: &F2Blocks) -> Result<()> {
        if part.populations() != self.populations.as_slice() {
            return Err(PopmixError::Internal(
                "a part of other populations".to_string(),
            ));
        }
        let fst = part
            .fst_blocks()
            .ok_or_else(|| PopmixError::Internal("a part without FST blocks".to_string()))?;
        let views = [
            part.block_sums().view(),
            part.block_raw_sums().view(),
            part.block_weights().view(),
            part.block_het_sums().view(),
            part.block_het_weights().view(),
            fst.num().view(),
            fst.den().view(),
        ];
        for (w, v) in self.arrays.iter_mut().zip(views) {
            w.push(v)?;
        }
        self.block_lengths.extend_from_slice(part.block_lengths());
        self.fst_block_lengths
            .extend_from_slice(fst.block_lengths());
        Ok(())
    }

    /// The SNPs of the blocks pushed so far (their lengths' sum).
    pub fn n_snps(&self) -> usize {
        self.block_lengths.iter().sum()
    }

    /// Finish the arrays and write `metadata.json`; returns the number of
    /// blocks.
    pub fn finish(
        self,
        source: Option<&str>,
        provenance: Option<serde_json::Value>,
    ) -> Result<usize> {
        for w in self.arrays {
            w.finish()?;
        }
        let n = self.populations.len();
        let metadata = F2Metadata {
            version: 4,
            populations: self
                .populations
                .iter()
                .map(|p| p.as_str().to_string())
                .collect(),
            n_blocks: self.block_lengths.len(),
            n_pairs: n * (n + 1) / 2,
            n_populations: n,
            block_lengths: self.block_lengths,
            fst_block_lengths: Some(self.fst_block_lengths),
            source: source.map(str::to_string),
            provenance,
        };
        write_metadata(&self.dir, &metadata)?;
        Ok(metadata.n_blocks)
    }
}

/// The f2 blocks of `populations` (in that order) from the f2 directory
/// `dir`, reading only their pairs' columns; as `read_f2_blocks(dir)`
/// followed by `subset_populations`, for directories of storage version 4.
pub fn read_f2_blocks_of(dir: &Path, populations: &[&str]) -> Result<F2Blocks> {
    let metadata = read_metadata(dir)?;
    if metadata.version != 4 || metadata.fst_block_lengths.is_none() {
        return crate::read_f2_blocks(dir)?.subset_populations(populations);
    }
    let all = metadata.populations.len();
    let idx: Vec<usize> = populations
        .iter()
        .map(|p| {
            metadata
                .populations
                .iter()
                .position(|q| q == p)
                .ok_or_else(|| PopmixError::PopulationNotFound((*p).to_string()))
        })
        .collect::<Result<_>>()?;
    // pair (i, j), i <= j, of all populations, in F2Blocks::pair_index order
    let pair = |i: usize, j: usize| {
        let (a, b) = if i <= j { (i, j) } else { (j, i) };
        a * all - a * (a + 1) / 2 + b
    };
    let pairs: Vec<usize> = (0..idx.len())
        .flat_map(|a| (a..idx.len()).map(move |b| (a, b)))
        .map(|(a, b)| pair(idx[a], idx[b]))
        .collect();
    let read =
        |name: &str, cols: &[usize]| read_npy_columns(&dir.join(format!("{name}.npy")), cols);
    let fst_lengths = metadata.fst_block_lengths.clone().unwrap_or_default();
    let arrays = F2Arrays {
        block_sums: read("block_sums", &pairs)?,
        block_raw_sums: read("block_raw_sums", &pairs)?,
        block_weights: read("block_weights", &pairs)?,
        block_het_sums: read("block_het_sums", &idx)?,
        block_het_weights: read("block_het_weights", &idx)?,
        populations: populations.iter().map(|p| PopulationId::from(*p)).collect(),
        block_lengths: metadata.block_lengths.clone(),
        fst: Some(FstBlocks::new(
            fst_lengths,
            read("block_fst_num", &pairs)?,
            read("block_fst_den", &pairs)?,
        )?),
    };
    let check = |a: &Array2<f64>, rows: usize, what: &str| {
        if a.nrows() == rows {
            Ok(())
        } else {
            Err(PopmixError::format(
                dir,
                format!("{what}: {} rows, the metadata says {rows}", a.nrows()),
            ))
        }
    };
    check(&arrays.block_sums, metadata.n_blocks, "block_sums")?;
    Ok(F2Blocks::from_arrays(arrays)?.with_provenance(metadata.provenance))
}
