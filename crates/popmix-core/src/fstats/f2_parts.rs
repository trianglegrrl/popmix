//! f2 blocks built a part at a time: runs of whole chromosomes, so that
//! only one part's frequency table is held (`from_parts`), and, for many
//! populations, so that only one part's blocks are held too
//! (`stream_parts`, TASKS.md F4): each part's blocks go to a sink, such as
//! an f2 directory written row by row. Blocks never span chromosomes, so the
//! parts' blocks are the whole's.

use ndarray::{Array2, Axis};

use super::blocks::{BlockCoordinate, BlockDefinition, CoordinateEvidence};
use super::f2_blocks::{F2Blocks, FstBlocks};
use super::snp_filter::F2Options;
use crate::genotype_matrix::AlleleFrequencies;
use crate::Result;

impl F2Blocks {
    /// As [`from_allele_frequencies_with_options`] over the concatenation of
    /// `parts`, the frequencies of consecutive runs of whole chromosomes in
    /// SNP order, built one part at a time so that only one part's table is
    /// held. Blocks never span chromosomes, so the parts' blocks are the
    /// whole's, in the coordinate (map or base pairs) chosen from
    /// `all_snps`, the evidence of every SNP the parts hold. None if the
    /// SNPs kept would choose another: the caller builds the whole at once.
    ///
    /// [`from_allele_frequencies_with_options`]: Self::from_allele_frequencies_with_options
    pub fn from_parts(
        parts: impl Iterator<Item = AlleleFrequencies>,
        all_snps: &CoordinateEvidence,
        block_def: &BlockDefinition,
        opts: &F2Options,
    ) -> Option<Self> {
        let guess = block_def.coordinate(all_snps);
        let mut evidence = (CoordinateEvidence::default(), CoordinateEvidence::default());
        let built: Vec<F2Blocks> = parts
            .map(|af| {
                Self::build(
                    &af,
                    block_def,
                    opts,
                    (Some(guess), Some(guess)),
                    &mut evidence,
                )
            })
            .collect();
        if block_def.coordinate(&evidence.0) != guess || block_def.coordinate(&evidence.1) != guess
        {
            return None;
        }
        let first = built.first()?;
        let rows = |f: fn(&F2Blocks) -> &Array2<f64>| {
            let views: Vec<_> = built.iter().map(|b| f(b).view()).collect();
            ndarray::concatenate(Axis(0), &views).expect("parts of the same populations")
        };
        let fsts: Vec<&FstBlocks> = built.iter().filter_map(|b| b.fst.as_ref()).collect();
        let fst_rows = |f: fn(&FstBlocks) -> &Array2<f64>| {
            let views: Vec<_> = fsts.iter().map(|b| f(b).view()).collect();
            ndarray::concatenate(Axis(0), &views).expect("parts of the same populations")
        };
        Some(F2Blocks {
            block_sums: rows(|b| &b.block_sums),
            block_raw_sums: rows(|b| &b.block_raw_sums),
            block_weights: rows(|b| &b.block_weights),
            block_het_sums: rows(|b| &b.block_het_sums),
            block_het_weights: rows(|b| &b.block_het_weights),
            populations: first.populations.clone(),
            n_blocks: built.iter().map(|b| b.n_blocks).sum(),
            block_lengths: built.iter().flat_map(|b| b.block_lengths.clone()).collect(),
            fst: Some(FstBlocks {
                block_lengths: fsts.iter().flat_map(|b| b.block_lengths.clone()).collect(),
                num: fst_rows(|b| &b.num),
                den: fst_rows(|b| &b.den),
            }),
            provenance: None,
        })
    }

    /// Build each of `parts` with f2 blocks in `coordinates.0` and FST blocks
    /// in `coordinates.1`, handing each part's blocks to `sink` in order,
    /// and return the evidence of the SNPs kept for f2 and for FST. The
    /// blocks are the whole's when `coordinates` are the coordinates that
    /// evidence chooses (`BlockDefinition::coordinate`); otherwise the
    /// caller builds again with those.
    pub fn stream_parts(
        parts: impl Iterator<Item = AlleleFrequencies>,
        coordinates: (BlockCoordinate, BlockCoordinate),
        block_def: &BlockDefinition,
        opts: &F2Options,
        sink: &mut dyn FnMut(F2Blocks) -> Result<()>,
    ) -> Result<(CoordinateEvidence, CoordinateEvidence)> {
        let mut evidence = (CoordinateEvidence::default(), CoordinateEvidence::default());
        for af in parts {
            let part = Self::build(
                &af,
                block_def,
                opts,
                (Some(coordinates.0), Some(coordinates.1)),
                &mut evidence,
            );
            sink(part)?;
        }
        Ok(evidence)
    }
}
