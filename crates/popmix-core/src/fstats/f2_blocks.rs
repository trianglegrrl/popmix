//! F2 block computation and storage.
//!
//! This module provides efficient computation and storage of f2 statistics
//! across all population pairs, organized by genomic blocks for jackknife resampling.
//!
//! The key data structure is `F2Blocks` which stores:
//! - f2 values for each block and population pair
//! - Block weights (number of valid SNPs per block per pair)
//! - Population names and indices

use ndarray::{linalg::general_mat_mul, Array2, Axis, ShapeBuilder};
use rayon::prelude::*;

use crate::fstats::blocks::{BlockCoordinate, BlockDefinition, CoordinateEvidence};
use crate::fstats::jackknife::block_jackknife;
use crate::genotype_matrix::AlleleFrequencies;
use crate::genotypes::PopulationId;

/// Precomputed F2 statistics organized by blocks.
///
/// This structure stores f2 values for all population pairs across genomic blocks,
/// enabling efficient jackknife resampling and computation of f3/f4 statistics.
#[derive(Debug, Clone)]
pub struct F2Blocks {
    /// F2 values (bias-corrected): [n_blocks, n_pairs]
    /// Each value is the sum of f2 for SNPs in that block for that pair.
    pub block_sums: Array2<f64>,

    /// Raw squared differences (not bias-corrected): [n_blocks, n_pairs]
    /// Used for FST denominator computation.
    pub block_raw_sums: Array2<f64>,

    /// Weights (SNP counts): [n_blocks, n_pairs]
    /// Number of valid SNPs contributing to each block/pair.
    pub block_weights: Array2<f64>,

    /// Heterozygosity sums: [n_blocks, n_pops]
    /// Sum of p*(1-p) for SNPs in each block for each population.
    /// Used for FST computation.
    pub block_het_sums: Array2<f64>,

    /// Heterozygosity weights: [n_blocks, n_pops]
    /// Number of valid SNPs contributing to each block/population.
    pub block_het_weights: Array2<f64>,

    /// Population names in order.
    pub populations: Vec<PopulationId>,

    /// Number of blocks.
    pub n_blocks: usize,

    /// Block lengths (number of SNPs per block).
    pub block_lengths: Vec<usize>,

    /// FST blocks on admixtools2's FST SNP set, or `None` if not available
    /// (f2 directories written before storage version 4).
    pub fst: Option<FstBlocks>,
}

/// Per-block FST numerator and denominator, on admixtools2's FST SNP set.
///
/// admixtools2 computes FST on all SNPs kept by `maxmiss`, monomorphic ones
/// included, with blocks built on those SNPs, so these blocks differ from
/// the f2 blocks. Per SNP, with allele counts n and frequencies p:
/// num = (p1 - p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1) and
/// den = num + p1(1-p1) n1/(n1-1) + p2(1-p2) n2/(n2-1). SNPs where these are
/// undefined (a population with one allele) are left out of the block means;
/// a block with none left is NaN for that pair. [`F2Blocks::fst`] leaves such
/// blocks out of that pair's jackknife, as admixtools2 `fst(f2_dir, pop1,
/// pop2)` does when it reads one pair (`read_f2(remove_na = TRUE)`).
#[derive(Debug, Clone)]
pub struct FstBlocks {
    /// SNPs per block.
    pub block_lengths: Vec<usize>,
    /// Block mean of the numerator: [n_blocks, n_pairs], pairs as `pair_index`.
    pub num: Array2<f64>,
    /// Block mean of the denominator: [n_blocks, n_pairs].
    pub den: Array2<f64>,
}

impl FstBlocks {
    /// FST blocks of `af`, in `coordinate` (None: as the kept SNPs decide),
    /// the kept SNPs observed into `evidence`.
    fn compute_in(
        af: &AlleleFrequencies,
        block_def: &BlockDefinition,
        max_missing: f64,
        coordinate: Option<BlockCoordinate>,
        evidence: &mut CoordinateEvidence,
    ) -> Self {
        let keep = kept_snps(
            af,
            &F2Options {
                max_missing,
                poly_only: false,
            },
        );
        let snps: Vec<_> = keep.iter().map(|&i| af.snps[i].clone()).collect();
        evidence.observe(&snps);
        let assignment = match coordinate {
            Some(c) => block_def.assign_blocks_in(&snps, c),
            None => block_def.assign_blocks(&snps),
        };
        let n_pops = af.populations.len();
        let n_pairs = n_pops * (n_pops + 1) / 2;

        let mut num = Array2::zeros((assignment.n_blocks, n_pairs));
        let mut den = Array2::zeros((assignment.n_blocks, n_pairs));
        num.axis_iter_mut(Axis(0))
            .into_par_iter()
            .zip(den.axis_iter_mut(Axis(0)))
            .enumerate()
            .for_each(|(b, (mut num, mut den))| {
                let (bn, bd) = fst_block(af, &keep[assignment.block_range(b)]);
                num.assign(&ndarray::ArrayView1::from(&bn));
                den.assign(&ndarray::ArrayView1::from(&bd));
            });
        FstBlocks {
            block_lengths: assignment.block_lengths,
            num,
            den,
        }
    }

    /// Per-block FST of pair `pair_idx`.
    pub fn block_values(&self, pair_idx: usize) -> Vec<f64> {
        (0..self.block_lengths.len())
            .map(|b| self.num[[b, pair_idx]] / self.den[[b, pair_idx]])
            .collect()
    }
}

/// Which SNPs go into the f2 blocks, as in admixtools2 `extract_f2`.
#[derive(Debug, Clone)]
pub struct F2Options {
    /// Drop SNPs missing in more than this fraction of the populations
    /// (admixtools2 `maxmiss`). Default: 0, so every population must have
    /// data, and all f-statistics are computed on the same SNPs.
    pub max_missing: f64,
    /// Use only SNPs whose frequency differs between populations
    /// (admixtools2 `poly_only = "f2"`). Default: true.
    pub poly_only: bool,
}

impl Default for F2Options {
    fn default() -> Self {
        F2Options {
            max_missing: 0.0,
            poly_only: true,
        }
    }
}

/// Whether a SNP is polymorphic across populations: some two observed
/// frequencies differ (missing populations are ignored).
fn is_polymorphic(freqs: ndarray::ArrayView1<f64>) -> bool {
    let mut observed = freqs.iter().filter(|f| f.is_finite());
    match observed.next() {
        Some(first) => observed.any(|f| f != first),
        None => false,
    }
}

/// Rows of `af` kept for the f2 blocks: SNPs missing in more than
/// `max_missing` of the populations are dropped and, if `poly_only`, so are
/// SNPs not polymorphic across the populations.
fn kept_snps(af: &AlleleFrequencies, opts: &F2Options) -> Vec<usize> {
    let n_pops = af.populations.len() as f64;
    (0..af.snps.len())
        .filter(|&i| {
            let missing = af.counts.row(i).iter().filter(|&&c| c == 0).count() as f64;
            missing / n_pops <= opts.max_missing
                && (!opts.poly_only || is_polymorphic(af.frequencies.row(i)))
        })
        .collect()
}

impl F2Blocks {
    /// Compute F2 blocks from allele frequencies with the admixtools2
    /// default SNP selection ([`F2Options::default`]).
    ///
    /// # Arguments
    /// * `af` - Allele frequencies for all populations
    /// * `block_def` - Block definition (size in Morgans)
    pub fn from_allele_frequencies(af: &AlleleFrequencies, block_def: &BlockDefinition) -> Self {
        Self::from_allele_frequencies_with_options(af, block_def, &F2Options::default())
    }

    /// Compute F2 blocks from allele frequencies.
    ///
    /// SNPs are selected first, then split into blocks, as in admixtools2:
    /// block boundaries depend on which SNPs are kept.
    pub fn from_allele_frequencies_with_options(
        af: &AlleleFrequencies,
        block_def: &BlockDefinition,
        opts: &F2Options,
    ) -> Self {
        let mut evidence = (CoordinateEvidence::default(), CoordinateEvidence::default());
        Self::build(af, block_def, opts, None, &mut evidence)
    }

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
            .map(|af| Self::build(&af, block_def, opts, Some(guess), &mut evidence))
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
        })
    }

    /// The blocks of `af` in `coordinate` (None: as its kept SNPs decide),
    /// the SNPs kept for f2 and for FST observed into `evidence`.
    fn build(
        af: &AlleleFrequencies,
        block_def: &BlockDefinition,
        opts: &F2Options,
        coordinate: Option<BlockCoordinate>,
        evidence: &mut (CoordinateEvidence, CoordinateEvidence),
    ) -> Self {
        let fst = Some(FstBlocks::compute_in(
            af,
            block_def,
            opts.max_missing,
            coordinate,
            &mut evidence.1,
        ));
        let keep = kept_snps(af, opts);
        let kept_snps: Vec<_> = keep.iter().map(|&i| af.snps[i].clone()).collect();
        evidence.0.observe(&kept_snps);
        let populations = af.populations.clone();
        let n_pops = populations.len();
        let n_pairs = n_pops * (n_pops + 1) / 2; // Including self-pairs for verification

        // Assign SNPs to blocks
        let block_assignment = match coordinate {
            Some(c) => block_def.assign_blocks_in(&kept_snps, c),
            None => block_def.assign_blocks(&kept_snps),
        };
        let n_blocks = block_assignment.n_blocks;

        if n_blocks == 0 {
            return F2Blocks {
                block_sums: Array2::zeros((0, n_pairs)),
                block_raw_sums: Array2::zeros((0, n_pairs)),
                block_weights: Array2::zeros((0, n_pairs)),
                block_het_sums: Array2::zeros((0, n_pops)),
                block_het_weights: Array2::zeros((0, n_pops)),
                populations,
                n_blocks: 0,
                block_lengths: Vec::new(),
                fst,
            };
        }

        // f2 and heterozygosity of each block, in parallel, each written
        // straight into its rows (no second copy of every block's sums)
        let mut block_sums = Array2::zeros((n_blocks, n_pairs));
        let mut block_raw_sums = Array2::zeros((n_blocks, n_pairs));
        let mut block_weights = Array2::zeros((n_blocks, n_pairs));
        let mut block_het_sums = Array2::zeros((n_blocks, n_pops));
        let mut block_het_weights = Array2::zeros((n_blocks, n_pops));
        block_sums
            .axis_iter_mut(Axis(0))
            .into_par_iter()
            .zip(block_raw_sums.axis_iter_mut(Axis(0)))
            .zip(block_weights.axis_iter_mut(Axis(0)))
            .zip(block_het_sums.axis_iter_mut(Axis(0)))
            .zip(block_het_weights.axis_iter_mut(Axis(0)))
            .enumerate()
            .for_each(|(b, ((((mut f2, mut raw), mut w), mut het), mut het_w))| {
                let r = f2_block(af, &keep[block_assignment.block_range(b)]);
                f2.assign(&ndarray::ArrayView1::from(&r.f2_sums));
                raw.assign(&ndarray::ArrayView1::from(&r.f2_raw_sums));
                w.assign(&ndarray::ArrayView1::from(&r.f2_weights));
                het.assign(&ndarray::ArrayView1::from(&r.het_sums));
                het_w.assign(&ndarray::ArrayView1::from(&r.het_weights));
            });

        // f2_from_geno(remove_na = TRUE): drop blocks where some pair (or a
        // population on its own) has no SNP with data, for every pair. Only
        // possible with max_missing > 0.
        let keep_blocks: Vec<usize> = (0..n_blocks)
            .filter(|&b| block_weights.row(b).iter().all(|&w| w > 0.0))
            .collect();
        let block_lengths: Vec<usize> = keep_blocks
            .iter()
            .map(|&b| block_assignment.block_lengths[b])
            .collect();

        // usually every block is kept: no copy then
        let select = |a: Array2<f64>| {
            if keep_blocks.len() == n_blocks {
                a
            } else {
                a.select(Axis(0), &keep_blocks)
            }
        };
        F2Blocks {
            block_sums: select(block_sums),
            block_raw_sums: select(block_raw_sums),
            block_weights: select(block_weights),
            block_het_sums: select(block_het_sums),
            block_het_weights: select(block_het_weights),
            populations,
            n_blocks: keep_blocks.len(),
            block_lengths,
            fst,
        }
    }

    /// Block lengths in SNPs, as jackknife weights.
    pub(crate) fn block_lengths_f64(&self) -> Vec<f64> {
        self.block_lengths.iter().map(|&l| l as f64).collect()
    }

    /// Get the number of populations.
    /// F2 blocks from each pair's per-block f2 mean and SNP count, as
    /// admixtools2 stores them (`f2_blocks` and its counts). `means[i][j]`
    /// and `counts[i][j]` hold the blocks of pair (i, j); only i <= j is
    /// read. Without the SNPs, there is nothing for FST: its sums are NaN.
    pub fn from_block_means(
        populations: Vec<PopulationId>,
        block_lengths: Vec<usize>,
        means: &[Vec<Vec<f64>>],
        counts: &[Vec<Vec<f64>>],
    ) -> crate::Result<F2Blocks> {
        let n = populations.len();
        let nb = block_lengths.len();
        let shape_ok = |a: &[Vec<Vec<f64>>]| {
            a.len() == n
                && a.iter()
                    .all(|r| r.len() == n && r.iter().all(|b| b.len() == nb))
        };
        if !shape_ok(means) || !shape_ok(counts) {
            return Err(crate::PopmixError::DimensionMismatch {
                context: "from_block_means: means and counts".to_string(),
                expected: format!("{n} x {n} x {nb}"),
                actual: "nested lists of other lengths".to_string(),
            });
        }
        let n_pairs = n * (n + 1) / 2;
        let mut sums = Array2::zeros((nb, n_pairs));
        let mut weights = Array2::zeros((nb, n_pairs));
        let mut pair = 0;
        for i in 0..n {
            for j in i..n {
                for b in 0..nb {
                    weights[[b, pair]] = counts[i][j][b];
                    sums[[b, pair]] = means[i][j][b] * counts[i][j][b];
                }
                pair += 1;
            }
        }
        Ok(F2Blocks {
            block_sums: sums,
            block_raw_sums: Array2::from_elem((nb, n_pairs), f64::NAN),
            block_weights: weights,
            block_het_sums: Array2::from_elem((nb, n), f64::NAN),
            block_het_weights: Array2::zeros((nb, n)),
            populations,
            n_blocks: nb,
            block_lengths,
            fst: None,
        })
    }

    /// These populations only, in this order (admixtools2
    /// `extract_f2_subset`); everything else, FST included, is kept.
    pub fn subset_populations(&self, pops: &[&str]) -> crate::Result<F2Blocks> {
        let idx: Vec<usize> = pops
            .iter()
            .map(|p| {
                self.pop_index(p)
                    .ok_or_else(|| crate::PopmixError::PopulationNotFound(p.to_string()))
            })
            .collect::<crate::Result<_>>()?;
        let pairs: Vec<usize> = (0..idx.len())
            .flat_map(|a| (a..idx.len()).map(move |b| (a, b)))
            .map(|(a, b)| self.pair_index(idx[a], idx[b]))
            .collect();
        Ok(F2Blocks {
            block_sums: self.block_sums.select(Axis(1), &pairs),
            block_raw_sums: self.block_raw_sums.select(Axis(1), &pairs),
            block_weights: self.block_weights.select(Axis(1), &pairs),
            block_het_sums: self.block_het_sums.select(Axis(1), &idx),
            block_het_weights: self.block_het_weights.select(Axis(1), &idx),
            populations: idx.iter().map(|&i| self.populations[i].clone()).collect(),
            n_blocks: self.n_blocks,
            block_lengths: self.block_lengths.clone(),
            fst: self.fst.as_ref().map(|f| FstBlocks {
                block_lengths: f.block_lengths.clone(),
                num: f.num.select(Axis(1), &pairs),
                den: f.den.select(Axis(1), &pairs),
            }),
        })
    }

    /// SNPs with data for both populations in each block.
    pub fn block_counts(&self, pop_a: &str, pop_b: &str) -> Option<Vec<f64>> {
        let (i, j) = (self.pop_index(pop_a)?, self.pop_index(pop_b)?);
        Some(self.block_weights.column(self.pair_index(i, j)).to_vec())
    }

    /// These blocks in this order, repeats allowed, as admixtools2
    /// `f2_blocks[, , blocks]` (its bootstrap and leave-out resamplings).
    /// FST blocks, on their own SNP set, are not carried over.
    pub fn select_blocks(&self, blocks: &[usize]) -> crate::Result<F2Blocks> {
        if blocks.is_empty() {
            return Err(crate::PopmixError::too_few("select_blocks", "block", 1, 0));
        }
        if let Some(&b) = blocks.iter().find(|&&b| b >= self.n_blocks) {
            return Err(crate::PopmixError::out_of_range("block", b, self.n_blocks));
        }
        Ok(F2Blocks {
            block_sums: self.block_sums.select(Axis(0), blocks),
            block_raw_sums: self.block_raw_sums.select(Axis(0), blocks),
            block_weights: self.block_weights.select(Axis(0), blocks),
            block_het_sums: self.block_het_sums.select(Axis(0), blocks),
            block_het_weights: self.block_het_weights.select(Axis(0), blocks),
            populations: self.populations.clone(),
            n_blocks: blocks.len(),
            block_lengths: blocks.iter().map(|&b| self.block_lengths[b]).collect(),
            fst: None,
        })
    }

    pub fn n_populations(&self) -> usize {
        self.populations.len()
    }

    /// Get the population index by name.
    pub fn pop_index(&self, name: &str) -> Option<usize> {
        self.populations.iter().position(|p| p.as_str() == name)
    }

    /// Get the pair index for two population indices.
    ///
    /// Uses triangular indexing: pair(i, j) where i <= j.
    #[inline]
    pub fn pair_index(&self, pop_a: usize, pop_b: usize) -> usize {
        let (i, j) = if pop_a <= pop_b {
            (pop_a, pop_b)
        } else {
            (pop_b, pop_a)
        };
        // Upper triangle, row by row: row i holds the n - i pairs (i, i..n)
        let n = self.populations.len();
        i * n - i * (i + 1) / 2 + j
    }

    /// The f2 of populations `i` and `j` (indices) in one block: the block's
    /// mean, NaN where the pair has no SNPs in it, and 0 for a population
    /// with itself, as admixtools2 sets the diagonal of its f2 blocks. (The
    /// stored self-pair holds the bias correction alone, not an f2.)
    #[inline]
    pub fn block_mean(&self, block: usize, i: usize, j: usize) -> f64 {
        if i == j {
            return 0.0;
        }
        let pair = self.pair_index(i, j);
        let weight = self.block_weights[[block, pair]];
        if weight > 0.0 {
            self.block_sums[[block, pair]] / weight
        } else {
            f64::NAN
        }
    }

    /// Get f2 estimate between two populations.
    ///
    /// Returns (estimate, std_err, n_snps).
    ///
    /// The estimate and standard error are the weighted block jackknife of
    /// admixtools2 `f2()`: blocks are weighted by their length in SNPs.
    ///
    /// For self-comparisons (pop_a == pop_b), returns 0.0 to match admixtools2.
    pub fn f2(&self, pop_a: &str, pop_b: &str) -> Option<(f64, f64, usize)> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;
        let pair_idx = self.pair_index(idx_a, idx_b);
        let total_snps: usize = (0..self.n_blocks)
            .map(|b| self.block_weights[[b, pair_idx]] as usize)
            .sum();

        // Self-F2 is defined as 0 (matching admixtools2)
        if idx_a == idx_b {
            return Some((0.0, 0.0, total_snps));
        }
        let (estimate, std_err) = self.f2_jackknife(idx_a, idx_b)?;
        Some((estimate, std_err, total_snps))
    }

    /// The weighted block jackknife of f2 for populations `i` != `j`. A block
    /// with no SNPs for the pair is NA in admixtools2.
    fn f2_jackknife(&self, i: usize, j: usize) -> Option<(f64, f64)> {
        let block_means: Vec<f64> = (0..self.n_blocks)
            .map(|b| self.block_mean(b, i, j))
            .collect();
        block_jackknife(&block_means, &self.block_lengths_f64())
    }

    /// Get f2 for a specific block.
    pub fn f2_block(&self, pop_a: &str, pop_b: &str, block_idx: usize) -> Option<(f64, f64)> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;
        let pair_idx = self.pair_index(idx_a, idx_b);

        let sum = self.block_sums[[block_idx, pair_idx]];
        let weight = self.block_weights[[block_idx, pair_idx]];

        Some((sum, weight))
    }

    /// Get block-level f2 values for a pair of populations.
    ///
    /// Returns f2 estimate for each block (block mean), or None if pair not found.
    /// Used for computing covariance in qpWave/qpAdm.
    /// For self-pairs (pop_a == pop_b), returns all zeros.
    pub fn get_block_f2_values(&self, pop_a: &str, pop_b: &str) -> Option<ndarray::Array1<f64>> {
        let idx_a = self.pop_index(pop_a)?;
        let idx_b = self.pop_index(pop_b)?;

        // Self-F2 is 0 for all blocks
        if idx_a == idx_b {
            return Some(ndarray::Array1::zeros(self.n_blocks));
        }

        let pair_idx = self.pair_index(idx_a, idx_b);

        let mut block_values = ndarray::Array1::zeros(self.n_blocks);

        for block_idx in 0..self.n_blocks {
            let sum = self.block_sums[[block_idx, pair_idx]];
            let weight = self.block_weights[[block_idx, pair_idx]];

            if weight > 0.0 {
                block_values[block_idx] = sum / weight;
            } else {
                // Use NaN for missing blocks
                block_values[block_idx] = f64::NAN;
            }
        }

        Some(block_values)
    }

    /// Get all pairwise f2 values as a matrix.
    ///
    /// Returns a symmetric matrix of the `f2()` estimates [n_pops, n_pops];
    /// NaN where a pair has no SNPs. Diagonal elements (self-F2) are 0.0 to
    /// match admixtools2.
    pub fn f2_matrix(&self) -> Array2<f64> {
        let n = self.populations.len();
        let mut matrix = Array2::zeros((n, n));
        for i in 0..n {
            for j in i + 1..n {
                let f2 = self
                    .f2_jackknife(i, j)
                    .map_or(f64::NAN, |(estimate, _)| estimate);
                matrix[[i, j]] = f2;
                matrix[[j, i]] = f2;
            }
        }
        matrix
    }
}

/// Result of computing F2 and heterozygosity for a block.
struct BlockResult {
    f2_sums: Vec<f64>,
    f2_raw_sums: Vec<f64>,
    f2_weights: Vec<f64>,
    het_sums: Vec<f64>,
    het_weights: Vec<f64>,
}

/// Per-SNP columns of a block, one per population, for the sums over SNPs
/// of products of two populations' columns: `aᵀb` is `pops x pops`.
fn columns(n_snps: usize, n_pops: usize) -> Array2<f64> {
    Array2::zeros((n_snps, n_pops).f())
}

/// `aᵀb`, `pops x pops`.
fn cross(a: &Array2<f64>, b: &Array2<f64>) -> Array2<f64> {
    let mut out = Array2::zeros((a.ncols(), b.ncols()));
    general_mat_mul(1.0, &a.t(), b, 0.0, &mut out);
    out
}

/// f2 sums, weights and heterozygosity of the SNPs `rows` of `af` (a block).
///
/// Per SNP and pair, f2 = (pa - pb)^2 - ca - cb with c = p(1 - p)/max(1,
/// n - 1), over the SNPs where both have data. With m the indicator of data,
/// p and q = p^2 and c zero where there is none, the block's sums over SNPs
/// are matrix products: f2 = Q + Qᵀ - C - Cᵀ - 2 PᵀP with Q = qᵀm and
/// C = cᵀm, and the weights are mᵀm. The raw sums (pa - pb)^2 are over SNPs
/// where both have at least two alleles; a population with itself is
/// computed per SNP, as these products need not cancel exactly.
fn f2_block(af: &AlleleFrequencies, rows: &[usize]) -> BlockResult {
    let n_pops = af.populations.len();
    let n_pairs = n_pops * (n_pops + 1) / 2;
    let (mut p, mut q, mut c, mut m) = (
        columns(rows.len(), n_pops),
        columns(rows.len(), n_pops),
        columns(rows.len(), n_pops),
        columns(rows.len(), n_pops),
    );
    let mut single = false;
    let mut het_sums = vec![0.0; n_pops];
    let mut het_weights = vec![0.0; n_pops];
    let mut self_f2 = vec![0.0; n_pops];
    let mut self_weight = vec![0.0; n_pops];
    for (r, &snp) in rows.iter().enumerate() {
        for a in 0..n_pops {
            let n = af.counts[[snp, a]];
            if n == 0 {
                continue;
            }
            let f = af.frequencies[[snp, a]];
            let h = f * (1.0 - f);
            let bias = h / (n as f64 - 1.0).max(1.0);
            p[[r, a]] = f;
            q[[r, a]] = f * f;
            c[[r, a]] = bias;
            m[[r, a]] = 1.0;
            self_f2[a] += 0.0 - bias - bias;
            self_weight[a] += 1.0;
            if n >= 2 {
                // admixtools2 uses p(1-p), not 2p(1-p)
                het_sums[a] += h;
                het_weights[a] += 1.0;
            } else {
                single = true;
            }
        }
    }
    let pp = cross(&p, &p);
    let qm = cross(&q, &m);
    let cm = cross(&c, &m);
    let mm = cross(&m, &m);
    // raw sums only over populations with at least two alleles
    let (raw_pp, raw_qm) = if single {
        let r = Array2::from_shape_fn(m.raw_dim().f(), |(i, a)| {
            (af.counts[[rows[i], a]] >= 2) as u8 as f64
        });
        let (pr, qr) = (&p * &r, &q * &r);
        (cross(&pr, &pr), cross(&qr, &r))
    } else {
        (pp.clone(), qm.clone())
    };

    let mut f2_sums = vec![0.0; n_pairs];
    let mut f2_raw_sums = vec![0.0; n_pairs];
    let mut f2_weights = vec![0.0; n_pairs];
    let mut pair = 0;
    for a in 0..n_pops {
        for b in a..n_pops {
            if a == b {
                f2_sums[pair] = self_f2[a];
                f2_weights[pair] = self_weight[a];
            } else {
                f2_sums[pair] =
                    qm[[a, b]] + qm[[b, a]] - cm[[a, b]] - cm[[b, a]] - 2.0 * pp[[a, b]];
                f2_raw_sums[pair] = raw_qm[[a, b]] + raw_qm[[b, a]] - 2.0 * raw_pp[[a, b]];
                f2_weights[pair] = mm[[a, b]];
            }
            pair += 1;
        }
    }
    BlockResult {
        f2_sums,
        f2_raw_sums,
        f2_weights,
        het_sums,
        het_weights,
    }
}

/// FST numerator and denominator block means of the SNPs `rows` of `af`
/// (see [`FstBlocks`]), as matrix products as in [`f2_block`]: a SNP counts
/// for a pair where both populations have at least two alleles (with one,
/// the terms are 0/0, which admixtools2 drops).
fn fst_block(af: &AlleleFrequencies, rows: &[usize]) -> (Vec<f64>, Vec<f64>) {
    let n_pops = af.populations.len();
    let n_pairs = n_pops * (n_pops + 1) / 2;
    let (mut p, mut q, mut c, mut d, mut m) = (
        columns(rows.len(), n_pops),
        columns(rows.len(), n_pops),
        columns(rows.len(), n_pops),
        columns(rows.len(), n_pops),
        columns(rows.len(), n_pops),
    );
    let mut self_num = vec![0.0; n_pops];
    let mut self_den = vec![0.0; n_pops];
    for (r, &snp) in rows.iter().enumerate() {
        for a in 0..n_pops {
            let n = af.counts[[snp, a]] as f64;
            if n < 2.0 {
                continue;
            }
            let f = af.frequencies[[snp, a]];
            let h = f * (1.0 - f);
            p[[r, a]] = f;
            q[[r, a]] = f * f;
            c[[r, a]] = h / (n - 1.0);
            d[[r, a]] = h * n / (n - 1.0);
            m[[r, a]] = 1.0;
            let num = 0.0 - h / (n - 1.0) - h / (n - 1.0);
            self_num[a] += num;
            self_den[a] += num + h * n / (n - 1.0) + h * n / (n - 1.0);
        }
    }
    let pp = cross(&p, &p);
    let qm = cross(&q, &m);
    let cm = cross(&c, &m);
    let dm = cross(&d, &m);
    let mm = cross(&m, &m);
    let mut num = vec![0.0; n_pairs];
    let mut den = vec![0.0; n_pairs];
    let mut pair = 0;
    for a in 0..n_pops {
        for b in a..n_pops {
            let (sn, sd) = if a == b {
                (self_num[a], self_den[a])
            } else {
                let sn = qm[[a, b]] + qm[[b, a]] - cm[[a, b]] - cm[[b, a]] - 2.0 * pp[[a, b]];
                (sn, sn + dm[[a, b]] + dm[[b, a]])
            };
            num[pair] = sn / mm[[a, b]];
            den[pair] = sd / mm[[a, b]];
            pair += 1;
        }
    }
    (num, den)
}

#[cfg(test)]
mod tests;
