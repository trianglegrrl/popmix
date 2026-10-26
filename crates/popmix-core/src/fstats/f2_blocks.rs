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

use crate::fstats::blocks::BlockDefinition;
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
    fn compute(af: &AlleleFrequencies, block_def: &BlockDefinition, max_missing: f64) -> Self {
        let keep = kept_snps(
            af,
            &F2Options {
                max_missing,
                poly_only: false,
            },
        );
        let snps: Vec<_> = keep.iter().map(|&i| af.snps[i].clone()).collect();
        let assignment = block_def.assign_blocks(&snps);
        let n_pops = af.populations.len();
        let n_pairs = n_pops * (n_pops + 1) / 2;

        let blocks: Vec<(Vec<f64>, Vec<f64>)> = (0..assignment.n_blocks)
            .into_par_iter()
            .map(|b| fst_block(af, &keep[assignment.block_range(b)]))
            .collect();

        let mut num = Array2::zeros((blocks.len(), n_pairs));
        let mut den = Array2::zeros((blocks.len(), n_pairs));
        for (b, (bn, bd)) in blocks.iter().enumerate() {
            for pair in 0..n_pairs {
                num[[b, pair]] = bn[pair];
                den[[b, pair]] = bd[pair];
            }
        }
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
        let fst = Some(FstBlocks::compute(af, block_def, opts.max_missing));
        let keep = kept_snps(af, opts);
        let kept_snps: Vec<_> = keep.iter().map(|&i| af.snps[i].clone()).collect();
        let populations = af.populations.clone();
        let n_pops = populations.len();
        let n_pairs = n_pops * (n_pops + 1) / 2; // Including self-pairs for verification

        // Assign SNPs to blocks
        let block_assignment = block_def.assign_blocks(&kept_snps);
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

        // Compute f2 and heterozygosity for each block in parallel
        let block_results: Vec<BlockResult> = (0..n_blocks)
            .into_par_iter()
            .map(|block_idx| f2_block(af, &keep[block_assignment.block_range(block_idx)]))
            .collect();

        // Assemble into arrays
        let mut block_sums = Array2::zeros((n_blocks, n_pairs));
        let mut block_raw_sums = Array2::zeros((n_blocks, n_pairs));
        let mut block_weights = Array2::zeros((n_blocks, n_pairs));
        let mut block_het_sums = Array2::zeros((n_blocks, n_pops));
        let mut block_het_weights = Array2::zeros((n_blocks, n_pops));

        for (block_idx, result) in block_results.into_iter().enumerate() {
            for pair_idx in 0..n_pairs {
                block_sums[[block_idx, pair_idx]] = result.f2_sums[pair_idx];
                block_raw_sums[[block_idx, pair_idx]] = result.f2_raw_sums[pair_idx];
                block_weights[[block_idx, pair_idx]] = result.f2_weights[pair_idx];
            }
            for (pop_idx, (het_sum, het_weight)) in result
                .het_sums
                .into_iter()
                .zip(result.het_weights)
                .enumerate()
            {
                block_het_sums[[block_idx, pop_idx]] = het_sum;
                block_het_weights[[block_idx, pop_idx]] = het_weight;
            }
        }

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

        F2Blocks {
            block_sums: block_sums.select(Axis(0), &keep_blocks),
            block_raw_sums: block_raw_sums.select(Axis(0), &keep_blocks),
            block_weights: block_weights.select(Axis(0), &keep_blocks),
            block_het_sums: block_het_sums.select(Axis(0), &keep_blocks),
            block_het_weights: block_het_weights.select(Axis(0), &keep_blocks),
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
            return Err(crate::PopmixError::Validation(format!(
                "from_block_means: means and counts must be {n} x {n} x {nb}"
            )));
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
                self.pop_index(p).ok_or_else(|| {
                    crate::PopmixError::Validation(format!(
                        "population {p:?} is not in the f2 blocks"
                    ))
                })
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
            return Err(crate::PopmixError::Validation(
                "select_blocks: no blocks selected".to_string(),
            ));
        }
        if let Some(&b) = blocks.iter().find(|&&b| b >= self.n_blocks) {
            return Err(crate::PopmixError::Validation(format!(
                "select_blocks: block {b} out of range ({} blocks)",
                self.n_blocks
            )));
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
mod tests {
    use super::*;
    use crate::fstats::f2::f2_single_snp;

    /// Compute f2 sums, weights, and heterozygosity for a block of SNPs.
    fn direct_f2_block(af: &AlleleFrequencies, rows: &[usize]) -> BlockResult {
        let n_pops = af.populations.len();
        let n_pairs = n_pops * (n_pops + 1) / 2;

        let mut f2_sums = vec![0.0; n_pairs];
        let mut f2_raw_sums = vec![0.0; n_pairs];
        let mut f2_weights = vec![0.0; n_pairs];
        let mut het_sums = vec![0.0; n_pops];
        let mut het_weights = vec![0.0; n_pops];

        for &snp_idx in rows {
            // Compute heterozygosity for each population: h = p * (1 - p)
            // Note: admixtools2 uses p*(1-p), not 2*p*(1-p)
            for pop_idx in 0..n_pops {
                let count = af.get_count(snp_idx, pop_idx);

                if let Some(freq) = af.get(snp_idx, pop_idx) {
                    if count >= 2 {
                        let het = freq * (1.0 - freq);
                        het_sums[pop_idx] += het;
                        het_weights[pop_idx] += 1.0;
                    }
                }
            }

            // For each population pair
            let mut pair_idx = 0;
            for i in 0..n_pops {
                for j in i..n_pops {
                    let freq_a = af.get(snp_idx, i);
                    let freq_b = af.get(snp_idx, j);
                    let count_a = af.get_count(snp_idx, i);
                    let count_b = af.get_count(snp_idx, j);

                    let (f2, weight) = f2_single_snp(freq_a, freq_b, count_a, count_b);

                    // Also compute raw squared difference (without bias correction)
                    // for FST denominator
                    if let (Some(pa), Some(pb)) = (freq_a, freq_b) {
                        if count_a >= 2 && count_b >= 2 {
                            let raw = (pa - pb).powi(2);
                            f2_raw_sums[pair_idx] += raw * weight;
                        }
                    }

                    f2_sums[pair_idx] += f2 * weight;
                    f2_weights[pair_idx] += weight;

                    pair_idx += 1;
                }
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

    use crate::genotype_matrix::GenotypeMatrix;
    use crate::genotypes::{Chromosome, Sample, Snp};
    use approx::assert_relative_eq;

    fn make_test_data() -> AlleleFrequencies {
        // Create test data with 3 populations, 20 SNPs across 2 cM
        let mut samples = Vec::new();
        for i in 0..10 {
            samples.push(Sample::new(format!("A{}", i), PopulationId::new("PopA")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("B{}", i), PopulationId::new("PopB")));
        }
        for i in 0..10 {
            samples.push(Sample::new(format!("C{}", i), PopulationId::new("PopC")));
        }

        let snps: Vec<Snp> = (0..20)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    i as f64 * 0.1, // 0.1 cM spacing
                    'A',
                    'G',
                )
            })
            .collect();

        // Create genotype data
        // PopA: mostly ref (0), PopB: mostly alt (2), PopC: mixed
        let mut data = ndarray::Array2::zeros((20, 30));
        for snp_idx in 0..20 {
            // PopA (columns 0-9): mostly 0
            for j in 0..10 {
                data[[snp_idx, j]] = if snp_idx % 5 == 0 { 1 } else { 0 };
            }
            // PopB (columns 10-19): mostly 2
            for j in 10..20 {
                data[[snp_idx, j]] = if snp_idx % 5 == 0 { 1 } else { 2 };
            }
            // PopC (columns 20-29): mixed
            for j in 20..30 {
                data[[snp_idx, j]] = ((snp_idx + j) % 3) as u8;
            }
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        AlleleFrequencies::from_genotypes(&geno)
    }

    /// Allele frequencies for PopA, PopB, PopC at SNPs on chromosome 1.
    /// `rows` holds (position in Morgans, [(freq, count); 3]); count 0 is missing.
    fn af_from(rows: &[(f64, [(f64, u32); 3])]) -> AlleleFrequencies {
        let n = rows.len();
        let mut frequencies = Array2::from_elem((n, 3), f64::NAN);
        let mut counts = Array2::zeros((n, 3));
        for (i, (_, pops)) in rows.iter().enumerate() {
            for (j, &(f, c)) in pops.iter().enumerate() {
                if c > 0 {
                    frequencies[[i, j]] = f;
                    counts[[i, j]] = c;
                }
            }
        }
        let snps = rows
            .iter()
            .enumerate()
            .map(|(i, (pos, _))| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1),
                    i as u32,
                    *pos,
                    'A',
                    'G',
                )
            })
            .collect();
        AlleleFrequencies {
            frequencies,
            counts,
            populations: ["PopA", "PopB", "PopC"].map(PopulationId::new).to_vec(),
            snps,
        }
    }

    /// Per-SNP bias-corrected f2: (p1 - p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1).
    fn snp_f2(a: (f64, u32), b: (f64, u32)) -> f64 {
        let corr = |(p, n): (f64, u32)| p * (1.0 - p) / (n as f64 - 1.0).max(1.0);
        (a.0 - b.0).powi(2) - corr(a) - corr(b)
    }

    /// Per-SNP FST numerator and denominator block means (see [`FstBlocks`]).
    fn direct_fst_block(af: &AlleleFrequencies, rows: &[usize]) -> (Vec<f64>, Vec<f64>) {
        let n_pops = af.populations.len();
        let n_pairs = n_pops * (n_pops + 1) / 2;
        let (mut num, mut den, mut cnt) =
            (vec![0.0; n_pairs], vec![0.0; n_pairs], vec![0.0; n_pairs]);
        for &snp in rows {
            let (p, n) = (af.frequencies.row(snp), af.counts.row(snp));
            let mut pair = 0;
            for i in 0..n_pops {
                for j in i..n_pops {
                    let (p1, p2) = (p[i], p[j]);
                    let (n1, n2) = (n[i] as f64, n[j] as f64);
                    let (h1, h2) = (p1 * (1.0 - p1), p2 * (1.0 - p2));
                    let nm = (p1 - p2).powi(2) - h1 / (n1 - 1.0) - h2 / (n2 - 1.0);
                    let dn = nm + h1 * n1 / (n1 - 1.0) + h2 * n2 / (n2 - 1.0);
                    if !nm.is_nan() && !dn.is_nan() {
                        num[pair] += nm;
                        den[pair] += dn;
                        cnt[pair] += 1.0;
                    }
                    pair += 1;
                }
            }
        }
        let mean = |v: Vec<f64>| v.iter().zip(&cnt).map(|(s, c)| s / c).collect();
        (mean(num), mean(den))
    }

    /// Random frequencies of `n_pops` populations, with missing data and
    /// populations of one allele (pseudo-haploid samples).
    fn random_af(n_snps: usize, n_pops: usize) -> AlleleFrequencies {
        let mut state = 2024u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut frequencies = Array2::from_elem((n_snps, n_pops), f64::NAN);
        let mut counts = Array2::zeros((n_snps, n_pops));
        for i in 0..n_snps {
            for j in 0..n_pops {
                let n = [0u32, 1, 2, 5, 10, 20, 20, 20][(next() % 8) as usize];
                if n > 0 {
                    counts[[i, j]] = n;
                    frequencies[[i, j]] = (next() % (n as u64 + 1)) as f64 / n as f64;
                }
            }
        }
        let snps = (0..n_snps)
            .map(|i| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1),
                    i as u32,
                    i as f64 * 1e-4,
                    'A',
                    'G',
                )
            })
            .collect();
        AlleleFrequencies {
            frequencies,
            counts,
            populations: (0..n_pops)
                .map(|j| PopulationId::new(format!("P{j}")))
                .collect(),
            snps,
        }
    }

    fn assert_close(got: &[f64], want: &[f64], what: &str) {
        assert_eq!(got.len(), want.len());
        for (k, (g, w)) in got.iter().zip(want).enumerate() {
            let ok = (g.is_nan() && w.is_nan()) || (g - w).abs() <= 1e-15 + 1e-12 * w.abs();
            assert!(ok, "{what} [{k}]: {g} != {w}");
        }
    }

    #[test]
    fn block_products_match_the_per_snp_sums() {
        let af = random_af(600, 9);
        let rows: Vec<usize> = (0..600).filter(|i| i % 7 != 3).collect();
        let (got, want) = (f2_block(&af, &rows), direct_f2_block(&af, &rows));
        assert_close(&got.f2_sums, &want.f2_sums, "f2");
        assert_close(&got.f2_raw_sums, &want.f2_raw_sums, "raw");
        assert_eq!(got.f2_weights, want.f2_weights);
        assert_close(&got.het_sums, &want.het_sums, "het");
        assert_eq!(got.het_weights, want.het_weights);
        let (num, den) = fst_block(&af, &rows);
        let (want_num, want_den) = direct_fst_block(&af, &rows);
        assert_close(&num, &want_num, "fst num");
        assert_close(&den, &want_den, "fst den");
    }

    #[test]
    fn blocks_do_not_depend_on_the_thread_count() {
        let af = random_af(4_000, 6);
        let opts = F2Options {
            max_missing: 1.0,
            poly_only: true,
        };
        let run = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    F2Blocks::from_allele_frequencies_with_options(
                        &af,
                        &BlockDefinition::new(0.05),
                        &opts,
                    )
                })
        };
        let one = run(1);
        for threads in [3, 8] {
            let many = run(threads);
            assert_eq!(one.block_sums, many.block_sums);
            assert_eq!(one.block_raw_sums, many.block_raw_sums);
            assert_eq!(
                one.fst.as_ref().unwrap().num,
                many.fst.as_ref().unwrap().num
            );
        }
    }

    #[test]
    fn test_snp_missing_in_any_population_is_dropped_for_every_pair() {
        // admixtools2 maxmiss = 0: SNP 1 is missing in PopC, so f2(A, B) must
        // not use it either.
        let rows = [
            (0.00, [(0.1, 10), (0.5, 10), (0.3, 10)]),
            (0.01, [(0.9, 10), (0.2, 10), (0.0, 0)]),
            (0.10, [(0.4, 10), (0.6, 10), (0.5, 10)]),
        ];
        let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
        assert_eq!(f2.block_lengths, vec![1, 1]);
        let (est, _, n_snps) = f2.f2("PopA", "PopB").unwrap();
        let expected =
            (snp_f2(rows[0].1[0], rows[0].1[1]) + snp_f2(rows[2].1[0], rows[2].1[1])) / 2.0;
        assert_relative_eq!(est, expected, max_relative = 1e-12);
        assert_eq!(n_snps, 2);
    }

    #[test]
    fn test_monomorphic_snps_are_not_used() {
        // admixtools2 poly_only = "f2": SNP 1 has the same frequency in every
        // population, so it is left out of the f2 blocks and block lengths.
        let rows = [
            (0.00, [(0.1, 10), (0.5, 10), (0.3, 10)]),
            (0.01, [(0.0, 10), (0.0, 10), (0.0, 10)]),
            (0.02, [(0.4, 10), (0.6, 10), (0.5, 10)]),
        ];
        let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
        assert_eq!(f2.block_lengths, vec![2]);
    }

    #[test]
    fn test_fst_blocks_keep_monomorphic_snps() {
        // admixtools2 computes FST on every SNP kept by maxmiss: the
        // monomorphic SNP 1 is in the FST block but not the f2 block, and the
        // SNP missing in PopB is in neither.
        let rows = [
            (0.00, [(0.1, 10), (0.5, 10), (0.3, 10)]),
            (0.01, [(0.2, 10), (0.2, 10), (0.2, 10)]),
            (0.02, [(0.4, 10), (0.0, 0), (0.5, 10)]),
            (0.03, [(0.4, 10), (0.6, 10), (0.5, 10)]),
        ];
        let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
        assert_eq!(f2.block_lengths, vec![2]);
        let fst = f2.fst.as_ref().unwrap();
        assert_eq!(fst.block_lengths, vec![3]);

        // Block FST(A, B) = mean(num) / mean(den) over SNPs 0, 1, 3
        let (num, den): (Vec<f64>, Vec<f64>) = [0, 1, 3]
            .iter()
            .map(|&i| {
                let ((p1, n1), (p2, n2)) = (rows[i].1[0], rows[i].1[1]);
                let (n1, n2) = (n1 as f64, n2 as f64);
                let (h1, h2) = (p1 * (1.0 - p1), p2 * (1.0 - p2));
                let num = (p1 - p2).powi(2) - h1 / (n1 - 1.0) - h2 / (n2 - 1.0);
                (num, num + h1 * n1 / (n1 - 1.0) + h2 * n2 / (n2 - 1.0))
            })
            .unzip();
        let expected = num.iter().sum::<f64>() / den.iter().sum::<f64>();
        let pair = f2.pair_index(0, 1);
        assert_relative_eq!(fst.block_values(pair)[0], expected, max_relative = 1e-12);
    }

    #[test]
    fn test_fst_unavailable_for_a_single_allele_population() {
        // One allele per SNP: p(1-p)/(n-1) is 0/0 everywhere, so PopA has no
        // FST blocks; the other pair is unaffected.
        let rows = [
            (0.00, [(0.0, 1), (0.5, 10), (0.3, 10)]),
            (0.10, [(1.0, 1), (0.2, 10), (0.4, 10)]),
        ];
        let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
        assert!(f2.fst("PopA", "PopB").is_none());
        assert!(f2.fst("PopB", "PopC").is_some());
    }

    #[test]
    fn test_blocks_without_data_for_some_pair_are_dropped() {
        // maxmiss = 1/3 keeps SNPs missing in one population. In the first
        // block PopA and PopC never have data at the same SNP, so admixtools2
        // f2_from_geno(remove_na = TRUE) drops the block for every pair.
        let rows = [
            (0.00, [(0.1, 10), (0.5, 10), (0.0, 0)]),
            (0.01, [(0.0, 0), (0.5, 10), (0.3, 10)]),
            (0.10, [(0.2, 10), (0.6, 10), (0.4, 10)]),
            (0.11, [(0.3, 10), (0.7, 10), (0.5, 10)]),
        ];
        let opts = F2Options {
            max_missing: 1.0 / 3.0,
            poly_only: true,
        };
        let f2 = F2Blocks::from_allele_frequencies_with_options(
            &af_from(&rows),
            &BlockDefinition::new(0.05),
            &opts,
        );
        assert_eq!(f2.block_lengths, vec![2]);
        assert_eq!(f2.n_blocks, 1);
        assert_eq!(f2.block_sums.nrows(), 1);
    }

    #[test]
    fn test_blocks_start_from_kept_snps() {
        // The SNP at 0.00 is dropped, so the first block starts at 0.03 and
        // spans 0.03..0.08 (admixtools2 builds blocks after filtering).
        let rows = [
            (0.00, [(0.1, 10), (0.0, 0), (0.3, 10)]),
            (0.03, [(0.1, 10), (0.5, 10), (0.3, 10)]),
            (0.06, [(0.2, 10), (0.5, 10), (0.3, 10)]),
            (0.09, [(0.3, 10), (0.5, 10), (0.3, 10)]),
        ];
        let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
        assert_eq!(f2.block_lengths, vec![2, 1]);
    }

    #[test]
    fn test_f2_estimate_and_se_weight_blocks_by_length() {
        // Two blocks of 3 and 1 SNPs. As in admixtools2 f2(): the estimate
        // is the length-weighted mean of block means, and the variance is
        // the weighted block jackknife mean over blocks of
        // (h - 1)(est - loo)^2, h = total / length.
        let pop = |p: f64| [(p, 10), (0.5, 10), (0.5, 10)];
        let rows = [
            (0.00, pop(0.1)),
            (0.01, pop(0.2)),
            (0.02, pop(0.3)),
            (0.10, pop(0.9)),
        ];
        let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
        assert_eq!(f2.block_lengths, vec![3, 1]);

        let snp = |i: usize| snp_f2(rows[i].1[0], rows[i].1[1]);
        let x1 = (snp(0) + snp(1) + snp(2)) / 3.0;
        let x2 = snp(3);
        let est = (3.0 * x1 + x2) / 4.0;
        // Leaving out block 1 leaves block 2, and vice versa.
        let var = ((4.0 / 3.0 - 1.0) * (est - x2).powi(2) + (4.0 - 1.0) * (est - x1).powi(2)) / 2.0;

        let (got_est, got_se, _) = f2.f2("PopA", "PopB").unwrap();
        assert_relative_eq!(got_est, est, max_relative = 1e-12);
        assert_relative_eq!(got_se, var.sqrt(), max_relative = 1e-12);
    }

    #[test]
    fn test_from_block_means_round_trips_the_block_values() {
        let af = make_test_data();
        let f2 = F2Blocks::from_allele_frequencies(&af, &BlockDefinition::new(0.005));
        let names: Vec<String> = f2.populations.iter().map(|p| p.to_string()).collect();
        let grid = |f: &dyn Fn(&str, &str) -> Vec<f64>| -> Vec<Vec<Vec<f64>>> {
            names
                .iter()
                .map(|a| names.iter().map(|b| f(a, b)).collect())
                .collect()
        };
        let means = grid(&|a, b| f2.get_block_f2_values(a, b).unwrap().to_vec());
        let counts = grid(&|a, b| f2.block_counts(a, b).unwrap());
        let back = F2Blocks::from_block_means(
            f2.populations.clone(),
            f2.block_lengths.clone(),
            &means,
            &counts,
        )
        .unwrap();
        assert_eq!(back.block_weights, f2.block_weights);
        for a in &names {
            for b in &names {
                let (x, y) = (
                    back.get_block_f2_values(a, b).unwrap(),
                    f2.get_block_f2_values(a, b).unwrap(),
                );
                for (u, v) in x.iter().zip(y.iter()) {
                    assert!(
                        (u - v).abs() <= 1e-15 * v.abs() || (u.is_nan() && v.is_nan()),
                        "{a} {b}"
                    );
                }
            }
        }
        assert_eq!(back.f2("PopA", "PopB"), f2.f2("PopA", "PopB"));
        assert!(
            F2Blocks::from_block_means(f2.populations.clone(), vec![1], &means, &counts).is_err()
        );
    }

    #[test]
    fn test_subset_populations_keeps_their_statistics() {
        let af = make_test_data();
        let f2 = F2Blocks::from_allele_frequencies(&af, &BlockDefinition::new(0.005));
        let sub = f2.subset_populations(&["PopC", "PopA"]).unwrap();
        assert_eq!(sub.n_populations(), 2);
        assert_eq!(sub.f2("PopA", "PopC"), f2.f2("PopA", "PopC"));
        assert_eq!(
            sub.get_block_f2_values("PopC", "PopC"),
            f2.get_block_f2_values("PopC", "PopC")
        );
        assert_eq!(sub.block_het_sums.column(0), f2.block_het_sums.column(2));
        assert!(f2.subset_populations(&["PopZ"]).is_err());
    }

    #[test]
    fn test_select_blocks_keeps_repeats_and_order() {
        let af = make_test_data();
        let f2 = F2Blocks::from_allele_frequencies(&af, &BlockDefinition::new(0.005));
        assert!(f2.n_blocks >= 3, "{}", f2.n_blocks);
        let all: Vec<usize> = (0..f2.n_blocks).collect();
        let same = f2.select_blocks(&all).unwrap();
        assert_eq!(same.block_sums, f2.block_sums);
        assert_eq!(same.block_lengths, f2.block_lengths);
        assert!(same.fst.is_none());

        let boot = f2.select_blocks(&[2, 0, 2]).unwrap();
        assert_eq!(boot.n_blocks, 3);
        assert_eq!(
            boot.block_lengths,
            vec![
                f2.block_lengths[2],
                f2.block_lengths[0],
                f2.block_lengths[2]
            ]
        );
        assert_eq!(boot.block_weights.row(0), f2.block_weights.row(2));
        assert_eq!(boot.block_het_sums.row(1), f2.block_het_sums.row(0));
        assert_eq!(boot.block_raw_sums.row(2), f2.block_raw_sums.row(2));

        assert!(f2.select_blocks(&[f2.n_blocks]).is_err());
        assert!(f2.select_blocks(&[]).is_err());
    }

    #[test]
    fn test_f2_blocks_creation() {
        let af = make_test_data();
        let block_def = BlockDefinition::new(0.01); // 1 cM blocks

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        assert_eq!(f2_blocks.n_populations(), 3);
        assert!(f2_blocks.n_blocks > 0);

        // Should have 6 pairs: AA, AB, AC, BB, BC, CC
        let n_pairs = 3 * 4 / 2; // n*(n+1)/2 = 6
        assert_eq!(f2_blocks.block_sums.ncols(), n_pairs);
    }

    #[test]
    fn test_f2_self_pair_near_zero() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f2(A,A) should be close to 0
        let (f2_aa, _, _) = f2_blocks.f2("PopA", "PopA").unwrap();
        assert!(
            f2_aa.abs() < 0.01,
            "f2(A,A) should be near 0, got {}",
            f2_aa
        );
    }

    #[test]
    fn test_f2_symmetry() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f2(A,B) should equal f2(B,A)
        let (f2_ab, _, _) = f2_blocks.f2("PopA", "PopB").unwrap();
        let (f2_ba, _, _) = f2_blocks.f2("PopB", "PopA").unwrap();
        assert_relative_eq!(f2_ab, f2_ba, epsilon = 1e-10);
    }

    #[test]
    fn test_f2_different_populations() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // f2(A,B) should be positive (different populations)
        let (f2_ab, std_err, n_snps) = f2_blocks.f2("PopA", "PopB").unwrap();
        assert!(f2_ab > 0.0, "f2(A,B) should be positive, got {}", f2_ab);
        assert!(std_err >= 0.0);
        assert!(n_snps > 0);
    }

    #[test]
    fn f2_matrix_holds_the_f2_estimates() {
        // Blocks of 1 and 3 SNPs: the length-weighted estimate is 0.4, the
        // plain mean of block means 0.3.
        let block = |a: f64, b: f64| {
            vec![
                vec![vec![0.0; 2], vec![a, b]],
                vec![vec![a, b], vec![0.0; 2]],
            ]
        };
        let f2 = F2Blocks::from_block_means(
            vec![PopulationId::new("A"), PopulationId::new("B")],
            vec![1, 3],
            &block(0.1, 0.5),
            &block(1.0, 3.0),
        )
        .unwrap();
        let (estimate, _, _) = f2.f2("A", "B").unwrap();
        assert_relative_eq!(estimate, 0.4, epsilon = 1e-12);
        assert_eq!(f2.f2_matrix()[[0, 1]], estimate);
        assert_eq!(f2.f2_matrix()[[1, 0]], estimate);
    }

    #[test]
    fn test_f2_matrix() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);
        let matrix = f2_blocks.f2_matrix();

        assert_eq!(matrix.shape(), &[3, 3]);

        // Print matrix for debugging
        println!("F2 matrix:\n{:?}", matrix);

        // Diagonal should be near 0 (but can be slightly negative due to bias correction)
        for i in 0..3 {
            assert!(
                matrix[[i, i]].abs() < 0.1,
                "Diagonal element {} = {} should be near 0",
                i,
                matrix[[i, i]]
            );
        }

        // Should be symmetric
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(matrix[[i, j]], matrix[[j, i]], epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_pair_index() {
        let af = make_test_data();
        let block_def = BlockDefinition::default();
        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        // Verify pair indexing is consistent
        // For 3 populations: pairs are (0,0), (0,1), (0,2), (1,1), (1,2), (2,2)
        // Indices should be 0, 1, 2, 3, 4, 5
        assert_eq!(f2_blocks.pair_index(0, 0), 0);
        assert_eq!(f2_blocks.pair_index(0, 1), 1);
        assert_eq!(f2_blocks.pair_index(0, 2), 2);
        assert_eq!(f2_blocks.pair_index(1, 1), 3);
        assert_eq!(f2_blocks.pair_index(1, 2), 4);
        assert_eq!(f2_blocks.pair_index(2, 2), 5);

        // Symmetric access
        assert_eq!(f2_blocks.pair_index(1, 0), f2_blocks.pair_index(0, 1));
        assert_eq!(f2_blocks.pair_index(2, 0), f2_blocks.pair_index(0, 2));
    }

    #[test]
    fn test_jackknife_std_err() {
        let af = make_test_data();
        let block_def = BlockDefinition::new(0.005); // Small blocks for more jackknife samples

        let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

        let (_, std_err, _) = f2_blocks.f2("PopA", "PopB").unwrap();

        // Standard error should be non-negative
        assert!(std_err >= 0.0);

        // With multiple blocks, should have positive standard error
        if f2_blocks.n_blocks > 2 {
            // May still be 0 if all blocks have identical estimates
        }
    }
}
