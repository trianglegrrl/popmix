//! f4-, D- and f3-statistics computed directly from genotypes.
//!
//! These are the statistics admixtools2's `qpdstat()` and `qp3pop()` report
//! when given a genotype prefix with `allsnps = TRUE` (their documentation),
//! which differ from the statistics computed from f2 blocks:
//!
//! - every SNP passed in is used (no polymorphism or missingness filter), and
//!   the jackknife blocks are formed over all of them;
//! - each statistic uses the SNPs at which all of its populations have data;
//! - allele frequencies count every sample as diploid (the genotype sum over
//!   twice the number of non-missing samples). Pseudo-haploid samples are
//!   taken into account only in the f3 bias correction and heterozygosity.
//!
//! D-statistics are `(a-b)(c-d) / ((a+b-2ab)(c+d-2cd))`, with numerator and
//! denominator averaged per block and the ratio taken of the leave-one-block-out
//! means. f3 is `(a-b)(a-c)` less the target's bias correction, normalized by
//! the target's heterozygosity (or, in outgroup mode, not normalized).
//!
//! The jackknife weights each block by its number of SNPs used. Each
//! statistic's leave-one-block-out value is the leave-one-out mean (for D and
//! f3, the ratio of the leave-one-out numerator and denominator means); its
//! total is the mean of those values weighted by 1 − m_k / n, and the
//! estimate and standard error are the weighted block jackknife's
//! ([`jackknife_from_loo`]).
//!
//! Provenance: the jackknife is written from the clean-room specifications
//! `docs/cleanroom/small-helpers.md` (section 2) and
//! `docs/cleanroom/jackknife.md` (Busing et al. 1999; Patterson et al.
//! 2012), and checked against admixtools2's black-box outputs.

use rayon::prelude::*;
use statrs::distribution::{ContinuousCDF, Normal};

use crate::error::{PopmixError, Result};
use crate::fstats::blocks::BlockDefinition;
use crate::fstats::jackknife::{jackknife_from_loo, leave_one_out, z_score};
use crate::genotype_matrix::GenotypeMatrix;
use crate::genotypes::Snp;

/// A source of genotype rows, one SNP at a time, in SNP order.
pub trait GenotypeRows {
    /// Fill `out` with the genotype codes (0, 1, 2; 3 for missing) of the
    /// selected samples at the next SNP.
    fn next_row(&mut self, out: &mut [u8]) -> Result<()>;

    /// Fill `out` with the next `out.len() / width` rows, `width` codes each,
    /// as that many calls of [`next_row`](Self::next_row) would; sources
    /// that can decode rows in parallel do.
    fn next_rows(&mut self, out: &mut [u8], width: usize) -> Result<()> {
        for row in out.chunks_exact_mut(width.max(1)) {
            self.next_row(row)?;
        }
        Ok(())
    }
}

/// Rows of an in-memory genotype matrix: the given SNPs, for the given samples.
pub struct MatrixRows<'a> {
    matrix: &'a GenotypeMatrix,
    snps: Vec<usize>,
    samples: Vec<usize>,
    next: usize,
}

impl<'a> MatrixRows<'a> {
    pub fn new(matrix: &'a GenotypeMatrix, snps: Vec<usize>, samples: Vec<usize>) -> Self {
        MatrixRows {
            matrix,
            snps,
            samples,
            next: 0,
        }
    }
}

impl GenotypeRows for MatrixRows<'_> {
    fn next_row(&mut self, out: &mut [u8]) -> Result<()> {
        let snp = *self.snps.get(self.next).ok_or_else(|| {
            PopmixError::Internal("read past the last SNP of the genotype matrix".to_string())
        })?;
        let row = self.matrix.snp_raw(snp);
        for (o, &s) in out.iter_mut().zip(&self.samples) {
            *o = row[s].min(3);
        }
        self.next += 1;
        Ok(())
    }
}

/// The samples that enter the statistics and their populations.
#[derive(Debug, Clone)]
pub struct GenoPopulations {
    /// Population names; statistics refer to populations by index here.
    pub names: Vec<String>,
    /// For each selected sample (in row order), its population index.
    pub sample_pop: Vec<usize>,
    /// For each selected sample, its ploidy (1 for pseudo-haploid, 2 diploid).
    pub ploidy: Vec<u8>,
}

/// f4 numerator alone, or the D-statistic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum F4Mode {
    F4,
    D,
}

/// Options for f3 from genotypes (admixtools2 `qp3pop` defaults).
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct F3GenoOptions {
    /// Subtract the target's sampling bias, `h / n` (`apply_corr`).
    pub apply_corr: bool,
    /// Leave f3 unnormalized (`outgroupmode = TRUE`); otherwise it is divided
    /// by the target's heterozygosity.
    pub outgroup_mode: bool,
}

impl Default for F3GenoOptions {
    fn default() -> Self {
        F3GenoOptions {
            apply_corr: true,
            outgroup_mode: false,
        }
    }
}

/// Per-block data for one statistic: the mean numerator and denominator over the block's SNPs
/// used, and their number. Means are NaN for blocks with no SNP used; the
/// denominator is NaN throughout for f4.
#[derive(Debug, Clone)]
pub struct GenoStatBlocks {
    pub num: Vec<f64>,
    pub den: Vec<f64>,
    pub n: Vec<u64>,
}

/// One statistic: jackknife estimate, standard error, z, two-sided p and the
/// total number of SNPs used, with its per-block data.
#[derive(Debug, Clone)]
pub struct GenoStat {
    pub est: f64,
    pub se: f64,
    pub z: f64,
    pub p: f64,
    pub n: u64,
    pub blocks: GenoStatBlocks,
}

/// Statistics from genotypes, with the block lengths (SNPs per block) used.
#[derive(Debug, Clone)]
pub struct GenoStats {
    pub block_lengths: Vec<usize>,
    pub stats: Vec<GenoStat>,
}

/// Per-population sums over the samples at each SNP of one block.
struct BlockTable {
    n_pops: usize,
    /// Allele frequency counting every sample as diploid; NaN without data.
    af: Vec<f64>,
    /// Ploidy-adjusted counts of the two alleles (admixtools2 `ref`, `alt`).
    ref_count: Vec<f64>,
    alt_count: Vec<f64>,
}

impl BlockTable {
    #[inline]
    fn idx(&self, snp: usize, pop: usize) -> usize {
        snp * self.n_pops + pop
    }

    /// The table of `n_snps` rows of genotype codes (`raw`, one row of the
    /// selected samples per SNP).
    fn from_rows(raw: &[u8], pops: &GenoPopulations, n_snps: usize) -> Self {
        let n_pops = pops.names.len();
        let mut table = BlockTable {
            n_pops,
            af: vec![f64::NAN; n_snps * n_pops],
            ref_count: vec![0.0; n_snps * n_pops],
            alt_count: vec![0.0; n_snps * n_pops],
        };
        let n_samples = pops.sample_pop.len();
        let mut sum = vec![0.0f64; n_pops];
        let mut count = vec![0.0f64; n_pops];
        for snp in 0..n_snps {
            let row = &raw[snp * n_samples..(snp + 1) * n_samples];
            sum.fill(0.0);
            count.fill(0.0);

            for ((&g, &p), &ploidy) in row.iter().zip(&pops.sample_pop).zip(&pops.ploidy) {
                if g > 2 {
                    continue;
                }
                let g = f64::from(g);
                sum[p] += g;
                count[p] += 1.0;
                let i = snp * n_pops + p;
                if ploidy == 1 {
                    table.ref_count[i] += g / 2.0;
                    table.alt_count[i] += (2.0 - g) / 2.0;
                } else {
                    table.ref_count[i] += g;
                    table.alt_count[i] += 2.0 - g;
                }
            }
            for p in 0..n_pops {
                table.af[snp * n_pops + p] = sum[p] / count[p] / 2.0;
            }
        }
        table
    }

    fn n_snps(&self) -> usize {
        self.af.len() / self.n_pops.max(1)
    }
}

/// Validate population indices and set up the jackknife blocks.
fn blocks_for(snps: &[Snp], block_size: f64) -> Vec<usize> {
    BlockDefinition::new(block_size)
        .assign_blocks(snps)
        .block_lengths
}

fn check_indices(pops: &GenoPopulations, stats: &[&[usize]]) -> Result<()> {
    let n = pops.names.len();
    if pops.sample_pop.len() != pops.ploidy.len() {
        return Err(PopmixError::DimensionMismatch {
            context: "sample ploidies".to_string(),
            expected: format!("{}, one per sample", pops.sample_pop.len()),
            actual: pops.ploidy.len().to_string(),
        });
    }
    if let Some(&bad) = pops.sample_pop.iter().find(|&&p| p >= n) {
        return Err(PopmixError::out_of_range("sample population index", bad, n));
    }
    for s in stats {
        if let Some(&bad) = s.iter().find(|&&p| p >= n) {
            return Err(PopmixError::out_of_range("population index", bad, n));
        }
    }
    Ok(())
}

/// Read the blocks one at a time and accumulate, for every statistic, the
/// per-block (numerator sum, denominator sum, SNP count) that `per_snp` gives.
fn accumulate<F>(
    rows: &mut dyn GenotypeRows,
    pops: &GenoPopulations,
    block_lengths: &[usize],
    n_stats: usize,
    per_block: F,
) -> Result<Vec<Vec<(f64, f64, u64)>>>
where
    F: Fn(&BlockTable, usize) -> (f64, f64, u64) + Sync,
{
    let n_samples = pops.sample_pop.len();
    let mut out = vec![Vec::with_capacity(block_lengths.len()); n_stats];
    // the rows of a batch of blocks are read in order, then each block is
    // summed on its own, the blocks in parallel
    for batch in block_lengths.chunks(BATCH_BLOCKS) {
        let mut raw: Vec<Vec<u8>> = Vec::with_capacity(batch.len());
        for &len in batch {
            let mut block = vec![0u8; len * n_samples];
            for row in block.chunks_exact_mut(n_samples.max(1)).take(len) {
                rows.next_row(row)?;
            }
            raw.push(block);
        }
        let sums: Vec<Vec<(f64, f64, u64)>> = raw
            .par_iter()
            .zip(batch)
            .map(|(block, &len)| {
                let table = BlockTable::from_rows(block, pops, len);
                (0..n_stats).map(|j| per_block(&table, j)).collect()
            })
            .collect();
        for block in sums {
            for (o, s) in out.iter_mut().zip(block) {
                o.push(s);
            }
        }
    }
    Ok(out)
}

/// Jackknife blocks read at a time by [`accumulate`].
const BATCH_BLOCKS: usize = 32;

/// f4- or D-statistics `(pop1, pop2; pop3, pop4)` from genotypes.
///
/// `rows` must yield exactly one row per SNP in `snps`, in order.
pub fn f4_from_genotypes(
    rows: &mut dyn GenotypeRows,
    snps: &[Snp],
    pops: &GenoPopulations,
    quads: &[[usize; 4]],
    mode: F4Mode,
    block_size: f64,
) -> Result<GenoStats> {
    f4_from_genotypes_requiring(rows, snps, pops, quads, mode, block_size, &[])
}

/// [`f4_from_genotypes`] on only the SNPs where every population of
/// `require` has data too (admixtools2 `allsnps = FALSE`: a model's SNPs).
pub fn f4_from_genotypes_requiring(
    rows: &mut dyn GenotypeRows,
    snps: &[Snp],
    pops: &GenoPopulations,
    quads: &[[usize; 4]],
    mode: F4Mode,
    block_size: f64,
    require: &[usize],
) -> Result<GenoStats> {
    let mut refs: Vec<&[usize]> = quads.iter().map(|q| &q[..]).collect();
    refs.push(require);
    check_indices(pops, &refs)?;
    let block_lengths = blocks_for(snps, block_size);
    let sums = accumulate(rows, pops, &block_lengths, quads.len(), |t, j| {
        let [i1, i2, i3, i4] = quads[j];
        let (mut num, mut den, mut n) = (0.0, 0.0, 0u64);
        for s in 0..t.n_snps() {
            if require.iter().any(|&p| !t.af[t.idx(s, p)].is_finite()) {
                continue;
            }
            let (a, b, c, d) = (
                t.af[t.idx(s, i1)],
                t.af[t.idx(s, i2)],
                t.af[t.idx(s, i3)],
                t.af[t.idx(s, i4)],
            );
            let x = (a - b) * (c - d);
            if x.is_finite() {
                num += x;
                n += 1;
                if mode == F4Mode::D {
                    den += (a + b - 2.0 * a * b) * (c + d - 2.0 * c * d);
                }
            }
        }
        (num, den, n)
    })?;

    let stats = sums
        .into_iter()
        .map(|blocks| {
            let n: Vec<u64> = blocks.iter().map(|b| b.2).collect();
            let num: Vec<f64> = blocks.iter().map(|b| b.0 / b.2 as f64).collect();
            let den: Vec<f64> = match mode {
                F4Mode::F4 => vec![f64::NAN; blocks.len()],
                F4Mode::D => blocks.iter().map(|b| b.1 / b.2 as f64).collect(),
            };
            let loo = match mode {
                F4Mode::F4 => block_loo(&num, &n),
                F4Mode::D => loo_ratio(&block_loo(&num, &n), &block_loo(&den, &n)),
            };
            summarize(&loo, GenoStatBlocks { num, den, n })
        })
        .collect();
    Ok(GenoStats {
        block_lengths,
        stats,
    })
}

/// f3-statistics `(pop1; pop2, pop3)` from genotypes, `pop1` the target.
///
/// `rows` must yield exactly one row per SNP in `snps`, in order.
pub fn f3_from_genotypes(
    rows: &mut dyn GenotypeRows,
    snps: &[Snp],
    pops: &GenoPopulations,
    triples: &[[usize; 3]],
    options: F3GenoOptions,
    block_size: f64,
) -> Result<GenoStats> {
    let refs: Vec<&[usize]> = triples.iter().map(|q| &q[..]).collect();
    check_indices(pops, &refs)?;
    let block_lengths = blocks_for(snps, block_size);
    let needs_h = options.apply_corr || !options.outgroup_mode;
    let sums = accumulate(rows, pops, &block_lengths, triples.len(), |t, j| {
        let [i1, i2, i3] = triples[j];
        let corrected = options.apply_corr && i1 != i2 && i1 != i3;
        let (mut num, mut den, mut n) = (0.0, 0.0, 0u64);
        for s in 0..t.n_snps() {
            let (a, b, c) = (t.af[t.idx(s, i1)], t.af[t.idx(s, i2)], t.af[t.idx(s, i3)]);
            let mut x = (a - b) * (a - c);
            if !x.is_finite() {
                continue;
            }
            n += 1;
            let mut h = f64::NAN;
            if needs_h {
                let (r, al) = (t.ref_count[t.idx(s, i1)], t.alt_count[t.idx(s, i1)]);
                let tot = r + al;
                h = r * al / (tot * (tot - 1.0));
                if corrected {
                    x -= h / tot;
                }
            }
            if !x.is_finite() {
                continue;
            }
            num += x;
            if options.outgroup_mode {
                den += 1.0;
            } else if h.is_finite() {
                den += 2.0 * h;
            }
        }
        (num, den, n)
    })?;

    let stats = sums
        .into_iter()
        .map(|blocks| {
            let n: Vec<u64> = blocks.iter().map(|b| b.2).collect();
            let num: Vec<f64> = blocks.iter().map(|b| b.0 / b.2 as f64).collect();
            let den: Vec<f64> = blocks.iter().map(|b| b.1 / b.2 as f64).collect();
            let loo = loo_ratio(&block_loo(&num, &n), &block_loo(&den, &n));
            summarize(&loo, GenoStatBlocks { num, den, n })
        })
        .collect();
    Ok(GenoStats {
        block_lengths,
        stats,
    })
}

/// Leave-one-block-out values of a block mean, each block weighted by its
/// number of SNPs used; a block without any gets a non-finite value.
fn block_loo(values: &[f64], n: &[u64]) -> Vec<f64> {
    let values: Vec<f64> = values
        .iter()
        .zip(n)
        .map(|(&v, &k)| if k > 0 { v } else { f64::NAN })
        .collect();
    leave_one_out(&values, &weights(n))
}

/// Leave-one-out values of a ratio: numerator's over denominator's.
fn loo_ratio(num: &[f64], den: &[f64]) -> Vec<f64> {
    num.iter().zip(den).map(|(a, b)| a / b).collect()
}

fn weights(n: &[u64]) -> Vec<f64> {
    n.iter().map(|&k| k as f64).collect()
}

/// Jackknife estimate, standard error, z and two-sided normal p from the
/// leave-one-out values; `n` counts the SNPs of the blocks with a finite
/// leave-one-out value. Everything is NaN when no block has one.
fn summarize(loo: &[f64], blocks: GenoStatBlocks) -> GenoStat {
    let (est, se) = jackknife_from_loo(loo, &weights(&blocks.n)).unwrap_or((f64::NAN, f64::NAN));
    let z = z_score(est, se);
    let p = 2.0 * Normal::standard().sf(z.abs());
    let n = loo
        .iter()
        .zip(&blocks.n)
        .filter(|(t, _)| t.is_finite())
        .map(|(_, &k)| k)
        .sum();
    GenoStat {
        est,
        se,
        z,
        p,
        n,
        blocks,
    }
}

#[cfg(test)]
mod tests;
