//! How far LD in the test population is correlated with LD in a reference,
//! as ALDER (`Alder::find_ld_corr_stops`, with `approx_ld_corr: YES`).
//!
//! For distance bins from 0.05 cM up to 2 cM, the LD (genotype covariance)
//! of each SNP pair in the bin is computed in both populations, and their
//! uncentered correlation, sum(xy) / sqrt(sum(x^2) sum(y^2)), with a jackknife
//! over chromosomes. The fit starts where the correlation has lost
//! significance (one-sided p > 0.05) for the second time. Pairs are visited
//! in rounds of every 64th SNP, and a bin stops early once the correlation is
//! clearly significant, as ALDER does to save time; the stop is taken at two
//! bin sizes (0.05 and 0.1 cM by default) and the larger kept.

use rayon::prelude::*;

use super::fit::Jackknifed;
use crate::pca::PackedGenotypes;

const SIGNIF_THRESH: f64 = 0.05;
const LIM_FAILURES: usize = 2;
const MAX_LD_DIST: f64 = 0.02;
const ROUNDS: usize = 6;

/// Genotype covariance of two SNPs over the samples called at both (n - 1
/// denominator); NaN with fewer than two.
fn ld(g: &PackedGenotypes, s1: usize, s2: usize) -> f64 {
    let (mut sx, mut sy, mut sxy, mut n) = (0i64, 0i64, 0i64, 0i64);
    for i in 0..g.n_samples() {
        let (x, y) = (g.get(s1, i), g.get(s2, i));
        if x <= 2 && y <= 2 {
            let (x, y) = (x as i64, y as i64);
            sx += x;
            sy += y;
            sxy += x * y;
            n += 1;
        }
    }
    if n <= 1 {
        f64::NAN
    } else {
        (sxy as f64 - (sx * sy) as f64 / n as f64) / (n - 1) as f64
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct Sums {
    xy: f64,
    x2: f64,
    y2: f64,
}

/// Jackknifed uncentered correlation over chromosomes.
fn cosine(per_chrom: &[Sums]) -> Jackknifed {
    let at = |skip: Option<usize>| {
        let (mut xy, mut x2, mut y2) = (0.0, 0.0, 0.0);
        for (c, s) in per_chrom.iter().enumerate() {
            if Some(c) != skip {
                xy += s.xy;
                x2 += s.x2;
                y2 += s.y2;
            }
        }
        xy / (x2 * y2).sqrt()
    };
    Jackknifed {
        value: at(None),
        reps: (0..per_chrom.len()).map(|c| at(Some(c))).collect(),
    }
}

fn not_significant(j: &Jackknifed) -> bool {
    let p = statrs::function::erf::erfc(j.value / j.se() / std::f64::consts::SQRT_2);
    // a NaN estimate is not significant either
    p > SIGNIF_THRESH || j.value <= 0.0 || j.value.is_nan()
}

/// The distance (Morgans) to which test and reference LD are correlated:
/// infinite if still correlated at 2 cM. `chroms` holds each chromosome's
/// SNP indices (rows of both genotype sets) in order, `pos` their genetic
/// positions.
pub fn ld_corr_stop(
    test: &PackedGenotypes,
    reference: &PackedGenotypes,
    chroms: &[Vec<usize>],
    pos: &[f64],
    binsize0: f64,
) -> f64 {
    let mut stop = f64::INFINITY;
    let mut binsize = binsize0;
    let stride = 1usize << ROUNDS;
    loop {
        let numbins = (MAX_LD_DIST / binsize) as usize;
        let mut failures = 0;
        let mut pass_stop = None;
        for b in 1..numbins {
            let (lo, hi) = (b as f64 * binsize, (b + 1) as f64 * binsize);
            let mut sums = vec![Sums::default(); chroms.len()];
            let mut checks_left = ROUNDS + 1;
            for offset in 0..stride {
                // each chromosome's sums in ALDER's order: pair by pair
                sums.par_iter_mut()
                    .zip(chroms.par_iter())
                    .for_each(|(acc, snps)| {
                        let mut s1 = offset;
                        while s1 < snps.len() {
                            let p1 = pos[snps[s1]];
                            let first = snps[s1..].partition_point(|&s| pos[s] < p1 + lo) + s1;
                            for &s2 in snps[first..].iter().take_while(|&&s| pos[s] < p1 + hi) {
                                let (x, y) = (ld(test, snps[s1], s2), ld(reference, snps[s1], s2));
                                if !x.is_nan() && !y.is_nan() {
                                    acc.xy += x * y;
                                    acc.x2 += x * x;
                                    acc.y2 += y * y;
                                }
                            }
                            s1 += stride;
                        }
                    });
                // after rounds 1, 2, 4, ..., 64: stop once clearly significant
                if (offset + 1) & offset == 0 {
                    let cos = cosine(&sums);
                    let p = statrs::function::erf::erfc(
                        cos.value / cos.se() / std::f64::consts::SQRT_2,
                    );
                    if p * (checks_left as f64) < SIGNIF_THRESH {
                        break;
                    }
                    checks_left -= 1;
                }
            }
            if not_significant(&cosine(&sums)) {
                failures += 1;
                if failures == LIM_FAILURES {
                    pass_stop = Some(b as f64 * binsize);
                    break;
                }
            }
        }
        match pass_stop {
            Some(d) => {
                if stop.is_infinite() || d > stop {
                    stop = d;
                }
            }
            None => stop = f64::INFINITY,
        }
        binsize *= 2.0;
        if binsize >= 0.002 {
            break;
        }
    }
    stop
}
