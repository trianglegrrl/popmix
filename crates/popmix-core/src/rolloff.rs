//! ROLLOFF (Moorjani et al. 2011; Patterson et al. 2012), as in ADMIXTOOLS
//! 7.0.2: the correlation between a signed LD statistic in an admixed
//! population and the allele frequency differences of two references, as a
//! function of genetic distance.
//!
//! For each pair of SNPs on a chromosome closer than `maxdis`, with weight
//! w = w_x w_y (w the frequency difference A - B, divided by sqrt of the
//! mean frequency in run mode 0) and both references in linkage equilibrium
//! there (LD chi-square below `chithresh`), the admixed population's LD
//! statistic z (with `zdipcorrmode`, its genotype correlation Fisher-
//! transformed, times sqrt(n - 3); otherwise a likelihood-ratio signed root)
//! is added to the pair's distance bin. Each bin reports the weighted
//! Z-score sum(w z) / sqrt(sum w^2), the regression sum(w z) / sum(w^2),
//! the correlation of w and z, and the pair count; also each of these with
//! one chromosome left out, for the jackknife.
//!
//! The date comes from fitting C + A exp(-m d) to the correlations, from
//! 0.5 cM (`fit_rolloff`), with a weighted jackknife over chromosomes, as
//! ADMIXTOOLS' rexpfit.r and dowtjack.

use rayon::prelude::*;

use crate::error::{PopmixError, Result};
use crate::genotypes::Snp;
use crate::pca::PackedGenotypes;

#[derive(Debug, Clone)]
pub struct RolloffOptions {
    /// Distance bin (Morgans).
    pub binsize: f64,
    /// Largest distance (Morgans).
    pub maxdis: f64,
    /// 0: weights (pA - pB) / sqrt(mean frequency); 1: pA - pB.
    pub runmode: u8,
    /// Fisher-transformed genotype correlation as the LD statistic.
    pub zdipcorr: bool,
    /// Leave out pairs in LD in either reference (chi-square at least this);
    /// 0 for none.
    pub chithresh: f64,
    /// Copies of each allele needed in the admixed population.
    pub mincount: usize,
    /// Samples of each reference called at both SNPs of a pair.
    pub minparentcount: usize,
}

impl Default for RolloffOptions {
    fn default() -> Self {
        RolloffOptions {
            binsize: 0.0005,
            maxdis: 0.05,
            runmode: 1,
            zdipcorr: true,
            chithresh: 6.0,
            mincount: 5,
            minparentcount: 10,
        }
    }
}

/// Sufficient statistics of (w, z) pairs (ADMIXTOOLS `CORR`).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct CorrSums {
    pub s0: f64,
    pub s1: f64,
    pub s2: f64,
    pub s11: f64,
    pub s12: f64,
    pub s22: f64,
}

impl CorrSums {
    fn add(&mut self, x1: f64, x2: f64, n: f64) {
        self.s0 += n;
        self.s1 += x1 * n;
        self.s2 += x2 * n;
        self.s11 += x1 * x1 * n;
        self.s12 += x1 * x2 * n;
        self.s22 += x2 * x2 * n;
    }

    fn merge(&mut self, o: &CorrSums) {
        self.s0 += o.s0;
        self.s1 += o.s1;
        self.s2 += o.s2;
        self.s11 += o.s11;
        self.s12 += o.s12;
        self.s22 += o.s22;
    }

    fn minus(&self, o: &CorrSums) -> CorrSums {
        CorrSums {
            s0: self.s0 - o.s0,
            s1: self.s1 - o.s1,
            s2: self.s2 - o.s2,
            s11: self.s11 - o.s11,
            s12: self.s12 - o.s12,
            s22: self.s22 - o.s22,
        }
    }

    /// ADMIXTOOLS `calccorr` (mean taken off unless `uncentered`): the
    /// correlation, and with `ztrans` the Fisher-transformed Z; None with no
    /// data (or, for `ztrans`, fewer than 4).
    fn corr(&self, uncentered: bool, ztrans: bool) -> Option<(f64, f64)> {
        if self.s0 < 0.5 {
            return None;
        }
        let yn = self.s0;
        let (m1, m2) = if uncentered {
            (0.0, 0.0)
        } else {
            (self.s1 / yn, self.s2 / yn)
        };
        let v11 = (self.s11 - yn * m1 * m1) / yn;
        let v12 = (self.s12 - yn * m1 * m2) / yn;
        let v22 = (self.s22 - yn * m2 * m2) / yn;
        let r = v12 / (v11 * v22 + 1e-20).sqrt();
        let mut z = yn.sqrt() * r;
        if ztrans {
            if yn < 4.0 {
                return None;
            }
            let y = r.clamp(-0.9, 0.9);
            z = (yn - 3.0).sqrt() * 0.5 * ((1.0 + y) / (1.0 - y)).ln();
        }
        Some((r, z))
    }
}

/// Sum of x ln y over 9 (or 4) cells (ADMIXTOOLS `vldot`).
fn vldot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b.ln()).sum()
}

/// Expected haplotype counts from a 3x3 genotype table given haplotype
/// frequencies p (ADMIXTOOLS `mk4from9`); returns the log-likelihood.
fn mk4from9(xd: &mut [f64; 4], xc: &[f64; 9], p: &[f64; 4]) -> f64 {
    let mut pp = [1e-20; 9];
    for e1 in 0..2 {
        for e2 in 0..2 {
            for f1 in 0..2 {
                for f2 in 0..2 {
                    pp[3 * (e1 + e2) + f1 + f2] += p[2 * e1 + f1] * p[2 * e2 + f2];
                }
            }
        }
    }
    let ylike = vldot(xc, &pp);
    *xd = [0.0; 4];
    for e1 in 0..2 {
        for e2 in 0..2 {
            for f1 in 0..2 {
                for f2 in 0..2 {
                    let (k1, k2, k3) = (2 * e1 + f1, 2 * e2 + f2, 3 * (e1 + e2) + f1 + f2);
                    let y = p[k1] * p[k2] / pp[k3];
                    xd[k1] += y * xc[k3];
                    xd[k2] += y * xc[k3];
                }
            }
        }
    }
    ylike
}

/// Marginal allele counts from haplotype counts (ADMIXTOOLS `mk2from4`).
fn mk2from4(xd: &[f64; 4]) -> ([f64; 2], [f64; 2]) {
    let mut xe = [0.0; 2];
    let mut xf = [0.0; 2];
    for e in 0..2 {
        for f in 0..2 {
            xe[e] += xd[2 * e + f];
            xf[f] += xd[2 * e + f];
        }
    }
    (xe, xf)
}

/// Likelihood-ratio chi-square for LD from a 3x3 genotype table, with the
/// linked and unlinked haplotype frequencies (ADMIXTOOLS `lddip`, whose
/// first loop stops at `first_tol`; `zdip0` uses 0.001 there); None with
/// fewer than 2.5 samples.
fn lddip_full(xc: &[f64; 9], first_tol: f64) -> Option<(f64, [f64; 4], [f64; 4])> {
    if xc.iter().sum::<f64>() < 2.5 {
        return None;
    }
    let mut p4 = [0.25; 4];
    let mut xd = [0.0; 4];
    let mut ybase = 0.0;
    let mut yl = 0.0;
    let mut ylast = -1e40;
    for iter in 1..=10 {
        let ylike0 = mk4from9(&mut xd, xc, &p4);
        yl = ylike0;
        let (xe, xf) = mk2from4(&xd);
        if iter == 1 {
            ybase = ylike0;
        }
        let ylike = ylike0 - ybase;
        let (se, sf) = (xe[0] + xe[1], xf[0] + xf[1]);
        let pe = [xe[0] / se, xe[1] / se];
        let pf = [xf[0] / sf, xf[1] / sf];
        for e in 0..2 {
            for f in 0..2 {
                p4[2 * e + f] = pe[e] * pf[f];
            }
        }
        if ylike < ylast + first_tol {
            break;
        }
        ylast = ylike;
    }
    let p4t = p4;
    let mut ylike = 0.0;
    ylast = -1e40;
    for _ in 1..=10 {
        ylike = mk4from9(&mut xd, xc, &p4) - yl;
        let s: f64 = xd.iter().sum();
        for k in 0..4 {
            p4[k] = xd[k] / s;
        }
        if ylike < ylast + 0.001 {
            break;
        }
        ylast = ylike;
    }
    Some((2.0 * ylike, p4, p4t))
}

/// The LD chi-square (ADMIXTOOLS `lddip`), -1 with too few samples.
fn lddip(xc: &[f64; 9]) -> f64 {
    lddip_full(xc, 0.0001).map_or(-1.0, |x| x.0)
}

/// The admixed population's LD statistic for a pair (ADMIXTOOLS `zdip`).
fn zdip(xc: &[f64; 9], zdipcorr: bool) -> f64 {
    if zdipcorr {
        let mut c = CorrSums::default();
        for a in 0..3 {
            for b in 0..3 {
                c.add(a as f64, b as f64, xc[3 * a + b]);
            }
        }
        return c.corr(false, true).map_or(0.0, |(_, z)| z);
    }
    // the signed root of the chi-square (ADMIXTOOLS `zdip0`)
    match lddip_full(xc, 0.001) {
        None => -1.0,
        Some((ychi, p4, p4t)) if ychi > 0.0 => {
            let z = ychi.sqrt();
            if p4[0] < p4t[0] {
                -z
            } else {
                z
            }
        }
        _ => 0.0,
    }
}

/// The 3x3 table of genotype pairs at two SNPs over the called samples.
fn table(g: &PackedGenotypes, s1: usize, s2: usize) -> [f64; 9] {
    let mut xc = [0.0; 9];
    for i in 0..g.n_samples() {
        let (e, f) = (g.get(s1, i), g.get(s2, i));
        if e <= 2 && f <= 2 {
            xc[3 * e as usize + f as usize] += 1.0;
        }
    }
    xc
}

/// The smaller allele count in a genotype-count vector (ADMIXTOOLS `cminor`).
fn minor(x: &[f64; 3]) -> f64 {
    let y1 = x[1] + 2.0 * x[2];
    let y2 = 2.0 * (x[0] + x[1] + x[2]) - y1;
    y1.min(y2)
}

/// ROLLOFF's output per distance bin.
#[derive(Debug, Clone, PartialEq)]
pub struct RolloffBins {
    /// Upper edge of each bin (Morgans), as ROLLOFF labels it.
    pub d: Vec<f64>,
    /// sum(w z) / sqrt(sum w^2): ROLLOFF's weighted Z-score.
    pub zscore: Vec<f64>,
    /// sum(w z) / sum(w^2) (uncentered regression).
    pub regression: Vec<f64>,
    /// Correlation of w and z.
    pub corr: Vec<f64>,
    /// SNP pairs.
    pub count: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct RolloffResult {
    /// Chromosome labels, in order.
    pub chromosomes: Vec<String>,
    pub bins: RolloffBins,
    /// Bins with each chromosome left out; their `zscore` holds, as
    /// ROLLOFF's jackknife files, the left-out chromosome's pair count.
    pub jackknife: Vec<RolloffBins>,
    /// SNPs with a weight, after the filters.
    pub n_snps_used: usize,
}

fn round(x: f64) -> f64 {
    // nint: nearest integer, halves away from zero
    x.round()
}

/// ROLLOFF for an admixed population with references A and B (each a set of
/// samples at the same SNPs, `snps`, autosomes sorted by position).
pub fn rolloff(
    admixed: &PackedGenotypes,
    ref_a: &PackedGenotypes,
    ref_b: &PackedGenotypes,
    snps: &[Snp],
    options: &RolloffOptions,
) -> Result<RolloffResult> {
    crate::options::check_bins("binsize", options.binsize, "maxdis", options.maxdis)?;
    let m = snps.len();
    if admixed.n_snps() != m || ref_a.n_snps() != m || ref_b.n_snps() != m {
        return Err(PopmixError::DimensionMismatch {
            context: "ROLLOFF SNPs".to_string(),
            expected: format!("{m} SNPs"),
            actual: format!(
                "{}, {}, {}",
                admixed.n_snps(),
                ref_a.n_snps(),
                ref_b.n_snps()
            ),
        });
    }
    if options.runmode > 1 {
        return Err(PopmixError::Validation(format!(
            "runmode {} is not supported",
            options.runmode
        )));
    }
    let counts = |g: &PackedGenotypes, s: usize| -> [f64; 3] {
        let mut x = [0.0; 3];
        for i in 0..g.n_samples() {
            let v = g.get(s, i);
            if v <= 2 {
                x[v as usize] += 1.0;
            }
        }
        x
    };
    let freq = |x: [f64; 3]| -> Option<f64> {
        let y2 = 2.0 * (x[0] + x[1] + x[2]);
        (y2 > 0.01).then(|| (2.0 * x[2] + x[1]) / y2)
    };
    let mincount = options.mincount.min(admixed.n_samples()) as f64;
    // weights, and the SNPs kept
    let weights: Vec<Option<f64>> = (0..m)
        .into_par_iter()
        .map(|s| {
            let pa = freq(counts(ref_a, s))?;
            let pb = freq(counts(ref_b, s))?;
            let diff = pa - pb;
            if diff.abs() < 0.001 {
                return None;
            }
            let mean = 0.5 * (pa + pb);
            if mean * (1.0 - mean) < 0.001 {
                return None;
            }
            let w = if options.runmode == 0 {
                diff / mean.sqrt()
            } else {
                diff
            };
            let x = counts(admixed, s);
            let (y1, y2) = (2.0 * x[2] + x[1], 2.0 * (x[0] + x[1] + x[2]));
            if round(y1) < mincount || round(y2 - y1) < mincount || w == 0.0 {
                return None;
            }
            Some(w)
        })
        .collect();
    let used: Vec<usize> = (0..m).filter(|&s| weights[s].is_some()).collect();
    if used.len() < 100 {
        return Err(PopmixError::Validation(format!(
            "too few SNPs: {}",
            used.len()
        )));
    }
    let numbins = (options.maxdis / options.binsize).round() as usize + 5;
    let out_bins = numbins - 5;

    // chromosomes, in order
    let mut chroms: Vec<(String, Vec<usize>)> = Vec::new();
    for &s in &used {
        let label = snps[s].chromosome.to_string();
        match chroms.last_mut() {
            Some((l, v)) if *l == label => v.push(s),
            _ => chroms.push((label, vec![s])),
        }
    }
    let per_chrom: Vec<Vec<CorrSums>> = chroms
        .par_iter()
        .map(|(_, idx)| {
            let mut bins = vec![CorrSums::default(); numbins];
            for (a, &x) in idx.iter().enumerate() {
                for &x2 in &idx[a + 1..] {
                    let dis = snps[x2].position_morgans - snps[x].position_morgans;
                    if dis >= options.maxdis {
                        break;
                    }
                    let bin = (dis / options.binsize) as usize;
                    let wt = weights[x].expect("used") * weights[x2].expect("used");
                    if wt == 0.0 {
                        continue;
                    }
                    if options.chithresh > 0.001 {
                        let ta = table(ref_a, x, x2);
                        if ta.iter().sum::<f64>() < options.minparentcount as f64
                            || lddip(&ta) >= options.chithresh
                        {
                            continue;
                        }
                        let tb = table(ref_b, x, x2);
                        if tb.iter().sum::<f64>() < options.minparentcount as f64
                            || lddip(&tb) >= options.chithresh
                        {
                            continue;
                        }
                    }
                    let xc = table(admixed, x, x2);
                    if xc.iter().sum::<f64>() < 4.0 {
                        continue;
                    }
                    let mut r = [0.0; 3];
                    let mut c = [0.0; 3];
                    for e in 0..3 {
                        for f in 0..3 {
                            r[e] += xc[3 * e + f];
                            c[f] += xc[3 * e + f];
                        }
                    }
                    if round(minor(&r)) < mincount || round(minor(&c)) < mincount {
                        continue;
                    }
                    let yy = zdip(&xc, options.zdipcorr);
                    bins[bin].add(wt, yy, 1.0);
                }
            }
            bins
        })
        .collect();

    let mut total = vec![CorrSums::default(); numbins];
    for c in &per_chrom {
        for (t, b) in total.iter_mut().zip(c) {
            t.merge(b);
        }
    }
    let summarize = |sums: &[CorrSums], first: &dyn Fn(usize, &CorrSums) -> f64| -> RolloffBins {
        let mut out = RolloffBins {
            d: Vec::with_capacity(out_bins),
            zscore: Vec::new(),
            regression: Vec::new(),
            corr: Vec::new(),
            count: Vec::new(),
        };
        for (k, s) in sums.iter().enumerate().take(out_bins) {
            out.d.push((k + 1) as f64 * options.binsize);
            out.zscore.push(first(k, s));
            out.regression.push(s.s12 / (s.s11 + 1e-20));
            out.corr.push(s.corr(true, false).map_or(0.0, |(r, _)| r));
            out.count.push(s.s0);
        }
        out
    };
    let bins = summarize(&total, &|_, s| s.s12 / (s.s11 + 1e-10).sqrt());
    let jackknife = per_chrom
        .iter()
        .map(|c| {
            let rest: Vec<CorrSums> = total.iter().zip(c).map(|(t, b)| t.minus(b)).collect();
            summarize(&rest, &|k, _| c[k].s0)
        })
        .collect();
    Ok(RolloffResult {
        chromosomes: chroms.into_iter().map(|(l, _)| l).collect(),
        bins,
        jackknife,
        n_snps_used: used.len(),
    })
}

/// A fit of C + A exp(-m d) to ROLLOFF correlations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RolloffFit {
    pub affine: f64,
    pub amplitude: f64,
    /// Generations.
    pub decay: f64,
    pub rss: f64,
}

/// Least squares C + A exp(-m d) over the bins with lval <= d <= hval
/// (cM, as rexpfit.r), m within (0, 500]: for each m, C and A in closed
/// form; m by a grid then golden-section search.
pub fn fit_rolloff(d_cm: &[f64], y: &[f64], lval: f64, hval: f64) -> RolloffFit {
    let pts: Vec<(f64, f64)> = d_cm
        .iter()
        .zip(y)
        .filter(|(d, _)| **d >= lval && **d <= hval)
        .map(|(d, y)| (d / 100.0, *y))
        .collect();
    let at = |m: f64| -> RolloffFit {
        let n = pts.len() as f64;
        let (mut se, mut see, mut sy, mut sey) = (0.0, 0.0, 0.0, 0.0);
        for &(d, y) in &pts {
            let e = (-m * d).exp();
            se += e;
            see += e * e;
            sy += y;
            sey += e * y;
        }
        let det = n * see - se * se;
        let a = (n * sey - se * sy) / det;
        let c = (sy - a * se) / n;
        let rss = pts
            .iter()
            .map(|&(d, y)| (y - c - a * (-m * d).exp()).powi(2))
            .sum();
        RolloffFit {
            affine: c,
            amplitude: a,
            decay: m,
            rss,
        }
    };
    // a fine logarithmic grid, then golden-section refinement around the best
    let grid: Vec<f64> = (0..=4000)
        .map(|k| 0.01 * (500.0f64 / 0.01).powf(k as f64 / 4000.0))
        .collect();
    let best = grid
        .iter()
        .enumerate()
        .min_by(|a, b| at(*a.1).rss.total_cmp(&at(*b.1).rss))
        .map(|(k, _)| k)
        .expect("grid");
    let (mut lo, mut hi) = (
        grid[best.saturating_sub(1)],
        grid[(best + 1).min(grid.len() - 1)],
    );
    let phi = (5f64.sqrt() - 1.0) / 2.0;
    let mut x1 = hi - phi * (hi - lo);
    let mut x2 = lo + phi * (hi - lo);
    for _ in 0..200 {
        if at(x1).rss < at(x2).rss {
            hi = x2;
            x2 = x1;
            x1 = hi - phi * (hi - lo);
        } else {
            lo = x1;
            x1 = x2;
            x2 = lo + phi * (hi - lo);
        }
        if hi - lo < 1e-12 * hi {
            break;
        }
    }
    at((lo + hi) / 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lddip_is_zero_for_unlinked_and_large_for_linked_tables() {
        // genotypes at two SNPs independent: product table
        let row = [25.0, 50.0, 25.0];
        let mut indep = [0.0; 9];
        for e in 0..3 {
            for f in 0..3 {
                indep[3 * e + f] = row[e] * row[f] / 100.0;
            }
        }
        assert!(lddip(&indep).abs() < 1e-6);
        // perfectly linked: only matching genotypes
        let linked = [25.0, 0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 25.0];
        assert!(lddip(&linked) > 50.0);
        assert_eq!(lddip(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]), -1.0);
    }

    #[test]
    fn zdip_has_the_sign_of_the_correlation() {
        let pos = [20.0, 5.0, 1.0, 5.0, 30.0, 5.0, 1.0, 5.0, 20.0];
        let neg = [1.0, 5.0, 20.0, 5.0, 30.0, 5.0, 20.0, 5.0, 1.0];
        assert!(zdip(&pos, true) > 0.0 && zdip(&neg, true) < 0.0);
        assert!(zdip(&pos, false) > 0.0 && zdip(&neg, false) < 0.0);
    }

    #[test]
    fn the_fit_recovers_an_exponential() {
        let d: Vec<f64> = (1..=200).map(|k| k as f64 * 0.1).collect();
        let y: Vec<f64> = d
            .iter()
            .map(|x| 0.002 + 0.03 * (-25.0 * x / 100.0).exp())
            .collect();
        let f = fit_rolloff(&d, &y, 0.5, 100.0);
        assert!((f.decay - 25.0).abs() < 1e-6, "{}", f.decay);
        assert!((f.amplitude - 0.03).abs() < 1e-9 && (f.affine - 0.002).abs() < 1e-9);
    }
}
