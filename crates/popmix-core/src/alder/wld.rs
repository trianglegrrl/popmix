//! Weighted LD curves by FFT, as ALDER (`Alder::run_chrom`).
//!
//! SNPs on a chromosome are binned by genetic distance from its first SNP;
//! a pair's distance is the difference of its bins. Every sum over SNP pairs
//! at a distance is a cross-correlation of per-bin sums, accumulated in the
//! frequency domain and transformed back once.
//!
//! - Two references: the sample covariance of the test population's
//!   genotypes at the two SNPs times both weights (the references' allele
//!   frequency difference). A pair with one SNP missing data enters with the
//!   missing genotypes replaced by that SNP's mean (pairs with both missing
//!   some data are left out).
//! - One reference: ALDER's unbiased "polyache" estimator of
//!   E[(x - pA_x)(y - pA_y)] products from the test population alone, over
//!   SNPs with complete data.
//!
//! Per chromosome the sums needed for the inter-chromosome affine term are
//! kept too.

use rustfft::{num_complex::Complex, FftPlanner};

use crate::pca::PackedGenotypes;

/// A chromosome's binned sums: weighted LD (sum over pairs) and pair counts
/// per distance bin, and its terms of the affine sum.
#[derive(Debug, Clone)]
pub struct ChromCurve {
    pub wld: Vec<f64>,
    pub count: Vec<f64>,
    pub affine: AffineData,
}

/// Per-chromosome sums for the affine term, as ALDER `AffineData`.
#[derive(Debug, Clone, Default)]
pub struct AffineData {
    /// SNPs with complete data in the test population (ALDER counts those
    /// with no weight too, ALDER-0001).
    pub count: f64,
    pub ws: f64,
    pub ss: f64,
    pub s2: f64,
    pub wg: Vec<f64>,
    pub gg: Vec<f64>,
    pub gs: Vec<f64>,
    /// gigj[i * n + j] for i < j
    pub gigj: Vec<f64>,
}

/// The test population's genotypes at one chromosome's SNPs.
pub struct ChromData<'a> {
    pub geno: &'a PackedGenotypes,
    /// SNP indices (rows of `geno`) of the chromosome, in order.
    pub snps: &'a [usize],
    /// Their bins.
    pub bins: &'a [usize],
    pub weights: &'a [f64],
    /// Per SNP: missing count, genotype sum and sum of squares over the test
    /// samples.
    pub missing: &'a [usize],
    pub sum: &'a [f64],
    pub sum2: &'a [f64],
    pub ignore: &'a [bool],
    /// Count SNPs whose weight is not a number in the affine pair count
    /// (ALDER-0001).
    pub count_ignored_in_affine: bool,
}

/// Frequency-domain accumulator of cross-correlations of real per-bin arrays.
struct Correlator {
    n: usize,
    fft: std::sync::Arc<dyn rustfft::Fft<f64>>,
    ifft: std::sync::Arc<dyn rustfft::Fft<f64>>,
    acc: Vec<Complex<f64>>,
}

impl Correlator {
    fn new(numbins_chrom: usize) -> Self {
        // as ALDER: the smallest power of two >= the bins, doubled
        let mut shift = 0;
        while (1usize << shift) < numbins_chrom {
            shift += 1;
        }
        let n = 1usize << (shift + 1);
        let mut planner = FftPlanner::new();
        Correlator {
            n,
            fft: planner.plan_fft_forward(n),
            ifft: planner.plan_fft_inverse(n),
            acc: vec![Complex::new(0.0, 0.0); n],
        }
    }

    fn transform(&self, x: &[f64]) -> Vec<Complex<f64>> {
        let mut z: Vec<Complex<f64>> = (0..self.n)
            .map(|b| Complex::new(*x.get(b).unwrap_or(&0.0), 0.0))
            .collect();
        self.fft.process(&mut z);
        z
    }

    /// Accumulate scale * conj(F f) * (F g).
    fn add(&mut self, f: &[f64], g: &[f64], scale: f64) {
        let zf = self.transform(f);
        let zg = self.transform(g);
        for ((a, x), y) in self.acc.iter_mut().zip(&zf).zip(&zg) {
            *a += x.conj() * y * scale;
        }
    }

    /// Accumulate scale * |F f|^2.
    fn add_self(&mut self, f: &[f64], scale: f64) {
        let zf = self.transform(f);
        for (a, x) in self.acc.iter_mut().zip(&zf) {
            *a += Complex::new(x.norm_sqr() * scale, 0.0);
        }
    }

    fn scale_all(&mut self, s: f64) {
        self.acc.iter_mut().for_each(|a| *a *= s);
    }

    /// The correlation at lags +d and -d added, for d < numbins.
    fn finish(mut self, numbins: usize, numbins_chrom: usize) -> Vec<f64> {
        self.ifft.process(&mut self.acc);
        let inv = 1.0 / self.n as f64;
        let mut out = vec![0.0; numbins];
        for (d, o) in out.iter_mut().enumerate().take(numbins.min(numbins_chrom)) {
            let back = if d == 0 { 0.0 } else { self.acc[self.n - d].re };
            *o = (self.acc[d].re + back) * inv;
        }
        out
    }
}

/// Weighted LD sums over one chromosome, for distance bins 0..numbins.
pub fn run_chrom(c: &ChromData, n: usize, num_refs: usize, numbins: usize) -> ChromCurve {
    let numbins_chrom = c.bins.last().map_or(0, |b| b + 1);
    let n_snps = c.snps.len();
    let mut affine = AffineData {
        count: (0..n_snps)
            .filter(|&s| c.missing[s] == 0 && (c.count_ignored_in_affine || !c.ignore[s]))
            .count() as f64,
        wg: vec![0.0; n],
        gg: vec![0.0; n],
        gs: vec![0.0; n],
        gigj: vec![0.0; n * n],
        ..Default::default()
    };
    let g = |s: usize, i: usize| c.geno.get(c.snps[s], i);
    let used = |s: usize| !c.ignore[s];
    let mut fx = vec![0.0; numbins_chrom];
    let mut gy = vec![0.0; numbins_chrom];
    let clear = |a: &mut Vec<f64>| a.iter_mut().for_each(|x| *x = 0.0);

    // pair counts
    let mut counts = Correlator::new(numbins_chrom);
    for s in (0..n_snps).filter(|&s| used(s)) {
        if c.missing[s] == 0 {
            fx[c.bins[s]] += 1.0;
            gy[c.bins[s]] += 0.5;
        } else {
            gy[c.bins[s]] += 1.0;
        }
    }
    counts.add(&fx, &gy, 1.0);
    let count = counts.finish(numbins, numbins_chrom);

    let mut acc = Correlator::new(numbins_chrom);
    let nf = n as f64;
    if num_refs == 2 {
        for i in 0..=n {
            clear(&mut fx);
            clear(&mut gy);
            for s in (0..n_snps).filter(|&s| used(s)) {
                let k = c.missing[s];
                let kf = k as f64;
                let complete = k == 0;
                let b = c.bins[s];
                let w = c.weights[s];
                let half = if complete { 2.0 } else { 1.0 };
                if i < n {
                    let gt = g(s, i);
                    let value = if gt > 2 {
                        c.sum[s] / (nf - kf)
                    } else {
                        gt as f64
                    };
                    if complete {
                        affine.wg[i] += gt as f64 * w;
                        fx[b] += gt as f64 * w;
                    }
                    gy[b] += value * w / (half * (nf - kf - 1.0));
                } else {
                    if complete {
                        affine.ws += c.sum[s] * w;
                        fx[b] += c.sum[s] * w;
                    }
                    gy[b] -= c.sum[s] * w / (half * (nf - kf) * (nf - kf - 1.0));
                }
            }
            acc.add(&fx, &gy, 1.0);
        }
    } else {
        let s0 = nf;
        let s0p2 = s0 * (s0 - 1.0);
        let s0p3 = s0p2 * (s0 - 2.0);
        let s0p4 = s0p3 * (s0 - 3.0);
        let fill = |f: &mut Vec<f64>, v: &dyn Fn(usize) -> f64| {
            f.iter_mut().for_each(|x| *x = 0.0);
            for s in (0..n_snps).filter(|&s| used(s)) {
                f[c.bins[s]] += v(s);
            }
        };
        let total = |f: &[f64]| f.iter().sum::<f64>();
        let gf = |s: usize, i: usize| g(s, i) as f64;

        fill(&mut fx, &|s| c.weights[s] * c.sum[s]);
        affine.ws = total(&fx);
        acc.add_self(&fx, -4.0 / s0p2);

        for i in 0..n {
            fill(&mut fx, &|s| c.weights[s] * gf(s, i));
            affine.wg[i] = total(&fx);
            acc.add_self(&fx, 4.0 * (s0p2 + s0) / (s0 * s0p2));
        }

        fill(&mut fx, &|s| c.weights[s] * c.sum[s]);
        fill(&mut gy, &|s| c.sum[s] * c.sum[s] - c.sum2[s]);
        acc.add(&fx, &gy, 4.0 / s0p3);

        for i in 0..n {
            fill(&mut fx, &|s| c.weights[s] * gf(s, i));
            fill(&mut gy, &|s| gf(s, i) * gf(s, i) - gf(s, i) * c.sum[s]);
            acc.add(&fx, &gy, 4.0 * (2.0 * s0p2 + s0p3) / (s0p2 * s0p3));
        }

        fill(&mut fx, &|s| 2.0 * c.sum2[s] - c.sum[s] * c.sum[s]);
        fill(&mut gy, &|s| c.sum[s] * c.sum[s]);
        affine.ss = total(&gy);
        acc.add(&fx, &gy, 1.0 / s0p4);

        fill(&mut fx, &|s| c.sum2[s]);
        affine.s2 = total(&fx);
        acc.add_self(&fx, -1.0 / s0p4);

        for i in 0..n {
            fill(&mut fx, &|s| {
                gf(s, i) * c.sum[s] - 2.0 * gf(s, i) * gf(s, i)
            });
            fill(&mut gy, &|s| gf(s, i) * c.sum[s]);
            affine.gs[i] = total(&gy);
            acc.add(&fx, &gy, (4.0 * s0p3 + s0p4) / (s0p3 * s0p4));
        }

        for i in 0..n {
            for j in i + 1..n {
                fill(&mut fx, &|s| gf(s, i) * gf(s, j));
                affine.gigj[i * n + j] = total(&fx);
                acc.add_self(&fx, -2.0 * (2.0 * s0p3 + s0p4) / (s0p3 * s0p4));
            }
        }

        for i in 0..n {
            fill(&mut fx, &|s| gf(s, i) * gf(s, i));
            affine.gg[i] = total(&fx);
            acc.add_self(&fx, (4.0 * s0p3 + s0p4) / (s0p3 * s0p4));
        }
        // 4 for the estimator, 2 for counting each pair both ways
        acc.scale_all(1.0 / 8.0);
    }
    ChromCurve {
        wld: acc.finish(numbins, numbins_chrom),
        count,
        affine,
    }
}

/// The inter-chromosome affine term over the given chromosomes' sums, and
/// its pair count (ALDER `compute_inter_chrom_affine`).
pub fn inter_chrom_affine(chroms: &[&AffineData], n: usize, num_refs: usize) -> (f64, f64) {
    let nf = n as f64;
    let s0p2 = nf * (nf - 1.0);
    let s0p3 = s0p2 * (nf - 2.0);
    let s0p4 = s0p3 * (nf - 3.0);
    let (mut aff, mut pairs) = (0.0, 0.0);
    for (c1, a) in chroms.iter().enumerate() {
        for (c2, b) in chroms.iter().enumerate() {
            if c1 == c2 {
                continue;
            }
            pairs += a.count * b.count;
            if num_refs == 2 {
                for i in 0..n {
                    aff += a.wg[i] * b.wg[i] / (nf - 1.0);
                }
                aff -= a.ws * b.ws / (nf * (nf - 1.0));
            } else {
                aff += a.ws * b.ws * -4.0 / s0p2;
                aff += a.ws * (b.ss - b.s2) * 4.0 / s0p3;
                aff += (2.0 * a.s2 - a.ss) * b.ss / s0p4;
                aff -= a.s2 * b.s2 / s0p4;
                for i in 0..n {
                    aff += a.wg[i] * b.wg[i] * 4.0 * (s0p2 + nf) / (nf * s0p2);
                    aff +=
                        a.wg[i] * (b.gg[i] - b.gs[i]) * 4.0 * (2.0 * s0p2 + s0p3) / (s0p2 * s0p3);
                    aff +=
                        (a.gs[i] - 2.0 * a.gg[i]) * b.gs[i] * (4.0 * s0p3 + s0p4) / (s0p3 * s0p4);
                    aff += a.gg[i] * b.gg[i] * (4.0 * s0p3 + s0p4) / (s0p3 * s0p4);
                    for j in i + 1..n {
                        aff += a.gigj[i * n + j] * b.gigj[i * n + j] * -2.0 * (2.0 * s0p3 + s0p4)
                            / (s0p3 * s0p4);
                    }
                }
            }
        }
    }
    if num_refs == 1 {
        aff /= 4.0;
    }
    (aff / pairs, pairs)
}
