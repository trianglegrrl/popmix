//! Fitting weighted LD curves, as ALDER: `A exp(-n d) + C` by least squares,
//! the decay rate n found by a grid search then a ternary search, and
//! standard errors by a jackknife over chromosomes.

use nalgebra::{DMatrix, DVector};

/// Bounds of the decay rate (generations).
pub const GEN_MIN: f64 = 2.0;
pub const GEN_MAX: f64 = 500.0;
const GEN_MULTIPLIER: f64 = 1.2;
const GEN_RESOLUTION: f64 = 0.005;

/// One fit of `amp_exp * exp(-decay d) + amp_aff`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DecayFitParams {
    pub decay: f64,
    pub amp_exp: f64,
    pub amp_aff: f64,
    /// Whether the decay is strictly inside its bounds (ALDER's success).
    pub success: bool,
}

/// Least squares A c = b by QR (as LAPACK dgels): the coefficients.
fn least_squares(a: DMatrix<f64>, b: &DVector<f64>) -> Option<DVector<f64>> {
    let qr = a.qr();
    let qtb = qr.q().transpose() * b;
    qr.r().solve_upper_triangular(&qtb)
}

/// Residual sum of squares of the least-squares fit at decay rate `gen`,
/// and the coefficients (exponential amplitude, affine term).
fn lsq(gen: f64, x: &[f64], y: &[f64], affine: Option<f64>) -> (f64, f64, f64) {
    let n = x.len();
    let e: Vec<f64> = x.iter().map(|d| (-gen * d).exp()).collect();
    match affine {
        Some(c) => {
            let b = DVector::from_fn(n, |i, _| y[i] - c);
            let amp =
                least_squares(DMatrix::from_fn(n, 1, |i, _| e[i]), &b).map_or(f64::NAN, |s| s[0]);
            let resid: f64 = (0..n).map(|i| (b[i] - amp * e[i]).powi(2)).sum();
            (resid, amp, c)
        }
        None => {
            let b = DVector::from_fn(n, |i, _| y[i]);
            let a = DMatrix::from_fn(n, 2, |i, j| if j == 0 { 1.0 } else { e[i] });
            let (c, amp) = least_squares(a, &b).map_or((f64::NAN, f64::NAN), |s| (s[0], s[1]));
            let resid: f64 = (0..n).map(|i| (b[i] - c - amp * e[i]).powi(2)).sum();
            (resid, amp, c)
        }
    }
}

/// Fit the curve (x in Morgans, y weighted LD); with `affine` the constant
/// term is fixed to it. ALDER `ExpFit::fit_decay`.
pub fn fit_decay(x: &[f64], y: &[f64], affine: Option<f64>) -> DecayFitParams {
    let mut best_norm = 1e9;
    let mut best_gen = GEN_MIN;
    let mut gen = GEN_MIN;
    while gen < GEN_MAX {
        let norm = lsq(gen, x, y, affine).0;
        if norm < best_norm {
            best_norm = norm;
            best_gen = gen;
        }
        gen *= GEN_MULTIPLIER;
    }
    let mut lo = GEN_MIN.max(best_gen / GEN_MULTIPLIER.powi(2));
    let mut hi = GEN_MAX.min(best_gen * GEN_MULTIPLIER.powi(2));
    while hi - lo > GEN_RESOLUTION {
        let g1 = (2.0 * lo + hi) / 3.0;
        let g2 = (lo + 2.0 * hi) / 3.0;
        if lsq(g1, x, y, affine).0 < lsq(g2, x, y, affine).0 {
            hi = g2;
        } else {
            lo = g1;
        }
    }
    let decay = (lo + hi) / 2.0;
    let (_, amp_exp, amp_aff) = lsq(decay, x, y, affine);
    DecayFitParams {
        decay,
        amp_exp,
        amp_aff,
        success: GEN_MIN + GEN_RESOLUTION < decay && decay < GEN_MAX - GEN_RESOLUTION,
    }
}

/// A quantity's full-data value and its leave-one-chromosome-out replicates.
#[derive(Debug, Clone, PartialEq)]
pub struct Jackknifed {
    pub value: f64,
    pub reps: Vec<f64>,
}

impl Jackknifed {
    /// ALDER's jackknife standard error: sqrt((n - 1) / n * sum (x - mean)^2)
    /// over the replicates; NaN if one is NaN, infinite if one is infinite.
    pub fn se(&self) -> f64 {
        let n = self.reps.len() as f64;
        if self.reps.iter().any(|x| x.is_nan()) {
            return f64::NAN;
        }
        if self.reps.iter().any(|x| x.is_infinite()) {
            return f64::INFINITY;
        }
        let s: f64 = self.reps.iter().sum();
        let s2: f64 = self.reps.iter().map(|x| x * x).sum();
        ((s2 - s * s / n) * (n - 1.0) / n).sqrt()
    }

    pub fn z(&self) -> f64 {
        self.value / self.se()
    }

    /// The difference from `other`, replicate by replicate.
    pub fn minus(&self, other: &Jackknifed) -> Jackknifed {
        Jackknifed {
            value: self.value - other.value,
            reps: self
                .reps
                .iter()
                .zip(&other.reps)
                .map(|(a, b)| a - b)
                .collect(),
        }
    }
}

/// A curve fit from one start distance, with its jackknife.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveFit {
    /// Fit start (Morgans).
    pub mindis: f64,
    pub maxdis: f64,
    pub inter_chrom_affine: bool,
    pub decay: Jackknifed,
    /// amp_exp + amp_aff / 2, as ALDER reports it.
    pub amp_tot: Jackknifed,
    pub amp_exp: Jackknifed,
    pub amp_aff: Jackknifed,
    /// Jackknife replicates that failed to fit (their values infinite).
    pub failures: usize,
}

impl CurveFit {
    /// ALDER's p-value for a curve: from the smaller of the decay and
    /// amplitude z-scores, times a multiple-testing factor.
    pub fn p_value(&self, correction: f64) -> f64 {
        let z = self.decay.z().min(self.amp_exp.z());
        statrs::function::erf::erfc(z / std::f64::consts::SQRT_2) * correction
    }

    pub fn is_significant(&self, correction: f64) -> bool {
        self.p_value(correction) < 0.05
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_exponential_curve_is_recovered() {
        let x: Vec<f64> = (10..1000).map(|b| b as f64 * 0.0005).collect();
        let y: Vec<f64> = x.iter().map(|d| 3e-4 * (-40.0 * d).exp() + 2e-6).collect();
        let free = fit_decay(&x, &y, None);
        assert!((free.decay - 40.0).abs() < 0.01 && free.success);
        assert!((free.amp_exp - 3e-4).abs() < 1e-8 && (free.amp_aff - 2e-6).abs() < 1e-9);
        let fixed = fit_decay(&x, &y, Some(2e-6));
        assert!((fixed.decay - 40.0).abs() < 0.01 && fixed.amp_aff == 2e-6);
    }

    #[test]
    fn a_decay_at_the_bound_is_a_failure() {
        let x: Vec<f64> = (10..100).map(|b| b as f64 * 0.0005).collect();
        let y: Vec<f64> = x.iter().map(|_| 1e-5).collect();
        assert!(!fit_decay(&x, &y, Some(0.0)).success);
    }

    #[test]
    fn jackknife_standard_error() {
        let j = Jackknifed {
            value: 1.0,
            reps: vec![0.9, 1.1, 1.0],
        };
        let want = ((0.01 + 0.01) * 2.0 / 3.0f64).sqrt();
        assert!((j.se() - want).abs() < 1e-12);
        let inf = Jackknifed {
            value: 1.0,
            reps: vec![1.0, f64::INFINITY],
        };
        assert_eq!(inf.se(), f64::INFINITY);
    }
}
