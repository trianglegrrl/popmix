//! The exponential fit of the LD decay, and the weighted block jackknife.
//!
//! The fit: correlations y_k of consecutive distance bins, from `lovalfit`
//! on, are fitted by A·r^(k+1) + K (affine) or A·r^(k+1); for a fixed r,
//! A and K are the least-squares solution, so the fit is a one-dimensional
//! minimisation of the mean squared error over 0 < r < 1. The number of
//! generations is −ln(r) / binsize. popmix finds the global minimum to full
//! precision: a log-spaced grid over the decay per bin, each local minimum
//! refined by golden-section search.
//!
//! The jackknife: the weighted delete-one-block jackknife of Busing, Meijer
//! & van der Leeden 1999, over chromosomes weighted by their SNP counts, as
//! DATES uses it.
//!
//! Provenance: written from the clean-room specification
//! `docs/cleanroom/dates.md` (sections 5 and 6), Narasimhan et al. 2019
//! (Science, supplementary methods of DATES), Chintalapati, Patterson &
//! Moorjani 2022 (eLife), Busing et al. 1999 (Statistics and Computing
//! 9:3-8), and DATES v753 run as a black box (`tests/golden/dates/v753/`).

use super::ld_decay::LdBin;
use crate::error::{PopmixError, Result};

/// Which bins to fit and how.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct DecayFitOptions {
    /// First distance fitted (cM): leading bins below it are skipped
    /// (DATES's `lovalfit`, in cM).
    pub min_distance_cm: f64,
    /// Last distance fitted (cM).
    pub max_distance_cm: f64,
    /// Fit an additive constant K (DATES's `afffit`).
    pub use_affine: bool,
    /// Width of a bin (cM); converts the decay per bin to generations.
    pub bin_size_cm: f64,
}

impl Default for DecayFitOptions {
    fn default() -> Self {
        DecayFitOptions {
            min_distance_cm: 0.45,
            max_distance_cm: 100.0,
            use_affine: true,
            bin_size_cm: 0.1,
        }
    }
}

/// The fitted decay A·exp(−n·d) + K.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct DecayFit {
    /// A·r: the fitted excess over K at the first fitted bin.
    pub amplitude: f64,
    /// n, in generations.
    pub decay_rate: f64,
    /// K (0 when not affine).
    pub intercept: f64,
    /// 1 − RSS / total sum of squares (0 if the correlations are constant).
    pub r_squared: f64,
    /// Jackknife standard error of n (NaN until `dates` fills it in).
    pub decay_rate_se: f64,
    /// The mean squared error at the fit.
    pub rss: f64,
}

/// Grid points of the decay per bin u = −ln r.
const GRID_POINTS: usize = 4000;
const U_MIN: f64 = 1e-7;
const U_MAX: f64 = 20.0;
const GOLDEN_ITERATIONS: usize = 200;
const GOLDEN_REL_WIDTH: f64 = 1e-15;

/// The least-squares fit at one value of r.
struct LinearFit {
    a: f64,
    k: f64,
    mse: f64,
}

/// The least-squares A (and K, if `affine`) for x_k = r^(k+1), and the mean
/// squared error; +∞ error when the affine system is singular.
fn solve_at(y: &[f64], r: f64, affine: bool) -> LinearFit {
    let m = y.len() as f64;
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
    let mut x = r;
    for &yk in y {
        sx += x;
        sy += yk;
        sxx += x * x;
        sxy += x * yk;
        x *= r;
    }
    let (a, k) = if affine {
        let det = m * sxx - sx * sx;
        if det.abs() < 1e-300 {
            return LinearFit {
                a: f64::NAN,
                k: f64::NAN,
                mse: f64::INFINITY,
            };
        }
        ((m * sxy - sx * sy) / det, (sxx * sy - sx * sxy) / det)
    } else if sxx <= 1e-300 {
        (0.0, 0.0)
    } else {
        (sxy / sxx, 0.0)
    };
    let mut ss = 0.0;
    let mut x = r;
    for &yk in y {
        let e = yk - a * x - k;
        ss += e * e;
        x *= r;
    }
    LinearFit { a, k, mse: ss / m }
}

/// Mean squared error at r, NaN taken as +∞ so comparisons stay ordered.
fn mse_at(y: &[f64], r: f64, affine: bool) -> f64 {
    let e = solve_at(y, r, affine).mse;
    if e.is_nan() {
        f64::INFINITY
    } else {
        e
    }
}

/// Golden-section search for the minimum of f on [lo, hi]; returns (r, f(r)).
fn golden_section(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> (f64, f64) {
    let inv_phi = (5f64.sqrt() - 1.0) / 2.0;
    let mut c = hi - inv_phi * (hi - lo);
    let mut d = lo + inv_phi * (hi - lo);
    let mut fc = f(c);
    let mut fd = f(d);
    for _ in 0..GOLDEN_ITERATIONS {
        if hi - lo <= GOLDEN_REL_WIDTH * hi.abs() {
            break;
        }
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - inv_phi * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + inv_phi * (hi - lo);
            fd = f(d);
        }
    }
    if fc <= fd {
        (c, fc)
    } else {
        (d, fd)
    }
}

/// The r in (0, 1) minimising the mean squared error of the fit to `y`.
fn minimise(y: &[f64], affine: bool) -> Option<f64> {
    let (ln_lo, ln_hi) = (U_MIN.ln(), U_MAX.ln());
    let step = (ln_hi - ln_lo) / (GRID_POINTS - 1) as f64;
    let grid: Vec<f64> = (0..GRID_POINTS)
        .map(|t| (-(ln_lo + t as f64 * step).exp()).exp())
        .collect();
    let errs: Vec<f64> = grid.iter().map(|&r| mse_at(y, r, affine)).collect();

    let mut best: Option<(f64, f64)> = None;
    for t in 0..GRID_POINTS {
        let left = if t == 0 { f64::INFINITY } else { errs[t - 1] };
        let right = if t + 1 == GRID_POINTS {
            f64::INFINITY
        } else {
            errs[t + 1]
        };
        if !(errs[t] <= left && errs[t] <= right) {
            continue;
        }
        // r decreases along the grid: the bracket runs from the next grid
        // point (smaller r) to the previous one (larger r).
        let lo = grid[(t + 1).min(GRID_POINTS - 1)];
        let hi = grid[t.saturating_sub(1)];
        let (r, e) = golden_section(|r| mse_at(y, r, affine), lo, hi);
        if !r.is_finite() || !e.is_finite() {
            continue;
        }
        if best.is_none_or(|(_, be)| e < be) {
            best = Some((r, e));
        }
    }
    best.map(|(r, _)| r)
}

/// Fit A·r^(k+1) + K to the correlations of the bins in the fitting range.
///
/// Bins are taken in the order given: leading bins with `distance_cm` below
/// `min_distance_cm` are skipped, then bins are taken while `distance_cm` is
/// at most `max_distance_cm`. At least 4 are needed.
pub fn fit_exponential(ld_bins: &[LdBin], options: &DecayFitOptions) -> Result<DecayFit> {
    let opts = options;
    let y: Vec<f64> = ld_bins
        .iter()
        .skip_while(|b| b.distance_cm < opts.min_distance_cm)
        .take_while(|b| b.distance_cm <= opts.max_distance_cm)
        .map(|b| b.correlation)
        .collect();
    if y.len() < 4 {
        return Err(PopmixError::too_few(
            "the exponential fit",
            "bins in the fitting range",
            4,
            y.len(),
        ));
    }
    if let Some(i) = y.iter().position(|v| !v.is_finite()) {
        return Err(PopmixError::numerical(
            "exponential fit",
            format!("non-finite correlation in fitted bin {i}"),
        ));
    }

    let r = minimise(&y, opts.use_affine).ok_or_else(|| {
        PopmixError::numerical("exponential fit", "no minimum of the fit error found")
    })?;
    let fit = solve_at(&y, r, opts.use_affine);

    let m = y.len() as f64;
    let mean = y.iter().sum::<f64>() / m;
    let ss_tot: f64 = y.iter().map(|v| (v - mean) * (v - mean)).sum();
    let r_squared = if ss_tot > 0.0 {
        1.0 - fit.mse * m / ss_tot
    } else {
        0.0
    };
    let step = opts.bin_size_cm / 100.0;
    Ok(DecayFit {
        amplitude: fit.a * r,
        decay_rate: -r.ln() / step,
        intercept: if opts.use_affine { fit.k } else { 0.0 },
        r_squared,
        decay_rate_se: f64::NAN,
        rss: fit.mse,
    })
}

/// Weighted delete-one-block jackknife (Busing et al. 1999).
///
/// `estimates` holds (weight, estimate with that block left out); blocks of
/// weight below 1e-6 are ignored, and at least two must remain. Returns the
/// bias-corrected estimate and its standard error.
pub fn weighted_jackknife(full_estimate: f64, estimates: &[(f64, f64)]) -> Result<(f64, f64)> {
    let blocks: Vec<(f64, f64)> = estimates
        .iter()
        .copied()
        .filter(|&(w, _)| w >= 1e-6)
        .collect();
    let g = blocks.len();
    if g <= 1 {
        return Err(PopmixError::too_few(
            "the weighted jackknife",
            "blocks with positive weight",
            2,
            g,
        ));
    }
    let total: f64 = blocks.iter().map(|&(w, _)| w).sum();
    let bias: f64 = blocks.iter().map(|&(_, t)| full_estimate - t).sum();
    let weighted_mean: f64 = blocks.iter().map(|&(w, t)| w * t).sum::<f64>() / total;
    let estimate = bias + weighted_mean;

    let mut var = 0.0;
    for &(w, t) in &blocks {
        let h = total / w;
        let pseudo = h * full_estimate - (h - 1.0) * t;
        let d = pseudo - estimate;
        var += d * d / (h - 1.0);
    }
    var /= g as f64;
    Ok((estimate, if var > 0.0 { var.sqrt() } else { 0.0 }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_synthetic_ld_bins(true_n: f64, amplitude: f64, intercept: f64) -> Vec<LdBin> {
        // DATES model: a(d_i) = A * r^i + K where r = exp(-n * step)
        let step_morgans = 0.001; // 0.1 cM bins
        let true_r = (-true_n * step_morgans).exp();

        let mut bins = Vec::new();
        // Generate bins from 0.5 to 30 cM
        for i in 5..=300 {
            let distance_cm = f64::from(i) * 0.1;
            let bin_idx = f64::from(i - 5);
            let true_ld = amplitude * true_r.powf(bin_idx) + intercept;

            bins.push(LdBin {
                distance_cm,
                weighted_ld: true_ld,
                z_score: true_ld * 100.0,
                correlation: true_ld,
                regression_coef: true_ld,
                n_pairs: 10000,
                std_err: amplitude * 0.01,
            });
        }
        bins
    }

    #[test]
    fn test_exponential_fit_recovery() {
        let true_n = 50.0;
        let true_a = 0.001;
        let true_k = 0.00001;

        let bins = create_synthetic_ld_bins(true_n, true_a, true_k);

        let opts = DecayFitOptions {
            min_distance_cm: 0.45,
            max_distance_cm: 30.0,
            ..Default::default()
        };

        let result = fit_exponential(&bins, &opts).unwrap();

        assert!(
            (result.decay_rate - true_n).abs() < 1e-6 * true_n,
            "Decay rate {} should equal true value {}",
            result.decay_rate,
            true_n
        );
        assert!((result.amplitude - true_a).abs() < 1e-9);
        assert!((result.intercept - true_k).abs() < 1e-12);
        assert!(result.r_squared > 0.999999);
    }

    #[test]
    fn test_exponential_fit_insufficient_data() {
        let bins = vec![LdBin {
            distance_cm: 1.0,
            weighted_ld: 0.001,
            z_score: 0.1,
            correlation: 0.1,
            regression_coef: 0.1,
            n_pairs: 100,
            std_err: 0.0001,
        }];

        let result = fit_exponential(&bins, &Default::default());
        assert!(result.is_err());
    }

    #[test]
    fn test_exponential_fit_different_decay_rates() {
        for true_n in [1.0, 20.0, 50.0, 100.0, 200.0] {
            let bins = create_synthetic_ld_bins(true_n, 0.001, 0.00001);
            let opts = DecayFitOptions {
                min_distance_cm: 0.45,
                max_distance_cm: 30.0,
                ..Default::default()
            };
            let result = fit_exponential(&bins, &opts).unwrap();

            assert!(
                (result.decay_rate - true_n).abs() < 1e-6 * true_n,
                "For true_n={}, got {}",
                true_n,
                result.decay_rate,
            );
        }
    }

    #[test]
    fn test_weighted_jackknife_equal_weights_is_ordinary_jackknife() {
        // With equal weights the weighted jackknife reduces to the ordinary
        // delete-one jackknife: est = g*full - (g-1)*mean(m),
        // var = (g-1)/g * sum((m - mean(m))^2).
        let full = 10.0;
        let m = [9.0, 10.5, 11.0, 9.5];
        let g = m.len() as f64;
        let estimates: Vec<(f64, f64)> = m.iter().map(|&x| (1.0, x)).collect();
        let (est, se) = weighted_jackknife(full, &estimates).unwrap();

        let mean_m = m.iter().sum::<f64>() / g;
        let expected_est = g * full - (g - 1.0) * mean_m;
        let expected_var = (g - 1.0) / g * m.iter().map(|x| (x - mean_m).powi(2)).sum::<f64>();
        assert!((est - expected_est).abs() < 1e-12);
        assert!((se - expected_var.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_weighted_jackknife_ignores_zero_weight_blocks() {
        let with_empty = [(0.0, 3.0), (100.0, 9.0), (50.0, 11.0), (0.0, 3.0)];
        let without = [(100.0, 9.0), (50.0, 11.0)];
        assert_eq!(
            weighted_jackknife(10.0, &with_empty).unwrap(),
            weighted_jackknife(10.0, &without).unwrap()
        );
        assert!(weighted_jackknife(10.0, &[(1.0, 9.0), (0.0, 1.0)]).is_err());
    }
}
