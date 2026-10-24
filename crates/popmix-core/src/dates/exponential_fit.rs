//! Exponential decay fitting for DATES algorithm.
//!
//! Fits the LD decay curve to estimate admixture time, reproducing DATES v753
//! `dates_expfit` (one exponential, affine mode) and `dowtjack`.
//!
//! DATES model: a_i = A × r^(i+1) + K
//!
//! Where:
//! - a_i = correlation at the i-th fitted bin (column 4 of the `.out` file)
//! - A = amplitude, K = affine constant, both fitted by unconstrained least
//!   squares for each r
//! - r = per-bin decay (0 < r < 1)
//!
//! The admixture time in generations is n = -ln(r) / step, where step is the
//! bin width in Morgans. DATES finds r with a seeded random search followed by
//! a GSL simplex; popmix instead minimises the same mean squared error to full
//! precision, which agrees with DATES to its printed precision.

use super::LdBin;
use crate::error::PopmixError;
use crate::Result;

/// Options for exponential fitting.
#[derive(Debug, Clone)]
pub struct FitOptions {
    /// Fit bins from the first one at or beyond this distance (cM).
    /// DATES `lovalfit`. Default: 0.45
    pub min_distance_cm: f64,
    /// Stop at the first bin beyond this distance (cM). Default: 100.0
    pub max_distance_cm: f64,
    /// Use affine model (with constant K). Default: true
    pub use_affine: bool,
    /// Bin size in cM. Default: 0.1
    pub bin_size_cm: f64,
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            min_distance_cm: 0.45,
            max_distance_cm: 100.0, // DATES default: maxdis = 1.0 Morgan = 100 cM
            use_affine: true,
            bin_size_cm: 0.1,
        }
    }
}

/// Result of exponential decay fitting.
#[derive(Debug, Clone)]
pub struct DecayFit {
    /// Amplitude at the first fitted bin (A × r).
    pub amplitude: f64,
    /// Decay rate in generations (n).
    pub decay_rate: f64,
    /// Affine constant (K), 0 if affine model not used.
    pub intercept: f64,
    /// R-squared goodness of fit.
    pub r_squared: f64,
    /// Standard error of decay rate (from jackknife, if run).
    pub decay_rate_se: f64,
    /// Mean squared error of the fit (DATES `scorit`).
    pub rss: f64,
}

/// Least-squares fit of A and K for a fixed r (DATES `scorit`).
struct LinearFit {
    a: f64,
    k: f64,
    mse: f64,
}

fn fit_linear(y: &[f64], r: f64, use_affine: bool) -> LinearFit {
    let n = y.len() as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut x = 1.0;
    for &yi in y {
        x *= r;
        sx += x;
        sy += yi;
        sxx += x * x;
        sxy += x * yi;
    }

    let (a, k) = if use_affine {
        let det = n * sxx - sx * sx;
        if det.abs() < 1e-300 {
            return LinearFit {
                a: 0.0,
                k: 0.0,
                mse: f64::INFINITY,
            };
        }
        ((n * sxy - sx * sy) / det, (sxx * sy - sx * sxy) / det)
    } else if sxx > 1e-300 {
        (sxy / sxx, 0.0)
    } else {
        (0.0, 0.0)
    };

    let mut sse = 0.0;
    let mut x = 1.0;
    for &yi in y {
        x *= r;
        let e = yi - (a * x + k);
        sse += e * e;
    }
    LinearFit { a, k, mse: sse / n }
}

/// Golden-section minimisation of f on [lo, hi].
fn golden_section(f: &impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    let invphi = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut x1 = hi - invphi * (hi - lo);
    let mut x2 = lo + invphi * (hi - lo);
    let mut f1 = f(x1);
    let mut f2 = f(x2);
    for _ in 0..200 {
        if hi - lo <= 1e-15 * hi.abs().max(1e-300) {
            break;
        }
        if f1 <= f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - invphi * (hi - lo);
            f1 = f(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + invphi * (hi - lo);
            f2 = f(x2);
        }
    }
    if f1 <= f2 {
        x1
    } else {
        x2
    }
}

/// Fit exponential decay model to LD bins as DATES `dates_expfit` does.
///
/// Bins are fitted from the first with distance ≥ `min_distance_cm` up to the
/// last with distance ≤ `max_distance_cm`, consecutively.
///
/// # Arguments
/// * `ld_bins` - LD decay bins from compute_ld_decay, in distance order
/// * `options` - Fitting options
///
/// # Returns
/// DecayFit with the fitted parameters
pub fn fit_exponential(ld_bins: &[LdBin], options: Option<FitOptions>) -> Result<DecayFit> {
    let opts = options.unwrap_or_default();

    let y: Vec<f64> = ld_bins
        .iter()
        .skip_while(|b| b.distance_cm < opts.min_distance_cm)
        .take_while(|b| b.distance_cm <= opts.max_distance_cm)
        .map(|b| b.correlation)
        .collect();

    if y.len() < 4 {
        return Err(PopmixError::InsufficientData {
            context: "exponential fit (bins in fitting range)".to_string(),
            required: 4,
            actual: y.len(),
        });
    }
    if let Some(bad) = y.iter().position(|v| !v.is_finite()) {
        return Err(PopmixError::InvalidValue {
            context: "exponential fit".to_string(),
            details: format!("non-finite correlation in fitted bin {bad}"),
        });
    }

    // Search u = -ln(r) on a log grid covering 1e-7 .. 20 per bin, then refine
    // every local minimum on the grid and keep the best.
    let mse = |r: f64| fit_linear(&y, r, opts.use_affine).mse;
    let n_grid = 4000;
    let (u_lo, u_hi) = (1e-7_f64.ln(), 20.0_f64.ln());
    let grid: Vec<f64> = (0..n_grid)
        .map(|k| (-(u_lo + (u_hi - u_lo) * k as f64 / (n_grid - 1) as f64).exp()).exp())
        .collect();
    let values: Vec<f64> = grid.iter().map(|&r| mse(r)).collect();

    let mut best_r = f64::NAN;
    let mut best_mse = f64::INFINITY;
    for k in 0..n_grid {
        let left = if k == 0 { f64::INFINITY } else { values[k - 1] };
        let right = if k + 1 == n_grid {
            f64::INFINITY
        } else {
            values[k + 1]
        };
        if values[k] > left || values[k] > right {
            continue;
        }
        // grid is in decreasing r
        let hi = if k == 0 { grid[0] } else { grid[k - 1] };
        let lo = if k + 1 == n_grid {
            grid[k]
        } else {
            grid[k + 1]
        };
        let r = golden_section(&mse, lo, hi);
        let m = mse(r);
        if m < best_mse {
            best_mse = m;
            best_r = r;
        }
    }

    if !best_r.is_finite() {
        return Err(PopmixError::InvalidValue {
            context: "exponential fit".to_string(),
            details: "no minimum of the fit error found".to_string(),
        });
    }

    let fit = fit_linear(&y, best_r, opts.use_affine);
    let step_morgans = opts.bin_size_cm / 100.0;
    let decay_rate = -best_r.ln() / step_morgans;

    let mean_y = y.iter().sum::<f64>() / y.len() as f64;
    let ss_tot: f64 = y.iter().map(|v| (v - mean_y).powi(2)).sum();
    let r_squared = if ss_tot > 0.0 {
        1.0 - fit.mse * y.len() as f64 / ss_tot
    } else {
        0.0
    };

    Ok(DecayFit {
        amplitude: fit.a * best_r,
        decay_rate,
        intercept: fit.k,
        r_squared,
        decay_rate_se: f64::NAN,
        rss: fit.mse,
    })
}

/// Weighted block jackknife (DATES `dowtjack` / nicklib `weightjack`).
///
/// # Arguments
/// * `full_estimate` - Estimate from all blocks
/// * `estimates` - (weight, estimate with that block removed) per block;
///   blocks with weight < 1e-6 are ignored
///
/// # Returns
/// (bias-corrected jackknife estimate, standard error)
pub fn weighted_jackknife(full_estimate: f64, estimates: &[(f64, f64)]) -> Result<(f64, f64)> {
    let blocks: Vec<(f64, f64)> = estimates
        .iter()
        .copied()
        .filter(|(w, _)| *w >= 1e-6)
        .collect();
    let g = blocks.len();
    if g <= 1 {
        return Err(PopmixError::InsufficientData {
            context: "weighted jackknife (blocks with positive weight)".to_string(),
            required: 2,
            actual: g,
        });
    }

    let total_weight: f64 = blocks.iter().map(|(w, _)| w).sum();
    let jackest = blocks.iter().map(|(_, m)| full_estimate - m).sum::<f64>()
        + blocks.iter().map(|(w, m)| w * m).sum::<f64>() / total_weight;

    let var = blocks
        .iter()
        .map(|(w, m)| {
            let h = total_weight / w;
            let tau = h * full_estimate - (h - 1.0) * m - jackest;
            tau * tau / (h - 1.0)
        })
        .sum::<f64>()
        / g as f64;

    Ok((jackest, if var > 0.0 { var.sqrt() } else { 0.0 }))
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
            let distance_cm = i as f64 * 0.1;
            let bin_idx = (i - 5) as f64;
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

        let opts = FitOptions {
            min_distance_cm: 0.45,
            max_distance_cm: 30.0,
            ..Default::default()
        };

        let result = fit_exponential(&bins, Some(opts)).unwrap();

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

        let result = fit_exponential(&bins, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_exponential_fit_different_decay_rates() {
        for true_n in [1.0, 20.0, 50.0, 100.0, 200.0] {
            let bins = create_synthetic_ld_bins(true_n, 0.001, 0.00001);
            let opts = FitOptions {
                min_distance_cm: 0.45,
                max_distance_cm: 30.0,
                ..Default::default()
            };
            let result = fit_exponential(&bins, Some(opts)).unwrap();

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
