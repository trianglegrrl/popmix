//! Exponential decay fitting for DATES algorithm.
//!
//! Fits the LD decay curve to estimate admixture time.
//!
//! DATES model: a(d_i) = A × r^i + K
//!
//! Where:
//! - a(d_i) = weighted LD at bin i
//! - A = amplitude
//! - r = per-bin decay rate (0 < r < 1)
//! - K = background LD constant
//! - i = bin index
//!
//! The decay rate in generations is: n = -log(r) / step
//! where step is the bin width in Morgans.

use super::LdBin;
use crate::Result;

/// Options for exponential fitting.
#[derive(Debug, Clone)]
pub struct FitOptions {
    /// Minimum distance for fitting (cM). Default: 0.45
    /// Excludes short-range LD that may not be due to admixture.
    pub min_distance_cm: f64,
    /// Maximum distance for fitting (cM). Default: 50.0
    pub max_distance_cm: f64,
    /// Use affine model (with intercept K). Default: true
    pub use_affine: bool,
    /// Bin size in cM. Default: 0.1
    pub bin_size_cm: f64,
    /// Random seed for initialization. Default: 77
    pub seed: u64,
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            min_distance_cm: 0.45,
            max_distance_cm: 100.0, // DATES default: maxdis = 1.0 Morgan = 100 cM
            use_affine: true,
            bin_size_cm: 0.1,
            seed: 77,
        }
    }
}

/// Result of exponential decay fitting.
#[derive(Debug, Clone)]
pub struct DecayFit {
    /// Amplitude parameter (A).
    pub amplitude: f64,
    /// Decay rate in generations (n).
    pub decay_rate: f64,
    /// Background constant (K), 0 if affine model not used.
    pub intercept: f64,
    /// R-squared goodness of fit.
    pub r_squared: f64,
    /// Standard error of decay rate (from jackknife, if provided).
    pub decay_rate_se: f64,
    /// Residual sum of squares.
    pub rss: f64,
}

/// Simple linear congruential generator for reproducible random numbers.
/// Matches DATES' use of drand48-style PRNG.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        SimpleRng { state: seed }
    }

    /// Returns a random f64 in [0, 1)
    fn next_f64(&mut self) -> f64 {
        // LCG parameters (same as drand48)
        self.state = self.state.wrapping_mul(0x5DEECE66D).wrapping_add(0xB) & 0xFFFFFFFFFFFF;
        (self.state as f64) / (0x1000000000000u64 as f64)
    }
}

/// Fit exponential decay model to LD bins using DATES algorithm.
///
/// Model: a(d_i) = A × r^i + K
///
/// Uses random restarts with annealing (matching DATES), then
/// refines with local optimization. Converts to generations using n = -log(r) / step.
///
/// # Arguments
/// * `ld_bins` - Vector of LD decay bins from compute_ld_decay
/// * `options` - Fitting options
///
/// # Returns
/// DecayFit with the fitted parameters
pub fn fit_exponential(ld_bins: &[LdBin], options: Option<FitOptions>) -> Result<DecayFit> {
    let opts = options.unwrap_or_default();

    // Filter bins by distance
    let valid_bins: Vec<&LdBin> = ld_bins
        .iter()
        .filter(|b| {
            b.distance_cm >= opts.min_distance_cm
                && b.distance_cm <= opts.max_distance_cm
                && b.correlation.is_finite()
                && b.n_pairs > 0
        })
        .collect();

    if valid_bins.len() < 3 {
        return Err(crate::error::PopmixError::InsufficientData {
            context: "exponential fitting".to_string(),
            required: 3,
            actual: valid_bins.len(),
        });
    }

    // Compute bin indices relative to the first valid bin
    // DATES uses bin indices, not distances
    let first_dist = valid_bins[0].distance_cm;
    let bin_indices: Vec<f64> = valid_bins
        .iter()
        .map(|b| (b.distance_cm - first_dist) / opts.bin_size_cm)
        .collect();
    // DATES fits correlation (column 4 in output), not weighted_ld (column 2)
    let ld_values: Vec<f64> = valid_bins.iter().map(|b| b.correlation).collect();

    // Random restart optimization (matching DATES fitexp function)
    // DATES uses inititer = 100 * 2^numexp = 200 for 1 exponential
    let init_iter = 200;
    let out_iter = 10;
    let inn_iter = init_iter / out_iter + out_iter;

    let mut rng = SimpleRng::new(opts.seed);
    let mut best_r = 0.5;
    let mut best_rss = f64::INFINITY;
    let mut mul = 1.0;

    for iter in 1..=init_iter {
        // Generate random r in [0, 1)
        let rand_r = rng.next_f64();

        // Annealing: blend random value with best so far
        let try_r = if best_rss < f64::INFINITY {
            (1.0 - mul) * best_r + mul * rand_r
        } else {
            rand_r
        };

        // Clamp to valid range
        let try_r = try_r.clamp(0.001, 0.9999);

        let (a, _k, rss) = fit_linear_params_dates(&bin_indices, &ld_values, try_r, opts.use_affine);

        if rss < best_rss && a >= 0.0 {
            best_rss = rss;
            best_r = try_r;
        }

        // Reduce search radius periodically (annealing)
        if iter % inn_iter == 0 {
            mul *= 0.5;
        }
    }

    // Local refinement with fine grid search around best value
    let search_min = (best_r - 0.02).max(0.001);
    let search_max = (best_r + 0.02).min(0.9999);
    let fine_step = 0.0001;

    let mut r = search_min;
    while r <= search_max {
        let (a, _k, rss) = fit_linear_params_dates(&bin_indices, &ld_values, r, opts.use_affine);
        if rss < best_rss && a >= 0.0 {
            best_rss = rss;
            best_r = r;
        }
        r += fine_step;
    }

    // Get final fit parameters
    let (best_a, best_k, _) = fit_linear_params_dates(&bin_indices, &ld_values, best_r, opts.use_affine);

    // Convert r to generations
    // n = -log(r) / step, where step is bin_size in Morgans
    let step_morgans = opts.bin_size_cm / 100.0;
    let decay_rate = if best_r > 0.0 && best_r < 1.0 {
        -best_r.ln() / step_morgans
    } else {
        f64::NAN
    };

    // Compute R-squared
    let mean_y: f64 = ld_values.iter().sum::<f64>() / ld_values.len() as f64;
    let ss_tot: f64 = ld_values.iter().map(|y| (y - mean_y).powi(2)).sum();

    let unweighted_rss: f64 = bin_indices
        .iter()
        .zip(ld_values.iter())
        .map(|(i, y)| {
            let pred = best_a * best_r.powf(*i) + best_k;
            (y - pred).powi(2)
        })
        .sum();

    let r_squared = if ss_tot > 0.0 {
        1.0 - unweighted_rss / ss_tot
    } else {
        0.0
    };

    Ok(DecayFit {
        amplitude: best_a,
        decay_rate,
        intercept: best_k,
        r_squared,
        decay_rate_se: f64::NAN,
        rss: best_rss,
    })
}

/// Fit linear parameters (A, K) for a given per-bin decay rate r.
///
/// For fixed r, the model is linear in A and K:
/// y = A × r^i + K
///
/// Returns (A, K, RSS)
fn fit_linear_params_dates(
    bin_indices: &[f64],
    ld_values: &[f64],
    r: f64,
    use_affine: bool,
) -> (f64, f64, f64) {
    let n_points = bin_indices.len();

    if use_affine {
        // Least squares for y = A*r^i + K
        // Design matrix: [r^i, 1]
        let mut sum_1 = 0.0;
        let mut sum_x = 0.0; // x = r^i
        let mut sum_y = 0.0;
        let mut sum_xx = 0.0;
        let mut sum_xy = 0.0;

        for j in 0..n_points {
            let x = r.powf(bin_indices[j]);
            let y = ld_values[j];

            sum_1 += 1.0;
            sum_x += x;
            sum_y += y;
            sum_xx += x * x;
            sum_xy += x * y;
        }

        // Solve normal equations
        let det = sum_1 * sum_xx - sum_x * sum_x;
        if det.abs() < 1e-15 {
            return (0.0, 0.0, f64::INFINITY);
        }

        let a = (sum_1 * sum_xy - sum_x * sum_y) / det;
        let k = (sum_xx * sum_y - sum_x * sum_xy) / det;

        // Compute RSS (mean squared error, like DATES)
        let rss: f64 = bin_indices
            .iter()
            .zip(ld_values.iter())
            .map(|(i, y)| {
                let pred = a * r.powf(*i) + k;
                (y - pred).powi(2)
            })
            .sum::<f64>()
            / n_points as f64;

        (a, k, rss)
    } else {
        // Simple exponential without intercept: y = A*r^i
        let mut sum_xx = 0.0;
        let mut sum_xy = 0.0;

        for j in 0..n_points {
            let x = r.powf(bin_indices[j]);
            let y = ld_values[j];

            sum_xx += x * x;
            sum_xy += x * y;
        }

        let a = if sum_xx > 1e-15 { sum_xy / sum_xx } else { 0.0 };

        // Compute RSS (mean squared error)
        let rss: f64 = bin_indices
            .iter()
            .zip(ld_values.iter())
            .map(|(i, y)| {
                let pred = a * r.powf(*i);
                (y - pred).powi(2)
            })
            .sum::<f64>()
            / n_points as f64;

        (a, 0.0, rss)
    }
}

/// Compute jackknife standard error for decay rate.
///
/// Drops each chromosome and refits to estimate variance.
///
/// # Arguments
/// * `ld_bins_by_chrom` - LD bins grouped by chromosome
/// * `options` - Fitting options
///
/// # Returns
/// (mean_estimate, standard_error)
pub fn jackknife_decay_rate(
    ld_bins_by_chrom: &[Vec<LdBin>],
    options: Option<FitOptions>,
) -> Result<(f64, f64)> {
    let opts = options.unwrap_or_default();
    let n_chroms = ld_bins_by_chrom.len();

    if n_chroms < 2 {
        return Err(crate::error::PopmixError::InsufficientData {
            context: "jackknife (need at least 2 chromosomes)".to_string(),
            required: 2,
            actual: n_chroms,
        });
    }

    // Fit with all data
    let all_bins: Vec<LdBin> = ld_bins_by_chrom.iter().flatten().cloned().collect();
    let full_fit = fit_exponential(&all_bins, Some(opts.clone()))?;

    // Leave-one-out estimates
    let mut jackknife_estimates = Vec::with_capacity(n_chroms);

    for skip_idx in 0..n_chroms {
        let subset: Vec<LdBin> = ld_bins_by_chrom
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != skip_idx)
            .flat_map(|(_, bins)| bins.iter().cloned())
            .collect();

        if let Ok(fit) = fit_exponential(&subset, Some(opts.clone())) {
            jackknife_estimates.push(fit.decay_rate);
        }
    }

    if jackknife_estimates.len() < 2 {
        return Ok((full_fit.decay_rate, f64::NAN));
    }

    // Compute jackknife mean and standard error
    let n = jackknife_estimates.len() as f64;
    let mean: f64 = jackknife_estimates.iter().sum::<f64>() / n;

    let variance: f64 = jackknife_estimates
        .iter()
        .map(|x| (x - mean).powi(2))
        .sum::<f64>()
        * (n - 1.0)
        / n;

    let se = variance.sqrt();

    Ok((mean, se))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_synthetic_ld_bins(true_n: f64, amplitude: f64, intercept: f64) -> Vec<LdBin> {
        // Generate synthetic LD decay data
        // DATES model: a(d_i) = A * r^i + K where r = exp(-n * step)
        let step_morgans = 0.001; // 0.1 cM bins
        let true_r = (-true_n * step_morgans).exp();

        let mut bins = Vec::new();
        // Generate bins from 0.5 to 30 cM
        for i in 5..=300 {
            let distance_cm = i as f64 * 0.1;
            let bin_idx = (distance_cm - 0.5) / 0.1;
            let true_ld = amplitude * true_r.powf(bin_idx) + intercept;

            bins.push(LdBin {
                distance_cm,
                weighted_ld: true_ld,
                z_score: true_ld * 100.0,
                correlation: true_ld,
                n_pairs: 10000,
                std_err: amplitude * 0.01,
            });
        }
        bins
    }

    #[test]
    fn test_exponential_fit_recovery() {
        // Generate data with known parameters
        let true_n = 50.0;
        let true_a = 0.001;
        let true_k = 0.00001;

        let bins = create_synthetic_ld_bins(true_n, true_a, true_k);

        let opts = FitOptions {
            min_distance_cm: 0.5,
            max_distance_cm: 30.0,
            ..Default::default()
        };

        let result = fit_exponential(&bins, Some(opts)).unwrap();

        // Should recover approximately the true parameters
        assert!(
            (result.decay_rate - true_n).abs() < 10.0,
            "Decay rate {} should be close to true value {}",
            result.decay_rate,
            true_n
        );
        assert!(
            result.r_squared > 0.95,
            "R-squared {} should be high for clean data",
            result.r_squared
        );
    }

    #[test]
    fn test_exponential_fit_insufficient_data() {
        let bins = vec![LdBin {
            distance_cm: 1.0,
            weighted_ld: 0.001,
            z_score: 0.1,
            correlation: 0.1,
            n_pairs: 100,
            std_err: 0.0001,
        }];

        let result = fit_exponential(&bins, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_exponential_fit_different_decay_rates() {
        // Test recovery across different decay rates
        for true_n in [20.0, 50.0, 100.0, 200.0] {
            let bins = create_synthetic_ld_bins(true_n, 0.001, 0.00001);
            let opts = FitOptions {
                min_distance_cm: 0.5,
                max_distance_cm: 30.0,
                ..Default::default()
            };
            let result = fit_exponential(&bins, Some(opts)).unwrap();

            assert!(
                (result.decay_rate - true_n).abs() < 15.0,
                "For true_n={}, got {} (diff={})",
                true_n,
                result.decay_rate,
                (result.decay_rate - true_n).abs()
            );
        }
    }
}
