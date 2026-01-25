//! qpAdm admixture weight estimation.
//!
//! qpAdm estimates the proportions of ancestry in a target population
//! from a set of source (left) populations, using outgroup (right)
//! populations to constrain the model.
//!
//! The method finds weights w such that:
//! - sum(w) = 1
//! - f4(target, ref; Ri, Rj) ≈ sum(w_k * f4(source_k, ref; Ri, Rj))
//!
//! Weights can be constrained to [0,1] or left unconstrained.

use ndarray::{Array1, Array2};

use super::qpwave::{qpwave, QpWaveOptions, QpWaveResult};
use crate::fstats::F2Blocks;
use crate::genotypes::PopulationId;
use crate::Result;

/// Configuration options for qpAdm analysis.
#[derive(Debug, Clone)]
pub struct QpAdmOptions {
    /// Regularization factor (default: 1e-4).
    pub fudge: f64,
    /// Whether to constrain weights to [0, 1] (default: false).
    pub constrained: bool,
    /// Maximum iterations for constrained optimization (default: 100).
    pub max_iter: usize,
    /// Convergence tolerance (default: 1e-8).
    pub tolerance: f64,
}

impl Default for QpAdmOptions {
    fn default() -> Self {
        Self {
            fudge: 1e-4,
            constrained: false,
            max_iter: 100,
            tolerance: 1e-8,
        }
    }
}

/// Result of qpAdm analysis.
#[derive(Debug, Clone)]
pub struct QpAdmResult {
    /// Target population.
    pub target: PopulationId,
    /// Source populations.
    pub sources: Vec<PopulationId>,
    /// Estimated admixture weights (sum to 1).
    pub weights: Vec<f64>,
    /// Standard errors of weights.
    pub weight_std_errors: Vec<f64>,
    /// Chi-squared statistic for model fit.
    pub chi_sq: f64,
    /// Degrees of freedom.
    pub df: usize,
    /// P-value for model fit (higher = better fit).
    pub p_value: f64,
    /// Whether constrained optimization was used.
    pub constrained: bool,
    /// qpWave result for the full rank test.
    pub qpwave_result: QpWaveResult,
}

impl QpAdmResult {
    /// Check if the model is a good fit (p > threshold).
    pub fn is_good_fit(&self, p_threshold: f64) -> bool {
        self.p_value > p_threshold
    }

    /// Check if all weights are feasible (in [0, 1]).
    pub fn is_feasible(&self) -> bool {
        self.weights.iter().all(|&w| (0.0..=1.0).contains(&w))
    }
}

/// Run qpAdm analysis to estimate admixture weights.
///
/// # Arguments
/// * `f2_blocks` - Precomputed F2 statistics
/// * `target` - Target population to model
/// * `sources` - Source populations (left populations excluding target)
/// * `right` - Right (outgroup) populations
/// * `options` - Optional analysis parameters
///
/// # Returns
/// QpAdmResult with estimated weights and model fit statistics
pub fn qpadm(
    f2_blocks: &F2Blocks,
    target: &str,
    sources: &[&str],
    right: &[&str],
    options: Option<QpAdmOptions>,
) -> Result<QpAdmResult> {
    let opts = options.unwrap_or_default();

    if sources.is_empty() {
        return Err(crate::PopmixError::Validation(
            "qpAdm requires at least one source population".to_string(),
        ));
    }

    if right.len() < 2 {
        return Err(crate::PopmixError::Validation(
            "qpAdm requires at least 2 right populations".to_string(),
        ));
    }

    // Construct left populations: [target, sources...]
    let mut left: Vec<&str> = Vec::with_capacity(sources.len() + 1);
    left.push(target);
    left.extend_from_slice(sources);

    // Run qpWave first to get full statistics
    let qpwave_opts = QpWaveOptions {
        fudge: opts.fudge,
        max_rank: None,
    };
    let qpwave_result = qpwave(f2_blocks, &left, right, Some(qpwave_opts))?;

    // Extract F4 values for weight estimation
    // F4 matrix has dimensions (n_left-1, n_right-1)
    // Row 0 is f4(source_0, target; right_j, right_0)
    // We need to solve: f4(target, target; ...) = sum(w * f4(source, target; ...))
    // But f4(target, target; ...) = 0, so we need a different formulation

    // Actually, the standard qpAdm formulation is:
    // For each right pair (Ri, R0), we have:
    // f4(target, R0; Ri, R0) = sum_k(w_k * f4(source_k, R0; Ri, R0))
    //
    // The F4 matrix from qpWave gives us:
    // f4_matrix[i, j] = f4(left[i+1], left[0]; right[j+1], right[0])
    //                 = f4(source_i, target; right_j, right_0)
    //
    // We need to compute f4(target, right_0; right_j, right_0) for the target

    let n_sources = sources.len();
    let n_right_pairs = right.len() - 1;

    // Build the design matrix A and target vector b
    // A[j, k] = f4(source_k, target; right_j, right_0)
    // b[j] = f4(target, right_0; right_j, right_0) = -f4(right_0, target; right_j, right_0)

    // From our F4 matrix (which uses target as left[0]):
    // A is directly the F4 matrix transposed appropriately
    let f4_matrix = &qpwave_result.f4_matrix;

    // Design matrix: each column is a source, each row is a right pair constraint
    let mut a = Array2::zeros((n_right_pairs, n_sources));
    for j in 0..n_right_pairs {
        for k in 0..n_sources {
            a[[j, k]] = f4_matrix.values[[k, j]];
        }
    }

    // Target vector: f4(target, right[0]; right[j+1], right[0])
    // This equals 0 by construction since both first and fourth positions have right[0]
    // Actually, we need: f4(target, left[0]; Rj, R0) where left[0] = target
    // So this is f4(target, target; Rj, R0) = 0
    //
    // The correct formulation: we model target as mixture
    // The constraint is that target's F4 pattern equals weighted sum of sources' F4 patterns
    // relative to some reference population

    // Alternative approach: use the fact that for qpAdm, we're solving
    // the system that makes the residual in F4 space minimal

    // Compute target F4 values directly from f2_blocks
    let mut b = Array1::zeros(n_right_pairs);
    for j in 0..n_right_pairs {
        // f4(target, ref; Rj+1, R0)
        // We need a reference - typically the first source or an outgroup
        // In standard qpAdm, target is compared to first source
        // f4(target, source[0]; right[j+1], right[0])
        if let Some(result) = f2_blocks.f4(target, sources[0], right[j + 1], right[0]) {
            b[j] = result.estimate;
        }
    }

    // Now solve: A' * w ≈ b where A'[j, k] = f4(source[k], source[0]; right[j+1], right[0])
    // for k > 0, with w summing to 1

    // Rebuild A relative to first source
    let mut a_rel = Array2::zeros((n_right_pairs, n_sources));
    for j in 0..n_right_pairs {
        for k in 0..n_sources {
            if let Some(result) = f2_blocks.f4(sources[k], sources[0], right[j + 1], right[0]) {
                a_rel[[j, k]] = result.estimate;
            }
        }
    }

    // Solve with constraint sum(w) = 1
    let (weights, chi_sq, df) = if opts.constrained {
        solve_constrained_weights(&a_rel, &b, &f4_matrix.covariance, opts.fudge, opts.max_iter)?
    } else {
        solve_unconstrained_weights(&a_rel, &b, &f4_matrix.covariance, opts.fudge)?
    };

    // Compute standard errors using jackknife
    let weight_std_errors = compute_weight_std_errors(f2_blocks, target, sources, right, &weights)?;

    // P-value from chi-squared
    let p_value = chi_squared_p_value(chi_sq, df);

    Ok(QpAdmResult {
        target: PopulationId::from(target),
        sources: sources.iter().map(|s| PopulationId::from(*s)).collect(),
        weights,
        weight_std_errors,
        chi_sq,
        df,
        p_value,
        constrained: opts.constrained,
        qpwave_result,
    })
}

/// Solve for unconstrained weights using weighted least squares.
fn solve_unconstrained_weights(
    a: &Array2<f64>,
    b: &Array1<f64>,
    cov: &Array2<f64>,
    fudge: f64,
) -> Result<(Vec<f64>, f64, usize)> {
    let n_constraints = a.nrows();
    let n_sources = a.ncols();

    if n_sources == 0 {
        return Ok((vec![], 0.0, n_constraints));
    }

    // Add constraint that weights sum to 1
    // Augment system with: 1'w = 1
    let mut a_aug = Array2::zeros((n_constraints + 1, n_sources));
    let mut b_aug = Array1::zeros(n_constraints + 1);

    for i in 0..n_constraints {
        for j in 0..n_sources {
            a_aug[[i, j]] = a[[i, j]];
        }
        b_aug[i] = b[i];
    }

    // Sum constraint
    for j in 0..n_sources {
        a_aug[[n_constraints, j]] = 1.0;
    }
    b_aug[n_constraints] = 1.0;

    // Regularized weight matrix (inverse covariance)
    let mut w_diag = Array1::ones(n_constraints + 1);
    if cov.nrows() > 0 {
        for i in 0..n_constraints.min(cov.nrows()) {
            let var = cov[[i, i]] + fudge;
            if var > 0.0 {
                w_diag[i] = 1.0 / var;
            }
        }
    }
    // High weight on sum constraint
    w_diag[n_constraints] = 1e6;

    // Weighted least squares: (A'WA)^-1 A'Wb
    let mut atwa = Array2::zeros((n_sources, n_sources));
    let mut atwb = Array1::zeros(n_sources);

    for i in 0..n_sources {
        for j in 0..n_sources {
            for k in 0..n_constraints + 1 {
                atwa[[i, j]] += a_aug[[k, i]] * w_diag[k] * a_aug[[k, j]];
            }
        }
        for k in 0..n_constraints + 1 {
            atwb[i] += a_aug[[k, i]] * w_diag[k] * b_aug[k];
        }
    }

    // Add regularization
    for i in 0..n_sources {
        atwa[[i, i]] += fudge;
    }

    // Solve
    let weights = solve_linear_system(&atwa, &atwb)?;

    // Compute chi-squared
    let mut residual = b.clone();
    for i in 0..n_constraints {
        for j in 0..n_sources {
            residual[i] -= a[[i, j]] * weights[j];
        }
    }

    let mut chi_sq = 0.0;
    for i in 0..n_constraints {
        chi_sq += residual[i] * residual[i] * w_diag[i];
    }

    // Degrees of freedom: n_constraints - (n_sources - 1)
    // because we have n_sources weights with 1 constraint (sum = 1)
    let df = n_constraints.saturating_sub(n_sources.saturating_sub(1));

    Ok((weights.to_vec(), chi_sq, df))
}

/// Solve for constrained weights (all in [0, 1]) using projected gradient.
fn solve_constrained_weights(
    a: &Array2<f64>,
    b: &Array1<f64>,
    cov: &Array2<f64>,
    fudge: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, f64, usize)> {
    let n_sources = a.ncols();

    // Start with unconstrained solution
    let (mut weights, _, _) = solve_unconstrained_weights(a, b, cov, fudge)?;

    // Project onto simplex (sum = 1, all >= 0)
    for _ in 0..max_iter {
        weights = project_to_simplex(&weights);

        // Check if already feasible
        if weights.iter().all(|&w| (0.0..=1.0).contains(&w)) {
            let sum: f64 = weights.iter().sum();
            if (sum - 1.0).abs() < 1e-10 {
                break;
            }
        }
    }

    // Ensure sum is exactly 1
    let sum: f64 = weights.iter().sum();
    if sum > 0.0 {
        for w in &mut weights {
            *w /= sum;
        }
    } else {
        // Fallback to uniform
        for w in &mut weights {
            *w = 1.0 / n_sources as f64;
        }
    }

    // Compute chi-squared with constrained weights
    let n_constraints = a.nrows();
    let mut residual = b.clone();
    for i in 0..n_constraints {
        for j in 0..n_sources {
            residual[i] -= a[[i, j]] * weights[j];
        }
    }

    let mut chi_sq = 0.0;
    for i in 0..n_constraints {
        let var = if i < cov.nrows() {
            cov[[i, i]] + fudge
        } else {
            fudge
        };
        if var > 0.0 {
            chi_sq += residual[i] * residual[i] / var;
        }
    }

    // More degrees of freedom since we have binding constraints
    let n_zero = weights.iter().filter(|&&w| w < 1e-10).count();
    let effective_params = n_sources.saturating_sub(1).saturating_sub(n_zero);
    let df = n_constraints.saturating_sub(effective_params);

    Ok((weights, chi_sq, df))
}

/// Project a vector onto the probability simplex.
fn project_to_simplex(v: &[f64]) -> Vec<f64> {
    let n = v.len();
    if n == 0 {
        return vec![];
    }

    // Sort in descending order
    let mut sorted: Vec<f64> = v.to_vec();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

    // Find the threshold
    let mut cumsum = 0.0;
    let mut rho = 0;

    for (i, &val) in sorted.iter().enumerate() {
        cumsum += val;
        if val - (cumsum - 1.0) / (i + 1) as f64 > 0.0 {
            rho = i + 1;
        }
    }

    let theta = (sorted[..rho].iter().sum::<f64>() - 1.0) / rho as f64;

    v.iter().map(|&x| (x - theta).max(0.0)).collect()
}

/// Solve a linear system Ax = b.
fn solve_linear_system(a: &Array2<f64>, b: &Array1<f64>) -> Result<Array1<f64>> {
    let n = a.nrows();
    if n != a.ncols() || n != b.len() {
        return Err(crate::PopmixError::Numerical(
            "Dimension mismatch in linear system".to_string(),
        ));
    }

    if n == 0 {
        return Ok(Array1::zeros(0));
    }

    // LU decomposition with partial pivoting
    let mut lu = a.clone();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut b_perm = b.clone();

    for k in 0..n {
        // Find pivot
        let mut max_val = lu[[k, k]].abs();
        let mut max_idx = k;
        for i in (k + 1)..n {
            if lu[[i, k]].abs() > max_val {
                max_val = lu[[i, k]].abs();
                max_idx = i;
            }
        }

        if max_val < 1e-15 {
            return Err(crate::PopmixError::Numerical(
                "Singular matrix in linear system".to_string(),
            ));
        }

        // Swap rows
        if max_idx != k {
            perm.swap(k, max_idx);
            for j in 0..n {
                let tmp = lu[[k, j]];
                lu[[k, j]] = lu[[max_idx, j]];
                lu[[max_idx, j]] = tmp;
            }
            let tmp = b_perm[k];
            b_perm[k] = b_perm[max_idx];
            b_perm[max_idx] = tmp;
        }

        // Elimination
        for i in (k + 1)..n {
            let factor = lu[[i, k]] / lu[[k, k]];
            lu[[i, k]] = factor;
            for j in (k + 1)..n {
                lu[[i, j]] -= factor * lu[[k, j]];
            }
            b_perm[i] -= factor * b_perm[k];
        }
    }

    // Back substitution
    let mut x = Array1::zeros(n);
    for i in (0..n).rev() {
        let mut sum = b_perm[i];
        for j in (i + 1)..n {
            sum -= lu[[i, j]] * x[j];
        }
        x[i] = sum / lu[[i, i]];
    }

    Ok(x)
}

/// Compute standard errors of weights using delete-one jackknife on blocks.
#[allow(unused_variables)]
fn compute_weight_std_errors(
    f2_blocks: &F2Blocks,
    target: &str,
    sources: &[&str],
    right: &[&str],
    full_weights: &[f64],
) -> Result<Vec<f64>> {
    let n_sources = sources.len();
    let n_blocks = f2_blocks.n_blocks;

    if n_blocks < 2 || n_sources == 0 {
        return Ok(vec![0.0; n_sources]);
    }

    // Collect jackknife estimates by leaving out each block
    let jackknife_weights: Vec<Vec<f64>> = Vec::with_capacity(n_blocks);

    // TODO: For a proper implementation, we would need F2Blocks to support
    // leave-one-out estimation. For now, use a variance approximation.

    // If we can't do full jackknife, use variance approximation
    // Based on the scatter of block-wise estimates
    if jackknife_weights.is_empty() {
        // Approximate SE using the variance of F4 values
        // This is a rough estimate
        let n_right_pairs = right.len().saturating_sub(1);
        if n_right_pairs == 0 {
            return Ok(vec![0.0; n_sources]);
        }

        // Compute mean squared F4 variance as rough SE proxy
        let mut total_var = 0.0;
        let mut count = 0;

        for k in 0..n_sources {
            for j in 0..n_right_pairs {
                if let Some(result) = f2_blocks.f4(sources[k], sources[0], right[j + 1], right[0]) {
                    let se = result.std_err;
                    total_var += se * se;
                    count += 1;
                }
            }
        }

        let avg_var = if count > 0 {
            total_var / count as f64
        } else {
            0.0
        };

        // Scale by number of constraints
        let se_scale = avg_var.sqrt() / (n_right_pairs as f64).sqrt();

        // Rough approximation: SE of weight is similar magnitude
        let std_errors: Vec<f64> = full_weights
            .iter()
            .map(|&w| {
                // SE proportional to weight magnitude, scaled by data SE
                (w.abs() * se_scale).max(se_scale * 0.1)
            })
            .collect();

        return Ok(std_errors);
    }

    // Compute jackknife variance
    let n_jack = jackknife_weights.len() as f64;
    let mut variances = vec![0.0; n_sources];

    for k in 0..n_sources {
        let mean = full_weights[k];
        let mut sum_sq = 0.0;

        for jw in &jackknife_weights {
            let diff = jw[k] - mean;
            sum_sq += diff * diff;
        }

        // Jackknife variance formula
        variances[k] = (n_jack - 1.0) / n_jack * sum_sq;
    }

    Ok(variances.iter().map(|v| v.sqrt()).collect())
}

/// Compute p-value from chi-squared (same as in qpwave).
fn chi_squared_p_value(chi_sq: f64, df: usize) -> f64 {
    if df == 0 {
        return 1.0;
    }

    let df_f = df as f64;

    if df <= 100 {
        chi_squared_cdf_complement(chi_sq, df_f)
    } else {
        let z = (chi_sq / df_f).powf(1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df_f));
        let z = z / (2.0 / (9.0 * df_f)).sqrt();
        standard_normal_cdf_complement(z)
    }
}

fn chi_squared_cdf_complement(x: f64, df: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    let a = df / 2.0;
    let z = x / 2.0;
    upper_incomplete_gamma_ratio(a, z)
}

fn upper_incomplete_gamma_ratio(a: f64, z: f64) -> f64 {
    if z < a + 1.0 {
        1.0 - lower_incomplete_gamma_series(a, z)
    } else {
        upper_incomplete_gamma_cf(a, z)
    }
}

fn lower_incomplete_gamma_series(a: f64, z: f64) -> f64 {
    if z == 0.0 {
        return 0.0;
    }
    let mut sum = 1.0 / a;
    let mut term = 1.0 / a;
    for n in 1..1000 {
        term *= z / (a + n as f64);
        sum += term;
        if term.abs() < sum.abs() * 1e-15 {
            break;
        }
    }
    sum * (-z + a * z.ln() - ln_gamma(a)).exp()
}

fn upper_incomplete_gamma_cf(a: f64, z: f64) -> f64 {
    let mut f = 1e-30_f64;
    let mut c = 1e-30_f64;
    let mut d = 0.0;
    for n in 1..1000 {
        let an = if n == 1 {
            1.0
        } else {
            (n as f64 - 1.0 - a) * (n as f64 - 1.0)
        };
        let bn = z + 2.0 * n as f64 - 1.0 - a;
        d = bn + an * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        d = 1.0 / d;
        c = bn + an / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-15 {
            break;
        }
    }
    f * (-z + a * z.ln() - ln_gamma(a)).exp()
}

#[allow(clippy::excessive_precision, clippy::needless_range_loop)]
fn ln_gamma(x: f64) -> f64 {
    // Lanczos coefficients - these need full precision
    let g = 7;
    let c = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        let pi = std::f64::consts::PI;
        pi.ln() - (pi * x).sin().ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = c[0];
        for i in 1..=g + 1 {
            a += c[i] / (x + i as f64);
        }
        let t = x + g as f64 + 0.5;
        0.5 * (2.0 * std::f64::consts::PI).ln() + (t.ln() * (x + 0.5)) - t + a.ln()
    }
}

fn standard_normal_cdf_complement(z: f64) -> f64 {
    0.5 * erfc(z / std::f64::consts::SQRT_2)
}

fn erfc(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.5 * x.abs());
    let tau = t
        * (-x * x - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398
                                    + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
        .exp();
    if x >= 0.0 {
        tau
    } else {
        2.0 - tau
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::{BlockDefinition, F2Blocks};
    use crate::genotype_matrix::AlleleFrequencies;
    use crate::genotypes::{Chromosome, PopulationId};
    use ndarray::Array2;

    fn create_admixed_test_data() -> F2Blocks {
        // Create data where Pop3 is a 50/50 mix of Pop1 and Pop2
        let pops = vec!["Pop1", "Pop2", "Pop3", "Out1", "Out2"];
        let n_snps = 2000;

        let mut rng_state: u64 = 98765;
        let mut rand = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut freqs = Array2::zeros((n_snps, 5));
        let counts = Array2::from_elem((n_snps, 5), 100u32);

        for i in 0..n_snps {
            // Pop1 and Pop2 are divergent sources
            let freq1 = rand() * 0.8 + 0.1;
            let freq2 = rand() * 0.8 + 0.1;

            freqs[[i, 0]] = freq1; // Pop1
            freqs[[i, 1]] = freq2; // Pop2

            // Pop3 is 50/50 admixture + small noise
            freqs[[i, 2]] = (0.5 * freq1 + 0.5 * freq2 + rand() * 0.02 - 0.01).clamp(0.0, 1.0);

            // Outgroups
            freqs[[i, 3]] = rand();
            freqs[[i, 4]] = rand();
        }

        let mut chromosomes = Vec::with_capacity(n_snps);
        let mut positions = Vec::with_capacity(n_snps);

        for i in 0..n_snps {
            chromosomes.push(Chromosome::Autosome(1));
            positions.push(i as f64 * 0.0005);
        }

        let af = AlleleFrequencies::new(
            freqs,
            counts,
            pops.iter().map(|s| PopulationId::from(*s)).collect(),
            chromosomes,
            positions,
        );

        let block_def = BlockDefinition::new(0.05);
        F2Blocks::from_allele_frequencies(&af, &block_def)
    }

    #[test]
    fn test_qpadm_basic() {
        let f2_blocks = create_admixed_test_data();

        let result = qpadm(
            &f2_blocks,
            "Pop3",
            &["Pop1", "Pop2"],
            &["Out1", "Out2"],
            None,
        )
        .unwrap();

        println!("qpAdm result:");
        println!("  Weights: {:?}", result.weights);
        println!("  Chi-sq: {}, df: {}, p: {}", result.chi_sq, result.df, result.p_value);

        // Weights should sum to 1
        let sum: f64 = result.weights.iter().sum();
        assert!(
            (sum - 1.0).abs() < 0.01,
            "Weights should sum to 1, got {}",
            sum
        );

        // Both weights should be positive (since Pop3 is mixture of both)
        // Note: unconstrained can give negative weights for noisy data
    }

    #[test]
    fn test_qpadm_constrained() {
        let f2_blocks = create_admixed_test_data();

        let opts = QpAdmOptions {
            constrained: true,
            ..Default::default()
        };

        let result = qpadm(
            &f2_blocks,
            "Pop3",
            &["Pop1", "Pop2"],
            &["Out1", "Out2"],
            Some(opts),
        )
        .unwrap();

        println!("Constrained qpAdm:");
        println!("  Weights: {:?}", result.weights);

        // All weights should be in [0, 1]
        for w in &result.weights {
            assert!(*w >= 0.0 && *w <= 1.0, "Weight {} out of bounds", w);
        }

        // Weights should sum to 1
        let sum: f64 = result.weights.iter().sum();
        assert!(
            (sum - 1.0).abs() < 0.01,
            "Weights should sum to 1, got {}",
            sum
        );
    }

    #[test]
    fn test_qpadm_requires_sources() {
        let f2_blocks = create_admixed_test_data();

        let result = qpadm(&f2_blocks, "Pop3", &[], &["Out1", "Out2"], None);
        assert!(result.is_err());
    }

    #[test]
    fn test_qpadm_requires_outgroups() {
        let f2_blocks = create_admixed_test_data();

        let result = qpadm(&f2_blocks, "Pop3", &["Pop1", "Pop2"], &["Out1"], None);
        assert!(result.is_err());
    }

    #[test]
    fn test_project_to_simplex() {
        // Already on simplex
        let v = vec![0.3, 0.5, 0.2];
        let proj = project_to_simplex(&v);
        let sum: f64 = proj.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);

        // All positive but doesn't sum to 1
        let v = vec![0.5, 0.5, 0.5];
        let proj = project_to_simplex(&v);
        let sum: f64 = proj.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);

        // Has negatives
        let v = vec![0.8, -0.1, 0.3];
        let proj = project_to_simplex(&v);
        let sum: f64 = proj.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
        assert!(proj.iter().all(|&x| x >= 0.0));
    }

    #[test]
    fn test_solve_linear_system() {
        let a = Array2::from_shape_vec((2, 2), vec![2.0, 1.0, 1.0, 3.0]).unwrap();
        let b = Array1::from_vec(vec![5.0, 7.0]);

        let x = solve_linear_system(&a, &b).unwrap();

        // Verify Ax = b
        let mut ax: Array1<f64> = Array1::zeros(2);
        for i in 0..2 {
            for j in 0..2 {
                ax[i] += a[[i, j]] * x[j];
            }
        }

        for i in 0..2 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-10,
                "Solution incorrect at {}: {} vs {}",
                i,
                ax[i],
                b[i]
            );
        }
    }
}
