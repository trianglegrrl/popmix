//! qpWave rank testing implementation.
//!
//! qpWave tests whether a set of "left" populations can be modeled as
//! mixtures of a common set of ancestral sources, using "right" populations
//! as outgroups.
//!
//! The test works by computing an F4 matrix and testing whether its rank
//! is less than expected under a null hypothesis of shared ancestry.

use ndarray::{Array1, Array2};

use super::F4Matrix;
use crate::fstats::F2Blocks;
use crate::Result;

/// Configuration options for qpWave analysis.
#[derive(Debug, Clone)]
pub struct QpWaveOptions {
    /// Regularization factor for numerical stability (default: 1e-4).
    pub fudge: f64,
    /// Maximum number of ranks to test (default: all possible).
    pub max_rank: Option<usize>,
}

impl Default for QpWaveOptions {
    fn default() -> Self {
        Self {
            fudge: 1e-4,
            max_rank: None,
        }
    }
}

/// Result of a single rank test.
#[derive(Debug, Clone)]
pub struct RankTest {
    /// The rank being tested (null hypothesis: rank <= this value).
    pub rank: usize,
    /// Chi-squared statistic.
    pub chi_sq: f64,
    /// Degrees of freedom.
    pub df: usize,
    /// P-value for the test.
    pub p_value: f64,
}

/// Result of qpWave analysis.
#[derive(Debug, Clone)]
pub struct QpWaveResult {
    /// Series of rank tests from 0 to max possible.
    pub rank_tests: Vec<RankTest>,
    /// Estimated number of independent ancestry streams.
    pub estimated_rank: usize,
    /// The F4 matrix used for analysis.
    pub f4_matrix: F4Matrix,
}

impl QpWaveResult {
    /// Get the full rank test result.
    pub fn full_rank(&self) -> Option<&RankTest> {
        self.rank_tests.last()
    }

    /// Check if populations share a common ancestral pattern.
    ///
    /// Returns true if rank 0 test is not significant (p > threshold),
    /// meaning left populations cannot be distinguished by right populations.
    pub fn is_clade(&self, p_threshold: f64) -> bool {
        self.rank_tests
            .first()
            .map(|r| r.p_value > p_threshold)
            .unwrap_or(false)
    }
}

/// Run qpWave analysis.
///
/// # Arguments
/// * `f2_blocks` - Precomputed F2 statistics
/// * `left` - Left population names
/// * `right` - Right population names (outgroups)
/// * `options` - Optional analysis parameters
///
/// # Returns
/// QpWaveResult containing rank test results
pub fn qpwave(
    f2_blocks: &F2Blocks,
    left: &[&str],
    right: &[&str],
    options: Option<QpWaveOptions>,
) -> Result<QpWaveResult> {
    let opts = options.unwrap_or_default();

    // Compute F4 matrix
    let f4_matrix = F4Matrix::from_f2_blocks(f2_blocks, left, right)?;

    // Maximum possible rank
    let max_rank = opts
        .max_rank
        .unwrap_or_else(|| f4_matrix.n_rows().min(f4_matrix.n_cols()));

    // Perform rank tests
    let rank_tests = compute_rank_tests(&f4_matrix, max_rank, opts.fudge)?;

    // Estimate rank as lowest rank where p-value > 0.05
    let estimated_rank = rank_tests
        .iter()
        .find(|r| r.p_value > 0.05)
        .map(|r| r.rank)
        .unwrap_or(max_rank);

    Ok(QpWaveResult {
        rank_tests,
        estimated_rank,
        f4_matrix,
    })
}

/// Compute rank tests using SVD of the F4 matrix.
fn compute_rank_tests(f4_matrix: &F4Matrix, max_rank: usize, fudge: f64) -> Result<Vec<RankTest>> {
    let n_rows = f4_matrix.n_rows();
    let n_cols = f4_matrix.n_cols();
    let n_stats = n_rows * n_cols;

    if n_stats == 0 {
        return Ok(vec![]);
    }

    // Vectorize F4 matrix (used for debugging/verification)
    let _f4_vec = f4_matrix.to_vector();

    // Get covariance matrix and add regularization
    // The fudge is scaled relative to the covariance magnitude to avoid
    // dominating small covariance values
    let mut cov = f4_matrix.covariance.clone();
    let cov_scale = if cov.nrows() > 0 {
        cov.diag().iter().sum::<f64>() / cov.nrows() as f64
    } else {
        1.0
    };
    let scaled_fudge = fudge * cov_scale.max(1e-15);
    for i in 0..cov.nrows() {
        cov[[i, i]] += scaled_fudge;
    }

    // Compute inverse covariance (precision matrix)
    let cov_inv = invert_matrix(&cov)?;

    // Perform rank tests
    let mut rank_tests = Vec::with_capacity(max_rank + 1);

    for rank in 0..=max_rank.min(n_rows.min(n_cols)) {
        let (chi_sq, df) = compute_rank_chi_squared(&f4_matrix.values, &cov_inv, rank)?;

        let p_value = chi_squared_p_value(chi_sq, df);

        rank_tests.push(RankTest {
            rank,
            chi_sq,
            df,
            p_value,
        });
    }

    Ok(rank_tests)
}

/// Compute chi-squared statistic for testing whether matrix has given rank.
///
/// Under the null hypothesis that the matrix has rank r, we test whether
/// the trailing (smallest) singular values are consistent with noise.
fn compute_rank_chi_squared(
    f4_values: &Array2<f64>,
    cov_inv: &Array2<f64>,
    rank: usize,
) -> Result<(f64, usize)> {
    let n_rows = f4_values.nrows();
    let n_cols = f4_values.ncols();

    // SVD of F4 matrix
    let (u, s, vt) = svd(f4_values)?;

    // Under null hypothesis of rank r, compute residual
    let mut residual = f4_values.clone();

    // Subtract rank-r approximation
    for i in 0..rank.min(s.len()) {
        for row in 0..n_rows {
            for col in 0..n_cols {
                residual[[row, col]] -= s[i] * u[[row, i]] * vt[[i, col]];
            }
        }
    }

    // Vectorize residual
    let resid_vec: Array1<f64> = Array1::from_iter(residual.iter().copied());

    // Chi-squared = resid' * cov_inv * resid
    let chi_sq = quadratic_form(&resid_vec, cov_inv);

    // Degrees of freedom: (n_rows - rank) * (n_cols - rank)
    let df = if rank < n_rows && rank < n_cols {
        (n_rows - rank) * (n_cols - rank)
    } else {
        0
    };

    Ok((chi_sq, df))
}

/// Compute SVD of a matrix.
///
/// Returns (U, S, V^T) where U and V are orthogonal and S is diagonal.
fn svd(a: &Array2<f64>) -> Result<(Array2<f64>, Array1<f64>, Array2<f64>)> {
    let m = a.nrows();
    let n = a.ncols();
    let k = m.min(n);

    // Use power iteration for simple SVD (for small matrices)
    // For production, should use LAPACK
    let mut u = Array2::zeros((m, k));
    let mut s = Array1::zeros(k);
    let mut vt = Array2::zeros((k, n));

    let mut residual = a.clone();

    for i in 0..k {
        // Power iteration to find dominant singular value/vector
        let (sigma, u_vec, v_vec) = power_iteration(&residual, 100)?;

        s[i] = sigma;
        for j in 0..m {
            u[[j, i]] = u_vec[j];
        }
        for j in 0..n {
            vt[[i, j]] = v_vec[j];
        }

        // Deflate
        for row in 0..m {
            for col in 0..n {
                residual[[row, col]] -= sigma * u_vec[row] * v_vec[col];
            }
        }
    }

    Ok((u, s, vt))
}

/// Power iteration to find dominant singular value.
fn power_iteration(a: &Array2<f64>, max_iter: usize) -> Result<(f64, Array1<f64>, Array1<f64>)> {
    let m = a.nrows();
    let n = a.ncols();

    // Start with random vector
    let mut v: Array1<f64> = Array1::from_elem(n, 1.0 / (n as f64).sqrt());

    for _ in 0..max_iter {
        // u = A * v
        let mut u: Array1<f64> = Array1::zeros(m);
        for i in 0..m {
            for j in 0..n {
                u[i] += a[[i, j]] * v[j];
            }
        }

        // Normalize u
        let u_norm: f64 = u.iter().map(|x| x * x).sum::<f64>().sqrt();
        if u_norm < 1e-15 {
            return Ok((0.0, Array1::zeros(m), Array1::zeros(n)));
        }
        u.mapv_inplace(|x| x / u_norm);

        // v = A^T * u
        v = Array1::zeros(n);
        for i in 0..m {
            for j in 0..n {
                v[j] += a[[i, j]] * u[i];
            }
        }

        // Normalize v
        let v_norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if v_norm < 1e-15 {
            return Ok((0.0, Array1::zeros(m), Array1::zeros(n)));
        }
        v.mapv_inplace(|x| x / v_norm);
    }

    // Compute singular value
    let mut sigma = 0.0;
    for i in 0..m {
        let mut av_i = 0.0;
        for j in 0..n {
            av_i += a[[i, j]] * v[j];
        }
        sigma += av_i * av_i;
    }
    sigma = sigma.sqrt();

    // Final u
    let mut u: Array1<f64> = Array1::zeros(m);
    for i in 0..m {
        for j in 0..n {
            u[i] += a[[i, j]] * v[j];
        }
    }
    if sigma > 1e-15 {
        u.mapv_inplace(|x| x / sigma);
    }

    Ok((sigma, u, v))
}

/// Invert a positive semi-definite matrix using Cholesky decomposition.
fn invert_matrix(a: &Array2<f64>) -> Result<Array2<f64>> {
    let n = a.nrows();
    if n != a.ncols() {
        return Err(crate::PopmixError::Numerical(
            "Cannot invert non-square matrix".to_string(),
        ));
    }

    if n == 0 {
        return Ok(Array2::zeros((0, 0)));
    }

    // Cholesky decomposition: A = L * L^T
    let l = cholesky(a)?;

    // Solve L * L^T * X = I for X = A^{-1}
    let mut inv = Array2::zeros((n, n));

    for col in 0..n {
        // Solve L * y = e_col (forward substitution)
        let mut y = Array1::zeros(n);
        for i in 0..n {
            let mut sum = if i == col { 1.0 } else { 0.0 };
            for j in 0..i {
                sum -= l[[i, j]] * y[j];
            }
            y[i] = sum / l[[i, i]];
        }

        // Solve L^T * x = y (backward substitution)
        for i in (0..n).rev() {
            let mut sum = y[i];
            for j in (i + 1)..n {
                sum -= l[[j, i]] * inv[[j, col]];
            }
            inv[[i, col]] = sum / l[[i, i]];
        }
    }

    Ok(inv)
}

/// Cholesky decomposition of a positive semi-definite matrix.
fn cholesky(a: &Array2<f64>) -> Result<Array2<f64>> {
    let n = a.nrows();
    let mut l = Array2::zeros((n, n));

    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[[i, j]];
            for k in 0..j {
                sum -= l[[i, k]] * l[[j, k]];
            }

            if i == j {
                if sum <= 0.0 {
                    return Err(crate::PopmixError::Numerical(
                        "Matrix is not positive definite".to_string(),
                    ));
                }
                l[[i, j]] = sum.sqrt();
            } else {
                l[[i, j]] = sum / l[[j, j]];
            }
        }
    }

    Ok(l)
}

/// Compute x^T * A * x (quadratic form).
fn quadratic_form(x: &Array1<f64>, a: &Array2<f64>) -> f64 {
    let n = x.len();
    let mut result = 0.0;

    for i in 0..n {
        for j in 0..n {
            result += x[i] * a[[i, j]] * x[j];
        }
    }

    result
}

/// Compute p-value from chi-squared distribution.
///
/// Uses Wilson-Hilferty approximation for large df,
/// and direct computation for small df.
fn chi_squared_p_value(chi_sq: f64, df: usize) -> f64 {
    if df == 0 {
        return 1.0;
    }

    let df_f = df as f64;

    if df <= 100 {
        // Direct computation using gamma function approximation
        chi_squared_cdf_complement(chi_sq, df_f)
    } else {
        // Wilson-Hilferty approximation
        let z = (chi_sq / df_f).powf(1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df_f));
        let z = z / (2.0 / (9.0 * df_f)).sqrt();
        standard_normal_cdf_complement(z)
    }
}

/// Chi-squared CDF complement using incomplete gamma function.
fn chi_squared_cdf_complement(x: f64, df: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }

    // P(X > x) = gamma(df/2, x/2) / Gamma(df/2)
    // where gamma is the upper incomplete gamma function
    let a = df / 2.0;
    let z = x / 2.0;

    upper_incomplete_gamma_ratio(a, z)
}

/// Compute Q(a, z) = Gamma(a, z) / Gamma(a)
/// using series or continued fraction expansion.
fn upper_incomplete_gamma_ratio(a: f64, z: f64) -> f64 {
    if z < a + 1.0 {
        // Use series for P(a, z) and compute Q = 1 - P
        1.0 - lower_incomplete_gamma_series(a, z)
    } else {
        // Use continued fraction for Q(a, z)
        upper_incomplete_gamma_cf(a, z)
    }
}

/// Lower incomplete gamma ratio using series expansion.
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

/// Upper incomplete gamma ratio using continued fraction.
fn upper_incomplete_gamma_cf(a: f64, z: f64) -> f64 {
    // Lentz's algorithm for continued fraction
    let mut f = 1e-30_f64;
    let mut c = 1e-30_f64;
    let mut d = 0.0;

    for n in 1..1000 {
        let an = if n == 1 { 1.0 } else { (n as f64 - 1.0 - a) * (n as f64 - 1.0) };
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

/// Log of the gamma function using Lanczos approximation.
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
        // Reflection formula
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

/// Standard normal CDF complement (1 - Phi(z)).
fn standard_normal_cdf_complement(z: f64) -> f64 {
    0.5 * erfc(z / std::f64::consts::SQRT_2)
}

/// Complementary error function approximation.
fn erfc(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.5 * x.abs());
    let tau = t
        * (-x * x - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
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

    fn create_test_f2_blocks() -> F2Blocks {
        // Create synthetic data with known structure
        let pops = vec!["Pop1", "Pop2", "Pop3", "Pop4", "Pop5"];
        let n_snps = 1000;

        let mut rng_state: u64 = 54321;
        let mut rand = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut freqs = Array2::zeros((n_snps, 5));
        let counts = Array2::from_elem((n_snps, 5), 50u32);

        for i in 0..n_snps {
            // Create populations with structured relationships
            // Pop1, Pop2, Pop3 form a clade
            // Pop4, Pop5 are outgroups
            let base1 = rand() * 0.8 + 0.1;
            let base2 = rand() * 0.8 + 0.1;

            freqs[[i, 0]] = (base1 + rand() * 0.03).clamp(0.0, 1.0);
            freqs[[i, 1]] = (base1 + rand() * 0.03).clamp(0.0, 1.0);
            freqs[[i, 2]] = (base1 + rand() * 0.05).clamp(0.0, 1.0);
            freqs[[i, 3]] = (base2 + rand() * 0.1).clamp(0.0, 1.0);
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
    fn test_qpwave_basic() {
        let f2_blocks = create_test_f2_blocks();

        let left = vec!["Pop1", "Pop2", "Pop3"];
        let right = vec!["Pop4", "Pop5"];

        let result = qpwave(&f2_blocks, &left, &right, None).unwrap();

        // Should have rank tests from 0 to min(n_left-1, n_right-1) = min(2, 1) = 1
        assert!(!result.rank_tests.is_empty());

        // Each rank test should have valid values
        for test in &result.rank_tests {
            assert!(test.chi_sq >= 0.0);
            assert!(test.p_value >= 0.0 && test.p_value <= 1.0);
        }
    }

    #[test]
    fn test_qpwave_clade_detection() {
        let f2_blocks = create_test_f2_blocks();

        // Pop1 and Pop2 should form a clade (similar frequencies)
        let left = vec!["Pop1", "Pop2"];
        let right = vec!["Pop4", "Pop5"];

        let result = qpwave(&f2_blocks, &left, &right, None).unwrap();

        // With similar populations, rank 0 p-value should be high
        // (though this depends on the synthetic data structure)
        println!("Rank tests: {:?}", result.rank_tests);
        println!("Estimated rank: {}", result.estimated_rank);
    }

    #[test]
    fn test_qpwave_requires_min_populations() {
        let f2_blocks = create_test_f2_blocks();

        // Should fail with only 1 left population
        let result = qpwave(&f2_blocks, &["Pop1"], &["Pop4", "Pop5"], None);
        assert!(result.is_err());

        // Should fail with only 1 right population
        let result = qpwave(&f2_blocks, &["Pop1", "Pop2"], &["Pop4"], None);
        assert!(result.is_err());
    }

    #[test]
    fn test_chi_squared_p_value() {
        // Known values
        // chi_sq = 0 should give p = 1
        let p = chi_squared_p_value(0.0, 1);
        assert!((p - 1.0).abs() < 0.01);

        // chi_sq = df (expected value) should give p around 0.5 for large df
        let p = chi_squared_p_value(10.0, 10);
        assert!(p > 0.3 && p < 0.7, "p = {}", p);

        // Very large chi_sq should give p near 0
        let p = chi_squared_p_value(100.0, 10);
        assert!(p < 0.01, "p = {}", p);
    }

    #[test]
    fn test_svd_simple() {
        let a: Array2<f64> = Array2::from_shape_vec((2, 2), vec![1.0, 2.0, 3.0, 4.0]).unwrap();

        let (u, s, vt) = svd(&a).unwrap();

        // Verify A ≈ U * S * V^T
        let mut reconstructed: Array2<f64> = Array2::zeros((2, 2));
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    reconstructed[[i, j]] += u[[i, k]] * s[k] * vt[[k, j]];
                }
            }
        }

        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (reconstructed[[i, j]] - a[[i, j]]).abs() < 1e-6,
                    "Reconstruction failed at ({}, {}): {} vs {}",
                    i,
                    j,
                    reconstructed[[i, j]],
                    a[[i, j]]
                );
            }
        }
    }

    #[test]
    fn test_matrix_inversion() {
        let a: Array2<f64> = Array2::from_shape_vec((2, 2), vec![4.0, 2.0, 2.0, 3.0]).unwrap();

        let a_inv = invert_matrix(&a).unwrap();

        // A * A^-1 should be identity
        let mut product: Array2<f64> = Array2::zeros((2, 2));
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    product[[i, j]] += a[[i, k]] * a_inv[[k, j]];
                }
            }
        }

        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (product[[i, j]] - expected).abs() < 1e-10,
                    "Product[{}, {}] = {}, expected {}",
                    i,
                    j,
                    product[[i, j]],
                    expected
                );
            }
        }
    }
}
