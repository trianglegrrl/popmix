//! Graph fitting: optimize edge lengths to match observed f-statistics.
//!
//! The fitting process minimizes the weighted residual between predicted
//! and observed f-statistics (primarily f2 values).
//!
//! Uses Levenberg-Marquardt algorithm for nonlinear least squares optimization.

use std::collections::HashMap;

use ndarray::{Array1, Array2};

use crate::fstats::F2Blocks;
use crate::{PopmixError, Result};

use super::AdmixGraph;

/// Options for graph fitting.
#[derive(Debug, Clone)]
pub struct FitOptions {
    /// Maximum iterations for optimization.
    pub max_iter: usize,
    /// Convergence tolerance for parameter changes.
    pub tolerance: f64,
    /// Tolerance for cost function improvement.
    pub cost_tolerance: f64,
    /// Minimum edge length (for numerical stability).
    pub min_edge_length: f64,
    /// Maximum edge length.
    pub max_edge_length: f64,
    /// Whether to use weighted residuals (by standard error).
    pub weighted: bool,
    /// Initial Levenberg-Marquardt damping parameter.
    pub initial_lambda: f64,
    /// Factor to increase lambda on bad step.
    pub lambda_increase: f64,
    /// Factor to decrease lambda on good step.
    pub lambda_decrease: f64,
}

impl Default for FitOptions {
    fn default() -> Self {
        Self {
            max_iter: 1000,
            tolerance: 1e-8,
            cost_tolerance: 1e-10,
            min_edge_length: 1e-10,
            max_edge_length: 1.0,
            weighted: true,
            initial_lambda: 0.001,
            lambda_increase: 10.0,
            lambda_decrease: 0.1,
        }
    }
}

/// Result of graph fitting.
#[derive(Debug, Clone)]
pub struct GraphFitResult {
    /// The fitted graph with optimized edge lengths.
    pub graph: AdmixGraph,
    /// Total score (weighted sum of squared residuals).
    pub score: f64,
    /// Root mean square error.
    pub rmse: f64,
    /// Number of f-statistic constraints used.
    pub n_constraints: usize,
    /// Individual residuals for each population pair.
    pub residuals: HashMap<(String, String), f64>,
    /// Number of iterations used.
    pub iterations: usize,
    /// Whether optimization converged.
    pub converged: bool,
}

impl GraphFitResult {
    /// Check if the fit is acceptable (low score).
    pub fn is_good_fit(&self, threshold: f64) -> bool {
        self.score < threshold
    }

    /// Get worst-fitting population pair.
    pub fn worst_pair(&self) -> Option<((String, String), f64)> {
        self.residuals
            .iter()
            .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .map(|(k, v)| (k.clone(), *v))
    }
}

/// Fit an admixture graph to observed f-statistics.
///
/// Optimizes edge lengths to minimize the discrepancy between
/// predicted f2 values (from graph structure) and observed f2 values.
///
/// # Arguments
/// * `graph` - The admixture graph to fit (topology is fixed)
/// * `f2_blocks` - Observed f2 statistics
/// * `options` - Fitting options
///
/// # Returns
/// GraphFitResult containing the fitted graph and fit quality metrics
pub fn fit_graph(
    graph: &AdmixGraph,
    f2_blocks: &F2Blocks,
    options: Option<FitOptions>,
) -> Result<GraphFitResult> {
    let opts = options.unwrap_or_default();

    // Get leaf populations that are present in both graph and data
    let graph_leaves = graph.leaf_names();
    let data_pops: Vec<String> = f2_blocks.populations.iter().map(|p| p.to_string()).collect();

    let common_pops: Vec<String> = graph_leaves
        .iter()
        .filter(|p| data_pops.contains(p))
        .cloned()
        .collect();

    if common_pops.len() < 2 {
        return Err(PopmixError::Validation(
            "Need at least 2 populations present in both graph and data".to_string(),
        ));
    }

    // Build list of population pairs to fit
    let mut pairs: Vec<(String, String)> = Vec::new();
    for i in 0..common_pops.len() {
        for j in (i + 1)..common_pops.len() {
            pairs.push((common_pops[i].clone(), common_pops[j].clone()));
        }
    }

    // Get observed f2 values and weights
    let mut observed: Vec<f64> = Vec::with_capacity(pairs.len());
    let mut weights: Vec<f64> = Vec::with_capacity(pairs.len());

    for (pop_a, pop_b) in &pairs {
        match f2_blocks.f2(pop_a, pop_b) {
            Some((f2, std_err, _)) => {
                observed.push(f2);
                weights.push(if opts.weighted && std_err > 0.0 {
                    1.0 / (std_err * std_err)
                } else {
                    1.0
                });
            }
            None => {
                return Err(PopmixError::Validation(format!(
                    "No f2 data for pair ({}, {})",
                    pop_a, pop_b
                )));
            }
        }
    }

    // Clone graph for optimization
    let mut fitted_graph = graph.clone();

    // Levenberg-Marquardt optimization
    let initial_weights = fitted_graph.get_edge_weights();
    let n_params = initial_weights.len();
    let n_residuals = pairs.len();

    if n_params == 0 {
        return Err(PopmixError::Validation(
            "Graph has no edges to optimize".to_string(),
        ));
    }

    let mut current_weights = initial_weights;
    let mut best_weights = current_weights.clone();

    // Compute initial residuals and cost
    let (mut residual_vec, mut current_cost) =
        compute_residuals(&fitted_graph, &pairs, &observed, &weights)?;
    let mut best_cost = current_cost;

    let mut converged = false;
    let mut iterations = 0;
    let mut lambda = opts.initial_lambda;

    for iter in 0..opts.max_iter {
        iterations = iter + 1;

        // Compute Jacobian matrix numerically (n_residuals x n_params)
        let jacobian = compute_jacobian(
            &mut fitted_graph,
            &mut current_weights,
            &pairs,
            &weights,
            n_residuals,
            n_params,
        )?;

        // Compute J^T * J and J^T * r
        let jtj = jacobian.t().dot(&jacobian);
        let jtr = jacobian.t().dot(&residual_vec);

        // Levenberg-Marquardt: solve (J^T J + lambda * diag(J^T J)) * delta = -J^T * r
        let mut a = jtj.clone();
        for i in 0..n_params {
            a[[i, i]] += lambda * jtj[[i, i]].max(1e-10);
        }

        // Solve linear system using Cholesky or fallback to simple iteration
        let delta = solve_linear_system(&a, &(-&jtr))?;

        // Compute candidate weights
        let mut candidate_weights = current_weights.clone();
        let mut max_change = 0.0_f64;
        for i in 0..n_params {
            let change = delta[i];
            candidate_weights[i] += change;
            // Clamp to valid range
            candidate_weights[i] = candidate_weights[i]
                .max(opts.min_edge_length)
                .min(opts.max_edge_length);
            max_change = max_change.max(change.abs());
        }

        // Evaluate candidate
        fitted_graph.set_edge_weights(&candidate_weights);
        let (new_residual_vec, new_cost) =
            compute_residuals(&fitted_graph, &pairs, &observed, &weights)?;

        // Accept or reject step
        if new_cost < current_cost {
            // Accept step, decrease lambda
            current_weights = candidate_weights;
            current_cost = new_cost;
            residual_vec = new_residual_vec;
            lambda *= opts.lambda_decrease;
            lambda = lambda.max(1e-10);

            if new_cost < best_cost {
                best_cost = new_cost;
                best_weights = current_weights.clone();
            }

            // Check convergence
            let cost_improvement = (current_cost - new_cost).abs() / current_cost.max(1e-15);
            if cost_improvement < opts.cost_tolerance && max_change < opts.tolerance {
                converged = true;
                break;
            }
        } else {
            // Reject step, increase lambda
            lambda *= opts.lambda_increase;
            fitted_graph.set_edge_weights(&current_weights);

            // If lambda gets too large, we're stuck
            if lambda > 1e10 {
                break;
            }
        }

        // Check for parameter convergence
        if max_change < opts.tolerance {
            converged = true;
            break;
        }
    }

    // Set best weights
    fitted_graph.set_edge_weights(&best_weights);

    // Compute final residuals
    let mut residuals = HashMap::new();
    for (pair, obs) in pairs.iter().zip(observed.iter()) {
        let pred = fitted_graph.expected_f2(&pair.0, &pair.1)?;
        residuals.insert(pair.clone(), obs - pred);
    }

    // Compute RMSE
    let rmse = (best_cost / pairs.len() as f64).sqrt();

    Ok(GraphFitResult {
        graph: fitted_graph,
        score: best_cost,
        rmse,
        n_constraints: pairs.len(),
        residuals,
        iterations,
        converged,
    })
}

/// Compute weighted sum of squared residuals.
#[allow(dead_code)]
fn compute_score(
    graph: &AdmixGraph,
    pairs: &[(String, String)],
    observed: &[f64],
    weights: &[f64],
) -> Result<f64> {
    let mut score = 0.0;

    for (i, (pop_a, pop_b)) in pairs.iter().enumerate() {
        let predicted = graph.expected_f2(pop_a, pop_b)?;
        let residual = observed[i] - predicted;
        score += weights[i] * residual * residual;
    }

    Ok(score)
}

/// Compute weighted residual vector and total cost.
fn compute_residuals(
    graph: &AdmixGraph,
    pairs: &[(String, String)],
    observed: &[f64],
    weights: &[f64],
) -> Result<(Array1<f64>, f64)> {
    let n = pairs.len();
    let mut residuals = Array1::zeros(n);
    let mut cost = 0.0;

    for (i, (pop_a, pop_b)) in pairs.iter().enumerate() {
        let predicted = graph.expected_f2(pop_a, pop_b)?;
        let residual = (observed[i] - predicted) * weights[i].sqrt();
        residuals[i] = residual;
        cost += residual * residual;
    }

    Ok((residuals, cost))
}

/// Compute Jacobian matrix numerically using central differences.
fn compute_jacobian(
    graph: &mut AdmixGraph,
    weights_vec: &mut [f64],
    pairs: &[(String, String)],
    fit_weights: &[f64],
    n_residuals: usize,
    n_params: usize,
) -> Result<Array2<f64>> {
    let mut jacobian = Array2::zeros((n_residuals, n_params));
    let eps = 1e-7;

    for j in 0..n_params {
        let original = weights_vec[j];

        // Forward
        weights_vec[j] = original + eps;
        graph.set_edge_weights(weights_vec);
        let mut residuals_plus = Vec::with_capacity(n_residuals);
        for (i, (pop_a, pop_b)) in pairs.iter().enumerate() {
            let predicted = graph.expected_f2(pop_a, pop_b)?;
            residuals_plus.push((pairs[i].0.clone(), predicted));
        }

        // Backward
        weights_vec[j] = original - eps;
        graph.set_edge_weights(weights_vec);
        let mut residuals_minus = Vec::with_capacity(n_residuals);
        for (pop_a, pop_b) in pairs.iter() {
            let predicted = graph.expected_f2(pop_a, pop_b)?;
            residuals_minus.push(predicted);
        }

        // Central difference for each residual
        for i in 0..n_residuals {
            // d(residual)/d(param) = d(observed - predicted)/d(param) = -d(predicted)/d(param)
            let dpred = (residuals_plus[i].1 - residuals_minus[i]) / (2.0 * eps);
            // Weighted residual: r_i = sqrt(w_i) * (obs_i - pred_i)
            // d(r_i)/d(param) = -sqrt(w_i) * d(pred_i)/d(param)
            jacobian[[i, j]] = -fit_weights[i].sqrt() * dpred;
        }

        // Restore
        weights_vec[j] = original;
    }

    graph.set_edge_weights(weights_vec);
    Ok(jacobian)
}

/// Solve linear system A * x = b using Cholesky decomposition.
/// Falls back to regularized solve if Cholesky fails.
fn solve_linear_system(a: &Array2<f64>, b: &Array1<f64>) -> Result<Array1<f64>> {
    let n = a.nrows();

    // Try Cholesky decomposition (A should be positive semi-definite)
    // Simple implementation without external dependencies
    let mut l = Array2::zeros((n, n));

    // Cholesky: A = L * L^T
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[[i, j]];
            for k in 0..j {
                sum -= l[[i, k]] * l[[j, k]];
            }

            if i == j {
                if sum <= 0.0 {
                    // Not positive definite, add regularization and retry
                    let mut a_reg = a.clone();
                    let diag_mean = (0..n).map(|k| a[[k, k]]).sum::<f64>() / n as f64;
                    for k in 0..n {
                        a_reg[[k, k]] += diag_mean * 0.001;
                    }
                    return solve_linear_system_simple(&a_reg, b);
                }
                l[[i, j]] = sum.sqrt();
            } else {
                l[[i, j]] = sum / l[[j, j]];
            }
        }
    }

    // Solve L * y = b (forward substitution)
    let mut y = Array1::zeros(n);
    for i in 0..n {
        let mut sum = b[i];
        for j in 0..i {
            sum -= l[[i, j]] * y[j];
        }
        y[i] = sum / l[[i, i]];
    }

    // Solve L^T * x = y (backward substitution)
    let mut x = Array1::zeros(n);
    for i in (0..n).rev() {
        let mut sum = y[i];
        for j in (i + 1)..n {
            sum -= l[[j, i]] * x[j];
        }
        x[i] = sum / l[[i, i]];
    }

    Ok(x)
}

/// Simple iterative solver as fallback.
fn solve_linear_system_simple(a: &Array2<f64>, b: &Array1<f64>) -> Result<Array1<f64>> {
    let n = a.nrows();
    let mut x = Array1::zeros(n);

    // Gauss-Seidel iteration
    for _ in 0..1000 {
        let mut max_change = 0.0_f64;
        for i in 0..n {
            let mut sum = b[i];
            for j in 0..n {
                if i != j {
                    sum -= a[[i, j]] * x[j];
                }
            }
            let new_val: f64 = sum / a[[i, i]];
            let diff: f64 = new_val - x[i];
            max_change = max_change.max(diff.abs());
            x[i] = new_val;
        }
        if max_change < 1e-12 {
            break;
        }
    }

    Ok(x)
}

/// Compute predicted f-statistics for all pairs in a graph.
#[allow(dead_code)]
pub fn compute_predicted_f2(
    graph: &AdmixGraph,
    populations: &[String],
) -> Result<HashMap<(String, String), f64>> {
    let mut result = HashMap::new();

    for i in 0..populations.len() {
        for j in (i + 1)..populations.len() {
            let f2 = graph.expected_f2(&populations[i], &populations[j])?;
            result.insert((populations[i].clone(), populations[j].clone()), f2);
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::{BlockDefinition, F2Blocks};
    use crate::genotype_matrix::AlleleFrequencies;
    use crate::genotypes::{Chromosome, PopulationId};
    use crate::qpgraph::graph::GraphBuilder;
    use ndarray::Array2;

    fn create_test_f2_blocks() -> F2Blocks {
        // Create synthetic data that matches a simple tree
        // A and B are sister groups, both diverged from root
        let pops = vec!["A", "B", "C"];
        let n_snps = 1000;

        let mut rng_state: u64 = 11111;
        let mut rand = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut freqs = Array2::zeros((n_snps, 3));
        let counts = Array2::from_elem((n_snps, 3), 100u32);

        for i in 0..n_snps {
            // Base frequency
            let base = rand() * 0.8 + 0.1;

            // A and B are close
            freqs[[i, 0]] = (base + rand() * 0.05).clamp(0.0, 1.0);
            freqs[[i, 1]] = (base + rand() * 0.05).clamp(0.0, 1.0);

            // C is more divergent
            freqs[[i, 2]] = (base + rand() * 0.15 - 0.075).clamp(0.0, 1.0);
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
    fn test_fit_simple_tree() {
        let f2_blocks = create_test_f2_blocks();

        // Create a simple tree topology
        let graph = GraphBuilder::new()
            .root("R")
            .internal("AB", "R", 0.01)
            .leaf("A", "AB", 0.01)
            .leaf("B", "AB", 0.01)
            .leaf("C", "R", 0.01)
            .build_unchecked();

        let result = fit_graph(&graph, &f2_blocks, None);
        assert!(result.is_ok(), "Fitting failed: {:?}", result.err());

        let result = result.unwrap();
        println!("Fit result: score={}, rmse={}", result.score, result.rmse);
        println!("Iterations: {}, converged: {}", result.iterations, result.converged);
    }

    #[test]
    fn test_fit_requires_matching_populations() {
        let f2_blocks = create_test_f2_blocks();

        // Graph with different population names
        let graph = GraphBuilder::new()
            .root("R")
            .leaf("X", "R", 0.1)
            .leaf("Y", "R", 0.1)
            .build_unchecked();

        let result = fit_graph(&graph, &f2_blocks, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_compute_predicted_f2() {
        let graph = GraphBuilder::new()
            .root("R")
            .leaf("A", "R", 0.1)
            .leaf("B", "R", 0.2)
            .build_unchecked();

        let populations = vec!["A".to_string(), "B".to_string()];
        let predicted = compute_predicted_f2(&graph, &populations).unwrap();

        assert_eq!(predicted.len(), 1);
        let f2_ab = predicted.get(&("A".to_string(), "B".to_string())).unwrap();
        assert!((f2_ab - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_fit_options() {
        let opts = FitOptions {
            max_iter: 100,
            tolerance: 1e-6,
            weighted: false,
            ..Default::default()
        };

        assert_eq!(opts.max_iter, 100);
        assert!(!opts.weighted);
    }
}
