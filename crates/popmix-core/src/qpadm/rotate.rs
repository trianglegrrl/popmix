//! qpAdm rotate: Systematic model testing across source combinations.
//!
//! This module implements qpAdm model rotation, which tests all combinations
//! of candidate source populations to find the best-fitting admixture models.

use crate::fstats::F2Blocks;
use crate::genotypes::PopulationId;
use crate::qpadm::qpadm_impl::{qpadm, QpAdmOptions};
use crate::Result;
use itertools::Itertools;

/// Result of a single admixture model from rotation.
#[derive(Debug, Clone)]
pub struct RotationModel {
    /// Source populations used in this model.
    pub sources: Vec<PopulationId>,
    /// Estimated admixture weights (sum to 1).
    pub weights: Vec<f64>,
    /// Standard errors of weights.
    pub weight_std_errors: Vec<f64>,
    /// P-value for model fit (higher = better fit).
    pub p_value: f64,
    /// Chi-squared statistic.
    pub chi_sq: f64,
    /// Degrees of freedom.
    pub df: usize,
    /// Whether all weights are in [0, 1] range.
    pub is_feasible: bool,
}

/// Result of qpAdm rotation analysis.
#[derive(Debug, Clone)]
pub struct QpAdmRotateResult {
    /// Target population being modeled.
    pub target: PopulationId,
    /// Candidate source populations.
    pub candidates: Vec<PopulationId>,
    /// Outgroup populations.
    pub outgroups: Vec<PopulationId>,
    /// All tested models, sorted by p_value descending.
    pub models: Vec<RotationModel>,
    /// Total number of combinations tested.
    pub n_combinations: usize,
    /// Number of feasible models (weights in [0,1]).
    pub n_feasible: usize,
}

impl QpAdmRotateResult {
    /// Get the best model (highest p-value).
    pub fn best_model(&self) -> Option<&RotationModel> {
        self.models.first()
    }

    /// Get models passing a p-value threshold.
    pub fn passing_models(&self, p_threshold: f64) -> Vec<&RotationModel> {
        self.models
            .iter()
            .filter(|m| m.p_value > p_threshold)
            .collect()
    }

    /// Get feasible models only.
    pub fn feasible_models(&self) -> Vec<&RotationModel> {
        self.models.iter().filter(|m| m.is_feasible).collect()
    }
}

/// Options for qpAdm rotation.
#[derive(Debug, Clone)]
pub struct QpAdmRotateOptions {
    /// Minimum number of sources to test (default: 1).
    pub min_sources: usize,
    /// Maximum number of sources to test (default: all candidates).
    pub max_sources: Option<usize>,
    /// P-value threshold for "passing" models (default: 0.05).
    pub p_threshold: f64,
    /// Only return feasible models (default: false).
    pub feasible_only: bool,
    /// Options passed to qpAdm for each model.
    pub qpadm_opts: QpAdmOptions,
}

impl Default for QpAdmRotateOptions {
    fn default() -> Self {
        Self {
            min_sources: 1,
            max_sources: None,
            p_threshold: 0.05,
            feasible_only: false,
            qpadm_opts: QpAdmOptions::default(),
        }
    }
}

/// Run qpAdm rotation to test all source combinations.
///
/// Tests all C(n, k) combinations of candidate sources for k from min_sources
/// to max_sources. Results are sorted by p-value (descending).
///
/// # Arguments
/// * `f2_blocks` - Pre-computed F2 blocks
/// * `target` - Target population to model
/// * `candidates` - Candidate source populations
/// * `outgroups` - Outgroup populations (right populations)
/// * `options` - Optional settings
///
/// # Returns
/// QpAdmRotateResult with all tested models sorted by p-value
pub fn qpadm_rotate(
    f2_blocks: &F2Blocks,
    target: &str,
    candidates: &[&str],
    outgroups: &[&str],
    options: Option<QpAdmRotateOptions>,
) -> Result<QpAdmRotateResult> {
    let opts = options.unwrap_or_default();
    let max_sources = opts.max_sources.unwrap_or(candidates.len());
    let min_sources = opts.min_sources.max(1);

    let mut models = Vec::new();
    let mut n_combinations = 0;

    // Test all combinations from min_sources to max_sources
    for k in min_sources..=max_sources {
        for combo in candidates.iter().combinations(k) {
            n_combinations += 1;

            let sources: Vec<&str> = combo.iter().map(|&&s| s).collect();

            // Run qpAdm for this combination
            match qpadm(
                f2_blocks,
                target,
                &sources,
                outgroups,
                Some(opts.qpadm_opts.clone()),
            ) {
                Ok(result) => {
                    let is_feasible = result.weights.iter().all(|&w| w >= 0.0 && w <= 1.0);

                    // Skip non-feasible if requested
                    if opts.feasible_only && !is_feasible {
                        continue;
                    }

                    models.push(RotationModel {
                        sources: result.sources,
                        weights: result.weights,
                        weight_std_errors: result.weight_std_errors,
                        p_value: result.p_value,
                        chi_sq: result.chi_sq,
                        df: result.df,
                        is_feasible,
                    });
                }
                Err(_) => {
                    // Skip failed models (e.g., singular matrix)
                    continue;
                }
            }
        }
    }

    // Sort by p-value descending (best fits first)
    models.sort_by(|a, b| {
        b.p_value
            .partial_cmp(&a.p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let n_feasible = models.iter().filter(|m| m.is_feasible).count();

    Ok(QpAdmRotateResult {
        target: PopulationId::new(target),
        candidates: candidates.iter().map(|&s| PopulationId::new(s)).collect(),
        outgroups: outgroups.iter().map(|&s| PopulationId::new(s)).collect(),
        models,
        n_combinations,
        n_feasible,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::blocks::BlockDefinition;
    use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
    use crate::genotypes::{Chromosome, Sample, Snp};
    use ndarray::Array2;

    /// Create test F2 blocks with multiple populations.
    fn create_test_f2_blocks() -> F2Blocks {
        // Create populations: Target, Source1, Source2, Source3, Out1, Out2
        let mut samples = Vec::new();
        for (pop_name, count) in [
            ("Target", 10),
            ("Source1", 10),
            ("Source2", 10),
            ("Source3", 10),
            ("Out1", 10),
            ("Out2", 10),
        ] {
            for i in 0..count {
                samples.push(Sample::new(
                    format!("{}_{}", pop_name, i),
                    PopulationId::new(pop_name),
                ));
            }
        }

        // Create SNPs across blocks
        let snps: Vec<Snp> = (0..100)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 1000) as u32,
                    i as f64 * 0.2, // 0.2 cM spacing
                    'A',
                    'G',
                )
            })
            .collect();

        let n_samples = samples.len();
        let n_snps = snps.len();
        let mut data = Array2::zeros((n_snps, n_samples));

        // Create frequency patterns
        for snp_idx in 0..n_snps {
            let base = (snp_idx * 7 % 100) as f64 / 100.0;
            for sample_idx in 0..n_samples {
                let pop_idx = sample_idx / 10;
                let freq = match pop_idx {
                    0 => 0.3 + base * 0.2,  // Target
                    1 => 0.2 + base * 0.2,  // Source1
                    2 => 0.5 + base * 0.2,  // Source2
                    3 => 0.4 + base * 0.2,  // Source3
                    4 => 0.1 + base * 0.1,  // Out1
                    5 => 0.15 + base * 0.1, // Out2
                    _ => 0.5,
                };
                let threshold = freq * 2.0;
                let noise = ((sample_idx * 13 + snp_idx * 17) % 100) as f64 / 100.0 * 2.0;
                data[[snp_idx, sample_idx]] = if noise < threshold { 2 } else { 0 };
            }
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        let af = AlleleFrequencies::from_genotypes(&geno);
        let block_def = BlockDefinition::new(0.05);
        F2Blocks::from_allele_frequencies(&af, &block_def)
    }

    #[test]
    fn test_qpadm_rotate_basic() {
        let f2_blocks = create_test_f2_blocks();

        let result = qpadm_rotate(
            &f2_blocks,
            "Target",
            &["Source1", "Source2", "Source3"],
            &["Out1", "Out2"],
            None,
        )
        .unwrap();

        // Should test all combinations: C(3,1) + C(3,2) + C(3,3) = 3 + 3 + 1 = 7
        assert_eq!(result.n_combinations, 7);
        assert!(!result.models.is_empty(), "Should have some valid models");

        // Models should be sorted by p-value descending
        for i in 1..result.models.len() {
            assert!(
                result.models[i - 1].p_value >= result.models[i].p_value,
                "Models should be sorted by p-value descending"
            );
        }
    }

    #[test]
    fn test_qpadm_rotate_min_max_sources() {
        let f2_blocks = create_test_f2_blocks();

        let opts = QpAdmRotateOptions {
            min_sources: 2,
            max_sources: Some(2),
            ..Default::default()
        };

        let result = qpadm_rotate(
            &f2_blocks,
            "Target",
            &["Source1", "Source2", "Source3"],
            &["Out1", "Out2"],
            Some(opts),
        )
        .unwrap();

        // Should only test 2-source combinations: C(3,2) = 3
        assert_eq!(result.n_combinations, 3);

        // All models should have exactly 2 sources
        for model in &result.models {
            assert_eq!(model.sources.len(), 2);
        }
    }

    #[test]
    fn test_qpadm_rotate_feasible_only() {
        let f2_blocks = create_test_f2_blocks();

        let opts = QpAdmRotateOptions {
            feasible_only: true,
            ..Default::default()
        };

        let result = qpadm_rotate(
            &f2_blocks,
            "Target",
            &["Source1", "Source2", "Source3"],
            &["Out1", "Out2"],
            Some(opts),
        )
        .unwrap();

        // All returned models should be feasible
        for model in &result.models {
            assert!(model.is_feasible, "All models should be feasible");
            for weight in &model.weights {
                assert!(
                    *weight >= 0.0 && *weight <= 1.0,
                    "Weights should be in [0, 1]"
                );
            }
        }
    }

    #[test]
    fn test_qpadm_rotate_helper_methods() {
        let f2_blocks = create_test_f2_blocks();

        let result = qpadm_rotate(
            &f2_blocks,
            "Target",
            &["Source1", "Source2"],
            &["Out1", "Out2"],
            None,
        )
        .unwrap();

        // Test best_model
        if let Some(best) = result.best_model() {
            assert!(best.p_value >= 0.0);
        }

        // Test passing_models
        let passing = result.passing_models(0.01);
        for model in passing {
            assert!(model.p_value > 0.01);
        }

        // Test feasible_models
        let feasible = result.feasible_models();
        for model in feasible {
            assert!(model.is_feasible);
        }
    }
}
