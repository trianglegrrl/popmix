//! # popmix-core
//!
//! Core algorithms for population genetics F-statistics and admixture analysis.
//!
//! This crate provides the computational foundation for popmix, implementing:
//! - F-statistics (f2, f3, f4, Fst)
//! - qpAdm/qpWave admixture weight estimation
//! - qpGraph admixture graph fitting
//! - DATES admixture dating via LD decay
//! - Jackknife and bootstrap resampling
//!
//! All algorithms are optimized for parallel execution using rayon.

pub mod dates;
pub mod error;
pub mod fstats;
pub mod genotype_matrix;
pub mod genotypes;
pub mod qpadm;
pub mod qpgraph;
pub mod threading;

// Re-exports for convenience
pub use error::{PopmixError, Result};
pub use fstats::{
    f2_estimate, f2_single_snp, BlockAssignment, BlockDefinition, F2Blocks, F2Options, F2Result,
    F3Result, F4Result,
};
pub use genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
pub use genotypes::{Chromosome, Genotype, Ploidy, PopulationId, Sample, Snp};
pub use qpadm::{
    f4ratio, qpadm, qpadm_rotate, qpwave, qpwave_pairs, Admixtools2Compat, F4Matrix,
    F4RatioOptions, F4RatioResult, PopDropRow, QpAdmOptions, QpAdmResult, QpAdmRotateOptions,
    QpAdmRotateResult, QpWaveOptions, QpWavePairResult, QpWavePairsResult, QpWaveResult, RankTest,
    RotationModel,
};
pub use qpgraph::{
    find_graphs, fit_graph, AdmixGraph, AdmixNode, EdgeType, FitOptions, GraphBuilder,
    GraphFitResult, GraphSearchOptions, GraphSearchResult, NodeId,
};
pub use threading::init_thread_pool;

// DATES re-exports
pub use dates::{
    compute_ld_decay, compute_weights, dates, fit_exponential, weighted_jackknife, DatesOptions,
    DatesResult, DecayFit, LdBin, LdOptions, WeightResult,
};
// Note: FitOptions for dates is available as dates::FitOptions to avoid collision with qpgraph::FitOptions
