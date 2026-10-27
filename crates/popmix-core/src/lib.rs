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
//!
//! ## Linear algebra
//!
//! Arrays in the public API are [`ndarray`]'s. Inside, small dense systems
//! (qpAdm's SVD and inverses, qpGraph's and ADMIXTURE's constrained solves)
//! use nalgebra, and PCA's one large symmetric eigendecomposition uses faer,
//! whose blocked solver is about four times nalgebra's speed at 2,000
//! samples (2.0 s against 8.0 s on one core, 2026-10-27), and the gap grows
//! as n^3. No hot loop converts between them. Neither appears in a public
//! signature (`tests/python/test_dependencies.py`).


pub mod admixture;
pub mod dates;
pub mod error;
pub mod fstats;
pub mod genotype_matrix;
pub mod genotypes;
pub mod kinship;
mod options;
pub mod pca;
#[cfg(feature = "python")]
pub mod python;
pub mod qpadm;
pub mod qpgraph;
pub mod rcompat;
pub mod rolloff;
pub mod sfs;
pub mod threading;
pub mod warnings;

// Re-exports for convenience
pub use admixture::{
    admixture, admixture_cv, admixture_loglik, admixture_projection, admixture_supervised,
    AdmixtureFit, AdmixtureOptions,
};
pub use error::{ErrorKind, IoResultExt, PopmixError, Result};
pub use fstats::{
    f2_estimate, f2_single_snp, f3_from_genotypes, f4_from_genotypes, BlockAssignment,
    BlockDefinition, CoordinateEvidence, F2Arrays, F2Blocks, F2Options, F2Result, F3GenoOptions,
    F3Result, F4Mode, F4Result, GenoPopulations, GenoStat, GenoStatBlocks, GenoStats, GenotypeRows,
    MatrixRows,
};
pub use genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
pub use genotypes::{Chromosome, Genotype, Ploidy, PopulationId, Sample, Sex, Snp};
pub use kinship::{
    king_robust, read_kinship, BitGenotypes, KingPair, ReadCompat, ReadNorm, ReadOptions,
    ReadResult,
};
pub use pca::{
    smartpca, EigensoftCompat, Outlier, PackedGenotypes, PcaOptions, PcaResult, SampleRole,
    TracyWidomRow,
};
pub use qpadm::{
    f4ratio, qpadm, qpadm_p, qpadm_rotate, qpwave, qpwave_pairs, Admixtools2Compat, F4Matrix,
    F4RatioOptions, F4RatioResult, PopDropRow, QpAdmOptions, QpAdmP, QpAdmResult,
    QpAdmRotateOptions, QpAdmRotateResult, QpWaveOptions, QpWavePairResult, QpWavePairsResult,
    QpWaveResult, RankTest, RotationModel,
};
pub use qpgraph::{
    all_graphs, find_graphs, graph_key, qpgraph, AdmixGraph, AdmixNode, EdgeType, GraphBuilder,
    GraphSearchMethod, GraphSearchOptions, GraphSearchResult, NodeId, QpGraphEdge, QpGraphF3,
    QpGraphOptions, QpGraphResult, QpGraphStart, ScoredGraph,
};
pub use rolloff::{fit_rolloff, rolloff, RolloffBins, RolloffFit, RolloffOptions, RolloffResult};
pub use sfs::{joint_sfs, projected_sfs, JointSfsCompat, JointSfsSnps, SfsSummary, Spectrum};
pub use threading::init_thread_pool;
pub use warnings::Warning;

// DATES re-exports
pub use dates::{
    compute_ld_decay, compute_weights, dates, fit_exponential, weighted_jackknife, DatesOptions,
    DatesResult, DecayFit, LdBin, LdOptions, WeightResult,
};
// The exponential fit's options are dates::DecayFitOptions
