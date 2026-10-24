//! DATES: Distribution of Ancestry Tracts of Evolutionary Signals
//!
//! This module implements admixture dating via LD decay estimation,
//! based on the method described in Narasimhan et al. (2018).
//!
//! # Algorithm Overview
//!
//! DATES estimates admixture timing by measuring the exponential decay of
//! linkage disequilibrium (LD) with genetic distance. After admixture,
//! LD decays each generation due to recombination.
//!
//! ## Core Mathematical Model
//!
//! ```text
//! E[a(d)] = A × exp(-n × d) + K
//! ```
//!
//! Where:
//! - `a(d)` = weighted LD covariance at genetic distance d (in Morgans)
//! - `n` = generations since admixture (the parameter we estimate)
//! - `A` = amplitude ∝ 2αβ × F₂(source1, source2)² (mixture proportions × divergence)
//! - `K` = background LD constant (from inter-chromosomal pairs)
//! - `d` = genetic distance in Morgans (cM / 100)
//!
//! ## Weighted LD Statistic
//!
//! For SNP pairs (x, y) at genetic distance d:
//!
//! ```text
//! a(d) = Σ [ Cov(geno_x, geno_y) × w(x) × w(y) ]
//! ```
//!
//! Where:
//! - `w(x) = p_source1(x) - p_source2(x)` (allele frequency difference at SNP x)
//! - Covariance computed in the test population
//! - Weights from reference/source populations

mod compat;
mod dating;
mod exponential_fit;
mod ld_decay;
mod weights;

pub use compat::DatesCompat;
pub use dating::{dates, DatesOptions, DatesResult, JackknifeBlock};
pub use exponential_fit::{fit_exponential, weighted_jackknife, DecayFit, FitOptions};
pub use ld_decay::{
    compute_ld_decay, compute_ld_decay_by_chromosome, ChromosomeLd, LdBin, LdDecay, LdOptions,
};
pub use weights::{compute_weights, WeightResult};
