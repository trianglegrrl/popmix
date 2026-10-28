//! qpAdm and qpWave admixture analysis module.
//!
//! This module implements:
//! - qpWave: Test admixture wave count between population sets
//! - qpAdm: Estimate admixture weights for a target population
//! - qpF4ratio: Estimate admixture proportions via f4 ratios
//!
//! Both methods are based on F4 statistics computed between "left" and "right"
//! population sets. The right populations serve as outgroups that differentiate
//! ancestry sources.

mod compat;
mod f4_matrix;
mod f4ratio;
mod fit;
mod geno;
mod qpadm_impl;
mod qpadm_p;
mod qpwave;
mod qpwave_pairs;
mod rotate;

pub use compat::Admixtools2Compat;
pub use f4_matrix::F4Matrix;
pub use f4ratio::{f4ratio, F4RatioOptions, F4RatioResult};
pub use fit::PopDropRow;
pub use geno::f4_matrix_from_genotypes;
pub use qpadm_impl::{qpadm, qpadm_on, QpAdmOptions, QpAdmResult};
pub use qpadm_p::{qpadm_p, QpAdmP};
pub use qpwave::{qpwave, qpwave_on, QpWaveOptions, QpWaveResult, RankTest};
pub use qpwave_pairs::{qpwave_pairs, QpWavePairResult, QpWavePairsResult};
pub use rotate::{qpadm_rotate, QpAdmRotateOptions, QpAdmRotateResult, RotationModel};
