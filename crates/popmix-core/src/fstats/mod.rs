//! F-statistics computation module.
//!
//! This module implements the core F-statistics used in population genetics:
//! - f2: Genetic drift between pairs of populations
//! - f3: Admixture test (3-population test)
//! - f4: D-statistic / ABBA-BABA test
//! - Fst: Fixation index
//!
//! All statistics are computed using blocked estimation with jackknife
//! resampling for standard error estimation.

pub mod blocks;
pub mod f2;
pub mod f2_blocks;
pub mod f3_f4;
pub mod jackknife;

pub use blocks::{BlockAssignment, BlockDefinition};
pub use f2::{f2_estimate, f2_single_snp, F2Result};
pub use f2_blocks::{F2Blocks, F2Options};
pub use jackknife::{block_jackknife, block_jackknife_cov};
pub use f3_f4::{F3Result, F4Result, FstResult};
