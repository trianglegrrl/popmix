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
pub mod f2_arrays;
pub mod f2_blocks;
mod f2_parts;
pub mod f3_f4;
pub mod geno_fstats;
pub mod jackknife;
mod one_allele;
mod provenance;
pub mod snp_filter;

pub use blocks::{BlockAssignment, BlockCoordinate, BlockDefinition, CoordinateEvidence};
pub use f2::{f2_estimate, f2_single_snp, F2Result};
pub use f2_arrays::F2Arrays;
pub use f2_blocks::{F2Blocks, FstBlocks};
pub use f3_f4::{F3Result, F4Result, FstResult};
pub use geno_fstats::{
    f3_from_genotypes, f4_from_genotypes, f4_from_genotypes_requiring, F3GenoOptions, F4Mode,
    GenoPopulations, GenoStat, GenoStatBlocks, GenoStats, GenotypeRows, MatrixRows,
};
pub use jackknife::{
    block_jackknife, block_jackknife_cov, jackknife_cov_from_loo, jackknife_from_loo,
    jackknife_with_total, leave_one_out, leave_one_out_nafix, p_from_z, z_score,
};
pub use snp_filter::F2Options;
