//! Principal component analysis, as EIGENSOFT's smartpca.

pub mod smartpca;
pub mod tracy_widom;
mod tw_table;

pub use smartpca::{smartpca, Outlier, PackedGenotypes, PcaOptions, PcaResult, SampleRole};
pub use tracy_widom::{tracy_widom_table, tw_tail, EigensoftCompat, TracyWidomRow};
