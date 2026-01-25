//! # popmix-io
//!
//! File format I/O for population genetics data.
//!
//! This crate handles reading and writing of:
//! - PLINK binary format (.bed, .bim, .fam)
//! - EIGENSTRAT text format (.geno, .snp, .ind)
//! - PACKEDANCESTRYMAP binary format (.geno, .snp, .ind)
//! - Precomputed F2 block storage
//! - NumPy .npy arrays

pub mod eigenstrat;
pub mod f2_storage;
pub mod npy;
pub mod packedancestrymap;
pub mod plink;
pub mod streaming;

// Re-exports
pub use eigenstrat::EigenstratReader;
pub use f2_storage::{read_f2_blocks, write_f2_blocks};
pub use npy::{read_npy_f64, write_npy_f64};
pub use packedancestrymap::{is_packed_format, PackedAncestryMapReader};
pub use plink::PlinkReader;
pub use streaming::{
    compute_allele_frequencies_packedancestrymap, detect_ploidy_packedancestrymap,
};
