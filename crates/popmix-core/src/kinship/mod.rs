//! Relatedness between samples: KING-robust kinship for diploid genotypes,
//! READ-style kinship for pseudo-haploid (ancient) ones.

pub mod king;
pub mod read;

pub use king::{king_robust, BitGenotypes, KingPair};
pub use read::{
    read_kinship, read_thresholds_2pow, FirstDegreeType, ReadCompat, ReadNorm, ReadOptions,
    ReadPair, ReadRelationship, ReadResult, READ_THRESHOLDS,
};
