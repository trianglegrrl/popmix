//! # popmix-io
//!
//! File format I/O for population genetics data.
//!
//! This crate handles reading and writing of:
//! - ANCESTRYMAP sparse format (.ancestrymapgeno, .snp, .ind)
//! - EIGENSTRAT text format (.geno, .snp, .ind)
//! - PED text format (.ped, .map, .pedind)
//! - PACKEDPED binary format (.bed, .bim, .fam) - aka PLINK
//! - PACKEDANCESTRYMAP binary format (.geno, .snp, .ind)
//! - Precomputed F2 block storage
//! - NumPy .npy arrays

use std::path::Path;

use popmix_core::error::{PopmixError, Result};

pub mod ancestrymap;
pub mod eigenstrat;
pub mod f2_storage;
pub mod merge;
pub mod npy;
pub mod packedancestrymap;
pub mod ped;
pub mod plink;
pub mod streaming;
pub mod vcf;
pub mod writers;

/// Supported genotype file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GenoFormat {
    /// ANCESTRYMAP sparse format (.ancestrymapgeno, .snp, .ind)
    /// One line per non-missing genotype: SNP_ID Sample_ID Dosage
    Ancestrymap,

    /// EIGENSTRAT text format (.geno, .snp, .ind)
    /// One line per SNP, one character per sample (0/1/2/9)
    Eigenstrat,

    /// PED text format (.ped, .map)
    /// One line per sample with allele pairs
    Ped,

    /// PACKEDPED binary format (.bed, .bim, .fam) - PLINK format
    /// SNP-major packed genotypes
    PackedPed,

    /// PACKEDANCESTRYMAP binary format (.geno, .snp, .ind)
    /// Packed genotypes with GENO header
    PackedAncestrymap,

    /// VCF format (.vcf, .vcf.gz)
    /// Variant Call Format, optionally gzip-compressed
    Vcf,
}

impl GenoFormat {
    /// Detect format from file extensions at given prefix.
    ///
    /// Checks for format-specific file extensions and, for formats
    /// that share extensions (EIGENSTRAT vs PACKEDANCESTRYMAP),
    /// examines file contents.
    ///
    /// For VCF files, pass the full path (including .vcf or .vcf.gz extension).
    pub fn detect(prefix: impl AsRef<Path>) -> Result<Self> {
        let prefix = prefix.as_ref();
        let prefix_str = prefix.to_string_lossy();

        // Check for VCF format (full path with extension)
        if prefix_str.ends_with(".vcf") || prefix_str.ends_with(".vcf.gz") {
            if Path::new(prefix).exists() {
                return Ok(GenoFormat::Vcf);
            }
        }

        // Check for VCF format with prefix (add .vcf or .vcf.gz)
        let vcf_path = format!("{}.vcf", prefix_str);
        let vcf_gz_path = format!("{}.vcf.gz", prefix_str);
        if Path::new(&vcf_gz_path).exists() {
            return Ok(GenoFormat::Vcf);
        }
        if Path::new(&vcf_path).exists() {
            return Ok(GenoFormat::Vcf);
        }

        // Check for PACKEDPED (PLINK) format: .bed file
        let bed_path = format!("{}.bed", prefix_str);
        if Path::new(&bed_path).exists() {
            return Ok(GenoFormat::PackedPed);
        }

        // Check for PED format: .ped file
        let ped_path = format!("{}.ped", prefix_str);
        if Path::new(&ped_path).exists() {
            return Ok(GenoFormat::Ped);
        }

        // Check for ANCESTRYMAP format: .ancestrymapgeno file
        let ancestrymap_path = format!("{}.ancestrymapgeno", prefix_str);
        if Path::new(&ancestrymap_path).exists() {
            return Ok(GenoFormat::Ancestrymap);
        }

        // Check for EIGENSTRAT or PACKEDANCESTRYMAP format: .geno file
        let geno_path = format!("{}.geno", prefix_str);
        if Path::new(&geno_path).exists() {
            // Distinguish by checking for GENO magic header
            if packedancestrymap::is_packed_format(&geno_path)? {
                return Ok(GenoFormat::PackedAncestrymap);
            } else {
                return Ok(GenoFormat::Eigenstrat);
            }
        }

        Err(PopmixError::FileNotFound {
            path: format!("{} (no recognized genotype files found)", prefix_str),
        })
    }

    /// Parse format from string name.
    ///
    /// Accepts various common names for each format.
    pub fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "ancestrymap" => Ok(GenoFormat::Ancestrymap),
            "eigenstrat" => Ok(GenoFormat::Eigenstrat),
            "ped" | "pedmap" => Ok(GenoFormat::Ped),
            "packedped" | "plink" | "bed" => Ok(GenoFormat::PackedPed),
            "packedancestrymap" => Ok(GenoFormat::PackedAncestrymap),
            "vcf" => Ok(GenoFormat::Vcf),
            _ => Err(PopmixError::InvalidFileFormat {
                file: s.to_string(),
                details: format!(
                    "Unknown format '{}'. Valid formats: ancestrymap, eigenstrat, ped, packedped, packedancestrymap, vcf",
                    s
                ),
            }),
        }
    }

    /// Get the canonical name for this format.
    pub fn name(&self) -> &'static str {
        match self {
            GenoFormat::Ancestrymap => "ANCESTRYMAP",
            GenoFormat::Eigenstrat => "EIGENSTRAT",
            GenoFormat::Ped => "PED",
            GenoFormat::PackedPed => "PACKEDPED",
            GenoFormat::PackedAncestrymap => "PACKEDANCESTRYMAP",
            GenoFormat::Vcf => "VCF",
        }
    }
}

// Re-exports
pub use ancestrymap::AncestryMapReader;
pub use eigenstrat::EigenstratReader;
pub use f2_storage::{read_f2_blocks, write_f2_blocks};
pub use npy::{read_npy_f64, write_npy_f64};
pub use packedancestrymap::{is_packed_format, PackedAncestryMapReader};
pub use ped::PedReader;
pub use plink::PlinkReader;
pub use streaming::{
    compute_allele_frequencies_packedancestrymap, detect_ploidy_packedancestrymap, PackedRows,
};
pub use vcf::VcfReader;
pub use writers::{
    write_ancestrymap, write_eigenstrat, write_packedancestrymap, write_packedped, write_ped,
    write_plink,
};
