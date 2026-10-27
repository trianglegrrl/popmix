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
pub mod dataset;
pub mod eigenstrat;
pub mod extract;
pub mod f2_storage;
pub mod merge;
pub mod npy;
pub mod packedancestrymap;
pub mod ped;
pub mod plink;
pub mod streaming;
pub mod tgeno;
pub mod vcf;
pub mod writers;

/// A genetic map position: a finite number. (convertf reads anything else
/// as 0, which would put the SNP at the start of its chromosome's map.)
pub(crate) fn parse_genetic_position(field: &str, path: &str, line: usize) -> Result<f64> {
    match field.parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(v),
        _ => Err(PopmixError::parse(
            path,
            line,
            format!("genetic position {field:?} is not a number"),
        )),
    }
}

/// Check that `range` is a range of SNP indices of a file with `n_snps` SNPs.
pub(crate) fn check_snp_range(range: &std::ops::Range<usize>, n_snps: usize) -> Result<()> {
    if range.start > range.end || range.end > n_snps {
        return Err(PopmixError::parameter(
            "SNP range",
            format!(
                "{}..{} is not within the {n_snps} SNPs of the file",
                range.start, range.end
            ),
        ));
    }
    Ok(())
}

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

    /// TGENO (.geno, .snp, .ind): packed genotypes stored individual by
    /// individual, with a TGENO header (AADR v66 and later)
    Tgeno,
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
        if (prefix_str.ends_with(".vcf") || prefix_str.ends_with(".vcf.gz"))
            && Path::new(prefix).exists()
        {
            return Ok(GenoFormat::Vcf);
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
            // Distinguish by checking for the TGENO or GENO magic header
            if tgeno::is_tgeno_format(&geno_path)? {
                return Ok(GenoFormat::Tgeno);
            }
            if packedancestrymap::is_packed_format(&geno_path)? {
                return Ok(GenoFormat::PackedAncestrymap);
            } else {
                return Ok(GenoFormat::Eigenstrat);
            }
        }

        Err(PopmixError::FileNotFound {
            path: format!("{} (no recognized genotype files found)", prefix_str).into(),
        })
    }

    /// Parse format from string name.
    ///
    /// Accepts various common names for each format.
    // returns popmix's Result, not FromStr's
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "ancestrymap" => Ok(GenoFormat::Ancestrymap),
            "eigenstrat" => Ok(GenoFormat::Eigenstrat),
            "ped" | "pedmap" => Ok(GenoFormat::Ped),
            "packedped" | "plink" | "bed" => Ok(GenoFormat::PackedPed),
            "packedancestrymap" => Ok(GenoFormat::PackedAncestrymap),
            "vcf" => Ok(GenoFormat::Vcf),
            "tgeno" | "transpose_packed" => Ok(GenoFormat::Tgeno),
            _ => Err(PopmixError::InvalidFileFormat {
                path: s.to_string().into(),
                details: format!(
                    "Unknown format '{}'. Valid formats: ancestrymap, eigenstrat, ped, packedped, packedancestrymap, tgeno, vcf",
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
            GenoFormat::Tgeno => "TGENO",
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
    compute_allele_frequencies_packedancestrymap, detect_ploidy_packedancestrymap,
    detect_ploidy_packedancestrymap_selected, PackedRows, SnpMajor, SnpMajorRows,
};
pub use vcf::{VcfReader, VcfSkipped};
pub use writers::{
    write_ancestrymap, write_ancestrymap_rows, write_eigenstrat, write_eigenstrat_rows,
    write_packedancestrymap, write_packedancestrymap_rows, write_packedped, write_packedped_rows,
    write_ped, write_ped_samples, write_plink, write_tgeno, write_tgeno_samples, SampleGenotypes,
};

/// Parse the non-empty lines of `content` with `parse(line number from 0,
/// trimmed line)`, in parallel: the file is cut at line ends into pieces of
/// about a megabyte, each piece's first line number counted first. The
/// results are in file order, and the error is the first bad line's.
pub(crate) fn parse_lines<T: Send>(
    content: &str,
    parse: impl Fn(usize, &str) -> Result<T> + Sync,
) -> Result<Vec<T>> {
    use rayon::prelude::*;
    const PIECE: usize = 1 << 20;
    let bytes = content.as_bytes();
    let mut starts = vec![0];
    while let Some(&last) = starts.last() {
        let at = last + PIECE;
        if at >= bytes.len() {
            break;
        }
        match bytes[at..].iter().position(|&b| b == b'\n') {
            Some(k) => starts.push(at + k + 1),
            None => break,
        }
    }
    let pieces: Vec<&str> = starts
        .iter()
        .zip(starts.iter().skip(1).chain([&bytes.len()]))
        .map(|(&a, &b)| &content[a..b])
        .collect();
    let lines_in: Vec<usize> = pieces
        .par_iter()
        .map(|p| p.bytes().filter(|&b| b == b'\n').count())
        .collect();
    let first_line: Vec<usize> = lines_in
        .iter()
        .scan(0, |n, &k| {
            let first = *n;
            *n += k;
            Some(first)
        })
        .collect();
    let parsed: Vec<Result<Vec<T>>> = pieces
        .par_iter()
        .zip(first_line)
        .map(|(piece, first)| {
            piece
                .lines()
                .enumerate()
                .map(|(k, line)| (first + k, line.trim()))
                .filter(|(_, line)| !line.is_empty())
                .map(|(n, line)| parse(n, line))
                .collect()
        })
        .collect();
    let mut out = Vec::new();
    for piece in parsed {
        out.extend(piece?);
    }
    Ok(out)
}

/// The first six whitespace-separated fields of a line, or the error naming
/// how many it has.
pub(crate) fn six_fields<'a>(path: &str, line_num: usize, line: &'a str) -> Result<[&'a str; 6]> {
    let mut fields = line.split_whitespace();
    let mut out = [""; 6];
    for (k, slot) in out.iter_mut().enumerate() {
        *slot = fields.next().ok_or_else(|| {
            PopmixError::parse(path, line_num + 1, format!("Expected 6 fields, got {k}"))
        })?;
    }
    Ok(out)
}

/// Genetic positions of PLINK .bim and .map files in Morgans.
///
/// EIGENSOFT and admixtools2 write and read them in Morgans (as do the
/// Reich lab's PLINK releases), PLINK 2 documents centimorgans; no
/// chromosome is 10 Morgans long, so positions beyond that mean the file is
/// in centimorgans, and they are divided by 100.
pub fn to_morgans(snps: &mut [popmix_core::Snp]) {
    let max = snps
        .iter()
        .map(|s| s.position_morgans.abs())
        .fold(0.0, f64::max);
    if max > 10.0 {
        for s in snps.iter_mut() {
            s.position_morgans /= 100.0;
        }
    }
}

#[cfg(test)]
mod morgan_tests {
    use super::*;
    use popmix_core::{Chromosome, Snp};

    #[test]
    fn test_centimorgans_are_detected() {
        let snp = |g: f64| Snp::new("s", Chromosome::Autosome(1), 1, g, 'A', 'G');
        let mut morgans = vec![snp(0.02013), snp(2.85)];
        to_morgans(&mut morgans);
        assert_eq!(morgans[1].position_morgans, 2.85);
        let mut cm = vec![snp(2.013), snp(285.0)];
        to_morgans(&mut cm);
        assert!(
            (cm[0].position_morgans - 0.02013).abs() < 1e-15
                && (cm[1].position_morgans - 2.85).abs() < 1e-15
        );
    }
}

#[cfg(test)]
mod parse_lines_tests {
    use super::*;

    #[test]
    fn lines_keep_their_order_and_numbers_across_pieces() {
        // about 3 MB: several pieces, with blank and CRLF lines among them
        let mut text = String::new();
        for i in 0..300_000 {
            text.push_str(&format!("{i} x\n"));
            if i % 1000 == 0 {
                text.push_str("  \r\n");
            }
        }
        let got = parse_lines(&text, |n, line| Ok((n, line.to_string()))).unwrap();
        assert_eq!(got.len(), 300_000);
        for (k, (n, line)) in got.iter().enumerate() {
            assert_eq!(line, &format!("{k} x"));
            // one blank line after line 0, 1000, 2000, ...
            assert_eq!(*n, k + k.div_ceil(1000));
        }
        let lines: Vec<&str> = text.lines().collect();
        let bad = lines.iter().position(|l| *l == "250000 x").unwrap();
        let err = parse_lines(&text, |n, line| {
            if line == "250000 x" || line == "260000 x" {
                Err(PopmixError::parse("f", n + 1, "bad"))
            } else {
                Ok(())
            }
        })
        .unwrap_err();
        assert!(err.to_string().contains(&format!("{}", bad + 1)), "{err}");
    }

    #[test]
    fn six_fields_or_how_many_there_are() {
        assert_eq!(
            six_fields("f", 0, "a b c d e f g").unwrap(),
            ["a", "b", "c", "d", "e", "f"]
        );
        let err = six_fields("f", 4, "a b c").unwrap_err().to_string();
        assert!(
            err.contains("Expected 6 fields, got 3") && err.contains('5'),
            "{err}"
        );
    }
}
