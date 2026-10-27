//! File conversion and merging for the bindings: `convertf` and `merge`,
//! streamed from inputs read in place where the formats allow, and the
//! mapped-input helpers the other readers share.

use std::collections::HashMap;

use popmix_core::python::{input_error, internal_error};
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;

use popmix_core::{GenotypeMatrix, Sample, Snp};
use popmix_io::tgeno::TgenoReader;
use popmix_io::{
    merge::merge_genotypes, write_ancestrymap, write_eigenstrat, write_packedancestrymap,
    write_packedped, write_ped, write_tgeno, AncestryMapReader, EigenstratReader, GenoFormat,
    PackedAncestryMapReader, PedReader, PlinkReader,
};

use crate::compute;
use crate::warnings::{open_vcf, with_warnings};

/// Convert genotype data between file formats.
///
/// Supported formats:
/// - ancestrymap: Sparse text format (.ancestrymapgeno, .snp, .ind)
/// - eigenstrat: Dense text format (.geno, .snp, .ind)
/// - packedancestrymap: Binary format (.geno, .snp, .ind)
/// - plink/packedped: Binary PLINK format (.bed, .bim, .fam)
/// - ped: Text PLINK format (.ped, .map, .pedind)
///
/// Args:
///     input_prefix: Path prefix for the input dataset
///     output_prefix: Path prefix for the output dataset
///     output_format: Output format name (see above)
///     input_format: Optional input format (auto-detected if not specified)
///
/// Returns:
///     Dictionary with conversion info (n_samples, n_snps, output_format)
#[pyfunction]
#[pyo3(signature = (input_prefix, output_prefix, output_format, *, input_format=None))]
pub(crate) fn convertf(
    py: Python<'_>,
    input_prefix: &str,
    output_prefix: &str,
    output_format: &str,
    input_format: Option<&str>,
) -> PyResult<HashMap<String, Py<PyAny>>> {
    let (n_samples, n_snps, in_format, out_format) = compute(py, || {
        convertf_files(input_prefix, output_prefix, output_format, input_format)
    })?;
    let mut result: HashMap<String, Py<PyAny>> = HashMap::new();
    result.insert("n_samples".to_string(), n_samples.into_py_any(py)?);
    result.insert("n_snps".to_string(), n_snps.into_py_any(py)?);
    result.insert(
        "input_format".to_string(),
        in_format.name().into_py_any(py)?,
    );
    result.insert(
        "output_format".to_string(),
        out_format.name().into_py_any(py)?,
    );
    Ok(result)
}

/// A genotype file mapped to be read in place, with its samples and SNPs.
pub(crate) type MappedInput = (popmix_io::SnpMajor, Vec<Sample>, Vec<Snp>);

/// A packed .geno, PLINK .bed or text EIGENSTRAT input, mapped to be read in
/// place, or a VCF read once into packed genotypes in memory; None for the
/// other formats (`popmix_io::dataset::open_mapped_as`).
pub(crate) fn mapped_input(prefix: &str, format: GenoFormat) -> PyResult<Option<MappedInput>> {
    let mapped = with_warnings(|warn| popmix_io::dataset::open_mapped_as(prefix, format, warn))?;
    Ok(mapped.map(|m| (m.file, m.samples, m.snps)))
}

/// Can `format` be written from a file read in place: one SNP at a time, or
/// (PED, TGENO) a block of samples at a time?
fn written_streamed(format: GenoFormat) -> bool {
    written_by_rows(format) || by_samples(format)
}

/// Is `format` stored sample by sample (written by [`write_samples`])?
fn by_samples(format: GenoFormat) -> bool {
    matches!(format, GenoFormat::Ped | GenoFormat::Tgeno)
}

/// Write PED or TGENO from genotypes read a block of samples at a time.
fn write_samples(
    format: GenoFormat,
    output_prefix: &str,
    samples: &[Sample],
    snps: &[Snp],
    genotypes: &dyn popmix_io::SampleGenotypes,
) -> PyResult<()> {
    if format == GenoFormat::Ped {
        popmix_io::write_ped_samples(output_prefix, samples, snps, genotypes)
    } else {
        popmix_io::write_tgeno_samples(output_prefix, samples, snps, genotypes)
    }
    .map_err(PyErr::from)
}

/// Can `format` be written one SNP at a time?
fn written_by_rows(format: GenoFormat) -> bool {
    matches!(
        format,
        GenoFormat::Eigenstrat
            | GenoFormat::Ancestrymap
            | GenoFormat::PackedAncestrymap
            | GenoFormat::PackedPed
    )
}

/// Write `rows` (one per SNP of `snps`) in a format of [`written_by_rows`].
fn write_rows(
    format: GenoFormat,
    output_prefix: &str,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn popmix_core::GenotypeRows,
) -> PyResult<()> {
    use popmix_io::{
        write_ancestrymap_rows, write_eigenstrat_rows, write_packedancestrymap_rows,
        write_packedped_rows,
    };
    let out = std::path::Path::new(output_prefix);
    match format {
        GenoFormat::Eigenstrat => write_eigenstrat_rows(out, samples, snps, rows),
        GenoFormat::Ancestrymap => write_ancestrymap_rows(out, samples, snps, rows),
        GenoFormat::PackedAncestrymap => write_packedancestrymap_rows(out, samples, snps, rows),
        GenoFormat::PackedPed => write_packedped_rows(out, samples, snps, rows),
        other => {
            return Err(internal_error(format!(
                "{} is not written SNP by SNP",
                other.name()
            )))
        }
    }
    .map_err(PyErr::from)
}

/// `convertf` one SNP at a time, when the input is read in place and the
/// output written SNP by SNP: the sample and SNP counts, or None.
fn convertf_streamed(
    input_prefix: &str,
    output_prefix: &str,
    in_format: GenoFormat,
    out_format: GenoFormat,
) -> PyResult<Option<(usize, usize)>> {
    if !written_streamed(out_format) {
        return Ok(None);
    }
    let Some((file, samples, snps)) = mapped_input(input_prefix, in_format)? else {
        return Ok(None);
    };
    if by_samples(out_format) {
        write_samples(out_format, output_prefix, &samples, &snps, &file)?;
    } else {
        let mut rows = popmix_io::SnpMajorRows::from_file(
            file,
            vec![true; snps.len()],
            (0..samples.len()).collect(),
        )?;
        write_rows(out_format, output_prefix, &samples, &snps, &mut rows)?;
    }
    Ok(Some((samples.len(), snps.len())))
}

/// `convertf`'s work: the sample and SNP counts and both formats.
fn convertf_files(
    input_prefix: &str,
    output_prefix: &str,
    output_format: &str,
    input_format: Option<&str>,
) -> PyResult<(usize, usize, GenoFormat, GenoFormat)> {
    use std::path::Path;

    // Detect or parse input format
    let in_format = match input_format {
        Some(fmt) => fmt.parse::<GenoFormat>()?,
        None => GenoFormat::detect(input_prefix)?,
    };

    // Parse output format
    let out_format = output_format.parse::<GenoFormat>()?;

    // From a packed .geno or a PLINK .bed to a format written SNP by SNP:
    // stream, without holding the genotypes in memory
    if let Some(done) = convertf_streamed(input_prefix, output_prefix, in_format, out_format)? {
        return Ok((done.0, done.1, in_format, out_format));
    }

    // Read the input data
    let geno = match in_format {
        GenoFormat::PackedPed => {
            let reader = PlinkReader::new(input_prefix)?;
            reader.read_all()?
        }
        GenoFormat::PackedAncestrymap => {
            let reader = PackedAncestryMapReader::new(input_prefix)?;
            reader.read_all()?
        }
        GenoFormat::Eigenstrat => {
            let reader = EigenstratReader::new(input_prefix)?;
            reader.read_all()?
        }
        GenoFormat::Ped => {
            let reader = PedReader::new(input_prefix)?;
            reader.read_all()?
        }
        GenoFormat::Ancestrymap => {
            let reader = AncestryMapReader::new(input_prefix)?;
            reader.read_all()?
        }
        GenoFormat::Vcf => {
            let reader = open_vcf(input_prefix)?;
            reader.read_all()?
        }
        GenoFormat::Tgeno => TgenoReader::new(input_prefix).and_then(|r| r.read_all())?,
    };

    let n_samples = geno.n_samples();
    let n_snps = geno.n_snps();

    // Write the output data
    let out_path = Path::new(output_prefix);
    match out_format {
        GenoFormat::Ancestrymap => {
            write_ancestrymap(out_path, &geno)?;
        }
        GenoFormat::Eigenstrat => {
            write_eigenstrat(out_path, &geno)?;
        }
        GenoFormat::PackedAncestrymap => {
            write_packedancestrymap(out_path, &geno)?;
        }
        GenoFormat::PackedPed => {
            write_packedped(out_path, &geno)?;
        }
        GenoFormat::Ped => {
            write_ped(out_path, &geno)?;
        }
        GenoFormat::Tgeno => {
            write_tgeno(out_path, &geno)?;
        }
        GenoFormat::Vcf => {
            return Err(input_error(
                "VCF output format is not supported. Use plink, eigenstrat, or packedancestrymap.",
            ));
        }
    }

    Ok((n_samples, n_snps, in_format, out_format))
}

/// Merge two genotype datasets.
///
/// Merges two genotype files by SNP intersection. SNPs are matched by
/// chromosome and position (and by name where several share one). Alleles
/// are reconciled as AdmixTools mergeit does: swapped alleles are recoded,
/// a SNP on the other strand is complemented, the unknown allele of a
/// monomorphic SNP is filled in, and A/T, C/G and third-allele SNPs are left
/// out and counted.
///
/// Args:
///     file_a: Path prefix for first dataset
///     file_b: Path prefix for second dataset
///     output: Path prefix for output dataset
///     output_format: Output format ("plink", "eigenstrat", "packedancestrymap", etc.)
///     input_format_a: Optional input format for file_a (auto-detected if not specified)
///     input_format_b: Optional input format for file_b (auto-detected if not specified)
///
/// Returns:
///     Dictionary with merge statistics:
///     - snps_only_a: SNPs only in file A
///     - snps_only_b: SNPs only in file B
///     - snps_identical: SNPs with identical alleles
///     - snps_flipped: SNPs with flipped alleles (recoded)
///     - snps_strand_ambiguous: SNPs skipped due to A/T or C/G ambiguity
///     - snps_mismatch: SNPs with mismatched alleles
///     - samples_total: Total samples in merged output
///     - snps_merged: Total SNPs in merged output
#[pyfunction]
#[pyo3(signature = (file_a, file_b, output, output_format, *, input_format_a=None, input_format_b=None))]
pub(crate) fn merge(
    py: Python<'_>,
    file_a: &str,
    file_b: &str,
    output: &str,
    output_format: &str,
    input_format_a: Option<&str>,
    input_format_b: Option<&str>,
) -> PyResult<HashMap<String, Py<PyAny>>> {
    use std::path::Path;
    let work = || -> PyResult<_> {
        // Helper to read genotypes from a prefix
        fn read_geno(prefix: &str, format: Option<&str>) -> PyResult<GenotypeMatrix> {
            let fmt = match format {
                Some(f) => f.parse::<GenoFormat>()?,
                None => GenoFormat::detect(prefix)?,
            };

            match fmt {
                GenoFormat::PackedPed => {
                    let reader = PlinkReader::new(prefix)?;
                    reader.read_all().map_err(PyErr::from)
                }
                GenoFormat::PackedAncestrymap => {
                    let reader = PackedAncestryMapReader::new(prefix)?;
                    reader.read_all().map_err(PyErr::from)
                }
                GenoFormat::Eigenstrat => {
                    let reader = EigenstratReader::new(prefix)?;
                    reader.read_all().map_err(PyErr::from)
                }
                GenoFormat::Ped => {
                    let reader = PedReader::new(prefix)?;
                    reader.read_all().map_err(PyErr::from)
                }
                GenoFormat::Ancestrymap => {
                    let reader = AncestryMapReader::new(prefix)?;
                    reader.read_all().map_err(PyErr::from)
                }
                GenoFormat::Vcf => {
                    let reader = open_vcf(prefix)?;
                    reader.read_all().map_err(PyErr::from)
                }
                GenoFormat::Tgeno => TgenoReader::new(prefix)
                    .and_then(|r| r.read_all())
                    .map_err(PyErr::from),
            }
        }

        let format_of = |prefix: &str, format: Option<&str>| match format {
            Some(f) => f.parse::<GenoFormat>().map_err(PyErr::from),
            None => GenoFormat::detect(prefix).map_err(PyErr::from),
        };
        let out_format = output_format.parse::<GenoFormat>()?;

        // Both inputs read in place and an output written SNP by SNP: stream,
        // without holding either file's genotypes in memory
        if written_streamed(out_format) {
            let in_a = mapped_input(file_a, format_of(file_a, input_format_a)?)?;
            let in_b = mapped_input(file_b, format_of(file_b, input_format_b)?)?;
            if let (Some((fa, samples_a, snps_a)), Some((fb, samples_b, snps_b))) = (in_a, in_b) {
                let (merged_snps, samples, stats) =
                    popmix_io::merge::plan_merge(&samples_a, &snps_a, &samples_b, &snps_b)?;
                let snps: Vec<Snp> = merged_snps.iter().map(|m| m.snp.clone()).collect();
                let mut rows = popmix_io::merge::MergedRows::new(&fa, &fb, &merged_snps);
                if by_samples(out_format) {
                    write_samples(out_format, output, &samples, &snps, &rows)?;
                } else {
                    write_rows(out_format, output, &samples, &snps, &mut rows)?;
                }
                return Ok((
                    snps_a.len(),
                    snps_b.len(),
                    samples_a.len(),
                    samples_b.len(),
                    snps.len(),
                    samples.len(),
                    stats,
                ));
            }
        }

        // Read both input files
        let geno_a = read_geno(file_a, input_format_a)?;
        let geno_b = read_geno(file_b, input_format_b)?;

        // Capture input stats before merging
        let snps_a = geno_a.n_snps();
        let snps_b = geno_b.n_snps();
        let samples_a = geno_a.n_samples();
        let samples_b = geno_b.n_samples();

        // Merge
        let (merged, stats) = merge_genotypes(&geno_a, &geno_b)?;

        // Write output
        let out_path = Path::new(output);

        match out_format {
            GenoFormat::Ancestrymap => {
                write_ancestrymap(out_path, &merged)?;
            }
            GenoFormat::Eigenstrat => {
                write_eigenstrat(out_path, &merged)?;
            }
            GenoFormat::PackedAncestrymap => {
                write_packedancestrymap(out_path, &merged)?;
            }
            GenoFormat::PackedPed => {
                write_packedped(out_path, &merged)?;
            }
            GenoFormat::Ped => {
                write_ped(out_path, &merged)?;
            }
            GenoFormat::Tgeno => {
                write_tgeno(out_path, &merged)?;
            }
            GenoFormat::Vcf => {
                return Err(input_error(
                "VCF output format is not supported. Use plink, eigenstrat, or packedancestrymap.",
            ));
            }
        }

        Ok((
            snps_a,
            snps_b,
            samples_a,
            samples_b,
            merged.n_snps(),
            merged.n_samples(),
            stats,
        ))
    };
    let (snps_a, snps_b, samples_a, samples_b, snps_merged, samples_merged, stats) =
        compute(py, work)?;

    // Return stats
    {
        let mut result: HashMap<String, Py<PyAny>> = HashMap::new();
        // Input file stats
        result.insert("snps_a".to_string(), snps_a.into_py_any(py)?);
        result.insert("snps_b".to_string(), snps_b.into_py_any(py)?);
        result.insert("samples_a".to_string(), samples_a.into_py_any(py)?);
        result.insert("samples_b".to_string(), samples_b.into_py_any(py)?);
        // Merge stats
        result.insert("snps_merged".to_string(), snps_merged.into_py_any(py)?);
        result.insert(
            "samples_merged".to_string(),
            samples_merged.into_py_any(py)?,
        );
        result.insert(
            "snps_identical".to_string(),
            stats.snps_identical.into_py_any(py)?,
        );
        result.insert(
            "flipped_snps".to_string(),
            stats.snps_flipped.into_py_any(py)?,
        );
        result.insert(
            "strand_ambiguous_snps".to_string(),
            stats.snps_strand_ambiguous.into_py_any(py)?,
        );
        result.insert(
            "snps_mismatch".to_string(),
            stats.snps_mismatch.into_py_any(py)?,
        );
        result.insert(
            "other_strand_snps".to_string(),
            stats.snps_other_strand.into_py_any(py)?,
        );
        result.insert(
            "duplicate_position_snps".to_string(),
            stats.snps_duplicate_position.into_py_any(py)?,
        );
        Ok(result)
    }
}
