//! File conversion and merging for the bindings: `convertf` and `merge`,
//! streamed from inputs read in place where the formats allow, and the
//! mapped-input helpers the other readers share.

use popmix_io::sample_subset::{Source, SubsetSamples};
use popmix_io::subset::Subset;
use std::collections::HashMap;

use popmix_core::python::internal_error;
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;

use popmix_core::{GenotypeMatrix, Sample, Snp};
use popmix_io::tgeno::TgenoReader;
use popmix_io::{
    write_ancestrymap, write_eigenstrat, write_packedancestrymap, write_packedped, write_ped,
    write_tgeno, AncestryMapReader, EigenstratReader, GenoFormat, PackedAncestryMapReader,
    PedReader, PlinkReader,
};

use crate::compute;
use crate::warnings::{open_vcf, with_warnings};

/// Convert genotype data between file formats.
///
/// When the input can be read in place (packed .geno, PLINK .bed, text
/// EIGENSTRAT, VCF) the output is written without holding every genotype
/// in memory; otherwise the input is read whole first.
///
/// Parameters
/// ----------
/// input_prefix
///     Path prefix of the input dataset (the full path of a VCF file).
/// output_prefix
///     Path prefix of the output dataset.
/// output_format
///     Output format: "ancestrymap", "eigenstrat", "packedancestrymap", "packedped" (or
///     "plink", "bed"), "ped" (or "pedmap") or "tgeno" (or
///     "transpose_packed"), in any case.
///     "vcf" raises ``InputError``.
/// input_format
///     Input format, by the same names or "vcf"; None (the default) to
///     detect it from the files.
/// keep_inds, exclude_inds
///     Only these sample IDs; never these.
/// relabel
///     A sample ID's new population (dict).
/// relabel_pops
///     A population's new name (dict).
/// populations
///     Only samples of these populations, as relabelled.
/// keep_snps, exclude_snps
///     Only these SNP IDs; never these.
/// chromosomes
///     Only SNPs on these chromosomes ("1"-"22", "X", "Y", "MT").
/// region
///     Only SNPs in ``(chromosome, start, end)``, base pairs, both ends
///     included.
///
/// Samples are kept or excluded by ID, relabelled by ID, then by
/// population, then chosen by population; an ID or population named that
/// is not in the dataset raises ``InputError``.
///
/// Returns
/// -------
/// dict
///     ``n_samples`` and ``n_snps`` (int) written, and ``input_format`` and
///     ``output_format`` (str, the canonical names, such as "PACKEDPED").
#[pyfunction]
#[pyo3(signature = (input_prefix, output_prefix, output_format, *, input_format=None, keep_inds=None,
                    exclude_inds=None, relabel=None, relabel_pops=None, populations=None, keep_snps=None,
                    exclude_snps=None, chromosomes=None, region=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn convertf(
    py: Python<'_>,
    input_prefix: &str,
    output_prefix: &str,
    output_format: &str,
    input_format: Option<&str>,
    keep_inds: Option<Vec<String>>,
    exclude_inds: Option<Vec<String>>,
    relabel: Option<HashMap<String, String>>,
    relabel_pops: Option<HashMap<String, String>>,
    populations: Option<Vec<String>>,
    keep_snps: Option<Vec<String>>,
    exclude_snps: Option<Vec<String>>,
    chromosomes: Option<Vec<String>>,
    region: Option<(String, u32, u32)>,
) -> PyResult<HashMap<String, Py<PyAny>>> {
    let subset = crate::subset::subset(
        keep_inds,
        exclude_inds,
        relabel,
        relabel_pops,
        populations,
        keep_snps,
        exclude_snps,
        chromosomes,
        region,
    );
    let (n_samples, n_snps, in_format, out_format) = compute(py, || {
        convertf_files(
            input_prefix,
            output_prefix,
            output_format,
            input_format,
            &subset,
        )
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
            | GenoFormat::Vcf
    )
}

/// Write a genotype matrix as VCF (`popmix_io::vcf_writer`).
fn write_vcf_matrix(output: &str, geno: &GenotypeMatrix) -> PyResult<()> {
    let mut rows = popmix_core::MatrixRows::new(
        geno,
        (0..geno.n_snps()).collect(),
        (0..geno.n_samples()).collect(),
    );
    popmix_io::vcf_writer::write_vcf_rows(output, geno.samples(), geno.snps(), &mut rows)
        .map_err(PyErr::from)
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
        GenoFormat::Vcf => {
            popmix_io::vcf_writer::write_vcf_rows(output_prefix, samples, snps, rows)
        }
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
    subset: &Subset,
) -> PyResult<Option<(usize, usize)>> {
    if !written_streamed(out_format) {
        return Ok(None);
    }
    if in_format == GenoFormat::Tgeno {
        return convertf_tgeno(
            output_prefix,
            TgenoReader::new(input_prefix)?,
            out_format,
            subset,
        )
        .map(Some);
    }
    let Some((file, samples, snps)) = mapped_input(input_prefix, in_format)? else {
        return Ok(None);
    };
    let (kept, samples) = subset.samples(&samples)?;
    let keep = subset.snps(&snps)?;
    let snps = kept_snps(&snps, &keep);
    if by_samples(out_format) {
        let part = SubsetSamples {
            source: Source::SnpMajor(&file),
            samples: kept,
            snps: kept_indices(&keep),
        };
        write_samples(out_format, output_prefix, &samples, &snps, &part)?;
        return Ok(Some((samples.len(), snps.len())));
    }
    let mut rows = popmix_io::SnpMajorRows::from_file(file, keep, kept)?;
    write_rows(out_format, output_prefix, &samples, &snps, &mut rows)?;
    Ok(Some((samples.len(), snps.len())))
}

/// The SNPs `keep` marks.
fn kept_snps(snps: &[Snp], keep: &[bool]) -> Vec<Snp> {
    snps.iter()
        .zip(keep)
        .filter(|(_, &k)| k)
        .map(|(s, _)| s.clone())
        .collect()
}

/// The indices of the SNPs `keep` marks.
fn kept_indices(keep: &[bool]) -> Vec<usize> {
    keep.iter()
        .enumerate()
        .filter(|(_, &k)| k)
        .map(|(i, _)| i)
        .collect()
}

/// [`convertf_streamed`] from TGENO, read in place a part at a time
/// (`popmix_io::tgeno_stream`): the sample and SNP counts.
fn convertf_tgeno(
    output_prefix: &str,
    reader: TgenoReader,
    out_format: GenoFormat,
    subset: &Subset,
) -> PyResult<(usize, usize)> {
    let (kept, samples) = subset.samples(reader.samples())?;
    let keep = subset.snps(reader.snps())?;
    let snps = kept_snps(reader.snps(), &keep);
    if by_samples(out_format) {
        let part = SubsetSamples {
            source: Source::Tgeno(&reader),
            samples: kept,
            snps: kept_indices(&keep),
        };
        write_samples(out_format, output_prefix, &samples, &snps, &part)?;
        return Ok((samples.len(), snps.len()));
    }
    let mut rows = popmix_io::TgenoRows::new(reader, keep, kept)?;
    write_rows(out_format, output_prefix, &samples, &snps, &mut rows)?;
    Ok((samples.len(), snps.len()))
}

/// `convertf`'s work: the sample and SNP counts and both formats.
fn convertf_files(
    input_prefix: &str,
    output_prefix: &str,
    output_format: &str,
    input_format: Option<&str>,
    subset: &Subset,
) -> PyResult<(usize, usize, GenoFormat, GenoFormat)> {
    use std::path::Path;

    // Detect or parse input format
    let in_format = match input_format {
        Some(fmt) => fmt.parse::<GenoFormat>()?,
        None => GenoFormat::detect(input_prefix)?,
    };

    // Parse output format
    let out_format = output_format.parse::<GenoFormat>()?;

    // From a file read in place (packed, PLINK .bed, text EIGENSTRAT,
    // TGENO, VCF) to any format: stream, without holding the genotypes in
    // memory. PED and ANCESTRYMAP input are read whole, below.
    if let Some(done) =
        convertf_streamed(input_prefix, output_prefix, in_format, out_format, subset)?
    {
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
    let geno = if subset.keeps_all_samples() && subset.keeps_all_snps() {
        geno
    } else {
        let (kept, samples) = subset.samples(geno.samples())?;
        let keep = subset.snps(geno.snps())?;
        popmix_io::subset::select(&geno, &kept, samples, &keep)?
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
        GenoFormat::Vcf => write_vcf_matrix(output_prefix, &geno)?,
    }

    Ok((n_samples, n_snps, in_format, out_format))
}

/// Merge two genotype datasets.
///
/// Keeps the SNPs in both files. SNPs are matched by chromosome and
/// position (and by name where several share one). Alleles are reconciled
/// as AdmixTools mergeit does: swapped alleles are recoded, a SNP on the
/// other strand is complemented, the unknown allele of a monomorphic SNP
/// is filled in, and A/T, C/G and third-allele SNPs are left out and
/// counted. The samples of A come first, then those of B.
///
/// Parameters
/// ----------
/// file_a
///     Path prefix of the first dataset.
/// file_b
///     Path prefix of the second dataset.
/// output
///     Path prefix of the merged dataset.
/// output_format
///     Output format: "ancestrymap", "eigenstrat", "packedancestrymap", "packedped" (or
///     "plink", "bed"), "ped" (or "pedmap") or "tgeno" (or
///     "transpose_packed"), in any case.
///     "vcf" raises ``InputError``.
/// input_format_a, input_format_b
///     Format of ``file_a`` / ``file_b``, by the same names or "vcf"; None
///     (the default) to detect it from the files.
/// strandcheck
///     True (the default) to reconcile strands: complemented alleles are
///     the same SNP on the other strand, and A/T and C/G SNPs, whose strand
///     the alleles cannot tell, are left out. False when both files are on
///     one strand (two releases of one panel), as mergeit ``strandcheck:
///     NO``: A/T and C/G SNPs are kept, and complemented alleles are a
///     mismatch.
///
/// Returns
/// -------
/// dict
///     Counts (int): ``snps_a``, ``snps_b``, ``samples_a``, ``samples_b``
///     (the inputs); ``snps_merged``, ``samples_merged`` (the output);
///     ``snps_identical`` (same alleles), ``flipped_snps`` (swapped alleles,
///     recoded), ``strand_ambiguous_snps`` (A/T or C/G, left out),
///     ``snps_mismatch`` (other alleles, left out), ``other_strand_snps``
///     (kept, complemented from file B) and ``duplicate_position_snps``
///     (SNPs of A left out because several SNPs of B share their position
///     and none has their name).
#[pyfunction]
#[pyo3(signature = (file_a, file_b, output, output_format, *, input_format_a=None, input_format_b=None, strandcheck=true))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn merge(
    py: Python<'_>,
    file_a: &str,
    file_b: &str,
    output: &str,
    output_format: &str,
    input_format_a: Option<&str>,
    input_format_b: Option<&str>,
    strandcheck: bool,
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
                let (merged_snps, samples, stats) = popmix_io::merge::plan_merge_on(
                    &samples_a,
                    &snps_a,
                    &samples_b,
                    &snps_b,
                    strandcheck,
                )?;
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
        let (merged, stats) = popmix_io::merge::merge_genotypes_on(&geno_a, &geno_b, strandcheck)?;

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
            GenoFormat::Vcf => write_vcf_matrix(output, &merged)?,
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
    crate::warnings::emit(|w| {
        popmix_core::warnings::check_panel_overlap(snps_a, snps_b, snps_a - stats.snps_only_a, w);
    })?;

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
