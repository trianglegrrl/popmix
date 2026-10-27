//! `extract_f2`: f2 blocks computed from a genotype file, with population
//! relabelling and ploidy detection.

use crate::f2::PyF2Blocks;
use crate::genotypes::read_genotypes_internal;
use crate::ploidy;
use crate::{compute, require_sample_populations};
use popmix_core::python::input_error;
use popmix_core::{AlleleFrequencies, BlockDefinition, F2Blocks, Sample};
use popmix_io::tgeno::{is_tgeno_format, TgenoReader};
use pyo3::prelude::*;
use std::collections::HashMap;

/// Extract F2 blocks from genotype data.
///
/// This is the main entry point for computing f-statistics.
/// Uses a streaming approach (like admixtools2) to avoid loading all genotype
/// data into memory. Reads SNPs one at a time and computes allele frequencies
/// incrementally.
///
/// Args:
///     prefix: Path prefix for the genotype dataset
///     populations: List of populations to compute F2 for (REQUIRED)
///     auto_only: If True (default), only use autosomal SNPs (chr 1-22)
///     block_size: Block size in Morgans (default: 0.05 = 5 cM)
///     maxmiss: Drop SNPs missing in more than this fraction of populations,
///         as admixtools2 ``maxmiss`` (default: 0, every population needs data)
///     poly_only: Use only SNPs whose frequency differs between populations,
///         as admixtools2 ``poly_only = "f2"`` (default: True)
///     inds, pops: Use only these samples, each in the population given in
///         ``pops`` (same length) instead of its ``.ind`` label, as
///         admixtools2 ``f2_from_geno(inds, pops)``
///     pseudohaploid: Which samples count one allele per SNP: "detect"
///         (default: those with no heterozygous call at any SNP),
///         "admixtools2" (no heterozygous call among the first 1,000 SNPs,
///         which also catches some diploid genomes; ADMIXTOOLS2-0007) or
///         "diploid" (none)
///
/// Returns:
///     PyF2Blocks object for computing f2, f3, f4 statistics
#[pyfunction]
#[pyo3(signature = (prefix, populations, auto_only=true, block_size=None, maxmiss=0.0, poly_only=true, inds=None, pops=None, pseudohaploid="detect"))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn extract_f2(
    py: Python<'_>,
    prefix: &str,
    populations: Vec<String>,
    auto_only: bool,
    block_size: Option<f64>,
    maxmiss: f64,
    poly_only: bool,
    inds: Option<Vec<String>>,
    pops: Option<Vec<String>>,
    pseudohaploid: &str,
) -> PyResult<PyF2Blocks> {
    compute(py, move || {
        use std::path::Path;

        let relabel = sample_labels(inds, pops)?;
        let rule = ploidy::ploidy_rule(pseudohaploid)?;

        if !(0.0..=1.0).contains(&maxmiss) {
            return Err(input_error(format!(
                "maxmiss must be between 0 and 1, got {maxmiss}"
            )));
        }
        let f2_options = popmix_core::F2Options {
            max_missing: maxmiss,
            poly_only,
        };

        // Create block definition
        let block_def = match block_size {
            Some(size) => BlockDefinition::new(size),
            None => BlockDefinition::default(),
        };

        // Packed .geno and PLINK .bed are read in place, SNP by SNP
        let geno_path = format!("{}.geno", prefix);
        if let Some(ploidy::Mapped {
            file,
            mut samples,
            snps,
        }) = ploidy::mapped(prefix)?
        {
            popmix_core::genotypes::check_sorted(&snps)?;
            if let Some(labels) = &relabel {
                apply_sample_labels(&mut samples, labels)?;
            }
            require_sample_populations(&samples, &populations)?;
            let pop_refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
            rule.apply_mapped(&file, &mut samples, &pop_refs);
            // a chromosome at a time: only one chromosome's frequency table
            // is held (blocks never span chromosomes)
            let mut runs: Vec<std::ops::Range<usize>> = Vec::new();
            for (i, snp) in snps.iter().enumerate() {
                match runs.last_mut() {
                    Some(r) if snps[r.start].chromosome == snp.chromosome => r.end = i + 1,
                    _ => runs.push(i..i + 1),
                }
            }
            let mut evidence = popmix_core::CoordinateEvidence::default();
            for snp in snps.iter().filter(|s| !auto_only || s.is_autosomal()) {
                evidence.observe(std::slice::from_ref(snp));
            }
            let parts = runs.iter().map(|r| {
                file.allele_frequencies_at(
                    &samples,
                    &snps[r.clone()],
                    r.start,
                    &pop_refs,
                    auto_only,
                )
            });
            let f2_blocks = match F2Blocks::from_parts(parts, &evidence, &block_def, &f2_options) {
                Some(blocks) => blocks,
                None => {
                    let af = file.allele_frequencies(&samples, &snps, &pop_refs, auto_only);
                    F2Blocks::from_allele_frequencies_with_options(&af, &block_def, &f2_options)
                }
            };
            return Ok(PyF2Blocks { inner: f2_blocks });
        }

        let mut geno = if Path::new(&geno_path).exists() && is_tgeno_format(&geno_path)? {
            // TGENO stores each sample's genotypes together: read only the
            // records of the requested populations' samples
            let reader = TgenoReader::new(prefix)?;
            let mut samples = reader.samples().to_vec();
            if let Some(labels) = &relabel {
                apply_sample_labels(&mut samples, labels)?;
            }
            require_sample_populations(&samples, &populations)?;
            let selected: Vec<usize> = (0..samples.len())
                .filter(|&i| {
                    populations
                        .iter()
                        .any(|p| p == samples[i].population.as_str())
                })
                .collect();
            let mut geno = reader.read_samples(&selected)?;
            for (sample, &i) in geno.samples_mut().iter_mut().zip(&selected) {
                *sample = samples[i].clone();
            }
            geno
        } else {
            // Fallback: load all genotypes for other formats (PLINK, text EIGENSTRAT)
            let mut geno = read_genotypes_internal(prefix)?;
            if let Some(labels) = &relabel {
                apply_sample_labels(geno.samples_mut(), labels)?;
            }
            require_sample_populations(geno.samples(), &populations)?;
            geno
        };
        popmix_core::genotypes::check_sorted(geno.snps())?;
        rule.apply(&mut geno);

        // Filter to autosomal if requested
        let geno = if auto_only {
            geno.filter_autosomal()?
        } else {
            geno
        };

        // Compute allele frequencies for specified populations
        let pop_refs: Vec<&str> = populations.iter().map(|s| s.as_str()).collect();
        let af = AlleleFrequencies::from_genotypes_for_populations(&geno, Some(&pop_refs));

        // Compute F2 blocks
        let f2_blocks =
            F2Blocks::from_allele_frequencies_with_options(&af, &block_def, &f2_options);

        Ok(PyF2Blocks { inner: f2_blocks })
    })
}

/// The population label of an excluded sample: no population is called this.
pub(crate) const EXCLUDED: &str = "\u{0}excluded";

/// admixtools2 `inds`/`pops`: each listed sample's population.
pub(crate) fn sample_labels(
    inds: Option<Vec<String>>,
    pops: Option<Vec<String>>,
) -> PyResult<Option<HashMap<String, String>>> {
    match (inds, pops) {
        (None, None) => Ok(None),
        (Some(inds), Some(pops)) if inds.len() == pops.len() => {
            let n = inds.len();
            let labels: HashMap<String, String> = inds.into_iter().zip(pops).collect();
            if labels.len() < n {
                return Err(input_error("a sample is listed twice in inds"));
            }
            Ok(Some(labels))
        }
        (Some(_), Some(_)) => Err(input_error(
            "inds and pops must have the same length: one population per sample",
        )),
        _ => Err(input_error("give inds and pops together")),
    }
}

/// Give the listed samples their populations and leave out the rest.
pub(crate) fn apply_sample_labels(
    samples: &mut [Sample],
    labels: &HashMap<String, String>,
) -> PyResult<()> {
    let known: std::collections::HashSet<&str> = samples.iter().map(|s| s.id.as_str()).collect();
    if let Some(missing) = labels.keys().find(|id| !known.contains(id.as_str())) {
        return Err(input_error(format!(
            "sample {missing:?} is not in the genotype data"
        )));
    }
    for s in samples.iter_mut() {
        s.population =
            popmix_core::PopulationId::from(labels.get(&s.id).map_or(EXCLUDED, |p| p.as_str()));
    }
    Ok(())
}
