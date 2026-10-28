//! f2 blocks from a genotype file of any format, as admixtools2's
//! `f2_from_geno`: the samples relabelled if asked, their ploidy found, and
//! their allele frequencies computed where the format allows without
//! loading the genotypes (a packed .geno, a PLINK .bed or a text EIGENSTRAT
//! .geno read in place; TGENO only for the samples used).

use std::collections::{HashMap, HashSet};
use std::ops::Range;
use std::path::Path;

use popmix_core::genotypes::{check_populations, check_sorted};
use popmix_core::warnings::check_genetic_map;
use popmix_core::{
    AlleleFrequencies, BlockDefinition, CoordinateEvidence, F2Blocks, F2Options, PopmixError,
    PopulationId, Result, Sample, Warning,
};

use crate::dataset::{open_mapped, read_genotypes, Mapped, PloidyRule};
use crate::tgeno::{is_tgeno_format, TgenoReader};

/// What [`extract_f2`] computes.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ExtractOptions {
    /// The populations whose f2 blocks are computed.
    pub populations: Vec<String>,
    /// Use only autosomal SNPs (chromosomes 1-22).
    pub auto_only: bool,
    /// The jackknife blocks.
    pub blocks: BlockDefinition,
    /// SNP filters (`maxmiss`, `poly_only`).
    pub f2: F2Options,
    /// admixtools2 `inds`/`pops`: sample ID to population. When given, only
    /// these samples are used, each in the population given.
    pub labels: Option<HashMap<String, String>>,
    /// How pseudo-haploid samples are found.
    pub ploidy: PloidyRule,
}

impl ExtractOptions {
    /// The defaults for `populations`: autosomes, 5 cM blocks, every
    /// population with data at each SNP, ploidy detected per sample.
    pub fn new(populations: Vec<String>) -> Self {
        ExtractOptions {
            populations,
            auto_only: true,
            blocks: BlockDefinition::default(),
            f2: F2Options::default(),
            labels: None,
            ploidy: PloidyRule::default(),
        }
    }
}

impl ExtractOptions {
    /// Whether only autosomal SNPs are read: `auto_only`, unless a list of
    /// SNPs to keep overrides it (admixtools2 `keepsnps`).
    pub(crate) fn autosomes_only(&self) -> bool {
        self.auto_only && self.f2.keep_snps.is_none()
    }
}

/// The population label of an excluded sample: no population is called this.
const EXCLUDED: &str = "\u{0}excluded";

/// admixtools2 `inds`/`pops` as sample labels: each listed sample's
/// population.
pub fn sample_labels(inds: Vec<String>, pops: Vec<String>) -> Result<HashMap<String, String>> {
    if inds.len() != pops.len() {
        return Err(PopmixError::parameter(
            "inds and pops",
            format!(
                "must have the same length, one population per sample: {} samples, {} \
                 populations",
                inds.len(),
                pops.len()
            ),
        ));
    }

    let n = inds.len();
    let labels: HashMap<String, String> = inds.into_iter().zip(pops).collect();
    if labels.len() < n {
        return Err(PopmixError::parameter("inds", "a sample is listed twice"));
    }
    Ok(labels)
}

/// Give the listed samples their populations and leave out the rest.
pub fn relabel_samples<S: std::hash::BuildHasher>(
    samples: &mut [Sample],
    labels: &HashMap<String, String, S>,
) -> Result<()> {
    let known: HashSet<&str> = samples.iter().map(|s| s.id.as_str()).collect();
    if let Some(missing) = labels.keys().find(|id| !known.contains(id.as_str())) {
        return Err(PopmixError::data(format!(
            "sample {missing:?} is not in the genotype data"
        )));
    }
    for s in samples.iter_mut() {
        s.population = PopulationId::from(labels.get(&s.id).map_or(EXCLUDED, |p| p.as_str()));
    }
    Ok(())
}

/// Relabel `samples` as `opts` asks and check that every population asked
/// for has samples.
fn prepare_samples(samples: &mut [Sample], opts: &ExtractOptions) -> Result<()> {
    if let Some(labels) = &opts.labels {
        relabel_samples(samples, labels)?;
    }
    let present: HashSet<&str> = samples.iter().map(|s| s.population.as_str()).collect();
    check_populations(
        opts.populations.iter().map(std::string::String::as_str),
        |p| present.contains(p),
    )
}

/// The f2 blocks of `opts.populations` in the dataset at `prefix`.
pub fn extract_f2(
    prefix: &str,
    opts: &ExtractOptions,
    warn: &mut dyn FnMut(Warning),
) -> Result<F2Blocks> {
    check_options(opts)?;
    let pops: Vec<&str> = opts
        .populations
        .iter()
        .map(std::string::String::as_str)
        .collect();
    if let Some(mapped) = open_mapped(prefix, warn)? {
        return from_mapped(prefix, mapped, opts, &pops, warn);
    }
    in_memory(prefix, opts, &pops, warn)
}

/// The SNP filters' bounds checked.
pub(crate) fn check_options(opts: &ExtractOptions) -> Result<()> {
    for (name, value, range) in [
        ("minmaf", opts.f2.min_maf, 0.0..=0.5),
        ("maxmaf", opts.f2.max_maf, 0.0..=0.5),
    ] {
        if !range.contains(&value) {
            return Err(PopmixError::parameter(
                name,
                format!("must be between 0 and 0.5, got {value}"),
            ));
        }
    }
    if !(0.0..=1.0).contains(&opts.f2.max_missing) {
        return Err(PopmixError::parameter(
            "maxmiss",
            format!("must be between 0 and 1, got {}", opts.f2.max_missing),
        ));
    }
    Ok(())
}

/// The f2 blocks of a dataset that cannot be read in place: TGENO (the
/// requested samples' records) or every genotype of the other formats.
pub(crate) fn in_memory(
    prefix: &str,
    opts: &ExtractOptions,
    pops: &[&str],
    warn: &mut dyn FnMut(Warning),
) -> Result<F2Blocks> {
    let geno_path = format!("{prefix}.geno");
    let mut geno = if Path::new(&geno_path).exists() && is_tgeno_format(&geno_path)? {
        // TGENO stores each sample's genotypes together: read only the
        // records of the requested populations' samples
        let reader = TgenoReader::new(prefix)?;
        let mut samples = reader.samples().to_vec();
        prepare_samples(&mut samples, opts)?;
        let selected: Vec<usize> = (0..samples.len())
            .filter(|&i| pops.contains(&samples[i].population.as_str()))
            .collect();
        let mut geno = reader.read_samples(&selected)?;
        for (sample, &i) in geno.samples_mut().iter_mut().zip(&selected) {
            *sample = samples[i].clone();
        }
        geno
    } else {
        // the other formats (PED, ANCESTRYMAP): every genotype loaded
        let mut geno = read_genotypes(prefix, warn)?;
        prepare_samples(geno.samples_mut(), opts)?;
        geno
    };
    check_sorted(geno.snps())?;
    opts.ploidy.apply(&mut geno);
    let panel = crate::provenance::Panel::of(geno.snps());
    let geno = if opts.autosomes_only() {
        geno.filter_autosomal()?
    } else {
        geno
    };
    let mut evidence = CoordinateEvidence::default();
    evidence.observe(geno.snps());
    check_genetic_map(&evidence, &opts.blocks, warn);
    let af = AlleleFrequencies::from_genotypes_for_populations(&geno, Some(pops));
    let blocks = F2Blocks::from_allele_frequencies_with_options(&af, &opts.blocks, &opts.f2);
    let provenance =
        crate::provenance::extraction(prefix, opts, geno.samples(), &panel, snps_used(&blocks));
    Ok(blocks.with_provenance(Some(provenance)))
}

/// The f2 blocks of a dataset read in place, a chromosome at a time: only
/// one chromosome's frequency table is held (blocks never span
/// chromosomes).
fn from_mapped(
    prefix: &str,
    mapped: Mapped,
    opts: &ExtractOptions,
    pops: &[&str],
    warn: &mut dyn FnMut(Warning),
) -> Result<F2Blocks> {
    let MappedRuns {
        file,
        samples,
        snps,
        runs,
        evidence,
    } = MappedRuns::prepare(mapped, opts, pops, warn)?;
    let auto_only = opts.autosomes_only();
    let parts = runs
        .iter()
        .map(|r| file.allele_frequencies_at(&samples, &snps[r.clone()], r.start, pops, auto_only));
    let blocks = match F2Blocks::from_parts(parts, &evidence, &opts.blocks, &opts.f2) {
        Some(blocks) => blocks,
        // the coordinates a block definition needs were only known from
        // all chromosomes together
        None => {
            let af = file.allele_frequencies(&samples, &snps, pops, auto_only);
            F2Blocks::from_allele_frequencies_with_options(&af, &opts.blocks, &opts.f2)
        }
    };
    let provenance = crate::provenance::extraction(
        prefix,
        opts,
        &samples,
        &crate::provenance::Panel::of(&snps),
        snps_used(&blocks),
    );
    Ok(blocks.with_provenance(Some(provenance)))
}

/// The SNPs `blocks` hold: their block lengths' sum.
fn snps_used(blocks: &F2Blocks) -> usize {
    blocks.block_lengths().iter().sum()
}

/// A dataset read in place, ready to be read a chromosome at a time: its
/// samples relabelled and with their ploidy, its chromosomes' runs of
/// SNPs, and the evidence (map or base pairs) of the SNPs used.
pub(crate) struct MappedRuns {
    pub(crate) file: crate::SnpMajor,
    pub(crate) samples: Vec<popmix_core::Sample>,
    pub(crate) snps: Vec<popmix_core::Snp>,
    pub(crate) runs: Vec<Range<usize>>,
    pub(crate) evidence: CoordinateEvidence,
}

impl MappedRuns {
    pub(crate) fn prepare(
        mapped: Mapped,
        opts: &ExtractOptions,
        pops: &[&str],
        warn: &mut dyn FnMut(Warning),
    ) -> Result<Self> {
        let Mapped {
            file,
            mut samples,
            snps,
        } = mapped;
        check_sorted(&snps)?;
        prepare_samples(&mut samples, opts)?;
        opts.ploidy.apply_mapped(&file, &mut samples, pops);
        let mut runs: Vec<Range<usize>> = Vec::new();
        for (i, snp) in snps.iter().enumerate() {
            match runs.last_mut() {
                Some(r) if snps[r.start].chromosome == snp.chromosome => r.end = i + 1,
                _ => runs.push(i..i + 1),
            }
        }
        let mut evidence = CoordinateEvidence::default();
        let auto_only = opts.autosomes_only();
        for snp in snps.iter().filter(|s| !auto_only || s.is_autosomal()) {
            evidence.observe(std::slice::from_ref(snp));
        }
        check_genetic_map(&evidence, &opts.blocks, warn);
        Ok(MappedRuns {
            file,
            samples,
            snps,
            runs,
            evidence,
        })
    }
}
