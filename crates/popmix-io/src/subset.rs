//! Choosing and relabelling samples and SNPs of a dataset (TASKS.md H2,
//! P3): for `convertf`, and the samples for f2 extraction.
//!
//! Samples are taken in this order: kept or excluded by ID, relabelled by
//! ID (a new population for a sample), relabelled by population (a new name
//! for every sample of one), then chosen by their final population. SNPs
//! are kept or excluded by ID, and chosen by chromosome and by a region of
//! one chromosome. Every ID and population named must be in the dataset, so
//! a typo is an error rather than a silent no-op.

use std::collections::{HashMap, HashSet};

use popmix_core::genotypes::Chromosome;
use popmix_core::{GenotypeMatrix, PopmixError, PopulationId, Result, Sample, Snp};

/// What to keep of a dataset, and how to relabel it. The default keeps
/// everything as it is.
#[derive(Debug, Clone, Default)]
#[non_exhaustive]
pub struct Subset {
    /// Only these sample IDs.
    pub keep_samples: Option<HashSet<String>>,
    /// Never these sample IDs.
    pub exclude_samples: Option<HashSet<String>>,
    /// A sample ID's new population.
    pub relabel_samples: HashMap<String, String>,
    /// A population's new name.
    pub relabel_populations: HashMap<String, String>,
    /// Only samples of these populations, by their final labels.
    pub populations: Option<Vec<String>>,
    /// Only these SNP IDs.
    pub keep_snps: Option<HashSet<String>>,
    /// Never these SNP IDs.
    pub exclude_snps: Option<HashSet<String>>,
    /// Only SNPs on these chromosomes.
    pub chromosomes: Option<HashSet<Chromosome>>,
    /// Only SNPs on this chromosome from this to this base pair, both
    /// included.
    pub region: Option<(Chromosome, u32, u32)>,
}

fn check_known<'a>(
    what: &str,
    names: impl IntoIterator<Item = &'a String>,
    known: &HashSet<&str>,
) -> Result<()> {
    let mut missing: Vec<&str> = names
        .into_iter()
        .map(String::as_str)
        .filter(|n| !known.contains(n))
        .collect();
    if missing.is_empty() {
        return Ok(());
    }
    missing.sort_unstable();
    let shown = missing
        .iter()
        .take(5)
        .copied()
        .collect::<Vec<_>>()
        .join(", ");
    Err(PopmixError::parameter(
        what,
        format!(
            "{} not in the dataset: {shown}{}",
            missing.len(),
            if missing.len() > 5 { ", ..." } else { "" }
        ),
    ))
}

impl Subset {
    /// Whether the subset keeps every sample as labelled.
    pub fn keeps_all_samples(&self) -> bool {
        self.keep_samples.is_none()
            && self.exclude_samples.is_none()
            && self.relabel_samples.is_empty()
            && self.relabel_populations.is_empty()
            && self.populations.is_none()
    }

    /// Whether the subset keeps every SNP.
    pub fn keeps_all_snps(&self) -> bool {
        self.keep_snps.is_none()
            && self.exclude_snps.is_none()
            && self.chromosomes.is_none()
            && self.region.is_none()
    }

    /// The samples kept, as indices into `samples` in their order, and as
    /// relabelled.
    pub fn samples(&self, samples: &[Sample]) -> Result<(Vec<usize>, Vec<Sample>)> {
        let ids: HashSet<&str> = samples.iter().map(|s| s.id.as_str()).collect();
        // populations as relabelled by ID, which relabel_populations renames
        let pops: HashSet<&str> = samples
            .iter()
            .map(|s| {
                self.relabel_samples
                    .get(&s.id)
                    .map_or(s.population.as_str(), String::as_str)
            })
            .collect();
        for (what, set) in [
            ("keep samples", &self.keep_samples),
            ("exclude samples", &self.exclude_samples),
        ] {
            if let Some(set) = set {
                check_known(what, set, &ids)?;
            }
        }
        check_known("relabel samples", self.relabel_samples.keys(), &ids)?;
        check_known(
            "relabel populations",
            self.relabel_populations.keys(),
            &pops,
        )?;

        let mut kept = Vec::new();
        let mut out = Vec::new();
        for (i, s) in samples.iter().enumerate() {
            if self
                .keep_samples
                .as_ref()
                .is_some_and(|k| !k.contains(&s.id))
                || self
                    .exclude_samples
                    .as_ref()
                    .is_some_and(|x| x.contains(&s.id))
            {
                continue;
            }
            let mut pop = self
                .relabel_samples
                .get(&s.id)
                .cloned()
                .unwrap_or_else(|| s.population.as_str().to_string());
            if let Some(new) = self.relabel_populations.get(&pop) {
                pop = new.clone();
            }
            let mut sample = s.clone();
            sample.population = PopulationId::from(pop.as_str());
            kept.push(i);
            out.push(sample);
        }
        if let Some(wanted) = &self.populations {
            let present: HashSet<&str> = out.iter().map(|s| s.population.as_str()).collect();
            check_known("populations", wanted, &present)?;
            let (k, o): (Vec<usize>, Vec<Sample>) = kept
                .into_iter()
                .zip(out)
                .filter(|(_, s)| wanted.iter().any(|p| p == s.population.as_str()))
                .unzip();
            kept = k;
            out = o;
        }
        if out.is_empty() {
            return Err(PopmixError::data("the sample selection keeps no samples"));
        }
        Ok((kept, out))
    }

    /// Which of `snps` are kept, one entry per SNP.
    pub fn snps(&self, snps: &[Snp]) -> Result<Vec<bool>> {
        let ids: HashSet<&str> = snps.iter().map(|s| s.id.as_str()).collect();
        for (what, set) in [
            ("keep SNPs", &self.keep_snps),
            ("exclude SNPs", &self.exclude_snps),
        ] {
            if let Some(set) = set {
                check_known(what, set, &ids)?;
            }
        }
        let keep: Vec<bool> = snps
            .iter()
            .map(|s| {
                self.keep_snps.as_ref().is_none_or(|k| k.contains(&s.id))
                    && !self
                        .exclude_snps
                        .as_ref()
                        .is_some_and(|x| x.contains(&s.id))
                    && self
                        .chromosomes
                        .as_ref()
                        .is_none_or(|c| c.contains(&s.chromosome))
                    && self.region.as_ref().is_none_or(|(c, from, to)| {
                        s.chromosome == *c && (*from..=*to).contains(&s.position_bp)
                    })
            })
            .collect();
        if !keep.iter().any(|&k| k) {
            return Err(PopmixError::data("the SNP selection keeps no SNPs"));
        }
        Ok(keep)
    }
}

/// The genotypes of `geno` at the SNPs `keep` marks, for the samples
/// `kept` (indices, in order), labelled as `samples`.
pub fn select(
    geno: &GenotypeMatrix,
    kept: &[usize],
    samples: Vec<Sample>,
    keep: &[bool],
) -> Result<GenotypeMatrix> {
    let rows: Vec<usize> = (0..keep.len()).filter(|&i| keep[i]).collect();
    let data = geno.data();
    let out =
        ndarray::Array2::from_shape_fn((rows.len(), kept.len()), |(r, c)| data[[rows[r], kept[c]]]);
    let snps = rows.iter().map(|&i| geno.snps()[i].clone()).collect();
    GenotypeMatrix::new(out, samples, snps)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn samples() -> Vec<Sample> {
        [("a", "P"), ("b", "P"), ("c", "Q"), ("d", "R")]
            .iter()
            .map(|(i, p)| Sample::new(*i, *p))
            .collect()
    }

    fn set(xs: &[&str]) -> HashSet<String> {
        xs.iter().map(ToString::to_string).collect()
    }

    #[test]
    fn samples_are_kept_relabelled_then_chosen_by_population() {
        let subset = Subset {
            exclude_samples: Some(set(&["b"])),
            relabel_samples: [("c".to_string(), "P".to_string())].into(),
            relabel_populations: [("P".to_string(), "New".to_string())].into(),
            populations: Some(vec!["New".to_string()]),
            ..Subset::default()
        };
        let (kept, out) = subset.samples(&samples()).unwrap();
        assert_eq!(kept, vec![0, 2]);
        assert!(out.iter().all(|s| s.population.as_str() == "New"));
    }

    #[test]
    fn an_unknown_id_is_an_error() {
        let subset = Subset {
            keep_samples: Some(set(&["a", "zz"])),
            ..Subset::default()
        };
        let err = subset.samples(&samples()).unwrap_err().to_string();
        assert!(err.contains("zz"), "{err}");
    }

    #[test]
    fn snps_by_region_are_inclusive() {
        let snps: Vec<Snp> = (0..5)
            .map(|i| {
                Snp::new(
                    format!("s{i}"),
                    Chromosome::Autosome(2),
                    100 * i,
                    0.0,
                    'A',
                    'G',
                )
            })
            .collect();
        let subset = Subset {
            region: Some((Chromosome::Autosome(2), 100, 300)),
            ..Subset::default()
        };
        assert_eq!(
            subset.snps(&snps).unwrap(),
            vec![false, true, true, true, false]
        );
    }
}
