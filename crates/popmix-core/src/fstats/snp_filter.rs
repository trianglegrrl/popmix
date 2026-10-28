//! Which SNPs go into f2 blocks, with the filters of admixtools2's
//! `extract_f2` (its documentation; the details below were measured by
//! running it): `maxmiss`, `minmaf`/`maxmaf`, `transitions`/`transversions`,
//! `keepsnps` and `poly_only`, plus a list of SNPs to leave out (TASKS.md
//! H1).
//!
//! - The minor allele frequency is that of the populations' alleles pooled:
//!   Σ pᵢnᵢ / Σ nᵢ over the populations with data, nᵢ their observed allele
//!   counts. SNPs below `min_maf` or above `max_maf` are dropped.
//! - Transitions are A<->G and C<->T; every other pair of bases is a
//!   transversion (A/T and C/G included). A SNP whose alleles are not both
//!   known bases is neither, and either filter drops it.
//! - `keep_snps` overrides `max_missing`, the MAF bounds and the
//!   transition/transversion choice (and, in `popmix_io::extract`, the
//!   autosome restriction), as admixtools2's `keepsnps` does; `poly_only`
//!   still applies, and so does `exclude_snps`.

use std::collections::HashSet;

use crate::genotype_matrix::AlleleFrequencies;
use crate::genotypes::Snp;

/// Which SNPs go into the f2 blocks, as in admixtools2 `extract_f2`.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct F2Options {
    /// Drop SNPs missing in more than this fraction of the populations
    /// (admixtools2 `maxmiss`). Default: 0, so every population must have
    /// data, and all f-statistics are computed on the same SNPs.
    pub max_missing: f64,
    /// Use only SNPs whose frequency differs between populations
    /// (admixtools2 `poly_only = "f2"`). Default: true.
    pub poly_only: bool,
    /// Drop SNPs with a pooled minor allele frequency below this
    /// (admixtools2 `minmaf`). Default: 0.
    pub min_maf: f64,
    /// Drop SNPs with a pooled minor allele frequency above this
    /// (admixtools2 `maxmaf`). Default: 0.5.
    pub max_maf: f64,
    /// Keep transitions (A<->G, C<->T); false for transversions only, as
    /// for ancient DNA without UDG treatment (admixtools2 `transitions`).
    pub transitions: bool,
    /// Keep transversions (admixtools2 `transversions`).
    pub transversions: bool,
    /// Only these SNP IDs, overriding the other filters but `poly_only`
    /// and `exclude_snps` (admixtools2 `keepsnps`).
    pub keep_snps: Option<HashSet<String>>,
    /// Never these SNP IDs.
    pub exclude_snps: Option<HashSet<String>>,
}

impl Default for F2Options {
    fn default() -> Self {
        F2Options {
            max_missing: 0.0,
            poly_only: true,
            min_maf: 0.0,
            max_maf: 0.5,
            transitions: true,
            transversions: true,
            keep_snps: None,
            exclude_snps: None,
        }
    }
}

/// Whether a SNP is polymorphic across populations: some two observed
/// frequencies differ (missing populations are ignored).
fn is_polymorphic(freqs: ndarray::ArrayView1<f64>) -> bool {
    let mut observed = freqs.iter().filter(|f| f.is_finite());
    match observed.next() {
        Some(first) => observed.any(|f| f != first),
        None => false,
    }
}

/// Whether a SNP's alleles are a transition (Some(true)), a transversion
/// (Some(false)), or not both known bases (None).
fn is_transition(snp: &Snp) -> Option<bool> {
    let base = |c: char| matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T');
    let (a, b) = (
        snp.allele_ref.to_ascii_uppercase(),
        snp.allele_alt.to_ascii_uppercase(),
    );
    if !base(a) || !base(b) {
        return None;
    }
    Some(matches!(
        (a, b),
        ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C')
    ))
}

impl F2Options {
    /// Whether the transition/transversion choice keeps `snp`.
    fn keeps_class(&self, snp: &Snp) -> bool {
        if self.transitions && self.transversions {
            return true;
        }
        match is_transition(snp) {
            Some(true) => self.transitions,
            Some(false) => self.transversions,
            None => false,
        }
    }

    /// Whether the MAF bounds keep row `i` of `af`.
    fn keeps_maf(&self, af: &AlleleFrequencies, i: usize) -> bool {
        if self.min_maf <= 0.0 && self.max_maf >= 0.5 {
            return true;
        }
        let (mut alleles, mut total) = (0.0, 0.0);
        for (&p, &n) in af.frequencies.row(i).iter().zip(af.counts.row(i)) {
            if n > 0 {
                alleles += p * f64::from(n);
                total += f64::from(n);
            }
        }
        let p = alleles / total;
        let maf = p.min(1.0 - p);
        maf >= self.min_maf && maf <= self.max_maf
    }
}

/// Rows of `af` kept by `opts`; SNPs not polymorphic across the
/// populations are dropped too if `poly_only` (passed separately: FST uses
/// monomorphic SNPs).
pub(crate) fn kept_snps(af: &AlleleFrequencies, opts: &F2Options, poly_only: bool) -> Vec<usize> {
    let n_pops = af.populations.len() as f64;
    (0..af.snps.len())
        .filter(|&i| {
            let snp = &af.snps[i];
            if opts
                .exclude_snps
                .as_ref()
                .is_some_and(|x| x.contains(&snp.id))
            {
                return false;
            }
            let selected = match &opts.keep_snps {
                Some(keep) => keep.contains(&snp.id),
                None => {
                    let missing = af.counts.row(i).iter().filter(|&&c| c == 0).count() as f64;
                    missing / n_pops <= opts.max_missing
                        && opts.keeps_class(snp)
                        && opts.keeps_maf(af, i)
                }
            };
            selected && (!poly_only || is_polymorphic(af.frequencies.row(i)))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::{Chromosome, PopulationId};
    use ndarray::array;

    /// Two populations, four SNPs: A/G, C/A, A/X and G/T.
    fn af() -> AlleleFrequencies {
        let snps = [('A', 'G'), ('C', 'A'), ('A', 'X'), ('G', 'T')]
            .iter()
            .enumerate()
            .map(|(i, &(r, a))| {
                Snp::new(
                    format!("s{i}"),
                    Chromosome::Autosome(1),
                    i as u32,
                    0.0,
                    r,
                    a,
                )
            })
            .collect();
        AlleleFrequencies::new(
            array![[0.1, 0.3], [0.5, 0.7], [0.2, 0.4], [0.9, 0.9]],
            array![[10, 30], [2, 2], [4, 4], [4, 0]],
            vec![PopulationId::from("A"), PopulationId::from("B")],
            snps,
        )
        .unwrap()
    }

    fn kept(opts: F2Options) -> Vec<usize> {
        kept_snps(
            &af(),
            &F2Options {
                max_missing: 1.0,
                ..opts
            },
            false,
        )
    }

    #[test]
    fn transversions_only_drops_transitions_and_unknown_bases() {
        let opts = F2Options {
            transitions: false,
            ..F2Options::default()
        };
        assert_eq!(kept(opts), vec![1, 3]);
    }

    #[test]
    fn transitions_only_drops_transversions_and_unknown_bases() {
        let opts = F2Options {
            transversions: false,
            ..F2Options::default()
        };
        assert_eq!(kept(opts), vec![0]);
    }

    #[test]
    fn maf_is_of_the_pooled_alleles() {
        // pooled: (0.1·10 + 0.3·30)/40 = 0.25; 0.6; 0.3; 0.9 (B missing), MAF 0.1
        let opts = F2Options {
            min_maf: 0.26,
            ..F2Options::default()
        };
        assert_eq!(kept(opts), vec![1, 2]);
        let opts = F2Options {
            max_maf: 0.25,
            ..F2Options::default()
        };
        assert_eq!(kept(opts), vec![0, 3]);
    }

    #[test]
    fn keep_overrides_the_filters_but_not_exclude() {
        let opts = F2Options {
            transitions: false,
            min_maf: 0.45,
            max_missing: 0.0,
            keep_snps: Some(["s0", "s3"].map(String::from).into()),
            exclude_snps: Some(["s3"].map(String::from).into()),
            ..F2Options::default()
        };
        assert_eq!(kept_snps(&af(), &opts, false), vec![0]);
    }
}
