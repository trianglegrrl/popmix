//! How samples are found to be pseudo-haploid (counted as one allele per
//! SNP), shared by every command that reads genotypes.

use pyo3::prelude::*;

use popmix_core::{GenotypeMatrix, Ploidy, Sample};
use popmix_io::detect_ploidy_packedancestrymap_selected;

use crate::to_py_err;

/// The `pseudohaploid` option.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PloidyRule {
    /// admixtools2 (ADMIXTOOLS2-0007): no heterozygous call among the first
    /// 1000 SNPs, which also catches diploid genomes whose first calls happen
    /// to be homozygous.
    First1000,
    /// No heterozygous call at any SNP (the default).
    Genome,
    /// Every sample diploid.
    Diploid,
}

impl PloidyRule {
    /// The SNPs to look at, of `n_snps`.
    pub(crate) fn snps_checked(self, n_snps: usize) -> usize {
        match self {
            PloidyRule::First1000 => 1000.min(n_snps),
            _ => n_snps,
        }
    }

    /// Set the ploidy of the samples of a loaded genotype matrix.
    pub(crate) fn apply(self, geno: &mut GenotypeMatrix) {
        match self {
            PloidyRule::Diploid => geno.set_all_diploid(),
            _ => geno.detect_ploidy(self.snps_checked(geno.n_snps())),
        }
    }

    /// Set the ploidy of the samples of a packed .geno file in `populations`
    /// (the others are left as they are), reading only as far as needed.
    pub(crate) fn apply_packed(
        self,
        geno_path: &str,
        samples: &mut [Sample],
        populations: &[&str],
        n_snps: usize,
    ) -> PyResult<()> {
        let selected: Vec<usize> = (0..samples.len())
            .filter(|&k| populations.contains(&samples[k].population.as_str()))
            .collect();
        if self == PloidyRule::Diploid {
            selected
                .iter()
                .for_each(|&k| samples[k].ploidy = Ploidy::Diploid);
            return Ok(());
        }
        detect_ploidy_packedancestrymap_selected(
            geno_path,
            samples,
            &selected,
            self.snps_checked(n_snps),
            n_snps,
        )
        .map_err(to_py_err)
    }
}

/// Parse the `pseudohaploid` option.
pub(crate) fn ploidy_rule(name: &str) -> PyResult<PloidyRule> {
    match name {
        "admixtools2" => Ok(PloidyRule::First1000),
        "detect" => Ok(PloidyRule::Genome),
        "diploid" => Ok(PloidyRule::Diploid),
        other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "pseudohaploid must be \"detect\", \"admixtools2\" or \"diploid\", got {other:?}"
        ))),
    }
}
