//! How samples are found to be pseudo-haploid (counted as one allele per
//! SNP), shared by every command that reads genotypes.

use popmix_core::python::input_error;
use pyo3::prelude::*;

use popmix_core::{GenotypeMatrix, Ploidy, Sample, Snp};
use popmix_io::{is_packed_format, GenoFormat, SnpMajor};

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

    /// Set the ploidy of the samples of a mapped genotype file in
    /// `populations` (the others are left as they are), reading only as far
    /// as needed.
    pub(crate) fn apply_mapped(
        self,
        file: &SnpMajor,
        samples: &mut [Sample],
        populations: &[&str],
    ) {
        let selected: Vec<usize> = (0..samples.len())
            .filter(|&k| populations.contains(&samples[k].population.as_str()))
            .collect();
        if self == PloidyRule::Diploid {
            selected
                .iter()
                .for_each(|&k| samples[k].ploidy = Ploidy::Diploid);
            return;
        }
        file.detect_ploidy(samples, &selected, self.snps_checked(file.n_snps()));
    }
}

/// A dataset whose genotypes can be read in place: a PACKEDANCESTRYMAP
/// .geno, a PLINK .bed or a text EIGENSTRAT .geno of even lines, with its
/// samples and SNPs.
pub(crate) struct Mapped {
    pub file: SnpMajor,
    pub samples: Vec<Sample>,
    pub snps: Vec<Snp>,
}

/// `prefix` as a [`Mapped`] dataset, or None if it is in another format (as
/// the readers choose: packed .geno first, then as `GenoFormat::detect`).
pub(crate) fn mapped(prefix: &str) -> PyResult<Option<Mapped>> {
    let geno_path = format!("{prefix}.geno");
    let format = if std::path::Path::new(&geno_path).exists() && is_packed_format(&geno_path)? {
        GenoFormat::PackedAncestrymap
    } else {
        // no recognised files: the reader that follows says so
        match GenoFormat::detect(prefix) {
            Ok(format) => format,
            Err(_) => return Ok(None),
        }
    };
    Ok(
        crate::convert::mapped_input(prefix, format)?.map(|(file, samples, snps)| Mapped {
            file,
            samples,
            snps,
        }),
    )
}

/// Parse the `pseudohaploid` option.
pub(crate) fn ploidy_rule(name: &str) -> PyResult<PloidyRule> {
    match name {
        "admixtools2" => Ok(PloidyRule::First1000),
        "detect" => Ok(PloidyRule::Genome),
        "diploid" => Ok(PloidyRule::Diploid),
        other => Err(input_error(format!(
            "pseudohaploid must be \"detect\", \"admixtools2\" or \"diploid\", got {other:?}"
        ))),
    }
}
