//! How samples are found to be pseudo-haploid (counted as one allele per
//! SNP), shared by every command that reads genotypes.

use pyo3::prelude::*;

use popmix_core::{GenotypeMatrix, Ploidy, Sample, Snp};
use popmix_io::{is_packed_format, SnpMajor};

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
/// .geno or a PLINK .bed, with its samples and SNPs.
pub(crate) struct Mapped {
    pub file: SnpMajor,
    pub samples: Vec<Sample>,
    pub snps: Vec<Snp>,
}

/// `prefix` as a [`Mapped`] dataset, or None if it is in another format (as
/// the readers choose: packed .geno, then VCF, then .bed).
pub(crate) fn mapped(prefix: &str) -> PyResult<Option<Mapped>> {
    use std::path::Path;
    let geno_path = format!("{prefix}.geno");
    if Path::new(&geno_path).exists() && is_packed_format(&geno_path).map_err(to_py_err)? {
        let samples =
            popmix_io::eigenstrat::read_ind_file(&format!("{prefix}.ind")).map_err(to_py_err)?;
        let snps =
            popmix_io::eigenstrat::read_snp_file(&format!("{prefix}.snp")).map_err(to_py_err)?;
        let file =
            SnpMajor::open_packed(&geno_path, samples.len(), snps.len()).map_err(to_py_err)?;
        return Ok(Some(Mapped {
            file,
            samples,
            snps,
        }));
    }
    let vcf = [".vcf", ".vcf.gz"]
        .iter()
        .any(|e| Path::new(&format!("{prefix}{e}")).exists())
        || ((prefix.ends_with(".vcf") || prefix.ends_with(".vcf.gz"))
            && Path::new(prefix).exists());
    let bed_path = format!("{prefix}.bed");
    if vcf || !Path::new(&bed_path).exists() {
        return Ok(None);
    }
    let reader = popmix_io::PlinkReader::new(prefix).map_err(to_py_err)?;
    let (samples, snps) = (reader.samples().to_vec(), reader.snps().to_vec());
    let file = SnpMajor::open_bed(&bed_path, samples.len(), snps.len()).map_err(to_py_err)?;
    Ok(Some(Mapped {
        file,
        samples,
        snps,
    }))
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
