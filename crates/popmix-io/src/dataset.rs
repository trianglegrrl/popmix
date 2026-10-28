//! A genotype dataset by its prefix, in whichever format it is: found,
//! opened to be read in place where the format allows, or read whole; and
//! how its samples are found to be pseudo-haploid.

use std::path::{Path, PathBuf};
use std::str::FromStr;

use popmix_core::{GenotypeMatrix, Ploidy, PopmixError, Result, Sample, Snp, Warning};

use crate::{
    eigenstrat, is_packed_format, tgeno, with_suffix, EigenstratReader, GenoFormat,
    PackedAncestryMapReader, PedReader, PlinkReader, SnpMajor, VcfReader,
};

/// The VCF file of `prefix`: the path itself, or with .vcf.gz or .vcf added.
// the comparisons are of a lower-cased name
#[allow(clippy::case_sensitive_file_extension_comparisons)]
pub fn vcf_path(prefix: impl AsRef<Path>) -> Option<PathBuf> {
    let prefix = prefix.as_ref();
    // lower-cased, so that a.VCF.GZ is a VCF too
    let name = prefix.as_os_str().to_string_lossy().to_lowercase();

    if (name.ends_with(".vcf") || name.ends_with(".vcf.gz")) && prefix.exists() {
        return Some(prefix.to_path_buf());
    }
    [".vcf.gz", ".vcf"]
        .iter()
        .map(|e| with_suffix(prefix, e))
        .find(|p| p.exists())
}

/// Open a VCF, warning about the records left out because they are not
/// biallelic SNPs.
pub fn open_vcf(path: impl AsRef<Path>, warn: &mut dyn FnMut(Warning)) -> Result<VcfReader> {
    let path = path.as_ref();
    let reader = VcfReader::new(path)?;
    let skipped = reader.skipped();
    if skipped.total() > 0 {
        warn(Warning::VcfRecordsSkipped {
            path: path.into(),
            multiallelic: skipped.multiallelic,
            not_snp: skipped.not_snp,
        });
    }
    Ok(reader)
}

/// The samples of the dataset at `prefix` (or a VCF's full path), from its
/// sample file alone (.fam, .ind, a TGENO header, or the VCF header).
pub fn read_samples(prefix: &str, warn: &mut dyn FnMut(Warning)) -> Result<Vec<Sample>> {
    if let Some(path) = vcf_path(prefix) {
        return Ok(open_vcf(&path, warn)?.samples().to_vec());
    }
    if Path::new(&format!("{prefix}.bed")).exists() {
        return crate::plink::read_fam_file(format!("{prefix}.fam"));
    }
    let geno_path = format!("{prefix}.geno");
    if Path::new(&geno_path).exists() {
        if tgeno::is_tgeno_format(&geno_path)? {
            return Ok(tgeno::TgenoReader::new(prefix)?.samples().to_vec());
        }
        return eigenstrat::read_ind_file(format!("{prefix}.ind"));
    }
    Err(PopmixError::FileNotFound {
        path: PathBuf::from(prefix),
    })
}

/// A dataset whose genotypes can be read in place: a PACKEDANCESTRYMAP
/// .geno, a PLINK .bed or a text EIGENSTRAT .geno of even lines (mapped), or
/// a VCF read once into packed genotypes in memory.
pub struct Mapped {
    pub file: SnpMajor,
    pub samples: Vec<Sample>,
    pub snps: Vec<Snp>,
}

/// `prefix`, in `format`, as a [`Mapped`] dataset; None for the formats
/// that cannot be read in place, and for text whose lines are not all one
/// character per sample.
pub fn open_mapped_as(
    prefix: impl AsRef<Path>,
    format: GenoFormat,
    warn: &mut dyn FnMut(Warning),
) -> Result<Option<Mapped>> {
    let prefix = prefix.as_ref();
    let eigenstrat_parts = || -> Result<(Vec<Sample>, Vec<Snp>)> {
        Ok((
            eigenstrat::read_ind_file(with_suffix(prefix, ".ind"))?,
            eigenstrat::read_snp_file(with_suffix(prefix, ".snp"))?,
        ))
    };
    let (file, samples, snps) = match format {
        GenoFormat::Vcf => match vcf_path(prefix) {
            Some(path) => {
                let (samples, snps, file) = open_vcf(&path, warn)?.into_snp_major();
                (file, samples, snps)
            }
            None => return Ok(None),
        },
        GenoFormat::Eigenstrat => {
            let (samples, snps) = eigenstrat_parts()?;
            match SnpMajor::open_text(with_suffix(prefix, ".geno"), samples.len(), snps.len())? {
                Some(file) => (file, samples, snps),
                None => return Ok(None),
            }
        }
        GenoFormat::PackedAncestrymap => {
            let (samples, snps) = eigenstrat_parts()?;
            let file =
                SnpMajor::open_packed(with_suffix(prefix, ".geno"), samples.len(), snps.len())?;
            (file, samples, snps)
        }
        GenoFormat::PackedPed => {
            let (samples, snps) = PlinkReader::new(prefix)?.into_parts();
            let file = SnpMajor::open_bed(with_suffix(prefix, ".bed"), samples.len(), snps.len())?;
            (file, samples, snps)
        }
        _ => return Ok(None),
    };
    Ok(Some(Mapped {
        file,
        samples,
        snps,
    }))
}

/// `prefix` as a [`Mapped`] dataset, its format found as the readers find
/// it (a packed .geno first, then [`GenoFormat::detect`]); None if it is in
/// another format, or in none (the reader that follows says which).
pub fn open_mapped(
    prefix: impl AsRef<Path>,
    warn: &mut dyn FnMut(Warning),
) -> Result<Option<Mapped>> {
    let prefix = prefix.as_ref();
    let geno_path = with_suffix(prefix, ".geno");
    let format = if Path::new(&geno_path).exists() && is_packed_format(&geno_path)? {
        GenoFormat::PackedAncestrymap
    } else {
        match GenoFormat::detect(prefix) {
            Ok(format) => format,
            Err(_) => return Ok(None),
        }
    };
    open_mapped_as(prefix, format, warn)
}

/// Every genotype of `prefix`, in any format, as stored (no ploidy set).
pub fn read_genotypes(
    prefix: impl AsRef<Path>,
    warn: &mut dyn FnMut(Warning),
) -> Result<GenotypeMatrix> {
    let prefix = prefix.as_ref();
    if let Some(path) = vcf_path(prefix) {
        return open_vcf(&path, warn)?.read_all();
    }
    if with_suffix(prefix, ".bed").exists() {
        return PlinkReader::new(prefix)?.read_all();
    }
    if with_suffix(prefix, ".ped").exists() {
        return PedReader::new(prefix)?.read_all();
    }
    let geno_path = with_suffix(prefix, ".geno");
    if Path::new(&geno_path).exists() {
        if tgeno::is_tgeno_format(&geno_path)? {
            return tgeno::TgenoReader::new(prefix)?.read_all();
        }
        if is_packed_format(&geno_path)? {
            return PackedAncestryMapReader::new(prefix)?.read_all();
        }
        return EigenstratReader::new(prefix)?.read_all();
    }
    Err(PopmixError::FileNotFound {
        path: format!(
            "{} (no genotype files: .vcf, .vcf.gz, .bed, .ped or .geno)",
            prefix.display()
        )
        .into(),
    })
}

/// The genotypes of the samples of `populations`, or more: for TGENO only
/// their records are read, for the other formats the whole file.
pub fn read_genotypes_of(
    prefix: impl AsRef<Path>,
    populations: &[String],
    warn: &mut dyn FnMut(Warning),
) -> Result<GenotypeMatrix> {
    let prefix = prefix.as_ref();
    let geno_path = with_suffix(prefix, ".geno");
    if Path::new(&geno_path).exists() && tgeno::is_tgeno_format(&geno_path)? {
        let reader = tgeno::TgenoReader::new(prefix)?;
        let selected: Vec<usize> = (0..reader.samples().len())
            .filter(|&i| {
                let pop = reader.samples()[i].population.as_str();
                populations.iter().any(|p| p == pop)
            })
            .collect();
        return reader.read_samples(&selected);
    }
    read_genotypes(prefix, warn)
}

/// How samples are found to be pseudo-haploid (counted as one allele per
/// SNP): the `pseudohaploid` option, `"detect"`, `"admixtools2"` or
/// `"diploid"`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum PloidyRule {
    /// `"admixtools2"` (ADMIXTOOLS2-0007): no heterozygous call among the
    /// first 1000 SNPs, which also catches diploid genomes whose first calls
    /// happen to be homozygous.
    First1000,
    /// `"detect"`, the default: no heterozygous call at any SNP.
    #[default]
    Genome,
    /// `"diploid"`: every sample diploid.
    Diploid,
}

impl FromStr for PloidyRule {
    type Err = PopmixError;

    fn from_str(name: &str) -> Result<Self> {
        match name {
            "admixtools2" => Ok(PloidyRule::First1000),
            "detect" => Ok(PloidyRule::Genome),
            "diploid" => Ok(PloidyRule::Diploid),
            other => Err(PopmixError::parameter(
                "pseudohaploid",
                format!("must be \"detect\", \"admixtools2\" or \"diploid\", got {other:?}"),
            )),
        }
    }
}

impl PloidyRule {
    /// The option value that names this rule.
    pub fn name(self) -> &'static str {
        match self {
            PloidyRule::First1000 => "admixtools2",
            PloidyRule::Genome => "detect",
            PloidyRule::Diploid => "diploid",
        }
    }

    /// The SNPs to look at, of `n_snps`.
    pub fn snps_checked(self, n_snps: usize) -> usize {
        match self {
            PloidyRule::First1000 => 1000.min(n_snps),
            _ => n_snps,
        }
    }

    /// Set the ploidy of the samples of a loaded genotype matrix.
    pub fn apply(self, geno: &mut GenotypeMatrix) {
        match self {
            PloidyRule::Diploid => geno.set_all_diploid(),
            _ => geno.detect_ploidy(self.snps_checked(geno.n_snps())),
        }
    }

    /// Set the ploidy of the samples of a mapped genotype file in
    /// `populations` (the others are left as they are), reading only as far
    /// as needed.
    pub fn apply_mapped(self, file: &SnpMajor, samples: &mut [Sample], populations: &[&str]) {
        let selected: Vec<usize> = (0..samples.len())
            .filter(|&k| populations.contains(&samples[k].population.as_str()))
            .collect();
        if self == PloidyRule::Diploid {
            for &k in &selected {
                samples[k].ploidy = Ploidy::Diploid;
            }
            return;
        }
        file.detect_ploidy(samples, &selected, self.snps_checked(file.n_snps()));
    }
}
