//! Things a caller should be told that are not errors: records left out of
//! an input, results that need care. Functions that can warn take a sink,
//! `warn: &mut dyn FnMut(Warning)`; the Python bindings raise each as a
//! Python warning.

use std::fmt;
use std::path::PathBuf;

use crate::fstats::blocks::NO_MAP_BLOCK_SIZE_BP;
use crate::fstats::{BlockCoordinate, BlockDefinition, CoordinateEvidence};

/// A warning.
///
/// New variants may be added in minor releases.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Warning {
    /// VCF records left out because popmix's formats hold only biallelic
    /// SNPs (two single-base alleles).
    VcfRecordsSkipped {
        path: PathBuf,
        multiallelic: usize,
        not_snp: usize,
    },
    /// Fewer jackknife blocks than [`MIN_BLOCKS`]: standard errors, z and p
    /// from so few are rough (TASKS.md G4).
    FewBlocks { blocks: usize },
    /// Too few jackknife blocks for the number of f4 statistics a qpAdm or
    /// qpWave chi-squared uses: it treats their jackknife covariance as
    /// known, which inflates it (TASKS.md G4).
    FewBlocksForF4 { blocks: usize, statistics: usize },
    /// A leave-one-chromosome-out jackknife over fewer than
    /// [`MIN_JACKKNIFE_CHROMOSOMES`] chromosomes (TASKS.md G4).
    FewChromosomes {
        analysis: String,
        chromosomes: usize,
    },
    /// No genetic map: jackknife blocks fall back to 2 Mb of base-pair
    /// position, as admixtools2 (TASKS.md G5).
    NoGeneticMap,
}

/// Warn when the SNPs seen by `evidence` have no genetic map, so that
/// `block_def`'s blocks fall back to 2 Mb of base-pair position (unless
/// base-pair blocks were asked for, with a block size of 100 or more).
pub fn check_genetic_map(
    evidence: &CoordinateEvidence,
    block_def: &BlockDefinition,
    warn: &mut dyn FnMut(Warning),
) {
    if block_def.block_size < 100.0
        && block_def.coordinate(evidence) == BlockCoordinate::BasePairs(NO_MAP_BLOCK_SIZE_BP)
    {
        warn(Warning::NoGeneticMap);
    }
}

/// Fewest jackknife blocks before [`Warning::FewBlocks`]. With g blocks
/// the SE's own relative error is about 1/sqrt(2(g − 1)): 16% at 20.
pub const MIN_BLOCKS: usize = 20;

/// Fewest chromosomes for a leave-one-chromosome-out jackknife before
/// [`Warning::FewChromosomes`] (the SE's relative error is 24% at 10).
pub const MIN_JACKKNIFE_CHROMOSOMES: usize = 10;

/// Largest tolerated inflation of a chi-squared from g blocks and p
/// statistics before [`Warning::FewBlocksForF4`].
const MAX_INFLATION: f64 = 1.1;

/// How much a chi-squared of `p` statistics is inflated on average when
/// their covariance is estimated from `g` jackknife blocks rather than
/// known: Hotelling's T² from g samples has mean p (g − 1) / (g − p − 2),
/// against p for a known covariance. The jackknife's blocks are not
/// independent samples, so this is a guide. Infinite when g ≤ p + 2.
pub fn chisq_inflation(blocks: usize, statistics: usize) -> f64 {
    if blocks > statistics + 2 {
        (blocks - 1) as f64 / (blocks - statistics - 2) as f64
    } else {
        f64::INFINITY
    }
}

/// Warn when `blocks` is below [`MIN_BLOCKS`].
pub fn check_blocks(blocks: usize, warn: &mut dyn FnMut(Warning)) {
    if blocks < MIN_BLOCKS {
        warn(Warning::FewBlocks { blocks });
    }
}

/// Warn when a chi-squared of `statistics` f4 statistics from `blocks`
/// jackknife blocks is inflated by more than 10%.
pub fn check_blocks_for_f4(blocks: usize, statistics: usize, warn: &mut dyn FnMut(Warning)) {
    if chisq_inflation(blocks, statistics) > MAX_INFLATION {
        warn(Warning::FewBlocksForF4 { blocks, statistics });
    }
}

/// Warn when a jackknife runs over fewer than
/// [`MIN_JACKKNIFE_CHROMOSOMES`] chromosomes.
pub fn check_chromosomes(analysis: &str, chromosomes: usize, warn: &mut dyn FnMut(Warning)) {
    if chromosomes < MIN_JACKKNIFE_CHROMOSOMES {
        warn(Warning::FewChromosomes {
            analysis: analysis.to_string(),
            chromosomes,
        });
    }
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Warning::VcfRecordsSkipped {
                path,
                multiallelic,
                not_snp,
            } => write!(
                f,
                "{}: skipped {} records that are not biallelic SNPs ({multiallelic} multiallelic, \
                 {not_snp} indels or symbolic alleles)",
                path.display(),
                multiallelic + not_snp
            ),
            Warning::FewBlocks { blocks } => write!(
                f,
                "only {blocks} jackknife blocks: standard errors, z and p-values from fewer than \
                 {MIN_BLOCKS} are rough (use more data or a smaller block size)"
            ),
            Warning::FewBlocksForF4 { blocks, statistics } => {
                let inflation = chisq_inflation(*blocks, *statistics);
                if inflation.is_finite() {
                    write!(
                        f,
                        "{blocks} jackknife blocks for {statistics} f4 statistics: chi-squared \
                         treats their jackknife covariance as known, which inflates it about \
                         {inflation:.2}-fold here, so p is too small (use more blocks or fewer \
                         populations)"
                    )
                } else {
                    write!(
                        f,
                        "{blocks} jackknife blocks for {statistics} f4 statistics: too few to \
                         estimate their covariance, so chi-squared and p mean nothing (use more \
                         blocks or fewer populations)"
                    )
                }
            }
            Warning::FewChromosomes {
                analysis,
                chromosomes,
            } => write!(
                f,
                "{analysis}: the jackknife over {chromosomes} chromosomes gives a rough standard \
                 error (fewer than {MIN_JACKKNIFE_CHROMOSOMES})"
            ),
            Warning::NoGeneticMap => write!(
                f,
                "no genetic map (every SNP at the same map position): jackknife blocks are 2 Mb of \
                 base-pair position, as admixtools2, shorter than the 5 cM default, so standard \
                 errors come out slightly small; add genetic positions to the .snp or .bim file"
            ),
        }
    }
}

/// A sink that drops every warning.
pub fn ignore(_: Warning) {}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect(f: impl FnOnce(&mut dyn FnMut(Warning))) -> Vec<Warning> {
        let mut out = Vec::new();
        f(&mut |w| out.push(w));
        out
    }

    #[test]
    fn few_blocks_warn_below_twenty() {
        assert_eq!(collect(|w| check_blocks(20, w)), vec![]);
        let warned = collect(|w| check_blocks(19, w));
        assert_eq!(warned, vec![Warning::FewBlocks { blocks: 19 }]);
        assert_eq!(
            warned[0].to_string(),
            "only 19 jackknife blocks: standard errors, z and p-values from fewer than 20 are \
             rough (use more data or a smaller block size)"
        );
    }

    #[test]
    fn the_inflation_is_hotellings() {
        // (g − 1) / (g − p − 2)
        assert_eq!(chisq_inflation(100, 6), 99.0 / 92.0);
        assert_eq!(chisq_inflation(8, 6), f64::INFINITY);
    }

    #[test]
    fn few_blocks_for_the_f4_statistics_warn_past_ten_percent() {
        // 6 statistics: 1.1 = (g − 1)/(g − 8) at g = 78
        assert_eq!(collect(|w| check_blocks_for_f4(79, 6, w)), vec![]);
        let warned = collect(|w| check_blocks_for_f4(77, 6, w));
        assert_eq!(
            warned[0].to_string(),
            "77 jackknife blocks for 6 f4 statistics: chi-squared treats their jackknife \
             covariance as known, which inflates it about 1.10-fold here, so p is too small (use \
             more blocks or fewer populations)"
        );
        let warned = collect(|w| check_blocks_for_f4(8, 6, w));
        assert_eq!(
            warned[0].to_string(),
            "8 jackknife blocks for 6 f4 statistics: too few to estimate their covariance, so \
             chi-squared and p mean nothing (use more blocks or fewer populations)"
        );
    }

    #[test]
    fn few_chromosomes_warn_below_ten() {
        assert_eq!(collect(|w| check_chromosomes("DATES", 10, w)), vec![]);
        let warned = collect(|w| check_chromosomes("DATES", 3, w));
        assert_eq!(
            warned[0].to_string(),
            "DATES: the jackknife over 3 chromosomes gives a rough standard error (fewer than 10)"
        );
    }

    fn snps(morgans: impl Fn(u32) -> f64) -> Vec<crate::genotypes::Snp> {
        (0..5u32)
            .map(|i| {
                crate::genotypes::Snp::new(
                    format!("rs{i}"),
                    crate::genotypes::Chromosome::Autosome(1),
                    1000 * (i + 1),
                    morgans(i),
                    'A',
                    'G',
                )
            })
            .collect()
    }

    fn map_warnings(s: &[crate::genotypes::Snp], block_size: f64) -> Vec<Warning> {
        let mut evidence = CoordinateEvidence::default();
        evidence.observe(s);
        collect(|w| check_genetic_map(&evidence, &BlockDefinition::new(block_size), w))
    }

    #[test]
    fn snps_without_a_map_warn() {
        let warned = map_warnings(&snps(|_| 0.0), 0.05);
        assert_eq!(warned, vec![Warning::NoGeneticMap]);
        assert!(warned[0].to_string().starts_with("no genetic map"));
    }

    #[test]
    fn a_map_or_base_pair_blocks_asked_for_do_not_warn() {
        assert_eq!(map_warnings(&snps(|i| f64::from(i) * 0.01), 0.05), vec![]);
        assert_eq!(map_warnings(&snps(|_| 0.0), 2e6), vec![]);
    }

    #[test]
    fn a_skipped_vcf_record_warning_says_how_many_and_why() {
        let w = Warning::VcfRecordsSkipped {
            path: "a.vcf".into(),
            multiallelic: 2,
            not_snp: 3,
        };
        assert_eq!(
            w.to_string(),
            "a.vcf: skipped 5 records that are not biallelic SNPs (2 multiallelic, 3 indels or \
             symbolic alleles)"
        );
    }
}
