//! The DATES analysis: weights, LD table, exponential fit and the
//! leave-one-chromosome-out weighted jackknife, giving the date of
//! admixture in generations and years.
//!
//! Provenance: written from the clean-room specification
//! `docs/cleanroom/dates.md` (section 7), Narasimhan et al. 2019 (Science,
//! supplementary methods of DATES), Chintalapati, Patterson & Moorjani 2022
//! (eLife), DATES's README and parameter-file documentation, and DATES v753
//! run as a black box (`tests/golden/dates/v753/`).

use rayon::prelude::*;

use super::compat::DatesCompat;
use super::exponential_fit::{fit_exponential, weighted_jackknife, DecayFit, FitOptions};
use super::ld_decay::{compute_ld_decay_by_chromosome, LdBin, LdOptions};
use super::weights::{compute_weights, WeightResult};
use crate::error::{PopmixError, Result};
use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
use crate::genotypes::PopulationId;
use crate::options::check_positive;

/// Options of the analysis.
#[derive(Debug, Clone)]
pub struct DatesOptions {
    /// The LD table (its `compat` is replaced by [`DatesOptions::compat`]).
    pub ld_options: LdOptions,
    /// The exponential fit.
    pub fit_options: FitOptions,
    /// Years per generation.
    pub generation_time: f64,
    /// Fewest observed alleles in each source for a SNP to be used.
    pub min_source_count: u32,
    /// Run the leave-one-chromosome-out jackknife.
    pub run_jackknife: bool,
    /// DATES v753 behaviours to reproduce.
    pub compat: DatesCompat,
}

impl Default for DatesOptions {
    fn default() -> Self {
        DatesOptions {
            ld_options: LdOptions::default(),
            fit_options: FitOptions::default(),
            generation_time: 29.0,
            min_source_count: 1,
            run_jackknife: true,
            compat: DatesCompat::CORRECTED,
        }
    }
}

/// One jackknife block: a chromosome left out.
#[derive(Debug, Clone)]
pub struct JackknifeBlock {
    /// The chromosome left out.
    pub chromosome: u8,
    /// Its weight (SNP count).
    pub weight: f64,
    /// The estimate without it (generations).
    pub generations: f64,
}

/// The result of [`dates`].
#[derive(Debug, Clone)]
pub struct DatesResult {
    pub test_pop: PopulationId,
    pub source1: PopulationId,
    pub source2: PopulationId,
    /// Jackknife estimate (or the full fit without the jackknife).
    pub generations: f64,
    /// Jackknife standard error (NaN without the jackknife).
    pub generations_se: f64,
    pub years_bp: f64,
    pub years_bp_se: f64,
    /// The LD table over every chromosome.
    pub ld_bins: Vec<LdBin>,
    /// The full-data fit.
    pub fit: DecayFit,
    /// SNPs with a valid weight.
    pub n_snps: usize,
    pub weights: WeightResult,
    /// Distinct autosomes among the SNPs.
    pub n_chromosomes: usize,
    pub jackknife_by_chromosome: Vec<JackknifeBlock>,
}

/// `x` rounded to `decimals` decimals, as `printf("%.*f")` prints it.
fn printed(x: f64, decimals: usize) -> f64 {
    format!("{x:.decimals$}").parse().unwrap_or(x)
}

/// The fit of `bins`, on their correlations rounded to 6 decimals if
/// `rounded`.
fn fit_bins(bins: &[LdBin], options: &FitOptions, rounded: bool) -> Result<DecayFit> {
    if rounded {
        let r: Vec<LdBin> = bins
            .iter()
            .map(|b| LdBin {
                correlation: printed(b.correlation, 6),
                ..b.clone()
            })
            .collect();
        fit_exponential(&r, Some(options.clone()))
    } else {
        fit_exponential(bins, Some(options.clone()))
    }
}

/// Date the admixture of `test_pop` between `source1` and `source2`.
///
/// `af` supplies the allele frequencies; when `None` they are computed from
/// `geno` for the three populations.
pub fn dates(
    geno: &GenotypeMatrix,
    af: Option<&AlleleFrequencies>,
    test_pop: &str,
    source1: &str,
    source2: &str,
    options: Option<DatesOptions>,
) -> Result<DatesResult> {
    let opts = options.unwrap_or_default();
    check_positive("generation_time", opts.generation_time)?;
    let computed;
    let af = match af {
        Some(af) => af,
        None => {
            computed = AlleleFrequencies::from_genotypes_counting(
                geno,
                Some(&[test_pop, source1, source2]),
                opts.compat.every_sample_diploid,
            );
            &computed
        }
    };

    let weights = compute_weights(af, source1, source2, Some(opts.min_source_count))?;
    if weights.n_valid == 0 {
        return Err(PopmixError::InsufficientData {
            context: "No SNPs with valid weights".to_string(),
            required: 100,
            actual: 0,
        });
    }
    let i1 = af
        .pop_index(source1)
        .ok_or_else(|| PopmixError::PopulationNotFound(source1.to_string()))?;
    let i2 = af
        .pop_index(source2)
        .ok_or_else(|| PopmixError::PopulationNotFound(source2.to_string()))?;
    let p1: Vec<f64> = (0..af.n_snps()).map(|i| af.frequencies[[i, i1]]).collect();
    let p2: Vec<f64> = (0..af.n_snps()).map(|i| af.frequencies[[i, i2]]).collect();

    let ld_options = LdOptions {
        compat: opts.compat,
        ..opts.ld_options.clone()
    };
    let ld = compute_ld_decay_by_chromosome(
        geno,
        af,
        test_pop,
        &weights.weights,
        &p1,
        &p2,
        Some(ld_options),
    )?;
    if ld.bins.is_empty() {
        return Err(PopmixError::InsufficientData {
            context: "No LD bins computed".to_string(),
            required: 10,
            actual: 0,
        });
    }

    let rounded = opts.compat.fit_printed_values;
    let mut fit = fit_bins(&ld.bins, &opts.fit_options, rounded)?;

    let (generations, generations_se, blocks) = if opts.run_jackknife {
        let blocks: Vec<JackknifeBlock> = ld
            .chromosomes
            .par_iter()
            .map(|c| {
                fit_bins(&c.leave_one_out, &opts.fit_options, rounded).map(|f| JackknifeBlock {
                    chromosome: c.chromosome,
                    weight: c.n_snps as f64,
                    generations: f.decay_rate,
                })
            })
            .collect::<Result<_>>()?;
        let at_precision = |x: f64| if rounded { printed(x, 3) } else { x };
        let estimates: Vec<(f64, f64)> = blocks
            .iter()
            .map(|b| (b.weight, at_precision(b.generations)))
            .collect();
        let (est, se) = weighted_jackknife(at_precision(fit.decay_rate), &estimates)?;
        fit.decay_rate_se = se;
        (est, se, blocks)
    } else {
        (fit.decay_rate, f64::NAN, Vec::new())
    };

    let autosomes: std::collections::BTreeSet<_> = geno
        .snps()
        .iter()
        .filter(|s| s.is_autosomal())
        .map(|s| &s.chromosome)
        .collect();

    Ok(DatesResult {
        test_pop: PopulationId::new(test_pop),
        source1: PopulationId::new(source1),
        source2: PopulationId::new(source2),
        generations,
        generations_se,
        years_bp: generations * opts.generation_time,
        years_bp_se: if generations_se.is_finite() {
            generations_se * opts.generation_time
        } else {
            f64::NAN
        },
        ld_bins: ld.bins,
        fit,
        n_snps: weights.n_valid,
        weights,
        n_chromosomes: autosomes.len(),
        jackknife_by_chromosome: blocks,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::{Chromosome, Sample, Snp};
    use ndarray::Array2;

    fn create_test_genotypes() -> GenotypeMatrix {
        // Create test data with 3 populations
        let mut samples = Vec::new();
        for i in 0..20 {
            samples.push(Sample::new(
                format!("test_{}", i),
                PopulationId::new("test"),
            ));
        }
        for i in 0..10 {
            samples.push(Sample::new(
                format!("src1_{}", i),
                PopulationId::new("source1"),
            ));
        }
        for i in 0..10 {
            samples.push(Sample::new(
                format!("src2_{}", i),
                PopulationId::new("source2"),
            ));
        }

        // Create SNPs across multiple chromosomes
        let mut snps = Vec::new();
        for chrom in 1..=3 {
            for i in 0..100 {
                snps.push(Snp::new(
                    format!("rs{}_{}", chrom, i),
                    Chromosome::Autosome(chrom),
                    (i * 10000) as u32,
                    i as f64 * 0.003, // 0.3 cM = 0.003 M spacing
                    'A',
                    'G',
                ));
            }
        }

        let n_samples = samples.len();
        let n_snps = snps.len();
        let mut data = Array2::zeros((n_snps, n_samples));

        // Create genotypes with ancestry-specific patterns
        for snp_idx in 0..n_snps {
            // Source1 pattern (higher frequency)
            for sample_idx in 20..30 {
                data[[snp_idx, sample_idx]] = if snp_idx % 3 == 0 { 2 } else { 1 };
            }
            // Source2 pattern (lower frequency)
            for sample_idx in 30..40 {
                data[[snp_idx, sample_idx]] = if snp_idx % 3 == 0 { 0 } else { 1 };
            }
            // Test population (admixed - mix of patterns with some LD)
            for sample_idx in 0..20 {
                let base = if (snp_idx + sample_idx) % 2 == 0 {
                    data[[snp_idx, 20 + sample_idx % 10]] // Like source1
                } else {
                    data[[snp_idx, 30 + sample_idx % 10]] // Like source2
                };
                data[[snp_idx, sample_idx]] = base;
            }
        }

        GenotypeMatrix::new(data, samples, snps).unwrap()
    }

    #[test]
    fn test_dates_basic() {
        let geno = create_test_genotypes();

        let mut opts = DatesOptions::default();
        opts.ld_options.min_pairs_per_bin = 10; // Lower threshold for test
        opts.run_jackknife = false; // Skip jackknife for faster test

        let result = dates(&geno, None, "test", "source1", "source2", Some(opts));

        // Should complete without error
        assert!(result.is_ok(), "DATES should complete: {:?}", result);

        let result = result.unwrap();

        // Check basic properties
        assert!(result.n_snps > 0);
        assert!(result.generations > 0.0);
        assert!(result.years_bp > 0.0);
        assert!(!result.ld_bins.is_empty());
    }

    /// The test genotypes with source2 of mixed ploidy: its first sample
    /// diploid (heterozygous calls), the others pseudohaploid (0 or 2 only),
    /// as with a diploid genome among pseudohaploid ancient samples.
    fn mixed_ploidy_genotypes() -> GenotypeMatrix {
        let base = create_test_genotypes();
        let mut data = Array2::zeros((base.n_snps(), base.n_samples()));
        for snp_idx in 0..base.n_snps() {
            data.row_mut(snp_idx).assign(&base.snp_raw(snp_idx));
        }
        for snp_idx in 0..data.nrows() {
            for sample_idx in 31..40 {
                if data[[snp_idx, sample_idx]] == 1 {
                    data[[snp_idx, sample_idx]] = if (snp_idx + sample_idx) % 3 == 0 {
                        2
                    } else {
                        0
                    };
                }
            }
        }
        let mut geno =
            GenotypeMatrix::new(data, base.samples().to_vec(), base.snps().to_vec()).unwrap();
        geno.detect_ploidy(geno.n_snps());
        assert_eq!(geno.samples()[30].ploidy, crate::genotypes::Ploidy::Diploid);
        assert_eq!(
            geno.samples()[31].ploidy,
            crate::genotypes::Ploidy::Pseudohaploid
        );
        geno
    }

    fn mixed_ploidy_options(compat: DatesCompat) -> DatesOptions {
        let mut opts = DatesOptions::default();
        opts.ld_options.min_pairs_per_bin = 10;
        opts.run_jackknife = false;
        opts.compat = compat;
        opts
    }

    fn ld_summary(r: &DatesResult) -> Vec<(usize, f64, f64)> {
        r.ld_bins
            .iter()
            .map(|b| (b.n_pairs, b.weighted_ld, b.correlation))
            .collect()
    }

    /// DATES-0004 isolation: with `every_sample_diploid` alone, a source of
    /// mixed ploidy gives exactly what it gives with every sample marked
    /// diploid, and differs from the default, which counts a pseudohaploid
    /// call as one allele.
    #[test]
    fn every_sample_diploid_flag_counts_pseudohaploid_calls_as_two_alleles() {
        let mixed = mixed_ploidy_genotypes();
        let mut diploid = mixed.clone();
        diploid.set_all_diploid();
        let flag = DatesCompat {
            every_sample_diploid: true,
            ..DatesCompat::CORRECTED
        };

        let with_flag = dates(
            &mixed,
            None,
            "test",
            "source1",
            "source2",
            Some(mixed_ploidy_options(flag)),
        )
        .unwrap();
        let as_diploid = dates(
            &diploid,
            None,
            "test",
            "source1",
            "source2",
            Some(mixed_ploidy_options(DatesCompat::CORRECTED)),
        )
        .unwrap();
        let default = dates(
            &mixed,
            None,
            "test",
            "source1",
            "source2",
            Some(mixed_ploidy_options(DatesCompat::CORRECTED)),
        )
        .unwrap();

        assert_eq!(
            with_flag.weights.valid_indices,
            as_diploid.weights.valid_indices
        );
        assert_eq!(
            format!("{:?}", with_flag.weights.weights),
            format!("{:?}", as_diploid.weights.weights)
        );
        assert_eq!(ld_summary(&with_flag), ld_summary(&as_diploid));
        assert_ne!(
            format!("{:?}", with_flag.weights.weights),
            format!("{:?}", default.weights.weights)
        );
        let (v753, corrected) = (DatesCompat::V753, DatesCompat::CORRECTED);
        assert!(v753.every_sample_diploid && !corrected.every_sample_diploid);
    }

    #[test]
    fn test_dates_missing_population() {
        let geno = create_test_genotypes();

        let result = dates(&geno, None, "nonexistent", "source1", "source2", None);
        assert!(result.is_err());
    }
}
