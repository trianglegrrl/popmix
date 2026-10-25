//! LD decay computation for DATES algorithm.
//!
//! Computes weighted linkage disequilibrium binned by genetic distance,
//! following DATES v753 (`dates.c`) in both of its modes. With
//! [`DatesCompat::V753`] the output is identical to DATES; by default two
//! mesh FFT quirks are corrected (see [`DatesCompat`]).
//!
//! - `qbin = 0` (DATES default): every SNP pair on a chromosome closer than
//!   `maxdis` is visited directly.
//! - `qbin > 0`: weighted residuals are placed on a genome-wide mesh of
//!   `binsize / qbin` Morgans and pair sums are taken at every mesh lag with
//!   FFT cross-correlation (`setqbins`, `ddadd`, `ddcorr` in DATES).
//!
//! For each individual in the test population:
//! 1. α = regression of (g/2 - p2) on (p1 - p2) over valid SNPs
//! 2. residual = g/2 - (α p1 + (1 - α) p2)
//! 3. y = residual × (p1 - p2)
//!
//! Pairs (y1, y2) are accumulated per chromosome into DATES CORR sums, so the
//! all-chromosome table and each leave-one-chromosome-out table (used by the
//! jackknife) come from the same sums.

use super::compat::DatesCompat;
use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
use crate::genotypes::{Chromosome, Snp};
use crate::threading::ordered_sum;
use crate::Result;
use rustfft::{num_complex::Complex, FftPlanner};

/// DATES `numchrom`: chromosomes outside 1..=22 are not used.
const NUM_CHROMOSOMES: u8 = 22;

/// Options for LD decay computation.
#[derive(Debug, Clone)]
pub struct LdOptions {
    /// Width of each distance bin in cM (DATES `binsize` × 100). Default: 0.1
    pub bin_width_cm: f64,
    /// Maximum distance to consider in cM (DATES `maxdis` × 100). Default: 100.0
    pub max_distance_cm: f64,
    /// Drop bins with fewer pairs than this. Default: 0, which reports every
    /// bin as DATES does.
    pub min_pairs_per_bin: usize,
    /// Mesh cells per distance bin (DATES `qbin`). 0 visits every SNP pair
    /// directly; > 0 uses the FFT mesh. Default: 0 (the DATES default).
    pub qbin: usize,
    /// DATES v753 mesh quirks to reproduce. [`crate::dates::dates`] sets this
    /// from `DatesOptions::compat`. Default: none.
    pub compat: DatesCompat,
}

impl LdOptions {
    /// Most mesh cells (`qbin` per bin) a run may ask for.
    const MAX_MESH_CELLS: f64 = 1e7;

    /// Positive bin width and maximum distance, at most
    /// [`crate::options::MAX_BINS`] bins, and a mesh of a size that fits.
    pub fn validate(&self) -> Result<()> {
        crate::options::check_bins(
            "bin_width_cm",
            self.bin_width_cm,
            "max_distance_cm",
            self.max_distance_cm,
        )?;
        let cells = self.qbin as f64 * self.max_distance_cm / self.bin_width_cm;
        if cells > Self::MAX_MESH_CELLS {
            return Err(crate::error::PopmixError::Validation(format!(
                "qbin {} makes {cells:.0} mesh cells; at most {:.0} are allowed",
                self.qbin,
                Self::MAX_MESH_CELLS
            )));
        }
        Ok(())
    }
}

impl Default for LdOptions {
    fn default() -> Self {
        LdOptions {
            bin_width_cm: 0.1,      // DATES binsize = 0.001 Morgan
            max_distance_cm: 100.0, // DATES default: maxdis = 1.0 Morgan = 100 cM
            min_pairs_per_bin: 0,
            qbin: 0,
            compat: DatesCompat::CORRECTED,
        }
    }
}

/// A single distance bin for LD decay, as written to the DATES `.out` file.
#[derive(Debug, Clone)]
pub struct LdBin {
    /// Bin distance in cM: (bin index + 1) × bin width.
    pub distance_cm: f64,
    /// Weighted LD covariance without mean subtraction (v12 = S12/S0).
    pub weighted_ld: f64,
    /// Z-score: sqrt(S0) × correlation.
    pub z_score: f64,
    /// Correlation without mean subtraction (v12 / sqrt(v11 × v22)). This is
    /// the column DATES fits.
    pub correlation: f64,
    /// Regression coefficient S12/S11.
    pub regression_coef: f64,
    /// Number of SNP pairs in this bin.
    pub n_pairs: usize,
    /// Approximate standard error of the correlation (Fisher, 1/sqrt(n - 3)).
    pub std_err: f64,
}

/// LD decay for all chromosomes together and with each chromosome left out.
#[derive(Debug, Clone)]
pub struct LdDecay {
    /// Bins computed from all chromosomes.
    pub bins: Vec<LdBin>,
    /// One entry per chromosome 1..=22 present in the data, in order.
    pub chromosomes: Vec<ChromosomeLd>,
}

/// Leave-one-chromosome-out LD decay for the jackknife.
#[derive(Debug, Clone)]
pub struct ChromosomeLd {
    /// Chromosome number (1..=22).
    pub chromosome: u8,
    /// Number of SNPs on this chromosome in the input, before any filtering.
    /// DATES uses this as the jackknife weight.
    pub n_snps: usize,
    /// Bins computed from all other chromosomes.
    pub leave_one_out: Vec<LdBin>,
}

/// Accumulator for correlation statistics (DATES CORR struct).
#[derive(Debug, Clone, Default)]
struct CorrBin {
    s0: f64,
    s1: f64,
    s2: f64,
    s11: f64,
    s12: f64,
    s22: f64,
}

impl CorrBin {
    /// DATES `addcorr`.
    fn add(&mut self, y1: f64, y2: f64) {
        self.s0 += 1.0;
        self.s1 += y1;
        self.s2 += y2;
        self.s11 += y1 * y1;
        self.s12 += y1 * y2;
        self.s22 += y2 * y2;
    }

    fn merge(&mut self, other: &CorrBin) {
        self.s0 += other.s0;
        self.s1 += other.s1;
        self.s2 += other.s2;
        self.s11 += other.s11;
        self.s12 += other.s12;
        self.s22 += other.s22;
    }

    fn minus(&self, other: &CorrBin) -> CorrBin {
        CorrBin {
            s0: self.s0 - other.s0,
            s1: self.s1 - other.s1,
            s2: self.s2 - other.s2,
            s11: self.s11 - other.s11,
            s12: self.s12 - other.s12,
            s22: self.s22 - other.s22,
        }
    }

    /// One `.out` row: DATES `dumpit` with `calccorr(corrpt, 1, NO)`, which
    /// does not subtract means. Bins with no pairs are reported as zeros.
    fn to_ld_bin(&self, distance_cm: f64) -> LdBin {
        let (weighted_ld, correlation, z_score) = if self.s0 < 0.5 {
            (0.0, 0.0, 0.0)
        } else {
            let n = self.s0;
            let v11 = self.s11 / n;
            let v12 = self.s12 / n;
            let v22 = self.s22 / n;
            let corr = v12 / (v11 * v22 + 1.0e-20).sqrt();
            (v12, corr, n.sqrt() * corr)
        };
        let n = self.s0;
        LdBin {
            distance_cm,
            weighted_ld,
            z_score,
            correlation,
            regression_coef: self.s12 / (self.s11 + 1.0e-20),
            n_pairs: n.round() as usize,
            std_err: if n > 3.0 {
                1.0 / (n - 3.0).sqrt()
            } else {
                f64::NAN
            },
        }
    }
}

/// Compute weighted LD decay by genetic distance using the DATES algorithm.
///
/// # Arguments
/// * `geno` - Genotype matrix (needed for individual genotypes)
/// * `af` - Allele frequencies (unused, kept for API compatibility)
/// * `test_pop` - Test (admixed) population
/// * `weights` - Per-SNP weights (p_source1 - p_source2), NaN if invalid
/// * `source1_freqs` - Allele frequencies in source1 population
/// * `source2_freqs` - Allele frequencies in source2 population
/// * `options` - Computation options
///
/// # Returns
/// The LD decay bins for all chromosomes together.
pub fn compute_ld_decay(
    geno: &GenotypeMatrix,
    af: &AlleleFrequencies,
    test_pop: &str,
    weights: &[f64],
    source1_freqs: &[f64],
    source2_freqs: &[f64],
    options: Option<LdOptions>,
) -> Result<Vec<LdBin>> {
    compute_ld_decay_by_chromosome(
        geno,
        af,
        test_pop,
        weights,
        source1_freqs,
        source2_freqs,
        options,
    )
    .map(|ld| ld.bins)
}

/// Compute weighted LD decay for all chromosomes and leaving out each one.
///
/// Arguments are as for [`compute_ld_decay`].
pub fn compute_ld_decay_by_chromosome(
    geno: &GenotypeMatrix,
    _af: &AlleleFrequencies,
    test_pop: &str,
    weights: &[f64],
    source1_freqs: &[f64],
    source2_freqs: &[f64],
    options: Option<LdOptions>,
) -> Result<LdDecay> {
    let opts = options.unwrap_or_default();
    opts.validate()?;

    let test_pop_id = crate::genotypes::PopulationId::new(test_pop);
    let test_indices = geno.population_indices(&test_pop_id);
    if test_indices.is_empty() {
        return Err(crate::error::PopmixError::PopulationNotFound(
            test_pop.to_string(),
        ));
    }

    let snps = geno.snps();
    let binsize = opts.bin_width_cm / 100.0;
    let maxdis = opts.max_distance_cm / 100.0;
    // DATES: numbins = nnint(maxdis/binsize) + 5, of which numbins - 5 are reported
    let n_bins = (maxdis / binsize).round() as usize;
    let n_alloc = n_bins + 5;

    // Chromosomes present in the input and their unfiltered SNP counts
    let mut n_snps_by_chrom = [0usize; NUM_CHROMOSOMES as usize + 1];
    for snp in snps {
        if let Some(c) = chromosome_number(&snp.chromosome) {
            n_snps_by_chrom[c as usize] += 1;
        }
    }
    let chroms: Vec<u8> = (1..=NUM_CHROMOSOMES)
        .filter(|&c| n_snps_by_chrom[c as usize] > 0)
        .collect();
    let mut chrom_slot = [usize::MAX; NUM_CHROMOSOMES as usize + 1];
    for (slot, &c) in chroms.iter().enumerate() {
        chrom_slot[c as usize] = slot;
    }

    // SNPs DATES does not ignore: |p1 - p2| >= 0.001 and p(1-p) >= 0.001 at the
    // mean source frequency. Kept in file order within each chromosome.
    let mut snps_by_chrom: Vec<Vec<usize>> = vec![Vec::new(); chroms.len()];
    for (idx, snp) in snps.iter().enumerate() {
        let Some(c) = chromosome_number(&snp.chromosome) else {
            continue;
        };
        let w = weights[idx];
        if !w.is_finite() || w.abs() < 0.001 {
            continue;
        }
        let mean_freq = (source1_freqs[idx] + source2_freqs[idx]) * 0.5;
        if mean_freq * (1.0 - mean_freq) < 0.001 {
            continue;
        }
        snps_by_chrom[chrom_slot[c as usize]].push(idx);
    }
    if let Some(snp) = snps_by_chrom
        .iter()
        .flatten()
        .map(|&i| &snps[i])
        .find(|s| !s.position_morgans.is_finite())
    {
        return Err(crate::error::PopmixError::Validation(format!(
            "SNP {} has no genetic position ({}); DATES needs one for every SNP",
            snp.id, snp.position_morgans
        )));
    }
    for chrom_snps in &mut snps_by_chrom {
        chrom_snps.sort_by(|&a, &b| {
            snps[a]
                .position_morgans
                .total_cmp(&snps[b].position_morgans)
        });
    }

    let residuals = |sample_idx: usize| {
        weighted_residuals(
            geno,
            &snps_by_chrom,
            sample_idx,
            weights,
            source1_freqs,
            source2_freqs,
        )
    };

    let chrom_bins: Vec<Vec<CorrBin>> = if opts.qbin == 0 {
        ordered_sum(
            test_indices.len(),
            || vec![vec![CorrBin::default(); n_alloc]; chroms.len()],
            |i| {
                exact_individual_bins(
                    snps,
                    &snps_by_chrom,
                    &residuals(test_indices[i]),
                    binsize,
                    maxdis,
                    n_alloc,
                )
            },
            |a, b| merge_chrom_bins(a, b),
        )
    } else {
        let qbin = opts.qbin;
        let cells = mesh_cells(snps, binsize, qbin);
        // DATES: diffmax = nnint(qbin * maxdis / binsize)
        let diffmax = (qbin as f64 * maxdis / binsize).round() as usize;
        // per-sample sums are added in sample order, whatever the thread count
        let lag_sums = ordered_sum(
            test_indices.len(),
            || vec![MeshLagSums::new(diffmax); chroms.len()],
            |i| {
                mesh_individual_lag_sums(
                    &snps_by_chrom,
                    &residuals(test_indices[i]),
                    &cells,
                    diffmax,
                    opts.compat,
                )
            },
            |a, b| {
                for (x, y) in a.iter_mut().zip(b.iter()) {
                    x.merge(y);
                }
            },
        );
        lag_sums
            .iter()
            .map(|sums| sums.to_bins(binsize, qbin, n_alloc))
            .collect()
    };

    let mut total = vec![CorrBin::default(); n_alloc];
    for bins in &chrom_bins {
        for (t, b) in total.iter_mut().zip(bins.iter()) {
            t.merge(b);
        }
    }

    let report = |bins: &[CorrBin]| -> Vec<LdBin> {
        bins.iter()
            .take(n_bins)
            .enumerate()
            .filter(|(_, b)| b.s0 >= opts.min_pairs_per_bin as f64)
            .map(|(k, b)| b.to_ld_bin((k as f64 + 1.0) * opts.bin_width_cm))
            .collect()
    };

    let chromosomes = chroms
        .iter()
        .zip(chrom_bins.iter())
        .map(|(&chromosome, bins)| {
            let loo: Vec<CorrBin> = total
                .iter()
                .zip(bins.iter())
                .map(|(t, b)| t.minus(b))
                .collect();
            ChromosomeLd {
                chromosome,
                n_snps: n_snps_by_chrom[chromosome as usize],
                leave_one_out: report(&loo),
            }
        })
        .collect();

    Ok(LdDecay {
        bins: report(&total),
        chromosomes,
    })
}

fn chromosome_number(chrom: &Chromosome) -> Option<u8> {
    match chrom {
        Chromosome::Autosome(n) if (1..=NUM_CHROMOSOMES).contains(n) => Some(*n),
        _ => None,
    }
}

fn merge_chrom_bins(a: &mut [Vec<CorrBin>], b: Vec<Vec<CorrBin>>) {
    for (xa, xb) in a.iter_mut().zip(b.iter()) {
        for (x, y) in xa.iter_mut().zip(xb.iter()) {
            x.merge(y);
        }
    }
}

/// Weighted residuals y = (g/2 - prediction) × (p1 - p2) for one individual,
/// as `(snp_idx, y)` per chromosome slot, skipping missing genotypes.
fn weighted_residuals(
    geno: &GenotypeMatrix,
    snps_by_chrom: &[Vec<usize>],
    sample_idx: usize,
    weights: &[f64],
    source1_freqs: &[f64],
    source2_freqs: &[f64],
) -> Vec<Vec<(usize, f64)>> {
    // α = dot(g - p2, p1 - p2) / dot(p1 - p2, p1 - p2)
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let observed: Vec<Vec<(usize, f64)>> = snps_by_chrom
        .iter()
        .map(|chrom_snps| {
            chrom_snps
                .iter()
                .filter_map(|&snp_idx| {
                    let g = geno.snp_raw(snp_idx)[sample_idx];
                    if g > 2 {
                        return None;
                    }
                    let geno_val = g as f64 / 2.0;
                    let w = weights[snp_idx];
                    numerator += (geno_val - source2_freqs[snp_idx]) * w;
                    denominator += w * w;
                    Some((snp_idx, geno_val))
                })
                .collect()
        })
        .collect();

    let alpha = if denominator > 1e-15 {
        numerator / denominator
    } else {
        0.5
    };

    observed
        .into_iter()
        .map(|chrom_obs| {
            chrom_obs
                .into_iter()
                .map(|(snp_idx, geno_val)| {
                    let prediction =
                        alpha * source1_freqs[snp_idx] + (1.0 - alpha) * source2_freqs[snp_idx];
                    (snp_idx, (geno_val - prediction) * weights[snp_idx])
                })
                .collect()
        })
        .collect()
}

/// DATES `qbin = 0` main loop: every pair closer than `maxdis` on a chromosome.
fn exact_individual_bins(
    snps: &[Snp],
    snps_by_chrom: &[Vec<usize>],
    residuals: &[Vec<(usize, f64)>],
    binsize: f64,
    maxdis: f64,
    n_alloc: usize,
) -> Vec<Vec<CorrBin>> {
    let mut chrom_bins = vec![vec![CorrBin::default(); n_alloc]; snps_by_chrom.len()];

    for (bins, chrom_res) in chrom_bins.iter_mut().zip(residuals.iter()) {
        for (i, &(snp_i, y1)) in chrom_res.iter().enumerate() {
            let pos_i = snps[snp_i].position_morgans;
            for &(snp_j, y2) in &chrom_res[i + 1..] {
                let dis = snps[snp_j].position_morgans - pos_i;
                if dis >= maxdis {
                    break;
                }
                // DATES: s = (int) (dis/binsize)
                let s = (dis / binsize) as usize;
                if s < n_alloc {
                    bins[s].add(y1, y2);
                }
            }
        }
    }

    chrom_bins
}

/// Gap (Morgans) put between chromosomes on the genome-wide mesh axis.
const CHROMOSOME_GAP: f64 = 5.0;

/// The mesh cell of every SNP, in the order given.
///
/// The SNPs are laid end to end on one axis starting at 0: each step to the
/// next SNP adds the difference of their genetic positions on the same
/// chromosome and [`CHROMOSOME_GAP`] when the chromosome changes. A SNP's
/// cell is floor(axis / (binsize / qbin)).
///
/// Provenance: written from the clean-room specification
/// `docs/cleanroom/dates-mesh.md`, after Narasimhan et al. 2019 (Science,
/// supplementary methods of DATES), and checked against DATES v753 run as a
/// black box.
fn mesh_cells(snps: &[Snp], binsize: f64, qbin: usize) -> Vec<usize> {
    let width = binsize / qbin as f64;
    let mut axis = 0.0;
    snps.iter()
        .enumerate()
        .map(|(i, snp)| {
            if i > 0 {
                let prev = &snps[i - 1];
                axis += if prev.chromosome == snp.chromosome {
                    snp.position_morgans - prev.position_morgans
                } else {
                    CHROMOSOME_GAP
                };
            }
            (axis / width).floor() as usize
        })
        .collect()
}

/// Per-lag sums over mesh cell pairs (i, i + d) on one chromosome, for lags
/// d = 0..=D. With n(i), Σy(i) and Σy²(i) the count, sum and sum of squares
/// of the values in cell i, each field holds at lag d the sum over i of the
/// product named.
///
/// Provenance: as [`mesh_cells`] (`docs/cleanroom/dates-mesh.md`).
#[derive(Debug, Clone)]
struct MeshLagSums {
    /// n(i) · n(i + d)
    pairs: Vec<f64>,
    /// n(i) · Σy(i + d)
    later: Vec<f64>,
    /// Σy(i) · n(i + d)
    earlier: Vec<f64>,
    /// Σy(i) · Σy(i + d)
    cross: Vec<f64>,
    /// n(i) · Σy²(i + d)
    later_sq: Vec<f64>,
    /// Σy²(i) · n(i + d)
    earlier_sq: Vec<f64>,
}

impl MeshLagSums {
    /// All zero, lags 0..=`max_lag`.
    fn new(max_lag: usize) -> Self {
        let z = vec![0.0; max_lag + 1];
        MeshLagSums {
            pairs: z.clone(),
            later: z.clone(),
            earlier: z.clone(),
            cross: z.clone(),
            later_sq: z.clone(),
            earlier_sq: z,
        }
    }

    fn merge(&mut self, other: &MeshLagSums) {
        let add = |a: &mut Vec<f64>, b: &[f64]| a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
        add(&mut self.pairs, &other.pairs);
        add(&mut self.later, &other.later);
        add(&mut self.earlier, &other.earlier);
        add(&mut self.cross, &other.cross);
        add(&mut self.later_sq, &other.later_sq);
        add(&mut self.earlier_sq, &other.earlier_sq);
    }

    /// Distance bins: lag d >= 1 with at least one pair (pair count >= 0.5,
    /// the FFT leaves rounding error) is at distance d · binsize / qbin and
    /// goes to bin trunc(distance / binsize) if that is below `n_bins`. The
    /// later SNP of a pair is y1 in [`CorrBin`].
    fn to_bins(&self, binsize: f64, qbin: usize, n_bins: usize) -> Vec<CorrBin> {
        let width = binsize / qbin as f64;
        let mut bins = vec![CorrBin::default(); n_bins];
        for d in 1..self.pairs.len() {
            if self.pairs[d] < 0.5 {
                continue;
            }
            let bin = (d as f64 * width / binsize) as usize;
            let Some(b) = bins.get_mut(bin) else {
                continue;
            };
            b.s0 += self.pairs[d];
            b.s1 += self.later[d];
            b.s2 += self.earlier[d];
            b.s12 += self.cross[d];
            b.s11 += self.later_sq[d];
            b.s22 += self.earlier_sq[d];
        }
        bins
    }
}

/// DATES mesh accumulation (`z0q`, `z1q`, `z2q`) and `ddadd` for one individual.
fn mesh_individual_lag_sums(
    snps_by_chrom: &[Vec<usize>],
    residuals: &[Vec<(usize, f64)>],
    cells: &[usize],
    diffmax: usize,
    compat: DatesCompat,
) -> Vec<MeshLagSums> {
    let mut planner = FftPlanner::<f64>::new();

    snps_by_chrom
        .iter()
        .zip(residuals.iter())
        .map(|(_, chrom_res)| {
            let mut sums = MeshLagSums::new(diffmax);
            let (Some(slo), Some(shi)) = (
                chrom_res.iter().map(|&(i, _)| cells[i]).min(),
                chrom_res.iter().map(|&(i, _)| cells[i]).max(),
            ) else {
                return sums;
            };
            let len = shi - slo + 1;
            if len <= 1 {
                return sums;
            }

            let mut z0 = vec![0.0; len];
            let mut z1 = vec![0.0; len];
            let mut z2 = vec![0.0; len];
            for &(snp_idx, y) in chrom_res {
                let s = cells[snp_idx] - slo;
                z0[s] += 1.0;
                z1[s] += y;
                z2[s] += y * y;
            }

            // DATES fftsubs.c FFT length: 2^(ceil(log2 m) + 1). This is only
            // long enough for lags up to n - m; see the output loop below.
            let n = 2f64.powf((len as f64).log2().ceil() + 1.0).round() as usize;
            let fft = planner.plan_fft_forward(n);
            let ifft = planner.plan_fft_inverse(n);
            let forward = |z: &[f64]| {
                let mut buf: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n];
                for (b, &v) in buf.iter_mut().zip(z) {
                    b.re = v;
                }
                fft.process(&mut buf);
                buf
            };
            let a0 = forward(&z0);
            let a1 = forward(&z1);
            let a2 = forward(&z2);

            // Circular correlation c[d] = Σ_i a[i] b[(i + d) mod n] = IFFT(conj(A) B)[d];
            // c[(n - d) mod n] gives b ahead of a.
            let cross = |a: &[Complex<f64>], b: &[Complex<f64>]| {
                let mut buf: Vec<Complex<f64>> =
                    a.iter().zip(b).map(|(x, y)| x.conj() * y).collect();
                ifft.process(&mut buf);
                let scale = 1.0 / n as f64;
                buf.into_iter().map(|c| c.re * scale).collect::<Vec<f64>>()
            };
            let c00 = cross(&a0, &a0);
            let c01 = cross(&a0, &a1);
            let c11 = cross(&a1, &a1);
            let c02 = cross(&a0, &a2);

            // DATES fftauto / fftcv2 output has two quirks that only matter
            // when diffmax > n - m (maxdis long relative to the chromosome's
            // mesh window); `compat` chooses whether to reproduce them:
            // - DATES-0001: the correlation is circular, so lags beyond n - m
            //   also include wrapped negative-lag (short-range) pairs. Real
            //   pairs have lags below m, and n >= 2m, so the corrected output
            //   just stops at lag m - 1.
            // - DATES-0002: fftauto scales only cout[0..outmax) by 1/n, so its
            //   last lag outmax = min(n - 1, diffmax) is left n times too large.
            let outmax = (n - 1).min(diffmax);
            let last_scale = if compat.unscaled_last_fft_lag {
                n as f64
            } else {
                1.0
            };
            let auto = |c: &[f64], d: usize| match d.cmp(&outmax) {
                std::cmp::Ordering::Less => c[d],
                std::cmp::Ordering::Equal => c[d] * last_scale,
                std::cmp::Ordering::Greater => 0.0,
            };
            let maxlag = if compat.circular_fft_lags {
                diffmax
            } else {
                diffmax.min(len - 1)
            };
            for d in 0..=maxlag {
                sums.pairs[d] = auto(&c00, d);
                sums.cross[d] = auto(&c11, d);
                if d < n {
                    sums.later[d] = c01[d];
                    sums.later_sq[d] = c02[d];
                }
                if d <= n {
                    let neg = (n - d) % n;
                    sums.earlier[d] = c01[neg];
                    sums.earlier_sq[d] = c02[neg];
                }
            }
            sums
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::{PopulationId, Sample, Snp};
    use ndarray::Array2;

    fn create_test_data() -> (GenotypeMatrix, AlleleFrequencies, Vec<f64>) {
        create_test_data_with_spacing(0.005)
    }

    fn create_test_data_with_spacing(
        spacing_morgans: f64,
    ) -> (GenotypeMatrix, AlleleFrequencies, Vec<f64>) {
        let mut samples = Vec::new();
        for i in 0..10 {
            samples.push(Sample::new(
                format!("test_{}", i),
                PopulationId::new("test"),
            ));
        }
        for i in 0..5 {
            samples.push(Sample::new(
                format!("src1_{}", i),
                PopulationId::new("source1"),
            ));
        }
        for i in 0..5 {
            samples.push(Sample::new(
                format!("src2_{}", i),
                PopulationId::new("source2"),
            ));
        }

        let snps: Vec<Snp> = (0..50)
            .map(|i| {
                Snp::new(
                    format!("rs{}", i),
                    Chromosome::Autosome(1),
                    (i * 10000) as u32,
                    i as f64 * spacing_morgans,
                    'A',
                    'G',
                )
            })
            .collect();

        let n_samples = samples.len();
        let n_snps = snps.len();
        let mut data = Array2::zeros((n_snps, n_samples));

        for snp_idx in 0..n_snps {
            for sample_idx in 0..n_samples {
                let base = ((snp_idx + sample_idx) % 3) as u8;
                data[[snp_idx, sample_idx]] = base;
            }
        }

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        let af = AlleleFrequencies::from_genotypes(&geno);

        let idx1 = af.pop_index("source1").unwrap();
        let idx2 = af.pop_index("source2").unwrap();
        let weights: Vec<f64> = (0..n_snps)
            .map(|snp_idx| af.frequencies[[snp_idx, idx1]] - af.frequencies[[snp_idx, idx2]])
            .collect();

        (geno, af, weights)
    }

    fn source_freqs(geno: &GenotypeMatrix, af: &AlleleFrequencies) -> (Vec<f64>, Vec<f64>) {
        let idx1 = af.pop_index("source1").unwrap();
        let idx2 = af.pop_index("source2").unwrap();
        let n_snps = geno.n_snps();
        (
            (0..n_snps).map(|i| af.frequencies[[i, idx1]]).collect(),
            (0..n_snps).map(|i| af.frequencies[[i, idx2]]).collect(),
        )
    }

    #[test]
    fn bins_do_not_depend_on_the_thread_count() {
        // 64 test samples with irregular genotypes, so the order in which
        // their sums are added shows in the last bits.
        let mut samples: Vec<Sample> = (0..64)
            .map(|i| Sample::new(format!("test_{i}"), PopulationId::new("test")))
            .collect();
        for (pop, n) in [("source1", 8), ("source2", 8)] {
            samples
                .extend((0..n).map(|i| Sample::new(format!("{pop}_{i}"), PopulationId::new(pop))));
        }
        let snps: Vec<Snp> = (0..200)
            .map(|i| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1 + (i / 100) as u8),
                    (i * 7919) as u32,
                    (i % 100) as f64 * 0.0031,
                    'A',
                    'G',
                )
            })
            .collect();
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let data = Array2::from_shape_fn((snps.len(), samples.len()), |_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % 3) as u8
        });
        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
        let af = AlleleFrequencies::from_genotypes(&geno);
        let (s1, s2) = source_freqs(&geno, &af);
        let weights: Vec<f64> = s1.iter().zip(&s2).map(|(a, b)| a - b).collect();
        for qbin in [0, 10] {
            let run = |threads| {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let opts = LdOptions {
                    bin_width_cm: 0.05,
                    max_distance_cm: 2.0,
                    min_pairs_per_bin: 1,
                    qbin,
                    ..Default::default()
                };
                let bins = pool
                    .install(|| {
                        compute_ld_decay(&geno, &af, "test", &weights, &s1, &s2, Some(opts))
                    })
                    .unwrap();
                bins.iter()
                    .map(|b| (b.weighted_ld.to_bits(), b.correlation.to_bits()))
                    .collect::<Vec<_>>()
            };
            let one = run(1);
            assert!(!one.is_empty());
            for threads in [2, 3, 8] {
                assert_eq!(run(threads), one, "qbin {qbin}, {threads} threads");
            }
        }
    }

    #[test]
    fn test_ld_decay_basic() {
        let (geno, af, weights) = create_test_data();
        let (source1_freqs, source2_freqs) = source_freqs(&geno, &af);

        let opts = LdOptions {
            bin_width_cm: 1.0,
            max_distance_cm: 25.0,
            min_pairs_per_bin: 1,
            qbin: 0,
            ..Default::default()
        };

        let result = compute_ld_decay(
            &geno,
            &af,
            "test",
            &weights,
            &source1_freqs,
            &source2_freqs,
            Some(opts),
        )
        .unwrap();

        assert!(!result.is_empty(), "Should have LD decay bins");

        for bin in &result {
            assert!(bin.distance_cm > 0.0);
            assert!(bin.n_pairs > 0);
        }
    }

    #[test]
    fn test_mesh_matches_exact_when_snps_sit_on_mesh_cells() {
        // SNPs exactly one bin (1/128 Morgan) apart with qbin = 2 fall in
        // distinct mesh cells at lags that are exact multiples of the cell
        // size, so both DATES modes see the same pairs in the same bins.
        let (geno, af, weights) = create_test_data_with_spacing(1.0 / 128.0);
        let (source1_freqs, source2_freqs) = source_freqs(&geno, &af);
        let run = |qbin| {
            let opts = LdOptions {
                bin_width_cm: 100.0 / 128.0,
                max_distance_cm: 10.0 * 100.0 / 128.0,
                min_pairs_per_bin: 0,
                qbin,
                ..Default::default()
            };
            compute_ld_decay(
                &geno,
                &af,
                "test",
                &weights,
                &source1_freqs,
                &source2_freqs,
                Some(opts),
            )
            .unwrap()
        };
        let exact = run(0);
        let mesh = run(2);
        assert_eq!(exact.len(), mesh.len());
        for (e, m) in exact.iter().zip(mesh.iter()) {
            assert_eq!(e.n_pairs, m.n_pairs, "pairs at {} cM", e.distance_cm);
            assert!((e.weighted_ld - m.weighted_ld).abs() < 1e-12);
            assert!((e.correlation - m.correlation).abs() < 1e-9);
        }
    }

    // Mesh geometry for create_test_data_with_spacing(1/128) with qbin = 1:
    // 50 SNPs in consecutive mesh cells span a window of m = 50 cells, so
    // DATES's FFT length is n = 2^(ceil(log2 50) + 1) = 128. Real pairs have
    // lags of at most 49 cells. A circular lag d > n - m = 78 also picks up
    // the pairs at negative lag n - d, from bin 79 on. fftauto leaves lag
    // n - 1 = 127 unscaled, and it holds the wrapped lag-1 pairs.
    const MESH_TEST_BINS: usize = 130;
    const FIRST_WRAPPED_BIN: usize = 79;
    const UNSCALED_BIN: usize = 127;

    fn run_mesh_test(qbin: usize, compat: DatesCompat) -> Vec<LdBin> {
        let (geno, af, weights) = create_test_data_with_spacing(1.0 / 128.0);
        let (source1_freqs, source2_freqs) = source_freqs(&geno, &af);
        let opts = LdOptions {
            bin_width_cm: 100.0 / 128.0,
            max_distance_cm: MESH_TEST_BINS as f64 * 100.0 / 128.0,
            min_pairs_per_bin: 0,
            qbin,
            compat,
        };
        compute_ld_decay(
            &geno,
            &af,
            "test",
            &weights,
            &source1_freqs,
            &source2_freqs,
            Some(opts),
        )
        .unwrap()
    }

    fn same_bin(a: &LdBin, b: &LdBin) -> bool {
        a.n_pairs == b.n_pairs
            && (a.weighted_ld - b.weighted_ld).abs() < 1e-12
            && (a.correlation - b.correlation).abs() < 1e-9
    }

    /// Indices of the bins where two LD tables differ.
    fn differing_bins(a: &[LdBin], b: &[LdBin]) -> Vec<usize> {
        assert_eq!(a.len(), b.len());
        (0..a.len()).filter(|&k| !same_bin(&a[k], &b[k])).collect()
    }

    #[test]
    fn test_corrected_mesh_matches_exact_past_the_fft_window() {
        // The oracle: on this grid the mesh sees exactly the pairs the exact
        // mode sees, at every distance, so the corrected mesh must agree with
        // it everywhere, including lags DATES wraps.
        let exact = run_mesh_test(0, DatesCompat::CORRECTED);
        let mesh = run_mesh_test(1, DatesCompat::CORRECTED);
        assert_eq!(differing_bins(&exact, &mesh), Vec::<usize>::new());
    }

    #[test]
    fn test_v753_mesh_reproduces_dates_quirks() {
        let exact = run_mesh_test(0, DatesCompat::CORRECTED);
        let mesh = run_mesh_test(1, DatesCompat::V753);
        let diff = differing_bins(&exact, &mesh);
        assert_eq!(diff.first(), Some(&FIRST_WRAPPED_BIN));
        assert!(diff.contains(&UNSCALED_BIN));
    }

    #[test]
    fn test_circular_fft_lags_flag_only_changes_wrapped_bins() {
        let only = DatesCompat {
            circular_fft_lags: true,
            ..DatesCompat::CORRECTED
        };
        let fixed = run_mesh_test(1, DatesCompat::CORRECTED);
        let quirk = run_mesh_test(1, only);
        let diff = differing_bins(&fixed, &quirk);
        assert!(!diff.is_empty());
        assert!(diff.iter().all(|&k| k >= FIRST_WRAPPED_BIN), "{diff:?}");
    }

    #[test]
    fn test_unscaled_last_lag_flag_only_changes_its_bin() {
        // Lag n - 1 has no real pairs, so the missing scaling only shows when
        // the circular quirk puts wrapped pairs there.
        let circular = DatesCompat {
            circular_fft_lags: true,
            ..DatesCompat::CORRECTED
        };
        let both = DatesCompat {
            unscaled_last_fft_lag: true,
            ..circular
        };
        let diff = differing_bins(&run_mesh_test(1, circular), &run_mesh_test(1, both));
        assert_eq!(diff, vec![UNSCALED_BIN]);
    }

    #[test]
    fn test_leave_one_out_with_single_chromosome_is_empty() {
        let (geno, af, weights) = create_test_data();
        let (source1_freqs, source2_freqs) = source_freqs(&geno, &af);
        let ld = compute_ld_decay_by_chromosome(
            &geno,
            &af,
            "test",
            &weights,
            &source1_freqs,
            &source2_freqs,
            None,
        )
        .unwrap();
        assert_eq!(ld.chromosomes.len(), 1);
        assert_eq!(ld.chromosomes[0].chromosome, 1);
        assert_eq!(ld.chromosomes[0].n_snps, 50);
        assert!(ld.chromosomes[0]
            .leave_one_out
            .iter()
            .all(|b| b.n_pairs == 0));
    }

    #[test]
    fn test_ld_decay_missing_population() {
        let (geno, af, weights) = create_test_data();

        let n_snps = geno.n_snps();
        let source1_freqs: Vec<f64> = vec![0.5; n_snps];
        let source2_freqs: Vec<f64> = vec![0.5; n_snps];

        let result = compute_ld_decay(
            &geno,
            &af,
            "nonexistent",
            &weights,
            &source1_freqs,
            &source2_freqs,
            None,
        );
        assert!(result.is_err());
    }

    #[test]
    fn a_snp_without_a_genetic_position_is_an_error() {
        let (geno, af, weights) = create_test_data();
        let (s1, s2) = source_freqs(&geno, &af);
        let mut snps = geno.snps().to_vec();
        snps[7].position_morgans = f64::NAN;
        let data = Array2::from_shape_fn((geno.n_snps(), geno.n_samples()), |(i, j)| {
            geno.get(i, j) as u8
        });
        let geno = GenotypeMatrix::new(data, geno.samples().to_vec(), snps).unwrap();
        let err = compute_ld_decay(&geno, &af, "test", &weights, &s1, &s2, None).unwrap_err();
        assert!(err.to_string().contains("rs7"), "{err}");
    }
}
