//! The weighted LD decay table of DATES.
//!
//! Each test individual's genotypes are regressed on the two source
//! frequencies (its admixture proportion α), and the residual at each SNP is
//! weighted by w = p1 − p2. For pairs of SNPs on one chromosome, the
//! products of these values are summed into bins of genetic distance over
//! all test individuals; each bin reports their correlation. Two modes, as in
//! DATES: every pair of SNPs closer than `maxdis` (`qbin` = 0), or a mesh of
//! cells `binsize / qbin` wide whose per-lag sums come from FFT
//! cross-correlations (`qbin` > 0). Leave-one-chromosome-out tables feed the
//! jackknife.
//!
//! Provenance: written from the clean-room specifications
//! `docs/cleanroom/dates.md` (section 4) and `docs/cleanroom/dates-mesh.md`
//! (whose `mesh_cells` and `MeshLagSums` are kept as written from it),
//! Narasimhan et al. 2019 (Science, supplementary methods of DATES),
//! Chintalapati, Patterson & Moorjani 2022 (eLife), DATES's parameter-file
//! documentation (binsize, maxdis, qbin) and DATES v753 run as a black box
//! (`tests/golden/dates/v753/`).

use std::collections::BTreeMap;

use rustfft::num_complex::Complex;
use rustfft::FftPlanner;

use super::compat::DatesCompat;
use super::weights::informative;
use crate::error::{PopmixError, Result};
use crate::genotype_matrix::{AlleleFrequencies, GenotypeMatrix};
use crate::genotypes::{Chromosome, PopulationId, Snp};
use crate::options::check_bins;
use crate::threading::ordered_sum;

/// Mesh axis advance (Morgans) from one chromosome to the next.
const CHROMOSOME_GAP: f64 = 5.0;
/// Most mesh cells (`qbin` · maxdis / binsize) a run may ask for.
const MAX_MESH_CELLS: f64 = 1e7;

/// Options of the LD table (DATES's binsize, maxdis and qbin, in cM).
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct LdOptions {
    /// Bin width (cM); DATES's `binsize` is this / 100 Morgans.
    pub bin_width_cm: f64,
    /// Largest pair distance (cM); DATES's `maxdis` is this / 100 Morgans.
    pub max_distance_cm: f64,
    /// Bins with fewer pairs are left out of the table (0 keeps all).
    pub min_pairs_per_bin: usize,
    /// Mesh cells per bin; 0 counts every pair exactly.
    pub qbin: usize,
    /// DATES v753 behaviours to reproduce (mesh mode only here).
    pub compat: DatesCompat,
}

impl Default for LdOptions {
    fn default() -> Self {
        LdOptions {
            bin_width_cm: 0.1,
            max_distance_cm: 100.0,
            min_pairs_per_bin: 0,
            qbin: 0,
            compat: DatesCompat::CORRECTED,
        }
    }
}

impl LdOptions {
    /// Bin width and maximum distance positive and not too many bins; the
    /// mesh not too fine.
    pub fn validate(&self) -> Result<()> {
        check_bins(
            "bin_width_cm",
            self.bin_width_cm,
            "max_distance_cm",
            self.max_distance_cm,
        )?;
        let cells = self.qbin as f64 * self.max_distance_cm / self.bin_width_cm;
        if cells > MAX_MESH_CELLS {
            return Err(PopmixError::parameter(
                "qbin",
                format!(
                    "{} makes {cells:.0} mesh cells; at most 10000000 are allowed",
                    self.qbin
                ),
            ));
        }
        Ok(())
    }
}

/// One distance bin of the LD table.
#[derive(Debug, Clone)]
pub struct LdBin {
    /// Upper edge of the bin (cM), as DATES prints it.
    pub distance_cm: f64,
    /// Mean of y1·y2 over the bin's pairs.
    pub weighted_ld: f64,
    /// sqrt(pairs) · correlation.
    pub z_score: f64,
    /// Σy1y2 / sqrt(Σy1² · Σy2²) (no mean removed).
    pub correlation: f64,
    /// Σy1y2 / Σy1².
    pub regression_coef: f64,
    /// Number of pairs.
    pub n_pairs: usize,
    /// 1 / sqrt(pairs − 3), NaN for 3 pairs or fewer.
    pub std_err: f64,
}

/// The LD table of all chromosomes and, per chromosome, with it left out.
#[derive(Debug, Clone)]
pub struct LdDecay {
    /// The table over every chromosome.
    pub bins: Vec<LdBin>,
    /// One per autosome present, in increasing order.
    pub chromosomes: Vec<ChromosomeLd>,
}

/// The table with one chromosome left out.
#[derive(Debug, Clone)]
pub struct ChromosomeLd {
    /// The chromosome left out.
    pub chromosome: u8,
    /// Its SNPs in the input, before any filtering (the jackknife weight).
    pub n_snps: usize,
    /// The table of every other chromosome.
    pub leave_one_out: Vec<LdBin>,
}

/// Sums over the pairs of one bin: count and Σy1, Σy2, Σy1², Σy1y2, Σy2².
#[derive(Debug, Clone, Copy, Default)]
struct CorrBin {
    s0: f64,
    s1: f64,
    s2: f64,
    s11: f64,
    s12: f64,
    s22: f64,
}

impl CorrBin {
    fn add(&mut self, o: &CorrBin) {
        self.s0 += o.s0;
        self.s1 += o.s1;
        self.s2 += o.s2;
        self.s11 += o.s11;
        self.s12 += o.s12;
        self.s22 += o.s22;
    }

    fn minus(&self, o: &CorrBin) -> CorrBin {
        CorrBin {
            s0: self.s0 - o.s0,
            s1: self.s1 - o.s1,
            s2: self.s2 - o.s2,
            s11: self.s11 - o.s11,
            s12: self.s12 - o.s12,
            s22: self.s22 - o.s22,
        }
    }

    fn to_ld_bin(self, k: usize, bin_width_cm: f64) -> LdBin {
        let (weighted_ld, correlation, z_score) = if self.s0 < 0.5 {
            (0.0, 0.0, 0.0)
        } else {
            let v11 = self.s11 / self.s0;
            let v12 = self.s12 / self.s0;
            let v22 = self.s22 / self.s0;
            let corr = v12 / (v11 * v22 + 1e-20).sqrt();
            (v12, corr, self.s0.sqrt() * corr)
        };
        LdBin {
            distance_cm: (k + 1) as f64 * bin_width_cm,
            weighted_ld,
            z_score,
            correlation,
            regression_coef: self.s12 / (self.s11 + 1e-20),
            n_pairs: self.s0.round() as usize,
            std_err: if self.s0 > 3.0 {
                1.0 / (self.s0 - 3.0).sqrt()
            } else {
                f64::NAN
            },
        }
    }
}

fn add_bins(acc: &mut [CorrBin], other: &[CorrBin]) {
    acc.iter_mut().zip(other).for_each(|(a, b)| a.add(b));
}

/// The bins reported: those with at least `min_pairs` pairs.
fn report(sums: &[CorrBin], bin_width_cm: f64, min_pairs: usize) -> Vec<LdBin> {
    sums.iter()
        .enumerate()
        .filter(|(_, b)| b.s0 >= min_pairs as f64)
        .map(|(k, b)| b.to_ld_bin(k, bin_width_cm))
        .collect()
}

/// An autosome taking part: its SNP count and the SNPs used, by position.
struct ChromSnps {
    chromosome: u8,
    n_snps: usize,
    used: Vec<usize>,
}

/// The binned weighted LD over all chromosomes; see
/// [`compute_ld_decay_by_chromosome`].
pub fn compute_ld_decay(
    geno: &GenotypeMatrix,
    af: &AlleleFrequencies,
    test_pop: &str,
    weights: &[f64],
    source1_freqs: &[f64],
    source2_freqs: &[f64],
    options: &LdOptions,
) -> Result<Vec<LdBin>> {
    Ok(compute_ld_decay_by_chromosome(
        geno,
        af,
        test_pop,
        weights,
        source1_freqs,
        source2_freqs,
        options,
    )?
    .bins)
}

/// The binned weighted LD of the test population, over every autosome and
/// with each left out in turn.
///
/// `weights`, `source1_freqs` and `source2_freqs` hold one value per SNP of
/// `geno`; SNPs on autosomes 1-22 with an informative weight are used. `af`
/// is not used (kept for compatibility).
pub fn compute_ld_decay_by_chromosome(
    geno: &GenotypeMatrix,
    _af: &AlleleFrequencies,
    test_pop: &str,
    weights: &[f64],
    source1_freqs: &[f64],
    source2_freqs: &[f64],
    options: &LdOptions,
) -> Result<LdDecay> {
    let opts = options;
    opts.validate()?;
    let samples = geno.population_indices(&PopulationId::new(test_pop));
    if samples.is_empty() {
        return Err(PopmixError::PopulationNotFound(test_pop.to_string()));
    }
    let n = geno.n_snps();
    if weights.len() != n || source1_freqs.len() != n || source2_freqs.len() != n {
        return Err(PopmixError::DimensionMismatch {
            context: "DATES weights and source frequencies".to_string(),
            expected: format!("{n} of each, one per SNP"),
            actual: format!(
                "{}, {} and {}",
                weights.len(),
                source1_freqs.len(),
                source2_freqs.len()
            ),
        });
    }
    let binsize = opts.bin_width_cm / 100.0;
    let maxdis = opts.max_distance_cm / 100.0;
    let n_bins = (maxdis / binsize).round() as usize;
    let snps = geno.snps();

    let chroms = chromosomes_used(snps, weights, source1_freqs, source2_freqs)?;

    let per_chrom: Vec<Vec<CorrBin>> = if opts.qbin == 0 {
        let part = |s: usize| {
            let values = individual_values(
                geno,
                samples[s],
                &chroms,
                weights,
                source1_freqs,
                source2_freqs,
            );
            values
                .iter()
                .map(|v| exact_bins(v, snps, binsize, maxdis, n_bins))
                .collect::<Vec<_>>()
        };
        ordered_sum(
            samples.len(),
            || vec![vec![CorrBin::default(); n_bins]; chroms.len()],
            part,
            |acc, p| acc.iter_mut().zip(&p).for_each(|(a, b)| add_bins(a, b)),
        )
    } else {
        let qbin = opts.qbin;
        let cells = mesh_cells(snps, binsize, qbin);
        let diffmax = (qbin as f64 * maxdis / binsize).round() as usize;
        let part = |s: usize| {
            let values = individual_values(
                geno,
                samples[s],
                &chroms,
                weights,
                source1_freqs,
                source2_freqs,
            );
            let mut fft = MeshFft::new();
            values
                .iter()
                .map(|v| fft.lag_sums(v, &cells, diffmax, opts.compat))
                .collect::<Vec<_>>()
        };
        let sums = ordered_sum(
            samples.len(),
            || vec![MeshLagSums::new(diffmax); chroms.len()],
            part,
            |acc, p| acc.iter_mut().zip(&p).for_each(|(a, b)| a.merge(b)),
        );
        sums.iter()
            .map(|s| s.to_bins(binsize, qbin, n_bins))
            .collect()
    };

    let mut total = vec![CorrBin::default(); n_bins];
    for c in &per_chrom {
        add_bins(&mut total, c);
    }
    let chromosomes = chroms
        .iter()
        .zip(&per_chrom)
        .map(|(c, sums)| {
            let loo: Vec<CorrBin> = total.iter().zip(sums).map(|(t, s)| t.minus(s)).collect();
            ChromosomeLd {
                chromosome: c.chromosome,
                n_snps: c.n_snps,
                leave_one_out: report(&loo, opts.bin_width_cm, opts.min_pairs_per_bin),
            }
        })
        .collect();
    Ok(LdDecay {
        bins: report(&total, opts.bin_width_cm, opts.min_pairs_per_bin),
        chromosomes,
    })
}

/// Autosomes 1-22 present, in increasing order, each with its SNP count and
/// its used SNPs sorted by genetic position (stable). A used SNP without a
/// finite genetic position is an error.
fn chromosomes_used(
    snps: &[Snp],
    weights: &[f64],
    source1_freqs: &[f64],
    source2_freqs: &[f64],
) -> Result<Vec<ChromSnps>> {
    let mut by_chrom: BTreeMap<u8, ChromSnps> = BTreeMap::new();
    for (i, snp) in snps.iter().enumerate() {
        let c = match snp.chromosome {
            Chromosome::Autosome(c) if (1..=22).contains(&c) => c,
            _ => continue,
        };
        let entry = by_chrom.entry(c).or_insert_with(|| ChromSnps {
            chromosome: c,
            n_snps: 0,
            used: Vec::new(),
        });
        entry.n_snps += 1;
        if informative(weights[i], source1_freqs[i], source2_freqs[i]) {
            entry.used.push(i);
        }
    }
    let mut chroms: Vec<ChromSnps> = by_chrom.into_values().collect();
    for c in &mut chroms {
        if let Some(&i) = c
            .used
            .iter()
            .find(|&&i| !snps[i].position_morgans.is_finite())
        {
            return Err(PopmixError::data(format!(
                "SNP {} has no genetic position ({}); DATES needs one for every SNP",
                snps[i].id, snps[i].position_morgans
            )));
        }
        c.used.sort_by(|&a, &b| {
            snps[a]
                .position_morgans
                .total_cmp(&snps[b].position_morgans)
        });
    }
    Ok(chroms)
}

/// One individual's weighted residuals: per chromosome, (SNP index, y) at
/// its called used SNPs, in position order.
fn individual_values(
    geno: &GenotypeMatrix,
    sample: usize,
    chroms: &[ChromSnps],
    weights: &[f64],
    p1: &[f64],
    p2: &[f64],
) -> Vec<Vec<(usize, f64)>> {
    let called = |i: usize| -> Option<f64> {
        let g = geno.snp_raw(i)[sample];
        (g <= 2).then(|| f64::from(g) / 2.0)
    };
    let (mut num, mut den) = (0.0, 0.0);
    for c in chroms {
        for &i in &c.used {
            if let Some(x) = called(i) {
                let w = weights[i];
                num += (x - p2[i]) * w;
                den += w * w;
            }
        }
    }
    let alpha = if den <= 1e-15 { 0.5 } else { num / den };
    chroms
        .iter()
        .map(|c| {
            c.used
                .iter()
                .filter_map(|&i| {
                    called(i).map(|x| {
                        (
                            i,
                            (x - (alpha * p1[i] + (1.0 - alpha) * p2[i])) * weights[i],
                        )
                    })
                })
                .collect()
        })
        .collect()
}

/// Exact mode: every pair (i before j) closer than `maxdis`, y1 the earlier.
fn exact_bins(
    values: &[(usize, f64)],
    snps: &[Snp],
    binsize: f64,
    maxdis: f64,
    n_bins: usize,
) -> Vec<CorrBin> {
    let mut bins = vec![CorrBin::default(); n_bins];
    let pos: Vec<f64> = values
        .iter()
        .map(|&(i, _)| snps[i].position_morgans)
        .collect();
    for a in 0..values.len() {
        let y1 = values[a].1;
        for b in a + 1..values.len() {
            let dis = pos[b] - pos[a];
            if dis >= maxdis {
                break;
            }
            let k = (dis / binsize) as usize;
            let Some(bin) = bins.get_mut(k) else {
                continue;
            };
            let y2 = values[b].1;
            bin.s0 += 1.0;
            bin.s1 += y1;
            bin.s2 += y2;
            bin.s11 += y1 * y1;
            bin.s12 += y1 * y2;
            bin.s22 += y2 * y2;
        }
    }
    bins
}

/// An FFT planner (which caches its plans), kept across one individual's
/// chromosomes.
struct MeshFft {
    planner: FftPlanner<f64>,
}

impl MeshFft {
    fn new() -> Self {
        MeshFft {
            planner: FftPlanner::new(),
        }
    }

    /// One individual's lag sums on one chromosome (lags 0..=`diffmax`).
    fn lag_sums(
        &mut self,
        values: &[(usize, f64)],
        cells: &[usize],
        diffmax: usize,
        compat: DatesCompat,
    ) -> MeshLagSums {
        let mut out = MeshLagSums::new(diffmax);
        let (Some(lo), Some(hi)) = (
            values.iter().map(|&(i, _)| cells[i]).min(),
            values.iter().map(|&(i, _)| cells[i]).max(),
        ) else {
            return out;
        };
        let m = hi - lo + 1;
        if m <= 1 {
            return out;
        }
        let n = 1usize << ((m as f64).log2().ceil().round() as u32 + 1);
        let zero = Complex::new(0.0, 0.0);
        let (mut cn, mut cy, mut cq) = (vec![zero; n], vec![zero; n], vec![zero; n]);
        for &(i, y) in values {
            let c = cells[i] - lo;
            cn[c].re += 1.0;
            cy[c].re += y;
            cq[c].re += y * y;
        }
        let fwd = self.planner.plan_fft_forward(n);
        let inv = self.planner.plan_fft_inverse(n);
        fwd.process(&mut cn);
        fwd.process(&mut cy);
        fwd.process(&mut cq);
        let corr = |a: &[Complex<f64>], b: &[Complex<f64>]| -> Vec<f64> {
            let mut v: Vec<Complex<f64>> = a.iter().zip(b).map(|(x, y)| x.conj() * y).collect();
            inv.process(&mut v);
            v.iter().map(|c| c.re).collect()
        };
        let c_nn = corr(&cn, &cn);
        let c_ny = corr(&cn, &cy);
        let c_yy = corr(&cy, &cy);
        let c_nq = corr(&cn, &cq);

        let scale = 1.0 / n as f64;
        let last = if compat.circular_fft_lags {
            diffmax.min(n - 1)
        } else {
            diffmax.min(m - 1)
        };
        for d in 0..=last {
            let back = (n - d) % n;
            out.pairs[d] = c_nn[d] * scale;
            out.later[d] = c_ny[d] * scale;
            out.earlier[d] = c_ny[back] * scale;
            out.cross[d] = c_yy[d] * scale;
            out.later_sq[d] = c_nq[d] * scale;
            out.earlier_sq[d] = c_nq[back] * scale;
        }
        let unscaled = (n - 1).min(diffmax);
        if compat.unscaled_last_fft_lag && unscaled <= last {
            out.pairs[unscaled] = c_nn[unscaled];
            out.cross[unscaled] = c_yy[unscaled];
        }
        out
    }
}

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

#[cfg(test)]
mod tests;
