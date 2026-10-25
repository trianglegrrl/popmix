//! Principal component analysis as EIGENSOFT's smartpca computes it.
//!
//! For the PCA samples, each SNP is centered on its mean genotype and scaled
//! by 1/sqrt(p(1-p)), p half the mean (missing genotypes are 0 after
//! centering); SNPs monomorphic in the PCA samples are dropped. The
//! sample-by-sample matrix X X' is scaled to trace n - 1 and its eigenvalues
//! and eigenvectors are the PCA. Outlier samples (more than `outlier_sigma`
//! standard deviations out on one of the top `outlier_components`
//! eigenvectors) are removed and the PCA recomputed, up to
//! `outlier_iterations` times.
//!
//! Every sample, PCA or projected, gets coordinates from the SNP loadings:
//! either the dot product of its normalized genotypes with them (missing
//! genotypes count as 0), or, with `lsq_project`, the least-squares fit of
//! its observed genotypes, rescaled so the PCA samples' fits match their dot
//! products on average (smartpca `lsqproject: YES`). For PCA samples the dot
//! product is the eigenvector itself.
//!
//! Genotypes are counted as smartpca reads them: pseudo-haploid samples
//! (coded 0 or 2) are not treated differently.

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use ndarray::{linalg::general_mat_mul, Array2};
use rayon::prelude::*;

use super::tracy_widom::{tracy_widom_table, EigensoftCompat, TracyWidomRow};
use crate::error::{PopmixError, Result};
use crate::fstats::GenotypeRows;

/// Genotypes of a set of samples, two bits each (0, 1, 2; 3 missing), SNP by SNP.
#[derive(Debug, Clone)]
pub struct PackedGenotypes {
    n_samples: usize,
    n_snps: usize,
    bytes_per_snp: usize,
    data: Vec<u8>,
}

impl PackedGenotypes {
    /// Read `n_snps` rows of `n_samples` genotypes.
    pub fn read(rows: &mut dyn GenotypeRows, n_samples: usize, n_snps: usize) -> Result<Self> {
        let bytes_per_snp = n_samples.div_ceil(4);
        let mut data = vec![0u8; bytes_per_snp * n_snps];
        let mut row = vec![0u8; n_samples];
        for snp in 0..n_snps {
            rows.next_row(&mut row)?;
            let out = &mut data[snp * bytes_per_snp..(snp + 1) * bytes_per_snp];
            for (k, &g) in row.iter().enumerate() {
                out[k / 4] |= g.min(3) << (2 * (k % 4));
            }
        }
        Ok(PackedGenotypes {
            n_samples,
            n_snps,
            bytes_per_snp,
            data,
        })
    }

    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    pub fn n_snps(&self) -> usize {
        self.n_snps
    }

    /// The given SNPs (by index), in that order.
    pub fn select_snps(&self, snps: &[usize]) -> PackedGenotypes {
        let mut data = Vec::with_capacity(snps.len() * self.bytes_per_snp);
        for &s in snps {
            data.extend_from_slice(
                &self.data[s * self.bytes_per_snp..(s + 1) * self.bytes_per_snp],
            );
        }
        PackedGenotypes {
            n_samples: self.n_samples,
            n_snps: snps.len(),
            bytes_per_snp: self.bytes_per_snp,
            data,
        }
    }

    /// The given samples (by index), in that order.
    pub fn select_samples(&self, samples: &[usize]) -> PackedGenotypes {
        let bytes_per_snp = samples.len().div_ceil(4);
        let mut data = vec![0u8; bytes_per_snp * self.n_snps];
        for snp in 0..self.n_snps {
            let out = &mut data[snp * bytes_per_snp..(snp + 1) * bytes_per_snp];
            for (k, &s) in samples.iter().enumerate() {
                out[k / 4] |= self.get(snp, s) << (2 * (k % 4));
            }
        }
        PackedGenotypes {
            n_samples: samples.len(),
            n_snps: self.n_snps,
            bytes_per_snp,
            data,
        }
    }

    /// Whether some sample has a genotype at the SNP.
    pub fn any_called(&self, snp: usize) -> bool {
        (0..self.n_samples).any(|s| self.get(snp, s) <= 2)
    }

    #[inline]
    pub fn get(&self, snp: usize, sample: usize) -> u8 {
        (self.data[snp * self.bytes_per_snp + sample / 4] >> (2 * (sample % 4))) & 3
    }
}

/// Whether a sample defines the PCA or is projected onto it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleRole {
    Pca,
    Project,
}

/// smartpca's parameters, with its defaults.
#[derive(Debug, Clone)]
pub struct PcaOptions {
    /// Eigenvectors to report (`numoutevec`).
    pub n_components: usize,
    /// Rounds of outlier removal (`numoutlieriter`); 0 for none.
    pub outlier_iterations: usize,
    /// Eigenvectors examined for outliers (`numoutlierevec`).
    pub outlier_components: usize,
    /// Outlier threshold in standard deviations (`outliersigmathresh`).
    pub outlier_sigma: f64,
    /// Project by least squares on observed genotypes (`lsqproject`).
    pub lsq_project: bool,
    /// Fewest eigenvalues for a Tracy-Widom statistic (`nostatslim`).
    pub min_tw_eigs: usize,
    pub compat: EigensoftCompat,
}

impl Default for PcaOptions {
    fn default() -> Self {
        PcaOptions {
            n_components: 10,
            outlier_iterations: 5,
            outlier_components: 10,
            outlier_sigma: 6.0,
            lsq_project: false,
            min_tw_eigs: 10,
            compat: EigensoftCompat::CORRECTED,
        }
    }
}

/// A sample removed as an outlier.
#[derive(Debug, Clone, PartialEq)]
pub struct Outlier {
    pub sample: usize,
    /// Round of outlier removal, from 1.
    pub iteration: usize,
    /// Eigenvector (from 0) on which it was first found out.
    pub component: usize,
    /// Its standardized score on that eigenvector.
    pub score: f64,
}

#[derive(Debug, Clone)]
pub struct PcaResult {
    /// Every eigenvalue of the scaled matrix, decreasing.
    pub eigenvalues: Vec<f64>,
    /// Samples with coordinates (index into the input samples) and their
    /// coordinates on the first `n_components` eigenvectors; outliers and
    /// samples without data are left out.
    pub samples: Vec<usize>,
    pub coords: Vec<Vec<f64>>,
    pub tracy_widom: Vec<TracyWidomRow>,
    pub outliers: Vec<Outlier>,
    /// Samples left out for lack of data (none valid, or, projected by least
    /// squares, no more observed SNPs than components).
    pub no_data: Vec<usize>,
    /// PCA samples in the final pass.
    pub n_pca_samples: usize,
    /// SNPs in the final PCA (polymorphic in its samples).
    pub n_snps_used: usize,
    /// SNPs dropped as monomorphic or without data, over all passes.
    pub n_snps_killed: usize,
}

/// One SNP's normalization over the PCA samples: the values of genotypes
/// 0, 1, 2 (missing is 0), or None if the SNP is dropped.
#[derive(Debug, Clone, Copy)]
struct ColumnNorm {
    values: [f64; 3],
    /// smartpca `xmean` and `xfancy`, used to normalize any sample's row.
    mean: f64,
    fancy: f64,
    alive: bool,
}

fn column_norm(geno: &PackedGenotypes, snp: usize, used: &[usize]) -> ColumnNorm {
    let (mut c0, mut c1, mut nmiss) = (0u64, 0u64, 0usize);
    for &s in used {
        let g = geno.get(snp, s);
        if g > 2 {
            nmiss += 1;
        } else {
            c0 += g as u64;
            c1 += 2 - g as u64;
        }
    }
    if nmiss == used.len() {
        return ColumnNorm {
            values: [0.0; 3],
            mean: 0.0,
            fancy: 0.0,
            alive: false,
        };
    }
    let ymean = c0 as f64 / (used.len() - nmiss) as f64;
    let p = 0.5 * ymean;
    let y = p * (1.0 - p);
    let fancy = if y > 0.0 { 1.0 / y.sqrt() } else { 1.0 };
    ColumnNorm {
        values: [-ymean * fancy, (1.0 - ymean) * fancy, (2.0 - ymean) * fancy],
        mean: ymean * fancy,
        fancy,
        alive: c0.min(c1) >= 1,
    }
}

/// The normalized genotype of any sample at a SNP (smartpca `fixxrow`).
#[inline]
fn row_value(g: u8, norm: &ColumnNorm) -> f64 {
    if g > 2 {
        0.0
    } else {
        g as f64 * norm.fancy - norm.mean
    }
}

const BLOCK: usize = 512;

/// SNPs per group in the parallel sums. Groups are fixed, and their partial
/// sums are added in group order, so results do not depend on the thread count.
const GROUP: usize = 16 * BLOCK;

/// `zero() + part(0) + part(1) + ... + part(n_groups - 1)`, in that order,
/// computing `threads` parts at a time.
fn ordered_sum<T: Send>(
    n_groups: usize,
    zero: impl Fn() -> T,
    part: impl Fn(usize) -> T + Sync,
    add: impl Fn(&mut T, T),
) -> T {
    let mut acc = zero();
    let batch = rayon::current_num_threads().max(1);
    for start in (0..n_groups).step_by(batch) {
        let parts: Vec<T> = (start..(start + batch).min(n_groups))
            .into_par_iter()
            .map(&part)
            .collect();
        for p in parts {
            add(&mut acc, p);
        }
    }
    acc
}

/// X X' over the given SNPs, X the normalized genotypes of `used`.
fn cross_product(
    geno: &PackedGenotypes,
    snps: &[usize],
    norms: &[ColumnNorm],
    used: &[usize],
) -> Array2<f64> {
    let n = used.len();
    let part = |group: usize| {
        let lo = group * GROUP;
        let hi = (lo + GROUP).min(snps.len());
        let mut acc = Array2::<f64>::zeros((n, n));
        for (chunk, chunk_norms) in snps[lo..hi].chunks(BLOCK).zip(norms[lo..hi].chunks(BLOCK)) {
            let live: Vec<(usize, &ColumnNorm)> = chunk
                .iter()
                .zip(chunk_norms)
                .filter(|(_, c)| c.alive)
                .map(|(&s, c)| (s, c))
                .collect();
            let mut x = Array2::<f64>::zeros((n, live.len()));
            for (col, (snp, norm)) in live.iter().enumerate() {
                for (row, &s) in used.iter().enumerate() {
                    let g = geno.get(*snp, s);
                    x[[row, col]] = if g > 2 { 0.0 } else { norm.values[g as usize] };
                }
            }
            general_mat_mul(1.0, &x, &x.t(), 1.0, &mut acc);
        }
        acc
    };
    ordered_sum(
        snps.len().div_ceil(GROUP),
        || Array2::<f64>::zeros((n, n)),
        part,
        |a, b| *a += &b,
    )
}

/// Eigenvalues (decreasing) and unit eigenvectors (columns, same order).
fn eigen(m: &Array2<f64>) -> (Vec<f64>, DMatrix<f64>) {
    let n = m.nrows();
    let dm = DMatrix::from_fn(n, n, |i, j| m[[i, j]]);
    let eig = SymmetricEigen::new(dm);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| eig.eigenvalues[b].total_cmp(&eig.eigenvalues[a]));
    let values = order.iter().map(|&k| eig.eigenvalues[k]).collect();
    let vectors = DMatrix::from_fn(n, n, |i, j| eig.eigenvectors[(i, order[j])]);
    (values, vectors)
}

/// smartpca `ridoutlier` (outliermode 0): rows more than `thresh` standard
/// deviations from the mean on any of the first `neigs` eigenvectors, with
/// the first eigenvector and score each was found on.
fn find_outliers(evecs: &DMatrix<f64>, neigs: usize, thresh: f64) -> Vec<(usize, usize, f64)> {
    let n = evecs.nrows();
    if n < 3 {
        return Vec::new();
    }
    let mut found: Vec<Option<(usize, f64)>> = vec![None; n];
    for i in 0..neigs {
        let w = evecs.column(i);
        let mean = w.sum() / n as f64;
        let sd = (w.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / n as f64).sqrt();
        for (j, f) in found.iter_mut().enumerate() {
            let z = (w[j] - mean) / sd;
            if z.abs() > thresh && f.is_none() {
                *f = Some((i, z));
            }
        }
    }
    found
        .into_iter()
        .enumerate()
        .filter_map(|(j, f)| f.map(|(i, z)| (j, i, z)))
        .collect()
}

/// PCA of the samples with role `Pca`, projecting the others, as smartpca.
///
/// `geno` holds only the SNPs to consider (smartpca uses autosomes 1-22).
pub fn smartpca(
    geno: &PackedGenotypes,
    roles: &[SampleRole],
    options: &PcaOptions,
) -> Result<PcaResult> {
    let n_samples = geno.n_samples();
    if roles.len() != n_samples {
        return Err(PopmixError::DimensionMismatch {
            context: "smartpca sample roles".to_string(),
            expected: format!("{n_samples} roles"),
            actual: format!("{}", roles.len()),
        });
    }

    // SNPs with at most one valid genotype in the PCA samples ("nodata"),
    // then samples with no valid genotype at any remaining SNP.
    let pca: Vec<usize> = (0..n_samples)
        .filter(|&s| roles[s] == SampleRole::Pca)
        .collect();
    let mut snp_alive: Vec<bool> = (0..geno.n_snps())
        .into_par_iter()
        .map(|snp| pca.iter().filter(|&&s| geno.get(snp, s) <= 2).count() > 1)
        .collect();
    let mut ignored = vec![false; n_samples];
    let mut no_data = Vec::new();
    for (s, ign) in ignored.iter_mut().enumerate() {
        if !(0..geno.n_snps()).any(|snp| snp_alive[snp] && geno.get(snp, s) <= 2) {
            *ign = true;
            no_data.push(s);
        }
    }

    let numoutiter = options.outlier_iterations.max(0) + 1;
    let mut outliers = Vec::new();
    let mut n_killed = 0usize;
    let mut pass = 0;
    let (lambda, evecs, used, snps, norms) = loop {
        pass += 1;
        let used: Vec<usize> = pca.iter().copied().filter(|&s| !ignored[s]).collect();
        if used.len() < 2 {
            return Err(PopmixError::Numerical(format!(
                "{} PCA samples with data",
                used.len()
            )));
        }
        // loadsnpx: SNPs not dropped, with a valid genotype in some sample
        let snps: Vec<usize> = (0..geno.n_snps())
            .filter(|&snp| {
                snp_alive[snp] && (0..n_samples).any(|s| !ignored[s] && geno.get(snp, s) <= 2)
            })
            .collect();
        let norms: Vec<ColumnNorm> = snps
            .par_iter()
            .map(|&snp| column_norm(geno, snp, &used))
            .collect();
        for (&snp, norm) in snps.iter().zip(&norms) {
            if !norm.alive {
                snp_alive[snp] = false;
                n_killed += 1;
            }
        }
        let mut xtx = cross_product(geno, &snps, &norms, &used);
        let trace: f64 =
            (0..used.len()).map(|i| xtx[[i, i]]).sum::<f64>() / (used.len() - 1) as f64;
        if !trace.is_finite() || trace <= 0.0 {
            return Err(PopmixError::Numerical(
                "the genotype matrix has zero trace (no data)".to_string(),
            ));
        }
        xtx.mapv_inplace(|x| x / trace);
        let (lambda, evecs) = eigen(&xtx);

        if pass >= numoutiter || options.outlier_iterations == 0 {
            break (lambda, evecs, used, snps, norms);
        }
        let neigs = options.outlier_components.min(used.len() - 1);
        let bad = find_outliers(&evecs, neigs, options.outlier_sigma);
        if bad.is_empty() {
            break (lambda, evecs, used, snps, norms);
        }
        for (row, component, score) in bad {
            ignored[used[row]] = true;
            outliers.push(Outlier {
                sample: used[row],
                iteration: pass,
                component,
                score,
            });
        }
    };

    let n_used = used.len();
    let k = options.n_components.min(n_used).min(snps.len());
    // SNP loadings: ffvecs[j][col] = sum over PCA samples of 10 v_j * x
    let v: Vec<Vec<f64>> = (0..n_used)
        .map(|row| (0..k).map(|j| 10.0 * evecs[(row, j)]).collect())
        .collect();
    let by_snp: Vec<Vec<f64>> = snps
        .par_iter()
        .zip(&norms)
        .map(|(&snp, norm)| {
            let mut acc = vec![0.0; k];
            if norm.alive {
                for (row, &s) in used.iter().enumerate() {
                    let g = geno.get(snp, s);
                    if g <= 2 {
                        let x = norm.values[g as usize];
                        for (a, vj) in acc.iter_mut().zip(&v[row]) {
                            *a += vj * x;
                        }
                    }
                }
            }
            acc
        })
        .collect();
    let flat: Vec<f64> = by_snp.into_iter().flatten().collect();

    // One pass over the SNPs for every sample left: the dot product of its
    // normalized genotypes with the loadings and, for least squares, the
    // right-hand side and the normal matrix over its missing SNPs.
    let mut kept: Vec<usize> = (0..n_samples).filter(|&s| !ignored[s]).collect();
    let sums = sample_sums(geno, &kept, &snps, &norms, &flat, k, options.lsq_project);
    let row_of_used: Vec<usize> = used
        .iter()
        .map(|s| kept.binary_search(s).expect("PCA sample kept"))
        .collect();
    let scale: Vec<f64> = (0..k)
        .map(|j| {
            let ss: f64 = row_of_used
                .iter()
                .map(|&r| sums.dot[r * k + j].powi(2))
                .sum();
            1.0 / ss.sqrt()
        })
        .collect();
    let mut bcoeffs: Vec<Vec<f64>> = (0..kept.len())
        .map(|r| (0..k).map(|j| sums.dot[r * k + j] * scale[j]).collect())
        .collect();

    let mut lsq: Vec<Option<Vec<f64>>> = vec![None; n_samples];
    if options.lsq_project {
        let all = normal_matrix(&flat, k, &mut (0..snps.len()));
        let fits: Vec<Option<Vec<f64>>> = (0..kept.len())
            .into_par_iter()
            .map(|r| {
                let observed = snps.len() - sums.n_missing[r];
                if observed <= k {
                    return None;
                }
                let mut co = DMatrix::<f64>::zeros(k, k);
                for a in 0..k {
                    for b in a..k {
                        let x = (all[(a, b)] - sums.missing[r * k * k + a * k + b])
                            * scale[a]
                            * scale[b];
                        co[(a, b)] = x;
                        co[(b, a)] = x;
                    }
                }
                let rr = DVector::from_fn(k, |j, _| sums.dot[r * k + j] * scale[j]);
                Some(co.cholesky()?.solve(&rr).iter().copied().collect())
            })
            .collect();
        for (r, fit) in fits.into_iter().enumerate() {
            let s = kept[r];
            match fit {
                Some(f) => lsq[s] = Some(f),
                None => {
                    ignored[s] = true;
                    no_data.push(s);
                }
            }
        }
        let mut r = 0;
        bcoeffs.retain(|_| {
            r += 1;
            !ignored[kept[r - 1]]
        });
        kept.retain(|&s| !ignored[s]);
    }

    let coords = if options.lsq_project {
        // rescale each component so the PCA samples' fits match their
        // projections: y = <a, b> / <a, a> over the PCA samples
        let pos = |s: usize| kept.binary_search(&s).ok();
        let factors: Vec<f64> = (0..k)
            .map(|j| {
                let (mut ab, mut aa) = (0.0, 0.0);
                for &s in &used {
                    let a = lsq[s].as_ref().map_or(0.0, |f| f[j]);
                    let b = pos(s).map_or(0.0, |p| bcoeffs[p][j]);
                    ab += a * b;
                    aa += a * a;
                }
                ab / aa
            })
            .collect();
        kept.iter()
            .map(|&s| {
                let a = lsq[s].as_ref().expect("fitted");
                a.iter().zip(&factors).map(|(x, f)| x * f).collect()
            })
            .collect()
    } else {
        bcoeffs
    };

    no_data.sort_unstable();
    Ok(PcaResult {
        tracy_widom: tracy_widom_table(&lambda, options.min_tw_eigs, options.compat),
        eigenvalues: lambda,
        samples: kept,
        coords,
        outliers,
        no_data,
        n_pca_samples: n_used,
        n_snps_used: norms.iter().filter(|c| c.alive).count(),
        n_snps_killed: n_killed,
    })
}

/// Per-sample sums over the SNPs (rows follow `samples`, `k` per sample):
/// the dot product of the normalized genotypes with the loadings (missing
/// genotypes are 0, so this is also the least-squares right-hand side), the
/// upper triangle of l l' summed over missing SNPs (for least squares), and
/// the number of missing SNPs. `flat` holds the loadings SNP by SNP.
struct SampleSums {
    dot: Vec<f64>,
    missing: Vec<f64>,
    n_missing: Vec<usize>,
}

fn sample_sums(
    geno: &PackedGenotypes,
    samples: &[usize],
    snps: &[usize],
    norms: &[ColumnNorm],
    flat: &[f64],
    k: usize,
    lsq: bool,
) -> SampleSums {
    let n = samples.len();
    let zero = || SampleSums {
        dot: vec![0.0; n * k],
        missing: if lsq {
            vec![0.0; n * k * k]
        } else {
            Vec::new()
        },
        n_missing: vec![0; n],
    };
    let part = |group: usize| {
        let cols = group * GROUP..((group + 1) * GROUP).min(snps.len());
        let mut acc = zero();
        for col in cols {
            let (snp, norm) = (snps[col], &norms[col]);
            let l = &flat[col * k..(col + 1) * k];
            for (r, &s) in samples.iter().enumerate() {
                let g = geno.get(snp, s);
                if g > 2 {
                    acc.n_missing[r] += 1;
                    if lsq {
                        let m = &mut acc.missing[r * k * k..(r + 1) * k * k];
                        for a in 0..k {
                            for b in a..k {
                                m[a * k + b] += l[a] * l[b];
                            }
                        }
                    }
                    continue;
                }
                let x = row_value(g, norm);
                if x != 0.0 {
                    for (d, lj) in acc.dot[r * k..(r + 1) * k].iter_mut().zip(l) {
                        *d += x * lj;
                    }
                }
            }
        }
        acc
    };
    ordered_sum(snps.len().div_ceil(GROUP), zero, part, |a, b| {
        a.dot.iter_mut().zip(&b.dot).for_each(|(x, y)| *x += y);
        a.missing
            .iter_mut()
            .zip(&b.missing)
            .for_each(|(x, y)| *x += y);
        a.n_missing
            .iter_mut()
            .zip(&b.n_missing)
            .for_each(|(x, y)| *x += y);
    })
}

/// Sum of l l' over the given SNP columns, upper triangle (`flat` holds the
/// loadings SNP by SNP).
fn normal_matrix(flat: &[f64], k: usize, cols: &mut dyn Iterator<Item = usize>) -> DMatrix<f64> {
    let mut co = DMatrix::<f64>::zeros(k, k);
    for col in cols {
        let l = &flat[col * k..(col + 1) * k];
        for a in 0..k {
            for b in a..k {
                co[(a, b)] += l[a] * l[b];
            }
        }
    }
    co
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fstats::MatrixRows;
    use crate::genotype_matrix::GenotypeMatrix;
    use crate::genotypes::{Chromosome, Sample, Snp};

    /// Two groups of samples whose allele frequencies differ at every SNP.
    fn two_groups(n_per: usize, n_snps: usize) -> PackedGenotypes {
        let samples: Vec<Sample> = (0..2 * n_per)
            .map(|i| Sample::new(format!("s{i}"), "P"))
            .collect();
        let snps: Vec<Snp> = (0..n_snps)
            .map(|i| {
                Snp::new(
                    format!("rs{i}"),
                    Chromosome::Autosome(1),
                    i as u32,
                    0.0,
                    'A',
                    'G',
                )
            })
            .collect();
        let data = Array2::from_shape_fn((n_snps, 2 * n_per), |(snp, s)| {
            let h = (snp * 31 + s * 17 + snp * s) % 7;
            let group = s / n_per;
            let p = if (snp + group).is_multiple_of(2) {
                0.8
            } else {
                0.2
            };
            let u = h as f64 / 7.0;
            if u < p * p {
                2
            } else if u < p * p + 2.0 * p * (1.0 - p) {
                1
            } else {
                0
            }
        });
        let m = GenotypeMatrix::new(data, samples, snps).unwrap();
        let mut rows = MatrixRows::new(&m, (0..n_snps).collect(), (0..2 * n_per).collect());
        PackedGenotypes::read(&mut rows, 2 * n_per, n_snps).unwrap()
    }

    #[test]
    fn packed_genotypes_round_trip() {
        let g = two_groups(5, 30);
        assert_eq!(g.n_samples(), 10);
        let total: usize = (0..30)
            .flat_map(|snp| (0..10).map(move |s| (snp, s)))
            .map(|(snp, s)| g.get(snp, s) as usize)
            .sum();
        assert!(total > 0);
    }

    #[test]
    fn the_first_component_separates_the_groups() {
        let g = two_groups(10, 400);
        let roles = vec![SampleRole::Pca; 20];
        let opts = PcaOptions {
            outlier_iterations: 0,
            ..Default::default()
        };
        let out = smartpca(&g, &roles, &opts).unwrap();
        // trace n - 1
        assert!((out.eigenvalues.iter().sum::<f64>() - 19.0).abs() < 1e-9);
        let pc1: Vec<f64> = out.coords.iter().map(|c| c[0]).collect();
        let (a, b) = pc1.split_at(10);
        assert!(a.iter().all(|&x| x.signum() == a[0].signum()));
        assert!(b.iter().all(|&x| x.signum() == -a[0].signum()));
        // PCA samples' coordinates are the unit eigenvectors
        assert!((pc1.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn results_do_not_depend_on_the_thread_count() {
        // enough SNPs for several parallel chunks in both sums
        let g = two_groups(10, 40_000);
        let mut roles = vec![SampleRole::Pca; 20];
        roles[3] = SampleRole::Project;
        let opts = PcaOptions {
            outlier_iterations: 0,
            ..Default::default()
        };
        let run = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| smartpca(&g, &roles, &opts).unwrap())
        };
        let one = run(1);
        for threads in [2, 3, 8] {
            let many = run(threads);
            assert_eq!(one.eigenvalues, many.eigenvalues, "{threads} threads");
            assert_eq!(one.coords, many.coords, "{threads} threads");
        }
    }

    #[test]
    fn projecting_a_pca_sample_reproduces_its_coordinates() {
        let g = two_groups(10, 400);
        let roles = vec![SampleRole::Pca; 20];
        let opts = PcaOptions {
            outlier_iterations: 0,
            lsq_project: true,
            ..Default::default()
        };
        let lsq = smartpca(&g, &roles, &opts).unwrap();
        let opts = PcaOptions {
            outlier_iterations: 0,
            ..Default::default()
        };
        let dot = smartpca(&g, &roles, &opts).unwrap();
        // without missing data the least-squares fit is the projection
        for (a, b) in lsq.coords.iter().zip(&dot.coords) {
            for (x, y) in a.iter().zip(b) {
                assert!((x - y).abs() < 1e-9, "{x} != {y}");
            }
        }
    }

    #[test]
    fn projected_samples_do_not_change_the_pca() {
        let g = two_groups(10, 400);
        let mut roles = vec![SampleRole::Pca; 20];
        let opts = PcaOptions {
            outlier_iterations: 0,
            ..Default::default()
        };
        let all = smartpca(&g, &roles, &opts).unwrap();
        roles[3] = SampleRole::Project;
        let projected = smartpca(&g, &roles, &opts).unwrap();
        assert_eq!(projected.n_pca_samples, 19);
        assert_eq!(projected.samples.len(), 20);
        assert!((projected.eigenvalues[0] - all.eigenvalues[0]).abs() > 1e-6);
        // the projected sample keeps its group's side of PC1
        let side = |r: &PcaResult, s: usize| {
            r.coords[r.samples.iter().position(|&x| x == s).unwrap()][0].signum()
        };
        assert_eq!(side(&projected, 3), side(&projected, 0));
    }

    #[test]
    fn selecting_snps_keeps_their_genotypes() {
        let g = two_groups(3, 20);
        let picked = g.select_samples(&[4, 1]);
        assert_eq!(picked.n_samples(), 2);
        assert!((0..20).all(|j| picked.get(j, 0) == g.get(j, 4) && picked.get(j, 1) == g.get(j, 1)));
        let sub = g.select_snps(&[5, 2]);
        assert_eq!(sub.n_snps(), 2);
        for s in 0..6 {
            assert_eq!(sub.get(0, s), g.get(5, s));
            assert_eq!(sub.get(1, s), g.get(2, s));
        }
        assert!(g.any_called(0));
    }

    #[test]
    fn mismatched_roles_are_an_error() {
        let g = two_groups(3, 20);
        assert!(smartpca(&g, &[SampleRole::Pca; 5], &PcaOptions::default()).is_err());
    }
}
