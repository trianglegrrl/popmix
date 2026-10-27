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

use nalgebra::{DMatrix, DVector};
use ndarray::{linalg::general_mat_mul, s, Array2, Axis, ShapeBuilder};
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
#[non_exhaustive]
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
#[non_exhaustive]
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
            c0 += u64::from(g);
            c1 += 2 - u64::from(g);
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
        f64::from(g) * norm.fancy - norm.mean
    }
}

/// Samples per side of a tile of X X': small enough to give every thread
/// several tiles, chosen from the number of samples alone.
fn tile_size(n: usize) -> usize {
    if n < 2048 {
        128
    } else {
        256
    }
}

/// Normalized genotypes held at once by `cross_product`, in values.
const PANEL_VALUES: usize = 1 << 22;

/// X X' over the given SNPs, X the normalized genotypes of `used`.
///
/// The live SNPs are taken a panel at a time; each panel is decoded into a
/// samples x SNPs matrix and added to the upper triangle of X X' tile by
/// tile, the tiles in parallel. Each tile adds its panels in SNP order, and
/// the panel and tile sizes depend only on the number of samples, so the
/// result does not depend on the thread count.
fn cross_product(
    geno: &PackedGenotypes,
    snps: &[usize],
    norms: &[ColumnNorm],
    used: &[usize],
) -> Array2<f64> {
    let n = used.len();
    let live: Vec<(usize, &ColumnNorm)> = snps
        .iter()
        .zip(norms)
        .filter(|(_, c)| c.alive)
        .map(|(&s, c)| (s, c))
        .collect();
    let panel = (PANEL_VALUES / n.max(1)).clamp(256, 4096);
    let tile = tile_size(n);
    let starts: Vec<usize> = (0..n).step_by(tile).collect();
    let mut tiles: Vec<((usize, usize), Array2<f64>)> = starts
        .iter()
        .enumerate()
        .flat_map(|(i, &a)| starts[i..].iter().map(move |&b| (a, b)))
        .map(|(a, b)| {
            let shape = ((a + tile).min(n) - a, (b + tile).min(n) - b);
            ((a, b), Array2::<f64>::zeros(shape))
        })
        .collect();
    // samples x SNPs, column-major so that each SNP decodes contiguously
    let mut x = Array2::<f64>::zeros((n, panel).f());
    for chunk in live.chunks(panel) {
        let width = chunk.len();
        x.slice_mut(s![.., ..width])
            .axis_iter_mut(Axis(1))
            .into_par_iter()
            .zip(chunk)
            .for_each(|(mut col, (snp, norm))| {
                for (v, &s) in col.iter_mut().zip(used) {
                    let g = geno.get(*snp, s);
                    *v = if g > 2 { 0.0 } else { norm.values[g as usize] };
                }
            });
        let x = x.slice(s![.., ..width]);
        tiles.par_iter_mut().for_each(|((a, b), acc)| {
            let rows = x.slice(s![*a..*a + acc.nrows(), ..]);
            let cols = x.slice(s![*b..*b + acc.ncols(), ..]);
            general_mat_mul(1.0, &rows, &cols.t(), 1.0, acc);
        });
    }
    let mut out = Array2::<f64>::zeros((n, n));
    for ((a, b), acc) in &tiles {
        for ((i, j), &v) in acc.indexed_iter() {
            out[[a + i, b + j]] = v;
            out[[b + j, a + i]] = v;
        }
    }
    // a diagonal tile is computed whole: keep its upper triangle
    for i in 0..n {
        for j in i + 1..n {
            out[[j, i]] = out[[i, j]];
        }
    }
    out
}

/// Eigenvalues (decreasing) and unit eigenvectors (columns, same order).
///
/// faer's solver runs sequentially, so the result does not depend on the
/// thread count.
fn eigen(m: &Array2<f64>) -> (Vec<f64>, DMatrix<f64>) {
    use faer::linalg::evd::{self_adjoint_evd, self_adjoint_evd_scratch, ComputeEigenvectors};
    use faer::{diag::Diag, dyn_stack::MemBuffer, dyn_stack::MemStack, Mat, Par};
    let n = m.nrows();
    let a = Mat::<f64>::from_fn(n, n, |i, j| m[[i, j]]);
    let mut s = Diag::<f64>::zeros(n);
    let mut u = Mat::<f64>::zeros(n, n);
    let scratch =
        self_adjoint_evd_scratch::<f64>(n, ComputeEigenvectors::Yes, Par::Seq, Default::default());
    self_adjoint_evd(
        a.as_ref(),
        s.as_mut(),
        Some(u.as_mut()),
        Par::Seq,
        MemStack::new(&mut MemBuffer::new(scratch)),
        Default::default(),
    )
    .expect("eigendecomposition of a finite symmetric matrix");
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| s[b].total_cmp(&s[a]));
    let values = order.iter().map(|&k| s[k]).collect();
    let vectors = DMatrix::from_fn(n, n, |i, j| u[(i, order[j])]);
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

    let numoutiter = options.outlier_iterations + 1;
    let mut outliers = Vec::new();
    let mut n_killed = 0usize;
    let mut pass = 0;
    let (lambda, evecs, used, snps, norms) = loop {
        pass += 1;
        let used: Vec<usize> = pca.iter().copied().filter(|&s| !ignored[s]).collect();
        if used.len() < 2 {
            return Err(PopmixError::too_few(
                "PCA",
                "samples with data",
                2,
                used.len(),
            ));
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
            return Err(PopmixError::data(
                "the genotype matrix has zero trace (no data)",
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
    // SNP loadings, ffvecs[j][col] = sum over PCA samples of 10 v_j * x,
    // then one pass over the SNPs for every sample left: the dot product of
    // its normalized genotypes with the loadings and, for least squares, the
    // right-hand side and the normal matrix over its missing SNPs.
    let mut kept: Vec<usize> = (0..n_samples).filter(|&s| !ignored[s]).collect();
    let mut v = Array2::<f64>::zeros((kept.len(), k));
    for (row, s) in used.iter().enumerate() {
        let r = kept.binary_search(s).expect("PCA sample kept");
        for j in 0..k {
            v[[r, j]] = 10.0 * evecs[(row, j)];
        }
    }
    let (flat, sums) = project(geno, &kept, &snps, &norms, &v, options.lsq_project);
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
/// the number of missing SNPs.
struct SampleSums {
    dot: Vec<f64>,
    missing: Vec<f64>,
    n_missing: Vec<usize>,
}

/// SNPs per piece of a panel in `project`, and samples per piece of its rows.
const PIECE: usize = 256;

/// The SNP loadings (SNP by SNP, `k` each; `v` holds 10 times the
/// eigenvectors on the rows of the PCA samples and 0 elsewhere) and the
/// per-sample sums over `samples`.
///
/// The live SNPs are taken a panel at a time and decoded once; the loadings
/// of a panel are its genotypes' product with `v`, and the dot products
/// gain the genotypes' product with those loadings. Pieces of SNPs and of
/// samples run in parallel, each adding its panels in SNP order, so the
/// result does not depend on the thread count.
fn project(
    geno: &PackedGenotypes,
    samples: &[usize],
    snps: &[usize],
    norms: &[ColumnNorm],
    v: &Array2<f64>,
    lsq: bool,
) -> (Vec<f64>, SampleSums) {
    let (n, k) = (samples.len(), v.ncols());
    let live: Vec<usize> = (0..snps.len()).filter(|&c| norms[c].alive).collect();
    let panel = (PANEL_VALUES / n.max(1)).clamp(PIECE, 16 * PIECE) / PIECE * PIECE;
    let mut flat = vec![0.0; snps.len() * k];
    let mut dot = Array2::<f64>::zeros((n, k));
    let mut x = Array2::<f64>::zeros((n, panel).f());
    let mut loads = Array2::<f64>::zeros((panel, k));
    for chunk in live.chunks(panel) {
        let width = chunk.len();
        x.slice_mut(s![.., ..width])
            .axis_iter_mut(Axis(1))
            .into_par_iter()
            .zip(chunk)
            .for_each(|(mut col, &c)| {
                for (val, &s) in col.iter_mut().zip(samples) {
                    *val = row_value(geno.get(snps[c], s), &norms[c]);
                }
            });
        let x = x.slice(s![.., ..width]);
        loads
            .slice_mut(s![..width, ..])
            .axis_chunks_iter_mut(Axis(0), PIECE)
            .into_par_iter()
            .enumerate()
            .for_each(|(p, mut out)| {
                let cols = x.slice(s![.., p * PIECE..p * PIECE + out.nrows()]);
                general_mat_mul(1.0, &cols.t(), v, 0.0, &mut out);
            });
        let loads = loads.slice(s![..width, ..]);
        for (row, &c) in loads.outer_iter().zip(chunk) {
            flat[c * k..(c + 1) * k].copy_from_slice(row.as_slice().expect("row-major"));
        }
        dot.axis_chunks_iter_mut(Axis(0), PIECE)
            .into_par_iter()
            .enumerate()
            .for_each(|(p, mut out)| {
                let rows = x.slice(s![p * PIECE..p * PIECE + out.nrows(), ..]);
                general_mat_mul(1.0, &rows, &loads, 1.0, &mut out);
            });
    }
    let (missing, n_missing) = if lsq {
        let per_sample: Vec<(Vec<f64>, usize)> = samples
            .par_iter()
            .map(|&s| {
                let mut m = vec![0.0; k * k];
                let mut count = 0;
                for (col, &snp) in snps.iter().enumerate() {
                    if geno.get(snp, s) > 2 {
                        count += 1;
                        let l = &flat[col * k..(col + 1) * k];
                        for a in 0..k {
                            for b in a..k {
                                m[a * k + b] += l[a] * l[b];
                            }
                        }
                    }
                }
                (m, count)
            })
            .collect();
        let count = per_sample.iter().map(|p| p.1).collect();
        (per_sample.into_iter().flat_map(|p| p.0).collect(), count)
    } else {
        (Vec::new(), vec![0; n])
    };
    let sums = SampleSums {
        dot: dot.iter().copied().collect(),
        missing,
        n_missing,
    };
    (flat, sums)
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
mod tests;
