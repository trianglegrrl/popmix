//! Shrinkage of projected samples' coordinates, and its correction
//! (TASKS.md G6).
//!
//! A PCA sample's coordinate on component k is the dot product of its
//! normalized genotypes with loadings that it helped make: with G the
//! scaled sample matrix (trace n − 1), λ_k and u_k its eigenpairs, the
//! coordinate is proportional to Σ_j G_ij u_jk / √λ_k. Its own term G_ii
//! holds its noise as well as its signal, and a sample projected from the
//! same population has only the signal part. So, to first order, projected
//! samples sit closer to the origin by the factor
//!
//!   ρ_k = 1 − Σ_i n_i u_ik² / λ_k,   n_i = G_ii − Σ_{s ∈ S} λ_s u_is²,
//!
//! where n_i is sample i's noise (its diagonal less the part the signal
//! components S explain). The correction divides projected samples'
//! coordinates by ρ_k. It is applied only on components that the
//! Tracy-Widom test finds significant (p < 0.05, from the top), the
//! components S: on a noise component ρ is near 0 and dividing by it would
//! only amplify noise.
//!
//! Provenance: derived here from the definition of the projection (the
//! problem is described in EIGENSOFT's documentation of `shrinkmode`,
//! whose method it does not describe and popmix does not reproduce), and
//! checked by the oracle EIGENSOFT's documentation states: two halves of one
//! population, one in the PCA and one projected, should coincide
//! (`tests/python/test_pca_shrinkage.py`).

use nalgebra::DMatrix;

use super::tracy_widom::TracyWidomRow;

/// Significance level of the Tracy-Widom test for a signal component.
const SIGNAL_P: f64 = 0.05;

/// The leading components whose Tracy-Widom p is below 0.05.
pub(crate) fn signal_components(tw: &[TracyWidomRow]) -> usize {
    tw.iter()
        .take_while(|row| row.p.is_some_and(|p| p < SIGNAL_P))
        .count()
}

/// ρ_k for each of the first `k` components: the factor by which projected
/// samples are shrunk, for the `signal` leading components, NaN beyond.
/// `diag` is the scaled matrix's diagonal and `evecs` its eigenvectors, one
/// row per PCA sample.
pub(crate) fn shrink_factors(
    diag: &[f64],
    lambda: &[f64],
    evecs: &DMatrix<f64>,
    k: usize,
    signal: usize,
) -> Vec<f64> {
    let noise: Vec<f64> = (0..diag.len())
        .map(|i| {
            diag[i]
                - (0..signal)
                    .map(|s| lambda[s] * evecs[(i, s)].powi(2))
                    .sum::<f64>()
        })
        .collect();
    (0..k)
        .map(|j| {
            if j >= signal {
                return f64::NAN;
            }
            let excess: f64 = (0..diag.len())
                .map(|i| noise[i] * evecs[(i, j)].powi(2))
                .sum();
            1.0 - excess / lambda[j]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_factor_is_one_less_the_weighted_noise_over_the_eigenvalue() {
        // two samples, one signal component with eigenvalue 1.5 on
        // (1, 1)/√2: each diagonal 1 holds 0.75 of signal, 0.25 of noise
        let diag = [1.0, 1.0];
        let lambda = [1.5, 0.5];
        let r = std::f64::consts::FRAC_1_SQRT_2;
        let evecs = DMatrix::from_row_slice(2, 2, &[r, r, r, -r]);
        let rho = shrink_factors(&diag, &lambda, &evecs, 2, 1);
        assert!((rho[0] - (1.0 - 0.25 / 1.5)).abs() < 1e-12, "{}", rho[0]);
        assert!(rho[1].is_nan());
    }

    #[test]
    fn signal_components_are_the_leading_significant_ones() {
        let row = |p: Option<f64>| TracyWidomRow {
            eigenvalue: 1.0,
            twstat: None,
            p,
            neff: None,
        };
        let tw = [
            row(Some(1e-9)),
            row(Some(0.01)),
            row(Some(0.3)),
            row(Some(0.001)),
        ];
        assert_eq!(signal_components(&tw), 2);
        assert_eq!(signal_components(&[row(None)]), 0);
    }
}
