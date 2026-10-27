//! Tracy-Widom statistics for PCA eigenvalues, as EIGENSOFT's smartpca and
//! twstats compute them (Patterson, Price and Reich 2006).
//!
//! For eigenvalues `lambda[i..m]` (decreasing, all positive), the effective
//! number of markers is estimated from their spread, the largest is centered
//! and scaled as Johnstone (2001), and the p-value is the right tail of the
//! Tracy-Widom distribution, interpolated from smartpca's table.

use super::tw_table::{TW_PDF, TW_TAIL, TW_X};

/// Which smartpca quirks to reproduce.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub struct EigensoftCompat {
    /// EIGENSOFT-0001: beyond the table (TW statistic above 8) smartpca's
    /// `gettw` returns the asymptotic density as the p-value (and the tail as
    /// the density). Corrected, the p-value is the asymptotic tail.
    pub tw_density_beyond_table: bool,
}

impl EigensoftCompat {
    /// Corrected behavior (the default).
    pub const CORRECTED: EigensoftCompat = EigensoftCompat {
        tw_density_beyond_table: false,
    };
    /// EIGENSOFT 7.2.1 smartpca exactly.
    pub const V7_2_1: EigensoftCompat = EigensoftCompat {
        tw_density_beyond_table: true,
    };
}

/// One row of smartpca's Tracy-Widom table.
#[derive(Debug, Clone, PartialEq)]
pub struct TracyWidomRow {
    pub eigenvalue: f64,
    /// TW statistic and p-value; None with fewer than `min_eigs` eigenvalues
    /// left.
    pub twstat: Option<f64>,
    pub p: Option<f64>,
    /// Effective number of markers.
    pub neff: Option<f64>,
}

/// Eigenvalues above this count towards the matrix rank (smartpca `numgtz`).
pub const RANK_THRESHOLD: f64 = 1e-6;

/// The Tracy-Widom table for decreasing eigenvalues `lambda`, one row per
/// positive eigenvalue. Rows with fewer than `min_eigs` eigenvalues from
/// there on have no statistic (smartpca `nostatslim`, at least 3).
pub fn tracy_widom_table(
    lambda: &[f64],
    min_eigs: usize,
    compat: EigensoftCompat,
) -> Vec<TracyWidomRow> {
    let m = lambda.iter().filter(|&&x| x > RANK_THRESHOLD).count();
    let min_eigs = min_eigs.max(3);
    (0..m)
        .map(|i| {
            let stat = tw_stat(&lambda[i..m], min_eigs, compat);
            TracyWidomRow {
                eigenvalue: lambda[i],
                twstat: stat.map(|s| s.0),
                p: stat.map(|s| s.1),
                neff: stat.map(|s| s.2),
            }
        })
        .collect()
}

/// (TW statistic, p-value, effective number of markers) for the largest of
/// `lambda` (smartpca `dotwcalc` with no given marker count).
fn tw_stat(lambda: &[f64], min_eigs: usize, compat: EigensoftCompat) -> Option<(f64, f64, f64)> {
    let m = lambda.len();
    if m < min_eigs {
        return None;
    }
    let lsum: f64 = lambda.iter().sum();
    if lsum <= 0.0 {
        return None;
    }
    let y = m as f64 / lsum;
    let evals: Vec<f64> = lambda.iter().map(|&l| l * y).collect();
    let top = (m * (m + 2)) as f64;
    let bot = evals.iter().map(|e| e * e).sum::<f64>() - m as f64;
    let zn = top / bot;
    let tw = tw_norm(evals[0] * zn, m as f64, zn);
    Some((tw, tw_tail(tw, compat), zn))
}

/// Center and scale an eigenvalue of a p x n Wishart matrix (Johnstone 2001).
fn tw_norm(lam: f64, p: f64, n: f64) -> f64 {
    if n < 0.0 || p < 0.0 {
        return -10.0;
    }
    if n < p {
        return tw_norm(lam, n, p);
    }
    let y1 = (n - 1.0).sqrt() + p.sqrt();
    let mu = y1 * y1;
    let y2 = 1.0 / (n - 1.0).sqrt() + 1.0 / p.sqrt();
    let phi = y1 * y2.powf(1.0 / 3.0);
    (lam - mu) / phi
}

const SQRT_PI: f64 = 1.0 / 0.564_189_583_547_756_3;

/// Asymptotic density and right tail beyond the table (Margetis-Edelman).
fn tw_density_asymptotic(tw: f64) -> f64 {
    if tw <= 0.0 {
        return 0.0;
    }
    (-0.25 * tw.ln() - 2.0 * tw.powf(1.5) / 3.0 - (SQRT_PI * 4.0).ln()).exp()
}

fn tw_tail_asymptotic(tw: f64) -> f64 {
    if tw <= 0.0 {
        return 1.0;
    }
    (-0.75 * tw.ln() - 2.0 * tw.powf(1.5) / 3.0 - (SQRT_PI * 4.0).ln()).exp()
}

/// smartpca's `firstgtx`: the first table index above `val`, by its own
/// recursive search.
fn first_above(tab: &[f64], lo: usize, hi: usize, val: f64) -> usize {
    if val >= tab[hi] {
        return hi + 1;
    }
    if val < tab[lo] {
        return lo;
    }
    let k = (lo + hi) / 2;
    if val <= tab[k] {
        first_above(tab, lo + 1, k, val)
    } else {
        first_above(tab, k, hi - 1, val)
    }
}

/// Right tail of the Tracy-Widom distribution: cubic Hermite interpolation
/// in smartpca's table (`gettw`), with the asymptotic form beyond it.
pub fn tw_tail(x: f64, compat: EigensoftCompat) -> f64 {
    let n = TW_X.len();
    let k = first_above(&TW_X, 0, n - 1, x);
    if k == 0 {
        return 1.0;
    }
    if k >= n {
        return if compat.tw_density_beyond_table {
            tw_density_asymptotic(x)
        } else {
            tw_tail_asymptotic(x)
        };
    }
    let (x0, x1) = (TW_X[k - 1], TW_X[k]);
    // tail and its derivative (minus the density) at both ends
    let (f0, f0p, f1, f1p) = (TW_TAIL[k - 1], -TW_PDF[k - 1], TW_TAIL[k], -TW_PDF[k]);
    let inc = x1 - x0;
    let t = (x - x0) / inc;
    let (a0, a1) = (f0, f0p * inc);
    let cc0 = f1 - (a0 + a1);
    let cc1 = f1p * inc - a1;
    let a2 = 3.0 * cc0 - cc1;
    let a3 = cc1 - 2.0 * cc0;
    ((a3 * t + a2) * t + a1) * t + a0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_tail_interpolates_the_table_exactly_at_its_points() {
        for (i, &x) in TW_X.iter().enumerate().take(TW_X.len() - 1) {
            let got = tw_tail(x, EigensoftCompat::CORRECTED);
            // the search puts a table point at the start of the next interval
            assert!(
                (got - TW_TAIL[i]).abs() < 1e-12,
                "x = {x}: {got} != {}",
                TW_TAIL[i]
            );
        }
    }

    #[test]
    fn the_tail_decreases() {
        let mut last = 1.0;
        let mut x = -8.5;
        while x < 12.0 {
            let t = tw_tail(x, EigensoftCompat::CORRECTED);
            // cubic interpolation overshoots by up to ~1e-11 where the table is flat
            assert!(t <= last + 1e-10, "tail rises at {x}: {t} > {last}");
            last = t;
            x += 0.013;
        }
    }

    #[test]
    fn beyond_the_table_smartpca_reports_the_density() {
        let x = 26.255226058550154;
        let corrected = tw_tail(x, EigensoftCompat::CORRECTED);
        let v721 = tw_tail(x, EigensoftCompat::V7_2_1);
        // smartpca 7.2.1 printed 6.9780333866491966e-41 here (family_all)
        assert!((v721 - 6.978_033_386_649_197e-41).abs() < 1e-52);
        // the tail is the density over sqrt(x)
        assert!((corrected / v721 - x.powf(-0.5)).abs() < 1e-12);
    }

    #[test]
    fn too_few_eigenvalues_have_no_statistic() {
        let lambda = [5.0, 2.0, 1.0, 1.0, 0.5];
        let rows = tracy_widom_table(&lambda, 10, EigensoftCompat::CORRECTED);
        assert_eq!(rows.len(), 5);
        assert!(rows.iter().all(|r| r.twstat.is_none()));
        let rows = tracy_widom_table(&lambda, 3, EigensoftCompat::CORRECTED);
        assert!(rows[2].twstat.is_some() && rows[3].twstat.is_none());
    }
}
