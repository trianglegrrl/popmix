//! The date in years: the sample's age plus generations × generation time
//! (TASKS.md G3).
//!
//! DATES dates admixture in generations before the admixed samples lived.
//! In years before present that is t = a + g T for a sample age a (years
//! BP), g generations and a generation time of T years. With independent
//! errors of SD σ_a, σ_g and σ_T the variance of the product g T is exactly
//! g² σ_T² + T² σ_g² + σ_g² σ_T² (Goodman 1960, "On the exact variance of
//! products", JASA 55), so Var t = σ_a² + g² σ_T² + T² σ_g² + σ_g² σ_T².

/// Years before present and their SE from generations `g` ± `g_se`, a
/// generation time `t` ± `t_se` and a sample age `age` ± `age_se`. The SE is
/// NaN when `g_se` is (no jackknife).
pub(crate) fn years_bp(g: f64, g_se: f64, t: f64, t_se: f64, age: f64, age_se: f64) -> (f64, f64) {
    let var =
        age_se * age_se + g * g * t_se * t_se + t * t * g_se * g_se + g_se * g_se * t_se * t_se;
    (age + g * t, var.sqrt())
}

#[cfg(test)]
mod tests {
    use super::years_bp;

    #[test]
    fn without_age_or_generation_time_uncertainty_years_are_generations_times_t() {
        assert_eq!(years_bp(100.0, 10.0, 29.0, 0.0, 0.0, 0.0), (2900.0, 290.0));
    }

    #[test]
    fn the_sample_age_is_added() {
        assert_eq!(
            years_bp(100.0, 10.0, 29.0, 0.0, 4500.0, 0.0),
            (7400.0, 290.0)
        );
    }

    #[test]
    fn independent_errors_add_in_variance() {
        // Var = 30² + 100²·2² + 29²·10² + 10²·2² = 900 + 40000 + 84100 + 400
        let (y, se) = years_bp(100.0, 10.0, 29.0, 2.0, 4500.0, 30.0);
        assert_eq!(y, 7400.0);
        assert!((se - 125_400.0_f64.sqrt()).abs() < 1e-9, "{se}");
    }

    #[test]
    fn no_generations_se_gives_no_years_se() {
        assert!(years_bp(100.0, f64::NAN, 29.0, 2.0, 0.0, 0.0).1.is_nan());
    }
}
