use super::*;
use approx::assert_relative_eq;

#[test]
fn untested_models_rank_after_tested_ones() {
    let mut p = [0.2, f64::NAN, 0.9, 0.01, f64::NAN, 0.5];
    p.sort_by(|a, b| by_p_descending(*a, *b));
    assert_eq!(&p[..4], &[0.9, 0.5, 0.2, 0.01]);
    assert!(p[4..].iter().all(|x| x.is_nan()));
}

#[test]
fn a_saturated_model_has_no_p_value() {
    assert!(chisq_p_value(0.3, 0).is_nan());
    assert!(chisq_p_value(0.0, 0).is_nan());
    // R pchisq for nested differences keeps its point mass at 0
    assert_eq!(chisq_upper_tail(0.3, 0), 0.0);
}

#[test]
fn test_rank_one_matrix_fits_rank_one_exactly() {
    // x = u v' has rank 1, so the rank-1 residual is (near) zero.
    let x = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 2.0, 4.0, 6.0]);
    let qinv = DMatrix::<f64>::identity(6, 6);
    let fit = fit_rank(&x, &qinv, 1, 0.0).unwrap();
    let (chisq, dof) = evaluate_fit(&x, &qinv, &fit, 1);
    assert!(chisq < 1e-20, "chisq {chisq}");
    assert_eq!(dof, 2);
}

#[test]
fn test_rank_zero_chisq_is_the_full_quadratic_form() {
    let x = DMatrix::from_row_slice(1, 2, &[0.5, -1.0]);
    let qinv = DMatrix::from_row_slice(2, 2, &[2.0, 0.5, 0.5, 1.0]);
    let fit = fit_rank(&x, &qinv, 0, 1e-4).unwrap();
    let (chisq, dof) = evaluate_fit(&x, &qinv, &fit, 0);
    // [0.5, -1] Q [0.5, -1]' = 0.5 - 0.5 + 1
    assert_relative_eq!(chisq, 1.0, max_relative = 1e-14);
    assert_eq!(dof, 2);
}

#[test]
fn test_fudge_scales_with_trace() {
    let cov = ndarray::arr2(&[[2.0, 0.0], [0.0, 6.0]]);
    let q = fudged_inverse(&cov, 0.5).unwrap();
    // diag + 0.5 * 8 = [6, 10]
    assert_relative_eq!(q[(0, 0)], 1.0 / 6.0, max_relative = 1e-14);
    assert_relative_eq!(q[(1, 1)], 1.0 / 10.0, max_relative = 1e-14);
}

#[test]
fn test_simplex_qp_matches_unconstrained_inside_and_clamps_outside() {
    let d = DMatrix::from_row_slice(2, 2, &[2.0, 0.5, 0.5, 1.0]);
    // Unconstrained optimum with sum 1 inside the box: same as the KKT solution
    let w = simplex_qp(&d, &DVector::from_vec(vec![1.0, 1.0])).unwrap();
    assert_relative_eq!(w.sum(), 1.0, max_relative = 1e-14);
    assert!(w.iter().all(|&v| (0.0..=1.0).contains(&v)));
    // A strong pull towards the first weight pins it at 1
    let w = simplex_qp(&d, &DVector::from_vec(vec![10.0, -10.0])).unwrap();
    assert_eq!(w.as_slice(), &[1.0, 0.0]);
}

#[test]
fn test_simplex_qp_needs_a_square_d_of_its_size() {
    let d = DMatrix::<f64>::identity(3, 3);
    assert!(simplex_qp(&d, &DVector::from_vec(vec![1.0, 1.0])).is_err());
}

#[test]
fn test_simplex_qp_of_many_weights_converges() {
    // 60 weights, the first pulled to 1
    let d = DMatrix::<f64>::identity(60, 60);
    let mut v = DVector::zeros(60);
    v[0] = 10.0;
    let w = simplex_qp(&d, &v).unwrap();
    assert_eq!(w[0], 1.0);
}

/// Uniform numbers in [0, 1) from a fixed seed.
fn uniforms(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed
        .wrapping_mul(2862933555777941757)
        .wrapping_add(3037000493);
    move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn qp_objective(d: &DMatrix<f64>, v: &DVector<f64>, w: &DVector<f64>) -> f64 {
    0.5 * w.dot(&(d * w)) - v.dot(w)
}

/// The minimum of the simplex QP by brute force: on every support S, the
/// optimum with w_S summing to 1 and the rest 0 (its KKT system), kept if
/// w_S >= 0. Supports whose KKT system is singular are skipped.
fn simplex_qp_by_brute_force(d: &DMatrix<f64>, v: &DVector<f64>) -> DVector<f64> {
    let m = v.len();
    let mut best: Option<(f64, DVector<f64>)> = None;
    for mask in 1..(1usize << m) {
        let s: Vec<usize> = (0..m).filter(|&i| mask >> i & 1 == 1).collect();
        let k = s.len();
        let mut kkt = DMatrix::zeros(k + 1, k + 1);
        let mut rhs = DVector::zeros(k + 1);
        for (a, &i) in s.iter().enumerate() {
            for (b, &j) in s.iter().enumerate() {
                kkt[(a, b)] = d[(i, j)];
            }
            kkt[(a, k)] = -1.0;
            kkt[(k, a)] = 1.0;
            rhs[a] = v[i];
        }
        rhs[k] = 1.0;
        let Some(sol) = kkt.lu().solve(&rhs) else {
            continue;
        };
        if (0..k).any(|a| sol[a].is_nan() || sol[a] < -1e-12) {
            continue;
        }
        let mut w = DVector::zeros(m);
        for (a, &i) in s.iter().enumerate() {
            w[i] = sol[a];
        }
        let f = qp_objective(d, v, &w);
        if best.as_ref().is_none_or(|(g, _)| f < *g) {
            best = Some((f, w));
        }
    }
    best.unwrap().1
}

#[test]
fn test_simplex_qp_matches_brute_force_on_positive_definite_problems() {
    for seed in 0..400 {
        let mut u = uniforms(seed);
        let m = 2 + (seed as usize % 4);
        let b = DMatrix::from_fn(m, m, |_, _| 2.0 * u() - 1.0);
        let d = b.transpose() * &b + DMatrix::identity(m, m) * 0.01;
        let v = DVector::from_fn(m, |_, _| 4.0 * u() - 2.0);
        let got = simplex_qp(&d, &v).unwrap();
        let want = simplex_qp_by_brute_force(&d, &v);
        assert!((&got - &want).amax() < 1e-9, "seed {seed}: {got} != {want}");
    }
}

#[test]
fn test_simplex_qp_matches_brute_force_on_qpadm_problems() {
    // qpAdm's: D = CᵀC, e = (0, ..., 0, 1), C = [Aᵀ; 1ᵀ] with A m x (m - 1)
    for seed in 0..400 {
        let mut u = uniforms(seed + 1000);
        let m = 2 + (seed as usize % 4);
        let mut c = DMatrix::from_fn(m, m, |_, _| 2.0 * u() - 1.0);
        c.row_mut(m - 1).fill(1.0);
        let mut e = DVector::zeros(m);
        e[m - 1] = 1.0;
        let (d, v) = (c.transpose() * &c, c.transpose() * e);
        let got = simplex_qp(&d, &v).unwrap();
        let want = simplex_qp_by_brute_force(&d, &v);
        assert!((&got - &want).amax() < 1e-9, "seed {seed}: {got} != {want}");
    }
}

#[test]
fn test_unconstrained_weights_are_normalised_to_sum_one() {
    // A = I (2 sources, rank 2): CᵀC w = Cᵀe gives w = (1/3, 1/3)
    let fit = RankFit {
        a: DMatrix::identity(2, 2),
        b: DMatrix::zeros(2, 1),
    };
    let w = weights_from_fit(&fit, false).unwrap();
    assert!(w.iter().all(|x| (x - 0.5).abs() < 1e-15), "{w:?}");
}

#[test]
fn test_model_failures_are_singular_or_infeasible_fits() {
    assert!(is_model_failure(&singular("x")));
    assert!(is_model_failure(&PopmixError::Infeasible {
        context: "x".to_string()
    }));
    assert!(!is_model_failure(&PopmixError::PopulationNotFound(
        "x".to_string()
    )));
}

#[test]
fn test_weight_std_errors_need_two_replicates() {
    let x = |a: f64| ndarray::arr2(&[[a, 0.5], [0.2, a]]);
    let qinv = DMatrix::<f64>::identity(4, 4);
    let none = qpadm_weight_std_errors(&[], &qinv, 1, 0.0, false).unwrap();
    assert_eq!(none.len(), 2);
    assert!(none.iter().all(|s| s.is_nan()));
    let two = qpadm_weight_std_errors(&[x(1.0), x(1.5)], &qinv, 1, 0.0, false).unwrap();
    assert!(two.iter().all(|s| s.is_finite()), "{two:?}");
}

#[test]
fn test_combinations_follow_r_combn_order() {
    assert_eq!(combinations(3, 2), vec![vec![0, 1], vec![0, 2], vec![1, 2]]);
    assert_eq!(combinations(3, 3), vec![vec![0, 1, 2]]);
    assert_eq!(combinations(2, 1), vec![vec![0], vec![1]]);
}

#[test]
fn test_chisq_upper_tail_with_zero_df() {
    // R: pchisq(-3, 0, lower.tail = FALSE) = 1; pchisq(3, 0, ...) = 0
    assert_eq!(chisq_upper_tail(-3.0, 0), 1.0);
    assert_eq!(chisq_upper_tail(3.0, 0), 0.0);
    assert_eq!(chisq_upper_tail(0.0, 0), 0.0);
}

fn row(pattern: &str, f4rank: usize, dof: usize, chi_sq: f64, feasible: bool) -> PopDropRow {
    PopDropRow {
        pattern: pattern.to_string(),
        wt: 2 - f4rank,
        dof,
        chi_sq,
        p_value: chisq_p_value(chi_sq, dof),
        f4rank,
        weights: Vec::new(),
        feasible,
        best: false,
        dofdiff: None,
        chisqdiff: None,
        p_nested: None,
    }
}

/// Three sources, in admixtools2 power_set order. Full model 000; of the
/// one-drop models 010 fits best but is infeasible, so the path takes
/// 001; of 001's children (011, 101) 101 fits best. 110 fits better
/// still but is not nested in 001.
fn three_source_rows() -> Vec<PopDropRow> {
    vec![
        row("011", 0, 5, 100.0, true),
        row("101", 0, 5, 50.0, true),
        row("110", 0, 5, 30.0, true),
        row("001", 1, 4, 12.0, true),
        row("010", 1, 4, 10.0, false),
        row("100", 1, 4, 15.0, true),
        row("000", 2, 3, 8.0, false),
    ]
}

/// (pattern, dofdiff, chisqdiff, p_nested) of each model on the path.
type Nested<'a> = (&'a str, Option<usize>, Option<f64>, Option<f64>);

fn nested_columns(rows: &[PopDropRow]) -> Vec<Nested<'_>> {
    rows.iter()
        .filter(|r| r.best)
        .map(|r| (r.pattern.as_str(), r.dofdiff, r.chisqdiff, r.p_nested))
        .collect()
}

#[test]
fn test_corrected_popdrop_path_is_nested_and_greedy() {
    // ADMIXTOOLS2-0001 oracle: each step drops one more source, taking
    // the feasible child with the lowest chi-squared, and is compared
    // with the model it came from (1 degree of freedom).
    let mut rows = three_source_rows();
    mark_nested_path(&mut rows, 2, Admixtools2Compat::CORRECTED);
    let mut got = nested_columns(&rows);
    got.sort_by(|a, b| a.0.cmp(b.0));
    assert_eq!(
        got,
        vec![
            ("000", None, None, None),
            ("001", Some(1), Some(4.0), Some(chisq_upper_tail(4.0, 1))),
            ("101", Some(1), Some(38.0), Some(chisq_upper_tail(38.0, 1))),
        ]
    );
}

#[test]
fn test_corrected_popdrop_path_takes_the_first_of_tied_children() {
    let mut rows = three_source_rows();
    // 001's children 011 and 101 fit equally well
    rows[1].chi_sq = 100.0;
    mark_nested_path(&mut rows, 2, Admixtools2Compat::CORRECTED);
    assert!(rows[0].best);
    assert!(!rows[1].best);
}

#[test]
fn test_corrected_popdrop_path_stops_without_a_feasible_child() {
    let mut rows = three_source_rows();
    for r in rows.iter_mut().filter(|r| r.f4rank == 1) {
        r.feasible = false;
    }
    mark_nested_path(&mut rows, 2, Admixtools2Compat::CORRECTED);
    assert_eq!(nested_columns(&rows), vec![("000", None, None, None)]);
}

#[test]
fn test_v2_0_10_popdrop_compares_the_one_drop_siblings() {
    // ADMIXTOOLS2-0001 compat: admixtools2 2.0.10 marks every model at
    // full rank - 1 and compares each with the next (same dof, so
    // dofdiff 0 and pchisq(x, 0) is 0 for x >= 0, 1 for x < 0).
    let mut rows = three_source_rows();
    mark_nested_path(&mut rows, 2, Admixtools2Compat::V2_0_10);
    assert_eq!(
        nested_columns(&rows),
        vec![
            ("001", Some(0), Some(2.0), Some(0.0)),
            ("010", Some(0), Some(-5.0), Some(1.0)),
            ("100", None, None, None),
        ]
    );
}

#[test]
fn test_popdrop_compat_flag_only_changes_the_nested_columns() {
    // ADMIXTOOLS2-0001 isolation.
    let x = DMatrix::from_fn(3, 5, |i, j| {
        ((i * 5 + j) as f64 * 1.7).sin() + 0.3 * i as f64
    });
    let qinv = DMatrix::from_fn(15, 15, |i, j| if i == j { 2.0 } else { 0.1 });
    let fixed = drop_pops(&x, &qinv, 1e-4, false, Admixtools2Compat::CORRECTED).unwrap();
    let compat = drop_pops(&x, &qinv, 1e-4, false, Admixtools2Compat::V2_0_10).unwrap();
    assert_eq!(fixed.len(), 7);
    for (a, b) in fixed.iter().zip(&compat) {
        assert_eq!(a.pattern, b.pattern);
        assert_eq!(
            (a.wt, a.dof, a.f4rank, a.feasible),
            (b.wt, b.dof, b.f4rank, b.feasible)
        );
        assert_eq!(a.chi_sq.to_bits(), b.chi_sq.to_bits());
        assert_eq!(a.p_value.to_bits(), b.p_value.to_bits());
        assert_eq!(a.weights, b.weights);
    }
    assert_ne!(nested_columns(&fixed), nested_columns(&compat));
}

#[test]
fn test_constrained_popdrop_with_three_sources() {
    // ADMIXTOOLS2-0002: admixtools2 qpadm(constrained = TRUE) stops here.
    // Every subset model is fitted at full rank, where the constrained
    // weights are well defined, and the chi-squareds do not depend on
    // the weights at all.
    let x = DMatrix::from_fn(3, 5, |i, j| {
        ((i * 5 + j) as f64 * 1.7).sin() + 0.3 * i as f64
    });
    let qinv = DMatrix::from_fn(15, 15, |i, j| if i == j { 2.0 } else { 0.1 });
    let free = drop_pops(&x, &qinv, 1e-4, false, Admixtools2Compat::CORRECTED).unwrap();
    let con = drop_pops(&x, &qinv, 1e-4, true, Admixtools2Compat::CORRECTED).unwrap();
    for (f, c) in free.iter().zip(&con) {
        assert_eq!(f.chi_sq.to_bits(), c.chi_sq.to_bits());
        assert!(
            c.feasible,
            "pat {}: constrained weights {:?}",
            c.pattern, c.weights
        );
        let sum: f64 = c.weights.iter().flatten().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }
}

#[test]
fn test_p_value_tail_is_exact() {
    // pchisq(129.6103737533115, 4, lower.tail = FALSE) in R
    assert_relative_eq!(
        chisq_p_value(129.6103737533115, 4),
        4.717642937981588e-27,
        max_relative = 1e-10
    );
}
