use super::*;
use crate::fstats::MatrixRows;
use crate::genotype_matrix::GenotypeMatrix;
use crate::genotypes::{Chromosome, Sample, Snp};
use ndarray::Array2;

/// Two source populations with different frequencies, and samples that
/// are mixtures of them with known proportions.
fn simulated(qs: &[f64], m: usize) -> (PackedGenotypes, Vec<f64>) {
    let mut rng = Rng(7);
    let p: Vec<[f64; 2]> = (0..m)
        .map(|_| [0.05 + 0.9 * rng.uniform(), 0.05 + 0.9 * rng.uniform()])
        .collect();
    let n = qs.len();
    let data = Array2::from_shape_fn((m, n), |(j, i)| {
        let f = qs[i] * p[j][0] + (1.0 - qs[i]) * p[j][1];
        // count of the non-P allele
        (0..2).filter(|_| rng.uniform() >= f).count() as u8
    });
    let samples: Vec<Sample> = (0..n).map(|i| Sample::new(format!("s{i}"), "P")).collect();
    let snps: Vec<Snp> = (0..m)
        .map(|j| {
            Snp::new(
                format!("rs{j}"),
                Chromosome::Autosome(1),
                j as u32,
                0.0,
                'A',
                'G',
            )
        })
        .collect();
    let mat = GenotypeMatrix::new(data, samples, snps).unwrap();
    let mut rows = MatrixRows::new(&mat, (0..m).collect(), (0..n).collect());
    let flat: Vec<f64> = p.iter().flat_map(|x| x.iter().copied()).collect();
    (PackedGenotypes::read(&mut rows, n, m).unwrap(), flat)
}

#[test]
fn unsupervised_recovers_the_proportions() {
    let qs = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.25, 0.75, 1.0, 0.0, 0.5];
    let (g, _) = simulated(&qs, 4000);
    let opts = AdmixtureOptions {
        tolerance: 1e-8,
        ..Default::default()
    };
    let fit = admixture(&g, 2, &opts).unwrap();
    assert!(fit.converged);
    // cluster labels are arbitrary
    let col = usize::from(fit.q[0] <= 0.5);
    for (i, &want) in qs.iter().enumerate() {
        assert!(
            (fit.q[i * 2 + col] - want).abs() < 0.05,
            "sample {i}: {} vs {want}",
            fit.q[i * 2 + col]
        );
    }
}

#[test]
fn the_loglikelihood_never_decreases_and_matches_its_definition() {
    let qs = [1.0, 0.0, 0.5, 0.3, 0.9, 0.1];
    let (g, _) = simulated(&qs, 500);
    let fit = admixture(&g, 2, &AdmixtureOptions::default()).unwrap();
    let again = admixture_loglik(&g, &fit.q, &fit.p, 2).unwrap();
    assert!((again - fit.loglik).abs() < 1e-9 * fit.loglik.abs());
    let (q0, p0) = random_start(6, 500, 2, 1);
    assert!(admixture_loglik(&g, &q0, &p0, 2).unwrap() < fit.loglik);
}

#[test]
fn the_terms_of_some_samples_are_their_terms_among_all() {
    let qs = [1.0, 0.0, 0.5, 0.3, 0.9, 0.1, 0.2, 0.6];
    let (g, _) = simulated(&qs, 3_000);
    let (q, p) = random_start(8, 3_000, 2, 3);
    let d = Data {
        geno: &g,
        held_out: None,
    };
    let (all, _, _) = q_terms(&d, &q, &p, 2, false);
    let some = sample_ll(&d, &q, &p, 2, &[6, 1, 3]);
    assert_eq!(some, vec![all[6], all[1], all[3]]);
}

#[test]
fn projection_with_the_true_frequencies_recovers_the_proportions() {
    let qs = [0.2, 0.8, 0.5, 1.0];
    let (g, p) = simulated(&qs, 20000);
    let opts = AdmixtureOptions {
        tolerance: 1e-10,
        ..Default::default()
    };
    let fit = admixture_projection(&g, &p, 2, &opts).unwrap();
    for (i, &want) in qs.iter().enumerate() {
        assert!(
            (fit.q[i * 2] - want).abs() < 0.02,
            "sample {i}: {} vs {want}",
            fit.q[i * 2]
        );
    }
    assert_eq!(fit.p.len(), p.len());
}

#[test]
fn supervised_samples_keep_their_population() {
    let qs = [1.0, 1.0, 0.0, 0.0, 0.4];
    let (g, _) = simulated(&qs, 2000);
    let labels = [Some(0), Some(0), Some(1), Some(1), None];
    let opts = AdmixtureOptions {
        tolerance: 1e-8,
        ..Default::default()
    };
    let fit = admixture_supervised(&g, 2, &labels, &opts).unwrap();
    // sample 0 is in population 0, sample 2 in population 1
    assert!(fit.q[0] > 0.9999 && fit.q[4] < 1e-4 && fit.q[5] > 0.9999);
    assert!((fit.q[8] - 0.4).abs() < 0.05, "{}", fit.q[8]);
}

#[test]
fn cross_validation_prefers_the_true_k() {
    let qs: Vec<f64> = (0..20).map(|i| if i < 10 { 1.0 } else { 0.0 }).collect();
    let (g, _) = simulated(&qs, 1500);
    let opts = AdmixtureOptions {
        tolerance: 1e-6,
        ..Default::default()
    };
    let cv = |k| admixture_cv(&g, &admixture(&g, k, &opts).unwrap(), 5, &opts).unwrap();
    assert!(cv(2) < cv(1));
}

#[test]
fn cross_validation_folds_are_checked() {
    let qs: Vec<f64> = (0..4).map(|i| if i < 2 { 1.0 } else { 0.0 }).collect();
    let (g, _) = simulated(&qs, 50);
    let opts = AdmixtureOptions::default();
    let fit = admixture(&g, 1, &opts).unwrap();
    for folds in [1, 257, 1000] {
        let err = admixture_cv(&g, &fit, folds, &opts)
            .unwrap_err()
            .to_string();
        assert!(err.contains("folds"), "{folds}: {err}");
    }
}

#[test]
fn one_logarithm_terms_are_the_definition() {
    for (a, b) in [(0.3, 0.7), (1e-5, 1.0 - 1e-5), (0.999, 0.001), (0.5, 0.25)] {
        for g in 0..3u8 {
            let h = 2.0 - f64::from(g);
            let want = h * f64::ln(a) + (2.0 - h) * f64::ln(b);
            // the product rounds by one ulp, which moves its logarithm
            // by about epsilon however small the logarithm is
            let tolerance = 4.0 * f64::EPSILON * (1.0 + want.abs());
            assert!((term(g, a, b) - want).abs() <= tolerance, "{a} {b} {g}");
        }
    }
}

#[test]
fn box_qp_solves_inside_and_on_the_bounds() {
    let a = [2.0, 0.5, 0.5, 1.0];
    // unconstrained minimum A^-1 g
    let x = box_qp(&a, &[1.0, 1.0], &[-10.0, -10.0], &[10.0, 10.0]);
    assert!(
        (2.0 * x[0] + 0.5 * x[1] - 1.0).abs() < 1e-12 && (0.5 * x[0] + x[1] - 1.0).abs() < 1e-12
    );
    let x = box_qp(&a, &[1.0, 1.0], &[-10.0, -10.0], &[10.0, 0.1]);
    assert_eq!(x[1], 0.1);
    assert!((x[0] - (1.0 - 0.05) / 2.0).abs() < 1e-12);
}

#[test]
fn simplex_qp_keeps_the_sum_and_the_bounds() {
    let a = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    // pull hard towards the first coordinate; the third may only drop by 0.1
    let x = simplex_qp(&a, &[3.0, 0.0, -3.0], &[-0.5, -0.5, -0.1]);
    assert!(x.iter().sum::<f64>().abs() < 1e-12);
    // x3 = -0.1 and then x2 = -0.5 hit their bounds; x1 takes the rest
    for (got, want) in x.iter().zip([0.6, -0.5, -0.1]) {
        assert!((got - want).abs() < 1e-12, "{x:?}");
    }
    let free = simplex_qp(&a, &[1.0, 0.0, -1.0], &[-5.0, -5.0, -5.0]);
    assert!((free[0] - 1.0).abs() < 1e-12 && free[1].abs() < 1e-12);
}

#[test]
fn bad_inputs_are_errors() {
    let (g, p) = simulated(&[0.5, 0.5], 10);
    assert!(admixture(&g, 0, &AdmixtureOptions::default()).is_err());
    assert!(admixture_projection(&g, &p[..4], 2, &AdmixtureOptions::default()).is_err());
    assert!(admixture_supervised(&g, 2, &[Some(2), None], &AdmixtureOptions::default()).is_err());
}
