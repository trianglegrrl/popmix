
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
fn the_cross_product_is_the_sum_over_snps() {
    // several tiles of samples and several panels of SNPs
    let (n, n_snps) = (2 * 150, 5_000);
    let g = two_groups(150, n_snps);
    let used: Vec<usize> = (0..n).filter(|s| s % 7 != 3).collect();
    let snps: Vec<usize> = (0..n_snps).collect();
    let norms: Vec<ColumnNorm> = snps.iter().map(|&j| column_norm(&g, j, &used)).collect();
    let got = cross_product(&g, &snps, &norms, &used);
    let x = |j: usize, s: usize| match g.get(j, s) {
        v if v > 2 || !norms[j].alive => 0.0,
        v => norms[j].values[v as usize],
    };
    for (a, &sa) in used.iter().enumerate().step_by(13) {
        for (b, &sb) in used.iter().enumerate() {
            let want: f64 = snps.iter().map(|&j| x(j, sa) * x(j, sb)).sum();
            assert!(
                (got[[a, b]] - want).abs() <= 1e-9 * want.abs().max(1.0),
                "{a},{b}"
            );
            assert_eq!(got[[a, b]], got[[b, a]]);
        }
    }
}

#[test]
fn eigen_returns_decreasing_unit_eigenvectors() {
    let n = 40;
    let m = Array2::from_shape_fn((n, n), |(i, j)| {
        1.0 / (1.0 + i.abs_diff(j) as f64) + if i == j { (i % 5) as f64 } else { 0.0 }
    });
    let (values, vectors) = eigen(&m);
    assert!(values.windows(2).all(|w| w[0] >= w[1]));
    for (k, &l) in values.iter().enumerate() {
        let v: Vec<f64> = (0..n).map(|i| vectors[(i, k)]).collect();
        assert!((v.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-12);
        for i in 0..n {
            let mv: f64 = (0..n).map(|j| m[[i, j]] * v[j]).sum();
            assert!((mv - l * v[i]).abs() < 1e-10, "eigenpair {k}");
        }
    }
}

#[test]
fn many_samples_give_the_same_results_at_every_thread_count() {
    let g = two_groups(150, 5_000);
    let roles = vec![SampleRole::Pca; 300];
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| smartpca(&g, &roles, &PcaOptions::default()).unwrap())
    };
    let one = run(1);
    for threads in [2, 5, 16] {
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
