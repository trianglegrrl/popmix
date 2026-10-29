use super::*;
use crate::genotypes::{Chromosome, Sample};
use ndarray::Array2;

const BLOCK: f64 = 0.0625;

/// Genotypes of `pops.len()` populations of two samples over `n` SNPs,
/// 1/64 Morgan apart on chromosome 1 (so 4 SNPs per block of BLOCK).
fn matrix(n: usize, pops: &[&str], geno: impl Fn(usize, usize) -> u8) -> GenotypeMatrix {
    let samples: Vec<Sample> = pops
        .iter()
        .flat_map(|p| (0..2).map(move |k| Sample::new(format!("{p}{k}"), *p)))
        .collect();
    let snps: Vec<Snp> = (0..n)
        .map(|i| {
            Snp::new(
                format!("rs{i}"),
                Chromosome::Autosome(1),
                i as u32 * 1000,
                i as f64 / 64.0,
                'A',
                'G',
            )
        })
        .collect();
    let data = Array2::from_shape_fn((n, samples.len()), |(s, k)| geno(s, k));
    GenotypeMatrix::new(data, samples, snps).unwrap()
}

fn populations(m: &GenotypeMatrix, names: &[&str], ploidy: u8) -> GenoPopulations {
    GenoPopulations {
        names: names.iter().map(std::string::ToString::to_string).collect(),
        sample_pop: m
            .samples()
            .iter()
            .map(|s| {
                names
                    .iter()
                    .position(|n| *n == s.population.as_str())
                    .unwrap()
            })
            .collect(),
        ploidy: vec![ploidy; m.n_samples()],
    }
}

fn pseudo_random(s: usize, k: usize) -> u8 {
    ((s * 7 + k * 13 + (s * k) % 5) % 3) as u8
}

fn rows(m: &GenotypeMatrix) -> MatrixRows<'_> {
    MatrixRows::new(m, (0..m.n_snps()).collect(), (0..m.n_samples()).collect())
}

#[test]
fn statistics_do_not_depend_on_the_thread_count() {
    // more blocks than one batch, with missing calls
    let names = ["A", "B", "C", "D"];
    let m = matrix(3_000, &names, |s, k| {
        if (s + k) % 11 == 0 {
            3
        } else {
            pseudo_random(s, k)
        }
    });
    let pops = populations(&m, &names, 2);
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                let f4 = f4_from_genotypes(
                    &mut rows(&m),
                    m.snps(),
                    &pops,
                    &[[0, 1, 2, 3], [0, 2, 1, 3]],
                    F4Mode::D,
                    BLOCK,
                )
                .unwrap();
                let f3 = f3_from_genotypes(
                    &mut rows(&m),
                    m.snps(),
                    &pops,
                    &[[0, 1, 2]],
                    F3GenoOptions::default(),
                    BLOCK,
                )
                .unwrap();
                let est = |g: &GenoStats| g.stats.iter().map(|s| (s.est, s.se)).collect::<Vec<_>>();
                (est(&f4), est(&f3))
            })
    };
    let one = run(1);
    for threads in [3, 8] {
        assert_eq!(run(threads), one, "{threads} threads");
    }
}

#[test]
fn f4_without_missing_data_is_the_mean_of_the_snp_products() {
    let names = ["A", "B", "C", "D"];
    let m = matrix(40, &names, pseudo_random);
    let pops = populations(&m, &names, 2);
    let out = f4_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &pops,
        &[[0, 1, 2, 3]],
        F4Mode::F4,
        BLOCK,
    )
    .unwrap();
    assert_eq!(out.block_lengths, vec![4; 10]);
    let af =
        |s: usize, p: usize| f64::from(pseudo_random(s, 2 * p) + pseudo_random(s, 2 * p + 1)) / 4.0;
    let mean = (0..40)
        .map(|s| (af(s, 0) - af(s, 1)) * (af(s, 2) - af(s, 3)))
        .sum::<f64>()
        / 40.0;
    // Equal blocks: the jackknife estimate is the plain mean.
    assert!((out.stats[0].est - mean).abs() < 1e-15);
    assert_eq!(out.stats[0].n, 40);
    assert!(out.stats[0].se > 0.0);
}

#[test]
fn swapping_pop1_and_pop2_negates_f4_and_d() {
    let names = ["A", "B", "C", "D"];
    let m = matrix(40, &names, pseudo_random);
    let pops = populations(&m, &names, 2);
    for mode in [F4Mode::F4, F4Mode::D] {
        let out = f4_from_genotypes(
            &mut rows(&m),
            m.snps(),
            &pops,
            &[[0, 1, 2, 3], [1, 0, 2, 3]],
            mode,
            BLOCK,
        )
        .unwrap();
        assert!((out.stats[0].est + out.stats[1].est).abs() < 1e-15);
        assert!((out.stats[0].se - out.stats[1].se).abs() < 1e-15);
    }
}

#[test]
fn d_is_bounded_by_one() {
    // A, C fixed for one allele, B, D for the other: every SNP is BABA.
    let names = ["A", "B", "C", "D"];
    let m = matrix(20, &names, |_, k| if k / 2 % 2 == 0 { 0 } else { 2 });
    let pops = populations(&m, &names, 2);
    let out = f4_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &pops,
        &[[0, 1, 2, 3]],
        F4Mode::D,
        BLOCK,
    )
    .unwrap();
    assert!((out.stats[0].est - 1.0).abs() < 1e-15);
}

#[test]
fn z_and_p_are_nan_without_an_se() {
    // every SNP is BABA, so D is 1 in every block and the SE is 0
    let names = ["A", "B", "C", "D"];
    let m = matrix(20, &names, |_, k| if k / 2 % 2 == 0 { 0 } else { 2 });
    let pops = populations(&m, &names, 2);
    let out = f4_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &pops,
        &[[0, 1, 2, 3]],
        F4Mode::D,
        BLOCK,
    )
    .unwrap();
    assert_eq!(out.stats[0].se, 0.0);
    assert!(out.stats[0].z.is_nan());
    assert!(out.stats[0].p.is_nan());
}

#[test]
fn missing_genotypes_drop_only_the_statistics_that_use_them() {
    let names = ["A", "B", "C", "D", "E"];
    // E (samples 8, 9) has no data at SNP 3
    let m = matrix(20, &names, |s, k| {
        if s == 3 && k >= 8 {
            3
        } else {
            pseudo_random(s, k)
        }
    });
    let pops = populations(&m, &names, 2);
    let out = f4_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &pops,
        &[[0, 1, 2, 3], [0, 1, 2, 4]],
        F4Mode::F4,
        BLOCK,
    )
    .unwrap();
    assert_eq!(out.stats[0].n, 20);
    assert_eq!(out.stats[1].n, 19);
    assert_eq!(out.stats[1].blocks.n, vec![3, 4, 4, 4, 4]);
}

#[test]
fn outgroup_f3_is_the_mean_of_the_corrected_products() {
    let names = ["O", "A", "B"];
    let m = matrix(20, &names, pseudo_random);
    let pops = populations(&m, &names, 2);
    let opts = F3GenoOptions {
        apply_corr: true,
        outgroup_mode: true,
    };
    let out = f3_from_genotypes(&mut rows(&m), m.snps(), &pops, &[[0, 1, 2]], opts, BLOCK).unwrap();
    let g = |s: usize, k: usize| f64::from(pseudo_random(s, k));
    let mean = (0..20)
        .map(|s| {
            let af = |p: usize| (g(s, 2 * p) + g(s, 2 * p + 1)) / 4.0;
            let r = g(s, 0) + g(s, 1);
            let h = r * (4.0 - r) / 12.0;
            (af(0) - af(1)) * (af(0) - af(2)) - h / 4.0
        })
        .sum::<f64>()
        / 20.0;
    assert!((out.stats[0].est - mean).abs() < 1e-15);
}

#[test]
fn pseudo_haploid_targets_use_one_allele_per_sample_in_the_correction() {
    let names = ["T", "A", "B"];
    let m = matrix(20, &names, |s, k| 2 * (pseudo_random(s, k) / 2));
    let opts = F3GenoOptions::default();
    let diploid = f3_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &populations(&m, &names, 2),
        &[[0, 1, 2]],
        opts,
        BLOCK,
    )
    .unwrap();
    let haploid = f3_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &populations(&m, &names, 1),
        &[[0, 1, 2]],
        opts,
        BLOCK,
    )
    .unwrap();
    // Same frequencies, different sample sizes: the correction differs.
    assert!((diploid.stats[0].est - haploid.stats[0].est).abs() > 1e-6);
}

#[test]
fn f3_of_a_single_pseudo_haploid_target_is_not_available() {
    let names = ["T", "A", "B"];
    let m = matrix(20, &names, |s, k| 2 * (pseudo_random(s, k) / 2));
    let mut pops = populations(&m, &names, 2);
    // T is one pseudo-haploid sample: its correction h / n is 0/0
    pops.sample_pop[1] = 1;
    pops.ploidy[0] = 1;
    let out = f3_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &pops,
        &[[0, 1, 2]],
        F3GenoOptions::default(),
        BLOCK,
    )
    .unwrap();
    assert!(out.stats[0].est.is_nan());
    assert_eq!(out.stats[0].n, 0);
    let uncorrected = F3GenoOptions {
        apply_corr: false,
        outgroup_mode: true,
    };
    let out = f3_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &pops,
        &[[0, 1, 2]],
        uncorrected,
        BLOCK,
    )
    .unwrap();
    assert!(out.stats[0].est.is_finite());
    assert_eq!(out.stats[0].n, 20);
}

#[test]
fn out_of_range_populations_are_an_error() {
    let names = ["A", "B", "C", "D"];
    let m = matrix(10, &names, pseudo_random);
    let pops = populations(&m, &names, 2);
    let err = f4_from_genotypes(
        &mut rows(&m),
        m.snps(),
        &pops,
        &[[0, 1, 2, 4]],
        F4Mode::F4,
        BLOCK,
    );
    assert!(err.is_err());
}

#[test]
fn uncorrected_outgroup_f3_is_the_mean_of_the_products() {
    let names = ["O", "A", "B"];
    let m = matrix(20, &names, pseudo_random);
    let pops = populations(&m, &names, 2);
    let opts = F3GenoOptions {
        apply_corr: false,
        outgroup_mode: true,
    };
    let out = f3_from_genotypes(&mut rows(&m), m.snps(), &pops, &[[0, 1, 2]], opts, BLOCK).unwrap();
    let g = |s: usize, k: usize| f64::from(pseudo_random(s, k));
    let mean = (0..20)
        .map(|s| {
            let af = |p: usize| (g(s, 2 * p) + g(s, 2 * p + 1)) / 4.0;
            (af(0) - af(1)) * (af(0) - af(2))
        })
        .sum::<f64>()
        / 20.0;
    assert!((out.stats[0].est - mean).abs() < 1e-15);
}

#[test]
fn uncorrected_f3_is_scaled_by_twice_the_target_heterozygosity() {
    // without the outgroup mode, each block's denominator is the mean of
    // 2h over its SNPs, h the target's unbiased heterozygosity
    let names = ["T", "A", "B"];
    let m = matrix(20, &names, pseudo_random);
    let pops = populations(&m, &names, 2);
    let opts = F3GenoOptions {
        apply_corr: false,
        outgroup_mode: false,
    };
    let out = f3_from_genotypes(&mut rows(&m), m.snps(), &pops, &[[0, 1, 2]], opts, BLOCK).unwrap();
    let two_h = |s: usize| {
        let r = f64::from(pseudo_random(s, 0) + pseudo_random(s, 1));
        2.0 * r * (4.0 - r) / 12.0
    };
    let g = |s: usize, k: usize| f64::from(pseudo_random(s, k));
    let product = |s: usize| {
        let af = |p: usize| (g(s, 2 * p) + g(s, 2 * p + 1)) / 4.0;
        (af(0) - af(1)) * (af(0) - af(2))
    };
    // 4 SNPs per block; the numerator is not corrected
    let block_means = |f: &dyn Fn(usize) -> f64| -> Vec<f64> {
        (0..5)
            .map(|b| (4 * b..4 * b + 4).map(f).sum::<f64>() / 4.0)
            .collect()
    };
    let close = |got: &[f64], want: &[f64]| {
        got.len() == want.len() && got.iter().zip(want).all(|(g, w)| (g - w).abs() < 1e-15)
    };
    let blocks = &out.stats[0].blocks;
    assert!(
        close(&blocks.den, &block_means(&two_h)),
        "den {:?}",
        blocks.den
    );
    assert!(
        close(&blocks.num, &block_means(&product)),
        "num {:?}",
        blocks.num
    );
}
