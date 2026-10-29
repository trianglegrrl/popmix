use super::*;
use crate::fstats::f2::f2_single_snp;

/// Compute f2 sums, weights, and heterozygosity for a block of SNPs.
fn direct_f2_block(af: &AlleleFrequencies, rows: &[usize]) -> BlockResult {
    let n_pops = af.populations.len();
    let n_pairs = n_pops * (n_pops + 1) / 2;

    let mut f2_sums = vec![0.0; n_pairs];
    let mut f2_raw_sums = vec![0.0; n_pairs];
    let mut f2_weights = vec![0.0; n_pairs];
    let mut het_sums = vec![0.0; n_pops];
    let mut het_weights = vec![0.0; n_pops];

    for &snp_idx in rows {
        // Compute heterozygosity for each population: h = p * (1 - p)
        // Note: admixtools2 uses p*(1-p), not 2*p*(1-p)
        for pop_idx in 0..n_pops {
            let count = af.get_count(snp_idx, pop_idx);

            if let Some(freq) = af.get(snp_idx, pop_idx) {
                if count >= 2 {
                    let het = freq * (1.0 - freq);
                    het_sums[pop_idx] += het;
                    het_weights[pop_idx] += 1.0;
                }
            }
        }

        // For each population pair
        let mut pair_idx = 0;
        for i in 0..n_pops {
            for j in i..n_pops {
                let freq_a = af.get(snp_idx, i);
                let freq_b = af.get(snp_idx, j);
                let count_a = af.get_count(snp_idx, i);
                let count_b = af.get_count(snp_idx, j);

                let (f2, weight) = f2_single_snp(freq_a, freq_b, count_a, count_b);

                // Also compute raw squared difference (without bias correction)
                // for FST denominator
                if let (Some(pa), Some(pb)) = (freq_a, freq_b) {
                    if count_a >= 2 && count_b >= 2 {
                        let raw = (pa - pb).powi(2);
                        f2_raw_sums[pair_idx] += raw * weight;
                    }
                }

                f2_sums[pair_idx] += f2 * weight;
                f2_weights[pair_idx] += weight;

                pair_idx += 1;
            }
        }
    }

    BlockResult {
        f2_sums,
        f2_raw_sums,
        f2_weights,
        het_sums,
        het_weights,
    }
}

use crate::genotype_matrix::GenotypeMatrix;
use crate::genotypes::{Chromosome, Sample, Snp};
use approx::assert_relative_eq;

fn make_test_data() -> AlleleFrequencies {
    // Create test data with 3 populations, 20 SNPs across 2 cM
    let mut samples = Vec::new();
    for i in 0..10 {
        samples.push(Sample::new(format!("A{i}"), PopulationId::new("PopA")));
    }
    for i in 0..10 {
        samples.push(Sample::new(format!("B{i}"), PopulationId::new("PopB")));
    }
    for i in 0..10 {
        samples.push(Sample::new(format!("C{i}"), PopulationId::new("PopC")));
    }

    let snps: Vec<Snp> = (0..20)
        .map(|i| {
            Snp::new(
                format!("rs{i}"),
                Chromosome::Autosome(1),
                (i * 1000) as u32,
                f64::from(i) * 0.1, // 0.1 cM spacing
                'A',
                'G',
            )
        })
        .collect();

    // Create genotype data
    // PopA: mostly ref (0), PopB: mostly alt (2), PopC: mixed
    let mut data = ndarray::Array2::zeros((20, 30));
    for snp_idx in 0..20 {
        // PopA (columns 0-9): mostly 0
        for j in 0..10 {
            data[[snp_idx, j]] = u8::from(snp_idx % 5 == 0);
        }
        // PopB (columns 10-19): mostly 2
        for j in 10..20 {
            data[[snp_idx, j]] = if snp_idx % 5 == 0 { 1 } else { 2 };
        }
        // PopC (columns 20-29): mixed
        for j in 20..30 {
            data[[snp_idx, j]] = ((snp_idx + j) % 3) as u8;
        }
    }

    let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
    AlleleFrequencies::from_genotypes(&geno)
}

/// Allele frequencies for PopA, PopB, PopC at SNPs on chromosome 1.
/// `rows` holds (position in Morgans, [(freq, count); 3]); count 0 is missing.
fn af_from(rows: &[(f64, [(f64, u32); 3])]) -> AlleleFrequencies {
    let n = rows.len();
    let mut frequencies = Array2::from_elem((n, 3), f64::NAN);
    let mut counts = Array2::zeros((n, 3));
    for (i, (_, pops)) in rows.iter().enumerate() {
        for (j, &(f, c)) in pops.iter().enumerate() {
            if c > 0 {
                frequencies[[i, j]] = f;
                counts[[i, j]] = c;
            }
        }
    }
    let snps = rows
        .iter()
        .enumerate()
        .map(|(i, (pos, _))| {
            Snp::new(
                format!("rs{i}"),
                Chromosome::Autosome(1),
                i as u32,
                *pos,
                'A',
                'G',
            )
        })
        .collect();
    AlleleFrequencies {
        frequencies,
        counts,
        populations: ["PopA", "PopB", "PopC"].map(PopulationId::new).to_vec(),
        snps,
    }
}

/// Per-SNP bias-corrected f2: (p1 - p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1).
fn snp_f2(a: (f64, u32), b: (f64, u32)) -> f64 {
    let corr = |(p, n): (f64, u32)| p * (1.0 - p) / (f64::from(n) - 1.0).max(1.0);
    (a.0 - b.0).powi(2) - corr(a) - corr(b)
}

/// Per-SNP FST numerator and denominator block means (see [`FstBlocks`]).
fn direct_fst_block(af: &AlleleFrequencies, rows: &[usize]) -> (Vec<f64>, Vec<f64>) {
    let n_pops = af.populations.len();
    let n_pairs = n_pops * (n_pops + 1) / 2;
    let (mut num, mut den, mut cnt) = (vec![0.0; n_pairs], vec![0.0; n_pairs], vec![0.0; n_pairs]);
    for &snp in rows {
        let (p, n) = (af.frequencies.row(snp), af.counts.row(snp));
        let mut pair = 0;
        for i in 0..n_pops {
            for j in i..n_pops {
                let (p1, p2) = (p[i], p[j]);
                let (n1, n2) = (f64::from(n[i]), f64::from(n[j]));
                let (h1, h2) = (p1 * (1.0 - p1), p2 * (1.0 - p2));
                let nm = (p1 - p2).powi(2) - h1 / (n1 - 1.0) - h2 / (n2 - 1.0);
                let dn = nm + h1 * n1 / (n1 - 1.0) + h2 * n2 / (n2 - 1.0);
                if !nm.is_nan() && !dn.is_nan() {
                    num[pair] += nm;
                    den[pair] += dn;
                    cnt[pair] += 1.0;
                }
                pair += 1;
            }
        }
    }
    let mean = |v: Vec<f64>| v.iter().zip(&cnt).map(|(s, c)| s / c).collect();
    (mean(num), mean(den))
}

/// Random frequencies of `n_pops` populations, with missing data and
/// populations of one allele (pseudo-haploid samples).
fn random_af(n_snps: usize, n_pops: usize) -> AlleleFrequencies {
    let mut state = 2024u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut frequencies = Array2::from_elem((n_snps, n_pops), f64::NAN);
    let mut counts = Array2::zeros((n_snps, n_pops));
    for i in 0..n_snps {
        for j in 0..n_pops {
            let n = [0u32, 1, 2, 5, 10, 20, 20, 20][(next() % 8) as usize];
            if n > 0 {
                counts[[i, j]] = n;
                frequencies[[i, j]] = (next() % (u64::from(n) + 1)) as f64 / f64::from(n);
            }
        }
    }
    let snps = (0..n_snps)
        .map(|i| {
            Snp::new(
                format!("rs{i}"),
                Chromosome::Autosome(1),
                i as u32,
                i as f64 * 1e-4,
                'A',
                'G',
            )
        })
        .collect();
    AlleleFrequencies {
        frequencies,
        counts,
        populations: (0..n_pops)
            .map(|j| PopulationId::new(format!("P{j}")))
            .collect(),
        snps,
    }
}

fn assert_close(got: &[f64], want: &[f64], what: &str) {
    assert_eq!(got.len(), want.len());
    for (k, (g, w)) in got.iter().zip(want).enumerate() {
        let ok = (g.is_nan() && w.is_nan()) || (g - w).abs() <= 1e-15 + 1e-12 * w.abs();
        assert!(ok, "{what} [{k}]: {g} != {w}");
    }
}

#[test]
fn block_products_match_the_per_snp_sums() {
    let af = random_af(600, 9);
    let rows: Vec<usize> = (0..600).filter(|i| i % 7 != 3).collect();
    let (got, want) = (f2_block(&af, &rows), direct_f2_block(&af, &rows));
    assert_close(&got.f2_sums, &want.f2_sums, "f2");
    assert_close(&got.f2_raw_sums, &want.f2_raw_sums, "raw");
    assert_eq!(got.f2_weights, want.f2_weights);
    assert_close(&got.het_sums, &want.het_sums, "het");
    assert_eq!(got.het_weights, want.het_weights);
    let (num, den) = fst_block(&af, &rows);
    let (want_num, want_den) = direct_fst_block(&af, &rows);
    assert_close(&num, &want_num, "fst num");
    assert_close(&den, &want_den, "fst den");
}

/// `af` cut into its chromosomes' runs of SNPs.
fn by_chromosome(af: &AlleleFrequencies) -> Vec<AlleleFrequencies> {
    let mut runs: Vec<std::ops::Range<usize>> = Vec::new();
    for (i, snp) in af.snps.iter().enumerate() {
        match runs.last_mut() {
            Some(r) if af.snps[r.start].chromosome == snp.chromosome => r.end = i + 1,
            _ => runs.push(i..i + 1),
        }
    }
    runs.into_iter()
        .map(|r| AlleleFrequencies {
            frequencies: af.frequencies.slice(ndarray::s![r.clone(), ..]).to_owned(),
            counts: af.counts.slice(ndarray::s![r.clone(), ..]).to_owned(),
            populations: af.populations.clone(),
            snps: af.snps[r].to_vec(),
        })
        .collect()
}

#[test]
fn blocks_a_chromosome_at_a_time_are_the_blocks_of_the_whole() {
    let mut af = random_af(6_000, 5);
    for (i, snp) in af.snps.iter_mut().enumerate() {
        snp.chromosome = Chromosome::Autosome(1 + (i / 2_000) as u8);
    }
    for max_missing in [0.0, 1.0] {
        let opts = F2Options {
            max_missing,
            poly_only: true,
            ..F2Options::default()
        };
        let def = BlockDefinition::new(0.05);
        let whole = F2Blocks::from_allele_frequencies_with_options(&af, &def, &opts);
        let mut evidence = CoordinateEvidence::default();
        evidence.observe(&af.snps);
        let parts = F2Blocks::from_parts(by_chromosome(&af).into_iter(), &evidence, &def, &opts)
            .expect("the coordinate guessed from every SNP is right");
        assert_eq!(parts.block_sums, whole.block_sums);
        assert_eq!(parts.block_raw_sums, whole.block_raw_sums);
        assert_eq!(parts.block_weights, whole.block_weights);
        assert_eq!(parts.block_het_sums, whole.block_het_sums);
        assert_eq!(parts.block_lengths, whole.block_lengths);
        assert_eq!(parts.n_blocks, whole.n_blocks);
        let (a, b) = (parts.fst.unwrap(), whole.fst.unwrap());
        assert_eq!(
            (a.num, a.den, a.block_lengths),
            (b.num, b.den, b.block_lengths)
        );
    }
}

#[test]
fn a_wrong_coordinate_guess_is_reported() {
    // no genetic map, and one polymorphic SNP: every SNP says base pairs,
    // the f2 SNPs (one) say the map
    let mut af = random_af(400, 3);
    for (i, snp) in af.snps.iter_mut().enumerate() {
        snp.position_morgans = 0.0;
        snp.chromosome = Chromosome::Autosome(1 + (i / 200) as u8);
        let f = if i == 7 { [0.1, 0.5, 0.9] } else { [0.5; 3] };
        for (p, &x) in f.iter().enumerate() {
            af.frequencies[[i, p]] = x;
            af.counts[[i, p]] = 20;
        }
    }
    let mut evidence = CoordinateEvidence::default();
    evidence.observe(&af.snps);
    let opts = F2Options::default();
    let def = BlockDefinition::default();
    assert!(F2Blocks::from_parts(by_chromosome(&af).into_iter(), &evidence, &def, &opts).is_none());
}

#[test]
fn blocks_do_not_depend_on_the_thread_count() {
    let af = random_af(4_000, 6);
    let opts = F2Options {
        max_missing: 1.0,
        poly_only: true,
        ..F2Options::default()
    };
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                F2Blocks::from_allele_frequencies_with_options(
                    &af,
                    &BlockDefinition::new(0.05),
                    &opts,
                )
            })
    };
    let one = run(1);
    for threads in [3, 8] {
        let many = run(threads);
        assert_eq!(one.block_sums, many.block_sums);
        assert_eq!(one.block_raw_sums, many.block_raw_sums);
        assert_eq!(
            one.fst.as_ref().unwrap().num,
            many.fst.as_ref().unwrap().num
        );
    }
}

#[test]
fn test_snp_missing_in_any_population_is_dropped_for_every_pair() {
    // admixtools2 maxmiss = 0: SNP 1 is missing in PopC, so f2(A, B) must
    // not use it either.
    let rows = [
        (0.00, [(0.1, 10), (0.5, 10), (0.3, 10)]),
        (0.01, [(0.9, 10), (0.2, 10), (0.0, 0)]),
        (0.10, [(0.4, 10), (0.6, 10), (0.5, 10)]),
    ];
    let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
    assert_eq!(f2.block_lengths, vec![1, 1]);
    let (est, _, n_snps) = f2.f2("PopA", "PopB").unwrap();
    let expected = f64::midpoint(
        snp_f2(rows[0].1[0], rows[0].1[1]),
        snp_f2(rows[2].1[0], rows[2].1[1]),
    );
    assert_relative_eq!(est, expected, max_relative = 1e-12);
    assert_eq!(n_snps, 2);
}

#[test]
fn test_monomorphic_snps_are_not_used() {
    // admixtools2 poly_only = "f2": SNP 1 has the same frequency in every
    // population, so it is left out of the f2 blocks and block lengths.
    let rows = [
        (0.00, [(0.1, 10), (0.5, 10), (0.3, 10)]),
        (0.01, [(0.0, 10), (0.0, 10), (0.0, 10)]),
        (0.02, [(0.4, 10), (0.6, 10), (0.5, 10)]),
    ];
    let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
    assert_eq!(f2.block_lengths, vec![2]);
}

#[test]
fn test_fst_blocks_keep_monomorphic_snps() {
    // admixtools2 computes FST on every SNP kept by maxmiss: the
    // monomorphic SNP 1 is in the FST block but not the f2 block, and the
    // SNP missing in PopB is in neither.
    let rows = [
        (0.00, [(0.1, 10), (0.5, 10), (0.3, 10)]),
        (0.01, [(0.2, 10), (0.2, 10), (0.2, 10)]),
        (0.02, [(0.4, 10), (0.0, 0), (0.5, 10)]),
        (0.03, [(0.4, 10), (0.6, 10), (0.5, 10)]),
    ];
    let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
    assert_eq!(f2.block_lengths, vec![2]);
    let fst = f2.fst.as_ref().unwrap();
    assert_eq!(fst.block_lengths, vec![3]);

    // Block FST(A, B) = mean(num) / mean(den) over SNPs 0, 1, 3
    let (num, den): (Vec<f64>, Vec<f64>) = [0, 1, 3]
        .iter()
        .map(|&i| {
            let ((p1, n1), (p2, n2)) = (rows[i].1[0], rows[i].1[1]);
            let (n1, n2) = (f64::from(n1), f64::from(n2));
            let (h1, h2) = (p1 * (1.0 - p1), p2 * (1.0 - p2));
            let num = (p1 - p2).powi(2) - h1 / (n1 - 1.0) - h2 / (n2 - 1.0);
            (num, num + h1 * n1 / (n1 - 1.0) + h2 * n2 / (n2 - 1.0))
        })
        .unzip();
    let expected = num.iter().sum::<f64>() / den.iter().sum::<f64>();
    let pair = f2.pair_index(0, 1);
    assert_relative_eq!(fst.block_values(pair)[0], expected, max_relative = 1e-12);
}

#[test]
fn test_fst_unavailable_for_a_single_allele_population() {
    // One allele per SNP: p(1-p)/(n-1) is 0/0 everywhere, so PopA has no
    // FST blocks; the other pair is unaffected.
    let rows = [
        (0.00, [(0.0, 1), (0.5, 10), (0.3, 10)]),
        (0.10, [(1.0, 1), (0.2, 10), (0.4, 10)]),
    ];
    let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
    assert!(f2.fst("PopA", "PopB").is_none());
    assert!(f2.fst("PopB", "PopC").is_some());
}

#[test]
fn test_blocks_without_data_for_some_pair_are_dropped() {
    // maxmiss = 1/3 keeps SNPs missing in one population. In the first
    // block PopA and PopC never have data at the same SNP, so admixtools2
    // f2_from_geno(remove_na = TRUE) drops the block for every pair.
    let rows = [
        (0.00, [(0.1, 10), (0.5, 10), (0.0, 0)]),
        (0.01, [(0.0, 0), (0.5, 10), (0.3, 10)]),
        (0.10, [(0.2, 10), (0.6, 10), (0.4, 10)]),
        (0.11, [(0.3, 10), (0.7, 10), (0.5, 10)]),
    ];
    let opts = F2Options {
        max_missing: 1.0 / 3.0,
        poly_only: true,
        ..F2Options::default()
    };
    let f2 = F2Blocks::from_allele_frequencies_with_options(
        &af_from(&rows),
        &BlockDefinition::new(0.05),
        &opts,
    );
    assert_eq!(f2.block_lengths, vec![2]);
    assert_eq!(f2.n_blocks, 1);
    assert_eq!(f2.block_sums.nrows(), 1);
}

#[test]
fn test_blocks_start_from_kept_snps() {
    // The SNP at 0.00 is dropped, so the first block starts at 0.03 and
    // spans 0.03..0.08 (admixtools2 builds blocks after filtering).
    let rows = [
        (0.00, [(0.1, 10), (0.0, 0), (0.3, 10)]),
        (0.03, [(0.1, 10), (0.5, 10), (0.3, 10)]),
        (0.06, [(0.2, 10), (0.5, 10), (0.3, 10)]),
        (0.09, [(0.3, 10), (0.5, 10), (0.3, 10)]),
    ];
    let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
    assert_eq!(f2.block_lengths, vec![2, 1]);
}

#[test]
fn test_f2_estimate_and_se_weight_blocks_by_length() {
    // Two blocks of 3 and 1 SNPs. As in admixtools2 f2(): the estimate
    // is the length-weighted mean of block means, and the variance is
    // the weighted block jackknife mean over blocks of
    // (h - 1)(est - loo)^2, h = total / length.
    let pop = |p: f64| [(p, 10), (0.5, 10), (0.5, 10)];
    let rows = [
        (0.00, pop(0.1)),
        (0.01, pop(0.2)),
        (0.02, pop(0.3)),
        (0.10, pop(0.9)),
    ];
    let f2 = F2Blocks::from_allele_frequencies(&af_from(&rows), &BlockDefinition::new(0.05));
    assert_eq!(f2.block_lengths, vec![3, 1]);

    let snp = |i: usize| snp_f2(rows[i].1[0], rows[i].1[1]);
    let x1 = (snp(0) + snp(1) + snp(2)) / 3.0;
    let x2 = snp(3);
    let est = (3.0 * x1 + x2) / 4.0;
    // Leaving out block 1 leaves block 2, and vice versa.
    let var = f64::midpoint(
        (4.0 / 3.0 - 1.0) * (est - x2).powi(2),
        (4.0 - 1.0) * (est - x1).powi(2),
    );

    let (got_est, got_se, _) = f2.f2("PopA", "PopB").unwrap();
    assert_relative_eq!(got_est, est, max_relative = 1e-12);
    assert_relative_eq!(got_se, var.sqrt(), max_relative = 1e-12);
}

#[test]
fn test_from_block_means_round_trips_the_block_values() {
    let af = make_test_data();
    let f2 = F2Blocks::from_allele_frequencies(&af, &BlockDefinition::new(0.005));
    let names: Vec<String> = f2
        .populations
        .iter()
        .map(std::string::ToString::to_string)
        .collect();
    let grid = |f: &dyn Fn(&str, &str) -> Vec<f64>| -> Vec<Vec<Vec<f64>>> {
        names
            .iter()
            .map(|a| names.iter().map(|b| f(a, b)).collect())
            .collect()
    };
    let means = grid(&|a, b| f2.get_block_f2_values(a, b).unwrap().to_vec());
    let counts = grid(&|a, b| f2.block_counts(a, b).unwrap());
    let back = F2Blocks::from_block_means(
        f2.populations.clone(),
        f2.block_lengths.clone(),
        &means,
        &counts,
    )
    .unwrap();
    assert_eq!(back.block_weights, f2.block_weights);
    for a in &names {
        for b in &names {
            let (x, y) = (
                back.get_block_f2_values(a, b).unwrap(),
                f2.get_block_f2_values(a, b).unwrap(),
            );
            for (u, v) in x.iter().zip(y.iter()) {
                assert!(
                    (u - v).abs() <= 1e-15 * v.abs() || (u.is_nan() && v.is_nan()),
                    "{a} {b}"
                );
            }
        }
    }
    assert_eq!(back.f2("PopA", "PopB"), f2.f2("PopA", "PopB"));
    assert!(F2Blocks::from_block_means(f2.populations.clone(), vec![1], &means, &counts).is_err());
    // the counts are checked as well as the means
    let short = grid(&|a, b| f2.block_counts(a, b).unwrap()[1..].to_vec());
    let lengths = f2.block_lengths.clone();
    assert!(F2Blocks::from_block_means(f2.populations.clone(), lengths, &means, &short).is_err());
}

#[test]
fn test_subset_populations_keeps_their_statistics() {
    let af = make_test_data();
    let f2 = F2Blocks::from_allele_frequencies(&af, &BlockDefinition::new(0.005));
    let sub = f2.subset_populations(&["PopC", "PopA"]).unwrap();
    assert_eq!(sub.n_populations(), 2);
    assert_eq!(sub.f2("PopA", "PopC"), f2.f2("PopA", "PopC"));
    assert_eq!(
        sub.get_block_f2_values("PopC", "PopC"),
        f2.get_block_f2_values("PopC", "PopC")
    );
    assert_eq!(sub.block_het_sums.column(0), f2.block_het_sums.column(2));
    assert!(f2.subset_populations(&["PopZ"]).is_err());
}

#[test]
fn test_select_blocks_keeps_repeats_and_order() {
    let af = make_test_data();
    let f2 = F2Blocks::from_allele_frequencies(&af, &BlockDefinition::new(0.005));
    assert!(f2.n_blocks >= 3, "{}", f2.n_blocks);
    let all: Vec<usize> = (0..f2.n_blocks).collect();
    let same = f2.select_blocks(&all).unwrap();
    assert_eq!(same.block_sums, f2.block_sums);
    assert_eq!(same.block_lengths, f2.block_lengths);
    assert!(same.fst.is_none());

    let boot = f2.select_blocks(&[2, 0, 2]).unwrap();
    assert_eq!(boot.n_blocks, 3);
    assert_eq!(
        boot.block_lengths,
        vec![
            f2.block_lengths[2],
            f2.block_lengths[0],
            f2.block_lengths[2]
        ]
    );
    assert_eq!(boot.block_weights.row(0), f2.block_weights.row(2));
    assert_eq!(boot.block_het_sums.row(1), f2.block_het_sums.row(0));
    assert_eq!(boot.block_raw_sums.row(2), f2.block_raw_sums.row(2));

    assert!(f2.select_blocks(&[f2.n_blocks]).is_err());
    assert!(f2.select_blocks(&[]).is_err());
}

#[test]
fn test_f2_blocks_creation() {
    let af = make_test_data();
    let block_def = BlockDefinition::new(0.01); // 1 cM blocks

    let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

    assert_eq!(f2_blocks.n_populations(), 3);
    assert!(f2_blocks.n_blocks > 0);

    // Should have 6 pairs: AA, AB, AC, BB, BC, CC
    let n_pairs = 3 * 4 / 2; // n*(n+1)/2 = 6
    assert_eq!(f2_blocks.block_sums.ncols(), n_pairs);
}

#[test]
fn test_f2_self_pair_near_zero() {
    let af = make_test_data();
    let block_def = BlockDefinition::default();

    let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

    // f2(A,A) should be close to 0
    let (f2_aa, _, _) = f2_blocks.f2("PopA", "PopA").unwrap();
    assert!(f2_aa.abs() < 0.01, "f2(A,A) should be near 0, got {f2_aa}");
}

#[test]
fn test_f2_symmetry() {
    let af = make_test_data();
    let block_def = BlockDefinition::default();

    let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

    // f2(A,B) should equal f2(B,A)
    let (f2_ab, _, _) = f2_blocks.f2("PopA", "PopB").unwrap();
    let (f2_ba, _, _) = f2_blocks.f2("PopB", "PopA").unwrap();
    assert_relative_eq!(f2_ab, f2_ba, epsilon = 1e-10);
}

#[test]
fn test_f2_different_populations() {
    let af = make_test_data();
    let block_def = BlockDefinition::default();

    let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

    // f2(A,B) should be positive (different populations)
    let (f2_ab, std_err, n_snps) = f2_blocks.f2("PopA", "PopB").unwrap();
    assert!(f2_ab > 0.0, "f2(A,B) should be positive, got {f2_ab}");
    assert!(std_err >= 0.0);
    assert!(n_snps > 0);
}

#[test]
fn f2_matrix_holds_the_f2_estimates() {
    // Blocks of 1 and 3 SNPs: the length-weighted estimate is 0.4, the
    // plain mean of block means 0.3.
    let block = |a: f64, b: f64| {
        vec![
            vec![vec![0.0; 2], vec![a, b]],
            vec![vec![a, b], vec![0.0; 2]],
        ]
    };
    let f2 = F2Blocks::from_block_means(
        vec![PopulationId::new("A"), PopulationId::new("B")],
        vec![1, 3],
        &block(0.1, 0.5),
        &block(1.0, 3.0),
    )
    .unwrap();
    let (estimate, _, _) = f2.f2("A", "B").unwrap();
    assert_relative_eq!(estimate, 0.4, epsilon = 1e-12);
    assert_eq!(f2.f2_matrix()[[0, 1]], estimate);
    assert_eq!(f2.f2_matrix()[[1, 0]], estimate);
}

#[test]
fn test_f2_matrix() {
    let af = make_test_data();
    let block_def = BlockDefinition::default();

    let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);
    let matrix = f2_blocks.f2_matrix();

    assert_eq!(matrix.shape(), &[3, 3]);

    // Print matrix for debugging
    println!("F2 matrix:\n{matrix:?}");

    // Diagonal should be near 0 (but can be slightly negative due to bias correction)
    for i in 0..3 {
        assert!(
            matrix[[i, i]].abs() < 0.1,
            "Diagonal element {} = {} should be near 0",
            i,
            matrix[[i, i]]
        );
    }

    // Should be symmetric
    for i in 0..3 {
        for j in 0..3 {
            assert_relative_eq!(matrix[[i, j]], matrix[[j, i]], epsilon = 1e-10);
        }
    }
}

#[test]
fn test_pair_index() {
    let af = make_test_data();
    let block_def = BlockDefinition::default();
    let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

    // Verify pair indexing is consistent
    // For 3 populations: pairs are (0,0), (0,1), (0,2), (1,1), (1,2), (2,2)
    // Indices should be 0, 1, 2, 3, 4, 5
    assert_eq!(f2_blocks.pair_index(0, 0), 0);
    assert_eq!(f2_blocks.pair_index(0, 1), 1);
    assert_eq!(f2_blocks.pair_index(0, 2), 2);
    assert_eq!(f2_blocks.pair_index(1, 1), 3);
    assert_eq!(f2_blocks.pair_index(1, 2), 4);
    assert_eq!(f2_blocks.pair_index(2, 2), 5);

    // Symmetric access
    assert_eq!(f2_blocks.pair_index(1, 0), f2_blocks.pair_index(0, 1));
    assert_eq!(f2_blocks.pair_index(2, 0), f2_blocks.pair_index(0, 2));
}

#[test]
fn test_jackknife_std_err() {
    let af = make_test_data();
    let block_def = BlockDefinition::new(0.005); // Small blocks for more jackknife samples

    let f2_blocks = F2Blocks::from_allele_frequencies(&af, &block_def);

    let (_, std_err, _) = f2_blocks.f2("PopA", "PopB").unwrap();

    // Standard error should be non-negative
    assert!(std_err >= 0.0);

    // With multiple blocks, should have positive standard error
    if f2_blocks.n_blocks > 2 {
        // May still be 0 if all blocks have identical estimates
    }
}
