use super::*;
use crate::fstats::MatrixRows;
use crate::genotype_matrix::GenotypeMatrix;
use crate::genotypes::{Chromosome, Sample, Snp};
use ndarray::Array2;

fn bits(geno: &Array2<u8>) -> BitGenotypes {
    let (n_snps, n) = geno.dim();
    let samples: Vec<Sample> = (0..n).map(|i| Sample::new(format!("s{i}"), "P")).collect();
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
    let m = GenotypeMatrix::new(geno.clone(), samples, snps).unwrap();
    let mut rows = MatrixRows::new(&m, (0..n_snps).collect(), (0..n).collect());
    BitGenotypes::read(&mut rows, n, n_snps).unwrap()
}

#[test]
fn overlapping_windows_share_their_boundary_snp() {
    let pos = [1, 2, 6, 7, 11];
    let chr = [1; 5];
    assert_eq!(
        windows(&pos, &chr, 5, ReadCompat::CORRECTED),
        vec![(0, 1), (2, 3), (4, 4)]
    );
    assert_eq!(
        windows(&pos, &chr, 5, ReadCompat::V2_01),
        vec![(0, 1), (1, 3), (3, 4)]
    );
}

#[test]
fn chromosome_boundaries_split_windows() {
    let pos = [1, 2, 3, 4];
    let chr = [1, 1, 2, 2];
    assert_eq!(
        windows(&pos, &chr, 100, ReadCompat::CORRECTED),
        vec![(0, 1), (2, 3)]
    );
    assert_eq!(windows(&pos, &chr, 100, ReadCompat::V2_01), vec![(0, 3)]);
}

#[test]
fn jackknife_of_equal_windows_is_zero_and_skipping_one_changes_it() {
    assert_eq!(jackknife_se(&[0.2, 0.2, 0.2], false), 0.0);
    let x = [0.1, 0.2, 0.4, 0.3];
    // leave-one-out means 0.3, 0.2667, 0.2, 0.2333: the usual jackknife
    let full = jackknife_se(&x, false);
    let sd = {
        let m = mean(&x);
        (x.iter().map(|v| (v - m).powi(2)).sum::<f64>() / 3.0).sqrt()
    };
    assert!((full - sd / 2.0).abs() < 1e-15);
    assert!((jackknife_se(&x, true) - full).abs() > 1e-3);
}

#[test]
fn identical_samples_are_identical_and_halves_are_counted() {
    // pseudo-haploid samples: 0 and 1 identical, 2 differs at half the SNPs,
    // 3 has a heterozygote
    let n = 4000;
    let geno = Array2::from_shape_fn((n, 4), |(s, k)| match k {
        0 | 1 => 2 * ((s / 3) % 2) as u8,
        2 => 2 * ((s / 7) % 2) as u8,
        _ => {
            if s == 0 {
                1
            } else {
                2 * ((s / 5) % 2) as u8
            }
        }
    });
    let pos: Vec<u64> = (0..n as u64).map(|s| s * 1000).collect();
    let opts = ReadOptions {
        norm: ReadNorm::Value(0.5),
        ..Default::default()
    };
    let r = read_kinship(&bits(&geno), &pos, &vec![1; n], &opts).unwrap();
    let p = |i, j| r.pairs.iter().find(|p| p.i == i && p.j == j).unwrap();
    assert_eq!(p(0, 1).p0, 0.0);
    assert_eq!(p(0, 1).relationship, ReadRelationship::IdenticalOrTwins);
    assert_eq!(p(0, 1).kinship, 1.0);
    // the heterozygote at SNP 0 counts half a mismatch
    let mism = (1..n).filter(|&s| geno[[s, 0]] != geno[[s, 3]]).count() as f64 + 0.5;
    assert!((p(0, 3).p0 - mism / n as f64).abs() < 1e-15);
}

#[test]
fn kinship_interval_is_around_the_kinship_unless_compat() {
    let n = 3000;
    let geno = Array2::from_shape_fn((n, 3), |(s, k)| 2 * (((s * (k + 3)) / 11) % 2) as u8);
    let pos: Vec<u64> = (0..n as u64).map(|s| s * 10_000).collect();
    let corrected = read_kinship(&bits(&geno), &pos, &vec![1; n], &ReadOptions::default()).unwrap();
    let p = &corrected.pairs[0];
    assert!(p.kinship_ci.0 < p.kinship && p.kinship < p.kinship_ci.1);
    let compat = ReadOptions {
        compat: ReadCompat::V2_01,
        ..Default::default()
    };
    let v = read_kinship(&bits(&geno), &pos, &vec![1; n], &compat).unwrap();
    let q = &v.pairs[0];
    assert!(q.kinship_ci.0 < q.normalized_p0 && q.normalized_p0 < q.kinship_ci.1);
}

// ---- fixtures from the specification's black-box runs (section 12) ----

/// Samples as strings of calls, one character per SNP: 0, 1, 2, or `.`
/// for missing.
fn from_strs(samples: &[&str]) -> BitGenotypes {
    let m = samples[0].len();
    let geno = Array2::from_shape_fn((m, samples.len()), |(s, k)| {
        match samples[k].as_bytes()[s] {
            b'0' => 0,
            b'1' => 1,
            b'2' => 2,
            _ => 3,
        }
    });
    bits(&geno)
}

/// Two pseudo-haploid samples, all 0, the second 2 at `mismatches`.
fn two_samples(m: usize, mismatches: &[usize]) -> BitGenotypes {
    let a = "0".repeat(m);
    let b: String = (0..m)
        .map(|k| if mismatches.contains(&k) { '2' } else { '0' })
        .collect();
    from_strs(&[&a, &b])
}

fn pair(r: &ReadResult, i: usize, j: usize) -> &ReadPair {
    r.pairs.iter().find(|p| p.i == i && p.j == j).unwrap()
}

fn assert_close(actual: f64, expected: f64, rtol: f64) {
    assert!(
        (actual - expected).abs() <= rtol * expected.abs(),
        "{actual} != {expected} (rtol {rtol})"
    );
}

/// READv2 with `--window_size w -n value --norm_value 0.5`.
fn v2_half(w: u64) -> ReadOptions {
    ReadOptions {
        window_size: w,
        norm: ReadNorm::Value(0.5),
        compat: ReadCompat::V2_01,
        ..Default::default()
    }
}

fn run(g: &BitGenotypes, pos: &[u64], opts: &ReadOptions) -> ReadResult {
    read_kinship(g, pos, &vec![1; pos.len()], opts).unwrap()
}

/// Random pseudo-haploid samples with missing calls and heterozygotes.
fn random_strs(n: usize, m: usize, seed: u64) -> Vec<String> {
    let mut state = seed;
    let mut next = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) % 20
    };
    (0..n)
        .map(|_| {
            (0..m)
                .map(|_| match next() {
                    0 => '.',
                    1 => '1',
                    v if v % 2 == 0 => '0',
                    _ => '2',
                })
                .collect()
        })
        .collect()
}

mod classes;
mod p0_and_se;
