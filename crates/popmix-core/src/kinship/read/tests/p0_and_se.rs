//! The P0, overlap and standard error fixtures of the specification's
//! section 12 (READv2 run as a black box), and direct counts.

use super::*;

#[test]
fn e6_windows_bin_absolute_positions_and_recurring_bins_split() {
    let pos = [99, 100, 101, 150, 199, 200, 250, 301];
    let plain = ReadCompat {
        windows_overlap: false,
        ..ReadCompat::V2_01
    };
    assert_eq!(
        windows(&pos, &[1; 8], 100, plain),
        vec![(0, 0), (1, 4), (5, 6), (7, 7)]
    );
    // bins 0, 0, 1, 0: three windows
    assert_eq!(
        windows(&[1, 2, 101, 3], &[1; 4], 100, ReadCompat::CORRECTED),
        vec![(0, 1), (2, 2), (3, 3)]
    );
    // the chromosome case: only the bin splits READv2's windows
    let pos = [10, 20, 30, 40, 150, 160, 50, 60];
    let chr = [1, 1, 2, 2, 2, 3, 3, 3];
    assert_eq!(
        windows(&pos, &chr, 100, plain),
        vec![(0, 3), (4, 5), (6, 7)]
    );
    assert_eq!(
        windows(&pos, &chr, 100, ReadCompat::V2_01),
        vec![(0, 3), (3, 5), (5, 7)]
    );
    assert_eq!(
        windows(&pos, &chr, 100, ReadCompat::CORRECTED),
        vec![(0, 1), (2, 3), (4, 4), (5, 5), (6, 7)]
    );
    assert!(windows(&[], &[], 100, ReadCompat::V2_01).is_empty());
}

#[test]
fn jackknife_edge_cases() {
    assert!(jackknife_se(&[], false).is_nan());
    assert!(jackknife_se(&[], true).is_nan());
    assert!(jackknife_se(&[0.3], false).is_nan());
    assert_eq!(jackknife_se(&[0.3], true), 0.0);
    assert_eq!(jackknife_se(&[0.1, 0.7], true), 0.0);
    // e2a: two left out of three windows, factor (g - 1) / g with g = 3
    assert_close(
        jackknife_se(&[0.25, 0.5, 1.0 / 3.0], true),
        0.07216878364870322,
        1e-14,
    );
}

const E2_POS: [u64; 12] = [1, 2, 3, 4, 11, 12, 13, 21, 22, 23, 24, 25];

#[test]
fn e2_jackknife_of_window_means_with_overlap_and_last_window_kept() {
    // (mismatching SNPs, P0, READv2's SE)
    let cases: [(&[usize], f64, f64); 4] = [
        (&[0, 4, 5, 11], 1.0 / 3.0, 0.07216878364870322),
        (&[4, 5, 11], 0.25, 0.14433756729740643),
        (&[3, 4, 5, 11], 1.0 / 3.0, 0.14433756729740643),
        (&[3], 1.0 / 12.0, 0.0),
    ];
    for (mism, p0, se) in cases {
        let r = run(&two_samples(12, mism), &E2_POS, &v2_half(10));
        let p = &r.pairs[0];
        assert_close(p.p0, p0, 1e-15);
        assert_close(p.p0_se, se, 1e-14);
        assert_eq!(p.overlap_snps, 12);
        assert_eq!(p.effective_snps, 6.0);
        assert_eq!(p.normalized_p0, p.p0 / 0.5);
        assert_eq!(p.normalized_p0_se, p.p0_se / 0.5);
    }
    // e2a CORRECTED (section 12): windows {1-4}, {11-13}, {21-25}
    let corrected = ReadOptions {
        compat: ReadCompat::CORRECTED,
        ..v2_half(10)
    };
    let r = run(&two_samples(12, &[0, 4, 5, 11]), &E2_POS, &corrected);
    assert_close(r.pairs[0].p0, 1.0 / 3.0, 1e-15);
    assert_close(r.pairs[0].p0_se, 0.14792807728549273, 1e-14);
}

#[test]
fn e4_a_mismatch_at_a_window_end_moves_the_se_and_in_the_last_window_not() {
    let pos: Vec<u64> = (0..4)
        .flat_map(|w| (0..3).map(move |k| 100 * w + 1 + k))
        .collect();
    let expected = [
        0.07856742013183861,
        0.07856742013183861,
        0.07081971546656644,
        0.05892556509887896,
        0.05892556509887896,
        0.05892556509887896,
        0.05892556509887896,
        0.05892556509887896,
        0.05892556509887896,
        0.0,
        0.0,
        0.0,
    ];
    for (k, se) in expected.into_iter().enumerate() {
        let r = run(&two_samples(12, &[k]), &pos, &v2_half(100));
        let p = &r.pairs[0];
        assert_close(p.p0, 1.0 / 12.0, 1e-15);
        assert_close(p.p0_se, se, 1e-14);
        if se == 0.0 {
            // se_n = 0: z inside Some, inf; the interval collapses
            assert_eq!(p.z_upper, Some(f64::INFINITY));
            assert_eq!(p.kinship_ci, (p.normalized_p0, p.normalized_p0));
        }
    }
    // CORRECTED, mismatch only in the last window: sd([0, 0, 0, 1/3]) / 2
    let corrected = ReadOptions {
        compat: ReadCompat::CORRECTED,
        ..v2_half(100)
    };
    let r = run(&two_samples(12, &[10]), &pos, &corrected);
    assert_close(r.pairs[0].p0_se, 1.0 / 12.0, 1e-14);
}

#[test]
fn e5_windows_of_unequal_size() {
    let sizes = [2u64, 3, 4, 5];
    let pos: Vec<u64> = sizes
        .iter()
        .enumerate()
        .flat_map(|(w, &s)| (0..s).map(move |k| 100 * w as u64 + 1 + k))
        .collect();
    let mut expected = vec![
        0.11785113019775792,
        0.10206207261596575,
        0.05892556509887896,
        0.05892556509887896,
        0.05400617248673217,
    ];
    expected.extend([0.04714045207910317; 4]);
    expected.extend([0.0; 5]);
    for (k, se) in expected.into_iter().enumerate() {
        let r = run(&two_samples(14, &[k]), &pos, &v2_half(100));
        assert_close(r.pairs[0].p0, 1.0 / 14.0, 1e-15);
        assert_close(r.pairs[0].p0_se, se, 1e-14);
    }
}

#[test]
fn e1_random_samples_match_readv2() {
    let g = from_strs(&[
        "02200002220020002220000202020022022000222002022000222202200000202202220220000002022000000222220022220020200002202222020220000002020020020202222220002220002002220020202202000020202220002200222202200220",
        "20222202000000002000202202022022002200200222022022002002220222020220200220022002220200022020220022222222022222020202000222222022000200022002002200200220022220220022022000000022202200222202202022002022",
        "20000222002022020222020220022222222202202200222022020220020020200020022000222000020000222222002000200022202222020220202020000202200002220200220222222020000202002000220020000000202000000020202220022000",
    ]);
    let pos: Vec<u64> = (0..200).map(|i| 1_000_000 * i + 1).collect();
    let opts = ReadOptions {
        compat: ReadCompat::V2_01,
        ..Default::default()
    };
    let r = run(&g, &pos, &opts);
    assert_eq!(r.norm_value, 0.52);
    let expected = [
        ((0, 1), 0.46, 0.03454620077970638, 0.8),
        ((0, 2), 0.52, 0.025480982009181707, 0.9),
        ((1, 2), 0.55, 0.03360398558469755, 1.0),
    ];
    for ((i, j), p0, se, q) in expected {
        let p = pair(&r, i, j);
        assert_eq!(p.p0, p0);
        assert_close(p.p0_se, se, 1e-13);
        assert_eq!(p.first_degree_windows_outside, q);
        assert_eq!(p.overlap_snps, 200);
        assert_eq!(p.effective_snps, 104.0);
        assert_eq!(p.first_degree_type, None);
    }
    let p = pair(&r, 0, 1);
    assert_eq!(p.relationship, ReadRelationship::SecondDegree);
    assert_close(p.z_upper.unwrap(), 0.3256508601839845, 1e-12);
    assert_close(p.z_lower.unwrap(), 1.085502867279946, 1e-12);
    assert_close(p.kinship, 0.11538461538461542, 1e-14);
    assert_close(p.kinship_ci.0, 0.7544027816764913, 1e-13);
    assert_close(p.kinship_ci.1, 1.0148279875542778, 1e-13);
    // unrelated pairs: no z_upper, z_lower to t1 (E < 3000)
    let p = pair(&r, 1, 2);
    assert_eq!(p.relationship, ReadRelationship::Unrelated);
    assert_eq!(p.z_upper, None);
    assert_close(p.z_lower.unwrap(), 2.343472020648077, 1e-12);
}

#[test]
fn e7_missing_calls_heterozygotes_and_a_window_with_nothing_in_common() {
    let g = from_strs(&[
        "220202202202022",
        "22....022202002",
        ".20000212020020",
        "202200021200222",
    ]);
    let pos: Vec<u64> = (0..5)
        .flat_map(|w| (0..3).map(move |k| 100 * w + 1 + k))
        .collect();
    let r = run(&g, &pos, &v2_half(100));
    let expected = [
        ((0, 1), 0.2727272727272727, 0.1571348402636772, 11),
        ((0, 2), 0.4642857142857143, 0.12103072956898178, 14),
        ((0, 3), 0.5, 0.08385254915624211, 15),
        ((1, 2), 0.65, 0.15590239111558085, 10),
        ((1, 3), 0.4090909090909091, 0.06874649261535878, 11),
        ((2, 3), 0.6428571428571429, 0.09164298936634488, 14),
    ];
    let order: Vec<(usize, usize)> = r.pairs.iter().map(|p| (p.i, p.j)).collect();
    assert_eq!(order, expected.map(|e| e.0).to_vec());
    for ((i, j), p0, se, overlap) in expected {
        let p = pair(&r, i, j);
        assert_eq!(p.p0, p0);
        assert_close(p.p0_se, se, 1e-13);
        assert_eq!(p.overlap_snps, overlap);
    }
    let p = pair(&r, 0, 2);
    assert_eq!(p.relationship, ReadRelationship::UnrelatedOrThirdDegree);
    assert_close(p.z_upper.unwrap(), 0.10143527811495602, 1e-12);
    assert_close(p.z_lower.unwrap(), 0.09221388919541482, 1e-12);
    assert_close(p.kinship_ci.0, 0.45413096866102004, 1e-13);
    assert_close(p.kinship_ci.1, 1.4030118884818372, 1e-13);
    let p = pair(&r, 0, 1);
    assert_eq!(p.relationship, ReadRelationship::IdenticalOrTwins);
    assert_eq!(p.z_lower, None);
    assert_close(p.z_upper.unwrap(), 0.25311208644745753, 1e-12);
    assert_eq!(pair(&r, 1, 3).first_degree_windows_outside, 0.0);
}

#[test]
fn e8_first_degree_windows_overlap_too_and_none_is_skipped() {
    let w = 20_000_000;
    let pos = [1, 2, w + 1, w + 2, 2 * w + 1, 2 * w + 2];
    let opts = v2_half(5_000_000);
    let r = run(&two_samples(6, &[1]), &pos, &opts);
    assert_close(r.pairs[0].first_degree_windows_outside, 2.0 / 3.0, 1e-15);
    assert_close(r.pairs[0].p0_se, 0.048112522432468816, 1e-13);
    // a mismatch at the first SNP of the last window: window P0s 0, 0, 1/3
    // (the last window holds SNPs 3-5), y = 0, 0, 2/3: one of three inside
    let r = run(&two_samples(6, &[4]), &pos, &opts);
    assert_close(r.pairs[0].first_degree_windows_outside, 2.0 / 3.0, 1e-15);
    // CORRECTED windows do not overlap: y = 1, 0, 0, all outside
    let corrected = ReadOptions {
        compat: ReadCompat::CORRECTED,
        ..opts
    };
    let r = run(&two_samples(6, &[1]), &pos, &corrected);
    assert_eq!(r.pairs[0].first_degree_windows_outside, 1.0);
}

#[test]
fn e16w_window_mode_counts_overlap_snps_twice_and_drops_empty_windows() {
    let g = from_strs(&[
        "022220202202000200022002200000222202022222202022000002200220",
        "220222200.......20000220022020202200202200000020202020222000",
        "022222020200202202002022202222000220002002000202002022222200",
    ]);
    let pos: Vec<u64> = (0..60).map(|i| 1_000_000 * i + 1).collect();
    let opts = ReadOptions {
        window_estimate: true,
        ..v2_half(5_000_000)
    };
    let r = run(&g, &pos, &opts);
    let expected = [
        ((0, 1), 0.5424242424242425, 0.06357861737797649, 62),
        ((0, 2), 0.5138888888888888, 0.06309342297297468, 71),
        ((1, 2), 0.4666666666666666, 0.05790971726347274, 62),
    ];
    for ((i, j), p0, se, overlap) in expected {
        let p = pair(&r, i, j);
        assert_close(p.p0, p0, 1e-14);
        assert_close(p.p0_se, se, 1e-13);
        assert_eq!(p.overlap_snps, overlap);
        assert_eq!(p.effective_snps, overlap as f64 * 0.5);
    }
    assert_close(
        pair(&r, 0, 2).first_degree_windows_outside,
        2.0 / 3.0,
        1e-15,
    );
    // without overlap, each SNP once
    let corrected = ReadOptions {
        compat: ReadCompat::CORRECTED,
        ..opts
    };
    assert_eq!(pair(&run(&g, &pos, &corrected), 0, 2).overlap_snps, 60);
}

#[test]
fn counts_across_words_match_a_direct_count() {
    let m = 300;
    let strs = random_strs(4, m, 3);
    let refs: Vec<&str> = strs.iter().map(String::as_str).collect();
    let pos: Vec<u64> = (0..m as u64).map(|i| 37 * i * i + 5).collect();
    let opts = ReadOptions {
        window_size: 50_000,
        norm: ReadNorm::Mean,
        ..Default::default()
    };
    let r = run(&from_strs(&refs), &pos, &opts);
    let wins = windows(&pos, &vec![1; m], 50_000, ReadCompat::CORRECTED);
    assert!(wins.len() > 5);
    for p in &r.pairs {
        let (a, b) = (refs[p.i].as_bytes(), refs[p.j].as_bytes());
        let d = |k: usize| -> Option<f64> {
            match (a[k], b[k]) {
                (b'.', _) | (_, b'.') => None,
                (b'1', _) | (_, b'1') => Some(0.5),
                (x, y) => Some(if x == y { 0.0 } else { 1.0 }),
            }
        };
        let p0_of = |s: usize, e: usize| {
            let v: Vec<f64> = (s..=e).filter_map(d).collect();
            (!v.is_empty()).then(|| v.iter().sum::<f64>() / v.len() as f64)
        };
        let n = (0..m).filter_map(d).count() as u64;
        assert_eq!(p.overlap_snps, n);
        assert_close(p.p0, p0_of(0, m - 1).unwrap(), 1e-15);
        let wp: Vec<f64> = wins.iter().filter_map(|&(s, e)| p0_of(s, e)).collect();
        assert_close(p.p0_se, jackknife_se(&wp, false), 1e-12);
    }
    let p0s: Vec<f64> = r.pairs.iter().map(|p| p.p0).collect();
    assert_eq!(r.norm_value, mean(&p0s));
}
