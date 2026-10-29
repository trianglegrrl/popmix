//! The classification, first-degree, normalisation and error fixtures of
//! the specification's section 12 (READv2 run as a black box).

use super::*;

/// e9: 6016 SNPs, norm 0.5, 47 k mismatches: normalized P0 = k / 64,
/// E = 3008.
fn e9(ks: &[usize], thresholds: [f64; 4]) -> ReadResult {
    let m = 6016;
    let mut cols = vec!["0".repeat(m)];
    for &k in ks {
        cols.push((0..m).map(|s| if s < 47 * k { '2' } else { '0' }).collect());
    }
    let refs: Vec<&str> = cols.iter().map(String::as_str).collect();
    let pos: Vec<u64> = (0..m as u64).map(|i| 16_000 * i + 1).collect();
    let opts = ReadOptions {
        thresholds,
        ..v2_half(5_000_000)
    };
    run(&from_strs(&refs), &pos, &opts)
}

#[test]
fn e9_classes_at_the_thresholds() {
    use ReadRelationship::*;
    let cases = [
        (39, IdenticalOrTwins),
        (40, FirstDegree),
        (51, FirstDegree),
        (52, SecondDegree),
        (57, SecondDegree),
        (58, ThirdDegree),
        (61, ThirdDegree),
        (62, Unrelated),
    ];
    let ks: Vec<usize> = cases.iter().map(|c| c.0).collect();
    let r = e9(&ks, READ_THRESHOLDS);
    for (c, (k, rel)) in cases.into_iter().enumerate() {
        let p = pair(&r, 0, c + 1);
        assert_eq!(p.normalized_p0, k as f64 / 64.0);
        assert_eq!(p.effective_snps, 3008.0);
        assert_eq!(p.relationship, rel, "k = {k}");
    }
    // unrelated: z_lower to t1, not t0
    let p = pair(&r, 0, 8);
    let x = 62.0 / 64.0;
    assert_eq!(p.z_lower, Some((x - 0.90625) / p.normalized_p0_se));
    // third degree: (t0 - x, x - t1)
    let p = pair(&r, 0, 6);
    let x = 58.0 / 64.0;
    assert_eq!(p.z_upper, Some((0.953125 - x) / p.normalized_p0_se));
    assert_eq!(p.z_lower, Some((x - 0.90625) / p.normalized_p0_se));
}

#[test]
fn e9p_2pow_thresholds() {
    let t = read_thresholds_2pow();
    let expected = [
        0.9558058261758408,
        0.9116116523516815,
        0.8232233047033631,
        0.6464466094067263,
    ];
    for (a, b) in t.iter().zip(expected) {
        assert_close(*a, b, 1e-15);
    }
    use ReadRelationship::*;
    // the 2^-6-step values that bracket each threshold
    let cases = [
        (41, IdenticalOrTwins),
        (42, FirstDegree),
        (52, FirstDegree),
        (53, SecondDegree),
        (58, SecondDegree),
        (59, ThirdDegree),
        (61, ThirdDegree),
        (62, Unrelated),
    ];
    let ks: Vec<usize> = cases.iter().map(|c| c.0).collect();
    let r = e9(&ks, t);
    for (c, (k, rel)) in cases.into_iter().enumerate() {
        assert_eq!(pair(&r, 0, c + 1).relationship, rel, "k = {k}");
    }
}

/// e12: `nw` 20 Mb windows of `per` SNPs, norm 0.5; for each sample, the
/// number of windows outside [t3, t1] (alternately y = 0.5 and 1, the rest
/// 0.75), never mismatching at a window's last SNP, and SNPs missing.
fn e12(nw: usize, per: usize, specs: &[(usize, usize)]) -> ReadResult {
    let m = nw * per;
    let mut cols = vec!["0".repeat(m)];
    for &(n_out, missing) in specs {
        let mut c = vec![b'0'; m];
        for w in 0..nw {
            let y = if w < n_out {
                if w % 2 == 0 {
                    0.5
                } else {
                    1.0
                }
            } else {
                0.75
            };
            let k = (y / 2.0 * per as f64).round() as usize;
            for s in 0..k {
                c[w * per + s] = b'2';
            }
        }
        // missing calls at SNPs that do not mismatch
        for s in 0..missing {
            c[per - 2 - s] = b'.';
        }
        cols.push(String::from_utf8(c).unwrap());
    }
    let refs: Vec<&str> = cols.iter().map(String::as_str).collect();
    let pos: Vec<u64> = (0..nw)
        .flat_map(|w| (0..per).map(move |i| 20_000_000 * w as u64 + 1 + 1000 * i as u64))
        .collect();
    run(&from_strs(&refs), &pos, &v2_half(5_000_000))
}

#[test]
fn e12_first_degree_type_from_the_share_of_windows_outside() {
    use FirstDegreeType::*;
    let r = e12(
        20,
        1000,
        &[
            (5, 0),
            (6, 0),
            (7, 0),
            (10, 0),
            (12, 0),
            (13, 0),
            (2, 2),
            (2, 1),
        ],
    );
    let expected = [
        (0.25, Some(ParentOffspring)),
        (0.3, Some(ParentOffspring)),
        (0.35, Some(Siblings)),
        (0.5, Some(Siblings)),
        (0.6, Some(NotApplicable)),
        (0.65, Some(NotApplicable)),
        // E = 9999 and 9999.5: not enough data
        (0.1, None),
        (0.1, None),
    ];
    for (c, (q, ty)) in expected.into_iter().enumerate() {
        let p = pair(&r, 0, c + 1);
        assert_eq!(p.relationship, ReadRelationship::FirstDegree);
        assert_close(p.first_degree_windows_outside, q, 1e-15);
        assert_eq!(p.first_degree_type, ty, "sample {}", c + 1);
    }
    assert_eq!(pair(&r, 0, 1).effective_snps, 10_000.0);
    assert_eq!(pair(&r, 0, 7).effective_snps, 9999.0);
    assert_eq!(pair(&r, 0, 8).effective_snps, 9999.5);
    let r = e12(
        40,
        500,
        &[(11, 0), (12, 0), (13, 0), (14, 0), (23, 0), (24, 0)],
    );
    let expected = [
        (0.275, Some(ParentOffspring)),
        (0.3, Some(ParentOffspring)),
        (0.325, None),
        (0.35, Some(Siblings)),
        (0.575, Some(Siblings)),
        (0.6, Some(NotApplicable)),
    ];
    for (c, (q, ty)) in expected.into_iter().enumerate() {
        let p = pair(&r, 0, c + 1);
        assert_close(p.first_degree_windows_outside, q, 1e-15);
        assert_eq!(p.first_degree_type, ty, "sample {}", c + 1);
    }
    // not first degree: no type, but the share is reported
    let p = pair(&r, 1, 2);
    assert_ne!(p.relationship, ReadRelationship::FirstDegree);
    assert_eq!(p.first_degree_type, None);
    assert!(p.first_degree_windows_outside > 0.0);
}

#[test]
fn e13_first_degree_window_boundaries() {
    let per = 64;
    let m = 4 * per;
    let pos: Vec<u64> = (0..4)
        .flat_map(|w| (0..per).map(move |i| 20_000_000 * w as u64 + 1 + 1000 * i as u64))
        .collect();
    let first_k =
        |k: usize| -> Vec<u8> { (0..m).map(|s| if s < k { b'2' } else { b'0' }).collect() };
    let mut cols: Vec<Vec<u8>> = vec![vec![b'0'; m]];
    for k in [19, 20, 21, 28, 29, 30] {
        cols.push(first_k(k));
    }
    let mut c = first_k(20);
    c[2 * per - 1..3 * per].fill(b'.'); // window 2 and its overlap SNP missing
    cols.push(c);
    let mut c = first_k(20);
    c[3 * per..4 * per].fill(b'.'); // window 3 missing but its overlap SNP
    cols.push(c);
    let strs: Vec<String> = cols
        .into_iter()
        .map(|c| String::from_utf8(c).unwrap())
        .collect();
    let refs: Vec<&str> = strs.iter().map(String::as_str).collect();
    let g = from_strs(&refs);
    let r = run(&g, &pos, &v2_half(5_000_000));
    // window 0 at y = k/32: 0.59375 out, 0.625 in, ..., 0.90625 in, 0.9375 out
    let expected = [1.0, 0.75, 0.75, 0.75, 0.75, 1.0, 2.0 / 3.0, 0.75];
    for (c, q) in expected.into_iter().enumerate() {
        assert_close(pair(&r, 0, c + 1).first_degree_windows_outside, q, 1e-15);
    }
    let pow = ReadOptions {
        thresholds: read_thresholds_2pow(),
        ..v2_half(5_000_000)
    };
    let r = run(&g, &pos, &pow);
    // 0.625 now outside; 0.65625 and 0.90625 inside
    let expected = [1.0, 1.0, 0.75, 0.75, 0.75, 1.0];
    for (c, q) in expected.into_iter().enumerate() {
        assert_eq!(pair(&r, 0, c + 1).first_degree_windows_outside, q);
    }
}

#[test]
fn e14_one_or_two_windows_give_se_zero_in_compat() {
    let g = from_strs(&["0000000000", "2200000000"]);
    let pos: Vec<u64> = (0..10).map(|i| 1000 * i + 1).collect();
    let r = run(&g, &pos, &v2_half(5_000_000));
    let p = &r.pairs[0];
    assert_eq!(p.p0_se, 0.0);
    assert_eq!(p.relationship, ReadRelationship::IdenticalOrTwins);
    assert_eq!(p.z_upper, Some(f64::INFINITY));
    assert_eq!(p.kinship_ci, (0.4, 0.4));
    // CORRECTED leaves the one window out: 0/0
    let corrected = ReadOptions {
        compat: ReadCompat::CORRECTED,
        ..v2_half(5_000_000)
    };
    assert!(run(&g, &pos, &corrected).pairs[0].p0_se.is_nan());
    // two windows
    let pos: Vec<u64> = (0..10).map(|i| 1_000_000 * i + 1).collect();
    assert_eq!(run(&g, &pos, &v2_half(5_000_000)).pairs[0].p0_se, 0.0);
}

#[test]
fn e16_third_degree_needs_3000_effective_snps() {
    let m = 6000;
    let pos: Vec<u64> = (0..m as u64).map(|i| 16_000 * i + 1).collect();
    let b: String = (0..m).map(|s| if s < 2812 { '2' } else { '0' }).collect();
    let mut c = b.clone().into_bytes();
    c[m - 1] = b'.';
    let c = String::from_utf8(c).unwrap();
    let r = run(
        &from_strs(&["0".repeat(m).as_str(), &b, &c]),
        &pos,
        &v2_half(5_000_000),
    );
    assert_eq!(pair(&r, 0, 1).effective_snps, 3000.0);
    assert_eq!(pair(&r, 0, 1).relationship, ReadRelationship::ThirdDegree);
    assert_eq!(pair(&r, 0, 2).effective_snps, 2999.5);
    assert_eq!(
        pair(&r, 0, 2).relationship,
        ReadRelationship::UnrelatedOrThirdDegree
    );
    // x = t0 exactly is third degree
    assert_eq!(
        classify(0.953125, 3000.0, &READ_THRESHOLDS),
        ReadRelationship::ThirdDegree
    );
    assert_eq!(
        classify(0.90625, 3000.0, &READ_THRESHOLDS),
        ReadRelationship::ThirdDegree
    );
    assert_eq!(
        classify(0.8125, 3000.0, &READ_THRESHOLDS),
        ReadRelationship::SecondDegree
    );
    assert_eq!(
        classify(0.625, 3000.0, &READ_THRESHOLDS),
        ReadRelationship::FirstDegree
    );
}

#[test]
fn e10_pairs_without_a_common_snp_are_left_out_of_output_and_norm() {
    let m = 400;
    let mut strs = random_strs(5, m, 11);
    // samples 3 and 4 never both called
    strs[3].replace_range(..200, &".".repeat(200));
    strs[4].replace_range(200.., &".".repeat(200));
    let refs: Vec<&str> = strs.iter().map(String::as_str).collect();
    let g = from_strs(&refs);
    let pos: Vec<u64> = (0..m as u64).map(|i| 250_000 * i + 1).collect();
    let run_norm = |norm| {
        run(
            &g,
            &pos,
            &ReadOptions {
                norm,
                ..Default::default()
            },
        )
    };
    let r = run_norm(ReadNorm::Median);
    assert_eq!(r.pairs.len(), 9);
    assert!(!r.pairs.iter().any(|p| (p.i, p.j) == (3, 4)));
    let mut p0s: Vec<f64> = r.pairs.iter().map(|p| p.p0).collect();
    let mean_p0 = mean(&p0s);
    p0s.sort_by(f64::total_cmp);
    assert_eq!(r.norm_value, p0s[4]);
    assert_eq!(run_norm(ReadNorm::Mean).norm_value, mean_p0);
    assert_eq!(run_norm(ReadNorm::Max).norm_value, p0s[8]);
    for p in &r.pairs {
        assert_eq!(p.normalized_p0, p.p0 / r.norm_value);
        assert_eq!(p.kinship, 1.0 - p.normalized_p0);
    }
}

#[test]
fn e15_median_of_an_even_count_and_identical_samples() {
    let m = 200;
    let strs = random_strs(4, m, 19);
    let refs: Vec<&str> = strs.iter().map(String::as_str).collect();
    let pos: Vec<u64> = (0..m as u64).map(|i| 250_000 * i + 1).collect();
    let r = run(&from_strs(&refs), &pos, &ReadOptions::default());
    let mut p0s: Vec<f64> = r.pairs.iter().map(|p| p.p0).collect();
    p0s.sort_by(f64::total_cmp);
    assert_eq!(r.norm_value, f64::midpoint(p0s[2], p0s[3]));
    let same = "0220".repeat(50);
    let g = from_strs(&[&same, &same, &same]);
    let err = read_kinship(&g, &pos, &vec![1; m], &ReadOptions::default()).unwrap_err();
    assert_eq!(
        err.to_string(),
        PopmixError::data("READ: the normalization value is 0: are most samples duplicates?")
            .to_string()
    );
}

#[test]
fn e17_too_few_samples_and_no_pair_with_a_common_snp() {
    let pos = [1, 2, 3, 4];
    let one = from_strs(&["0202"]);
    assert!(matches!(
        read_kinship(&one, &pos, &[1; 4], &ReadOptions::default()),
        Err(PopmixError::TooFew {
            required: 2,
            actual: 1,
            ..
        })
    ));
    let apart = from_strs(&["02..", "..20"]);
    let r = read_kinship(
        &apart,
        &pos,
        &[1; 4],
        &ReadOptions {
            norm: ReadNorm::Value(0.5),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(r.norm_value, 0.5);
    assert!(r.pairs.is_empty());
    for norm in [ReadNorm::Median, ReadNorm::Mean, ReadNorm::Max] {
        let opts = ReadOptions {
            norm,
            ..Default::default()
        };
        assert!(matches!(
            read_kinship(&apart, &pos, &[1; 4], &opts),
            Err(PopmixError::InvalidData { .. })
        ));
    }
}

#[test]
fn invalid_input_is_refused_in_order() {
    let g = from_strs(&["0202", "0220"]);
    let pos = [1, 2, 3, 4];
    let ok = ReadOptions::default();
    let dim = |p: &[u64], c: &[u32], what: &str| match read_kinship(&g, p, c, &ok) {
        Err(PopmixError::DimensionMismatch {
            context,
            expected,
            actual,
        }) => {
            assert_eq!(context, format!("READ kinship {what}"));
            assert_eq!(expected, "4");
            assert_eq!(actual, (p.len() + c.len() - 4).to_string());
        }
        other => panic!("{other:?}"),
    };
    dim(&pos[..3], &[1; 4], "positions");
    dim(&pos, &[1; 5], "chroms");
    let param = |opts: ReadOptions| match read_kinship(&g, &pos, &[1; 4], &opts) {
        Err(PopmixError::InvalidParameter { name, details }) => (name, details),
        other => panic!("{other:?}"),
    };
    let (name, details) = param(ReadOptions {
        window_size: 0,
        first_degree_window_size: 0,
        ..Default::default()
    });
    assert_eq!(name, "window_size");
    assert_eq!(details, "must be at least 1 base pair");
    let (name, _) = param(ReadOptions {
        first_degree_window_size: 0,
        ..Default::default()
    });
    assert_eq!(name, "first_degree_window_size");
    for v in [0.0, -0.5, f64::NAN, f64::INFINITY] {
        let (name, details) = param(ReadOptions {
            norm: ReadNorm::Value(v),
            ..Default::default()
        });
        assert_eq!(name, "norm_value");
        assert_eq!(details, format!("must be a positive number, got {v}"));
    }
    for t in [
        [0.9, 0.9, 0.8, 0.6],
        [0.6, 0.8, 0.9, 0.95],
        [f64::NAN, 0.9, 0.8, 0.6],
    ] {
        let (name, _) = param(ReadOptions {
            thresholds: t,
            ..Default::default()
        });
        assert_eq!(name, "thresholds");
    }
}

#[test]
fn defaults_and_labels() {
    let o = ReadOptions::default();
    assert_eq!(o.window_size, 5_000_000);
    assert_eq!(o.first_degree_window_size, 20_000_000);
    assert_eq!(o.norm, ReadNorm::Median);
    assert!(!o.window_estimate);
    assert_eq!(o.thresholds, READ_THRESHOLDS);
    assert_eq!(o.compat, ReadCompat::CORRECTED);
    assert_eq!(ReadCompat::default(), ReadCompat::CORRECTED);
    assert_eq!(
        ReadCompat::V2_01,
        ReadCompat {
            windows_overlap: true,
            jackknife_skips_last_window: true,
            ci_of_normalized_p0: true,
            windows_ignore_chromosomes: true,
        }
    );
    use ReadRelationship::*;
    let labels = [
        (Unrelated, "Unrelated"),
        (
            UnrelatedOrThirdDegree,
            "Unrelated/Consistent with Third Degree",
        ),
        (ThirdDegree, "Third Degree"),
        (SecondDegree, "Second Degree"),
        (FirstDegree, "First Degree"),
        (IdenticalOrTwins, "IdenticalTwins/SameIndividual"),
    ];
    for (r, l) in labels {
        assert_eq!(r.label(), l);
    }
    assert_eq!(FirstDegreeType::ParentOffspring.label(), "Parent-offspring");
    assert_eq!(FirstDegreeType::Siblings.label(), "Siblings");
    assert_eq!(FirstDegreeType::NotApplicable.label(), "NotApplicable");
}
