
use super::*;
use crate::genotypes::{PopulationId, Sample, Snp};
use ndarray::Array2;

fn create_test_data() -> (GenotypeMatrix, AlleleFrequencies, Vec<f64>) {
    create_test_data_with_spacing(0.005)
}

fn create_test_data_with_spacing(
    spacing_morgans: f64,
) -> (GenotypeMatrix, AlleleFrequencies, Vec<f64>) {
    let mut samples = Vec::new();
    for i in 0..10 {
        samples.push(Sample::new(
            format!("test_{}", i),
            PopulationId::new("test"),
        ));
    }
    for i in 0..5 {
        samples.push(Sample::new(
            format!("src1_{}", i),
            PopulationId::new("source1"),
        ));
    }
    for i in 0..5 {
        samples.push(Sample::new(
            format!("src2_{}", i),
            PopulationId::new("source2"),
        ));
    }

    let snps: Vec<Snp> = (0..50)
        .map(|i| {
            Snp::new(
                format!("rs{}", i),
                Chromosome::Autosome(1),
                (i * 10000) as u32,
                i as f64 * spacing_morgans,
                'A',
                'G',
            )
        })
        .collect();

    let n_samples = samples.len();
    let n_snps = snps.len();
    let mut data = Array2::zeros((n_snps, n_samples));

    for snp_idx in 0..n_snps {
        for sample_idx in 0..n_samples {
            let base = ((snp_idx + sample_idx) % 3) as u8;
            data[[snp_idx, sample_idx]] = base;
        }
    }

    let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
    let af = AlleleFrequencies::from_genotypes(&geno);

    let idx1 = af.pop_index("source1").unwrap();
    let idx2 = af.pop_index("source2").unwrap();
    let weights: Vec<f64> = (0..n_snps)
        .map(|snp_idx| af.frequencies[[snp_idx, idx1]] - af.frequencies[[snp_idx, idx2]])
        .collect();

    (geno, af, weights)
}

fn source_freqs(geno: &GenotypeMatrix, af: &AlleleFrequencies) -> (Vec<f64>, Vec<f64>) {
    let idx1 = af.pop_index("source1").unwrap();
    let idx2 = af.pop_index("source2").unwrap();
    let n_snps = geno.n_snps();
    (
        (0..n_snps).map(|i| af.frequencies[[i, idx1]]).collect(),
        (0..n_snps).map(|i| af.frequencies[[i, idx2]]).collect(),
    )
}

#[test]
fn bins_do_not_depend_on_the_thread_count() {
    // 64 test samples with irregular genotypes, so the order in which
    // their sums are added shows in the last bits.
    let mut samples: Vec<Sample> = (0..64)
        .map(|i| Sample::new(format!("test_{i}"), PopulationId::new("test")))
        .collect();
    for (pop, n) in [("source1", 8), ("source2", 8)] {
        samples.extend((0..n).map(|i| Sample::new(format!("{pop}_{i}"), PopulationId::new(pop))));
    }
    let snps: Vec<Snp> = (0..200)
        .map(|i| {
            Snp::new(
                format!("rs{i}"),
                Chromosome::Autosome(1 + (i / 100) as u8),
                (i * 7919) as u32,
                (i % 100) as f64 * 0.0031,
                'A',
                'G',
            )
        })
        .collect();
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let data = Array2::from_shape_fn((snps.len(), samples.len()), |_| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % 3) as u8
    });
    let geno = GenotypeMatrix::new(data, samples, snps).unwrap();
    let af = AlleleFrequencies::from_genotypes(&geno);
    let (s1, s2) = source_freqs(&geno, &af);
    let weights: Vec<f64> = s1.iter().zip(&s2).map(|(a, b)| a - b).collect();
    for qbin in [0, 10] {
        let run = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let opts = LdOptions {
                bin_width_cm: 0.05,
                max_distance_cm: 2.0,
                min_pairs_per_bin: 1,
                qbin,
                ..Default::default()
            };
            let bins = pool
                .install(|| compute_ld_decay(&geno, &af, "test", &weights, &s1, &s2, Some(opts)))
                .unwrap();
            bins.iter()
                .map(|b| (b.weighted_ld.to_bits(), b.correlation.to_bits()))
                .collect::<Vec<_>>()
        };
        let one = run(1);
        assert!(!one.is_empty());
        for threads in [2, 3, 8] {
            assert_eq!(run(threads), one, "qbin {qbin}, {threads} threads");
        }
    }
}

#[test]
fn test_ld_decay_basic() {
    let (geno, af, weights) = create_test_data();
    let (source1_freqs, source2_freqs) = source_freqs(&geno, &af);

    let opts = LdOptions {
        bin_width_cm: 1.0,
        max_distance_cm: 25.0,
        min_pairs_per_bin: 1,
        qbin: 0,
        ..Default::default()
    };

    let result = compute_ld_decay(
        &geno,
        &af,
        "test",
        &weights,
        &source1_freqs,
        &source2_freqs,
        Some(opts),
    )
    .unwrap();

    assert!(!result.is_empty(), "Should have LD decay bins");

    for bin in &result {
        assert!(bin.distance_cm > 0.0);
        assert!(bin.n_pairs > 0);
    }
}

#[test]
fn test_mesh_matches_exact_when_snps_sit_on_mesh_cells() {
    // SNPs exactly one bin (1/128 Morgan) apart with qbin = 2 fall in
    // distinct mesh cells at lags that are exact multiples of the cell
    // size, so both DATES modes see the same pairs in the same bins.
    let (geno, af, weights) = create_test_data_with_spacing(1.0 / 128.0);
    let (source1_freqs, source2_freqs) = source_freqs(&geno, &af);
    let run = |qbin| {
        let opts = LdOptions {
            bin_width_cm: 100.0 / 128.0,
            max_distance_cm: 10.0 * 100.0 / 128.0,
            min_pairs_per_bin: 0,
            qbin,
            ..Default::default()
        };
        compute_ld_decay(
            &geno,
            &af,
            "test",
            &weights,
            &source1_freqs,
            &source2_freqs,
            Some(opts),
        )
        .unwrap()
    };
    let exact = run(0);
    let mesh = run(2);
    assert_eq!(exact.len(), mesh.len());
    for (e, m) in exact.iter().zip(mesh.iter()) {
        assert_eq!(e.n_pairs, m.n_pairs, "pairs at {} cM", e.distance_cm);
        assert!((e.weighted_ld - m.weighted_ld).abs() < 1e-12);
        assert!((e.correlation - m.correlation).abs() < 1e-9);
    }
}

// Mesh geometry for create_test_data_with_spacing(1/128) with qbin = 1:
// 50 SNPs in consecutive mesh cells span a window of m = 50 cells, so
// DATES's FFT length is n = 2^(ceil(log2 50) + 1) = 128. Real pairs have
// lags of at most 49 cells. A circular lag d > n - m = 78 also picks up
// the pairs at negative lag n - d, from bin 79 on. DATES v753 leaves lag
// n - 1 = 127 unscaled (DATES-0002), and it holds the wrapped lag-1 pairs.
const MESH_TEST_BINS: usize = 130;
const FIRST_WRAPPED_BIN: usize = 79;
const UNSCALED_BIN: usize = 127;

fn run_mesh_test(qbin: usize, compat: DatesCompat) -> Vec<LdBin> {
    let (geno, af, weights) = create_test_data_with_spacing(1.0 / 128.0);
    let (source1_freqs, source2_freqs) = source_freqs(&geno, &af);
    let opts = LdOptions {
        bin_width_cm: 100.0 / 128.0,
        max_distance_cm: MESH_TEST_BINS as f64 * 100.0 / 128.0,
        min_pairs_per_bin: 0,
        qbin,
        compat,
    };
    compute_ld_decay(
        &geno,
        &af,
        "test",
        &weights,
        &source1_freqs,
        &source2_freqs,
        Some(opts),
    )
    .unwrap()
}

fn same_bin(a: &LdBin, b: &LdBin) -> bool {
    a.n_pairs == b.n_pairs
        && (a.weighted_ld - b.weighted_ld).abs() < 1e-12
        && (a.correlation - b.correlation).abs() < 1e-9
}

/// Indices of the bins where two LD tables differ.
fn differing_bins(a: &[LdBin], b: &[LdBin]) -> Vec<usize> {
    assert_eq!(a.len(), b.len());
    (0..a.len()).filter(|&k| !same_bin(&a[k], &b[k])).collect()
}

#[test]
fn test_corrected_mesh_matches_exact_past_the_fft_window() {
    // The oracle: on this grid the mesh sees exactly the pairs the exact
    // mode sees, at every distance, so the corrected mesh must agree with
    // it everywhere, including lags DATES wraps.
    let exact = run_mesh_test(0, DatesCompat::CORRECTED);
    let mesh = run_mesh_test(1, DatesCompat::CORRECTED);
    assert_eq!(differing_bins(&exact, &mesh), Vec::<usize>::new());
}

#[test]
fn test_v753_mesh_reproduces_dates_quirks() {
    let exact = run_mesh_test(0, DatesCompat::CORRECTED);
    let mesh = run_mesh_test(1, DatesCompat::V753);
    let diff = differing_bins(&exact, &mesh);
    assert_eq!(diff.first(), Some(&FIRST_WRAPPED_BIN));
    assert!(diff.contains(&UNSCALED_BIN));
}

#[test]
fn test_circular_fft_lags_flag_only_changes_wrapped_bins() {
    let only = DatesCompat {
        circular_fft_lags: true,
        ..DatesCompat::CORRECTED
    };
    let fixed = run_mesh_test(1, DatesCompat::CORRECTED);
    let quirk = run_mesh_test(1, only);
    let diff = differing_bins(&fixed, &quirk);
    assert!(!diff.is_empty());
    assert!(diff.iter().all(|&k| k >= FIRST_WRAPPED_BIN), "{diff:?}");
}

#[test]
fn test_unscaled_last_lag_flag_only_changes_its_bin() {
    // Lag n - 1 has no real pairs, so the missing scaling only shows when
    // the circular quirk puts wrapped pairs there.
    let circular = DatesCompat {
        circular_fft_lags: true,
        ..DatesCompat::CORRECTED
    };
    let both = DatesCompat {
        unscaled_last_fft_lag: true,
        ..circular
    };
    let diff = differing_bins(&run_mesh_test(1, circular), &run_mesh_test(1, both));
    assert_eq!(diff, vec![UNSCALED_BIN]);
}

#[test]
fn test_leave_one_out_with_single_chromosome_is_empty() {
    let (geno, af, weights) = create_test_data();
    let (source1_freqs, source2_freqs) = source_freqs(&geno, &af);
    let ld = compute_ld_decay_by_chromosome(
        &geno,
        &af,
        "test",
        &weights,
        &source1_freqs,
        &source2_freqs,
        None,
    )
    .unwrap();
    assert_eq!(ld.chromosomes.len(), 1);
    assert_eq!(ld.chromosomes[0].chromosome, 1);
    assert_eq!(ld.chromosomes[0].n_snps, 50);
    assert!(ld.chromosomes[0]
        .leave_one_out
        .iter()
        .all(|b| b.n_pairs == 0));
}

#[test]
fn test_ld_decay_missing_population() {
    let (geno, af, weights) = create_test_data();

    let n_snps = geno.n_snps();
    let source1_freqs: Vec<f64> = vec![0.5; n_snps];
    let source2_freqs: Vec<f64> = vec![0.5; n_snps];

    let result = compute_ld_decay(
        &geno,
        &af,
        "nonexistent",
        &weights,
        &source1_freqs,
        &source2_freqs,
        None,
    );
    assert!(result.is_err());
}

#[test]
fn a_snp_without_a_genetic_position_is_an_error() {
    let (geno, af, weights) = create_test_data();
    let (s1, s2) = source_freqs(&geno, &af);
    let mut snps = geno.snps().to_vec();
    snps[7].position_morgans = f64::NAN;
    let data = Array2::from_shape_fn((geno.n_snps(), geno.n_samples()), |(i, j)| {
        geno.get(i, j) as u8
    });
    let geno = GenotypeMatrix::new(data, geno.samples().to_vec(), snps).unwrap();
    let err = compute_ld_decay(&geno, &af, "test", &weights, &s1, &s2, None).unwrap_err();
    assert!(err.to_string().contains("rs7"), "{err}");
}
