
use super::*;
use popmix_core::genotypes::{Chromosome, Ploidy, PopulationId};

fn make_snp(chr: u8, pos: u32, ref_allele: char, alt_allele: char) -> Snp {
    Snp {
        id: format!("snp_{}_{}", chr, pos),
        chromosome: Chromosome::Autosome(chr),
        position_morgans: pos as f64 / 1_000_000.0,
        position_bp: pos,
        allele_ref: ref_allele,
        allele_alt: alt_allele,
    }
}

fn make_sample(id: &str, pop: &str) -> Sample {
    Sample {
        id: id.to_string(),
        population: PopulationId(pop.to_string()),
        ploidy: Ploidy::Diploid,
        sex: Default::default(),
    }
}

// =========================================================================
// TDD RED PHASE: Tests that should initially fail
// =========================================================================

#[test]
fn test_merge_snp_intersection() {
    // File A: SNPs at pos 100, 200, 300 on chr 1
    let snps_a = vec![
        make_snp(1, 100, 'A', 'G'),
        make_snp(1, 200, 'C', 'T'),
        make_snp(1, 300, 'G', 'A'),
    ];

    // File B: SNPs at pos 200, 300, 400 on chr 1
    let snps_b = vec![
        make_snp(1, 200, 'C', 'T'),
        make_snp(1, 300, 'G', 'A'),
        make_snp(1, 400, 'T', 'C'),
    ];

    let (merged, stats) = build_snp_intersection(&snps_a, &snps_b);

    // Should only keep SNPs at pos 200, 300 (intersection)
    assert_eq!(merged.len(), 2);
    assert_eq!(merged[0].snp.position_bp, 200);
    assert_eq!(merged[1].snp.position_bp, 300);

    // Stats should show 1 SNP only in A, 1 only in B
    assert_eq!(stats.snps_only_a, 1); // pos 100
    assert_eq!(stats.snps_only_b, 1); // pos 400
    assert_eq!(stats.snps_identical, 2);
}

#[test]
fn test_merge_matching_alleles() {
    // Both files have A/G at position 100 on chr 1
    let snp_a = make_snp(1, 100, 'A', 'G');
    let snp_b = make_snp(1, 100, 'A', 'G');

    let result = match_alleles(&snp_a, &snp_b);
    assert_eq!(result, AlleleMatch::Identical);
}

#[test]
fn test_merge_flipped_alleles() {
    // File A: REF=A, ALT=G
    // File B: REF=G, ALT=A (flipped)
    let snp_a = make_snp(1, 100, 'A', 'G');
    let snp_b = make_snp(1, 100, 'G', 'A');

    let result = match_alleles(&snp_a, &snp_b);
    assert_eq!(result, AlleleMatch::Flipped);
}

#[test]
fn test_merge_rejects_strand_ambiguous() {
    // A/T SNP - cannot reliably resolve strand
    let snp_a = make_snp(1, 100, 'A', 'T');
    let snp_b = make_snp(1, 100, 'A', 'T');

    let result = match_alleles(&snp_a, &snp_b);
    assert_eq!(result, AlleleMatch::StrandAmbiguous);

    // C/G SNP - also strand ambiguous
    let snp_c = make_snp(1, 200, 'C', 'G');
    let snp_d = make_snp(1, 200, 'C', 'G');

    let result2 = match_alleles(&snp_c, &snp_d);
    assert_eq!(result2, AlleleMatch::StrandAmbiguous);
}

#[test]
fn test_merge_combines_samples() {
    // File A: 3 samples
    let samples_a = vec![
        make_sample("A1", "PopA"),
        make_sample("A2", "PopA"),
        make_sample("A3", "PopA"),
    ];

    // File B: 5 samples (different IDs)
    let samples_b = vec![
        make_sample("B1", "PopB"),
        make_sample("B2", "PopB"),
        make_sample("B3", "PopB"),
        make_sample("B4", "PopB"),
        make_sample("B5", "PopB"),
    ];

    // Check no duplicates
    let duplicates = find_duplicate_samples(&samples_a, &samples_b);
    assert!(duplicates.is_empty());

    // Combined should have 8 samples
    let mut combined = samples_a.clone();
    combined.extend(samples_b.iter().cloned());
    assert_eq!(combined.len(), 8);
}

#[test]
fn test_merge_duplicate_sample_error() {
    // Same sample ID "shared" in both files
    let samples_a = vec![
        make_sample("unique_a", "PopA"),
        make_sample("shared", "PopA"),
    ];

    let samples_b = vec![
        make_sample("unique_b", "PopB"),
        make_sample("shared", "PopB"),
    ];

    let duplicates = find_duplicate_samples(&samples_a, &samples_b);
    assert_eq!(duplicates.len(), 1);
    assert_eq!(duplicates[0], "shared");
}

#[test]
fn test_recode_flipped_genotypes() {
    // HomRef (0) should become HomAlt (2)
    assert_eq!(recode_flipped(Genotype::HomRef), Genotype::HomAlt);

    // HomAlt (2) should become HomRef (0)
    assert_eq!(recode_flipped(Genotype::HomAlt), Genotype::HomRef);

    // Het (1) stays Het (1)
    assert_eq!(recode_flipped(Genotype::Het), Genotype::Het);

    // Missing stays Missing
    assert_eq!(recode_flipped(Genotype::Missing), Genotype::Missing);
}

#[test]
fn test_allele_mismatch() {
    // A third allele (A/G vs A/C)
    let snp_a = make_snp(1, 100, 'A', 'G');
    let snp_b = make_snp(1, 100, 'A', 'C');

    let result = match_alleles(&snp_a, &snp_b);
    assert_eq!(result, AlleleMatch::Mismatch);
}

#[test]
fn alleles_on_the_other_strand_are_matched() {
    let a = make_snp(1, 100, 'A', 'G');
    let same = reconcile_alleles(&a, &make_snp(1, 100, 'T', 'C'));
    assert_eq!(
        (same.allele_match, same.other_strand),
        (AlleleMatch::Identical, true)
    );
    let swapped = reconcile_alleles(&a, &make_snp(1, 100, 'C', 'T'));
    assert_eq!(
        (swapped.allele_match, swapped.other_strand),
        (AlleleMatch::Flipped, true)
    );
    assert_eq!(swapped.alleles, ('A', 'G'));
}

#[test]
fn an_unknown_allele_is_taken_from_the_other_file() {
    // PLINK writes 0, EIGENSTRAT X, VCF . for the allele of a monomorphic SNP
    let b = reconcile_alleles(&make_snp(1, 100, 'A', 'G'), &make_snp(1, 100, 'G', '0'));
    assert_eq!(
        (b.allele_match, b.alleles),
        (AlleleMatch::Flipped, ('A', 'G'))
    );
    let a = reconcile_alleles(&make_snp(1, 100, 'C', 'X'), &make_snp(1, 100, 'T', 'C'));
    assert_eq!(
        (a.allele_match, a.alleles),
        (AlleleMatch::Flipped, ('C', 'T'))
    );
    let both = reconcile_alleles(&make_snp(1, 100, 'C', '.'), &make_snp(1, 100, 'C', 'X'));
    assert_eq!(both.allele_match, AlleleMatch::Identical);
    // C on one side and G on the other is a C/G SNP
    let cg = reconcile_alleles(&make_snp(1, 100, 'C', 'X'), &make_snp(1, 100, 'G', '0'));
    assert_eq!(cg.allele_match, AlleleMatch::StrandAmbiguous);
}

#[test]
fn snps_sharing_a_position_are_matched_by_id() {
    let named = |id: &str, pos, r, a| Snp {
        id: id.to_string(),
        ..make_snp(1, pos, r, a)
    };
    let snps_a = vec![named("rs7", 700, 'T', 'C'), named("rs9", 900, 'A', 'G')];
    let snps_b = vec![
        named("rs7x", 700, 'A', 'G'),
        named("rs7", 700, 'C', 'T'),
        named("rs9a", 900, 'A', 'G'),
        named("rs9b", 900, 'A', 'G'),
    ];
    let (merged, stats) = build_snp_intersection(&snps_a, &snps_b);
    assert_eq!(merged.len(), 1);
    assert_eq!(
        (merged[0].file_b_idx, merged[0].allele_match),
        (1, AlleleMatch::Flipped)
    );
    // two candidates for rs9 and neither has its name: left out, counted
    assert_eq!(stats.snps_duplicate_position, 1);
    assert_eq!(stats.snps_only_b, 1);
}

#[test]
fn test_is_strand_ambiguous() {
    // A/T and T/A are strand ambiguous
    assert!(is_strand_ambiguous('A', 'T'));
    assert!(is_strand_ambiguous('T', 'A'));

    // C/G and G/C are strand ambiguous
    assert!(is_strand_ambiguous('C', 'G'));
    assert!(is_strand_ambiguous('G', 'C'));

    // Other combinations are NOT strand ambiguous
    assert!(!is_strand_ambiguous('A', 'G'));
    assert!(!is_strand_ambiguous('A', 'C'));
    assert!(!is_strand_ambiguous('G', 'T'));
    assert!(!is_strand_ambiguous('C', 'T'));
}

#[test]
fn test_merge_genotypes_full() {
    // Create two small genotype matrices and merge them

    // File A: 2 samples, 3 SNPs
    let samples_a = vec![make_sample("A1", "PopA"), make_sample("A2", "PopA")];
    let snps_a = vec![
        make_snp(1, 100, 'A', 'G'), // Only in A
        make_snp(1, 200, 'C', 'T'), // In both, identical
        make_snp(1, 300, 'G', 'A'), // In both, flipped
    ];
    // Genotypes: [n_snps, n_samples]
    // SNP1: A1=HomRef(0), A2=Het(1)
    // SNP2: A1=Het(1), A2=HomAlt(2)
    // SNP3: A1=HomAlt(2), A2=Missing(3)
    let data_a = ndarray::array![
        [0u8, 1], // SNP at 100
        [1u8, 2], // SNP at 200
        [2u8, 3], // SNP at 300
    ];
    let matrix_a = GenotypeMatrix::new(data_a, samples_a, snps_a).unwrap();

    // File B: 3 samples, 3 SNPs
    let samples_b = vec![
        make_sample("B1", "PopB"),
        make_sample("B2", "PopB"),
        make_sample("B3", "PopB"),
    ];
    let snps_b = vec![
        make_snp(1, 200, 'C', 'T'), // In both, identical
        make_snp(1, 300, 'A', 'G'), // In both, flipped (G/A -> A/G)
        make_snp(1, 400, 'T', 'C'), // Only in B
    ];
    // Genotypes:
    // SNP1 (200): B1=HomRef(0), B2=Het(1), B3=HomAlt(2)
    // SNP2 (300): B1=HomRef(0), B2=Het(1), B3=HomAlt(2) <- will be recoded
    // SNP3 (400): B1=Het(1), B2=HomRef(0), B3=Missing(3)
    let data_b = ndarray::array![
        [0u8, 1, 2], // SNP at 200
        [0u8, 1, 2], // SNP at 300 (needs recoding: 0->2, 2->0)
        [1u8, 0, 3], // SNP at 400
    ];
    let matrix_b = GenotypeMatrix::new(data_b, samples_b, snps_b).unwrap();

    // Merge
    let (merged, stats) = merge_genotypes(&matrix_a, &matrix_b).unwrap();

    // Check stats
    assert_eq!(stats.snps_only_a, 1); // SNP at 100
    assert_eq!(stats.snps_only_b, 1); // SNP at 400
    assert_eq!(stats.snps_identical, 1); // SNP at 200
    assert_eq!(stats.snps_flipped, 1); // SNP at 300
    assert_eq!(stats.samples_total, 5); // 2 + 3

    // Check merged matrix dimensions
    assert_eq!(merged.n_snps(), 2); // Only intersection (200, 300)
    assert_eq!(merged.n_samples(), 5); // A1, A2, B1, B2, B3

    // Check SNP positions
    assert_eq!(merged.snps()[0].position_bp, 200);
    assert_eq!(merged.snps()[1].position_bp, 300);

    // Check samples
    assert_eq!(merged.samples()[0].id, "A1");
    assert_eq!(merged.samples()[1].id, "A2");
    assert_eq!(merged.samples()[2].id, "B1");
    assert_eq!(merged.samples()[3].id, "B2");
    assert_eq!(merged.samples()[4].id, "B3");

    // Check genotypes at SNP 200 (identical alleles, no recoding)
    // Original A: [1, 2], Original B: [0, 1, 2]
    assert_eq!(merged.get(0, 0), Genotype::Het); // A1
    assert_eq!(merged.get(0, 1), Genotype::HomAlt); // A2
    assert_eq!(merged.get(0, 2), Genotype::HomRef); // B1
    assert_eq!(merged.get(0, 3), Genotype::Het); // B2
    assert_eq!(merged.get(0, 4), Genotype::HomAlt); // B3

    // Check genotypes at SNP 300 (flipped alleles, B needs recoding)
    // Original A: [2, 3] (HomAlt, Missing)
    // Original B: [0, 1, 2] -> recoded to [2, 1, 0] (HomAlt, Het, HomRef)
    assert_eq!(merged.get(1, 0), Genotype::HomAlt); // A1
    assert_eq!(merged.get(1, 1), Genotype::Missing); // A2
    assert_eq!(merged.get(1, 2), Genotype::HomAlt); // B1 (was 0, recoded to 2)
    assert_eq!(merged.get(1, 3), Genotype::Het); // B2 (was 1, stays 1)
    assert_eq!(merged.get(1, 4), Genotype::HomRef); // B3 (was 2, recoded to 0)
}

#[test]
fn merged_rows_read_in_place_equal_the_merged_matrix() {
    use crate::writers::{write_packedancestrymap, write_packedped};

    // A: 5 samples (a record of 2 bytes, padded), B: 3 samples, B's SNPs
    // in another order, one swapped, one only in each file
    let samples_a: Vec<Sample> = (0..5)
        .map(|i| make_sample(&format!("A{i}"), "PopA"))
        .collect();
    let snps_a = vec![
        make_snp(1, 100, 'A', 'G'),
        make_snp(1, 200, 'C', 'T'),
        make_snp(1, 300, 'G', 'A'),
        make_snp(2, 50, 'A', 'C'),
    ];
    let data_a = ndarray::array![
        [0u8, 1, 2, 3, 0],
        [1, 2, 3, 0, 1],
        [2, 3, 0, 1, 2],
        [3, 0, 1, 2, 3],
    ];
    let samples_b: Vec<Sample> = (0..3)
        .map(|i| make_sample(&format!("B{i}"), "PopB"))
        .collect();
    let snps_b = vec![
        make_snp(2, 50, 'A', 'C'),
        make_snp(1, 300, 'A', 'G'),
        make_snp(1, 400, 'T', 'C'),
        make_snp(1, 200, 'C', 'T'),
    ];
    let data_b = ndarray::array![[0u8, 1, 2], [0, 2, 3], [1, 0, 3], [2, 1, 0]];
    let a = GenotypeMatrix::new(data_a, samples_a, snps_a).unwrap();
    let b = GenotypeMatrix::new(data_b, samples_b, snps_b).unwrap();
    let (merged, _) = merge_genotypes(&a, &b).unwrap();

    let dir = tempfile::tempdir().unwrap();
    let (pa, pb) = (dir.path().join("a"), dir.path().join("b"));
    write_packedancestrymap(&pa, &a).unwrap();
    write_packedped(&pb, &b).unwrap();
    let fa = SnpMajor::open_packed(&format!("{}.geno", pa.display()), 5, 4).unwrap();
    let fb = SnpMajor::open_bed(&format!("{}.bed", pb.display()), 3, 4).unwrap();
    let (snps, _) = build_snp_intersection(a.snps(), b.snps());
    let mut rows = MergedRows::new(&fa, &fb, &snps);
    let mut row = vec![9u8; 8];
    for s in 0..merged.n_snps() {
        rows.next_row(&mut row).unwrap();
        let expected: Vec<u8> = (0..8).map(|i| merged.get(s, i) as u8).collect();
        assert_eq!(row, expected, "SNP {s}");
    }
    assert!(rows.next_row(&mut row).is_err());

    // and a block of samples at a time, across the two files
    use crate::writers::SampleGenotypes;
    let n = merged.n_snps();
    let mut block = vec![9u8; 4 * n];
    rows.fill_samples(3..7, n, &mut block);
    for (k, s) in (3..7).enumerate() {
        let expected: Vec<u8> = (0..n).map(|j| merged.get(j, s) as u8).collect();
        assert_eq!(&block[k * n..(k + 1) * n], &expected[..], "sample {s}");
    }
}

#[test]
fn test_merge_duplicate_sample_error_full() {
    // Create matrices with duplicate sample IDs
    let samples_a = vec![make_sample("shared", "PopA")];
    let snps_a = vec![make_snp(1, 100, 'A', 'G')];
    let data_a = ndarray::array![[0u8]];
    let matrix_a = GenotypeMatrix::new(data_a, samples_a, snps_a).unwrap();

    let samples_b = vec![make_sample("shared", "PopB")];
    let snps_b = vec![make_snp(1, 100, 'A', 'G')];
    let data_b = ndarray::array![[1u8]];
    let matrix_b = GenotypeMatrix::new(data_b, samples_b, snps_b).unwrap();

    // Should error due to duplicate sample
    let result = merge_genotypes(&matrix_a, &matrix_b);
    assert!(result.is_err());
    assert!(result
        .unwrap_err()
        .to_string()
        .contains("Duplicate sample IDs"));
}

#[test]
fn test_merge_no_overlap_error() {
    // Create matrices with no overlapping SNPs
    let samples_a = vec![make_sample("A1", "PopA")];
    let snps_a = vec![make_snp(1, 100, 'A', 'G')];
    let data_a = ndarray::array![[0u8]];
    let matrix_a = GenotypeMatrix::new(data_a, samples_a, snps_a).unwrap();

    let samples_b = vec![make_sample("B1", "PopB")];
    let snps_b = vec![make_snp(1, 200, 'C', 'T')]; // Different position
    let data_b = ndarray::array![[1u8]];
    let matrix_b = GenotypeMatrix::new(data_b, samples_b, snps_b).unwrap();

    // Should error due to no overlapping SNPs
    let result = merge_genotypes(&matrix_a, &matrix_b);
    assert!(result.is_err());
    assert!(result
        .unwrap_err()
        .to_string()
        .contains("No overlapping SNPs"));
}
