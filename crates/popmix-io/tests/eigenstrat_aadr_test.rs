//! Integration tests for EIGENSTRAT reader with text format data.
//!
//! Note: AADR is distributed as packed binary (PACKEDANCESTRYMAP, TGENO), not text EIGENSTRAT.
//! These tests use a small text EIGENSTRAT fixture.

use popmix_io::EigenstratReader;

const TEST_PREFIX: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/fixtures/small_eigenstrat"
);

/// Test that we can open and read basic metadata from EIGENSTRAT dataset.
#[test]
fn test_open_eigenstrat_dataset() {
    let reader = EigenstratReader::new(TEST_PREFIX).expect("Failed to open EIGENSTRAT dataset");

    // Verify expected dimensions
    assert_eq!(reader.num_samples(), 3, "Expected 3 samples");
    assert_eq!(reader.num_snps(), 5, "Expected 5 SNPs");

    // Check samples
    let samples = reader.samples();
    assert_eq!(samples[0].id, "IND1");
    assert_eq!(samples[0].population.to_string(), "PopA");
    assert_eq!(samples[1].id, "IND2");
    assert_eq!(samples[2].id, "IND3");
    assert_eq!(samples[2].population.to_string(), "PopB");

    // Check SNPs
    let snps = reader.snps();
    assert_eq!(snps[0].id, "rs1");
    assert_eq!(snps[4].id, "rs5");
}

/// Test reading all SNPs from EIGENSTRAT dataset.
#[test]
fn test_read_snps_eigenstrat() {
    let reader = EigenstratReader::new(TEST_PREFIX).expect("Failed to open EIGENSTRAT dataset");

    // Read all SNPs
    let geno_matrix = reader
        .read_snp_range(0..5)
        .expect("Failed to read SNP range");

    assert_eq!(geno_matrix.n_snps(), 5);
    assert_eq!(geno_matrix.n_samples(), 3);

    // Check specific genotypes from our test file:
    // Row 1: 012 -> HomRef, Het, HomAlt
    // Row 2: 120 -> Het, HomAlt, HomRef
    // etc.

    use popmix_core::Genotype;

    // SNP 0 (row "012"): IND1=0, IND2=1, IND3=2
    assert_eq!(geno_matrix.get(0, 0), Genotype::HomRef);
    assert_eq!(geno_matrix.get(0, 1), Genotype::Het);
    assert_eq!(geno_matrix.get(0, 2), Genotype::HomAlt);

    // SNP 1 (row "120"): IND1=1, IND2=2, IND3=0
    assert_eq!(geno_matrix.get(1, 0), Genotype::Het);
    assert_eq!(geno_matrix.get(1, 1), Genotype::HomAlt);
    assert_eq!(geno_matrix.get(1, 2), Genotype::HomRef);

    // SNP 4 (row "222"): All HomAlt
    assert_eq!(geno_matrix.get(4, 0), Genotype::HomAlt);
    assert_eq!(geno_matrix.get(4, 1), Genotype::HomAlt);
    assert_eq!(geno_matrix.get(4, 2), Genotype::HomAlt);
}

/// Test ploidy detection on EIGENSTRAT dataset.
#[test]
fn test_detect_ploidy_eigenstrat() {
    let mut reader = EigenstratReader::new(TEST_PREFIX).expect("Failed to open EIGENSTRAT dataset");

    // Detect ploidy based on all SNPs
    reader.detect_ploidy(5).expect("Failed to detect ploidy");

    // All samples should have ploidy set
    for sample in reader.samples() {
        // Since we have heterozygous calls (1), samples should be diploid
        assert_eq!(sample.ploidy, popmix_core::Ploidy::Diploid);
    }
}
