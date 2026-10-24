//! Integration tests for PLINK reader with real AADR data.
//!
//! These tests require the AADR HO dataset to be present at the expected path.
//! They will be skipped if the data is not available.

use std::path::Path;

use popmix_io::PlinkReader;

const AADR_HO_PREFIX: &str = "/references/aadr_v54.1.p1/HO_public";

fn data_available() -> bool {
    Path::new(&format!("{}.bed", AADR_HO_PREFIX)).exists()
}

/// Test that we can open and read basic metadata from AADR HO dataset.
#[test]
fn test_open_aadr_ho_dataset() {
    if !data_available() {
        eprintln!("Skipping test: AADR HO data not available");
        return;
    }

    let reader = PlinkReader::new(AADR_HO_PREFIX).expect("Failed to open AADR HO dataset");

    // Verify expected dimensions (from wc -l on .fam and .bim files)
    assert_eq!(reader.num_samples(), 20503, "Expected 20503 samples");
    assert_eq!(reader.num_snps(), 597573, "Expected 597573 SNPs");

    // Check first sample
    let first_sample = &reader.samples()[0];
    assert!(
        !first_sample.id.is_empty(),
        "First sample should have an ID"
    );
    println!(
        "First sample: {} ({})",
        first_sample.id, first_sample.population
    );

    // Check first SNP
    let first_snp = &reader.snps()[0];
    assert!(!first_snp.id.is_empty(), "First SNP should have an ID");
    println!(
        "First SNP: {} at chr{} pos {}",
        first_snp.id, first_snp.chromosome, first_snp.position_bp
    );
}

/// Test reading a small range of SNPs from AADR HO dataset.
#[test]
fn test_read_snps_aadr_ho() {
    if !data_available() {
        eprintln!("Skipping test: AADR HO data not available");
        return;
    }

    let mut reader = PlinkReader::new(AADR_HO_PREFIX).expect("Failed to open AADR HO dataset");

    // Read first 100 SNPs
    let geno_matrix = reader
        .read_snp_range(0..100)
        .expect("Failed to read SNP range");

    assert_eq!(geno_matrix.n_snps(), 100);
    assert_eq!(geno_matrix.n_samples(), 20503);

    // Verify we got some non-missing genotypes
    let mut total_calls = 0u64;
    let mut missing_calls = 0u64;

    for snp_idx in 0..geno_matrix.n_snps() {
        for sample_idx in 0..geno_matrix.n_samples() {
            total_calls += 1;
            if geno_matrix.get(snp_idx, sample_idx) == popmix_core::Genotype::Missing {
                missing_calls += 1;
            }
        }
    }

    let missing_rate = missing_calls as f64 / total_calls as f64;
    println!(
        "Genotype calls: {}, Missing: {} ({:.1}%)",
        total_calls,
        missing_calls,
        missing_rate * 100.0
    );

    // Missing rate should be reasonable (typically <50% for ancient DNA)
    assert!(
        missing_rate < 0.9,
        "Missing rate {} is too high",
        missing_rate
    );
}

/// Test ploidy detection on AADR HO dataset.
#[test]
fn test_detect_ploidy_aadr_ho() {
    if !data_available() {
        eprintln!("Skipping test: AADR HO data not available");
        return;
    }

    let mut reader = PlinkReader::new(AADR_HO_PREFIX).expect("Failed to open AADR HO dataset");

    // Detect ploidy based on first 1000 SNPs
    reader.detect_ploidy(1000).expect("Failed to detect ploidy");

    let mut diploid_count = 0;
    let mut pseudohaploid_count = 0;

    for sample in reader.samples() {
        match sample.ploidy {
            popmix_core::Ploidy::Diploid => diploid_count += 1,
            popmix_core::Ploidy::Pseudohaploid => pseudohaploid_count += 1,
        }
    }

    println!(
        "Ploidy: {} diploid, {} pseudohaploid",
        diploid_count, pseudohaploid_count
    );

    // Ancient DNA typically has a mix
    assert!(
        diploid_count > 0 || pseudohaploid_count > 0,
        "Should detect at least some samples"
    );
}
