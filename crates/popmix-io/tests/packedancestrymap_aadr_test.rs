//! Integration tests for PACKEDANCESTRYMAP reader with real AADR data.
//!
//! These tests require the AADR v62 dataset to be present at the expected path.
//! They will be skipped if the data is not available.

use std::path::Path;

use popmix_io::PackedAncestryMapReader;

const AADR_PREFIX: &str = "/references/aadr/v62.0_1240k_public";

fn data_available() -> bool {
    Path::new(&format!("{}.geno", AADR_PREFIX)).exists()
}

/// Test that we can open and read basic metadata from AADR v62 PACKEDANCESTRYMAP dataset.
#[test]
fn test_open_aadr_packedancestrymap_dataset() {
    if !data_available() {
        eprintln!("Skipping test: AADR PACKEDANCESTRYMAP data not available");
        return;
    }
    let reader = PackedAncestryMapReader::new(AADR_PREFIX)
        .expect("Failed to open AADR PACKEDANCESTRYMAP dataset");

    // Verify expected dimensions
    assert_eq!(reader.num_samples(), 17629, "Expected 17629 samples");
    assert_eq!(reader.num_snps(), 1233013, "Expected 1233013 SNPs");

    // Check first sample
    let first_sample = &reader.samples()[0];
    assert!(!first_sample.id.is_empty(), "First sample should have an ID");
    println!("First sample: {} ({})", first_sample.id, first_sample.population);

    // Check first SNP
    let first_snp = &reader.snps()[0];
    assert!(!first_snp.id.is_empty(), "First SNP should have an ID");
    println!("First SNP: {} at chr{} pos {}", first_snp.id, first_snp.chromosome, first_snp.position_bp);
}

/// Test reading a small range of SNPs from AADR PACKEDANCESTRYMAP dataset.
#[test]
fn test_read_snps_aadr_packedancestrymap() {
    if !data_available() {
        eprintln!("Skipping test: AADR PACKEDANCESTRYMAP data not available");
        return;
    }

    let mut reader = PackedAncestryMapReader::new(AADR_PREFIX)
        .expect("Failed to open AADR PACKEDANCESTRYMAP dataset");

    // Read first 100 SNPs
    let geno_matrix = reader.read_snp_range(0..100)
        .expect("Failed to read SNP range");

    assert_eq!(geno_matrix.n_snps(), 100);
    assert_eq!(geno_matrix.n_samples(), 17629);

    // Verify we got some non-missing genotypes
    let mut total_calls = 0u64;
    let mut missing_calls = 0u64;
    let mut geno_counts = [0u64; 4]; // HomRef, Het, HomAlt, Missing

    for snp_idx in 0..geno_matrix.n_snps() {
        for sample_idx in 0..geno_matrix.n_samples() {
            total_calls += 1;
            let geno = geno_matrix.get(snp_idx, sample_idx);
            match geno {
                popmix_core::Genotype::HomRef => geno_counts[0] += 1,
                popmix_core::Genotype::Het => geno_counts[1] += 1,
                popmix_core::Genotype::HomAlt => geno_counts[2] += 1,
                popmix_core::Genotype::Missing => {
                    geno_counts[3] += 1;
                    missing_calls += 1;
                }
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
    println!(
        "Genotype distribution: HomRef={}, Het={}, HomAlt={}, Missing={}",
        geno_counts[0], geno_counts[1], geno_counts[2], geno_counts[3]
    );

    // Missing rate should be reasonable (typically <90% for ancient DNA with coverage issues)
    assert!(
        missing_rate < 0.95,
        "Missing rate {} is too high",
        missing_rate
    );

    // Should have some non-missing genotypes of each type
    assert!(geno_counts[0] > 0, "Should have some HomRef");
    assert!(geno_counts[2] > 0, "Should have some HomAlt");
}

/// Test ploidy detection on AADR PACKEDANCESTRYMAP dataset.
#[test]
fn test_detect_ploidy_aadr_packedancestrymap() {
    if !data_available() {
        eprintln!("Skipping test: AADR PACKEDANCESTRYMAP data not available");
        return;
    }

    let mut reader = PackedAncestryMapReader::new(AADR_PREFIX)
        .expect("Failed to open AADR PACKEDANCESTRYMAP dataset");

    // Detect ploidy based on first 1000 SNPs
    reader.detect_ploidy(1000)
        .expect("Failed to detect ploidy");

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

    // Should detect at least some samples
    assert!(diploid_count > 0 || pseudohaploid_count > 0, "Should detect at least some samples");
}

/// Test reading a single SNP.
#[test]
fn test_read_single_snp_aadr_packedancestrymap() {
    if !data_available() {
        eprintln!("Skipping test: AADR PACKEDANCESTRYMAP data not available");
        return;
    }

    let mut reader = PackedAncestryMapReader::new(AADR_PREFIX)
        .expect("Failed to open AADR PACKEDANCESTRYMAP dataset");

    // Read SNP at index 1000
    let genotypes = reader.read_snp(1000)
        .expect("Failed to read single SNP");

    assert_eq!(genotypes.len(), 17629);

    // Count genotypes
    let mut counts = [0usize; 4];
    for g in &genotypes {
        match g {
            popmix_core::Genotype::HomRef => counts[0] += 1,
            popmix_core::Genotype::Het => counts[1] += 1,
            popmix_core::Genotype::HomAlt => counts[2] += 1,
            popmix_core::Genotype::Missing => counts[3] += 1,
        }
    }

    println!(
        "SNP 1000: HomRef={}, Het={}, HomAlt={}, Missing={}",
        counts[0], counts[1], counts[2], counts[3]
    );
}
