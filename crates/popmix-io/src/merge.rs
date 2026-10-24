//! Genotype file merge algorithm.
//!
//! Merges two genotype datasets by SNP intersection, handling allele matching
//! and genotype recoding when alleles are flipped.
//!
//! Based on PLINK 2.0 merge implementation by Shaun Purcell & Christopher Chang.
//! Reference: https://github.com/chrchang/plink-ng
//! License: GPLv3

use std::collections::{HashMap, HashSet};

use popmix_core::error::{PopmixError, Result};
use popmix_core::genotype_matrix::GenotypeMatrix;
use popmix_core::genotypes::{Chromosome, Genotype, Sample, Snp};

/// Key for matching SNPs between files (chromosome + position)
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct SnpKey {
    pub chromosome: Chromosome,
    pub position: u32,
}

impl SnpKey {
    pub fn from_snp(snp: &Snp) -> Self {
        Self {
            chromosome: snp.chromosome.clone(),
            position: snp.position_bp,
        }
    }
}

/// Result of comparing alleles between two SNPs
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleMatch {
    /// REF/ALT match exactly
    Identical,
    /// REF/ALT are swapped (need genotype recode: 0↔2)
    Flipped,
    /// A/T or C/G - cannot reliably resolve strand
    StrandAmbiguous,
    /// Different alleles entirely
    Mismatch,
}

/// Information about a merged SNP
#[derive(Debug, Clone)]
pub struct MergedSnp {
    /// The SNP metadata (from file A)
    pub snp: Snp,
    /// Index in file A's SNP list
    pub file_a_idx: usize,
    /// Index in file B's SNP list
    pub file_b_idx: usize,
    /// How alleles match between files
    pub allele_match: AlleleMatch,
}

/// Statistics from a merge operation
#[derive(Debug, Clone, Default)]
pub struct MergeStats {
    /// SNPs only in file A
    pub snps_only_a: usize,
    /// SNPs only in file B
    pub snps_only_b: usize,
    /// SNPs in both files with identical alleles
    pub snps_identical: usize,
    /// SNPs in both files with flipped alleles
    pub snps_flipped: usize,
    /// SNPs skipped due to strand ambiguity
    pub snps_strand_ambiguous: usize,
    /// SNPs with mismatched alleles
    pub snps_mismatch: usize,
    /// Total samples in merged output
    pub samples_total: usize,
}

/// Check if a SNP has strand-ambiguous alleles (A/T or C/G)
pub fn is_strand_ambiguous(ref_allele: char, alt_allele: char) -> bool {
    matches!(
        (ref_allele, alt_allele),
        ('A', 'T') | ('T', 'A') | ('C', 'G') | ('G', 'C')
    )
}

/// Match alleles between two SNPs at the same position
pub fn match_alleles(snp_a: &Snp, snp_b: &Snp) -> AlleleMatch {
    let (ref_a, alt_a) = (snp_a.allele_ref, snp_a.allele_alt);
    let (ref_b, alt_b) = (snp_b.allele_ref, snp_b.allele_alt);

    // Check for strand ambiguity (A/T, T/A, C/G, G/C)
    if is_strand_ambiguous(ref_a, alt_a) {
        return AlleleMatch::StrandAmbiguous;
    }

    // Exact match
    if ref_a == ref_b && alt_a == alt_b {
        return AlleleMatch::Identical;
    }

    // Flipped (need to recode genotypes: 0↔2)
    if ref_a == alt_b && alt_a == ref_b {
        return AlleleMatch::Flipped;
    }

    AlleleMatch::Mismatch
}

/// Recode genotype for flipped alleles (0↔2, 1→1, missing→missing)
pub fn recode_flipped(geno: Genotype) -> Genotype {
    match geno {
        Genotype::HomRef => Genotype::HomAlt,
        Genotype::HomAlt => Genotype::HomRef,
        Genotype::Het => Genotype::Het,
        Genotype::Missing => Genotype::Missing,
    }
}

/// Build SNP intersection between two datasets
pub fn build_snp_intersection(snps_a: &[Snp], snps_b: &[Snp]) -> (Vec<MergedSnp>, MergeStats) {
    let mut stats = MergeStats::default();

    // Build index of SNPs in file B by (chromosome, position)
    let mut b_index: HashMap<SnpKey, (usize, &Snp)> = HashMap::new();
    for (idx, snp) in snps_b.iter().enumerate() {
        let key = SnpKey::from_snp(snp);
        b_index.insert(key, (idx, snp));
    }

    let mut merged_snps = Vec::new();
    let mut matched_b_indices: HashSet<usize> = HashSet::new();

    // Find matching SNPs from file A
    for (a_idx, snp_a) in snps_a.iter().enumerate() {
        let key = SnpKey::from_snp(snp_a);

        if let Some(&(b_idx, snp_b)) = b_index.get(&key) {
            let allele_match = match_alleles(snp_a, snp_b);

            match allele_match {
                AlleleMatch::Identical => {
                    stats.snps_identical += 1;
                    merged_snps.push(MergedSnp {
                        snp: snp_a.clone(),
                        file_a_idx: a_idx,
                        file_b_idx: b_idx,
                        allele_match,
                    });
                }
                AlleleMatch::Flipped => {
                    stats.snps_flipped += 1;
                    merged_snps.push(MergedSnp {
                        snp: snp_a.clone(),
                        file_a_idx: a_idx,
                        file_b_idx: b_idx,
                        allele_match,
                    });
                }
                AlleleMatch::StrandAmbiguous => {
                    stats.snps_strand_ambiguous += 1;
                    // Skip this SNP
                }
                AlleleMatch::Mismatch => {
                    stats.snps_mismatch += 1;
                    // Skip this SNP
                }
            }
            matched_b_indices.insert(b_idx);
        } else {
            stats.snps_only_a += 1;
        }
    }

    // Count SNPs only in file B
    stats.snps_only_b = snps_b.len() - matched_b_indices.len();

    (merged_snps, stats)
}

/// Check for duplicate sample IDs between two sample lists
pub fn find_duplicate_samples(samples_a: &[Sample], samples_b: &[Sample]) -> Vec<String> {
    let ids_a: HashSet<&str> = samples_a.iter().map(|s| s.id.as_str()).collect();
    samples_b
        .iter()
        .filter(|s| ids_a.contains(s.id.as_str()))
        .map(|s| s.id.clone())
        .collect()
}

/// Merge two genotype matrices
///
/// Returns the merged matrix and statistics about the merge.
pub fn merge_genotypes(
    matrix_a: &GenotypeMatrix,
    matrix_b: &GenotypeMatrix,
    allow_strand_ambiguous: bool,
) -> Result<(GenotypeMatrix, MergeStats)> {
    // Check for duplicate samples
    let duplicates = find_duplicate_samples(matrix_a.samples(), matrix_b.samples());
    if !duplicates.is_empty() {
        return Err(PopmixError::Validation(format!(
            "Duplicate sample IDs found in both files: {}",
            duplicates.join(", ")
        )));
    }

    // Build SNP intersection
    let (merged_snps, mut stats) = build_snp_intersection(matrix_a.snps(), matrix_b.snps());

    if merged_snps.is_empty() {
        return Err(PopmixError::Validation(
            "No overlapping SNPs found between the two files".to_string(),
        ));
    }

    // Merge samples
    let mut merged_samples: Vec<Sample> = matrix_a.samples().to_vec();
    merged_samples.extend(matrix_b.samples().iter().cloned());
    stats.samples_total = merged_samples.len();

    // Build merged genotype matrix
    let n_snps = merged_snps.len();
    let n_samples_a = matrix_a.samples().len();
    let n_samples_b = matrix_b.samples().len();
    let n_samples = n_samples_a + n_samples_b;

    let mut merged_data = ndarray::Array2::<u8>::zeros((n_snps, n_samples));

    for (merged_idx, merged_snp) in merged_snps.iter().enumerate() {
        // Copy genotypes from file A
        for sample_idx in 0..n_samples_a {
            let geno = matrix_a.get(merged_snp.file_a_idx, sample_idx);
            merged_data[[merged_idx, sample_idx]] = geno as u8;
        }

        // Copy genotypes from file B (with potential recoding)
        for sample_idx in 0..n_samples_b {
            let mut geno = matrix_b.get(merged_snp.file_b_idx, sample_idx);

            // Recode if alleles are flipped
            if merged_snp.allele_match == AlleleMatch::Flipped {
                geno = recode_flipped(geno);
            }

            merged_data[[merged_idx, n_samples_a + sample_idx]] = geno as u8;
        }
    }

    // Extract merged SNP metadata
    let merged_snp_metadata: Vec<Snp> = merged_snps.iter().map(|ms| ms.snp.clone()).collect();

    Ok((
        GenotypeMatrix::new(merged_data, merged_samples, merged_snp_metadata)?,
        stats,
    ))
}

#[cfg(test)]
mod tests {
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
        // Different alleles entirely (A/G vs C/T)
        let snp_a = make_snp(1, 100, 'A', 'G');
        let snp_b = make_snp(1, 100, 'C', 'T');

        let result = match_alleles(&snp_a, &snp_b);
        assert_eq!(result, AlleleMatch::Mismatch);
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
        let (merged, stats) = merge_genotypes(&matrix_a, &matrix_b, false).unwrap();

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
        let result = merge_genotypes(&matrix_a, &matrix_b, false);
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
        let result = merge_genotypes(&matrix_a, &matrix_b, false);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("No overlapping SNPs"));
    }
}
