//! Genotype file merge algorithm.
//!
//! Merges two genotype datasets by SNP intersection, handling allele matching
//! and genotype recoding when alleles are flipped.
//!
//! SNPs are paired by position, and by name where several share one. Alleles
//! are reconciled as AdmixTools `mergeit` does (checked black-box against it,
//! tests/golden/formats/generate_merge.sh): a SNP reported on the other strand
//! is complemented, an A/T or C/G SNP is left out, a third allele is left out,
//! and the unknown allele of a monomorphic SNP is taken from the other file.

use std::collections::{HashMap, HashSet};

use popmix_core::error::{PopmixError, Result};
use popmix_core::genotype_matrix::GenotypeMatrix;
use popmix_core::genotypes::{Chromosome, Genotype, Sample, Snp};
use popmix_core::GenotypeRows;

use crate::streaming::SnpMajor;

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
    /// Of the SNPs kept, those reported on the other strand in file B
    pub snps_other_strand: usize,
    /// SNPs of file A left out because several SNPs of file B share their
    /// position and none has their name
    pub snps_duplicate_position: usize,
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

/// Is this the allele written for the unknown second allele of a monomorphic
/// SNP (PLINK `0`, EIGENSTRAT `X`, VCF `.`)?
fn is_unknown(allele: char) -> bool {
    matches!(allele, '0' | 'X' | '.' | 'N')
}

fn complement(allele: char) -> char {
    match allele {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// How a SNP of file B lines up with the same SNP of file A.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Reconciled {
    pub allele_match: AlleleMatch,
    /// File B reports the SNP on the other strand.
    pub other_strand: bool,
    /// The merged SNP's alleles, in file A's order, unknowns filled from B.
    pub alleles: (char, char),
}

/// Line up the alleles of `snp_b` with those of `snp_a`: as they are, swapped,
/// then complemented and either. An unknown allele matches any.
pub fn reconcile_alleles(snp_a: &Snp, snp_b: &Snp) -> Reconciled {
    let (ra, aa) = (snp_a.allele_ref, snp_a.allele_alt);
    let fits = |x: char, y: char| x == y || is_unknown(x) || is_unknown(y);
    let fill = |x: char, y: char| if is_unknown(x) { y } else { x };
    for other_strand in [false, true] {
        let flip = |c: char| if other_strand { complement(c) } else { c };
        let (rb, ab) = (flip(snp_b.allele_ref), flip(snp_b.allele_alt));
        for (allele_match, (r, a)) in [
            (AlleleMatch::Identical, (rb, ab)),
            (AlleleMatch::Flipped, (ab, rb)),
        ] {
            if fits(ra, r) && fits(aa, a) {
                let alleles = (fill(ra, r), fill(aa, a));
                let allele_match = if is_strand_ambiguous(alleles.0, alleles.1) {
                    AlleleMatch::StrandAmbiguous
                } else {
                    allele_match
                };
                return Reconciled {
                    allele_match,
                    other_strand,
                    alleles,
                };
            }
        }
    }
    let allele_match = if is_strand_ambiguous(ra, aa) {
        AlleleMatch::StrandAmbiguous
    } else {
        AlleleMatch::Mismatch
    };
    Reconciled {
        allele_match,
        other_strand: false,
        alleles: (ra, aa),
    }
}

/// Match alleles between two SNPs at the same position
pub fn match_alleles(snp_a: &Snp, snp_b: &Snp) -> AlleleMatch {
    reconcile_alleles(snp_a, snp_b).allele_match
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

    // SNPs of file B by (chromosome, position), in file order
    let mut b_index: HashMap<SnpKey, Vec<usize>> = HashMap::new();
    for (idx, snp) in snps_b.iter().enumerate() {
        b_index.entry(SnpKey::from_snp(snp)).or_default().push(idx);
    }

    let mut merged_snps = Vec::new();
    // SNPs of file B paired with, or considered for, a SNP of file A
    let mut seen_b: HashSet<usize> = HashSet::new();

    for (a_idx, snp_a) in snps_a.iter().enumerate() {
        let Some(candidates) = b_index.get(&SnpKey::from_snp(snp_a)) else {
            stats.snps_only_a += 1;
            continue;
        };
        let unused: Vec<usize> = candidates
            .iter()
            .copied()
            .filter(|b| !seen_b.contains(b))
            .collect();
        let b_idx = match unused.iter().find(|&&b| snps_b[b].id == snp_a.id) {
            Some(&b) => b,
            None if unused.len() == 1 => unused[0],
            None if unused.is_empty() => {
                stats.snps_only_a += 1;
                continue;
            }
            None => {
                stats.snps_duplicate_position += 1;
                seen_b.extend(unused);
                continue;
            }
        };
        seen_b.insert(b_idx);

        let reconciled = reconcile_alleles(snp_a, &snps_b[b_idx]);
        match reconciled.allele_match {
            AlleleMatch::Identical => stats.snps_identical += 1,
            AlleleMatch::Flipped => stats.snps_flipped += 1,
            AlleleMatch::StrandAmbiguous => {
                stats.snps_strand_ambiguous += 1;
                continue;
            }
            AlleleMatch::Mismatch => {
                stats.snps_mismatch += 1;
                continue;
            }
        }
        if reconciled.other_strand {
            stats.snps_other_strand += 1;
        }
        let (allele_ref, allele_alt) = reconciled.alleles;
        merged_snps.push(MergedSnp {
            snp: Snp {
                allele_ref,
                allele_alt,
                ..snp_a.clone()
            },
            file_a_idx: a_idx,
            file_b_idx: b_idx,
            allele_match: reconciled.allele_match,
        });
    }

    stats.snps_only_b = snps_b.len() - seen_b.len();

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

/// What a merge of files A and B keeps: the SNPs both have and can be
/// reconciled, and the samples of A then B, with the merge's counts. Fails if
/// a sample is in both files or no SNP is.
pub fn plan_merge(
    samples_a: &[Sample],
    snps_a: &[Snp],
    samples_b: &[Sample],
    snps_b: &[Snp],
) -> Result<(Vec<MergedSnp>, Vec<Sample>, MergeStats)> {
    let duplicates = find_duplicate_samples(samples_a, samples_b);
    if !duplicates.is_empty() {
        return Err(PopmixError::Validation(format!(
            "Duplicate sample IDs found in both files: {}",
            duplicates.join(", ")
        )));
    }
    let (merged_snps, mut stats) = build_snp_intersection(snps_a, snps_b);
    if merged_snps.is_empty() {
        return Err(PopmixError::Validation(
            "No overlapping SNPs found between the two files".to_string(),
        ));
    }
    let mut samples = samples_a.to_vec();
    samples.extend(samples_b.iter().cloned());
    stats.samples_total = samples.len();
    Ok((merged_snps, samples, stats))
}

/// Merge two genotype matrices
///
/// Returns the merged matrix and statistics about the merge. Strand-ambiguous
/// SNPs (A/T, C/G) are left out and counted in the statistics.
pub fn merge_genotypes(
    matrix_a: &GenotypeMatrix,
    matrix_b: &GenotypeMatrix,
) -> Result<(GenotypeMatrix, MergeStats)> {
    let (merged_snps, merged_samples, stats) = plan_merge(
        matrix_a.samples(),
        matrix_a.snps(),
        matrix_b.samples(),
        matrix_b.snps(),
    )?;

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

/// The merged genotypes, read in place from two SNP-major files one SNP at a
/// time: for each SNP of `snps`, the samples of file A then those of file B,
/// B's genotypes recoded where its alleles are swapped.
pub struct MergedRows<'a> {
    a: &'a SnpMajor,
    b: &'a SnpMajor,
    snps: &'a [MergedSnp],
    next: usize,
}

impl<'a> MergedRows<'a> {
    pub fn new(a: &'a SnpMajor, b: &'a SnpMajor, snps: &'a [MergedSnp]) -> Self {
        MergedRows {
            a,
            b,
            snps,
            next: 0,
        }
    }
}

impl GenotypeRows for MergedRows<'_> {
    fn next_row(&mut self, out: &mut [u8]) -> Result<()> {
        let snp = self
            .snps
            .get(self.next)
            .ok_or_else(|| PopmixError::Validation("read past the last merged SNP".to_string()))?;
        self.next += 1;
        let (na, nb) = (self.a.n_samples(), self.b.n_samples());
        let (out_a, out_b) = out[..na + nb].split_at_mut(na);
        let record = self.a.record(snp.file_a_idx);
        for (s, o) in out_a.iter_mut().enumerate() {
            *o = self.a.code(record, s);
        }
        let record = self.b.record(snp.file_b_idx);
        let flipped = snp.allele_match == AlleleMatch::Flipped;
        for (s, o) in out_b.iter_mut().enumerate() {
            let g = self.b.code(record, s);
            // 0 <-> 2; heterozygous and missing (1, 3) stay
            *o = if flipped && g != 1 && g != 3 {
                2 - g
            } else {
                g
            };
        }
        Ok(())
    }
}

impl crate::writers::SampleGenotypes for MergedRows<'_> {
    fn fill_samples(&self, samples: std::ops::Range<usize>, n_snps: usize, out: &mut [u8]) {
        let na = self.a.n_samples();
        let split = samples.end.clamp(samples.start, na.max(samples.start));
        let (out_a, out_b) = out.split_at_mut((split - samples.start) * n_snps);
        let snps = self.snps;
        if split > samples.start {
            self.a
                .fill_samples_of(samples.start..split, n_snps, out_a, &|j| {
                    (snps[j].file_a_idx, false)
                });
        }
        if samples.end > split {
            self.b
                .fill_samples_of(split - na..samples.end - na, n_snps, out_b, &|j| {
                    (
                        snps[j].file_b_idx,
                        snps[j].allele_match == AlleleMatch::Flipped,
                    )
                });
        }
    }
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
}
