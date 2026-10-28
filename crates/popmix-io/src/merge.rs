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
    reconcile_alleles_on(snp_a, snp_b, true)
}

/// [`reconcile_alleles`], or with `strand_check` false, taking both files to
/// be on one strand (mergeit `strandcheck: NO`, TASKS.md H5): alleles are
/// matched as they are or swapped, never complemented, so an A/T or C/G SNP
/// is kept like any other and a complemented one is a mismatch.
pub fn reconcile_alleles_on(snp_a: &Snp, snp_b: &Snp, strand_check: bool) -> Reconciled {
    let (ra, aa) = (snp_a.allele_ref, snp_a.allele_alt);
    let ambiguous = |r: char, a: char| strand_check && is_strand_ambiguous(r, a);
    let strands: &[bool] = if strand_check {
        &[false, true]
    } else {
        &[false]
    };
    let fits = |x: char, y: char| x == y || is_unknown(x) || is_unknown(y);
    let fill = |x: char, y: char| if is_unknown(x) { y } else { x };
    for &other_strand in strands {
        let flip = |c: char| if other_strand { complement(c) } else { c };
        let (rb, ab) = (flip(snp_b.allele_ref), flip(snp_b.allele_alt));
        for (allele_match, (r, a)) in [
            (AlleleMatch::Identical, (rb, ab)),
            (AlleleMatch::Flipped, (ab, rb)),
        ] {
            if fits(ra, r) && fits(aa, a) {
                let alleles = (fill(ra, r), fill(aa, a));
                let allele_match = if ambiguous(alleles.0, alleles.1) {
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
    let allele_match = if ambiguous(ra, aa) {
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
    build_snp_intersection_on(snps_a, snps_b, true)
}

/// [`build_snp_intersection`], with alleles reconciled as
/// [`reconcile_alleles_on`] with `strand_check`.
pub fn build_snp_intersection_on(
    snps_a: &[Snp],
    snps_b: &[Snp],
    strand_check: bool,
) -> (Vec<MergedSnp>, MergeStats) {
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

        let reconciled = reconcile_alleles_on(snp_a, &snps_b[b_idx], strand_check);
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
    plan_merge_on(samples_a, snps_a, samples_b, snps_b, true)
}

/// [`plan_merge`], with alleles reconciled as [`reconcile_alleles_on`] with
/// `strand_check`.
pub fn plan_merge_on(
    samples_a: &[Sample],
    snps_a: &[Snp],
    samples_b: &[Sample],
    snps_b: &[Snp],
    strand_check: bool,
) -> Result<(Vec<MergedSnp>, Vec<Sample>, MergeStats)> {
    let duplicates = find_duplicate_samples(samples_a, samples_b);
    if !duplicates.is_empty() {
        return Err(PopmixError::data(format!(
            "Duplicate sample IDs found in both files: {}",
            duplicates.join(", ")
        )));
    }
    let (merged_snps, mut stats) = build_snp_intersection_on(snps_a, snps_b, strand_check);
    if merged_snps.is_empty() {
        return Err(PopmixError::data(
            "No overlapping SNPs found between the two files",
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
    merge_genotypes_on(matrix_a, matrix_b, true)
}

/// [`merge_genotypes`], with alleles reconciled as [`reconcile_alleles_on`]
/// with `strand_check`.
pub fn merge_genotypes_on(
    matrix_a: &GenotypeMatrix,
    matrix_b: &GenotypeMatrix,
    strand_check: bool,
) -> Result<(GenotypeMatrix, MergeStats)> {
    let (merged_snps, merged_samples, stats) = plan_merge_on(
        matrix_a.samples(),
        matrix_a.snps(),
        matrix_b.samples(),
        matrix_b.snps(),
        strand_check,
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
            .ok_or_else(|| PopmixError::Internal("read past the last merged SNP".to_string()))?;

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
mod tests;
