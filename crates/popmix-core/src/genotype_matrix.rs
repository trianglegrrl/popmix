//! Genotype matrix storage and operations.
//!
//! This module provides efficient storage for genotype data, supporting
//! both in-memory and memory-mapped access patterns.

use ndarray::{Array2, ArrayView1, Axis};
use rayon::prelude::*;

use crate::error::{PopmixError, Result};
use crate::genotypes::{Chromosome, Genotype, Ploidy, PopulationId, Sample, Snp};

/// A matrix of genotypes with sample and SNP metadata.
///
/// Dimensions: [n_snps, n_samples] (SNP-major order for efficient iteration)
#[derive(Debug, Clone)]
pub struct GenotypeMatrix {
    /// Genotype data as u8 (0, 1, 2, 3 for missing)
    /// Shape: [n_snps, n_samples]
    data: Array2<u8>,
    /// Sample metadata
    samples: Vec<Sample>,
    /// SNP metadata
    snps: Vec<Snp>,
}

impl GenotypeMatrix {
    /// Create a new genotype matrix from raw data.
    pub fn new(data: Array2<u8>, samples: Vec<Sample>, snps: Vec<Snp>) -> Result<Self> {
        let (n_snps, n_samples) = data.dim();

        if n_samples != samples.len() {
            return Err(PopmixError::DimensionMismatch {
                context: "GenotypeMatrix creation".to_string(),
                expected: format!("{} samples", samples.len()),
                actual: format!("{} columns", n_samples),
            });
        }

        if n_snps != snps.len() {
            return Err(PopmixError::DimensionMismatch {
                context: "GenotypeMatrix creation".to_string(),
                expected: format!("{} SNPs", snps.len()),
                actual: format!("{} rows", n_snps),
            });
        }

        Ok(GenotypeMatrix {
            data,
            samples,
            snps,
        })
    }

    /// Get the number of samples.
    #[inline]
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }

    /// Get the number of SNPs.
    #[inline]
    pub fn n_snps(&self) -> usize {
        self.snps.len()
    }

    /// Get sample metadata.
    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }

    /// The samples, to relabel their populations.
    pub fn samples_mut(&mut self) -> &mut [Sample] {
        &mut self.samples
    }

    /// Get SNP metadata.
    pub fn snps(&self) -> &[Snp] {
        &self.snps
    }

    /// Get genotype at a specific position.
    #[inline]
    pub fn get(&self, snp_idx: usize, sample_idx: usize) -> Genotype {
        match self.data[[snp_idx, sample_idx]] {
            0 => Genotype::HomRef,
            1 => Genotype::Het,
            2 => Genotype::HomAlt,
            _ => Genotype::Missing,
        }
    }

    /// Get all genotypes for a single SNP.
    pub fn snp_genotypes(&self, snp_idx: usize) -> Vec<Genotype> {
        self.data
            .row(snp_idx)
            .iter()
            .map(|&g| match g {
                0 => Genotype::HomRef,
                1 => Genotype::Het,
                2 => Genotype::HomAlt,
                _ => Genotype::Missing,
            })
            .collect()
    }

    /// Get raw genotype row as u8 slice (0, 1, 2, 3=missing).
    pub fn snp_raw(&self, snp_idx: usize) -> ArrayView1<'_, u8> {
        self.data.row(snp_idx)
    }

    /// Detect ploidy for each sample based on heterozygosity.
    ///
    /// Samples with no heterozygous calls in the first `n_snps_check` SNPs
    /// are marked as pseudohaploid.
    pub fn detect_ploidy(&mut self, n_snps_check: usize) {
        let n_check = n_snps_check.min(self.n_snps());

        for (sample_idx, sample) in self.samples.iter_mut().enumerate() {
            let has_het = (0..n_check).any(|snp_idx| self.data[[snp_idx, sample_idx]] == 1);

            sample.ploidy = if has_het {
                Ploidy::Diploid
            } else {
                Ploidy::Pseudohaploid
            };
        }
    }

    /// Mark every sample diploid.
    pub fn set_all_diploid(&mut self) {
        self.samples
            .iter_mut()
            .for_each(|s| s.ploidy = Ploidy::Diploid);
    }

    /// Get unique population IDs in order of first appearance.
    pub fn populations(&self) -> Vec<PopulationId> {
        let mut seen = std::collections::HashSet::new();
        let mut pops = Vec::new();

        for sample in &self.samples {
            if seen.insert(sample.population.clone()) {
                pops.push(sample.population.clone());
            }
        }

        pops
    }

    /// Get sample indices for a specific population.
    pub fn population_indices(&self, pop: &PopulationId) -> Vec<usize> {
        self.samples
            .iter()
            .enumerate()
            .filter(|(_, s)| &s.population == pop)
            .map(|(i, _)| i)
            .collect()
    }

    /// Filter to only autosomal SNPs.
    pub fn filter_autosomal(&self) -> Result<Self> {
        let auto_indices: Vec<usize> = self
            .snps
            .iter()
            .enumerate()
            .filter(|(_, snp)| snp.is_autosomal())
            .map(|(i, _)| i)
            .collect();

        if auto_indices.is_empty() {
            return Err(PopmixError::NoSnpsRemaining);
        }

        let new_data = self.data.select(Axis(0), &auto_indices);
        let new_snps: Vec<Snp> = auto_indices.iter().map(|&i| self.snps[i].clone()).collect();

        Ok(GenotypeMatrix {
            data: new_data.to_owned(),
            samples: self.samples.clone(),
            snps: new_snps,
        })
    }

    /// Compute allele frequency for a single SNP across specified samples.
    ///
    /// Returns (frequency, count) where count is number of observed alleles.
    /// A pseudohaploid sample's call (0 or 2) is one allele, as in admixtools2
    /// and the streaming reader: it adds `geno / (3 - ploidy)` copies.
    pub fn allele_frequency(&self, snp_idx: usize, sample_indices: &[usize]) -> (f64, u32) {
        self.allele_frequency_counting(snp_idx, sample_indices, false)
    }

    /// [`Self::allele_frequency`], or with `all_diploid` every sample counted
    /// as diploid whatever its ploidy (a call of 2 is two copies), as DATES
    /// counts (DATES-0004).
    pub fn allele_frequency_counting(
        &self,
        snp_idx: usize,
        sample_indices: &[usize],
        all_diploid: bool,
    ) -> (f64, u32) {
        let mut sum = 0.0f64;
        let mut count = 0u32;

        for &sample_idx in sample_indices {
            let geno = self.data[[snp_idx, sample_idx]];
            if geno <= 2 {
                // Not missing
                let ploidy = if all_diploid {
                    2
                } else {
                    self.samples[sample_idx].ploidy.effective_ploidy() as u32
                };
                sum += geno as f64 / (3 - ploidy) as f64;
                count += ploidy;
            }
        }

        if count == 0 {
            (f64::NAN, 0)
        } else {
            (sum / count as f64, count)
        }
    }
}

/// Allele frequency matrix for all populations.
#[derive(Debug, Clone)]
pub struct AlleleFrequencies {
    /// Allele frequencies: [n_snps, n_pops]
    /// NaN for missing (no observations)
    pub frequencies: Array2<f64>,
    /// Observation counts: [n_snps, n_pops]
    /// Number of observed alleles (not individuals)
    pub counts: Array2<u32>,
    /// Population names (column labels)
    pub populations: Vec<PopulationId>,
    /// SNP metadata (row labels)
    pub snps: Vec<Snp>,
}

impl AlleleFrequencies {
    /// Create AlleleFrequencies directly from precomputed data.
    ///
    /// Used primarily for testing and when frequencies are computed externally.
    pub fn new(
        frequencies: Array2<f64>,
        counts: Array2<u32>,
        populations: Vec<PopulationId>,
        chromosomes: Vec<Chromosome>,
        positions: Vec<f64>,
    ) -> Self {
        // Create SNPs from chromosomes and positions
        let snps: Vec<Snp> = chromosomes
            .into_iter()
            .zip(positions)
            .enumerate()
            .map(|(i, (chr, pos))| {
                use crate::genotypes::Snp;
                Snp::new(format!("snp{}", i), chr, i as u32 * 1000, pos, 'A', 'T')
            })
            .collect();

        AlleleFrequencies {
            frequencies,
            counts,
            populations,
            snps,
        }
    }

    /// Compute allele frequencies from a genotype matrix.
    ///
    /// Frequencies are computed per population, with ploidy-aware counting.
    pub fn from_genotypes(geno: &GenotypeMatrix) -> Self {
        Self::from_genotypes_for_populations(geno, None)
    }

    /// Compute allele frequencies for specific populations only.
    ///
    /// This is much more efficient when you only need a few populations
    /// from a large dataset.
    ///
    /// If `population_filter` is None, computes for all populations.
    pub fn from_genotypes_for_populations(
        geno: &GenotypeMatrix,
        population_filter: Option<&[&str]>,
    ) -> Self {
        Self::from_genotypes_counting(geno, population_filter, false)
    }

    /// [`Self::from_genotypes_for_populations`], or with `all_diploid` every
    /// sample counted as diploid (see
    /// [`GenotypeMatrix::allele_frequency_counting`]).
    pub fn from_genotypes_counting(
        geno: &GenotypeMatrix,
        population_filter: Option<&[&str]>,
        all_diploid: bool,
    ) -> Self {
        let all_populations = geno.populations();

        // Filter to requested populations
        let populations: Vec<PopulationId> = match population_filter {
            Some(filter) => {
                let filter_set: std::collections::HashSet<&str> = filter.iter().copied().collect();
                all_populations
                    .iter()
                    .filter(|p| filter_set.contains(p.as_str()))
                    .cloned()
                    .collect()
            }
            None => all_populations,
        };

        let n_pops = populations.len();
        let n_snps = geno.n_snps();

        // Pre-compute population indices
        let pop_indices: Vec<Vec<usize>> = populations
            .iter()
            .map(|pop| geno.population_indices(pop))
            .collect();

        // Compute frequencies in parallel across SNPs
        let results: Vec<(Vec<f64>, Vec<u32>)> = (0..n_snps)
            .into_par_iter()
            .map(|snp_idx| {
                let mut freqs = Vec::with_capacity(n_pops);
                let mut counts = Vec::with_capacity(n_pops);

                for indices in &pop_indices {
                    let (freq, count) =
                        geno.allele_frequency_counting(snp_idx, indices, all_diploid);
                    freqs.push(freq);
                    counts.push(count);
                }

                (freqs, counts)
            })
            .collect();

        // Assemble into arrays
        let mut frequencies = Array2::zeros((n_snps, n_pops));
        let mut counts = Array2::zeros((n_snps, n_pops));

        for (snp_idx, (freqs, cnts)) in results.into_iter().enumerate() {
            for (pop_idx, (freq, count)) in freqs.into_iter().zip(cnts).enumerate() {
                frequencies[[snp_idx, pop_idx]] = freq;
                counts[[snp_idx, pop_idx]] = count;
            }
        }

        AlleleFrequencies {
            frequencies,
            counts,
            populations,
            snps: geno.snps().to_vec(),
        }
    }

    /// Get the number of SNPs.
    #[inline]
    pub fn n_snps(&self) -> usize {
        self.snps.len()
    }

    /// Get the number of populations.
    #[inline]
    pub fn n_populations(&self) -> usize {
        self.populations.len()
    }

    /// Get population index by name.
    pub fn pop_index(&self, pop: &str) -> Option<usize> {
        self.populations.iter().position(|p| p.as_str() == pop)
    }

    /// Get allele frequency for a specific SNP and population.
    #[inline]
    pub fn get(&self, snp_idx: usize, pop_idx: usize) -> Option<f64> {
        let freq = self.frequencies[[snp_idx, pop_idx]];
        if freq.is_nan() {
            None
        } else {
            Some(freq)
        }
    }

    /// Get observation count for a specific SNP and population.
    #[inline]
    pub fn get_count(&self, snp_idx: usize, pop_idx: usize) -> u32 {
        self.counts[[snp_idx, pop_idx]]
    }

    /// Filter to only include specified populations.
    ///
    /// Returns a new AlleleFrequencies with only the specified populations,
    /// in that order; an unknown or repeated population is an error.
    pub fn filter_populations(&self, populations: &[&str]) -> Result<Self> {
        crate::genotypes::check_populations(populations.iter().copied(), |p| {
            self.pop_index(p).is_some()
        })?;
        let indices: Vec<(usize, usize)> = populations
            .iter()
            .enumerate()
            .filter_map(|(new_idx, pop)| self.pop_index(pop).map(|old_idx| (new_idx, old_idx)))
            .collect();

        let n_snps = self.n_snps();
        let n_pops = indices.len();

        // Extract frequency and count columns
        let mut frequencies = Array2::zeros((n_snps, n_pops));
        let mut counts = Array2::zeros((n_snps, n_pops));
        let mut new_populations = Vec::with_capacity(n_pops);

        for (new_idx, old_idx) in &indices {
            frequencies
                .slice_mut(ndarray::s![.., *new_idx])
                .assign(&self.frequencies.slice(ndarray::s![.., *old_idx]));
            counts
                .slice_mut(ndarray::s![.., *new_idx])
                .assign(&self.counts.slice(ndarray::s![.., *old_idx]));
            new_populations.push(self.populations[*old_idx].clone());
        }

        Ok(AlleleFrequencies {
            frequencies,
            counts,
            populations: new_populations,
            snps: self.snps.clone(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::Chromosome;
    use ndarray::array;

    fn make_test_samples() -> Vec<Sample> {
        vec![
            Sample::new("IND1", "PopA"),
            Sample::new("IND2", "PopA"),
            Sample::new("IND3", "PopB"),
            Sample::new("IND4", "PopB"),
        ]
    }

    fn make_test_snps() -> Vec<Snp> {
        vec![
            Snp::new("rs1", Chromosome::Autosome(1), 100, 0.01, 'A', 'G'),
            Snp::new("rs2", Chromosome::Autosome(1), 200, 0.02, 'C', 'T'),
            Snp::new("rs3", Chromosome::X, 300, 0.03, 'A', 'T'),
        ]
    }

    #[test]
    fn a_pseudohaploid_call_counts_one_allele() {
        // PopA: IND1 diploid het (1 of 2), IND2 pseudohaploid "2" (1 of 1)
        let data = array![[1, 2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
        let mut samples = make_test_samples();
        samples[1].ploidy = Ploidy::Pseudohaploid;
        let geno = GenotypeMatrix::new(data, samples, make_test_snps()).unwrap();
        let (freq, count) = geno.allele_frequency(0, &[0, 1]);
        assert_eq!(count, 3);
        assert!((freq - 2.0 / 3.0).abs() < 1e-15, "{freq}");
    }

    #[test]
    fn test_genotype_matrix_creation() {
        let data = array![[0, 1, 2, 3], [1, 1, 0, 2], [2, 0, 1, 1]];
        let samples = make_test_samples();
        let snps = make_test_snps();

        let geno = GenotypeMatrix::new(data, samples, snps).unwrap();

        assert_eq!(geno.n_samples(), 4);
        assert_eq!(geno.n_snps(), 3);
    }

    #[test]
    fn test_genotype_access() {
        let data = array![[0, 1, 2, 3], [1, 1, 0, 2], [2, 0, 1, 1]];
        let geno = GenotypeMatrix::new(data, make_test_samples(), make_test_snps()).unwrap();

        assert_eq!(geno.get(0, 0), Genotype::HomRef);
        assert_eq!(geno.get(0, 1), Genotype::Het);
        assert_eq!(geno.get(0, 2), Genotype::HomAlt);
        assert_eq!(geno.get(0, 3), Genotype::Missing);
    }

    #[test]
    fn test_population_indices() {
        let data = array![[0, 1, 2, 0]];
        let snps = vec![Snp::new(
            "rs1",
            Chromosome::Autosome(1),
            100,
            0.01,
            'A',
            'G',
        )];
        let geno = GenotypeMatrix::new(data, make_test_samples(), snps).unwrap();

        let pop_a = geno.population_indices(&PopulationId::new("PopA"));
        let pop_b = geno.population_indices(&PopulationId::new("PopB"));

        assert_eq!(pop_a, vec![0, 1]);
        assert_eq!(pop_b, vec![2, 3]);
    }

    #[test]
    fn test_allele_frequency() {
        // SNP with genotypes: 0, 2, 0, 2 (all diploid)
        // PopA (indices 0,1): (0 + 2) / 4 = 0.5
        // PopB (indices 2,3): (0 + 2) / 4 = 0.5
        let data = array![[0, 2, 0, 2]];
        let snps = vec![Snp::new(
            "rs1",
            Chromosome::Autosome(1),
            100,
            0.01,
            'A',
            'G',
        )];
        let geno = GenotypeMatrix::new(data, make_test_samples(), snps).unwrap();

        let (freq, count) = geno.allele_frequency(0, &[0, 1]);
        assert!((freq - 0.5).abs() < 1e-10);
        assert_eq!(count, 4); // 2 diploid individuals = 4 alleles
    }

    #[test]
    fn test_allele_frequencies_from_genotypes() {
        let data = array![[0, 2, 0, 2], [1, 1, 2, 2]];
        let snps = vec![
            Snp::new("rs1", Chromosome::Autosome(1), 100, 0.01, 'A', 'G'),
            Snp::new("rs2", Chromosome::Autosome(1), 200, 0.02, 'C', 'T'),
        ];
        let geno = GenotypeMatrix::new(data, make_test_samples(), snps).unwrap();

        let af = AlleleFrequencies::from_genotypes(&geno);

        assert_eq!(af.n_snps(), 2);
        assert_eq!(af.n_populations(), 2);

        // SNP 0, PopA: (0 + 2) / 4 = 0.5
        assert!((af.get(0, 0).unwrap() - 0.5).abs() < 1e-10);

        // SNP 1, PopB: (2 + 2) / 4 = 1.0
        assert!((af.get(1, 1).unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_filter_autosomal() {
        let data = array![[0, 1, 2, 0], [1, 1, 0, 2], [2, 0, 1, 1]];
        let geno = GenotypeMatrix::new(data, make_test_samples(), make_test_snps()).unwrap();

        let auto = geno.filter_autosomal().unwrap();

        assert_eq!(auto.n_snps(), 2); // Only rs1 and rs2 are autosomal
        assert!(auto.snps().iter().all(|s| s.is_autosomal()));
    }
}
