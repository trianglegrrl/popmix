//! Core genotype data types.
//!
//! This module defines the fundamental types for representing genetic data:
//! - `Genotype`: Individual genotype values (0, 1, 2, or missing)
//! - `Sample`: Individual sample with population assignment
//! - `Snp`: SNP metadata (position, alleles)
//! - `Ploidy`: Diploid vs pseudohaploid
//! - `Chromosome`: Chromosome identifier

use std::fmt;

/// A genotype value representing allele dosage.
///
/// Values:
/// - `HomRef` (0): Homozygous reference (0/0)
/// - `Het` (1): Heterozygous (0/1)
/// - `HomAlt` (2): Homozygous alternate (1/1)
/// - `Missing`: No data available
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Genotype {
    HomRef = 0,
    Het = 1,
    HomAlt = 2,
    Missing = 3,
}

impl Genotype {
    /// Create a genotype from PLINK 2-bit encoding.
    ///
    /// PLINK encoding:
    /// - 00 = Homozygous for first allele (typically HomAlt = 2)
    /// - 01 = Missing
    /// - 10 = Heterozygous
    /// - 11 = Homozygous for second allele (typically HomRef = 0)
    #[inline]
    pub fn from_plink_bits(bits: u8) -> Self {
        match bits & 0b11 {
            0b00 => Genotype::HomAlt,
            0b01 => Genotype::Missing,
            0b10 => Genotype::Het,
            0b11 => Genotype::HomRef,
            _ => unreachable!(),
        }
    }

    /// Create a genotype from EIGENSTRAT character encoding.
    ///
    /// EIGENSTRAT encoding:
    /// - '0' = Homozygous reference
    /// - '1' = Heterozygous
    /// - '2' = Homozygous alternate
    /// - '9' = Missing
    #[inline]
    pub fn from_eigenstrat_char(c: u8) -> Option<Self> {
        match c {
            b'0' => Some(Genotype::HomRef),
            b'1' => Some(Genotype::Het),
            b'2' => Some(Genotype::HomAlt),
            b'9' => Some(Genotype::Missing),
            _ => None,
        }
    }

    /// Convert to allele dosage (0.0, 1.0, 2.0) or None if missing.
    #[inline]
    pub fn to_dosage(self) -> Option<f64> {
        match self {
            Genotype::HomRef => Some(0.0),
            Genotype::Het => Some(1.0),
            Genotype::HomAlt => Some(2.0),
            Genotype::Missing => None,
        }
    }

    /// Check if this genotype is missing.
    #[inline]
    pub fn is_missing(self) -> bool {
        matches!(self, Genotype::Missing)
    }

    /// Check if this genotype is heterozygous.
    #[inline]
    pub fn is_het(self) -> bool {
        matches!(self, Genotype::Het)
    }

    /// Convert to EIGENSTRAT character encoding.
    ///
    /// Returns '0', '1', '2', or '9' for missing.
    #[inline]
    pub fn to_eigenstrat_char(self) -> u8 {
        match self {
            Genotype::HomRef => b'0',
            Genotype::Het => b'1',
            Genotype::HomAlt => b'2',
            Genotype::Missing => b'9',
        }
    }

    /// Convert to PLINK 2-bit encoding.
    ///
    /// Returns the 2-bit value for PLINK .bed format.
    /// Note: This is the inverse of `from_plink_bits`.
    #[inline]
    pub fn to_plink_bits(self) -> u8 {
        match self {
            Genotype::HomRef => 0b11,
            Genotype::Het => 0b10,
            Genotype::HomAlt => 0b00,
            Genotype::Missing => 0b01,
        }
    }

    /// Convert to PACKEDANCESTRYMAP 2-bit encoding.
    ///
    /// Encoding: 00=HomRef, 01=Het, 10=HomAlt, 11=Missing
    /// Note: This differs from PLINK encoding!
    #[inline]
    pub fn to_packed_bits(self) -> u8 {
        match self {
            Genotype::HomRef => 0b00,
            Genotype::Het => 0b01,
            Genotype::HomAlt => 0b10,
            Genotype::Missing => 0b11,
        }
    }

    /// Create a genotype from PACKEDANCESTRYMAP 2-bit encoding.
    ///
    /// Encoding: 00=HomRef, 01=Het, 10=HomAlt, 11=Missing
    #[inline]
    pub fn from_packed_bits(bits: u8) -> Self {
        match bits & 0b11 {
            0b00 => Genotype::HomRef,
            0b01 => Genotype::Het,
            0b10 => Genotype::HomAlt,
            0b11 => Genotype::Missing,
            _ => unreachable!(),
        }
    }
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Genotype::HomRef => write!(f, "0"),
            Genotype::Het => write!(f, "1"),
            Genotype::HomAlt => write!(f, "2"),
            Genotype::Missing => write!(f, "NA"),
        }
    }
}

/// Ploidy of a sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Ploidy {
    /// Diploid sample (can have heterozygous genotypes).
    Diploid,
    /// Pseudohaploid sample (no heterozygous genotypes observed).
    Pseudohaploid,
}

impl Ploidy {
    /// Get the effective ploidy for allele frequency calculations.
    ///
    /// Returns 2 for diploid, 1 for pseudohaploid.
    #[inline]
    pub fn effective_ploidy(self) -> u8 {
        match self {
            Ploidy::Diploid => 2,
            Ploidy::Pseudohaploid => 1,
        }
    }
}

/// Population identifier.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PopulationId(pub String);

impl PopulationId {
    pub fn new(id: impl Into<String>) -> Self {
        PopulationId(id.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for PopulationId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<&str> for PopulationId {
    fn from(s: &str) -> Self {
        PopulationId(s.to_string())
    }
}

impl From<String> for PopulationId {
    fn from(s: String) -> Self {
        PopulationId(s)
    }
}

/// Sample (individual) metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct Sample {
    /// Unique identifier for this sample.
    pub id: String,
    /// Population this sample belongs to.
    pub population: PopulationId,
    /// Ploidy of this sample.
    pub ploidy: Ploidy,
}

impl Sample {
    /// Create a new sample with default diploid ploidy.
    pub fn new(id: impl Into<String>, population: impl Into<PopulationId>) -> Self {
        Sample {
            id: id.into(),
            population: population.into(),
            ploidy: Ploidy::Diploid,
        }
    }

    /// Create a new sample with specified ploidy.
    pub fn with_ploidy(
        id: impl Into<String>,
        population: impl Into<PopulationId>,
        ploidy: Ploidy,
    ) -> Self {
        Sample {
            id: id.into(),
            population: population.into(),
            ploidy,
        }
    }
}

/// Chromosome identifier.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Chromosome {
    /// Autosomal chromosome (1-22).
    Autosome(u8),
    /// X chromosome.
    X,
    /// Y chromosome.
    Y,
    /// Mitochondrial DNA.
    MT,
    /// Unknown or non-standard chromosome.
    Other(String),
}

impl Chromosome {
    /// Parse chromosome from string.
    pub fn parse(s: &str) -> Self {
        let s = s.trim();
        // Handle common prefixes
        let s = s.strip_prefix("chr").unwrap_or(s);

        match s {
            "X" | "x" | "23" => Chromosome::X,
            "Y" | "y" | "24" => Chromosome::Y,
            "MT" | "Mt" | "mt" | "M" | "chrM" | "26" => Chromosome::MT,
            _ => {
                if let Ok(n) = s.parse::<u8>() {
                    if (1..=22).contains(&n) {
                        return Chromosome::Autosome(n);
                    }
                }
                Chromosome::Other(s.to_string())
            }
        }
    }

    /// Check if this is an autosomal chromosome.
    #[inline]
    pub fn is_autosome(&self) -> bool {
        matches!(self, Chromosome::Autosome(_))
    }
}

impl fmt::Display for Chromosome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Chromosome::Autosome(n) => write!(f, "{}", n),
            Chromosome::X => write!(f, "X"),
            Chromosome::Y => write!(f, "Y"),
            Chromosome::MT => write!(f, "MT"),
            Chromosome::Other(s) => write!(f, "{}", s),
        }
    }
}

/// SNP (Single Nucleotide Polymorphism) metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct Snp {
    /// SNP identifier (e.g., rs number).
    pub id: String,
    /// Chromosome.
    pub chromosome: Chromosome,
    /// Position in base pairs.
    pub position_bp: u32,
    /// Position in Morgans (genetic distance).
    /// Stored in Morgans to match DATES binning precision.
    pub position_morgans: f64,
    /// Reference allele.
    pub allele_ref: char,
    /// Alternate allele.
    pub allele_alt: char,
}

impl Snp {
    /// Create a new SNP.
    pub fn new(
        id: impl Into<String>,
        chromosome: Chromosome,
        position_bp: u32,
        position_morgans: f64,
        allele_ref: char,
        allele_alt: char,
    ) -> Self {
        Snp {
            id: id.into(),
            chromosome,
            position_bp,
            position_morgans,
            allele_ref,
            allele_alt,
        }
    }

    /// Check if this SNP is on an autosomal chromosome.
    #[inline]
    pub fn is_autosomal(&self) -> bool {
        self.chromosome.is_autosome()
    }
}

/// Check a list of population names: each at most once, and each one
/// `is_known` accepts. A typo or a population in two roles (qpAdm's
/// target among the right, say) is an error naming them, never an empty
/// result.
pub fn check_populations<'a>(
    requested: impl IntoIterator<Item = &'a str>,
    is_known: impl Fn(&str) -> bool,
) -> crate::Result<()> {
    let mut seen = std::collections::HashSet::new();
    let mut unknown = Vec::new();
    for pop in requested {
        if !seen.insert(pop) {
            return Err(crate::PopmixError::DuplicatePopulation(pop.to_string()));
        }
        if !is_known(pop) && !unknown.contains(&pop) {
            unknown.push(pop);
        }
    }
    if unknown.is_empty() {
        Ok(())
    } else {
        Err(crate::PopmixError::PopulationNotFound(unknown.join("', '")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn population_lists_name_each_population_once_and_only_known_ones() {
        let known = |p: &str| ["A", "B", "C"].contains(&p);
        assert!(check_populations(["A", "B"], known).is_ok());
        let dup = check_populations(["A", "B", "A"], known)
            .unwrap_err()
            .to_string();
        assert!(dup.contains("'A'"), "{dup}");
        let unknown = check_populations(["A", "X", "Y"], known)
            .unwrap_err()
            .to_string();
        assert!(unknown.contains("X") && unknown.contains("Y"), "{unknown}");
    }

    #[test]
    fn test_genotype_from_plink_bits() {
        assert_eq!(Genotype::from_plink_bits(0b00), Genotype::HomAlt);
        assert_eq!(Genotype::from_plink_bits(0b01), Genotype::Missing);
        assert_eq!(Genotype::from_plink_bits(0b10), Genotype::Het);
        assert_eq!(Genotype::from_plink_bits(0b11), Genotype::HomRef);
    }

    #[test]
    fn test_genotype_from_eigenstrat_char() {
        assert_eq!(Genotype::from_eigenstrat_char(b'0'), Some(Genotype::HomRef));
        assert_eq!(Genotype::from_eigenstrat_char(b'1'), Some(Genotype::Het));
        assert_eq!(Genotype::from_eigenstrat_char(b'2'), Some(Genotype::HomAlt));
        assert_eq!(
            Genotype::from_eigenstrat_char(b'9'),
            Some(Genotype::Missing)
        );
        assert_eq!(Genotype::from_eigenstrat_char(b'X'), None);
    }

    #[test]
    fn test_genotype_to_dosage() {
        assert_eq!(Genotype::HomRef.to_dosage(), Some(0.0));
        assert_eq!(Genotype::Het.to_dosage(), Some(1.0));
        assert_eq!(Genotype::HomAlt.to_dosage(), Some(2.0));
        assert_eq!(Genotype::Missing.to_dosage(), None);
    }

    #[test]
    fn test_chromosome_parse() {
        assert_eq!(Chromosome::parse("1"), Chromosome::Autosome(1));
        assert_eq!(Chromosome::parse("22"), Chromosome::Autosome(22));
        assert_eq!(Chromosome::parse("chr1"), Chromosome::Autosome(1));
        assert_eq!(Chromosome::parse("X"), Chromosome::X);
        assert_eq!(Chromosome::parse("chrX"), Chromosome::X);
        assert_eq!(Chromosome::parse("Y"), Chromosome::Y);
        assert_eq!(Chromosome::parse("MT"), Chromosome::MT);
        assert_eq!(Chromosome::parse("chrM"), Chromosome::MT);
        assert!(matches!(Chromosome::parse("unknown"), Chromosome::Other(_)));
    }

    #[test]
    fn test_sample_creation() {
        let sample = Sample::new("IND001", "French");
        assert_eq!(sample.id, "IND001");
        assert_eq!(sample.population.as_str(), "French");
        assert_eq!(sample.ploidy, Ploidy::Diploid);
    }

    #[test]
    fn test_snp_is_autosomal() {
        let auto_snp = Snp::new("rs1", Chromosome::Autosome(1), 100, 0.01, 'A', 'G');
        assert!(auto_snp.is_autosomal());

        let x_snp = Snp::new("rsX", Chromosome::X, 100, 0.01, 'A', 'G');
        assert!(!x_snp.is_autosomal());
    }

    #[test]
    fn test_genotype_to_eigenstrat_char() {
        assert_eq!(Genotype::HomRef.to_eigenstrat_char(), b'0');
        assert_eq!(Genotype::Het.to_eigenstrat_char(), b'1');
        assert_eq!(Genotype::HomAlt.to_eigenstrat_char(), b'2');
        assert_eq!(Genotype::Missing.to_eigenstrat_char(), b'9');
    }

    #[test]
    fn test_genotype_to_plink_bits() {
        // Verify round-trip
        assert_eq!(
            Genotype::from_plink_bits(Genotype::HomRef.to_plink_bits()),
            Genotype::HomRef
        );
        assert_eq!(
            Genotype::from_plink_bits(Genotype::Het.to_plink_bits()),
            Genotype::Het
        );
        assert_eq!(
            Genotype::from_plink_bits(Genotype::HomAlt.to_plink_bits()),
            Genotype::HomAlt
        );
        assert_eq!(
            Genotype::from_plink_bits(Genotype::Missing.to_plink_bits()),
            Genotype::Missing
        );
    }

    #[test]
    fn test_genotype_packed_bits_roundtrip() {
        // Verify round-trip for PACKEDANCESTRYMAP encoding
        assert_eq!(
            Genotype::from_packed_bits(Genotype::HomRef.to_packed_bits()),
            Genotype::HomRef
        );
        assert_eq!(
            Genotype::from_packed_bits(Genotype::Het.to_packed_bits()),
            Genotype::Het
        );
        assert_eq!(
            Genotype::from_packed_bits(Genotype::HomAlt.to_packed_bits()),
            Genotype::HomAlt
        );
        assert_eq!(
            Genotype::from_packed_bits(Genotype::Missing.to_packed_bits()),
            Genotype::Missing
        );
    }

    #[test]
    fn test_packed_bits_encoding() {
        // Verify specific encoding values
        assert_eq!(Genotype::HomRef.to_packed_bits(), 0b00);
        assert_eq!(Genotype::Het.to_packed_bits(), 0b01);
        assert_eq!(Genotype::HomAlt.to_packed_bits(), 0b10);
        assert_eq!(Genotype::Missing.to_packed_bits(), 0b11);
    }
}
