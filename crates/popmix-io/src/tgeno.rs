//! TGENO: the packed genotypes of PACKEDANCESTRYMAP stored individual by
//! individual (AdmixTools `transpose_packed`; AADR from v66 on).
//!
//! Layout, determined black-box from files AdmixTools 8 convertf writes:
//! - Header: 48 bytes, "TGENO {n_samples} {n_snps} {ind_hash} {snp_hash}"
//!   padded with zeros; the hashes are those of the packed format
//!   ([`crate::packedancestrymap::header_hash`]).
//! - Then one record per individual, in .ind order, of ceil(n_snps / 4)
//!   bytes: 2 bits per SNP, first SNP in the high bits; 0, 1 and 2 copies of
//!   the .snp file's first allele, 3 missing.
//!
//! An individual's genotypes are contiguous, so reading a few individuals
//! reads only their records.

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;
use ndarray::Array2;

use popmix_core::error::{PopmixError, Result};
use popmix_core::genotype_matrix::GenotypeMatrix;
use popmix_core::genotypes::{Sample, Snp};

/// Length of the header, whatever the record length.
pub const HEADER_LEN: usize = 48;

/// Does this .geno file start with the TGENO magic word?
pub fn is_tgeno_format(geno_path: &str) -> Result<bool> {
    use std::io::Read;
    let mut file = File::open(geno_path).map_err(|e| PopmixError::io_with_path(e, geno_path))?;
    let mut magic = [0u8; 5];
    match file.read_exact(&mut magic) {
        Ok(()) => Ok(&magic == b"TGENO"),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(PopmixError::io_with_path(e, geno_path)),
    }
}

/// Reader for a TGENO dataset (.geno, .snp, .ind).
pub struct TgenoReader {
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    record_len: usize,
    geno: Mmap,
}

impl TgenoReader {
    /// Open a TGENO dataset by prefix: the .ind and .snp files are those of
    /// EIGENSTRAT. The header's counts and the file size are checked.
    pub fn new(prefix: impl AsRef<Path>) -> Result<Self> {
        let prefix = prefix.as_ref().to_string_lossy().to_string();
        let samples = crate::eigenstrat::read_ind_file(&format!("{prefix}.ind"))?;
        let snps = crate::eigenstrat::read_snp_file(&format!("{prefix}.snp"))?;
        let geno_path = format!("{prefix}.geno");
        if !is_tgeno_format(&geno_path)? {
            return Err(PopmixError::InvalidFileFormat {
                file: geno_path,
                details: "expected a TGENO header".to_string(),
            });
        }
        crate::packedancestrymap::check_header_counts(
            &geno_path,
            "TGENO",
            samples.len(),
            snps.len(),
        )?;
        let file = File::open(&geno_path).map_err(|e| PopmixError::io_with_path(e, &geno_path))?;
        // SAFETY: the map is only read; like every reader of these files,
        // popmix assumes no other process truncates the file while it is open.
        let geno =
            unsafe { Mmap::map(&file) }.map_err(|e| PopmixError::io_with_path(e, &geno_path))?;
        let record_len = snps.len().div_ceil(4);
        let expected = HEADER_LEN + samples.len() * record_len;
        if geno.len() != expected {
            return Err(PopmixError::UnexpectedFileSize {
                file: geno_path,
                expected: expected as u64,
                actual: geno.len() as u64,
            });
        }
        Ok(TgenoReader {
            samples,
            snps,
            record_len,
            geno,
        })
    }

    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }

    pub fn snps(&self) -> &[Snp] {
        &self.snps
    }

    /// The genotypes of the samples at `indices` (in that order), for every
    /// SNP; only their records are read.
    pub fn read_samples(&self, indices: &[usize]) -> Result<GenotypeMatrix> {
        self.read(indices, 0..self.snps.len())
    }

    /// Every sample's genotypes for a range of SNPs.
    pub fn read_snp_range(&self, range: std::ops::Range<usize>) -> Result<GenotypeMatrix> {
        let all: Vec<usize> = (0..self.samples.len()).collect();
        self.read(&all, range)
    }

    /// The samples at `indices` (in that order) and the SNPs in `snps`.
    pub fn read(&self, indices: &[usize], snps: std::ops::Range<usize>) -> Result<GenotypeMatrix> {
        crate::check_snp_range(&snps, self.snps.len())?;
        if let Some(&bad) = indices.iter().find(|&&i| i >= self.samples.len()) {
            return Err(PopmixError::Validation(format!(
                "sample {bad} out of range for {} samples",
                self.samples.len()
            )));
        }
        let mut data = Array2::<u8>::zeros((snps.len(), indices.len()));
        for (col, &sample) in indices.iter().enumerate() {
            let start = HEADER_LEN + sample * self.record_len;
            let record = &self.geno[start..start + self.record_len];
            for (row, snp) in snps.clone().enumerate() {
                data[[row, col]] = (record[snp / 4] >> (6 - 2 * (snp % 4))) & 0b11;
            }
        }
        let samples = indices.iter().map(|&i| self.samples[i].clone()).collect();
        GenotypeMatrix::new(data, samples, self.snps[snps].to_vec())
    }

    /// Every sample's genotypes.
    pub fn read_all(&self) -> Result<GenotypeMatrix> {
        let all: Vec<usize> = (0..self.samples.len()).collect();
        self.read_samples(&all)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_tgeno"
    );
    const SAME_AS_TEXT: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_eigenstrat"
    );

    #[test]
    fn reads_what_convertf_wrote() {
        // tests/fixtures/small_tgeno.geno is AdmixTools 8 convertf's TGENO of
        // small_eigenstrat
        assert!(is_tgeno_format(&format!("{FIXTURE}.geno")).unwrap());
        let tgeno = TgenoReader::new(FIXTURE).unwrap().read_all().unwrap();
        let text = crate::eigenstrat::EigenstratReader::new(SAME_AS_TEXT)
            .unwrap()
            .read_all()
            .unwrap();
        assert_eq!((tgeno.n_snps(), tgeno.n_samples()), (5, 3));
        for snp in 0..5 {
            for sample in 0..3 {
                assert_eq!(
                    tgeno.get(snp, sample),
                    text.get(snp, sample),
                    "{snp} {sample}"
                );
            }
        }
    }

    #[test]
    fn reads_a_subset_in_the_order_asked() {
        let reader = TgenoReader::new(FIXTURE).unwrap();
        let all = reader.read_all().unwrap();
        let two = reader.read_samples(&[2, 0]).unwrap();
        assert_eq!(two.samples()[0].id, "IND3");
        for snp in 0..5 {
            assert_eq!(two.get(snp, 0), all.get(snp, 2));
            assert_eq!(two.get(snp, 1), all.get(snp, 0));
        }
        assert!(reader.read_samples(&[3]).is_err());
        let range = reader.read_snp_range(1..4).unwrap();
        assert_eq!(range.snps()[0].id, "rs2");
        for (row, snp) in (1..4).enumerate() {
            assert_eq!(range.get(row, 1), all.get(snp, 1));
        }
        assert!(reader.read_snp_range(2..9).is_err());
    }

    #[test]
    fn written_files_are_byte_for_byte_convertf_s() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        let matrix = TgenoReader::new(FIXTURE).unwrap().read_all().unwrap();
        crate::writers::write_tgeno(&prefix, &matrix).unwrap();
        let ours = std::fs::read(prefix.with_extension("geno")).unwrap();
        let theirs = std::fs::read(format!("{FIXTURE}.geno")).unwrap();
        // the records are identical; the header differs only in spacing
        assert_eq!(ours[HEADER_LEN..], theirs[HEADER_LEN..]);
        let fields = |b: &[u8]| -> Vec<String> {
            String::from_utf8_lossy(&b[..HEADER_LEN])
                .trim_end_matches('\0')
                .split_whitespace()
                .map(String::from)
                .collect()
        };
        assert_eq!(fields(&ours), fields(&theirs));
        let back = TgenoReader::new(&prefix).unwrap().read_all().unwrap();
        for snp in 0..5 {
            for sample in 0..3 {
                assert_eq!(back.get(snp, sample), matrix.get(snp, sample));
            }
        }
    }

    #[test]
    fn a_packed_file_is_not_tgeno() {
        let packed = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/fixtures/small_eigenstrat.geno"
        );
        assert!(!is_tgeno_format(packed).unwrap());
    }

    #[test]
    fn header_counts_and_size_are_checked() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        for ext in ["geno", "snp", "ind"] {
            std::fs::copy(format!("{FIXTURE}.{ext}"), prefix.with_extension(ext)).unwrap();
        }
        let ind = prefix.with_extension("ind");
        let lines: Vec<String> = std::fs::read_to_string(&ind)
            .unwrap()
            .lines()
            .map(String::from)
            .collect();
        std::fs::write(&ind, lines[1..].join("\n") + "\n").unwrap();
        let err = TgenoReader::new(&prefix).err().unwrap().to_string();
        assert!(err.contains("header gives 3 samples"), "{err}");

        std::fs::write(&ind, lines.join("\n") + "\n").unwrap();
        let geno = prefix.with_extension("geno");
        let mut bytes = std::fs::read(&geno).unwrap();
        bytes.pop();
        std::fs::write(&geno, bytes).unwrap();
        assert!(TgenoReader::new(&prefix).is_err());
    }
}
