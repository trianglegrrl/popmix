//! VCF (Variant Call Format) reader.
//!
//! VCF is a text format for storing genetic variants.
//! This module reads VCF files and converts them to genotype matrices.
//!
//! Implemented from the VCF v4.2 specification; which records are read
//! follows PLINK 2's documented `--snps-only just-acgt --max-alleles 2`.
//! The file is read once: the reading thread decompresses and splits lines,
//! batches of them are parsed in parallel, and the genotypes are kept packed
//! (2 bits each) rather than a byte per genotype.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use ndarray::Array2;

use crate::streaming::SnpMajor;
use popmix_core::{
    Chromosome, Genotype, GenotypeMatrix, Ploidy, PopmixError, PopulationId, Result, Sample, Snp,
};

/// Reader for VCF (Variant Call Format) files.
///
/// Supports both plain text (.vcf) and gzip-compressed (.vcf.gz) files.
/// Only SNPs with at most two alleles are read (a single-base REF, and a
/// single-base ALT or none, `.`); other records are skipped and counted
/// ([`VcfReader::skipped`]), as with PLINK 2's `--snps-only just-acgt
/// --max-alleles 2`. A site with no ALT is kept, monomorphic, with `.` as
/// its counted allele, as PLINK 2 and convertf write it.
pub struct VcfReader {
    path: String,
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    skipped: VcfSkipped,
    /// Every SNP's genotypes, packed as PACKEDANCESTRYMAP packs them
    /// (`samples.len().div_ceil(4)` bytes each).
    records: Vec<u8>,
}

/// Records a VCF reader skipped, by reason.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct VcfSkipped {
    /// More than one ALT allele.
    pub multiallelic: usize,
    /// An indel or a symbolic allele (not one base each).
    pub not_snp: usize,
}

impl VcfSkipped {
    /// Records skipped in all.
    pub fn total(&self) -> usize {
        self.multiallelic + self.not_snp
    }
}

/// Why a record is not a biallelic SNP, or None if it is one.
fn skip_reason(ref_allele: &str, alt: &str) -> Option<fn(&mut VcfSkipped)> {
    let base = |a: &str| {
        a.len() == 1
            && matches!(
                a.as_bytes()[0].to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T'
            )
    };
    if alt.contains(',') {
        Some(|k| k.multiallelic += 1)
    } else if !base(ref_allele) || !(base(alt) || alt == ".") {
        Some(|k| k.not_snp += 1)
    } else {
        None
    }
}

impl VcfReader {
    /// Open a VCF file for reading.
    ///
    /// The path can be either:
    /// - A plain text VCF file (.vcf)
    /// - A gzip-compressed VCF file (.vcf.gz)
    pub fn new(path: impl AsRef<Path>) -> Result<Self> {
        let path_str = path.as_ref().to_string_lossy().to_string();
        let is_gzipped = path_str.ends_with(".gz");
        let scan = scan_vcf(&path_str, is_gzipped)?;
        Ok(VcfReader {
            path: path_str,
            samples: scan.samples,
            snps: scan.snps,
            skipped: scan.skipped,
            records: scan.records,
        })
    }

    /// Get the number of samples.
    pub fn num_samples(&self) -> usize {
        self.samples.len()
    }

    /// Get the number of SNPs.
    pub fn num_snps(&self) -> usize {
        self.snps.len()
    }

    /// Get sample metadata.
    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }

    /// Get SNP metadata.
    pub fn snps(&self) -> &[Snp] {
        &self.snps
    }

    /// Records skipped because they are not biallelic SNPs.
    pub fn skipped(&self) -> VcfSkipped {
        self.skipped
    }

    /// Read all genotypes into a GenotypeMatrix.
    pub fn read_all(&self) -> Result<GenotypeMatrix> {
        use rayon::prelude::*;
        let n = self.samples.len();
        let bytes = n.div_ceil(4);
        let mut data = Array2::<u8>::zeros((self.snps.len(), n));
        data.axis_iter_mut(ndarray::Axis(0))
            .into_par_iter()
            .zip(self.records.par_chunks(bytes.max(1)))
            .for_each(|(mut row, record)| {
                for (s, g) in row.iter_mut().enumerate() {
                    *g = (record[s / 4] >> (6 - 2 * (s % 4))) & 0b11;
                }
            });
        GenotypeMatrix::new(data, self.samples.clone(), self.snps.clone())
    }

    /// The samples, the SNPs, and the genotypes as a [`SnpMajor`] file held
    /// in memory, to be read as a mapped one is.
    pub fn into_snp_major(self) -> (Vec<Sample>, Vec<Snp>, SnpMajor) {
        let file = SnpMajor::from_packed_records(
            self.path,
            self.records,
            self.samples.len(),
            self.snps.len(),
        );
        (self.samples, self.snps, file)
    }
}

/// Parse a VCF GT field into a Genotype.
///
/// Supports both phased (|) and unphased (/) separators.
///
/// VCF GT encoding:
/// - 0/0 or 0|0 = Homozygous reference
/// - 0/1 or 0|1 or 1/0 or 1|0 = Heterozygous
/// - 1/1 or 1|1 = Homozygous alternate
/// - ./. or .|. = Missing
///
/// Returns None if the GT field cannot be parsed.
pub fn parse_vcf_gt(gt: &str) -> Option<Genotype> {
    // Handle missing
    if gt.starts_with('.') {
        return Some(Genotype::Missing);
    }

    // Split on / or |
    let parts: Vec<&str> = if gt.contains('|') {
        gt.split('|').collect()
    } else if gt.contains('/') {
        gt.split('/').collect()
    } else {
        // A haploid call (male chrX, say) is homozygous, as plink2 reads it
        return match gt {
            "0" => Some(Genotype::HomRef),
            "1" => Some(Genotype::HomAlt),
            _ => None,
        };
    };

    if parts.len() != 2 {
        return None;
    }

    // Parse allele indices
    let a1 = parts[0].parse::<u8>().ok()?;
    let a2 = parts[1].parse::<u8>().ok()?;

    // Convert to genotype (count of ALT allele)
    let dosage = a1 + a2;
    match dosage {
        0 => Some(Genotype::HomRef), // 0/0
        1 => Some(Genotype::Het),    // 0/1 or 1/0
        2 => Some(Genotype::HomAlt), // 1/1
        _ => None,                   // Multi-allelic not supported
    }
}

/// One data line of a VCF, parsed.
enum Line {
    /// A biallelic SNP and its packed genotypes.
    Snp(Snp, Vec<u8>),
    /// A record that is not a biallelic SNP, counted by reason.
    Skipped(fn(&mut VcfSkipped)),
    /// A blank or `#` line.
    Blank,
}

/// Lines decompressed and split off by the reading thread, parsed a batch at
/// a time in parallel.
const BATCH_BYTES: usize = 16 << 20;

/// The genotype code (0, 1, 2; 3 missing) of a sample's GT, as
/// [`parse_vcf_gt`] reads it, without allocating for the usual `a|b`.
#[inline]
fn gt_code(gt: &[u8]) -> u8 {
    if let [a @ (b'0' | b'1'), b'|' | b'/', b @ (b'0' | b'1')] = gt {
        return (a - b'0') + (b - b'0');
    }
    std::str::from_utf8(gt)
        .ok()
        .and_then(parse_vcf_gt)
        .unwrap_or(Genotype::Missing) as u8
}

/// Parse data line `number` (1-based, for errors) of a VCF of `n_samples`
/// samples, packing its genotypes as PACKEDANCESTRYMAP does (the first
/// sample in a byte's high bits).
fn parse_line(path: &str, number: usize, line: &[u8], n_samples: usize) -> Result<Line> {
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    if line.is_empty() || line[0] == b'#' {
        return Ok(Line::Blank);
    }
    let text = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
    let mut fields = line.split(|&c| c == b'\t');
    let mut fixed: [&[u8]; 9] = [&[]; 9];
    let mut n_fixed = 0;
    for slot in fixed.iter_mut() {
        match fields.next() {
            Some(f) => {
                *slot = f;
                n_fixed += 1;
            }
            None => break,
        }
    }
    if n_fixed < 8 {
        return Err(PopmixError::parse(
            path,
            number,
            format!("VCF data line must have at least 8 columns, got {n_fixed}"),
        ));
    }
    let (chrom, pos_field) = (text(fixed[0]), text(fixed[1]));
    let pos = pos_field
        .parse::<u32>()
        .map_err(|_| PopmixError::parse(path, number, format!("Invalid position: {pos_field}")))?;
    let (ref_allele, alt) = (text(fixed[3]), text(fixed[4]));
    // every data line holds a column per sample
    let samples: Vec<&[u8]> = fields.take(n_samples).collect();
    if n_fixed < 9 || samples.len() < n_samples {
        return Err(PopmixError::parse(
            path,
            number,
            format!(
                "Expected {} columns, got {}",
                9 + n_samples,
                n_fixed + samples.len()
            ),
        ));
    }
    if let Some(count) = skip_reason(&ref_allele, &alt) {
        return Ok(Line::Skipped(count));
    }
    let gt_idx = fixed[8]
        .split(|&c| c == b':')
        .position(|f| f == b"GT")
        .ok_or_else(|| PopmixError::parse(path, number, "No GT field in FORMAT column"))?;
    let mut record = vec![0u8; n_samples.div_ceil(4)];
    for (s, field) in samples.iter().enumerate() {
        let code = match field.split(|&c| c == b':').nth(gt_idx) {
            Some(gt) => gt_code(gt),
            None => Genotype::Missing as u8,
        };
        record[s / 4] |= code << (6 - 2 * (s % 4));
    }
    let id = if fixed[2] == b"." {
        format!("{chrom}:{pos}")
    } else {
        text(fixed[2])
    };
    let ref_base = ref_allele.chars().next().unwrap_or('N');
    let alt_base = alt.chars().next().unwrap_or('N');
    // VCF doesn't have genetic position - use 0.0. Genotypes count ALT
    // copies, and the counted allele is listed first (EIGENSTRAT .snp
    // column 5, PLINK A1), as PLINK 2 and convertf list it.
    let snp = Snp::new(id, Chromosome::parse(&chrom), pos, 0.0, alt_base, ref_base);
    Ok(Line::Snp(snp, record))
}

/// A VCF read in one pass: its samples, SNPs, skipped records, and every
/// SNP's genotypes packed in `records`, `n_samples.div_ceil(4)` bytes each.
struct Scan {
    samples: Vec<Sample>,
    snps: Vec<Snp>,
    skipped: VcfSkipped,
    records: Vec<u8>,
}

/// Read a VCF once: the header on this thread, then the data lines,
/// decompressed and split into batches by a reading thread, parsed in
/// parallel and kept in file order.
fn scan_vcf(path: &str, is_gzipped: bool) -> Result<Scan> {
    use rayon::prelude::*;
    use std::sync::mpsc::sync_channel;

    let file = File::open(path).map_err(|e| PopmixError::io(e, path))?;
    let mut reader: Box<dyn BufRead + Send> = if is_gzipped {
        // bgzip (the usual .vcf.gz) is a series of gzip members
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };

    // the header, to the #CHROM line
    let mut samples = Vec::new();
    let mut number = 0;
    let mut line = Vec::new();
    loop {
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .map_err(|e| PopmixError::io(e, path))?;
        if n == 0 {
            break;
        }
        number += 1;
        let text = String::from_utf8_lossy(&line);
        let text = text.trim_end_matches(['\n', '\r']);
        if text.starts_with("##") || text.is_empty() {
            continue;
        }
        if text.starts_with("#CHROM") {
            let fields: Vec<&str> = text.split('\t').collect();
            if fields.len() < 10 {
                return Err(PopmixError::parse(
                    path,
                    number,
                    "VCF header must have at least 10 columns (including sample columns)",
                ));
            }
            // VCF doesn't have population info - use sample ID as population
            samples = fields[9..]
                .iter()
                .map(|name| {
                    Sample::with_ploidy(
                        name.to_string(),
                        PopulationId::new(name.to_string()),
                        Ploidy::Diploid,
                    )
                })
                .collect();
        }
        break;
    }
    if samples.is_empty() {
        return Err(PopmixError::parse(
            path,
            0,
            "No samples found in VCF header",
        ));
    }
    let n_samples = samples.len();

    let mut snps = Vec::new();
    let mut skipped = VcfSkipped::default();
    let mut records = Vec::new();
    // batches of whole lines: the bytes, where each line ends, the number
    // of the first
    type Batch = std::io::Result<(Vec<u8>, Vec<usize>, usize)>;
    let (tx, rx) = sync_channel::<Batch>(2);
    let first = number + 1;
    std::thread::scope(|scope| -> Result<()> {
        scope.spawn(move || {
            let mut next = first;
            loop {
                let mut buf = Vec::with_capacity(BATCH_BYTES + (1 << 16));
                let mut ends = Vec::new();
                while buf.len() < BATCH_BYTES {
                    match reader.read_until(b'\n', &mut buf) {
                        Ok(0) => break,
                        Ok(_) => ends.push(buf.len()),
                        Err(e) => {
                            let _ = tx.send(Err(e));
                            return;
                        }
                    }
                }
                if ends.is_empty() {
                    return;
                }
                let n = ends.len();
                if tx.send(Ok((buf, ends, next))).is_err() {
                    return;
                }
                next += n;
            }
        });
        for batch in rx {
            let (buf, ends, start) = batch.map_err(|e| PopmixError::io(e, path))?;
            let parsed: Vec<Result<Line>> = (0..ends.len())
                .into_par_iter()
                .map(|i| {
                    let lo = if i == 0 { 0 } else { ends[i - 1] };
                    let line = &buf[lo..ends[i]];
                    let line = line.strip_suffix(b"\n").unwrap_or(line);
                    parse_line(path, start + i, line, n_samples)
                })
                .collect();
            for line in parsed {
                match line? {
                    Line::Snp(snp, record) => {
                        snps.push(snp);
                        records.extend_from_slice(&record);
                    }
                    Line::Skipped(count) => count(&mut skipped),
                    Line::Blank => {}
                }
            }
        }
        Ok(())
    })?;

    Ok(Scan {
        samples,
        snps,
        skipped,
        records,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Fixture paths relative to workspace root
    const VCF_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_vcf.vcf"
    );
    const VCF_MISSING_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_vcf_missing.vcf"
    );
    const VCF_PHASED_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_vcf_phased.vcf"
    );
    const EIGENSTRAT_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_eigenstrat"
    );

    /// Every genotype of `reader`'s in-memory records, decoded.
    fn decoded(reader: VcfReader) -> Vec<Vec<u8>> {
        let (_, snps, file) = reader.into_snp_major();
        (0..snps.len())
            .map(|i| {
                let r = file.record(i);
                (0..file.n_samples()).map(|s| file.code(r, s)).collect()
            })
            .collect()
    }

    #[test]
    fn records_held_in_memory_equal_read_all() {
        for path in [
            VCF_FIXTURE,
            VCF_MISSING_FIXTURE,
            VCF_PHASED_FIXTURE,
            concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/fixtures/small_vcf_haploid.vcf"
            ),
            concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/fixtures/small_vcf_multiallelic.vcf"
            ),
        ] {
            let geno = VcfReader::new(path).unwrap().read_all().unwrap();
            let want: Vec<Vec<u8>> = (0..geno.n_snps())
                .map(|i| geno.snp_raw(i).to_vec())
                .collect();
            assert_eq!(decoded(VcfReader::new(path).unwrap()), want, "{path}");
        }
    }

    #[test]
    fn odd_calls_crlf_lines_and_gzip_read_as_the_rules_say() {
        let dir = tempfile::tempdir().unwrap();
        let text = "##fileformat=VCFv4.2\r\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\te\r\n\
1\t10\t.\tA\tG\t.\t.\t.\tGT:DP\t0|1:3\t1/1\t./.:0\t0/.\t1\r\n\
1\t20\trs2\tA\tAT\t.\t.\t.\tGT\t0|1\t1/1\t0/0\t0/0\t0\r\n\
1\t30\trs3\tC\tT\t.\t.\t.\tDP:GT\t3:0/0\t4:1|0\t5\t6:10/0\t7:0\r\n";
        let vcf = dir.path().join("x.vcf");
        std::fs::write(&vcf, text).unwrap();
        let gz = dir.path().join("x.vcf.gz");
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut enc, text.as_bytes()).unwrap();
        std::fs::write(&gz, enc.finish().unwrap()).unwrap();
        for path in [&vcf, &gz] {
            let reader = VcfReader::new(path).unwrap();
            assert_eq!(reader.skipped().not_snp, 1);
            let ids: Vec<&str> = reader.snps().iter().map(|s| s.id.as_str()).collect();
            assert_eq!(ids, ["1:10", "rs3"]);
            // 0/. and 10/0 are not calls popmix reads: missing; a lone 0 is HomRef
            assert_eq!(
                decoded(reader),
                vec![vec![1, 2, 3, 3, 2], vec![0, 1, 3, 3, 0]]
            );
        }
    }

    // TDD RED: Test GT parsing
    #[test]
    fn test_parse_vcf_gt_homref() {
        assert_eq!(parse_vcf_gt("0/0"), Some(Genotype::HomRef));
        assert_eq!(parse_vcf_gt("0|0"), Some(Genotype::HomRef));
    }

    #[test]
    fn test_parse_vcf_gt_het() {
        assert_eq!(parse_vcf_gt("0/1"), Some(Genotype::Het));
        assert_eq!(parse_vcf_gt("1/0"), Some(Genotype::Het));
        assert_eq!(parse_vcf_gt("0|1"), Some(Genotype::Het));
        assert_eq!(parse_vcf_gt("1|0"), Some(Genotype::Het));
    }

    #[test]
    fn test_parse_vcf_gt_homalt() {
        assert_eq!(parse_vcf_gt("1/1"), Some(Genotype::HomAlt));
        assert_eq!(parse_vcf_gt("1|1"), Some(Genotype::HomAlt));
    }

    #[test]
    fn test_parse_vcf_gt_missing() {
        assert_eq!(parse_vcf_gt("./."), Some(Genotype::Missing));
        assert_eq!(parse_vcf_gt(".|."), Some(Genotype::Missing));
        assert_eq!(parse_vcf_gt("."), Some(Genotype::Missing));
    }

    #[test]
    fn test_parse_vcf_gt_invalid() {
        assert_eq!(parse_vcf_gt(""), None);
        assert_eq!(parse_vcf_gt("invalid"), None);
    }

    #[test]
    fn a_haploid_call_is_homozygous() {
        // as plink2 reads one (male chrX, for instance)
        assert_eq!(parse_vcf_gt("0"), Some(Genotype::HomRef));
        assert_eq!(parse_vcf_gt("1"), Some(Genotype::HomAlt));
        assert_eq!(parse_vcf_gt("2"), None);
    }

    // TDD RED: Test VCF reader with fixture
    #[test]
    fn test_vcf_reader_samples() {
        let reader = VcfReader::new(VCF_FIXTURE).unwrap();
        assert_eq!(reader.num_samples(), 3);

        let sample_ids: Vec<&str> = reader.samples().iter().map(|s| s.id.as_str()).collect();
        assert_eq!(sample_ids, vec!["IND1", "IND2", "IND3"]);
    }

    #[test]
    fn test_bgzipped_vcf_reads_every_member() {
        // bgzip writes a series of gzip members; split the fixture across two
        use flate2::{write::GzEncoder, Compression};
        use std::io::Write;
        let text = std::fs::read_to_string(VCF_FIXTURE).unwrap();
        let cut = text
            .match_indices('\n')
            .nth(text.lines().count() - 3)
            .unwrap()
            .0
            + 1;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.vcf.gz");
        let mut bytes = Vec::new();
        for part in [&text[..cut], &text[cut..]] {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(part.as_bytes()).unwrap();
            bytes.extend(enc.finish().unwrap());
        }
        std::fs::write(&path, bytes).unwrap();
        let reader = VcfReader::new(path.to_str().unwrap()).unwrap();
        assert_eq!(reader.num_snps(), 5);
        let a = reader.read_all().unwrap();
        let b = VcfReader::new(VCF_FIXTURE).unwrap().read_all().unwrap();
        assert_eq!(a.n_snps(), b.n_snps());
        assert!((0..a.n_snps()).all(|j| a.snp_raw(j) == b.snp_raw(j)));
    }

    const VCF_MULTIALLELIC_FIXTURE: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/fixtures/small_vcf_multiallelic.vcf"
    );

    #[test]
    fn only_biallelic_snps_are_read() {
        // rs2 has two ALT alleles, rs3 and rs7 are indels, rs5 is a
        // symbolic allele: skipped and counted, not miscoded. rs4 has no ALT:
        // kept, monomorphic, as PLINK 2 keeps it
        let reader = VcfReader::new(VCF_MULTIALLELIC_FIXTURE).unwrap();
        let ids: Vec<&str> = reader.snps().iter().map(|s| s.id.as_str()).collect();
        assert_eq!(ids, ["rs1", "rs4", "rs6"]);
        assert_eq!(
            reader.skipped(),
            VcfSkipped {
                multiallelic: 1,
                not_snp: 3
            }
        );
        let m = reader.read_all().unwrap();
        assert_eq!(m.n_snps(), 3);
        assert_eq!(reader.snps()[1].allele_ref, '.');
        // rs6: C>T with 1/1, 0/1, 0/0, counted on ALT
        assert_eq!(m.get(2, 0), Genotype::HomAlt);
        assert_eq!(m.get(2, 1), Genotype::Het);
        assert_eq!(m.get(2, 2), Genotype::HomRef);
    }

    #[test]
    fn test_vcf_reader_snps() {
        let reader = VcfReader::new(VCF_FIXTURE).unwrap();
        assert_eq!(reader.num_snps(), 5);

        // Check first SNP
        let snp = &reader.snps()[0];
        assert_eq!(snp.id, "rs1");
        assert_eq!(snp.chromosome, Chromosome::Autosome(1));
        assert_eq!(snp.position_bp, 100);
        // REF A, ALT G: the counted allele (listed first, as PLINK 2 and
        // convertf list it after importing a VCF) is ALT
        assert_eq!(snp.allele_ref, 'G');
        assert_eq!(snp.allele_alt, 'A');
    }

    #[test]
    fn test_vcf_reader_genotypes() {
        let reader = VcfReader::new(VCF_FIXTURE).unwrap();
        let matrix = reader.read_all().unwrap();

        // Check dimensions
        assert_eq!(matrix.n_snps(), 5);
        assert_eq!(matrix.n_samples(), 3);

        // Verify genotypes match EIGENSTRAT fixture
        // EIGENSTRAT .geno row 1: 012 → IND1=HomRef, IND2=Het, IND3=HomAlt
        // VCF row 1: IND1=0/0 (HomRef), IND2=0/1 (Het), IND3=1/1 (HomAlt)
        assert_eq!(matrix.get(0, 0), Genotype::HomRef); // IND1, rs1
        assert_eq!(matrix.get(0, 1), Genotype::Het); // IND2, rs1
        assert_eq!(matrix.get(0, 2), Genotype::HomAlt); // IND3, rs1

        // EIGENSTRAT .geno row 2: 120 → IND1=Het, IND2=HomAlt, IND3=HomRef
        // VCF row 2: IND1=0/1 (Het), IND2=1/1 (HomAlt), IND3=0/0 (HomRef)
        assert_eq!(matrix.get(1, 0), Genotype::Het); // IND1, rs2
        assert_eq!(matrix.get(1, 1), Genotype::HomAlt); // IND2, rs2
        assert_eq!(matrix.get(1, 2), Genotype::HomRef); // IND3, rs2
    }

    #[test]
    fn test_vcf_reader_matches_eigenstrat() {
        // Read VCF
        let vcf_reader = VcfReader::new(VCF_FIXTURE).unwrap();
        let vcf_matrix = vcf_reader.read_all().unwrap();

        // Read EIGENSTRAT
        let eigenstrat_reader =
            crate::eigenstrat::EigenstratReader::new(EIGENSTRAT_FIXTURE).unwrap();
        let eigenstrat_matrix = eigenstrat_reader.read_all().unwrap();

        // Verify dimensions match
        assert_eq!(vcf_matrix.n_snps(), eigenstrat_matrix.n_snps());
        assert_eq!(vcf_matrix.n_samples(), eigenstrat_matrix.n_samples());

        // A genotype value counts copies of the allele listed first (EIGENSTRAT
        // .snp column 5); compare copies of the same named allele, so a
        // reversed allele order cannot pass. The fixture's .snp lists each
        // VCF ALT allele first, as convertf does for an imported VCF.
        let copies = |snp: &popmix_core::genotypes::Snp, g: Genotype, allele: char| match g {
            Genotype::Missing => None,
            _ if snp.allele_ref == allele => Some(g as u8),
            _ if snp.allele_alt == allele => Some(2 - g as u8),
            _ => panic!("{} has no allele {allele}", snp.id),
        };
        for snp_idx in 0..vcf_matrix.n_snps() {
            for sample_idx in 0..vcf_matrix.n_samples() {
                let allele = eigenstrat_matrix.snps()[snp_idx].allele_ref;
                let vcf = copies(
                    &vcf_matrix.snps()[snp_idx],
                    vcf_matrix.get(snp_idx, sample_idx),
                    allele,
                );
                let eig = copies(
                    &eigenstrat_matrix.snps()[snp_idx],
                    eigenstrat_matrix.get(snp_idx, sample_idx),
                    allele,
                );
                assert_eq!(
                    vcf, eig,
                    "copies of {allele} at SNP {snp_idx}, sample {sample_idx}"
                );
            }
        }
    }

    #[test]
    fn test_vcf_reader_missing_data() {
        let reader = VcfReader::new(VCF_MISSING_FIXTURE).unwrap();
        let matrix = reader.read_all().unwrap();

        // Check missing genotypes
        assert_eq!(matrix.get(0, 0), Genotype::Missing); // IND1, rs1
        assert_eq!(matrix.get(1, 1), Genotype::Missing); // IND2, rs2
        assert_eq!(matrix.get(2, 2), Genotype::Missing); // IND3, rs3
    }

    #[test]
    fn test_vcf_reader_phased() {
        let reader = VcfReader::new(VCF_PHASED_FIXTURE).unwrap();
        let matrix = reader.read_all().unwrap();

        // Phased genotypes should be parsed correctly (same data as unphased, with | instead of /)
        assert_eq!(matrix.get(0, 0), Genotype::HomRef); // 0|0
        assert_eq!(matrix.get(0, 1), Genotype::Het); // 0|1
        assert_eq!(matrix.get(0, 2), Genotype::HomAlt); // 1|1
    }
}
