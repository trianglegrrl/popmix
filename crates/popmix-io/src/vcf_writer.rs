//! VCF output (TASKS.md P1): one record per SNP, GT for every sample, as
//! PLINK 2's `--export vcf` writes it from the same genotypes
//! (`tests/golden/formats/vcf_out`).
//!
//! A genotype counts copies of the SNP's first allele (EIGENSTRAT `.snp`
//! column 5, PLINK A1), which is VCF's ALT; the second is REF. On
//! chromosome X a male's homozygous call is haploid (`0` or `1`) and a
//! heterozygous one stays `0/1`, as PLINK 2 writes them. An unknown allele
//! (`0`, `X`, `.`) is `.` as ALT and `N` as REF, where PLINK 2 writes the
//! code as it is, which the VCF specification does not allow.

use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use popmix_core::genotypes::{Chromosome, Sex};
use popmix_core::{GenotypeRows, IoResultExt, Result, Sample, Snp};

/// The VCF file `output` names: itself when it ends in .vcf or .vcf.gz
/// (gzip-compressed), or `output.vcf`.
// the comparisons are of a lower-cased name
#[allow(clippy::case_sensitive_file_extension_comparisons)]
pub fn vcf_output_path(output: &str) -> PathBuf {
    let lower = output.to_lowercase();
    if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
        PathBuf::from(output)
    } else {
        PathBuf::from(format!("{output}.vcf"))
    }
}

fn is_unknown(allele: char) -> bool {
    matches!(allele, '0' | 'X' | 'x' | '.' | 'N' | 'n' | '?')
}

/// The GT of genotype code `g` (copies of ALT; 3 missing).
fn gt(g: u8, haploid: bool) -> &'static str {
    match (g, haploid) {
        (0, false) => "0/0",
        (1, _) => "0/1",
        (2, false) => "1/1",
        (0, true) => "0",
        (2, true) => "1",
        (_, false) => "./.",
        (_, true) => ".",
    }
}

/// Write `rows` (one per SNP of `snps`, every sample) as the VCF `output`
/// names (see [`vcf_output_path`]).
pub fn write_vcf_rows(
    output: &str,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn GenotypeRows,
) -> Result<()> {
    let path = vcf_output_path(output);
    let file = File::create(&path).with_path(&path)?;
    let mut out: Box<dyn Write> = if path.to_string_lossy().to_lowercase().ends_with(".gz") {
        Box::new(BufWriter::new(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::default(),
        )))
    } else {
        Box::new(BufWriter::new(file))
    };
    write_vcf_to(&mut out, samples, snps, rows).with_path(&path)?;
    out.flush().with_path(&path)?;
    Ok(())
}

fn write_vcf_to(
    out: &mut dyn Write,
    samples: &[Sample],
    snps: &[Snp],
    rows: &mut dyn GenotypeRows,
) -> std::io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.3")?;
    writeln!(out, "##source=popmix {}", env!("CARGO_PKG_VERSION"))?;
    let mut contigs: Vec<String> = Vec::new();
    for snp in snps {
        let c = snp.chromosome.to_string();
        if !contigs.contains(&c) {
            contigs.push(c);
        }
    }
    for c in &contigs {
        writeln!(out, "##contig=<ID={c}>")?;
    }
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for s in samples {
        write!(out, "\t{}", s.id)?;
    }
    writeln!(out)?;
    let male: Vec<bool> = samples.iter().map(|s| s.sex == Sex::Male).collect();
    let mut row = vec![0u8; samples.len()];
    let mut line = String::new();
    for snp in snps {
        rows.next_row(&mut row).map_err(std::io::Error::other)?;
        let alt = if is_unknown(snp.allele_ref) {
            '.'
        } else {
            snp.allele_ref
        };
        let reference = if is_unknown(snp.allele_alt) {
            'N'
        } else {
            snp.allele_alt
        };
        let on_x = snp.chromosome == Chromosome::X;
        line.clear();
        let _ = write!(
            line,
            "{}\t{}\t{}\t{reference}\t{alt}\t.\t.\t.\tGT",
            snp.chromosome, snp.position_bp, snp.id
        );
        for (&g, &m) in row.iter().zip(&male) {
            line.push('\t');
            line.push_str(gt(g, on_x && m));
        }
        line.push('\n');
        out.write_all(line.as_bytes())?;
    }
    Ok(())
}
