use super::*;
use crate::{write_tgeno, PackedAncestryMapReader};
use popmix_core::GenotypeMatrix;

const SMALL: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/golden/dates/family_chr20_22"
);

/// The family data (50 samples, 21,175 SNPs) as TGENO, and in memory.
fn family_tgeno() -> (tempfile::TempDir, String, GenotypeMatrix) {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("t").display().to_string();
    let geno = PackedAncestryMapReader::new(SMALL)
        .unwrap()
        .read_all()
        .unwrap();
    write_tgeno(&prefix, &geno).unwrap();
    (dir, prefix, geno)
}

#[test]
fn rows_are_the_genotypes_whatever_the_chunk() {
    let (_dir, prefix, geno) = family_tgeno();
    let n_snps = geno.n_snps();
    // every third SNP, and samples out of order
    let keep: Vec<bool> = (0..n_snps).map(|i| i % 3 != 1).collect();
    let samples = vec![49, 0, 7, 7, 20];
    for chunk in [1, 3, 4, 1000, n_snps + 5] {
        let reader = TgenoReader::new(&prefix).unwrap();
        let mut rows = TgenoRows::with_chunk(reader, keep.clone(), samples.clone(), chunk).unwrap();
        let mut row = vec![0u8; samples.len()];
        for snp in (0..n_snps).filter(|&i| keep[i]) {
            rows.next_row(&mut row).unwrap();
            let want: Vec<u8> = samples.iter().map(|&s| geno.get(snp, s) as u8).collect();
            assert_eq!(row, want, "chunk {chunk}, SNP {snp}");
        }
        assert!(rows.next_row(&mut row).is_err(), "chunk {chunk}");
    }
}

#[test]
fn batches_of_rows_are_the_rows_one_at_a_time() {
    let (_dir, prefix, geno) = family_tgeno();
    let keep = vec![true; geno.n_snps()];
    let samples: Vec<usize> = (0..geno.n_samples()).collect();
    let width = samples.len();
    let reader = TgenoReader::new(&prefix).unwrap();
    let mut rows = TgenoRows::with_chunk(reader, keep, samples, 100).unwrap();
    let mut batch = vec![0u8; 250 * width];
    rows.next_rows(&mut batch, width).unwrap();
    for (snp, row) in batch.chunks(width).enumerate() {
        let want: Vec<u8> = (0..width).map(|s| geno.get(snp, s) as u8).collect();
        assert_eq!(row, want.as_slice(), "SNP {snp}");
    }
}

#[test]
fn blocks_of_samples_are_the_genotypes() {
    let (_dir, prefix, geno) = family_tgeno();
    let reader = TgenoReader::new(&prefix).unwrap();
    let n_snps = geno.n_snps();
    let mut out = vec![0u8; 3 * n_snps];
    reader.fill_samples(10..13, n_snps, &mut out);
    for (k, s) in (10..13).enumerate() {
        let want: Vec<u8> = (0..n_snps).map(|snp| geno.get(snp, s) as u8).collect();
        assert_eq!(
            &out[k * n_snps..(k + 1) * n_snps],
            want.as_slice(),
            "sample {s}"
        );
    }
}

#[test]
fn a_sample_out_of_range_is_an_error() {
    let (_dir, prefix, geno) = family_tgeno();
    let reader = TgenoReader::new(&prefix).unwrap();
    assert!(TgenoRows::new(reader, vec![true; geno.n_snps()], vec![50]).is_err());
}
