use super::*;
use crate::eigenstrat::{read_ind_file, read_snp_file};
use crate::PackedAncestryMapReader;

const SMALL: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/golden/dates/family_chr20_22"
);

/// 54 samples: each record is padded to 48 bytes, not ceil(54 / 4) = 14.
#[test]
fn small_files_stream_with_48_byte_records() {
    let geno_path = format!("{SMALL}.geno");
    let mut samples = read_ind_file(&format!("{SMALL}.ind")).unwrap();
    let snps = read_snp_file(&format!("{SMALL}.snp")).unwrap();
    assert_eq!(samples.len(), 54);
    detect_ploidy_packedancestrymap(&geno_path, &mut samples, 1000, snps.len()).unwrap();
    let pops = ["French", "Yoruba"];
    let streamed =
        compute_allele_frequencies_packedancestrymap(&geno_path, &samples, &snps, &pops, false)
            .unwrap();

    let geno = PackedAncestryMapReader::new(SMALL)
        .unwrap()
        .read_all()
        .unwrap();
    let in_memory = AlleleFrequencies::from_genotypes_for_populations(&geno, Some(&pops));
    assert_eq!(streamed.n_snps(), in_memory.n_snps());
    for (col, pop) in pops.iter().enumerate() {
        let other = in_memory.pop_index(pop).unwrap();
        for snp in 0..streamed.n_snps() {
            let (a, b) = (
                streamed.frequencies()[[snp, col]],
                in_memory.frequencies()[[snp, other]],
            );
            assert!(
                a == b || (a.is_nan() && b.is_nan()),
                "{pop} SNP {snp}: {a} != {b}"
            );
        }
    }
}

#[test]
fn frequencies_do_not_depend_on_the_thread_count() {
    let geno_path = format!("{SMALL}.geno");
    let mut samples = read_ind_file(&format!("{SMALL}.ind")).unwrap();
    let snps = read_snp_file(&format!("{SMALL}.snp")).unwrap();
    detect_ploidy_packedancestrymap(&geno_path, &mut samples, 1000, snps.len()).unwrap();
    let pops = ["Yoruba", "French", "Han"];
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                compute_allele_frequencies_packedancestrymap(
                    &geno_path, &samples, &snps, &pops, true,
                )
                .unwrap()
            })
    };
    let one = run(1);
    for threads in [2, 7] {
        let many = run(threads);
        assert_eq!(one.snps().len(), many.snps().len());
        assert_eq!(one.counts(), many.counts());
        assert!(one
            .frequencies()
            .iter()
            .zip(many.frequencies())
            .all(|(a, b)| a.to_bits() == b.to_bits()));
    }
}

/// The family data as a PLINK .bed in a fresh directory, and its prefix.
fn family_bed() -> (tempfile::TempDir, String) {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("family").to_string_lossy().to_string();
    let geno = PackedAncestryMapReader::new(SMALL)
        .unwrap()
        .read_all()
        .unwrap();
    crate::writers::write_packedped(&prefix, &geno).unwrap();
    (dir, prefix)
}

#[test]
fn readers_are_shared_between_threads() {
    use rayon::prelude::*;
    fn sync<T: Sync>(_: &T) {}
    let (_dir, bed) = family_bed();
    let plink = crate::PlinkReader::new(&bed).unwrap();
    let packed = PackedAncestryMapReader::new(SMALL).unwrap();
    sync(&plink);
    sync(&packed);
    let n = plink.num_snps();
    let parts: Vec<(Vec<u8>, Vec<u8>)> = (0..4)
        .into_par_iter()
        .map(|k| {
            let r = k * n / 4..(k + 1) * n / 4;
            let a = plink.read_snp_range(r.clone()).unwrap();
            let b = packed.read_snp_range(r).unwrap();
            (a.snp_raw(0).to_vec(), b.snp_raw(0).to_vec())
        })
        .collect();
    for (a, b) in parts {
        assert_eq!(a, b);
    }
}

#[test]
fn bed_rows_match_the_in_memory_genotypes() {
    use popmix_core::GenotypeRows;
    let (_dir, bed) = family_bed();
    let geno = crate::PlinkReader::new(&bed).unwrap().read_all().unwrap();
    let keep: Vec<bool> = (0..geno.n_snps()).map(|i| i % 3 != 1).collect();
    let samples = vec![53, 0, 7, 20];
    let mut rows =
        SnpMajorRows::open_bed(&format!("{bed}.bed"), 54, keep.clone(), samples.clone()).unwrap();
    let mut row = vec![0u8; samples.len()];
    for snp in (0..geno.n_snps()).filter(|&i| keep[i]) {
        rows.next_row(&mut row).unwrap();
        let want: Vec<u8> = samples.iter().map(|&s| geno.snp_raw(snp)[s]).collect();
        assert_eq!(row, want, "SNP {snp}");
    }
    assert!(rows.next_row(&mut row).is_err());
}

#[test]
fn batches_of_rows_are_the_rows_one_at_a_time() {
    use popmix_core::GenotypeRows;
    let (_dir, bed) = family_bed();
    let n_snps = crate::PlinkReader::new(&bed).unwrap().num_snps();
    let keep: Vec<bool> = (0..n_snps).map(|i| i % 5 != 2).collect();
    let samples = vec![53, 0, 7, 20, 31];
    let open = || {
        SnpMajorRows::open_bed(&format!("{bed}.bed"), 54, keep.clone(), samples.clone()).unwrap()
    };
    let (mut one, mut many) = (open(), open());
    let kept = keep.iter().filter(|&&k| k).count();
    let mut batch = vec![0u8; 1000 * samples.len()];
    let mut row = vec![0u8; samples.len()];
    let mut done = 0;
    while done < kept {
        let n = (kept - done).min(1000);
        many.next_rows(&mut batch[..n * samples.len()], samples.len())
            .unwrap();
        for r in 0..n {
            one.next_row(&mut row).unwrap();
            assert_eq!(&batch[r * samples.len()..(r + 1) * samples.len()], &row[..]);
        }
        done += n;
    }
    assert!(many
        .next_rows(&mut batch[..samples.len()], samples.len())
        .is_err());
}

#[test]
fn bed_frequencies_and_ploidy_match_the_in_memory_ones() {
    let (_dir, bed) = family_bed();
    let mut geno = crate::PlinkReader::new(&bed).unwrap().read_all().unwrap();
    geno.detect_ploidy(geno.n_snps());
    let mut samples = geno.samples().to_vec();
    samples.iter_mut().for_each(|s| s.ploidy = Ploidy::Diploid);
    let file = SnpMajor::open_bed(&format!("{bed}.bed"), 54, geno.n_snps()).unwrap();
    let all: Vec<usize> = (0..54).collect();
    file.detect_ploidy(&mut samples, &all, geno.n_snps());
    let pops = ["French", "Simulation", "Yoruba"];
    let ploidy = |s: &[Sample]| s.iter().map(|x| x.ploidy).collect::<Vec<_>>();
    assert_eq!(ploidy(&samples), ploidy(geno.samples()));
    let streamed = file.allele_frequencies(&samples, geno.snps(), &pops, false);
    let in_memory = AlleleFrequencies::from_genotypes_for_populations(&geno, Some(&pops));
    for (col, pop) in pops.iter().enumerate() {
        let other = in_memory.pop_index(pop).unwrap();
        assert_eq!(
            streamed.counts().column(col),
            in_memory.counts().column(other)
        );
        let (a, b) = (
            streamed.frequencies().column(col),
            in_memory.frequencies().column(other),
        );
        assert!(a
            .iter()
            .zip(b)
            .all(|(x, y)| x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan())));
    }
}

/// The family data as text EIGENSTRAT in a fresh directory, and its prefix.
fn family_text() -> (tempfile::TempDir, String) {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("family").to_string_lossy().to_string();
    let geno = PackedAncestryMapReader::new(SMALL)
        .unwrap()
        .read_all()
        .unwrap();
    crate::writers::write_eigenstrat(&prefix, &geno).unwrap();
    (dir, prefix)
}

#[test]
fn text_geno_is_read_in_place() {
    use popmix_core::GenotypeRows;
    let (_dir, text) = family_text();
    let mut geno = crate::EigenstratReader::new(&text)
        .unwrap()
        .read_all()
        .unwrap();
    let path = format!("{text}.geno");
    let file = SnpMajor::open_text(&path, 54, geno.n_snps())
        .unwrap()
        .unwrap();
    let keep: Vec<bool> = (0..geno.n_snps()).map(|i| i % 3 != 1).collect();
    let samples = vec![53, 0, 7, 20];
    let mut rows = SnpMajorRows::from_file(
        SnpMajor::open_text(&path, 54, geno.n_snps())
            .unwrap()
            .unwrap(),
        keep.clone(),
        samples.clone(),
    )
    .unwrap();
    let mut row = vec![0u8; samples.len()];
    for snp in (0..geno.n_snps()).filter(|&i| keep[i]) {
        rows.next_row(&mut row).unwrap();
        let want: Vec<u8> = samples.iter().map(|&s| geno.snp_raw(snp)[s]).collect();
        assert_eq!(row, want, "SNP {snp}");
    }

    geno.detect_ploidy(geno.n_snps());
    let mut samples = geno.samples().to_vec();
    samples.iter_mut().for_each(|s| s.ploidy = Ploidy::Diploid);
    let all: Vec<usize> = (0..54).collect();
    file.detect_ploidy(&mut samples, &all, geno.n_snps());
    let ploidy = |s: &[Sample]| s.iter().map(|x| x.ploidy).collect::<Vec<_>>();
    assert_eq!(ploidy(&samples), ploidy(geno.samples()));
    let pops = ["French", "Yoruba"];
    let streamed = file.allele_frequencies(&samples, geno.snps(), &pops, false);
    let in_memory = AlleleFrequencies::from_genotypes_for_populations(&geno, Some(&pops));
    for (col, pop) in pops.iter().enumerate() {
        let other = in_memory.pop_index(pop).unwrap();
        assert_eq!(
            streamed.counts().column(col),
            in_memory.counts().column(other)
        );
    }
}

#[test]
fn text_geno_with_uneven_lines_is_left_to_the_reader() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("x.geno");
    // CRLF line ends: the right size for neither layout
    std::fs::write(&path, b"012\r\n929\r\n").unwrap();
    assert!(SnpMajor::open_text(path.to_str().unwrap(), 3, 2)
        .unwrap()
        .is_none());
    // the right size, but a line one short and the next one long
    std::fs::write(&path, b"01\n2929\n").unwrap();
    assert!(SnpMajor::open_text(path.to_str().unwrap(), 3, 2)
        .unwrap()
        .is_none());
    std::fs::write(&path, b"012\n929\n").unwrap();
    let file = SnpMajor::open_text(path.to_str().unwrap(), 3, 2)
        .unwrap()
        .unwrap();
    let r = file.record(1);
    assert_eq!(
        (0..3).map(|s| file.code(r, s)).collect::<Vec<_>>(),
        vec![3, 2, 3]
    );
}

#[test]
fn a_truncated_file_is_an_error() {
    let dir = std::env::temp_dir().join(format!("popmix-streaming-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let geno_path = dir.join("short.geno");
    let full = std::fs::read(format!("{SMALL}.geno")).unwrap();
    std::fs::write(&geno_path, &full[..full.len() - 10]).unwrap();
    let samples = read_ind_file(&format!("{SMALL}.ind")).unwrap();
    let snps = read_snp_file(&format!("{SMALL}.snp")).unwrap();
    let err = compute_allele_frequencies_packedancestrymap(
        geno_path.to_str().unwrap(),
        &samples,
        &snps,
        &["French"],
        false,
    );
    assert!(err.is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn packed_rows_match_the_in_memory_genotypes() {
    use popmix_core::GenotypeRows;
    let geno = PackedAncestryMapReader::new(SMALL)
        .unwrap()
        .read_all()
        .unwrap();
    let keep: Vec<bool> = (0..geno.n_snps()).map(|i| i % 3 != 1).collect();
    let samples = vec![53, 0, 7, 20];
    let mut rows =
        PackedRows::open(&format!("{SMALL}.geno"), 54, keep.clone(), samples.clone()).unwrap();
    let mut row = vec![0u8; samples.len()];
    for snp in (0..geno.n_snps()).filter(|&i| keep[i]) {
        rows.next_row(&mut row).unwrap();
        let want: Vec<u8> = samples.iter().map(|&s| geno.snp_raw(snp)[s]).collect();
        assert_eq!(row, want, "SNP {snp}");
    }
    assert!(rows.next_row(&mut row).is_err());
}

#[test]
fn selected_ploidy_matches_the_full_detection() {
    let geno_path = format!("{SMALL}.geno");
    let mut all = read_ind_file(&format!("{SMALL}.ind")).unwrap();
    let n_snps = read_snp_file(&format!("{SMALL}.snp")).unwrap().len();
    detect_ploidy_packedancestrymap(&geno_path, &mut all, n_snps, n_snps).unwrap();
    let mut some = read_ind_file(&format!("{SMALL}.ind")).unwrap();
    let selected = [0, 5, 53];
    detect_ploidy_packedancestrymap_selected(&geno_path, &mut some, &selected, n_snps, n_snps)
        .unwrap();
    for &s in &selected {
        assert_eq!(some[s].ploidy, all[s].ploidy);
    }
}
