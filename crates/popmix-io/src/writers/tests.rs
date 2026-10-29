use super::*;
use ndarray::array;
use popmix_core::genotypes::{Chromosome, Sex};
use tempfile::TempDir;

#[test]
fn ped_and_tgeno_written_by_sample_blocks_are_the_matrix_writes() {
    let family = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/golden/dates/family_chr20_22"
    );
    let geno = crate::PackedAncestryMapReader::new(family)
        .unwrap()
        .read_all()
        .unwrap();
    let file =
        crate::SnpMajor::open_packed(format!("{family}.geno"), geno.n_samples(), geno.n_snps())
            .unwrap();
    let dir = TempDir::new().unwrap();
    let p = |name: &str| dir.path().join(name).to_string_lossy().to_string();
    write_ped(p("matrix"), &geno).unwrap();
    write_tgeno(p("matrix_t"), &geno).unwrap();
    // 50 samples in blocks of 7: the last block is short
    for block in [7, 1000] {
        let (ped, tgeno) = (p(&format!("file{block}")), p(&format!("file_t{block}")));
        write_ped_blocks(&ped, geno.samples(), geno.snps(), &file, block).unwrap();
        write_tgeno_blocks(&tgeno, geno.samples(), geno.snps(), &file, block).unwrap();
        for ext in ["ped", "map", "pedind"] {
            assert_eq!(
                std::fs::read(format!("{}.{ext}", p("matrix"))).unwrap(),
                std::fs::read(format!("{ped}.{ext}")).unwrap(),
                "{ext}, blocks of {block}"
            );
        }
        for ext in ["geno", "snp", "ind"] {
            assert_eq!(
                std::fs::read(format!("{}.{ext}", p("matrix_t"))).unwrap(),
                std::fs::read(format!("{tgeno}.{ext}")).unwrap(),
                "{ext}, blocks of {block}"
            );
        }
    }
}

#[test]
fn streamed_writes_are_the_matrix_writes() {
    let family = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/golden/dates/family_chr20_22"
    );
    let reader = crate::PackedAncestryMapReader::new(family).unwrap();
    let geno = reader.read_all().unwrap();
    let (samples, snps) = (geno.samples().to_vec(), geno.snps().to_vec());
    let dir = TempDir::new().unwrap();
    type Whole = fn(&Path, &GenotypeMatrix) -> Result<()>;
    type Rows = fn(&Path, &[Sample], &[Snp], &mut dyn GenotypeRows) -> Result<()>;
    let formats: [(&str, Whole, Rows, &[&str]); 4] = [
        (
            "eigenstrat",
            |p, m| write_eigenstrat(p, m),
            |p, a, b, r| write_eigenstrat_rows(p, a, b, r),
            &["geno", "snp", "ind"],
        ),
        (
            "ancestrymap",
            |p, m| write_ancestrymap(p, m),
            |p, a, b, r| write_ancestrymap_rows(p, a, b, r),
            &["ancestrymapgeno", "snp", "ind"],
        ),
        (
            "packed",
            |p, m| write_packedancestrymap(p, m),
            |p, a, b, r| write_packedancestrymap_rows(p, a, b, r),
            &["geno", "snp", "ind"],
        ),
        (
            "bed",
            |p, m| write_packedped(p, m),
            |p, a, b, r| write_packedped_rows(p, a, b, r),
            &["bed", "bim", "fam"],
        ),
    ];
    for (name, whole, rows, exts) in formats {
        let (a, b) = (
            dir.path().join(format!("{name}_whole")),
            dir.path().join(format!("{name}_rows")),
        );
        whole(&a, &geno).unwrap();
        let mut source = crate::SnpMajorRows::open(
            format!("{family}.geno"),
            samples.len(),
            vec![true; snps.len()],
            (0..samples.len()).collect(),
        )
        .unwrap();
        rows(&b, &samples, &snps, &mut source).unwrap();
        for ext in exts {
            let read = |p: &Path| std::fs::read(format!("{}.{ext}", p.display())).unwrap();
            assert!(read(&a) == read(&b), "{name}: .{ext} differs");
        }
    }
}

fn make_test_samples() -> Vec<Sample> {
    vec![
        Sample::new("IND1", "PopA"),
        Sample::new("IND2", "PopA"),
        Sample::new("IND3", "PopB"),
    ]
}

fn make_test_snps() -> Vec<Snp> {
    vec![
        Snp::new("rs1", Chromosome::Autosome(1), 100, 0.0001, 'A', 'G'),
        Snp::new("rs2", Chromosome::Autosome(1), 200, 0.0002, 'C', 'T'),
        Snp::new("rs3", Chromosome::Autosome(2), 300, 0.0003, 'A', 'T'),
    ]
}

fn make_test_matrix() -> GenotypeMatrix {
    // 3 SNPs x 3 samples
    // SNP0: 0, 1, 2
    // SNP1: 1, 2, 0
    // SNP2: 2, 0, 1
    let data = array![[0, 1, 2], [1, 2, 0], [2, 0, 1]];
    GenotypeMatrix::new(data, make_test_samples(), make_test_snps()).unwrap()
}

#[test]
fn test_eigenstrat_roundtrip() {
    use crate::eigenstrat::EigenstratReader;

    let tmp = TempDir::new().unwrap();
    let prefix = tmp.path().join("test");
    let matrix = make_test_matrix();

    // Write
    write_eigenstrat(&prefix, &matrix).unwrap();

    // Read back
    let reader = EigenstratReader::new(&prefix).unwrap();
    let reloaded = reader.read_all().unwrap();

    // Compare
    assert_eq!(reloaded.n_samples(), matrix.n_samples());
    assert_eq!(reloaded.n_snps(), matrix.n_snps());

    for snp_idx in 0..matrix.n_snps() {
        for sample_idx in 0..matrix.n_samples() {
            assert_eq!(
                reloaded.get(snp_idx, sample_idx),
                matrix.get(snp_idx, sample_idx),
                "Mismatch at SNP {snp_idx} sample {sample_idx}"
            );
        }
    }
}

#[test]
fn test_packedancestrymap_roundtrip() {
    use crate::packedancestrymap::PackedAncestryMapReader;

    let tmp = TempDir::new().unwrap();
    let prefix = tmp.path().join("test");
    let matrix = make_test_matrix();

    // Write
    write_packedancestrymap(&prefix, &matrix).unwrap();

    // Read back
    let reader = PackedAncestryMapReader::new(&prefix).unwrap();
    let reloaded = reader.read_all().unwrap();

    // Compare
    assert_eq!(reloaded.n_samples(), matrix.n_samples());
    assert_eq!(reloaded.n_snps(), matrix.n_snps());

    for snp_idx in 0..matrix.n_snps() {
        for sample_idx in 0..matrix.n_samples() {
            assert_eq!(
                reloaded.get(snp_idx, sample_idx),
                matrix.get(snp_idx, sample_idx),
                "Mismatch at SNP {snp_idx} sample {sample_idx}"
            );
        }
    }
}

#[test]
fn test_packedped_roundtrip() {
    use crate::plink::PlinkReader;

    let tmp = TempDir::new().unwrap();
    let prefix = tmp.path().join("test");
    let matrix = make_test_matrix();

    // Write
    write_packedped(&prefix, &matrix).unwrap();

    // Read back
    let reader = PlinkReader::new(&prefix).unwrap();
    let reloaded = reader.read_all().unwrap();

    // Compare
    assert_eq!(reloaded.n_samples(), matrix.n_samples());
    assert_eq!(reloaded.n_snps(), matrix.n_snps());

    for snp_idx in 0..matrix.n_snps() {
        for sample_idx in 0..matrix.n_samples() {
            assert_eq!(
                reloaded.get(snp_idx, sample_idx),
                matrix.get(snp_idx, sample_idx),
                "Mismatch at SNP {snp_idx} sample {sample_idx}"
            );
        }
    }
}

fn sexed_matrix() -> GenotypeMatrix {
    let samples = vec![
        Sample::new("IND1", "PopA").with_sex(Sex::Male),
        Sample::new("IND2", "PopA").with_sex(Sex::Female),
        Sample::new("IND3", "PopB"),
    ];
    let data = array![[0, 1, 2], [1, 2, 0], [2, 0, 1]];
    GenotypeMatrix::new(data, samples, make_test_snps()).unwrap()
}

fn sexes(matrix: &GenotypeMatrix) -> Vec<Sex> {
    matrix.samples().iter().map(|s| s.sex).collect()
}

#[test]
fn test_sex_survives_ind_fam_and_ped() {
    use crate::eigenstrat::EigenstratReader;
    use crate::plink::PlinkReader;

    let tmp = TempDir::new().unwrap();
    let want = vec![Sex::Male, Sex::Female, Sex::Unknown];

    let eigen = tmp.path().join("eigen");
    write_eigenstrat(&eigen, &sexed_matrix()).unwrap();
    let ind = std::fs::read_to_string(eigen.with_extension("ind")).unwrap();
    assert_eq!(
        ind,
        format!(
            "{:>20} M {:>10}\n{:>20} F {:>10}\n{:>20} U {:>10}\n",
            "IND1", "PopA", "IND2", "PopA", "IND3", "PopB"
        )
    );
    let from_ind = EigenstratReader::new(&eigen).unwrap().read_all().unwrap();
    assert_eq!(sexes(&from_ind), want);

    let plink = tmp.path().join("plink");
    write_packedped(&plink, &from_ind).unwrap();
    let fam = std::fs::read_to_string(plink.with_extension("fam")).unwrap();
    let fam_sex: Vec<&str> = fam.lines().map(|l| l.split('\t').nth(4).unwrap()).collect();
    assert_eq!(fam_sex, ["1", "2", "0"]);
    let from_fam = PlinkReader::new(&plink).unwrap().read_all().unwrap();
    assert_eq!(sexes(&from_fam), want);

    let ped = tmp.path().join("ped");
    write_ped(&ped, &from_fam).unwrap();
    let from_ped = crate::ped::PedReader::new(&ped)
        .unwrap()
        .read_all()
        .unwrap();
    assert_eq!(sexes(&from_ped), want);
}

#[test]
fn test_eigenstrat_to_packedped_conversion() {
    use crate::eigenstrat::EigenstratReader;
    use crate::plink::PlinkReader;

    let tmp = TempDir::new().unwrap();
    let eigen_prefix = tmp.path().join("eigen");
    let plink_prefix = tmp.path().join("plink");
    let matrix = make_test_matrix();

    // Write EIGENSTRAT
    write_eigenstrat(&eigen_prefix, &matrix).unwrap();

    // Read EIGENSTRAT
    let reader = EigenstratReader::new(&eigen_prefix).unwrap();
    let eigen_matrix = reader.read_all().unwrap();

    // Write PACKEDPED
    write_packedped(&plink_prefix, &eigen_matrix).unwrap();

    // Read PACKEDPED
    let reader = PlinkReader::new(&plink_prefix).unwrap();
    let plink_matrix = reader.read_all().unwrap();

    // Compare genotypes
    for snp_idx in 0..matrix.n_snps() {
        for sample_idx in 0..matrix.n_samples() {
            assert_eq!(
                plink_matrix.get(snp_idx, sample_idx),
                matrix.get(snp_idx, sample_idx),
                "Mismatch at SNP {snp_idx} sample {sample_idx}"
            );
        }
    }
}

#[test]
fn test_packedped_to_packedancestrymap_conversion() {
    use crate::packedancestrymap::PackedAncestryMapReader;
    use crate::plink::PlinkReader;

    let tmp = TempDir::new().unwrap();
    let plink_prefix = tmp.path().join("plink");
    let packed_prefix = tmp.path().join("packed");
    let matrix = make_test_matrix();

    // Write PACKEDPED
    write_packedped(&plink_prefix, &matrix).unwrap();

    // Read PACKEDPED
    let reader = PlinkReader::new(&plink_prefix).unwrap();
    let plink_matrix = reader.read_all().unwrap();

    // Write PACKEDANCESTRYMAP
    write_packedancestrymap(&packed_prefix, &plink_matrix).unwrap();

    // Read PACKEDANCESTRYMAP
    let reader = PackedAncestryMapReader::new(&packed_prefix).unwrap();
    let packed_matrix = reader.read_all().unwrap();

    // Compare genotypes
    for snp_idx in 0..matrix.n_snps() {
        for sample_idx in 0..matrix.n_samples() {
            assert_eq!(
                packed_matrix.get(snp_idx, sample_idx),
                matrix.get(snp_idx, sample_idx),
                "Mismatch at SNP {snp_idx} sample {sample_idx}"
            );
        }
    }
}

#[test]
fn test_with_missing_genotypes() {
    use crate::eigenstrat::EigenstratReader;

    let tmp = TempDir::new().unwrap();
    let prefix = tmp.path().join("test");

    // Include missing genotypes (value 3)
    let data = array![[0, 1, 3], [3, 2, 0], [2, 3, 1]];
    let matrix = GenotypeMatrix::new(data, make_test_samples(), make_test_snps()).unwrap();

    // Write and read back
    write_eigenstrat(&prefix, &matrix).unwrap();
    let reader = EigenstratReader::new(&prefix).unwrap();
    let reloaded = reader.read_all().unwrap();

    // Verify missing genotypes are preserved
    assert_eq!(reloaded.get(0, 2), Genotype::Missing);
    assert_eq!(reloaded.get(1, 0), Genotype::Missing);
    assert_eq!(reloaded.get(2, 1), Genotype::Missing);
}

#[test]
fn test_ancestrymap_write() {
    let tmp = TempDir::new().unwrap();
    let prefix = tmp.path().join("test");

    // Include missing genotypes
    let data = array![[0, 1, 3], [3, 2, 0], [2, 3, 1]];
    let matrix = GenotypeMatrix::new(data, make_test_samples(), make_test_snps()).unwrap();

    // Write ANCESTRYMAP
    write_ancestrymap(&prefix, &matrix).unwrap();

    // Every genotype, the EIGENSTRAT code, -1 missing, as convertf writes it
    let geno_content =
        std::fs::read_to_string(format!("{}.ancestrymapgeno", prefix.display())).unwrap();
    let lines: Vec<&str> = geno_content.lines().collect();
    assert_eq!(lines.len(), 9);
    let codes: Vec<&str> = lines
        .iter()
        .map(|l| l.split_whitespace().nth(2).unwrap())
        .collect();
    assert_eq!(codes, ["0", "1", "-1", "-1", "2", "0", "2", "-1", "1"]);
}

#[test]
fn test_ped_write() {
    let tmp = TempDir::new().unwrap();
    let prefix = tmp.path().join("test");
    let matrix = make_test_matrix();

    // Write PED
    write_ped(&prefix, &matrix).unwrap();

    // Verify .ped file exists and has correct structure
    let ped_content = std::fs::read_to_string(format!("{}.ped", prefix.display())).unwrap();
    let lines: Vec<&str> = ped_content.lines().collect();

    // Should have 3 lines (one per sample)
    assert_eq!(lines.len(), 3);

    for line in &lines {
        let parts: Vec<&str> = line.split('\t').collect();
        // 6 info columns + 2 allele columns per SNP (3 SNPs = 6 allele columns)
        assert_eq!(parts.len(), 6 + 3 * 2);
    }
    // the first sample's codes 0, 1, 2 count copies of each SNP's first
    // allele (A/G, C/T, A/T)
    let first: Vec<&str> = lines[0].split('\t').skip(6).collect();
    assert_eq!(first, ["G", "G", "C", "T", "A", "A"]);

    // Verify .map file
    let map_content = std::fs::read_to_string(format!("{}.map", prefix.display())).unwrap();
    let map_lines: Vec<&str> = map_content.lines().collect();
    assert_eq!(map_lines.len(), 3);
}

/// .snp and .ind as AdmixTools convertf writes them (AADR's own layout):
/// right-aligned columns, X and Y as 23 and 24, a zero genetic position as 0.
#[test]
fn snp_and_ind_lines_are_laid_out_as_admixtools_writes_them() {
    let snps = vec![
        Snp::new(
            "rs3094315",
            Chromosome::Autosome(1),
            752566,
            0.020130,
            'G',
            'A',
        ),
        Snp::new("rs2562132", Chromosome::Autosome(16), 84170, 0.0, 'C', 'T'),
        Snp::new("rs4129148", Chromosome::X, 990180, 0.071651, 'C', 'G'),
        Snp::new("M288", Chromosome::Y, 2649694, 0.0, 'G', 'T'),
    ];
    let samples = vec![
        Sample::new("Loschbour.AG", "Luxembourg_Loschbour_Mesolithic"),
        Sample::new("I1838.AG", "Spain_C"),
    ];
    let matrix =
        GenotypeMatrix::new(array![[0, 1], [1, 1], [2, 1], [0, 1]], samples, snps).unwrap();
    let tmp = TempDir::new().unwrap();
    let prefix = tmp.path().join("out");
    write_eigenstrat(&prefix, &matrix).unwrap();
    let snp = std::fs::read_to_string(prefix.with_extension("snp")).unwrap();
    assert_eq!(
        snp,
        "           rs3094315     1        0.020130          752566 G A\n\
         \x20          rs2562132    16               0           84170 C T\n\
         \x20          rs4129148    23        0.071651          990180 C G\n\
         \x20               M288    24               0         2649694 G T\n"
    );
    let ind = std::fs::read_to_string(prefix.with_extension("ind")).unwrap();
    assert_eq!(
        ind,
        "        Loschbour.AG U Luxembourg_Loschbour_Mesolithic\n            I1838.AG U    Spain_C\n"
    );
}
