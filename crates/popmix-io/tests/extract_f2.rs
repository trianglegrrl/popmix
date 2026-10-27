//! `popmix_io::extract::extract_f2`, the f2 pipeline the Python bindings
//! call (TASKS.md I4): every format gives the same blocks, samples are
//! relabelled as admixtools2's `inds`/`pops`, and bad requests are errors.

use std::path::Path;

use popmix_core::{warnings, F2Blocks, GenotypeMatrix, PopmixError};
use popmix_io::dataset::{read_genotypes, PloidyRule};
use popmix_io::extract::{extract_f2, ExtractOptions};
use popmix_io::{write_eigenstrat, write_packedped, write_ped, write_tgeno};

const POPS: [&str; 3] = ["French", "Yoruba", "Simulation.1"];

fn family() -> String {
    format!(
        "{}/../../tests/golden/dates/family_chr20_22",
        env!("CARGO_MANIFEST_DIR")
    )
}

fn options(pops: &[&str]) -> ExtractOptions {
    ExtractOptions::new(pops.iter().map(|p| p.to_string()).collect())
}

fn extract(prefix: &str, opts: &ExtractOptions) -> popmix_core::Result<F2Blocks> {
    extract_f2(prefix, opts, &mut warnings::ignore)
}

fn assert_same(got: &F2Blocks, want: &F2Blocks, what: &str) {
    assert_eq!(got.populations, want.populations, "{what}");
    assert_eq!(got.block_lengths, want.block_lengths, "{what}");
    assert_eq!(got.block_weights, want.block_weights, "{what}");
    for (g, w) in got.block_sums.iter().zip(&want.block_sums) {
        assert!((g - w).abs() <= 1e-12 * w.abs(), "{what}: {g} != {w}");
    }
}

#[test]
fn every_format_gives_the_same_f2_blocks() {
    let packed = extract(&family(), &options(&POPS)).unwrap();
    assert_eq!(packed.populations.len(), 3);
    let geno = read_genotypes(&family(), &mut warnings::ignore).unwrap();
    let dir = tempfile::tempdir().unwrap();
    type Write = fn(&Path, &GenotypeMatrix) -> popmix_core::Result<()>;
    let writers: [(&str, Write); 4] = [
        ("text EIGENSTRAT", |p, g| write_eigenstrat(p, g)),
        ("PLINK", |p, g| write_packedped(p, g)),
        ("TGENO", |p, g| write_tgeno(p, g)),
        ("PED", |p, g| write_ped(p, g)),
    ];
    for (name, write) in writers {
        let prefix = dir.path().join(name.replace(' ', "_"));
        write(&prefix, &geno).unwrap();
        let got = extract(prefix.to_str().unwrap(), &options(&POPS)).unwrap();
        assert_same(&got, &packed, name);
    }
}

#[test]
fn samples_take_the_populations_given_and_the_others_are_left_out() {
    let geno = read_genotypes(&family(), &mut warnings::ignore).unwrap();
    let french: Vec<String> = geno
        .samples()
        .iter()
        .filter(|s| s.population.as_str() == "French")
        .map(|s| s.id.clone())
        .collect();
    // the French samples split in two populations of their own
    let half = french.len() / 2;
    let labels = french
        .iter()
        .enumerate()
        .map(|(i, id)| (id.clone(), if i < half { "A" } else { "B" }.to_string()))
        .collect();
    let mut opts = options(&["A", "B"]);
    opts.labels = Some(labels);
    let got = extract(&family(), &opts).unwrap();
    assert_eq!(got.populations.len(), 2);
    // Yoruba is no longer there: only the listed samples are kept
    let mut opts = options(&["A", "Yoruba"]);
    opts.labels = Some([(french[0].clone(), "A".to_string())].into());
    let err = extract(&family(), &opts).unwrap_err();
    assert!(
        matches!(err, PopmixError::PopulationNotFound(ref p) if p == "Yoruba"),
        "{err}"
    );
}

#[test]
fn an_unknown_population_or_sample_is_an_error() {
    let err = extract(&family(), &options(&["French", "Nowhere"])).unwrap_err();
    assert!(
        matches!(err, PopmixError::PopulationNotFound(ref p) if p == "Nowhere"),
        "{err}"
    );
    let mut opts = options(&["A"]);
    opts.labels = Some([("no-such-sample".to_string(), "A".to_string())].into());
    let err = extract(&family(), &opts).unwrap_err();
    assert!(err.to_string().contains("no-such-sample"), "{err}");
}

#[test]
fn the_ploidy_rule_is_one_of_three_names() {
    assert_eq!("detect".parse::<PloidyRule>().unwrap(), PloidyRule::Genome);
    assert_eq!(
        "admixtools2".parse::<PloidyRule>().unwrap(),
        PloidyRule::First1000
    );
    assert_eq!(
        "diploid".parse::<PloidyRule>().unwrap(),
        PloidyRule::Diploid
    );
    let err = "haploid".parse::<PloidyRule>().unwrap_err();
    assert!(matches!(err, PopmixError::InvalidParameter { .. }), "{err}");
}
