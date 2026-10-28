//! The provenance `extract_f2` records with its f2 blocks (TASKS.md H3):
//! the popmix version, the source, every extraction parameter, the SNP
//! counts and each sample's ploidy, saved in an f2 directory's
//! `metadata.json`.

use popmix_core::genotypes::{Ploidy, Sample};
use popmix_core::F2Blocks;
use serde_json::{json, Value};

use crate::extract::ExtractOptions;

/// The provenance of `blocks`, extracted from `prefix` with `opts`;
/// `samples` are the dataset's samples after relabelling and ploidy
/// detection, of which those in `opts.populations` are recorded.
pub(crate) fn extraction(
    prefix: &str,
    opts: &ExtractOptions,
    samples: &[Sample],
    n_snps_read: usize,
    blocks: &F2Blocks,
) -> Value {
    let samples: Vec<Value> = samples
        .iter()
        .filter(|s| opts.populations.iter().any(|p| p == s.population.as_str()))
        .map(|s| {
            json!({
                "id": s.id,
                "population": s.population.as_str(),
                "ploidy": match s.ploidy {
                    Ploidy::Diploid => "diploid",
                    Ploidy::Pseudohaploid => "pseudohaploid",
                },
            })
        })
        .collect();
    json!({
        "popmix_version": env!("CARGO_PKG_VERSION"),
        "source": prefix,
        "parameters": {
            "populations": opts.populations,
            "auto_only": opts.auto_only,
            "block_size": opts.blocks.block_size,
            "maxmiss": opts.f2.max_missing,
            "poly_only": opts.f2.poly_only,
            "pseudohaploid": opts.ploidy.name(),
            "relabelled_samples": opts.labels.as_ref().map_or(0, std::collections::HashMap::len),
        },
        "n_snps_read": n_snps_read,
        "n_snps_used": blocks.block_lengths().iter().sum::<usize>(),
        "samples": samples,
    })
}
