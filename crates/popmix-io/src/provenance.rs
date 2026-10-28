//! The provenance `extract_f2` records with its f2 blocks (TASKS.md H3):
//! the popmix version, the source, every extraction parameter, the SNP
//! counts and each sample's ploidy, saved in an f2 directory's
//! `metadata.json`; and the SNP panel it read (TASKS.md H9), as its size
//! and a fingerprint of its positions, since f-statistics depend on how the
//! SNPs were ascertained.

use popmix_core::genotypes::{Ploidy, Sample, Snp};
use popmix_core::F2Blocks;
use serde_json::{json, Value};

use crate::extract::ExtractOptions;

/// The provenance of `blocks`, extracted from `prefix` with `opts`;
/// `samples` are the dataset's samples after relabelling and ploidy
/// detection, of which those in `opts.populations` are recorded, and
/// `panel` the SNPs read.
pub(crate) fn extraction(
    prefix: &str,
    opts: &ExtractOptions,
    samples: &[Sample],
    panel: &Panel,
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
        "n_snps_read": panel.snps,
        "snp_panel": {"snps": panel.snps, "fingerprint": panel.fingerprint},
        "n_snps_used": blocks.block_lengths().iter().sum::<usize>(),
        "samples": samples,
    })
}

/// A SNP panel: its size, and a 64-bit FNV-1a of the SNPs' chromosomes and
/// positions, in order, as 16 hex digits (the same for two datasets on one
/// panel whatever their SNP names or samples).
pub(crate) struct Panel {
    snps: usize,
    fingerprint: String,
}

impl Panel {
    pub(crate) fn of(snps: &[Snp]) -> Panel {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for snp in snps {
            for byte in format!("{}:{}\n", snp.chromosome, snp.position_bp).bytes() {
                h ^= u64::from(byte);
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        Panel {
            snps: snps.len(),
            fingerprint: format!("{h:016x}"),
        }
    }
}
