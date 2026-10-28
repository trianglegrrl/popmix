//! f2 blocks extracted straight into an f2 directory (TASKS.md F4), as
//! admixtools2's `extract_f2(pref, outdir)`: from a dataset read in place,
//! one chromosome's frequency table and blocks are held at a time, so the
//! memory no longer grows with every block of every pair. The directory is
//! the same bytes as `write_f2_blocks` of `extract_f2`'s blocks.

use std::path::Path;

use popmix_core::{BlockDefinition, F2Blocks, PopulationId, Result, Warning};

use crate::dataset::open_mapped;
use crate::extract::{check_options, in_memory, ExtractOptions, MappedRuns};
use crate::f2_stream::F2DirWriter;
use crate::write_f2_blocks;

/// Extract the f2 blocks of `opts.populations` from `prefix` into `dir`,
/// recording `source`; returns the number of blocks.
pub fn extract_f2_to_dir(
    prefix: &str,
    opts: &ExtractOptions,
    dir: &Path,
    source: Option<&str>,
    warn: &mut dyn FnMut(Warning),
) -> Result<usize> {
    check_options(opts)?;
    let pops: Vec<&str> = opts.populations.iter().map(String::as_str).collect();
    let Some(mapped) = open_mapped(prefix, warn)? else {
        let blocks = in_memory(prefix, opts, &pops, warn)?;
        write_f2_blocks(&blocks, dir, source)?;
        return Ok(blocks.n_blocks());
    };
    let m = MappedRuns::prepare(mapped, opts, &pops, warn)?;
    let auto_only = opts.autosomes_only();
    let populations: Vec<PopulationId> = pops.iter().map(|p| PopulationId::from(*p)).collect();
    let def: &BlockDefinition = &opts.blocks;
    let guess = def.coordinate(&m.evidence);
    let mut coordinates = (guess, guess);
    // the coordinates the kept SNPs choose are known at the end: when they
    // are not the guess, a second pass in them gives the whole's blocks
    for _ in 0..2 {
        let mut writer = F2DirWriter::create(dir, &populations)?;
        let parts = m.runs.iter().map(|r| {
            m.file
                .allele_frequencies_at(&m.samples, &m.snps[r.clone()], r.start, &pops, auto_only)
        });
        let evidence = F2Blocks::stream_parts(parts, coordinates, def, &opts.f2, &mut |part| {
            writer.push(&part)
        })?;
        let chosen = (def.coordinate(&evidence.0), def.coordinate(&evidence.1));
        if chosen == coordinates {
            let provenance = crate::provenance::extraction(
                prefix,
                opts,
                &m.samples,
                &crate::provenance::Panel::of(&m.snps),
                writer.n_snps(),
            );
            return writer.finish(source, Some(provenance));
        }
        coordinates = chosen;
    }
    Err(popmix_core::PopmixError::Internal(
        "the block coordinates did not settle in two passes".to_string(),
    ))
}
