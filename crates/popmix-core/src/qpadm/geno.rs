//! The qpWave/qpAdm f4 matrix from genotypes, as admixtools2's
//! `qpwave(prefix)` and `qpadm(prefix)` (TASKS.md H11): each f4 is computed
//! from allele frequencies at the SNPs, block by block, rather than from f2
//! blocks.
//!
//! - `allsnps = false` (admixtools2's default): the SNPs where every
//!   population of the model has data, monomorphic ones included.
//! - `allsnps = true`: for each f4, the SNPs where its four populations have
//!   data, so different statistics use different SNPs.
//!
//! Each block's statistics are weighted by the most SNPs any of them has
//! there; that weight, and the f4 estimates as block means, were found by
//! running admixtools2 (see [`F4Matrix::from_block_means`]).

use super::f4_matrix::F4Matrix;
use crate::fstats::{f4_from_genotypes_requiring, F4Mode, GenoPopulations, GenotypeRows};
use crate::genotypes::Snp;
use crate::Result;

/// The f4 matrix f4(left[i+1], left[0]; right[j+1], right[0]) from
/// genotypes, `left` and `right` indexing the populations of `pops`, with
/// each block's weight (the most SNPs any statistic has there).
pub fn f4_matrix_from_genotypes(
    rows: &mut dyn GenotypeRows,
    snps: &[Snp],
    pops: &GenoPopulations,
    left: &[usize],
    right: &[usize],
    block_size: f64,
    allsnps: bool,
) -> Result<(F4Matrix, Vec<f64>)> {
    if left.len() < 2 || right.len() < 2 {
        return Err(crate::PopmixError::too_few(
            "qpWave/qpAdm",
            "left and right populations each",
            2,
            left.len().min(right.len()),
        ));
    }
    let quads: Vec<[usize; 4]> = left[1..]
        .iter()
        .flat_map(|&l| right[1..].iter().map(move |&r| [l, left[0], r, right[0]]))
        .collect();
    let model: Vec<usize> = left.iter().chain(right).copied().collect();
    let require: &[usize] = if allsnps { &[] } else { &model };
    let stats =
        f4_from_genotypes_requiring(rows, snps, pops, &quads, F4Mode::F4, block_size, require)?;
    let n_blocks = stats.block_lengths.len();
    let weights: Vec<f64> = (0..n_blocks)
        .map(|k| stats.stats.iter().map(|s| s.blocks.n[k]).max().unwrap_or(0) as f64)
        .collect();
    let means: Vec<Vec<f64>> = stats.stats.iter().map(|s| s.blocks.num.clone()).collect();
    let name = |i: usize| pops.names[i].as_str();
    let left_names: Vec<&str> = left.iter().map(|&i| name(i)).collect();
    let right_names: Vec<&str> = right.iter().map(|&i| name(i)).collect();
    let matrix = F4Matrix::from_block_means(&left_names, &right_names, &means, &weights)?;
    Ok((matrix, weights))
}
