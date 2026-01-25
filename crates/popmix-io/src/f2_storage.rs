//! F2 block storage for precomputed statistics.
//!
//! This module handles reading and writing of precomputed F2 statistics,
//! which can be reused across multiple analyses.
//!
//! ## File Format
//!
//! F2 blocks are stored in a directory with:
//! - `metadata.json`: Human-readable metadata (populations, block info)
//! - `block_sums.npy`: NumPy format 2D array of F2 sums [n_blocks, n_pairs]
//! - `block_weights.npy`: NumPy format 2D array of weights [n_blocks, n_pairs]
//! - `block_het_sums.npy`: NumPy format 2D array of heterozygosity sums [n_blocks, n_pops] (v2+)
//! - `block_het_weights.npy`: NumPy format 2D array of het weights [n_blocks, n_pops] (v2+)

use std::fs::{self, File};
use std::io::{BufReader, BufWriter};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::npy::{read_npy_f64, write_npy_f64};
use popmix_core::fstats::F2Blocks;
use popmix_core::genotypes::PopulationId;
use popmix_core::{PopmixError, Result};

/// Metadata structure for F2 block storage.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct F2Metadata {
    /// Format version for future compatibility.
    /// v1: F2 only
    /// v2: F2 + heterozygosity (for FST)
    version: u32,
    /// Population names in order.
    populations: Vec<String>,
    /// Number of blocks.
    n_blocks: usize,
    /// Number of population pairs.
    n_pairs: usize,
    /// Number of populations (v2+, inferred from populations.len() for v1).
    #[serde(default)]
    n_populations: usize,
    /// Number of SNPs per block.
    block_lengths: Vec<usize>,
    /// Optional source file path.
    #[serde(skip_serializing_if = "Option::is_none")]
    source: Option<String>,
}

/// Read F2 blocks from a directory.
///
/// The directory must contain:
/// - `metadata.json`: Population names, block lengths, version info
/// - `block_sums.npy`: NumPy format array of F2 sums
/// - `block_weights.npy`: NumPy format array of weights
/// - `block_het_sums.npy`: NumPy format array of het sums (v2+, optional for v1)
/// - `block_het_weights.npy`: NumPy format array of het weights (v2+, optional for v1)
pub fn read_f2_blocks(dir: &Path) -> Result<F2Blocks> {
    // Read metadata
    let metadata_path = dir.join("metadata.json");
    let file = File::open(&metadata_path).map_err(|e| {
        PopmixError::Io(std::io::Error::new(
            e.kind(),
            format!("Failed to open {}: {}", metadata_path.display(), e),
        ))
    })?;
    let reader = BufReader::new(file);
    let metadata: F2Metadata = serde_json::from_reader(reader).map_err(|e| {
        PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Failed to parse {}: {}", metadata_path.display(), e),
        ))
    })?;

    // Validate version
    if metadata.version < 1 || metadata.version > 3 {
        return Err(PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Unsupported F2 storage version {} in {}",
                metadata.version,
                dir.display()
            ),
        )));
    }

    // Read F2 arrays
    let block_sums = read_npy_f64(&dir.join("block_sums.npy"))?;
    let block_weights = read_npy_f64(&dir.join("block_weights.npy"))?;

    // Validate F2 array shapes
    if block_sums.nrows() != metadata.n_blocks || block_sums.ncols() != metadata.n_pairs {
        return Err(PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "block_sums shape {:?} doesn't match metadata ({}, {})",
                block_sums.shape(),
                metadata.n_blocks,
                metadata.n_pairs
            ),
        )));
    }

    if block_weights.shape() != block_sums.shape() {
        return Err(PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "block_weights shape {:?} doesn't match block_sums {:?}",
                block_weights.shape(),
                block_sums.shape()
            ),
        )));
    }

    // Read heterozygosity arrays (v2+) or create zeros for v1
    let n_pops = metadata.populations.len();
    let (block_het_sums, block_het_weights) = if metadata.version >= 2 {
        let het_sums = read_npy_f64(&dir.join("block_het_sums.npy"))?;
        let het_weights = read_npy_f64(&dir.join("block_het_weights.npy"))?;

        // Validate heterozygosity array shapes
        if het_sums.nrows() != metadata.n_blocks || het_sums.ncols() != n_pops {
            return Err(PopmixError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "block_het_sums shape {:?} doesn't match metadata ({}, {})",
                    het_sums.shape(),
                    metadata.n_blocks,
                    n_pops
                ),
            )));
        }

        if het_weights.shape() != het_sums.shape() {
            return Err(PopmixError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "block_het_weights shape {:?} doesn't match block_het_sums {:?}",
                    het_weights.shape(),
                    het_sums.shape()
                ),
            )));
        }

        (het_sums, het_weights)
    } else {
        // v1 format: create zero arrays for heterozygosity (FST will return None)
        (
            ndarray::Array2::zeros((metadata.n_blocks, n_pops)),
            ndarray::Array2::zeros((metadata.n_blocks, n_pops)),
        )
    };

    // Read raw F2 arrays (v3+) or create zeros for older versions
    let block_raw_sums = if metadata.version >= 3 {
        read_npy_f64(&dir.join("block_raw_sums.npy"))?
    } else {
        // For older versions, use block_sums as approximation (won't give exact FST)
        block_sums.clone()
    };

    // Convert population names
    let populations: Vec<PopulationId> = metadata
        .populations
        .iter()
        .map(|s| PopulationId::from(s.as_str()))
        .collect();

    Ok(F2Blocks {
        block_sums,
        block_raw_sums,
        block_weights,
        block_het_sums,
        block_het_weights,
        populations,
        n_blocks: metadata.n_blocks,
        block_lengths: metadata.block_lengths,
    })
}

/// Write F2 blocks to a directory.
///
/// Creates:
/// - `metadata.json`: Population names, block lengths, version info
/// - `block_sums.npy`: NumPy format array of F2 sums
/// - `block_weights.npy`: NumPy format array of weights
/// - `block_het_sums.npy`: NumPy format array of heterozygosity sums
/// - `block_het_weights.npy`: NumPy format array of heterozygosity weights
///
/// # Arguments
/// * `f2_blocks` - The F2 blocks to save
/// * `dir` - Directory path to write to (will be created if needed)
/// * `source` - Optional source file path for metadata
pub fn write_f2_blocks(f2_blocks: &F2Blocks, dir: &Path, source: Option<&str>) -> Result<()> {
    // Create directory if needed
    fs::create_dir_all(dir).map_err(|e| {
        PopmixError::Io(std::io::Error::new(
            e.kind(),
            format!("Failed to create directory {}: {}", dir.display(), e),
        ))
    })?;

    // Calculate n_pairs and n_populations
    let n_pairs = f2_blocks.block_sums.ncols();
    let n_populations = f2_blocks.populations.len();

    // Build metadata (version 3 includes raw sums for FST)
    let metadata = F2Metadata {
        version: 3,
        populations: f2_blocks
            .populations
            .iter()
            .map(|p| p.as_str().to_string())
            .collect(),
        n_blocks: f2_blocks.n_blocks,
        n_pairs,
        n_populations,
        block_lengths: f2_blocks.block_lengths.clone(),
        source: source.map(|s| s.to_string()),
    };

    // Write metadata
    let metadata_path = dir.join("metadata.json");
    let file = File::create(&metadata_path).map_err(|e| {
        PopmixError::Io(std::io::Error::new(
            e.kind(),
            format!("Failed to create {}: {}", metadata_path.display(), e),
        ))
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &metadata).map_err(|e| {
        PopmixError::Io(std::io::Error::other(format!(
            "Failed to write {}: {}",
            metadata_path.display(),
            e
        )))
    })?;

    // Write F2 arrays
    write_npy_f64(&dir.join("block_sums.npy"), &f2_blocks.block_sums)?;
    write_npy_f64(&dir.join("block_raw_sums.npy"), &f2_blocks.block_raw_sums)?;
    write_npy_f64(&dir.join("block_weights.npy"), &f2_blocks.block_weights)?;

    // Write heterozygosity arrays (v2+)
    write_npy_f64(&dir.join("block_het_sums.npy"), &f2_blocks.block_het_sums)?;
    write_npy_f64(
        &dir.join("block_het_weights.npy"),
        &f2_blocks.block_het_weights,
    )?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;
    use tempfile::tempdir;

    /// Create test F2Blocks with known values.
    fn create_test_f2_blocks() -> F2Blocks {
        // 3 populations = 6 pairs (including self-pairs)
        // 4 blocks
        let n_blocks = 4;
        let n_pairs = 6;
        let n_pops = 3;

        let mut block_sums = Array2::zeros((n_blocks, n_pairs));
        let mut block_raw_sums = Array2::zeros((n_blocks, n_pairs));
        let mut block_weights = Array2::zeros((n_blocks, n_pairs));
        let mut block_het_sums = Array2::zeros((n_blocks, n_pops));
        let mut block_het_weights = Array2::zeros((n_blocks, n_pops));

        // Fill with recognizable test values
        for i in 0..n_blocks {
            for j in 0..n_pairs {
                block_sums[[i, j]] = (i * 10 + j) as f64 * 0.001;
                // Raw sums are slightly larger (before bias correction)
                block_raw_sums[[i, j]] = (i * 10 + j) as f64 * 0.0012;
                block_weights[[i, j]] = (100 + i * n_pairs + j) as f64;
            }
            for j in 0..n_pops {
                // Heterozygosity values between 0 and 0.5
                block_het_sums[[i, j]] = (i * 3 + j) as f64 * 0.02;
                block_het_weights[[i, j]] = (100 + i * n_pops + j) as f64;
            }
        }

        let populations = vec![
            PopulationId::from("PopA"),
            PopulationId::from("PopB"),
            PopulationId::from("PopC"),
        ];

        let block_lengths = vec![100, 120, 110, 95];

        F2Blocks {
            block_sums,
            block_raw_sums,
            block_weights,
            block_het_sums,
            block_het_weights,
            populations,
            n_blocks,
            block_lengths,
        }
    }

    #[test]
    fn test_write_read_roundtrip() {
        let f2_blocks = create_test_f2_blocks();
        let dir = tempdir().unwrap();

        // Write
        write_f2_blocks(&f2_blocks, dir.path(), Some("/test/source")).unwrap();

        // Read back
        let loaded = read_f2_blocks(dir.path()).unwrap();

        // Verify all fields match
        assert_eq!(loaded.populations.len(), f2_blocks.populations.len());
        for (a, b) in loaded.populations.iter().zip(f2_blocks.populations.iter()) {
            assert_eq!(a.as_str(), b.as_str());
        }
        assert_eq!(loaded.n_blocks, f2_blocks.n_blocks);
        assert_eq!(loaded.block_lengths, f2_blocks.block_lengths);

        // F2 arrays should be exactly equal (no floating point drift in binary format)
        assert_eq!(loaded.block_sums.shape(), f2_blocks.block_sums.shape());
        assert_eq!(loaded.block_weights.shape(), f2_blocks.block_weights.shape());

        for i in 0..f2_blocks.n_blocks {
            for j in 0..6 {
                assert_eq!(
                    loaded.block_sums[[i, j]],
                    f2_blocks.block_sums[[i, j]],
                    "block_sums mismatch at [{}, {}]",
                    i,
                    j
                );
                assert_eq!(
                    loaded.block_raw_sums[[i, j]],
                    f2_blocks.block_raw_sums[[i, j]],
                    "block_raw_sums mismatch at [{}, {}]",
                    i,
                    j
                );
                assert_eq!(
                    loaded.block_weights[[i, j]],
                    f2_blocks.block_weights[[i, j]],
                    "block_weights mismatch at [{}, {}]",
                    i,
                    j
                );
            }
        }

        // Heterozygosity arrays should also match
        assert_eq!(
            loaded.block_het_sums.shape(),
            f2_blocks.block_het_sums.shape()
        );
        assert_eq!(
            loaded.block_het_weights.shape(),
            f2_blocks.block_het_weights.shape()
        );

        for i in 0..f2_blocks.n_blocks {
            for j in 0..3 {
                assert_eq!(
                    loaded.block_het_sums[[i, j]],
                    f2_blocks.block_het_sums[[i, j]],
                    "block_het_sums mismatch at [{}, {}]",
                    i,
                    j
                );
                assert_eq!(
                    loaded.block_het_weights[[i, j]],
                    f2_blocks.block_het_weights[[i, j]],
                    "block_het_weights mismatch at [{}, {}]",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_metadata_json_valid() {
        let f2_blocks = create_test_f2_blocks();
        let dir = tempdir().unwrap();

        write_f2_blocks(&f2_blocks, dir.path(), None).unwrap();

        // Read and parse metadata
        let metadata_path = dir.path().join("metadata.json");
        let content = std::fs::read_to_string(&metadata_path).unwrap();
        let meta: serde_json::Value = serde_json::from_str(&content).unwrap();

        // Check required fields
        assert_eq!(meta["version"], 3); // v3 includes raw sums for FST
        assert!(meta["populations"].is_array());
        assert_eq!(meta["populations"].as_array().unwrap().len(), 3);
        assert_eq!(meta["n_blocks"], 4);
        assert_eq!(meta["n_pairs"], 6);
        assert_eq!(meta["n_populations"], 3);
        assert!(meta["block_lengths"].is_array());
    }

    #[test]
    fn test_read_missing_directory_errors() {
        let result = read_f2_blocks(Path::new("/nonexistent/path/that/does/not/exist"));
        assert!(result.is_err());
    }

    #[test]
    fn test_read_corrupted_metadata_errors() {
        let dir = tempdir().unwrap();
        std::fs::write(dir.path().join("metadata.json"), "invalid json{{{").unwrap();

        let result = read_f2_blocks(dir.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_files_created() {
        let f2_blocks = create_test_f2_blocks();
        let dir = tempdir().unwrap();

        write_f2_blocks(&f2_blocks, dir.path(), None).unwrap();

        // All six files should exist (v3 format)
        assert!(dir.path().join("metadata.json").exists());
        assert!(dir.path().join("block_sums.npy").exists());
        assert!(dir.path().join("block_raw_sums.npy").exists());
        assert!(dir.path().join("block_weights.npy").exists());
        assert!(dir.path().join("block_het_sums.npy").exists());
        assert!(dir.path().join("block_het_weights.npy").exists());
    }
}
