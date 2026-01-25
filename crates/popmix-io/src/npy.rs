//! NumPy .npy file format reader/writer for 2D f64 arrays.
//!
//! This module implements the NumPy v1.0 .npy format for storing 2D arrays
//! of f64 values. The format is simple and well-documented:
//!
//! ## Format Specification (v1.0)
//!
//! - Magic number: `\x93NUMPY` (6 bytes)
//! - Version: `\x01\x00` (major=1, minor=0)
//! - Header length: 2 bytes, little-endian uint16
//! - Header: Python dict literal as ASCII, padded with spaces to 64-byte alignment
//! - Data: Raw binary data in C order (row-major)
//!
//! The header dict contains:
//! - `'descr'`: Data type descriptor (e.g., `'<f8'` for little-endian float64)
//! - `'fortran_order'`: Whether data is in Fortran order (always `False` for us)
//! - `'shape'`: Tuple of dimensions (e.g., `(100, 6)`)

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use ndarray::Array2;

use popmix_core::{PopmixError, Result};

/// Magic bytes that identify a .npy file.
const NPY_MAGIC: &[u8] = b"\x93NUMPY";

/// Write a 2D f64 array to a .npy file.
pub fn write_npy_f64(path: &Path, array: &Array2<f64>) -> Result<()> {
    let file = File::create(path).map_err(|e| {
        PopmixError::Io(std::io::Error::new(
            e.kind(),
            format!("Failed to create {}: {}", path.display(), e),
        ))
    })?;
    let mut writer = BufWriter::new(file);

    // Write magic number
    writer.write_all(NPY_MAGIC)?;

    // Write version (1.0)
    writer.write_u8(1)?;
    writer.write_u8(0)?;

    // Build header dict
    let (rows, cols) = (array.nrows(), array.ncols());
    let header = format!(
        "{{'descr': '<f8', 'fortran_order': False, 'shape': ({}, {}), }}",
        rows, cols
    );

    // Calculate padding to align data to 64 bytes
    // Header length includes the header string plus newline
    let header_len = header.len() + 1; // +1 for newline
    let total_header_size = 10 + header_len; // 6 magic + 2 version + 2 len + header + newline
    let padding_needed = (64 - (total_header_size % 64)) % 64;
    let padded_header_len = header_len + padding_needed;

    // Write header length (little-endian u16)
    writer.write_u16::<LittleEndian>(padded_header_len as u16)?;

    // Write header with padding
    writer.write_all(header.as_bytes())?;
    for _ in 0..padding_needed {
        writer.write_u8(b' ')?;
    }
    writer.write_u8(b'\n')?;

    // Write data in C order (row-major)
    for row in array.rows() {
        for &val in row.iter() {
            writer.write_f64::<LittleEndian>(val)?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Read a 2D f64 array from a .npy file.
pub fn read_npy_f64(path: &Path) -> Result<Array2<f64>> {
    let file = File::open(path).map_err(|e| {
        PopmixError::Io(std::io::Error::new(
            e.kind(),
            format!("Failed to open {}: {}", path.display(), e),
        ))
    })?;
    let mut reader = BufReader::new(file);

    // Read and verify magic number
    let mut magic = [0u8; 6];
    reader.read_exact(&mut magic)?;
    if magic != NPY_MAGIC {
        return Err(PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Invalid NPY magic number in {}", path.display()),
        )));
    }

    // Read version
    let major = reader.read_u8()?;
    let minor = reader.read_u8()?;

    // Read header length based on version
    let header_len = if major == 1 {
        reader.read_u16::<LittleEndian>()? as usize
    } else if major == 2 {
        reader.read_u32::<LittleEndian>()? as usize
    } else {
        return Err(PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Unsupported NPY version {}.{} in {}",
                major,
                minor,
                path.display()
            ),
        )));
    };

    // Read header
    let mut header_bytes = vec![0u8; header_len];
    reader.read_exact(&mut header_bytes)?;
    let header = String::from_utf8_lossy(&header_bytes);

    // Parse shape from header
    let (rows, cols) = parse_shape(&header).ok_or_else(|| {
        PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Failed to parse shape from NPY header in {}", path.display()),
        ))
    })?;

    // Verify dtype is f64
    if !header.contains("'<f8'") && !header.contains("'<f8'") {
        return Err(PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Expected dtype '<f8' (float64) in {}, got: {}",
                path.display(),
                header.trim()
            ),
        )));
    }

    // Read data
    let mut data = Vec::with_capacity(rows * cols);
    for _ in 0..(rows * cols) {
        data.push(reader.read_f64::<LittleEndian>()?);
    }

    // Convert to Array2
    Array2::from_shape_vec((rows, cols), data).map_err(|e| {
        PopmixError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Failed to create array from NPY data: {}", e),
        ))
    })
}

/// Parse shape tuple from NPY header string.
fn parse_shape(header: &str) -> Option<(usize, usize)> {
    // Find 'shape': (rows, cols)
    let shape_start = header.find("'shape':")? + 8;
    let after_shape = &header[shape_start..];

    // Find the opening paren
    let paren_start = after_shape.find('(')? + 1;
    let paren_end = after_shape.find(')')?;
    let shape_str = &after_shape[paren_start..paren_end];

    // Parse the two numbers
    let parts: Vec<&str> = shape_str.split(',').collect();
    if parts.len() < 2 {
        return None;
    }

    let rows: usize = parts[0].trim().parse().ok()?;
    let cols: usize = parts[1].trim().parse().ok()?;

    Some((rows, cols))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_write_read_roundtrip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.npy");

        // Create test array
        let array = Array2::from_shape_fn((5, 3), |(i, j)| (i * 10 + j) as f64);

        // Write and read back
        write_npy_f64(&path, &array).unwrap();
        let loaded = read_npy_f64(&path).unwrap();

        assert_eq!(loaded.shape(), array.shape());
        for i in 0..5 {
            for j in 0..3 {
                assert_eq!(loaded[[i, j]], array[[i, j]]);
            }
        }
    }

    #[test]
    fn test_preserves_exact_values() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.npy");

        // Use values that could have precision issues
        let array = Array2::from_shape_fn((3, 2), |(i, j)| {
            std::f64::consts::PI * (i as f64 + 0.1 * j as f64)
        });

        write_npy_f64(&path, &array).unwrap();
        let loaded = read_npy_f64(&path).unwrap();

        // Should be bit-for-bit identical
        for i in 0..3 {
            for j in 0..2 {
                assert_eq!(loaded[[i, j]].to_bits(), array[[i, j]].to_bits());
            }
        }
    }

    #[test]
    fn test_read_nonexistent_errors() {
        let result = read_npy_f64(Path::new("/nonexistent/file.npy"));
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_shape() {
        let header = "{'descr': '<f8', 'fortran_order': False, 'shape': (100, 6), }";
        let (rows, cols) = parse_shape(header).unwrap();
        assert_eq!(rows, 100);
        assert_eq!(cols, 6);
    }
}
