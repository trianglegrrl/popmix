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

use std::path::Path;

use ndarray::Array2;

use popmix_core::{IoResultExt, PopmixError, Result};

/// Magic bytes that identify a .npy file.
const NPY_MAGIC: &[u8] = b"\x93NUMPY";

/// Write a 2D f64 array to a .npy file.
pub fn write_npy_f64(path: &Path, array: &Array2<f64>) -> Result<()> {
    let (rows, cols) = (array.nrows(), array.ncols());
    let header = format!("{{'descr': '<f8', 'fortran_order': False, 'shape': ({rows}, {cols}), }}");
    // the header, then a newline, padded with spaces so that the data
    // starts at a multiple of 64 bytes (6 magic + 2 version + 2 length)
    let header_len = header.len() + 1;
    let padding = (64 - (10 + header_len) % 64) % 64;
    let mut bytes = Vec::with_capacity(10 + header_len + padding + 8 * rows * cols);
    bytes.extend_from_slice(NPY_MAGIC);
    bytes.extend_from_slice(&[1, 0]);
    bytes.extend_from_slice(&((header_len + padding) as u16).to_le_bytes());
    bytes.extend_from_slice(header.as_bytes());
    bytes.resize(bytes.len() + padding, b' ');
    bytes.push(b'\n');
    // data in C order (row-major)
    for &val in array {
        bytes.extend_from_slice(&val.to_le_bytes());
    }
    std::fs::write(path, bytes).with_path(path)
}

/// Read a 2D f64 array from a .npy file.
pub fn read_npy_f64(path: &Path) -> Result<Array2<f64>> {
    let bytes = std::fs::read(path).with_path(path)?;
    let bad = |details: String| PopmixError::format(path, details);
    let take = |at: usize, n: usize| {
        bytes
            .get(at..at + n)
            .ok_or_else(|| bad(format!("truncated: {} bytes", bytes.len())))
    };

    if bytes.get(..6) != Some(NPY_MAGIC) {
        return Err(bad("not a .npy file (no NPY magic number)".to_string()));
    }
    let head = take(0, 10)?;
    let (major, minor) = (head[6], head[7]);
    // the header length: 2 bytes in version 1, 4 in version 2
    let (header_len, header_at) = match major {
        1 => (u16::from_le_bytes([head[8], head[9]]) as usize, 10),
        2 => {
            let b = take(8, 4)?;
            (u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize, 12)
        }
        _ => return Err(bad(format!("unsupported NPY version {major}.{minor}"))),
    };
    let header = String::from_utf8_lossy(take(header_at, header_len)?);

    let (rows, cols) =
        parse_shape(&header).ok_or_else(|| bad("no 2-D shape in the NPY header".to_string()))?;
    if !header.contains("'<f8'") {
        return Err(bad(format!(
            "expected dtype '<f8' (float64), got: {}",
            header.trim()
        )));
    }

    let data = take(header_at + header_len, 8 * rows * cols)?
        .as_chunks::<8>()
        .0
        .iter()
        .map(|&b| f64::from_le_bytes(b))
        .collect();

    Array2::from_shape_vec((rows, cols), data).map_err(|e| bad(e.to_string()))
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
    fn a_missing_file_is_file_not_found() {
        let err = read_npy_f64(Path::new("/nonexistent/file.npy")).unwrap_err();
        assert!(
            matches!(err, PopmixError::FileNotFound { ref path } if path == Path::new("/nonexistent/file.npy"))
        );
    }

    #[test]
    fn a_damaged_file_is_an_invalid_format_naming_the_file() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("x.npy");
        write_npy_f64(&path, &Array2::from_elem((4, 3), 1.5)).unwrap();
        let whole = std::fs::read(&path).unwrap();
        // not .npy at all, a header cut short, the data cut short
        for bytes in [
            &b"PK\x03\x04 a zip file"[..],
            &whole[..8],
            &whole[..whole.len() - 1],
        ] {
            std::fs::write(&path, bytes).unwrap();
            let err = read_npy_f64(&path).unwrap_err();
            assert!(
                matches!(err, PopmixError::InvalidFileFormat { path: ref p, .. } if *p == path),
                "{err}"
            );
        }
    }

    #[test]
    fn the_layout_is_numpys() {
        // numpy.save of a 2 x 3 float64 array: header padded so data starts at 128
        let dir = tempdir().unwrap();
        let path = dir.path().join("x.npy");
        write_npy_f64(
            &path,
            &Array2::from_shape_fn((2, 3), |(i, j)| (3 * i + j) as f64),
        )
        .unwrap();
        let bytes = std::fs::read(&path).unwrap();
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        assert_eq!(&bytes[..8], b"\x93NUMPY\x01\x00");
        assert_eq!((10 + header_len) % 64, 0);
        assert_eq!(bytes[10 + header_len - 1], b'\n');
        assert_eq!(bytes.len(), 10 + header_len + 6 * 8);
        assert_eq!(
            &bytes[10 + header_len + 8..10 + header_len + 16],
            &1.0f64.to_le_bytes()
        );
    }

    #[test]
    fn test_parse_shape() {
        let header = "{'descr': '<f8', 'fortran_order': False, 'shape': (100, 6), }";
        let (rows, cols) = parse_shape(header).unwrap();
        assert_eq!(rows, 100);
        assert_eq!(cols, 6);
    }
}
