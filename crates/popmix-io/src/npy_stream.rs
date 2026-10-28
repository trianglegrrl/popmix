//! .npy arrays written a few rows at a time, and read a few columns at a
//! time (TASKS.md F4): f2 blocks of many populations, written by chromosome
//! and read back for the populations a command needs, without holding every
//! pair of every block.

use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use ndarray::{Array2, ArrayView2};

use crate::npy::{data_layout, header_bytes};
use popmix_core::{IoResultExt, PopmixError, Result};

/// A C-order float64 .npy file of `cols` columns written row by row; its
/// header is rewritten with the row count when finished. The file is the
/// same bytes as `write_npy_f64` of the whole array.
pub struct NpyRowWriter {
    path: PathBuf,
    out: BufWriter<File>,
    cols: usize,
    rows: usize,
    header_len: usize,
}

impl NpyRowWriter {
    pub fn create(path: &Path, cols: usize) -> Result<Self> {
        let header = header_bytes(0, cols);
        let mut out = BufWriter::new(File::create(path).with_path(path)?);
        out.write_all(&header).with_path(path)?;
        Ok(NpyRowWriter {
            path: path.to_path_buf(),
            out,
            cols,
            rows: 0,
            header_len: header.len(),
        })
    }

    /// Append the rows of `rows`, which must have this file's columns.
    pub fn push(&mut self, rows: ArrayView2<f64>) -> Result<()> {
        if rows.ncols() != self.cols {
            return Err(PopmixError::DimensionMismatch {
                context: format!("rows of {}", self.path.display()),
                expected: format!("{} columns", self.cols),
                actual: rows.ncols().to_string(),
            });
        }
        let mut bytes = Vec::with_capacity(8 * rows.len());
        for &v in rows {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        self.out.write_all(&bytes).with_path(&self.path)?;
        self.rows += rows.nrows();
        Ok(())
    }

    /// Write the final header.
    pub fn finish(self) -> Result<()> {
        let NpyRowWriter {
            path,
            out,
            cols,
            rows,
            header_len,
        } = self;
        let mut file = out
            .into_inner()
            .map_err(std::io::IntoInnerError::into_error)
            .with_path(&path)?;
        let header = header_bytes(rows, cols);
        if header.len() == header_len {
            file.seek(SeekFrom::Start(0)).with_path(&path)?;
            file.write_all(&header).with_path(&path)?;
            return file.flush().with_path(&path);
        }
        // the row count moved the header past a 64-byte boundary: rewrite
        drop(file);
        let bytes = std::fs::read(&path).with_path(&path)?;
        let mut whole = header;
        whole.extend_from_slice(&bytes[header_len..]);
        std::fs::write(&path, whole).with_path(&path)
    }
}

/// The columns `cols` of a C-order float64 .npy file, read through a
/// memory map so that only their values are touched.
pub fn read_npy_columns(path: &Path, cols: &[usize]) -> Result<Array2<f64>> {
    let file = File::open(path).with_path(path)?;
    // SAFETY: the file is opened read-only and not changed while mapped
    let map = unsafe { memmap2::Mmap::map(&file) }.with_path(path)?;
    let (at, rows, ncols) = data_layout(path, &map)?;
    if map.len() < at + 8 * rows * ncols {
        return Err(PopmixError::format(
            path,
            format!("truncated: {} bytes", map.len()),
        ));
    }
    if let Some(&c) = cols.iter().find(|&&c| c >= ncols) {
        return Err(PopmixError::out_of_range("column", c, ncols));
    }
    Ok(Array2::from_shape_fn((rows, cols.len()), |(r, k)| {
        let o = at + 8 * (r * ncols + cols[k]);
        f64::from_le_bytes(map[o..o + 8].try_into().expect("8 bytes"))
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::npy::{read_npy_f64, write_npy_f64};
    use ndarray::s;

    #[test]
    fn rows_written_in_pieces_are_the_whole_arrays_bytes() {
        let dir = tempfile::tempdir().unwrap();
        let a = Array2::from_shape_fn((7, 5), |(i, j)| (i * 10 + j) as f64 + 0.5);
        write_npy_f64(&dir.path().join("whole.npy"), &a).unwrap();
        let mut w = NpyRowWriter::create(&dir.path().join("rows.npy"), 5).unwrap();
        w.push(a.slice(s![..3, ..])).unwrap();
        w.push(a.slice(s![3.., ..])).unwrap();
        w.finish().unwrap();
        let whole = std::fs::read(dir.path().join("whole.npy")).unwrap();
        assert_eq!(std::fs::read(dir.path().join("rows.npy")).unwrap(), whole);
        assert_eq!(read_npy_f64(&dir.path().join("rows.npy")).unwrap(), a);
    }

    #[test]
    fn columns_are_read_alone() {
        let dir = tempfile::tempdir().unwrap();
        let a = Array2::from_shape_fn((4, 6), |(i, j)| (i * 10 + j) as f64);
        write_npy_f64(&dir.path().join("a.npy"), &a).unwrap();
        let got = read_npy_columns(&dir.path().join("a.npy"), &[4, 1]).unwrap();
        assert_eq!(
            got,
            ndarray::stack![ndarray::Axis(1), a.column(4), a.column(1)]
        );
    }
}
