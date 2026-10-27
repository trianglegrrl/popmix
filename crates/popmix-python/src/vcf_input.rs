//! VCF input for the bindings: finding a prefix's VCF file and opening it,
//! with a Python warning for the records popmix's formats cannot hold.

use popmix_io::VcfReader;
use std::path::Path;

use pyo3::prelude::*;

/// Open a VCF, warning (UserWarning) about the records skipped because they
/// are not biallelic SNPs: popmix's formats hold two single-base alleles.
pub(crate) fn open_vcf(path: &str) -> PyResult<VcfReader> {
    let reader = VcfReader::new(path)?;
    let k = reader.skipped();
    if k.total() > 0 {
        Python::attach(|py| {
            PyErr::warn(
                py,
                &py.get_type::<pyo3::exceptions::PyUserWarning>(),
                &std::ffi::CString::new(format!(
                    "{path}: skipped {} records that are not biallelic SNPs ({} multiallelic, {} indels or symbolic alleles)",
                    k.total(),
                    k.multiallelic,
                    k.not_snp
                ))?,
                1,
            )
        })?;
    }
    Ok(reader)
}

/// The VCF file of `prefix`: the path itself, or with .vcf.gz or .vcf added.
pub(crate) fn vcf_path(prefix: &str) -> Option<String> {
    if (prefix.ends_with(".vcf") || prefix.ends_with(".vcf.gz")) && Path::new(prefix).exists() {
        return Some(prefix.to_string());
    }
    [".vcf.gz", ".vcf"]
        .iter()
        .map(|e| format!("{prefix}{e}"))
        .find(|p| Path::new(p).exists())
}
