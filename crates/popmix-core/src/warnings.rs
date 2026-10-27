//! Things a caller should be told that are not errors: records left out of
//! an input, results that need care. Functions that can warn take a sink,
//! `warn: &mut dyn FnMut(Warning)`; the Python bindings raise each as a
//! Python warning.

use std::fmt;
use std::path::PathBuf;

/// A warning.
///
/// New variants may be added in minor releases.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Warning {
    /// VCF records left out because popmix's formats hold only biallelic
    /// SNPs (two single-base alleles).
    VcfRecordsSkipped {
        path: PathBuf,
        multiallelic: usize,
        not_snp: usize,
    },
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Warning::VcfRecordsSkipped {
                path,
                multiallelic,
                not_snp,
            } => write!(
                f,
                "{}: skipped {} records that are not biallelic SNPs ({multiallelic} multiallelic, \
                 {not_snp} indels or symbolic alleles)",
                path.display(),
                multiallelic + not_snp
            ),
        }
    }
}

/// A sink that drops every warning.
pub fn ignore(_: Warning) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_skipped_vcf_record_warning_says_how_many_and_why() {
        let w = Warning::VcfRecordsSkipped {
            path: "a.vcf".into(),
            multiallelic: 2,
            not_snp: 3,
        };
        assert_eq!(
            w.to_string(),
            "a.vcf: skipped 5 records that are not biallelic SNPs (2 multiallelic, 3 indels or \
             symbolic alleles)"
        );
    }
}
