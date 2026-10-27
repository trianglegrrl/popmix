//! Checks on numeric options shared by the LD-decay methods.

use crate::error::{PopmixError, Result};

/// Most distance bins a run may ask for. Each bin holds a few sums per
/// chromosome; past this a mistyped bin size would exhaust memory.
pub(crate) const MAX_BINS: f64 = 1e6;

/// `value` must be a finite number greater than 0.
pub(crate) fn check_positive(name: &str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(PopmixError::parameter(
            name,
            format!("must be a finite number greater than 0, not {value}"),
        ))
    }
}

/// `value` must be a finite number at least 0.
pub(crate) fn check_non_negative(name: &str, value: f64) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(PopmixError::parameter(
            name,
            format!("must be a finite number at least 0, not {value}"),
        ))
    }
}

/// Bins of width `binsize` out to `maxdis`: both positive, and at most
/// [`MAX_BINS`] of them.
pub(crate) fn check_bins(bin_name: &str, binsize: f64, max_name: &str, maxdis: f64) -> Result<()> {
    check_positive(bin_name, binsize)?;
    check_positive(max_name, maxdis)?;
    if maxdis / binsize > MAX_BINS {
        return Err(PopmixError::parameter(
            bin_name,
            format!(
                "{binsize} makes {:.0} bins out to {max_name} {maxdis}; at most {MAX_BINS:.0} are allowed",
                maxdis / binsize
            ),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_must_be_positive_and_not_too_many() {
        assert!(check_bins("binsize", 0.001, "maxdis", 0.5).is_ok());
        for (binsize, maxdis) in [
            (0.0, 0.5),
            (-0.1, 0.5),
            (f64::NAN, 0.5),
            (f64::INFINITY, 0.5),
            (1e-12, 0.5),
        ] {
            let err = check_bins("binsize", binsize, "maxdis", maxdis)
                .unwrap_err()
                .to_string();
            assert!(err.contains("binsize"), "{binsize}: {err}");
        }
        for maxdis in [0.0, f64::NAN, -1.0] {
            let err = check_bins("binsize", 0.001, "maxdis", maxdis)
                .unwrap_err()
                .to_string();
            assert!(err.contains("maxdis"), "{maxdis}: {err}");
        }
    }
}
