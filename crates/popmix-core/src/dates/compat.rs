//! The DATES v753 behaviours popmix corrects by default, each behind a flag
//! (DEVIATIONS.md, DATES-0001 to DATES-0004).
//!
//! Provenance: written from the clean-room specification
//! `docs/cleanroom/dates.md` (section 8) and `DEVIATIONS.md`, whose evidence
//! comes from DATES v753 run as a black box.

/// Which DATES v753 behaviours to reproduce. [`DatesCompat::CORRECTED`]
/// (the default) reproduces none; [`DatesCompat::V753`] reproduces all four
/// and so matches DATES's own output.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct DatesCompat {
    /// DATES-0001: in mesh mode, lags past the window of real pairs take the
    /// circular cross-correlation's wrapped values (short-range pairs
    /// counted again at long distance).
    pub circular_fft_lags: bool,
    /// DATES-0002: in mesh mode, the pair and cross sums at lag
    /// min(n − 1, diffmax) are left without the 1/n FFT scaling.
    pub unscaled_last_fft_lag: bool,
    /// DATES-0003: fit the correlations rounded to 6 decimals, and run the
    /// jackknife on estimates rounded to 3, as DATES prints them.
    pub fit_printed_values: bool,
    /// DATES-0004: count every sample as diploid in the source frequencies,
    /// so a pseudo-haploid call is two copies of one allele.
    pub every_sample_diploid: bool,
}

impl DatesCompat {
    /// popmix's corrected behaviour (no DATES behaviour reproduced).
    pub const CORRECTED: DatesCompat = DatesCompat {
        circular_fft_lags: false,
        unscaled_last_fft_lag: false,
        fit_printed_values: false,
        every_sample_diploid: false,
    };

    /// Every DATES v753 behaviour reproduced.
    pub const V753: DatesCompat = DatesCompat {
        circular_fft_lags: true,
        unscaled_last_fft_lag: true,
        fit_printed_values: true,
        every_sample_diploid: true,
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn corrected_is_the_default_and_v753_sets_every_flag() {
        assert_eq!(DatesCompat::CORRECTED, DatesCompat::default());
        let v = DatesCompat::V753;
        assert!(
            v.circular_fft_lags
                && v.unscaled_last_fft_lag
                && v.fit_printed_values
                && v.every_sample_diploid
        );
    }
}
