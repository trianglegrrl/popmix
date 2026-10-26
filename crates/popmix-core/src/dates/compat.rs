//! DATES v753 behaviours that popmix corrects by default.
//!
//! popmix reproduces DATES v753 exactly so it can be checked against it, but
//! a few DATES behaviours are bugs or artefacts that change results. popmix
//! fixes these by default; [`DatesCompat::V753`] turns every fix back off, so
//! the golden tests can still require exact agreement with DATES. Each flag
//! can also be set on its own, so a test can show that a fix changes only
//! what it is meant to.
//!
//! Every flag is listed, with its evidence and tests, in DEVIATIONS.md.

/// Which DATES v753 quirks to reproduce. Default: none (all corrected).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DatesCompat {
    /// DATES-0001. The mesh FFT is circular, and DATES sizes it so that lags
    /// beyond n - m (n = FFT length, m = the chromosome's mesh window) also
    /// pick up wrapped short-range pairs. Corrected: lags with no real pairs
    /// are zero.
    pub circular_fft_lags: bool,
    /// DATES-0002. `fftauto` scales its output lags by 1/n except the last,
    /// outmax = min(n - 1, diffmax). Corrected: every lag is scaled.
    pub unscaled_last_fft_lag: bool,
    /// DATES-0003. DATES fits the correlations as printed in its `.out` file
    /// (%12.6f) and runs the jackknife on estimates as printed in its `.jin`
    /// file (%9.3f). Corrected: full precision throughout.
    pub fit_printed_values: bool,
    /// DATES-0004. DATES counts every sample as diploid: a pseudohaploid
    /// call (0 or 2) is two copies of one allele, so in a source population
    /// that mixes pseudohaploid and diploid samples each pseudohaploid
    /// sample weighs as much as a diploid one in its allele frequency.
    /// Corrected: a pseudohaploid call is one allele, as everywhere else in
    /// popmix. Populations of one ploidy give the same frequencies either way.
    pub every_sample_diploid: bool,
}

impl DatesCompat {
    /// popmix behaviour: every quirk corrected.
    pub const CORRECTED: DatesCompat = DatesCompat {
        circular_fft_lags: false,
        unscaled_last_fft_lag: false,
        fit_printed_values: false,
        every_sample_diploid: false,
    };

    /// Reproduce DATES v753 exactly.
    pub const V753: DatesCompat = DatesCompat {
        circular_fft_lags: true,
        unscaled_last_fft_lag: true,
        fit_printed_values: true,
        every_sample_diploid: true,
    };
}
