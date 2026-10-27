//! DATES: dating admixture from the decay of weighted ancestry covariance.
//!
//! An admixed (test) population inherits, from an admixture `n` generations
//! ago, ancestry tracts whose covariance between two SNPs decays as
//! exp(−n·d) with their genetic distance d. Weighting each SNP by the
//! difference of the two source frequencies, correlating the residual
//! genotypes of each test individual over pairs of SNPs, binning by distance
//! and fitting A·exp(−n·d) + K gives the number of generations; a
//! leave-one-chromosome-out weighted jackknife gives its standard error.
//!
//! Stages: `weights` (SNP weights), `ld_decay` (the binned weighted LD
//! table, by every pair or on an FFT mesh), `exponential_fit` (the fit and
//! the weighted jackknife), `dating` (the whole analysis), `years` (the date
//! in years) and `compat` (the DATES v753 behaviours popmix corrects by
//! default).
//!
//! Provenance: written from the clean-room specification
//! `docs/cleanroom/dates.md` (and `docs/cleanroom/dates-mesh.md`),
//! Narasimhan et al. 2019 (Science 365, eaat7487, supplementary methods),
//! Chintalapati, Patterson & Moorjani 2022 (eLife 11:e77625), Busing, Meijer
//! & van der Leeden 1999 (the weighted delete-m jackknife), DATES's README and
//! parameter-file documentation, and the output of DATES v753 run as a black
//! box (`tests/golden/dates/v753/`).

mod compat;
mod dating;
mod exponential_fit;
mod ld_decay;
mod weights;
mod years;

pub use compat::DatesCompat;
pub use dating::{dates, DatesOptions, DatesResult, JackknifeBlock};
pub use exponential_fit::{fit_exponential, weighted_jackknife, DecayFit, DecayFitOptions};
pub use ld_decay::{
    compute_ld_decay, compute_ld_decay_by_chromosome, ChromosomeLd, LdBin, LdDecay, LdOptions,
};
pub use weights::{compute_weights, WeightResult};
