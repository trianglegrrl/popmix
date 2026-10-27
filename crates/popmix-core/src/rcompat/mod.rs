//! Numerical routines that reproduce R's, for results that match tools
//! written in R: its random number generator and `optim`'s L-BFGS-B.
//!
//! Public, and behaving as R 4.5 does to the last bit (checked against R by
//! `tests/rcompat_golden.rs`); a change in R that these follow is a change
//! of behaviour, listed in the release notes.

mod lbfgsb;
mod rng;

pub use lbfgsb::{optim_lbfgsb, Gradient, LbfgsbOptions, LbfgsbResult};
pub use rng::RRng;
