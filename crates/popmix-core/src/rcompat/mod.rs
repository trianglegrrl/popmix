//! Numerical routines that reproduce R's, for results that match tools
//! written in R: its random number generator and `optim`'s L-BFGS-B.

mod lbfgsb;
mod rng;

pub use lbfgsb::{optim_lbfgsb, Gradient, LbfgsbOptions, LbfgsbResult};
pub use rng::RRng;
