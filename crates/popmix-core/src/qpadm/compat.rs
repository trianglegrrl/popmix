//! admixtools2 2.0.10 behaviours that popmix corrects by default.
//!
//! As with [`crate::dates::DatesCompat`], popmix reproduces admixtools2
//! exactly so it can be checked against it, fixes its bugs by default, and
//! [`Admixtools2Compat::V2_0_10`] turns every fix back off for the golden
//! tests. Each flag can be set on its own.
//!
//! Every flag is listed, with its evidence and tests, in DEVIATIONS.md.

/// Which admixtools2 2.0.10 quirks to reproduce. Default: none (all corrected).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct Admixtools2Compat {
    /// ADMIXTOOLS2-0001. qpAdm `popdrop` marks every model that drops
    /// exactly one source as `best` and compares those siblings with each
    /// other, instead of walking a nested path. Corrected: the path starts at
    /// the full model and each step drops the one further source that leaves
    /// the best-fitting feasible model; each step is compared with the model
    /// it came from.
    pub popdrop_sibling_path: bool,
    /// ADMIXTOOLS2-0009. qpGraph's best admixture weights are where R's
    /// finite-difference L-BFGS-B stops, which is short of a minimum at a
    /// zero-bound kink. Corrected: the best start's weights are polished by
    /// a line search on each weight (`qpgraph::polish`).
    pub qpgraph_stall_at_kink: bool,
}

impl Admixtools2Compat {
    /// popmix behaviour: every quirk corrected.
    pub const CORRECTED: Admixtools2Compat = Admixtools2Compat {
        popdrop_sibling_path: false,
        qpgraph_stall_at_kink: false,
    };

    /// Reproduce admixtools2 2.0.10 exactly.
    pub const V2_0_10: Admixtools2Compat = Admixtools2Compat {
        popdrop_sibling_path: true,
        qpgraph_stall_at_kink: true,
    };
}
