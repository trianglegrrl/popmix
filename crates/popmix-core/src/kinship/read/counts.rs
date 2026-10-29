//! Per-pair mismatch counts from the bit sets of `BitGenotypes`.
//!
//! For a pair and a 64-SNP word: called in both = (het|hom0|hom2) of one
//! sample & the same of the other; full mismatches = opposite homozygotes;
//! half mismatches = called in both & a heterozygote in either. Twice the
//! mismatch sum, 2D, and the SNP count N are then integers, so every P0 is
//! one exact division and does not depend on summation order.

use crate::kinship::king::BitGenotypes;

/// A window as the words it covers, each with the mask of its SNPs.
pub(super) type Segments = Vec<(usize, u64)>;

/// Twice the mismatch sum and the number of SNPs called in both.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct Count {
    pub d2: u64,
    pub n: u64,
}

impl Count {
    /// P0 = D / N; `None` when no SNP is called in both.
    pub fn p0(self) -> Option<f64> {
        (self.n > 0).then(|| self.d2 as f64 / (2 * self.n) as f64)
    }
}

/// The word segments of each inclusive SNP-index window.
pub(super) fn segments(windows: &[(usize, usize)]) -> Vec<Segments> {
    windows
        .iter()
        .map(|&(a, b)| {
            let (wa, wb) = (a / 64, b / 64);
            (wa..=wb)
                .map(|w| {
                    let lo = if w == wa { a % 64 } else { 0 };
                    let hi = if w == wb { b % 64 } else { 63 };
                    let upper = if hi == 63 {
                        u64::MAX
                    } else {
                        (1u64 << (hi + 1)) - 1
                    };
                    (w, upper & (u64::MAX << lo))
                })
                .collect()
        })
        .collect()
}

/// The bit sets of one pair.
pub(super) struct PairBits<'a> {
    a: (&'a [u64], &'a [u64], &'a [u64]),
    b: (&'a [u64], &'a [u64], &'a [u64]),
}

impl<'a> PairBits<'a> {
    pub fn new(g: &'a BitGenotypes, i: usize, j: usize) -> Self {
        PairBits {
            a: g.masks(i),
            b: g.masks(j),
        }
    }

    /// The counts of word `w` restricted to `mask`.
    #[inline]
    fn word(&self, w: usize, mask: u64) -> Count {
        let (het_a, h0_a, h2_a) = (self.a.0[w], self.a.1[w], self.a.2[w]);
        let (het_b, h0_b, h2_b) = (self.b.0[w], self.b.1[w], self.b.2[w]);
        let called = (het_a | h0_a | h2_a) & (het_b | h0_b | h2_b) & mask;
        let full = ((h0_a & h2_b) | (h2_a & h0_b)) & mask;
        let half = called & (het_a | het_b);
        Count {
            d2: 2 * u64::from(full.count_ones()) + u64::from(half.count_ones()),
            n: u64::from(called.count_ones()),
        }
    }

    /// The counts over all SNPs, each once.
    pub fn total(&self) -> Count {
        (0..self.a.0.len()).fold(Count::default(), |acc, w| add(acc, self.word(w, u64::MAX)))
    }

    /// The counts of each window.
    pub fn windows(&self, segs: &[Segments]) -> Vec<Count> {
        segs.iter()
            .map(|s| {
                s.iter()
                    .fold(Count::default(), |acc, &(w, m)| add(acc, self.word(w, m)))
            })
            .collect()
    }
}

fn add(a: Count, b: Count) -> Count {
    Count {
        d2: a.d2 + b.d2,
        n: a.n + b.n,
    }
}
