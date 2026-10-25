//! R's default uniform random numbers (`RNGkind()`: "Mersenne-Twister",
//! "Inversion", "Rejection"), so random starts and resamples match tools
//! written in R draw for draw.
//!
//! Sources: Matsumoto & Nishimura 1998 (ACM TOMACS 8:3) and their reference
//! implementation `mt19937ar.c` (2002, BSD-3; see NOTICE) for the generator
//! and its tempering; R's documentation (`?Random`, `?set.seed`, `?runif`,
//! `?sample`, and the R 3.6.0 NEWS entry on `sample.kind = "Rejection"`) and
//! R 4.5.2 run as a black box (`.Random.seed` after `set.seed`, draws of
//! `runif` and `sample`) for the seeding, the scaling and the index
//! sampling; `docs/cleanroom/small-helpers.md` section 3.
//!
//! - Seeding: the seed, as an unsigned 32-bit integer, is scrambled by 50
//!   steps of the congruential map s -> 69069 s + 1 (mod 2^32); 625 further
//!   steps give R's state words, the first of which is a position counter
//!   (reset to 624, so the first draw regenerates the state) and the rest
//!   the 624-word MT19937 state.
//! - A draw is the tempered 32-bit output times 2^-32, kept strictly inside
//!   (0, 1).
//! - `unif_index(n)`: with b = ceil(log2 n), 16-bit chunks floor(65536 u)
//!   are concatenated (most significant first) until more than b bits are
//!   held, i.e. b / 16 + 1 draws; the low b bits are the candidate, redrawn
//!   until it is below n.

/// Words of MT19937 state.
const N: usize = 624;
/// The middle offset of the recurrence.
const M: usize = 397;
/// The twist matrix's last row.
const TWIST: u32 = 0x9908_b0df;
/// 2^-32.
const TWO_POW_M32: f64 = 2.328_306_436_538_696_3e-10;
/// Half of 1 / (2^32 - 1): the nearest a draw may come to 0 or 1.
const EDGE: f64 = 0.5 / 4_294_967_295.0;

/// R's Mersenne-Twister stream after `set.seed(seed)`.
#[derive(Debug, Clone)]
pub struct RRng {
    state: [u32; N],
    /// Index of the next state word to temper; `N` means regenerate first.
    pos: usize,
}

/// One step of the seeding congruential generator.
fn lcg(s: u32) -> u32 {
    s.wrapping_mul(69069).wrapping_add(1)
}

impl RRng {
    /// The stream R gives after `set.seed(seed)`.
    pub fn set_seed(seed: i32) -> Self {
        let mut s = seed as u32;
        for _ in 0..50 {
            s = lcg(s);
        }
        // The first word R fills is its position counter, reset to N.
        s = lcg(s);
        let mut state = [0u32; N];
        for word in state.iter_mut() {
            s = lcg(s);
            *word = s;
        }
        RRng { state, pos: N }
    }

    /// Regenerate all N state words (the MT19937 recurrence).
    fn twist(&mut self) {
        for i in 0..N {
            let y = (self.state[i] & 0x8000_0000) | (self.state[(i + 1) % N] & 0x7fff_ffff);
            let mut next = self.state[(i + M) % N] ^ (y >> 1);
            if y & 1 == 1 {
                next ^= TWIST;
            }
            self.state[i] = next;
        }
        self.pos = 0;
    }

    /// The next tempered 32-bit output.
    fn next_u32(&mut self) -> u32 {
        if self.pos >= N {
            self.twist();
        }
        let mut y = self.state[self.pos];
        self.pos += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// One uniform draw on (0, 1).
    pub fn unif_rand(&mut self) -> f64 {
        let u = f64::from(self.next_u32()) * TWO_POW_M32;
        if u <= 0.0 {
            EDGE
        } else if u >= 1.0 {
            1.0 - EDGE
        } else {
            u
        }
    }

    /// `n` draws: R's `runif(n)`.
    pub fn runif(&mut self, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.unif_rand()).collect()
    }

    /// A uniform index on `0..n` by rejection sampling, as R's `sample`
    /// draws one (0-based). Panics if `n` is 0.
    pub fn unif_index(&mut self, n: usize) -> usize {
        assert!(n > 0, "unif_index: n must be positive");
        // b = ceil(log2(n)) bits.
        let bits = usize::BITS - (n - 1).leading_zeros();
        let mask: u128 = (1u128 << bits) - 1;
        loop {
            let mut acc: u128 = 0;
            let mut held = 0;
            while held <= bits {
                let chunk = (self.unif_rand() * 65536.0).floor() as u128;
                acc = (acc << 16) | chunk;
                held += 16;
            }
            let candidate = (acc & mask) as usize;
            if candidate < n {
                return candidate;
            }
        }
    }

    /// `size` independent indices on `0..n`: R's
    /// `sample.int(n, size, replace = TRUE)`, 0-based.
    pub fn sample_replace(&mut self, n: usize, size: usize) -> Vec<usize> {
        (0..size).map(|_| self.unif_index(n)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_runif_after_set_seed_matches_r() {
        // R 4.5.2: set.seed(s); sprintf("%.17g", runif(5))
        let cases: [(i32, [f64; 2]); 3] = [
            (1, [0.265_508_663_142_1, 0.372_123_899_636_790_16]),
            (42, [0.914_806_043_496_355_4, 0.937_075_413_297_861_8]),
            (-7, [0.326_494_855_573_400_86, 0.539_861_117_023_974_7]),
        ];
        for (seed, want) in cases {
            let got = RRng::set_seed(seed).runif(2);
            assert_eq!(got, want, "seed {seed}");
        }
    }

    #[test]
    fn test_state_regenerates_after_624_draws() {
        let x = RRng::set_seed(1).runif(626);
        assert_eq!(x[623], 0.136_487_586_190_924_05);
        assert_eq!(x[624], 0.324_865_140_486_508_6);
        assert_eq!(x[625], 0.621_076_293_755_322_7);
    }

    #[test]
    fn test_unif_index_draws_until_more_than_b_bits() {
        // R 4.5.2: set.seed(1); sample(n, 1, TRUE); runif(1) is the k-th
        // draw of set.seed(1); runif(4): k = 2 for n = 1, 2 (one draw used),
        // k = 3 for n = 65535, 65536, 65537 (two draws; b = 16 needs two).
        for (n, want, used) in [
            (1, 0, 1),
            (2, 0, 1),
            (65535, 24387, 2),
            (65536, 24387, 2),
            (65537, 24387, 2),
        ] {
            let mut rng = RRng::set_seed(1);
            assert_eq!(rng.unif_index(n), want, "n {n}");
            let next = rng.unif_rand();
            assert_eq!(next, RRng::set_seed(1).runif(used + 1)[used], "n {n}");
        }
    }
}
