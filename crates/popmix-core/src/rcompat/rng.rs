//! R's default random number generator: Mersenne Twister (MT19937) seeded
//! by `set.seed`, with `runif`'s mapping to (0, 1).
//!
//! `set.seed(seed)` scrambles the seed with 50 steps of the linear
//! congruential generator seed = 69069 seed + 1, then fills R's 625-word
//! state (the twister's position followed by its 624 words) with further
//! steps of it, and sets the position to 624 so the first draw regenerates
//! the state. Each uniform is a tempered 32-bit word times 2^-32, moved
//! into the open interval (0, 1) if it falls on 0 or 1. `sample` draws
//! indices by rejection from as many random bits as the range needs, taken
//! 16 at a time from uniforms (R 3.6.0 and later).

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;
const I2_32M1: f64 = 2.328_306_437_080_797e-10; // 1 / (2^32 - 1)

/// R's Mersenne Twister, as after `set.seed(seed)`.
#[derive(Debug, Clone)]
pub struct RRng {
    mt: [u32; N],
    mti: usize,
}

impl RRng {
    /// The generator after R's `set.seed(seed)` (default kind).
    pub fn set_seed(seed: i32) -> Self {
        let mut s = seed as u32;
        for _ in 0..50 {
            s = s.wrapping_mul(69069).wrapping_add(1);
        }
        // R's first state word is the position; it is overwritten with 624
        s = s.wrapping_mul(69069).wrapping_add(1);
        let mut mt = [0u32; N];
        for w in mt.iter_mut() {
            s = s.wrapping_mul(69069).wrapping_add(1);
            *w = s;
        }
        Self { mt, mti: N }
    }

    fn genrand(&mut self) -> f64 {
        let mag01 = [0u32, MATRIX_A];
        if self.mti >= N {
            let mt = &mut self.mt;
            for kk in 0..N - M {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            }
            for kk in N - M..N - 1 {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M - N] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            }
            let y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            self.mti = 0;
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y as f64 * 2.328_306_436_538_696_3e-10 // 2^-32
    }

    /// One draw of R's `unif_rand()`.
    pub fn unif_rand(&mut self) -> f64 {
        let v = self.genrand();
        // R's fixup: keep the value in the open interval (0, 1)
        if v <= 0.0 {
            0.5 * I2_32M1
        } else if 1.0 - v <= 0.0 {
            1.0 - 0.5 * I2_32M1
        } else {
            v
        }
    }

    /// R's `runif(n)`.
    pub fn runif(&mut self, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.unif_rand()).collect()
    }

    /// `bits` random bits from 16-bit chunks of uniforms (R's `rbits`).
    fn rbits(&mut self, bits: u32) -> f64 {
        let mut v: i64 = 0;
        let mut n = 0;
        while n <= bits {
            let v1 = (self.unif_rand() * 65536.0).floor() as i64;
            v = 65536 * v + v1;
            n += 16;
        }
        (v & ((1i64 << bits) - 1)) as f64
    }

    /// A uniform index in 0..n (R's `R_unif_index`, the default "Rejection"
    /// sample kind since R 3.6.0).
    pub fn unif_index(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        let dn = n as f64;
        let bits = dn.log2().ceil() as u32;
        loop {
            let dv = self.rbits(bits);
            if dn > dv {
                return dv as usize;
            }
        }
    }

    /// R's `sample.int(n, size, replace = TRUE)`, 0-based.
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
}
