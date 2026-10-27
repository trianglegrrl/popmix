//! The limited-memory BFGS matrices: the stored corrections and the
//! products and factorizations made from them (Byrd et al. 1995, section 3).

use super::*;

/// The limited-memory BFGS matrix in compact form: S and Y (n x m, a
/// circular buffer of columns from `head`), S'Y and S'S (m x m), the
/// Cholesky factor of theta S'S + L D^-1 L' (`wt`), and the factored middle
/// matrix of the subspace problem (`wn`, with its pieces kept in `wn1`).
pub(super) struct Memory {
    pub(super) n: usize,
    pub(super) m: usize,
    pub(super) ws: Vec<f64>,
    pub(super) wy: Vec<f64>,
    pub(super) sy: Vec<f64>,
    pub(super) ss: Vec<f64>,
    pub(super) wt: Vec<f64>,
    pub(super) wn: Vec<f64>,
    pub(super) wn1: Vec<f64>,
    pub(super) theta: f64,
    pub(super) col: usize,
    pub(super) head: usize,
    pub(super) itail: usize,
    pub(super) iupdat: usize,
    pub(super) updatd: bool,
}

impl Memory {
    pub(super) fn new(n: usize, m: usize) -> Self {
        Memory {
            n,
            m,
            ws: vec![0.0; n * m],
            wy: vec![0.0; n * m],
            sy: vec![0.0; m * m],
            ss: vec![0.0; m * m],
            wt: vec![0.0; m * m],
            wn: vec![0.0; 4 * m * m],
            wn1: vec![0.0; 4 * m * m],
            theta: 1.0,
            col: 0,
            head: 1,
            itail: 0,
            iupdat: 0,
            updatd: false,
        }
    }

    /// Forget the corrections (after a numerical failure).
    pub(super) fn reset(&mut self) {
        self.col = 0;
        self.head = 1;
        self.theta = 1.0;
        self.iupdat = 0;
        self.updatd = false;
    }

    #[inline]
    pub(super) fn wy(&self, i: usize, j: usize) -> f64 {
        self.wy[at(i, j, self.n)]
    }

    #[inline]
    pub(super) fn ws(&self, i: usize, j: usize) -> f64 {
        self.ws[at(i, j, self.n)]
    }
}

/// p = M v for the 2col x 2col middle matrix M of the compact form, using
/// S'Y and the factor in `wt`. Returns a nonzero info when `wt` is singular.
pub(super) fn bmv(mem: &Memory, v: &[f64], p: &mut [f64]) -> usize {
    let (m, col) = (mem.m, mem.col);
    if col == 0 {
        return 0;
    }
    let sy = |i: usize, j: usize| mem.sy[at(i, j, m)];
    p[col] = v[col];
    for i in 2..=col {
        let i2 = col + i;
        let mut sum = 0.0;
        for k in 1..i {
            sum += sy(i, k) * v[k - 1] / sy(k, k);
        }
        p[i2 - 1] = v[i2 - 1] + sum;
    }
    let info = dtrsl(&mem.wt, 0, m, col, &mut p[col..2 * col], true);
    if info != 0 {
        return info;
    }
    for i in 1..=col {
        p[i - 1] = v[i - 1] / sy(i, i).sqrt();
    }
    let info = dtrsl(&mem.wt, 0, m, col, &mut p[col..2 * col], false);
    if info != 0 {
        return info;
    }
    for i in 1..=col {
        p[i - 1] = -p[i - 1] / sy(i, i).sqrt();
    }
    for i in 1..=col {
        let mut sum = 0.0;
        for k in (i + 1)..=col {
            sum += sy(k, i) * p[col + k - 1] / sy(i, i);
        }
        p[i - 1] += sum;
    }
    0
}

/// Form and factor the middle matrix K of the subspace problem,
///   K = [-D - Y'ZZ'Y/theta     L_a' - R_z'  ]
///       [L_a - R_z            theta S'AA'S ],
/// updating its pieces in `wn1` for the variables that entered or left the
/// free set. Returns a nonzero info if a factorisation fails.
pub(super) fn formk(mem: &mut Memory, fs: &FreeSet) -> i32 {
    let (n, m, col, head, theta) = (mem.n, mem.m, mem.col, mem.head, mem.theta);
    let (nsub, nenter, ileave) = (fs.nfree, fs.nenter, fs.ileave);
    let (ind, indx2) = (&fs.index, &fs.indx2);
    let m2 = 2 * m;
    let w1 = |i: usize, j: usize| at(i, j, m2);
    let wyv = |mem: &Memory, i: usize, j: usize| mem.wy[at(i, j, n)];
    let wsv = |mem: &Memory, i: usize, j: usize| mem.ws[at(i, j, n)];

    let upcl = if mem.updatd {
        if mem.iupdat > m {
            // shift the old part of wn1 up-left
            for jy in 1..m {
                let js = m + jy;
                for k in 0..(m - jy) {
                    mem.wn1[w1(jy + k, jy)] = mem.wn1[w1(jy + 1 + k, jy + 1)];
                }
                for k in 0..(m - jy) {
                    mem.wn1[w1(js + k, js)] = mem.wn1[w1(js + 1 + k, js + 1)];
                }
                for k in 0..(m - 1) {
                    mem.wn1[w1(m + 1 + k, jy)] = mem.wn1[w1(m + 2 + k, jy + 1)];
                }
            }
        }
        // new row and column of the blocks from the free (pbegin..pend)
        // and active (dbegin..dend) variables
        let (pbegin, pend, dbegin, dend) = (1, nsub, nsub + 1, n);
        let iy = col;
        let is = m + col;
        let mut ipntr = head + col - 1;
        if ipntr > m {
            ipntr -= m;
        }
        let mut jpntr = head;
        for jy in 1..=col {
            let js = m + jy;
            let mut temp1 = 0.0;
            let mut temp2 = 0.0;
            let mut temp3 = 0.0;
            for k in pbegin..=pend {
                let k1 = ind[k - 1];
                temp1 += wyv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            for k in dbegin..=dend {
                let k1 = ind[k - 1];
                temp2 += wsv(mem, k1, ipntr) * wsv(mem, k1, jpntr);
                temp3 += wsv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            mem.wn1[w1(iy, jy)] = temp1;
            mem.wn1[w1(is, js)] = temp2;
            mem.wn1[w1(is, jy)] = temp3;
            jpntr = next(jpntr, m);
        }
        let jy = col;
        let mut jpntr = head + col - 1;
        if jpntr > m {
            jpntr -= m;
        }
        let mut ipntr = head;
        for i in 1..=col {
            let is = m + i;
            let mut temp3 = 0.0;
            for k in pbegin..=pend {
                let k1 = ind[k - 1];
                temp3 += wsv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            ipntr = next(ipntr, m);
            mem.wn1[w1(is, jy)] = temp3;
        }
        col - 1
    } else {
        col
    };

    // the old part, for the variables that entered or left
    let mut ipntr = head;
    for iy in 1..=upcl {
        let is = m + iy;
        let mut jpntr = head;
        for jy in 1..=iy {
            let js = m + jy;
            let mut temp1 = 0.0;
            let mut temp2 = 0.0;
            let mut temp3 = 0.0;
            let mut temp4 = 0.0;
            for k in 1..=nenter {
                let k1 = indx2[k - 1];
                temp1 += wyv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
                temp2 += wsv(mem, k1, ipntr) * wsv(mem, k1, jpntr);
            }
            for k in ileave..=n {
                let k1 = indx2[k - 1];
                temp3 += wyv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
                temp4 += wsv(mem, k1, ipntr) * wsv(mem, k1, jpntr);
            }
            mem.wn1[w1(iy, jy)] = mem.wn1[w1(iy, jy)] + temp1 - temp3;
            mem.wn1[w1(is, js)] = mem.wn1[w1(is, js)] - temp2 + temp4;
            jpntr = next(jpntr, m);
        }
        ipntr = next(ipntr, m);
    }
    let mut ipntr = head;
    for is in (m + 1)..=(m + upcl) {
        let mut jpntr = head;
        for jy in 1..=upcl {
            let mut temp1 = 0.0;
            let mut temp3 = 0.0;
            for k in 1..=nenter {
                let k1 = indx2[k - 1];
                temp1 += wsv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            for k in ileave..=n {
                let k1 = indx2[k - 1];
                temp3 += wsv(mem, k1, ipntr) * wyv(mem, k1, jpntr);
            }
            if is <= jy + m {
                mem.wn1[w1(is, jy)] = mem.wn1[w1(is, jy)] + temp1 - temp3;
            } else {
                mem.wn1[w1(is, jy)] = mem.wn1[w1(is, jy)] - temp1 + temp3;
            }
            jpntr = next(jpntr, m);
        }
        ipntr = next(ipntr, m);
    }

    // the upper triangle of K, in wn
    for iy in 1..=col {
        let is = col + iy;
        let is1 = m + iy;
        for jy in 1..=iy {
            let js = col + jy;
            let js1 = m + jy;
            mem.wn[w1(jy, iy)] = mem.wn1[w1(iy, jy)] / theta;
            mem.wn[w1(js, is)] = mem.wn1[w1(is1, js1)] * theta;
        }
        for jy in 1..iy {
            mem.wn[w1(jy, is)] = -mem.wn1[w1(is1, jy)];
        }
        for jy in iy..=col {
            mem.wn[w1(jy, is)] = mem.wn1[w1(is1, jy)];
        }
        mem.wn[w1(iy, iy)] += mem.sy[at(iy, iy, m)];
    }

    // LEL' factorisation: factor the (1,1) block, solve for the (1,2)
    // block, update and factor the (2,2) block
    if dpofa(&mut mem.wn, 0, m2, col) != 0 {
        return -1;
    }
    let col2 = 2 * col;
    let mut b = vec![0.0; col];
    for js in (col + 1)..=col2 {
        for i in 1..=col {
            b[i - 1] = mem.wn[w1(i, js)];
        }
        dtrsl(&mem.wn, 0, m2, col, &mut b, true);
        for i in 1..=col {
            mem.wn[w1(i, js)] = b[i - 1];
        }
    }
    for is in (col + 1)..=col2 {
        for js in is..=col2 {
            let mut dot = 0.0;
            for k in 1..=col {
                dot += mem.wn[w1(k, is)] * mem.wn[w1(k, js)];
            }
            mem.wn[w1(is, js)] += dot;
        }
    }
    if dpofa(&mut mem.wn, w1(col + 1, col + 1), m2, col) != 0 {
        return -2;
    }
    0
}

/// Form and factor T = theta S'S + L D^-1 L' (upper triangle) into `wt`.
pub(super) fn formt(mem: &mut Memory) -> i32 {
    let (m, col, theta) = (mem.m, mem.col, mem.theta);
    let sy = |mem: &Memory, i: usize, j: usize| mem.sy[at(i, j, m)];
    for j in 1..=col {
        mem.wt[at(1, j, m)] = theta * mem.ss[at(1, j, m)];
    }
    for i in 2..=col {
        for j in i..=col {
            let k1 = i.min(j) - 1;
            let mut ddum = 0.0;
            for k in 1..=k1 {
                ddum += sy(mem, i, k) * sy(mem, j, k) / sy(mem, k, k);
            }
            mem.wt[at(i, j, m)] = ddum + theta * mem.ss[at(i, j, m)];
        }
    }
    if dpofa(&mut mem.wt, 0, m, col) != 0 {
        return -3;
    }
    0
}

/// Add the correction pair (s, y) = (d, r) to the memory.
pub(super) fn matupd(mem: &mut Memory, d: &[f64], r: &[f64], rr: f64, dr: f64, stp: f64, dtd: f64) {
    let (n, m) = (mem.n, mem.m);
    if mem.iupdat <= m {
        mem.col = mem.iupdat;
        mem.itail = (mem.head + mem.iupdat - 2) % m + 1;
    } else {
        mem.itail = next(mem.itail, m);
        mem.head = next(mem.head, m);
    }
    let it = mem.itail;
    mem.ws[at(1, it, n)..at(1, it, n) + n].copy_from_slice(d);
    mem.wy[at(1, it, n)..at(1, it, n) + n].copy_from_slice(r);
    mem.theta = rr / dr;
    let col = mem.col;
    if mem.iupdat > m {
        // drop the oldest pair from S'S and S'Y
        for j in 1..col {
            for k in 0..j {
                mem.ss[at(1 + k, j, m)] = mem.ss[at(2 + k, j + 1, m)];
            }
            for k in 0..(col - j) {
                mem.sy[at(j + k, j, m)] = mem.sy[at(j + 1 + k, j + 1, m)];
            }
        }
    }
    let mut pointr = mem.head;
    for j in 1..col {
        let wy = &mem.wy[at(1, pointr, n)..at(1, pointr, n) + n];
        let ws = &mem.ws[at(1, pointr, n)..at(1, pointr, n) + n];
        let syv = ddot(d, wy);
        let ssv = ddot(ws, d);
        mem.sy[at(col, j, m)] = syv;
        mem.ss[at(j, col, m)] = ssv;
        pointr = next(pointr, m);
    }
    mem.ss[at(col, col, m)] = if stp == 1.0 { dtd } else { stp * stp * dtd };
    mem.sy[at(col, col, m)] = dr;
}
