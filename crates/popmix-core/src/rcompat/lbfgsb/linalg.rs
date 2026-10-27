//! Dense linear algebra L-BFGS-B needs: column-major indexing, dot
//! products, LINPACK's triangular solve and Cholesky factorization.

#[inline]
pub(super) fn at(i: usize, j: usize, ld: usize) -> usize {
    (i - 1) + (j - 1) * ld
}

#[inline]
pub(super) fn next(p: usize, m: usize) -> usize {
    p % m + 1
}

pub(super) fn ddot(x: &[f64], y: &[f64]) -> f64 {
    let mut s = 0.0;
    for (a, b) in x.iter().zip(y) {
        s += a * b;
    }
    s
}

/// Solve with an upper triangular matrix `t` (column-major, leading
/// dimension `ldt`, its (1, 1) at offset `off`), of order `n`: T x = b
/// (`transpose` false) or T' x = b. Returns the index of a zero diagonal
/// element, or 0.
pub(super) fn dtrsl(
    t: &[f64],
    off: usize,
    ldt: usize,
    n: usize,
    b: &mut [f64],
    transpose: bool,
) -> usize {
    let tt = |i: usize, j: usize| t[off + at(i, j, ldt)];
    for info in 1..=n {
        if tt(info, info) == 0.0 {
            return info;
        }
    }
    if n == 0 {
        return 0;
    }
    if transpose {
        b[0] /= tt(1, 1);
        for j in 2..=n {
            let mut s = 0.0;
            for k in 1..j {
                s += tt(k, j) * b[k - 1];
            }
            b[j - 1] -= s;
            b[j - 1] /= tt(j, j);
        }
    } else {
        b[n - 1] /= tt(n, n);
        for jj in 2..=n {
            let j = n - jj + 1;
            let temp = -b[j];
            for i in 1..=j {
                b[i - 1] += temp * tt(i, j + 1);
            }
            b[j - 1] /= tt(j, j);
        }
    }
    0
}

/// Cholesky factorisation A = R'R of the upper triangle of a symmetric
/// matrix (order `n`, at offset `off`, leading dimension `lda`), in place.
/// Returns 0, or the order of the leading minor that is not positive
/// definite.
pub(super) fn dpofa(a: &mut [f64], off: usize, lda: usize, n: usize) -> usize {
    let ix = |i: usize, j: usize| off + at(i, j, lda);
    for j in 1..=n {
        let mut s = 0.0;
        for k in 1..j {
            let mut dot = 0.0;
            for i in 1..k {
                dot += a[ix(i, k)] * a[ix(i, j)];
            }
            let mut t = a[ix(k, j)] - dot;
            t /= a[ix(k, k)];
            a[ix(k, j)] = t;
            s += t * t;
        }
        s = a[ix(j, j)] - s;
        if s <= 0.0 {
            return j;
        }
        a[ix(j, j)] = s.sqrt();
    }
    0
}
