// Copyright (c) 2008, Colorado School of Mines,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or
// without modification, are permitted provided that the following
// conditions are met:
//
//     *  Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     *  Redistributions in binary form must reproduce the above
//        copyright notice, this list of conditions and the following
//        disclaimer in the documentation and/or other materials provided
//        with the distribution.
//     *  Neither the name of the Colorado School of Mines nor the names of
//        its contributors may be used to endorse or promote products
//        derived from this software without specific prior written permission.
//
// Warranty Disclaimer:
// THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Rust port for OxWave, a Rust port of CWP/SU.

//! Single precision QR routines from cwp/lib/sqr.c, adapted from LINPACK FORTRAN.
//! Credits: Dave Hale, Colorado School of Mines, 12/29/89 (sqrdc, sqrsl)
//! and 12/31/89 (sqrst).

use oxwave_sblas::*;

// Preserve the C ABS macro, including negative zero.
fn abs(x: f32) -> f32 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

// Swap column contents, not their slice descriptors.
fn swap_columns(x: &mut [&mut [f32]], n: i32, a: usize, b: usize) {
    if a < b {
        let (left, right) = x.split_at_mut(b);
        sswap(n, left[a], 1, right[0], 1);
    } else {
        let (left, right) = x.split_at_mut(a);
        sswap(n, right[0], 1, left[b], 1);
    }
}

/// Compute a Householder QR decomposition, optionally with column pivoting.
///
/// C prototype:
/// `void sqrdc(float **x, int n, int p, float *qraux, int *jpvt, float *work, int job);`
///
/// `x[column][row]` contains at least `p` columns of at least `n` rows.
/// The upper triangle becomes R; the lower triangle and `qraux` encode Q.
/// For nonzero `job`, positive/zero/negative input `jpvt` entries mark
/// initial/free/final columns, and output entries are original column indices.
/// `qraux` and `work` need `p` elements; `jpvt` is unused when `job == 0`.
/// Requires nonnegative n,p and sufficient storage for every referenced slice.
pub fn sqrdc(
    x: &mut [&mut [f32]],
    n: i32,
    p: i32,
    qraux: &mut [f32],
    jpvt: &mut [i32],
    work: &mut [f32],
    job: i32,
) {
    let mut pl = 0i32;
    let mut pu = -1i32;
    if job != 0 {
        // In the C zero-based LINPACK adaptation, marking column zero
        // final produces -0 == 0; preserve that behavior.
        for j in 0..p {
            let ji = j as usize;
            let swapj = jpvt[ji] > 0;
            let negj = jpvt[ji] < 0;
            jpvt[ji] = j;
            if negj {
                jpvt[ji] = -j;
            }
            if swapj {
                if j != pl {
                    swap_columns(x, n, pl as usize, ji);
                }
                jpvt[ji] = jpvt[pl as usize];
                jpvt[pl as usize] = j;
                pl += 1;
            }
        }
        pu = p - 1;
        for j in (0..p).rev() {
            let ji = j as usize;
            if jpvt[ji] < 0 {
                jpvt[ji] = -jpvt[ji];
                if j != pu {
                    swap_columns(x, n, pu as usize, ji);
                    jpvt.swap(pu as usize, ji);
                }
                pu -= 1;
            }
        }
    }
    // Compute norms of free columns.
    for j in pl..=pu {
        let j = j as usize;
        qraux[j] = snrm2(n, x[j], 1);
        work[j] = qraux[j];
    }
    // Perform the Householder reduction.
    let lup = n.min(p);
    for l in 0..lup {
        let li = l as usize;
        if l >= pl && l < pu {
            let mut maxnrm = 0.0f32;
            let mut maxj = l;
            for j in l..=pu {
                if qraux[j as usize] > maxnrm {
                    maxnrm = qraux[j as usize];
                    maxj = j;
                }
            }
            if maxj != l {
                let maxj = maxj as usize;
                swap_columns(x, n, li, maxj);
                qraux[maxj] = qraux[li];
                work[maxj] = work[li];
                jpvt.swap(maxj, li);
            }
        }
        qraux[li] = 0.0;
        if l != n - 1 {
            let mut nrmxl = snrm2(n - l, &x[li][li..], 1);
            if nrmxl != 0.0 {
                if x[li][li] != 0.0 {
                    nrmxl = if x[li][li] > 0.0 {
                        abs(nrmxl)
                    } else {
                        -abs(nrmxl)
                    };
                }
                sscal(n - l, (1.0 / nrmxl as f64) as f32, &mut x[li][li..], 1);
                x[li][li] = (x[li][li] as f64 + 1.0) as f32;
                // Apply transformation to remaining columns, updating norms.
                for j in l + 1..p {
                    let ji = j as usize;
                    let mut t = -sdot(n - l, &x[li][li..], 1, &x[ji][li..], 1) / x[li][li];
                    let (left, right) = x.split_at_mut(ji);
                    saxpy(n - l, t, &left[li][li..], 1, &mut right[0][li..], 1);
                    if j >= pl && j <= pu && qraux[ji] != 0.0 {
                        let mut tt = abs(x[ji][li]) / qraux[ji];
                        tt = (1.0 - (tt * tt) as f64) as f32;
                        tt = if tt > 0.0 { tt } else { 0.0 };
                        t = tt;
                        let ttt = qraux[ji] / work[ji];
                        tt = (1.0 + 0.05 * tt as f64 * ttt as f64 * ttt as f64) as f32;
                        if tt != 1.0 {
                            qraux[ji] = (qraux[ji] as f64 * (t as f64).sqrt()) as f32;
                        } else {
                            qraux[ji] = snrm2(n - l - 1, &x[ji][li + 1..], 1);
                            work[ji] = qraux[ji];
                        }
                    }
                }
                qraux[li] = x[li][li];
                x[li][li] = -nrmxl;
            }
        }
    }
}

/// Use a QR decomposition for transforms, projections, and least squares.
///
/// C prototype:
/// `void sqrsl(float **x, int n, int k, float *qraux, float *y, float *qy, float *qty, float *b, float *rsd, float *xb, int job, int *info);`
///
/// `x[column][row]` and `qraux` come from `sqrdc`, with `0 <= k <= min(n,p)`
/// and `n >= 1`. The decimal digits ABCDE of `job` request Qy, Q'y, b,
/// residuals, and fitted values respectively; C, D, or E also requests Q'y.
/// Requested vector outputs need `n` elements, except `b`, which needs `k`
/// (at least one if the C special case `k == 0` requests b). Unrequested
/// outputs may be empty slices. Rust requires separate output storage;
/// the C routine also permits certain aliased buffers.
///
/// `info` follows the implementation: a zero diagonal reports its zero-based
/// column index, except the `min(k,n-1) == 0` special case reports 1.
/// A singular solve can leave b partly updated, as in C.
pub fn sqrsl(
    x: &mut [&mut [f32]],
    n: i32,
    k: i32,
    qraux: &[f32],
    y: &[f32],
    qy: &mut [f32],
    qty: &mut [f32],
    b: &mut [f32],
    rsd: &mut [f32],
    xb: &mut [f32],
    job: i32,
    info: &mut i32,
) {
    *info = 0;
    let cqy = job / 10000 != 0;
    let cqty = job % 10000 != 0;
    let cb = (job % 1000) / 100 != 0;
    let cr = (job % 100) / 10 != 0;
    let cxb = job % 10 != 0;
    let ju = k.min(n - 1);
    // The original takes this path for k == 0 as well as n == 1.
    if ju == 0 {
        if cqy {
            qy[0] = y[0];
        }
        if cqty {
            qty[0] = y[0];
        }
        if cxb {
            xb[0] = y[0];
        }
        if cb {
            if x[0][0] == 0.0 {
                *info = 1;
            } else {
                b[0] = y[0] / x[0][0];
            }
        }
        if cr {
            rsd[0] = 0.0;
        }
        return;
    }
    if cqy {
        scopy(n, y, 1, qy, 1);
    }
    if cqty {
        scopy(n, y, 1, qty, 1);
    }
    if cqy {
        for j in (0..ju).rev() {
            let ji = j as usize;
            if qraux[ji] != 0.0 {
                let temp = x[ji][ji];
                x[ji][ji] = qraux[ji];
                let t = -sdot(n - j, &x[ji][ji..], 1, &qy[ji..], 1) / x[ji][ji];
                saxpy(n - j, t, &x[ji][ji..], 1, &mut qy[ji..], 1);
                x[ji][ji] = temp;
            }
        }
    }
    if cqty {
        for j in 0..ju {
            let ji = j as usize;
            if qraux[ji] != 0.0 {
                let temp = x[ji][ji];
                x[ji][ji] = qraux[ji];
                let t = -sdot(n - j, &x[ji][ji..], 1, &qty[ji..], 1) / x[ji][ji];
                saxpy(n - j, t, &x[ji][ji..], 1, &mut qty[ji..], 1);
                x[ji][ji] = temp;
            }
        }
    }
    if cb {
        scopy(k, qty, 1, b, 1);
    }
    if cxb {
        scopy(k, qty, 1, xb, 1);
    }
    if cr && k < n {
        scopy(n - k, &qty[k as usize..], 1, &mut rsd[k as usize..], 1);
    }
    if cxb && k < n {
        for i in k..n {
            xb[i as usize] = 0.0;
        }
    }
    if cr {
        for i in 0..k {
            rsd[i as usize] = 0.0;
        }
    }
    if cb {
        for j in (0..k).rev() {
            let ji = j as usize;
            if x[ji][ji] == 0.0 {
                *info = j;
                break;
            }
            b[ji] /= x[ji][ji];
            if j != 0 {
                let t = -b[ji];
                saxpy(j, t, x[ji], 1, b, 1);
            }
        }
    }
    if cr || cxb {
        for j in (0..ju).rev() {
            let ji = j as usize;
            if qraux[ji] != 0.0 {
                let temp = x[ji][ji];
                x[ji][ji] = qraux[ji];
                if cr {
                    let t = -sdot(n - j, &x[ji][ji..], 1, &rsd[ji..], 1) / x[ji][ji];
                    saxpy(n - j, t, &x[ji][ji..], 1, &mut rsd[ji..], 1);
                }
                if cxb {
                    let t = -sdot(n - j, &x[ji][ji..], 1, &xb[ji..], 1) / x[ji][ji];
                    saxpy(n - j, t, &x[ji][ji..], 1, &mut xb[ji..], 1);
                }
                x[ji][ji] = temp;
            }
        }
    }
}

/// Solve an underdetermined or overdetermined least squares problem with tolerance.
///
/// C prototype:
/// `void sqrst(float **x, int n, int p, float *y, float tol, float *b, float *rsd, int *k, int *jpvt, float *qraux, float *work);`
///
/// `x[column][row]` has `p` columns of `n` rows and becomes its pivoted QR
/// decomposition. `y` and `rsd` need `n` elements; `b`, `jpvt`, `qraux`, and
/// `work` need `p`. `k` is the number of selected columns. Unused coefficients
/// become zero. The original strict diagonal tolerance test is retained:
/// even `tol == 0` rejects an exactly zero diagonal, and when no columns are
/// selected `rsd` is left untouched. Requires n,p >= 1.
pub fn sqrst(
    x: &mut [&mut [f32]],
    n: i32,
    p: i32,
    y: &[f32],
    tol: f32,
    b: &mut [f32],
    rsd: &mut [f32],
    k: &mut i32,
    jpvt: &mut [i32],
    qraux: &mut [f32],
    work: &mut [f32],
) {
    for j in 0..p {
        jpvt[j as usize] = 0;
    }
    sqrdc(x, n, p, qraux, jpvt, work, 1);
    let m = n.min(p);
    let mut kk = 0;
    for j in 0..m {
        if !(abs(x[j as usize][j as usize]) > tol * abs(x[0][0])) {
            break;
        }
        kk += 1;
    }
    if kk > 0 {
        // C aliases qty and rsd. Separate qty storage preserves the same
        // arithmetic and residual values while obeying Rust borrowing.
        let mut qty = vec![0.0; n as usize];
        let mut info = 0;
        sqrsl(
            x,
            n,
            kk,
            qraux,
            y,
            &mut [],
            &mut qty,
            work,
            rsd,
            &mut [],
            110,
            &mut info,
        );
    }
    for j in 0..kk {
        b[jpvt[j as usize] as usize] = work[j as usize];
    }
    for j in kk..p {
        b[jpvt[j as usize] as usize] = 0.0;
    }
    *k = kk;
}
