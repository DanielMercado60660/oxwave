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

//! Single precision LU routines ported from `cwp/lib/sge.c` (LINPACK FORTRAN adaptations).
//! Credits: Dave Hale, Colorado School of Mines, 10/01/89.
//!
//! Matrices use the C layout `a[column][row]`. Each column must contain at
//! least `n` elements, and `n` must be positive, as required by the C routines.
//! Pivot indices are zero based. Extra columns and trailing slice elements
//! are left untouched.

use oxwave_sblas::*;

// C's ABS macro preserves negative zero and the sign of NaNs.
fn abs(x: f32) -> f32 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

/// Gaussian elimination with partial pivoting, overwriting `a` with LU factors.
/// `info` is the last zero pivot index, or -1 if all pivots are nonzero.
///
/// C prototype: `void sgefa (float **a, int n, int *ipvt, int *info);`
pub fn sgefa(a: &mut [&mut [f32]], n: i32, ipvt: &mut [i32], info: &mut i32) {
    assert!(n > 0, "sgefa requires n >= 1");
    let n = n as usize;
    *info = -1;
    let nm1 = n - 1;
    for k in 0..nm1 {
        let kp1 = k + 1;

        // Find l = pivot index.
        let l = k + isamax((n - k) as i32, &a[k][k..], 1) as usize;
        ipvt[k] = l as i32;

        // A zero pivot implies this column is already triangularized.
        if a[k][l] == 0.0 {
            *info = k as i32;
            continue;
        }

        if l != k {
            let t = a[k][l];
            a[k][l] = a[k][k];
            a[k][k] = t;
        }

        // The unsuffixed -1.0 in C makes the division double precision.
        let t = (-1.0f64 / a[k][k] as f64) as f32;
        sscal((n - k - 1) as i32, t, &mut a[k][k + 1..], 1);

        // Row elimination with column indexing.
        for j in kp1..n {
            let t = a[j][l];
            if l != k {
                a[j][l] = a[j][k];
                a[j][k] = t;
            }
            let (left, right) = a.split_at_mut(j);
            saxpy(
                (n - k - 1) as i32,
                t,
                &left[k][k + 1..],
                1,
                &mut right[0][k + 1..],
                1,
            );
        }
    }
    ipvt[n - 1] = (n - 1) as i32;
    if a[n - 1][n - 1] == 0.0 {
        *info = (n - 1) as i32;
    }
}

/// Factor `a` and estimate its reciprocal 1-norm condition number in `rcond`.
/// `z` is workspace and receives the C estimator's normalized vector.
///
/// C prototype: `void sgeco (float **a, int n, int *ipvt, float *rcond, float *z);`
pub fn sgeco(a: &mut [&mut [f32]], n: i32, ipvt: &mut [i32], rcond: &mut f32, z: &mut [f32]) {
    assert!(n > 0, "sgeco requires n >= 1");
    let ni = n;
    let n = n as usize;

    // Compute the 1-norm of a.
    let mut anorm = 0.0f32;
    for j in 0..n {
        let t = sasum(ni, a[j], 1);
        anorm = if t > anorm { t } else { anorm };
    }

    let mut info = 0;
    sgefa(a, ni, ipvt, &mut info);

    // Solve U'w = e, choosing e to maximize local growth. Rescale often
    // to avoid overflow, following LINPACK's condition estimator.
    let mut ek = 1.0f32;
    for j in 0..n {
        z[j] = 0.0;
    }
    for k in 0..n {
        if z[k] != 0.0 {
            ek = if z[k] > 0.0 { -abs(ek) } else { abs(ek) };
        }
        if abs(ek - z[k]) > abs(a[k][k]) {
            let s = abs(a[k][k]) / abs(ek - z[k]);
            sscal(ni, s, z, 1);
            ek *= s;
        }
        let mut wk = ek - z[k];
        let mut wkm = -ek - z[k];
        let mut s = abs(wk);
        let mut sm = abs(wkm);
        if a[k][k] == 0.0 {
            wk = 1.0;
            wkm = 1.0;
        } else {
            wk /= a[k][k];
            wkm /= a[k][k];
        }
        let kp1 = k + 1;
        if kp1 < n {
            for j in kp1..n {
                let t = z[j] + wkm * a[j][k];
                sm += abs(t);
                z[j] += wk * a[j][k];
                s += abs(z[j]);
            }
            if s < sm {
                let t = wkm - wk;
                wk = wkm;
                for j in kp1..n {
                    z[j] += t * a[j][k];
                }
            }
        }
        z[k] = wk;
    }
    let s = (1.0f64 / sasum(ni, z, 1) as f64) as f32;
    sscal(ni, s, z, 1);

    // Solve L'y = w.
    for k in (0..n).rev() {
        if k < n - 1 {
            z[k] += sdot((n - k - 1) as i32, &a[k][k + 1..], 1, &z[k + 1..], 1);
        }
        if abs(z[k]) as f64 > 1.0 {
            let s = (1.0f64 / abs(z[k]) as f64) as f32;
            sscal(ni, s, z, 1);
        }
        let l = ipvt[k] as usize;
        let t = z[l];
        z[l] = z[k];
        z[k] = t;
    }
    let s = (1.0f64 / sasum(ni, z, 1) as f64) as f32;
    sscal(ni, s, z, 1);

    let mut ynorm = 1.0f32;

    // Solve Lv = y.
    for k in 0..n {
        let l = ipvt[k] as usize;
        let t = z[l];
        z[l] = z[k];
        z[k] = t;
        if k < n - 1 {
            saxpy((n - k - 1) as i32, t, &a[k][k + 1..], 1, &mut z[k + 1..], 1);
        }
        if abs(z[k]) as f64 > 1.0 {
            let s = (1.0f64 / abs(z[k]) as f64) as f32;
            sscal(ni, s, z, 1);
            ynorm *= s;
        }
    }
    let s = (1.0f64 / sasum(ni, z, 1) as f64) as f32;
    sscal(ni, s, z, 1);
    ynorm *= s;

    // Solve Uz = v.
    for k in (0..n).rev() {
        if abs(z[k]) > abs(a[k][k]) {
            let s = abs(a[k][k]) / abs(z[k]);
            sscal(ni, s, z, 1);
            ynorm *= s;
        }
        if a[k][k] != 0.0 {
            z[k] /= a[k][k];
        } else {
            z[k] = 1.0;
        }
        let t = -z[k];
        saxpy(k as i32, t, a[k], 1, z, 1);
    }

    // Make znorm = 1.0.
    let s = (1.0f64 / sasum(ni, z, 1) as f64) as f32;
    sscal(ni, s, z, 1);
    ynorm *= s;

    *rcond = if anorm != 0.0 { ynorm / anorm } else { 0.0 };
}

/// Solve Ax = b (`job == 0`) or A'x = b (any other `job`) using the
/// factors and pivots produced by `sgefa` or `sgeco`, overwriting `b`.
/// As in C, zero pivots may produce infinities or NaNs.
///
/// C prototype: `void sgesl (float **a, int n, int *ipvt, float *b, int job);`
pub fn sgesl(a: &[&[f32]], n: i32, ipvt: &[i32], b: &mut [f32], job: i32) {
    assert!(n > 0, "sgesl requires n >= 1");
    let n = n as usize;
    let nm1 = n - 1;
    if job == 0 {
        // First solve Ly = b.
        for k in 0..nm1 {
            let l = ipvt[k] as usize;
            let t = b[l];
            if l != k {
                b[l] = b[k];
                b[k] = t;
            }
            saxpy((n - k - 1) as i32, t, &a[k][k + 1..], 1, &mut b[k + 1..], 1);
        }

        // Now solve Ux = y.
        for k in (0..n).rev() {
            b[k] /= a[k][k];
            let t = -b[k];
            saxpy(k as i32, t, a[k], 1, b, 1);
        }
    } else {
        // First solve U'y = b.
        for k in 0..n {
            let t = sdot(k as i32, a[k], 1, b, 1);
            b[k] = (b[k] - t) / a[k][k];
        }

        // Now solve L'x = y.
        for k in (0..n - 1).rev() {
            b[k] += sdot((n - k - 1) as i32, &a[k][k + 1..], 1, &b[k + 1..], 1);
            let l = ipvt[k] as usize;
            if l != k {
                let t = b[l];
                b[l] = b[k];
                b[k] = t;
            }
        }
    }
}
