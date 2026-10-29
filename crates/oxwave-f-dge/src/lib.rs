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

//! Port of cwp/lib/dge.c, adapted from LINPACK FORTRAN.
//! Author: Dave Hale, Colorado School of Mines, 10/01/89.

use oxwave_f_dblas::{dasum, daxpy, ddot, dscal, idamax};

fn abs(x: f64) -> f64 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

/// Factors a positive-order matrix given as a slice of columns.
/// C prototype: `void dgefa(double **a, int n, int *ipvt, int *info);`
pub fn dgefa(a: &mut [&mut [f64]], n: i32, ipvt: &mut [i32], info: &mut i32) {
    let n = n as usize;
    *info = -1;
    for k in 0..n - 1 {
        let l = k + idamax((n - k) as i32, &a[k][k..], 1) as usize;
        ipvt[k] = l as i32;
        if a[k][l] == 0.0 {
            *info = k as i32;
            continue;
        }
        if l != k {
            a[k].swap(l, k);
        }
        let t = -1.0 / a[k][k];
        dscal((n - k - 1) as i32, t, &mut a[k][k + 1..], 1);
        for j in k + 1..n {
            let t = a[j][l];
            if l != k {
                a[j][l] = a[j][k];
                a[j][k] = t;
            }
            let (left, right) = a.split_at_mut(j);
            daxpy(
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

/// Factors a positive-order column-oriented matrix and estimates its reciprocal condition number.
/// C prototype: `void dgeco(double **a, int n, int *ipvt, double *rcond, double *z);`
pub fn dgeco(a: &mut [&mut [f64]], n: i32, ipvt: &mut [i32], rcond: &mut f64, z: &mut [f64]) {
    let nn = n as usize;
    let mut anorm = 0.0;
    for col in &a[..nn] {
        let t = dasum(n, col, 1);
        anorm = if t > anorm { t } else { anorm };
    }
    let mut info = 0;
    dgefa(a, n, ipvt, &mut info);
    let mut ek = 1.0;
    for v in &mut z[..nn] {
        *v = 0.0;
    }
    for k in 0..nn {
        if z[k] != 0.0 {
            ek = if z[k] > 0.0 { -abs(ek) } else { abs(ek) };
        }
        if abs(ek - z[k]) > abs(a[k][k]) {
            let s = abs(a[k][k]) / abs(ek - z[k]);
            dscal(n, s, z, 1);
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
        if k + 1 < nn {
            for j in k + 1..nn {
                let t = z[j] + wkm * a[j][k];
                sm += abs(t);
                z[j] += wk * a[j][k];
                s += abs(z[j]);
            }
            if s < sm {
                let t = wkm - wk;
                wk = wkm;
                for j in k + 1..nn {
                    z[j] += t * a[j][k];
                }
            }
        }
        z[k] = wk;
    }
    let mut s = 1.0 / dasum(n, z, 1);
    dscal(n, s, z, 1);
    for k in (0..nn).rev() {
        if k < nn - 1 {
            z[k] += ddot((nn - k - 1) as i32, &a[k][k + 1..], 1, &z[k + 1..], 1);
        }
        if abs(z[k]) > 1.0 {
            s = 1.0 / abs(z[k]);
            dscal(n, s, z, 1);
        }
        let l = ipvt[k] as usize;
        z.swap(l, k);
    }
    s = 1.0 / dasum(n, z, 1);
    dscal(n, s, z, 1);
    let mut ynorm = 1.0;
    for k in 0..nn {
        let l = ipvt[k] as usize;
        let t = z[l];
        z.swap(l, k);
        if k < nn - 1 {
            daxpy(
                (nn - k - 1) as i32,
                t,
                &a[k][k + 1..],
                1,
                &mut z[k + 1..],
                1,
            );
        }
        if abs(z[k]) > 1.0 {
            s = 1.0 / abs(z[k]);
            dscal(n, s, z, 1);
            ynorm *= s;
        }
    }
    s = 1.0 / dasum(n, z, 1);
    dscal(n, s, z, 1);
    ynorm *= s;
    for k in (0..nn).rev() {
        if abs(z[k]) > abs(a[k][k]) {
            s = abs(a[k][k]) / abs(z[k]);
            dscal(n, s, z, 1);
            ynorm *= s;
        }
        if a[k][k] != 0.0 {
            z[k] /= a[k][k];
        } else {
            z[k] = 1.0;
        }
        let t = -z[k];
        daxpy(k as i32, t, a[k], 1, z, 1);
    }
    s = 1.0 / dasum(n, z, 1);
    dscal(n, s, z, 1);
    ynorm *= s;
    *rcond = if anorm != 0.0 { ynorm / anorm } else { 0.0 };
}

/// Solves a factored positive-order system, or its transpose when `job != 0`.
/// C prototype: `void dgesl(double **a, int n, int *ipvt, double *b, int job);`
pub fn dgesl(a: &[&[f64]], n: i32, ipvt: &[i32], b: &mut [f64], job: i32) {
    let n = n as usize;
    if job == 0 {
        for k in 0..n - 1 {
            let l = ipvt[k] as usize;
            let t = b[l];
            if l != k {
                b[l] = b[k];
                b[k] = t;
            }
            daxpy((n - k - 1) as i32, t, &a[k][k + 1..], 1, &mut b[k + 1..], 1);
        }
        for k in (0..n).rev() {
            b[k] /= a[k][k];
            let t = -b[k];
            daxpy(k as i32, t, a[k], 1, b, 1);
        }
    } else {
        for k in 0..n {
            let t = ddot(k as i32, a[k], 1, b, 1);
            b[k] = (b[k] - t) / a[k][k];
        }
        for k in (0..n - 1).rev() {
            b[k] += ddot((n - k - 1) as i32, &a[k][k + 1..], 1, &b[k + 1..], 1);
            let l = ipvt[k] as usize;
            if l != k {
                b.swap(l, k);
            }
        }
    }
}
