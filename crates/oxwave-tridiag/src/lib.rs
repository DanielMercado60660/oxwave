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

//! Tridiagonal solvers ported from `cwp/lib/tridiagonal.c`.
//! Credits: tridif/tridid: Dave Hale, Colorado School of Mines, 10/03/89;
//! tripd/tripp: Zhenyue Liu, Colorado School of Mines, 1992-1993.

/// Solve a tridiagonal system in double precision without pivoting.
///
/// C prototype:
/// `void tridid(int n, double a[], double b[], double c[], double r[], double u[]);`
///
/// The dimension is `b.len()`. All other slices must contain at least that
/// many elements; `a[0]` and `c[n-1]` are ignored. Overwrites the first `n`
/// elements of `u`. Zero elimination pivots propagate IEEE infinities/NaNs.
/// Panics if the system is empty or a slice is too short.
pub fn tridid(a: &[f64], b: &[f64], c: &[f64], r: &[f64], u: &mut [f64]) {
    let n = b.len();
    assert!(n > 0 && a.len() >= n && c.len() >= n && r.len() >= n && u.len() >= n);
    let mut w = vec![0.0_f64; n];
    let mut t = b[0];
    u[0] = r[0] / t;
    for j in 1..n {
        w[j] = c[j - 1] / t;
        t = b[j] - a[j] * w[j];
        u[j] = (r[j] - a[j] * u[j - 1]) / t;
    }
    for j in (0..n - 1).rev() {
        u[j] -= w[j + 1] * u[j + 1];
    }
}

/// Solve a tridiagonal system in single precision without pivoting.
///
/// C prototype:
/// `void tridif(int n, float a[], float b[], float c[], float r[], float u[]);`
///
/// The dimension is `b.len()`. All other slices must contain at least that
/// many elements; `a[0]` and `c[n-1]` are ignored. Overwrites the first `n`
/// elements of `u`. Zero elimination pivots propagate IEEE infinities/NaNs.
/// Panics if the system is empty or a slice is too short.
pub fn tridif(a: &[f32], b: &[f32], c: &[f32], r: &[f32], u: &mut [f32]) {
    let n = b.len();
    assert!(n > 0 && a.len() >= n && c.len() >= n && r.len() >= n && u.len() >= n);
    let mut w = vec![0.0_f32; n];
    let mut t = b[0];
    u[0] = r[0] / t;
    for j in 1..n {
        w[j] = c[j - 1] / t;
        t = b[j] - a[j] * w[j];
        u[j] = (r[j] - a[j] * u[j - 1]) / t;
    }
    for j in (0..n - 1).rev() {
        u[j] -= w[j + 1] * u[j + 1];
    }
}

/// Solve a symmetric positive definite tridiagonal system in place.
///
/// C prototype: `void tripd(float *d, float *e, float *b, int n);`
///
/// The dimension is `d.len()`. `e` must contain at least `n-1` elements and
/// `b` at least `n`. The decomposition overwrites `d` and the first `n-1`
/// elements of `e`; the first `n` elements of `b` become the solution.
/// Panics if the system is empty or a slice is too short.
pub fn tripd(d: &mut [f32], e: &mut [f32], b: &mut [f32]) {
    let n = d.len();
    assert!(n > 0 && e.len() >= n - 1 && b.len() >= n);

    // Decomposition.
    for k in 1..n {
        let temp = e[k - 1];
        e[k - 1] = temp / d[k - 1];
        d[k] -= temp * e[k - 1];
    }

    // Substitution.
    for k in 1..n {
        b[k] -= e[k - 1] * b[k - 1];
    }
    b[n - 1] /= d[n - 1];
    for k in (1..n).rev() {
        b[k - 1] = b[k - 1] / d[k - 1] - e[k - 1] * b[k];
    }
}

// The original exch takes both arguments BY VALUE. Its assignments affect
// only local copies. Retain that no-op rather than fixing the C pivoting bug.
fn exch(_x: f32, _y: f32) {}

/// Solve an unsymmetric tridiagonal system in place, retaining C's behavior.
///
/// C prototype: `void tripp(int n, float *d, float *e, float *c, float *b);`
///
/// The dimension is `d.len()`. All other slices must contain at least `n`
/// elements. `c[0]` is ignored on input; the subdiagonal starts at `c[1]`.
/// `d`, `e`, and `c` are overwritten during elimination, and the first `n`
/// elements of `b` become the solution. The C pivot helper is a no-op, which
/// is preserved here. A zero pivot writes `coefficient matrix is singular!`
/// to stderr and calculation continues with IEEE infinities/NaNs.
///
/// For `n == 1`, only `b[0] /= d[0]` is performed: C additionally accesses
/// index -1 of `b`, `d`, and `e`, outside its documented arrays.
/// Panics if the system is empty or a slice is too short.
pub fn tripp(d: &mut [f32], e: &mut [f32], c: &mut [f32], b: &mut [f32]) {
    let n = d.len();
    assert!(n > 0 && e.len() >= n && c.len() >= n && b.len() >= n);

    // Elimination.
    for k in 0..n - 1 {
        c[k] = 0.0;
        // ABS in cwp.h is a comparison and conditional negation.
        let abs_d = if d[k] < 0.0 { -d[k] } else { d[k] };
        let abs_c = if c[k + 1] < 0.0 { -c[k + 1] } else { c[k + 1] };
        if abs_d < abs_c {
            exch(d[k], c[k + 1]);
            exch(e[k], d[k + 1]);
            exch(c[k], e[k + 1]);
            exch(b[k], b[k + 1]);
        }
        if d[k] == 0.0 {
            singular();
        }
        let temp = c[k + 1] / d[k];
        d[k + 1] -= temp * e[k];
        e[k + 1] -= temp * c[k];
        b[k + 1] -= temp * b[k];
    }

    // Substitution.
    if d[n - 1] == 0.0 {
        singular();
    }
    b[n - 1] /= d[n - 1];
    if n == 1 {
        return;
    }
    b[n - 2] = (b[n - 2] - b[n - 1] * e[n - 2]) / d[n - 2];
    for k in (0..n - 2).rev() {
        b[k] = (b[k] - b[k + 1] * e[k] - b[k + 2] * c[k]) / d[k];
    }
}

fn singular() {
    use std::io::Write;
    let _ = std::io::stderr().write_all(b"coefficient matrix is singular!\n");
}
