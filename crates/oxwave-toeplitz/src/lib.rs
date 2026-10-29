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

//! Port of cwp/lib/stoep.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

/// Solves a symmetric Toeplitz system using the top row `r` and right-hand side `g`.
///
/// C prototype: `void stoepd (int n, double r[], double g[], double f[], double a[]);`
///
/// The dimension is `r.len()`, which must be nonzero. The other slices must
/// contain at least that many elements. If `r[0]` is zero, outputs are unchanged.
/// `a` receives the auxiliary solution used by the recurrence. Elements beyond
/// the dimension are unchanged. Singular intermediate systems can produce
/// infinities or NaNs, just as in C.
///
/// # Panics
/// Panics if `r` is empty or another slice is shorter than `r`.
pub fn stoepd(r: &[f64], g: &[f64], f: &mut [f64], a: &mut [f64]) {
    let n = r.len();
    assert!(n > 0 && g.len() >= n && f.len() >= n && a.len() >= n);
    if r[0] == 0.0 {
        return;
    }
    a[0] = 1.0;
    let mut v = r[0];
    f[0] = g[0] / r[0];
    for j in 1..n {
        // Solve Ra=v as in Claerbout, FGDP, p. 57.
        a[j] = 0.0;
        f[j] = 0.0;
        let mut e = 0.0;
        for i in 0..j {
            e += a[i] * r[j - i];
        }
        let mut c = e / v;
        v -= c * e;
        for i in 0..=j / 2 {
            let bot = a[j - i] - c * a[i];
            a[i] -= c * a[j - i];
            a[j - i] = bot;
        }
        // Use a and v above to get f[i], i = 0, 1, ..., j.
        let mut w = 0.0;
        for i in 0..j {
            w += f[i] * r[j - i];
        }
        c = (w - g[j]) / v;
        for i in 0..=j {
            f[i] -= c * a[j - i];
        }
    }
}

/// Solves a symmetric Toeplitz system using the top row `r` and right-hand side `g`.
///
/// C prototype: `void stoepf (int n, float r[], float g[], float f[], float a[]);`
///
/// The dimension is `r.len()`, which must be nonzero. The other slices must
/// contain at least that many elements. If `r[0]` is zero, outputs are unchanged.
/// `a` receives the auxiliary solution used by the recurrence. Elements beyond
/// the dimension are unchanged. Singular intermediate systems can produce
/// infinities or NaNs, just as in C.
///
/// # Panics
/// Panics if `r` is empty or another slice is shorter than `r`.
pub fn stoepf(r: &[f32], g: &[f32], f: &mut [f32], a: &mut [f32]) {
    let n = r.len();
    assert!(n > 0 && g.len() >= n && f.len() >= n && a.len() >= n);
    if r[0] == 0.0 {
        return;
    }
    a[0] = 1.0;
    let mut v = r[0];
    f[0] = g[0] / r[0];
    for j in 1..n {
        // Solve Ra=v as in Claerbout, FGDP, p. 57.
        a[j] = 0.0;
        f[j] = 0.0;
        let mut e = 0.0;
        for i in 0..j {
            e += a[i] * r[j - i];
        }
        let mut c = e / v;
        v -= c * e;
        for i in 0..=j / 2 {
            let bot = a[j - i] - c * a[i];
            a[i] -= c * a[j - i];
            a[j - i] = bot;
        }
        // Use a and v above to get f[i], i = 0, 1, ..., j.
        let mut w = 0.0;
        for i in 0..j {
            w += f[i] * r[j - i];
        }
        c = (w - g[j]) / v;
        for i in 0..=j {
            f[i] -= c * a[j - i];
        }
    }
}
