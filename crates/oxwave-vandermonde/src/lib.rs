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

//! Port of cwp/lib/vandermonde.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.
//! Adapted from Golub and Van Loan, *Matrix Computations* (1983), Algorithm 5.6-2.

/// Solve a Vandermonde system in double precision.
///
/// C prototype: `void vanded (int n, double v[], double b[], double x[]);`
/// The dimension is `v.len()`; `b` and `x` must contain at least that many elements.
/// Repeated nodes cause division by zero, as in C.
pub fn vanded(v: &[f64], b: &[f64], x: &mut [f64]) {
    let n = v.len();
    assert!(b.len() >= n && x.len() >= n);
    for i in 0..n {
        x[i] = b[i];
    }
    for i in 0..n.saturating_sub(1) {
        for j in (i + 1..n).rev() {
            x[j] -= v[i] * x[j - 1];
        }
    }
    for i in (1..n).rev() {
        for j in i..n {
            x[j] /= v[j] - v[j - i];
        }
        for j in i..n {
            x[j - 1] -= x[j];
        }
    }
}

/// Solve a Vandermonde system in single precision.
///
/// C prototype: `void vandef (int n, float v[], float b[], float x[]);`
/// The dimension is `v.len()`; `b` and `x` must contain at least that many elements.
/// Repeated nodes cause division by zero, as in C.
pub fn vandef(v: &[f32], b: &[f32], x: &mut [f32]) {
    let n = v.len();
    assert!(b.len() >= n && x.len() >= n);
    for i in 0..n {
        x[i] = b[i];
    }
    for i in 0..n.saturating_sub(1) {
        for j in (i + 1..n).rev() {
            x[j] -= v[i] * x[j - 1];
        }
    }
    for i in (1..n).rev() {
        for j in i..n {
            x[j] /= v[j] - v[j - i];
        }
        for j in i..n {
            x[j - 1] -= x[j];
        }
    }
}
