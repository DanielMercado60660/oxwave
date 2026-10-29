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

//! Port of cwp/lib/linregress.c.
//! Author: Chris Liner, UTulsa, 11/16/03.

/// Compute slope, intercept, correlation coefficient, and standard error.
///
/// C prototype: `void linear_regression(float *y, float *x, int n, float coeff[4]);`
///
/// `coeff` receives the slope, intercept, absolute correlation coefficient,
/// and standard error, in that order. A zero correlation denominator yields
/// the C routine's sentinel value of 999.
///
/// The slices must have equal, nonzero lengths representable by a C `int`.
/// As in C, degenerate data can yield NaNs or infinities.
pub fn linear_regression(y: &[f32], x: &[f32], coeff: &mut [f32; 4]) {
    assert!(!x.is_empty() && x.len() == y.len() && x.len() <= i32::MAX as usize);
    let n = x.len() as i32;

    // Initialize sums from element zero, as in the C routine.
    let mut sx = x[0];
    let mut sx2 = x[0] * x[0];
    let mut sy = y[0];
    let mut sy2 = y[0] * y[0];
    let mut sxy = x[0] * y[0];
    for i in 1..x.len() {
        sx += x[i];
        sx2 += x[i] * x[i];
        sy += y[i];
        sy2 += y[i] * y[i];
        sxy += x[i] * y[i];
    }

    // Slope.
    let mut num = n as f32 * sxy - sx * sy;
    let mut den = n as f32 * sx2 - sx * sx;
    let a = num / den;
    coeff[0] = a;

    // Intercept.
    let b = (sy - a * sx) / n as f32;
    coeff[1] = b;

    // Correlation coefficient; preserve the C ABS macro's signed-zero behavior.
    num = n as f32 * sxy - sx * sy;
    if num < 0.0 {
        num = -num;
    }
    den = (n as f32 * sx2 - sx * sx) * (n as f32 * sy2 - sy * sy);
    // C sqrt takes a double, then the assignment narrows back to float.
    den = (den as f64).sqrt() as f32;
    let cc = if den != 0.0 { num / den } else { 999.0 };
    coeff[2] = cc;

    // Standard error. Division takes place in float before calling sqrt.
    let mut tmp = 0.0f32;
    for i in 0..x.len() {
        tmp += (y[i] - (a * x[i] + b)) * (y[i] - (a * x[i] + b));
    }
    let se = ((tmp / (n - 2) as f32) as f64).sqrt() as f32;
    coeff[3] = se;
}
