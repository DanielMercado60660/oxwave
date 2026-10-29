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

//! Phase unwrapping, ported from cwp/lib/unwrapphase.c.
//! Author: John Stockwell, CWP, 2010.
//! The private differentiator follows cwp/lib/differentiate.c;
//! Author: John Stockwell, CWP, 1994.

use oxwave_cwp::cwp::PI;
use oxwave_linear_regression::linear_regression;

/// Unwrap phase by integrating its frequency derivative (Oppenheim, 1975).
///
/// C prototype:
/// `void oppenheim_unwrap_phase(int n, int trend, int zeromean, float df, float *xr, float *xi, float *phase);`
///
/// All slices must have the same nonzero length, representable by a C `int`.
/// `df` is the frequency sampling interval. Nonzero `trend` removes the fitted
/// slope; nonzero `zeromean` also removes the intercept when `trend` is set.
/// Degenerate data retains C's NaN/infinity behavior. For one sample, the
/// result is zero without detrending and NaN with detrending.
pub fn oppenheim_unwrap_phase(
    trend: i32,
    zeromean: i32,
    df: f32,
    xr: &[f32],
    xi: &[f32],
    phase: &mut [f32],
) {
    let n = phase.len();
    assert!(n > 0 && n <= i32::MAX as usize && xr.len() == n && xi.len() == n);

    let mut phaseprime = vec![0.0f32; n];
    let mut xrprime = vec![0.0f32; n];
    let mut xiprime = vec![0.0f32; n];
    let mut freq = vec![0.0f32; n];
    phase.fill(0.0);

    phase[0] = 0.0;
    xrprime[0] = 0.0;
    xiprime[0] = 0.0;
    phaseprime[0] = 0.0;

    differentiate(df, xr, &mut xrprime);
    differentiate(df, xi, &mut xiprime);

    for i in 0..n {
        // These expressions are evaluated in float before assignment to double.
        let denominator = (xr[i] * xr[i] + xi[i] * xi[i]) as f64;
        let numerator = (xr[i] * xiprime[i] - xrprime[i] * xi[i]) as f64;

        // CLOSETO(denominator, 0.0) uses double arithmetic here. Its relative
        // tolerance is zero; keep the negated comparison for NaN behavior.
        if !(denominator.abs() <= 0.0) {
            phaseprime[i] = (numerator / denominator) as f32;
        } else {
            phaseprime[i] = 0.0;
        }

        // The product is float; division by 2.0 and addition are double.
        if i > 0 {
            phase[i] = (phase[i - 1] as f64
                + (df * (phaseprime[i] + phaseprime[i - 1])) as f64 / 2.0)
                as f32;
        }
        freq[i] = i as f32;
    }

    let mut coeff = [0.0f32; 4];
    linear_regression(phase, &freq, &mut coeff);

    if trend != 0 {
        for i in 0..n {
            if zeromean != 0 {
                phase[i] = phase[i] - (coeff[0] * i as f32 + coeff[1]);
            } else {
                phase[i] = phase[i] - coeff[0] * i as f32;
            }
        }
    }
}

/// Unwrap phase by correcting jumps in consecutive phase differences.
///
/// C prototype:
/// `void simple_unwrap_phase(int n, int trend, int zeromean, float w, float *phase);`
///
/// The slice must have nonzero length representable by a C `int`. Jumps at
/// least `PI/w` are replaced by the average of neighboring differences (or
/// the preceding difference at the last sample). Nonzero `trend` removes
/// the fitted slope, and nonzero `zeromean` then also removes the intercept.
/// As in C, `w == 0` writes `wrapping parameter is zero` to stderr and exits
/// with status 1. Detrending a single sample produces NaN, as in C.
pub fn simple_unwrap_phase(trend: i32, zeromean: i32, w: f32, phase: &mut [f32]) {
    if w == 0.0 {
        oxwave_core::fprintf!(oxwave_core::cio::stderr(), "wrapping parameter is zero");
        oxwave_core::exit(1);
    }
    let pibyw = (PI / w as f64) as f32;
    let n = phase.len();
    assert!(n > 0 && n <= i32::MAX as usize);

    let mut dphase = vec![0.0f32; n];
    let mut temp = vec![0.0f32; n];
    let mut freq = vec![0.0f32; n];
    temp[0] = phase[0];
    dphase[0] = 0.0;
    freq[0] = 0.0;

    for i in 1..n {
        dphase[i] = phase[i] - phase[i - 1];
        if i < n - 1 {
            dphase[i + 1] = phase[i + 1] - phase[i - 1];
        }

        let jump = dphase[i] - dphase[i - 1];
        // C's ABS macro preserves negative zero and leaves NaN unchanged.
        let abs_jump = if jump < 0.0 { -jump } else { jump };
        if abs_jump >= pibyw {
            if i < n - 1 {
                dphase[i] = ((dphase[i - 1] + dphase[i + 1]) as f64 / 2.0) as f32;
            } else {
                dphase[i] = dphase[i - 1];
            }
        }
        temp[i] = temp[i - 1] + dphase[i];
        freq[i] = i as f32;
    }

    if trend == 0 {
        for i in 0..n {
            phase[i] = temp[i];
        }
    } else {
        let mut coeff = [0.0f32; 4];
        linear_regression(&temp, &freq, &mut coeff);
        for i in 0..n {
            if zeromean != 0 {
                phase[i] = temp[i] - (coeff[0] * i as f32 + coeff[1]);
            } else {
                phase[i] = temp[i] - coeff[0] * i as f32;
            }
        }
    }
}

// Leading/centered/lagging differences from cwp/lib/differentiate.c.
fn differentiate(h: f32, f: &[f32], fprime: &mut [f32]) {
    let n = f.len();
    let h2 = 2.0f32 * h;
    if n == 1 {
        // C reads f[1] and f[-1] here. The resulting derivative is never
        // integrated when n == 1, so omit those invalid reads. phase[0] and
        // its detrended result are identical regardless of these values.
        return;
    }
    fprime[0] = (f[1] - f[0]) / h;
    for i in 1..n - 1 {
        fprime[i] = (f[i + 1] - f[i - 1]) / h2;
    }
    fprime[n - 1] = (f[n - 1] - f[n - 2]) / h;
}
