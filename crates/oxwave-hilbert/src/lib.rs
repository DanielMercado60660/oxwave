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

//! Port of cwp/lib/hilbert.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

use std::sync::OnceLock;

// cwp.h defines PI with this precision.
const PI: f64 = 3.141592653589793;

/// Compute the Hilbert transform using the SU 61-sample Hamming-windowed filter.
///
/// C prototype: `void hilbert (int n, float x[], float y[]);`
///
/// The input and output slices must have the same length. Samples outside the
/// input are treated as zero, as in the original convolution.
pub fn hilbert(x: &[f32], y: &mut [f32]) {
    assert_eq!(x.len(), y.len());
    static H: OnceLock<[f32; 61]> = OnceLock::new();
    let h = H.get_or_init(|| {
        let mut h = [0.0f32; 61];
        for i in 1..=30 {
            let taper = (0.54 + 0.46 * (PI * i as f64 / 30.0).cos()) as f32;
            h[30 + i] = (taper as f64 * (-(i as f32 % 2.0) as f64 * 2.0
                / (PI * i as f64))) as f32;
            h[30 - i] = -h[30 + i];
        }
        h
    });
    // The C convolution swaps its inputs when the signal is shorter than the
    // filter. Accumulate in that same input order.
    let n = x.len();
    // In the rolling-on region the optimized C convolution evaluates paired
    // output samples with the filter indices descending.
    let rolling_pairs_end = if n > 31 { ((n.min(61) - 31) / 2) * 2 } else { 0 };
    for i in 0..n {
        let mut sum = 0.0f32;
        if n < 61 {
            let low = i.saturating_sub(30);
            let high = n.min(i.saturating_add(31));
            if i < rolling_pairs_end && n > 30 {
                for j in (low..high).rev() {
                    sum += x[j] * h[(i + 30) - j];
                }
            } else {
                for j in low..high {
                    sum += x[j] * h[(i + 30) - j];
                }
            }
        } else {
            let low = (i + 31).saturating_sub(n);
            let high = 61.min(i + 31);
            if i < rolling_pairs_end {
                for j in (low..high).rev() {
                    sum += h[j] * x[i + 30 - j];
                }
            } else {
                for j in low..high {
                    sum += h[j] * x[i + 30 - j];
                }
            }
        }
        y[i] = sum;
    }
}
