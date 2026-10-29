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

//! Port of cwp/lib/intcub.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

use oxwave_index_search::xindex;

/// Evaluate a piecewise cubic or one of its first three derivatives.
///
/// C prototype: `void intcub(int ideriv, int nin, float xin[], float ydin[][4], int nout, float xout[], float yout[])`.
///
/// `xin` must be nonempty and monotonic, and `ydin` must contain one row
/// per input coordinate. Each row contains the function value and its first
/// three derivatives. Queries outside the input range use the endpoint's
/// cubic polynomial. Results are written to the first `xout.len()` elements
/// of `yout`. Other derivative orders produce zero.
///
/// # Panics
/// Panics if the required input or output slices are too short, or if the
/// input length does not fit in a C `int`.
pub fn intcub(ideriv: i32, xin: &[f32], ydin: &[[f32; 4]], xout: &[f32], yout: &mut [f32]) {
    assert!(yout.len() >= xout.len());
    if !(0..=3).contains(&ideriv) {
        for y in &mut yout[..xout.len()] {
            *y = 0.0;
        }
        return;
    }
    if xout.is_empty() {
        return;
    }
    assert!(!xin.is_empty() && xin.len() <= i32::MAX as usize && ydin.len() >= xin.len());
    // C caches this search hint in a static variable. A per-call hint gives
    // the same interval for finite monotonic inputs and permits concurrent use.
    let mut idx = 0;
    for iout in 0..xout.len() {
        xindex(xin, xout[iout], &mut idx);
        let i = idx as usize;
        let delx = xout[iout] - xin[i];
        let y = ydin[i];
        yout[iout] = match ideriv {
            // O2 and O6 are unsuffixed C literals: their enclosing polynomial
            // is evaluated in double after the float subtraction above.
            0 => {
                (y[0] as f64
                    + delx as f64
                        * (y[1] as f64
                            + delx as f64
                                * (y[2] as f64 * 0.5000000
                                    + delx as f64 * (y[3] as f64 * 0.1666667))))
                    as f32
            }
            1 => {
                (y[1] as f64
                    + delx as f64 * (y[2] as f64 + delx as f64 * (y[3] as f64 * 0.5000000)))
                    as f32
            }
            2 => y[2] + delx * y[3],
            _ => y[3],
        };
    }
}
