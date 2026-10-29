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

//! Ports of cwp/lib/intlin.c and intlinc.c (Dave Hale, 1989; complex
//! adaptation Hans Ecke, 2002), and su/lib/linterpd.c (no credits listed).
use oxwave_cwp::complex::{cadd, crmul, csub, Complex};
use oxwave_index_search::xindex;

/// Upper-bound search. C: `int bhighd(double *all, int last, double iguy)`.
pub fn bhighd(all: &[f64], iguy: f64) -> i32 {
    assert!(all.len() <= i32::MAX as usize);
    let (mut low, mut high) = (0, all.len());
    while low < high {
        let mid = low + (high - low) / 2;
        if iguy >= all[mid] {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    low as i32
}
/// Upper-bound search. C: `int bhighi(int *all, int last, int iguy)`.
pub fn bhighi(all: &[i32], iguy: i32) -> i32 {
    assert!(all.len() <= i32::MAX as usize);
    let (mut low, mut high) = (0, all.len());
    while low < high {
        let mid = low + (high - low) / 2;
        if iguy >= all[mid] {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    low as i32
}
/// Interpolate on a nonempty monotonic table.
/// C: `void intlin(int nin, float xin[], float yin[], float yinl, float yinr, int nout, float xout[], float yout[])`.
pub fn intlin(xin: &[f32], yin: &[f32], yinl: f32, yinr: f32, xout: &[f32], yout: &mut [f32]) {
    assert!(!xin.is_empty() && xin.len() == yin.len() && xout.len() == yout.len());
    let last = xin.len() - 1;
    let increasing = xin[0] <= xin[last];
    let mut idx = 0;
    for j in 0..xout.len() {
        let x = xout[j];
        if if increasing { x < xin[0] } else { x > xin[0] } {
            yout[j] = yinl;
        } else if if increasing {
            x > xin[last]
        } else {
            x < xin[last]
        } {
            yout[j] = yinr;
        } else if x == xin[last] || xin.len() == 1 {
            yout[j] = yin[last];
        } else {
            xindex(xin, x, &mut idx);
            let i = idx as usize;
            yout[j] = yin[i] + (x - xin[i]) * (yin[i + 1] - yin[i]) / (xin[i + 1] - xin[i]);
        }
    }
}
/// Interpolate complex values on a nonempty monotonic table.
/// C: `void intlinc(int nin, float xin[], complex yin[], complex yinl, complex yinr, int nout, float xout[], complex yout[])`.
pub fn intlinc(
    xin: &[f32],
    yin: &[Complex],
    yinl: Complex,
    yinr: Complex,
    xout: &[f32],
    yout: &mut [Complex],
) {
    assert!(!xin.is_empty() && xin.len() == yin.len() && xout.len() == yout.len());
    let last = xin.len() - 1;
    let increasing = xin[0] <= xin[last];
    let mut idx = 0;
    for j in 0..xout.len() {
        let x = xout[j];
        if if increasing { x < xin[0] } else { x > xin[0] } {
            yout[j] = yinl;
        } else if if increasing {
            x > xin[last]
        } else {
            x < xin[last]
        } {
            yout[j] = yinr;
        } else if x == xin[last] || xin.len() == 1 {
            yout[j] = yin[last];
        } else {
            xindex(xin, x, &mut idx);
            let i = idx as usize;
            yout[j] = cadd(
                yin[i],
                crmul(
                    csub(yin[i + 1], yin[i]),
                    (x - xin[i]) / (xin[i + 1] - xin[i]),
                ),
            );
        }
    }
}
/// Interpolate in double precision. Modes 0, 1, 2, 3 extrapolate neither, both, lower, or upper end.
/// C: `void linterpd(double offset, double *offs, double *tims, int nto, int mgtextr, double *time)`.
/// A NaN offset with at least two entries exits via `err!`; C would access outside the arrays.
pub fn linterpd(offset: f64, offs: &[f64], tims: &[f64], mgtextr: i32, time: &mut f64) {
    assert!(!offs.is_empty() && offs.len() == tims.len() && offs.len() <= i32::MAX as usize);
    let nto = offs.len();
    if nto < 2 {
        *time = tims[0];
        return;
    }
    // NaN bypasses both endpoint checks and makes bhighd return zero.
    // C then reads tims[-1] and offs[-1]; report this undefined input cleanly.
    if offset.is_nan() {
        oxwave_core::err!("linterpd: offset is NaN");
    }
    let n;
    if offset <= offs[0] {
        if mgtextr == 0 || mgtextr == 3 {
            *time = tims[0];
            return;
        }
        n = 1;
    } else if offset >= offs[nto - 1] {
        if mgtextr == 0 || mgtextr == 2 {
            *time = tims[nto - 1];
            return;
        }
        n = nto - 1;
    } else {
        n = bhighd(offs, offset) as usize;
    }
    let wo = (offs[n] - offset) / (offs[n] - offs[n - 1]);
    *time = wo * tims[n - 1] + (1.0 - wo) * tims[n];
}
