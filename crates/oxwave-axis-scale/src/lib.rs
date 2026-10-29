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

//! Port of cwp/lib/scaxis.c.
//! Author: Dave Hale, Colorado School of Mines, 01/13/89.

fn checked_int(value: f32) -> i32 {
    let truncated = (value as f64).trunc();
    if !truncated.is_finite() || truncated < i32::MIN as f64 || truncated > i32::MAX as f64 {
        oxwave_core::err!("scaxis: invalid or unrepresentable axis scale");
    }
    truncated as i32
}

/// Compute a readable scale, updating the tick count, interval and first tick.
///
/// C prototype: `void scaxis(float x1, float x2, int *nxnum, float *dxnum, float *fxnum);`
///
/// Exits with an error if the scale is unrepresentable or tick adjustment
/// cannot advance in `f32` arithmetic.
pub fn scaxis(x1: f32, x2: f32, nxnum: &mut i32, dxnum: &mut f32, fxnum: &mut f32) {
    let rdint = [1.0f32, 2.0, 5.0, 10.0];
    if x1 == x2 {
        *nxnum = 1;
        *dxnum = 1.0;
        *fxnum = x1;
        return;
    }
    let xmin = if x1 < x2 { x1 } else { x2 };
    let xmax = if x1 > x2 { x1 } else { x2 };
    let mut n = *nxnum;
    n = if 2 > n { 2 } else { n };
    let a = (xmax - xmin) / ((n - 1) as f32);
    if !a.is_finite() || a <= 0.0 {
        oxwave_core::err!("scaxis: invalid or unrepresentable axis scale");
    }
    let mut iloga = (a as f64).log10() as i32;
    if (a as f64) < 1.0 { iloga = iloga.wrapping_sub(1); }
    let b = ((a as f64) / 10.0f64.powf(iloga as f64)) as f32;
    let mut i = 0;
    while i < 3 && (b as f64) >= ((rdint[i] * rdint[i + 1]) as f64).sqrt() {
        i += 1;
    }
    let d = ((rdint[i] as f64) * 10.0f64.powf((iloga as f32) as f64)) as f32;
    if !d.is_finite() || d <= 0.0 {
        oxwave_core::err!("scaxis: invalid or unrepresentable axis scale");
    }
    let mut f = (checked_int(xmin / d) as f32) * d - d;
    let eps = (0.0001f64 * ((xmax - xmin) as f64)) as f32;
    while f < xmin - eps {
        let next = f + d;
        if next == f {
            oxwave_core::err!("scaxis: tick increment does not advance");
        }
        f = next;
    }
    n = checked_int((xmax + eps - f) / d).checked_add(1).unwrap_or_else(|| {
        oxwave_core::err!("scaxis: invalid or unrepresentable axis scale")
    });
    *nxnum = n;
    *dxnum = d;
    *fxnum = f;
}
