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

//! Port of cwp/lib/xcor.c and its convolution calculation.
//! Credits: Dave Hale, Colorado School of Mines, 11/23/91.

/// Cross-correlate the input slices, writing samples beginning at `ifz`.
/// Slice lengths supply the C lengths `lx`, `ly`, and `lz`.
///
/// C prototype: `void xcor(int lx, int ifx, float *x, int ly, int ify,
/// float *y, int lz, int ifz, float *z);`
pub fn xcor(ifx: i32, x: &[f32], ify: i32, y: &[f32], ifz: i32, z: &mut [f32]) {
    let mut xr = vec![0.0f32; x.len()];
    for i in 0..x.len() {
        xr[i] = x[x.len() - 1 - i];
    }
    let fx = 1i64 - ifx as i64 - x.len() as i64;
    convolve(&xr, fx, y, ify as i64, z, ifz as i64);
}

// C swaps operands to sum the shorter one. Its paired rolling-on path
// accumulates backwards, and its unrolled short middle starts with a product.
fn convolve<'a>(mut x: &'a [f32], mut ifx: i64, mut y: &'a [f32], mut ify: i64, z: &mut [f32], ifz: i64) {
    if x.len() > y.len() {
        std::mem::swap(&mut x, &mut y);
        std::mem::swap(&mut ifx, &mut ify);
    }
    let ilx = ifx + x.len() as i64 - 1;
    let ily = ify + y.len() as i64 - 1;
    let rolling_start = ifz.max(ifx + ify);
    let rolling_end = (ify + ilx - 1).min(ifz + z.len() as i64 - 1);
    let mut paired_start = rolling_start;
    if (rolling_start - ify - ifx) % 2 != 0 {
        paired_start += 1;
    }
    let paired_end = paired_start + ((rolling_end - paired_start + 1).max(0) / 2) * 2;
    let middle_start = ifz.max(ify + ilx);
    let middle_end = (ily + ifx).min(ifz + z.len() as i64 - 1);
    let middle_paired_end = middle_start + ((middle_end - middle_start + 1).max(0) / 2) * 2;
    for (offset, out) in z.iter_mut().enumerate() {
        let i = ifz + offset as i64;
        let jlow = (i - ily).max(ifx);
        let jhigh = (i - ify).min(ilx);
        let mut sum = 0.0f32;
        if x.len() > 30 && i >= paired_start && i < paired_end {
            for j in (jlow..=jhigh).rev() {
                sum += x[(j - ifx) as usize] * y[(i - j - ify) as usize];
            }
        } else if !x.is_empty() && x.len() <= 30 && i >= middle_start && i < middle_paired_end {
            sum = x[0] * y[(i - ifx - ify) as usize];
            for j in ifx + 1..=jhigh {
                sum += x[(j - ifx) as usize] * y[(i - j - ify) as usize];
            }
        } else {
            for j in jlow..=jhigh {
                sum += x[(j - ifx) as usize] * y[(i - j - ify) as usize];
            }
        }
        *out = sum;
    }
}
