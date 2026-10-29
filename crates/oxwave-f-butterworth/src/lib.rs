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

//! Port of cwp/lib/butterworth.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

const PI: f64 = 3.141592653589793;

/// Design a Butterworth filter from pass and stop band constraints.
///
/// C: `void bfdesign(float fpass, float apass, float fstop, float astop, int *npoles, float *f3db);`
/// Frequencies must lie between 0 and 0.5 and differ; amplitudes must
/// satisfy `0 < astop < apass < 1`. Designs a low-pass filter when
/// `fpass < fstop`, otherwise a high-pass filter.
pub fn bfdesign(fpass: f32, apass: f32, fstop: f32, astop: f32, npoles: &mut i32, f3db: &mut f32) {
    // Warp frequencies according to the bilinear transform.
    let wpass = (2.0 * (PI * fpass as f64).tan()) as f32;
    let wstop = (2.0 * (PI * fstop as f64).tan()) as f32;
    let fnpoles;
    let w3db;
    // C multiplies the float amplitudes and divides the warped frequencies
    // before promoting those intermediate results to double.
    if fstop > fpass {
        fnpoles = (((1.0 / (apass * apass) as f64 - 1.0) / (1.0 / (astop * astop) as f64 - 1.0))
            .ln()
            / ((wpass / wstop) as f64).powf(2.0).ln()) as f32;
        w3db =
            (wpass as f64 / (1.0 / (apass * apass) as f64 - 1.0).powf(0.5 / fnpoles as f64)) as f32;
    } else {
        fnpoles = (((1.0 / (apass * apass) as f64 - 1.0) / (1.0 / (astop * astop) as f64 - 1.0))
            .ln()
            / ((wstop / wpass) as f64).powf(2.0).ln()) as f32;
        w3db =
            (wpass as f64 * (1.0 / (apass * apass) as f64 - 1.0).powf(0.5 / fnpoles as f64)) as f32;
    }
    // Determine the integer number of poles, then unwarp the -3 dB frequency.
    *npoles = 1i32.wrapping_add(fnpoles as i32);
    *f3db = ((0.5 * w3db as f64).atan() / PI) as f32;
}

/// Apply a high-pass Butterworth filter. `q` must have at least `p.len()` elements.
///
/// C: `void bfhighpass(int npoles, float f3db, int n, float p[], float q[]);`
/// C requires `npoles >= 0` and `0 <= f3db <= 0.5`; Nyquist is `0.5`.
/// Like C's `fabs`, negative `f3db`
/// values are treated as their absolute values. Each call starts with zero
/// filter history. Only the first `p.len()` output elements are written.
pub fn bfhighpass(npoles: i32, f3db: f32, p: &[f32], q: &mut [f32]) {
    assert!(q.len() >= p.len());
    let r = (2.0 * (PI * (f3db as f64).abs()).tan()) as f32;
    if npoles % 2 != 0 {
        let scale = (r as f64 + 2.0) as f32;
        let a = (2.0 / scale as f64) as f32;
        let b1 = ((r as f64 - 2.0) / scale as f64) as f32;
        let mut pj = 0.0f32;
        let mut qjm1 = 0.0f32;
        for j in 0..p.len() {
            let pjm1 = pj;
            pj = p[j];
            q[j] = a * (pj - pjm1) - b1 * qjm1;
            qjm1 = q[j];
        }
    } else {
        for j in 0..p.len() {
            q[j] = p[j];
        }
    }
    for jpair in 0..npoles / 2 {
        let theta = (PI * (2 * jpair + 1) as f64 / npoles.wrapping_mul(2) as f64) as f32;
        let scale = (4.0 + 4.0 * r as f64 * (theta as f64).sin() + (r * r) as f64) as f32;
        let a = (4.0 / scale as f64) as f32;
        let b1 = ((2.0 * r as f64 * r as f64 - 8.0) / scale as f64) as f32;
        let b2 =
            ((4.0 - 4.0 * r as f64 * (theta as f64).sin() + (r * r) as f64) / scale as f64) as f32;
        let mut pjm1 = 0.0f32;
        let mut pj = 0.0f32;
        let mut qjm2 = 0.0f32;
        let mut qjm1 = 0.0f32;
        for j in 0..p.len() {
            let pjm2 = pjm1;
            pjm1 = pj;
            pj = q[j];
            // The 2.0 literal promotes the numerator to double, while
            // each feedback product is evaluated in float before promotion.
            q[j] = (a as f64 * (pj as f64 - 2.0 * pjm1 as f64 + pjm2 as f64)
                - (b1 * qjm1) as f64
                - (b2 * qjm2) as f64) as f32;
            qjm2 = qjm1;
            qjm1 = q[j];
        }
    }
}

/// Apply a low-pass Butterworth filter. `q` must have at least `p.len()` elements.
///
/// C: `void bflowpass(int npoles, float f3db, int n, float p[], float q[]);`
/// C requires `npoles >= 0` and `0 <= f3db <= 0.5`; Nyquist is `0.5`.
/// Like C's `fabs`, negative `f3db`
/// values are treated as their absolute values. Each call starts with zero
/// filter history. Only the first `p.len()` output elements are written.
pub fn bflowpass(npoles: i32, f3db: f32, p: &[f32], q: &mut [f32]) {
    assert!(q.len() >= p.len());
    let r = (2.0 * (PI * (f3db as f64).abs()).tan()) as f32;
    if npoles % 2 != 0 {
        let scale = (r as f64 + 2.0) as f32;
        let a = r / scale;
        let b1 = ((r as f64 - 2.0) / scale as f64) as f32;
        let mut pj = 0.0f32;
        let mut qjm1 = 0.0f32;
        for j in 0..p.len() {
            let pjm1 = pj;
            pj = p[j];
            q[j] = a * (pj + pjm1) - b1 * qjm1;
            qjm1 = q[j];
        }
    } else {
        for j in 0..p.len() {
            q[j] = p[j];
        }
    }
    for jpair in 0..npoles / 2 {
        let theta = (PI * (2 * jpair + 1) as f64 / npoles.wrapping_mul(2) as f64) as f32;
        let scale = (4.0 + 4.0 * r as f64 * (theta as f64).sin() + (r * r) as f64) as f32;
        let a = r * r / scale;
        let b1 = ((2.0 * r as f64 * r as f64 - 8.0) / scale as f64) as f32;
        let b2 =
            ((4.0 - 4.0 * r as f64 * (theta as f64).sin() + (r * r) as f64) / scale as f64) as f32;
        let mut pjm1 = 0.0f32;
        let mut pj = 0.0f32;
        let mut qjm2 = 0.0f32;
        let mut qjm1 = 0.0f32;
        for j in 0..p.len() {
            let pjm2 = pjm1;
            pjm1 = pj;
            pj = q[j];
            // The 2.0 literal promotes the numerator to double, while
            // each feedback product is evaluated in float before promotion.
            q[j] = (a as f64 * (pj as f64 + 2.0 * pjm1 as f64 + pjm2 as f64)
                - (b1 * qjm1) as f64
                - (b2 * qjm2) as f64) as f32;
            qjm2 = qjm1;
            qjm1 = q[j];
        }
    }
}
