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
//
// Portions of this file are a Rust translation of code that is
// Copyright (c) 2007 by the Society of Exploration Geophysicists.
// For more information, go to http://software.seg.org/2007/0004 .
// You must read and accept usage terms at:
// http://software.seg.org/disclaimer.txt before use.
//
// Original code: Werner M. Heigl, Apache Corporation, E&P Technology, 2006-2007,
// accompanying: Heigl, W. M., 2007, Computing Gaussian derivative waveforms of
// any order: Geophysics, 72(4), H39-H42, doi:10.1190/1.2716624.
// The SEG-derived code is marked below between `SEG-BEGIN` and `SEG-END` comments.
//
// Modification notice (SEG terms of use, item 4):
//   Nature:      translated from C to Rust; same algorithm, statement order and
//                floating-point precision as the C; no functional change intended.
//   Extent:      the SEG-BEGIN/SEG-END blocks in this file.
//   Purpose:     OxWave, a Rust port of CWP/SU (Seismic Un*x).
//   When:        October 2026.
//   Who:         the OxWave project, maintained by Daniel Mercado;
//                translation produced with AI coding agents and checked against the
//                original C by differential tests.
//   Affiliation: OxWave project (independent open-source project).
// Keep this notice in all redistributions.

//! Port of cwp/lib/waveforms.c. Authors: Tong Fei, Ken Larner,
//! Nils Maercklin, and Werner M. Heigl; Gaussian derivative code credited to
//! Dominique Rousset.
const PI: f64 = 3.141592653589793;

/// C prototype: `void ricker1_wavelet(int nt, float dt, float fpeak, float *wavelet)`.
pub fn ricker1_wavelet(dt: f32, fpeak: f32, wavelet: &mut [f32]) {
    let t0 = (1.0 / fpeak as f64) as f32;
    for (it, w) in wavelet.iter_mut().enumerate() {
        let t1 = it as f32 * dt;
        let d = (t1 - t0) as f64;
        *w = ((-PI * PI * fpeak as f64 * fpeak as f64 * d * d).exp()
            * (1.0 - 2.0 * PI * PI * fpeak as f64 * fpeak as f64 * d * d)) as f32;
    }
}
/// C prototype: `void ricker2_wavelet(int hlw, float dt, float period, float ampl, float distort, float *wavelet)`.
/// The output contains at least `2 * hlw - 1` samples and `hlw` is positive.
pub fn ricker2_wavelet(
    hlw: usize,
    dt: f32,
    period: f32,
    ampl: f32,
    distort: f32,
    wavelet: &mut [f32],
) {
    assert!(hlw > 0 && wavelet.len() >= 2 * hlw - 1);
    let mut t = 0.0f32;
    for j in 0..hlw {
        let tnorm = t / period;
        let xsq = ((2.5 * tnorm as f64) * (2.5 * tnorm as f64)) as f32;
        let wsym = (ampl as f64 * (1.0 - 2.0 * xsq as f64) * (-(xsq as f64)).exp()) as f32;
        wavelet[hlw + j - 1] = (wsym as f64 * (1.0 - 2.0 * distort as f64 * tnorm as f64)) as f32;
        wavelet[hlw - j - 1] = (wsym as f64 * (1.0 + 2.0 * distort as f64 * tnorm as f64)) as f32;
        t += dt;
    }
}
/// C prototype: `void akb_wavelet(int nt, float dt, float fpeak, float *wavelet)`.
pub fn akb_wavelet(dt: f32, fpeak: f32, wavelet: &mut [f32]) {
    let t0 = (1.8 / fpeak as f64) as f32;
    for (it, w) in wavelet.iter_mut().enumerate() {
        let t1 = it as f32 * dt;
        let d = (t1 - t0) as f64;
        *w = (-60.0 * d * (-2.0 * fpeak as f64 * fpeak as f64 * d * d).exp()) as f32;
    }
}
/// C prototype: `void spike_wavelet(int nt, int tindex, float *wavelet)`.
/// The index must be within the slice.
pub fn spike_wavelet(tindex: usize, wavelet: &mut [f32]) {
    for w in wavelet.iter_mut() {
        *w = 0.0;
    }
    wavelet[tindex] = 1.0;
}
/// C prototype: `void unit_wavelet(int nt, float *wavelet)`.
pub fn unit_wavelet(wavelet: &mut [f32]) {
    for w in wavelet {
        *w = 1.0;
    }
}
/// C prototype: `void zero_wavelet(int nt, float *wavelet)`.
pub fn zero_wavelet(wavelet: &mut [f32]) {
    for w in wavelet {
        *w = 0.0;
    }
}
/// C prototype: `void berlage_wavelet(int nt, float dt, float fpeak, float ampl, float tn, float decay, float ipa, float *wavelet)`.
pub fn berlage_wavelet(
    dt: f32,
    fpeak: f32,
    ampl: f32,
    tn: f32,
    decay: f32,
    ipa: f32,
    wavelet: &mut [f32],
) {
    for (it, w) in wavelet.iter_mut().enumerate() {
        let t = dt * it as f32;
        *w = (ampl as f64
            * (t as f64).powf(tn as f64)
            * ((-decay * t) as f64).exp()
            * (2.0 * PI * fpeak as f64 * t as f64 + ipa as f64).cos()) as f32;
    }
}
/// C prototype: `void gaussian_wavelet(int nt, float dt, float fpeak, float *wavelet)`.
pub fn gaussian_wavelet(dt: f32, fpeak: f32, wavelet: &mut [f32]) {
    let t0 = (1.0 / fpeak as f64) as f32;
    let s = (1.0 / (2.0f64.sqrt() * PI * fpeak as f64)) as f32;
    for (it, w) in wavelet.iter_mut().enumerate() {
        let d = dt * it as f32 - t0;
        *w = (1.0 / (s as f64 * (2.0 * PI).sqrt())
            * ((-(d * d) as f64) / (2.0 * s as f64 * s as f64)).exp()) as f32;
    }
}
/// C prototype: `void gaussderiv_wavelet(int nt, float dt, float fpeak, float *wavelet)`.
pub fn gaussderiv_wavelet(dt: f32, fpeak: f32, wavelet: &mut [f32]) {
    let t0 = (1.0 / fpeak as f64) as f32;
    let s = (1.0 / (2.0f64.sqrt() * PI * fpeak as f64)) as f32;
    for (it, w) in wavelet.iter_mut().enumerate() {
        let d = dt * it as f32 - t0;
        *w = (-1.0 / ((s * s * s) as f64 * (2.0 * PI).sqrt())
            * d as f64
            * ((-(d * d) as f64) / (2.0 * s as f64 * s as f64)).exp()) as f32;
    }
}
// SEG-BEGIN waveforms.c:deriv_n_gauss
/// Compute a normalized nth Gaussian derivative.
/// C prototype: `void deriv_n_gauss(double dt, int nt, double t0, float fpeak, int n, double *w, int sign, int verbose)`.
/// The output slice must be nonempty.
///
/// Notes:
/// Copyright (c) 2007 by the Society of Exploration Geophysicists.
/// For more information, go to http://software.seg.org/2007/0004 .
/// You must read and accept usage terms at:
/// http://software.seg.org/disclaimer.txt before use.
/// Author: Werner M. Heigl, Apache Corporation, E&P Technology, November 2006.
pub fn deriv_n_gauss(dt: f64, t0: f64, fpeak: f32, n: i32, w: &mut [f64], sign: i32, verbose: i32) {
    let nt = w.len();
    assert!(nt > 0);
    let mut t = vec![0.0; nt];
    let mut h = vec![0.0; nt];
    let mut h0 = vec![0.0; nt];
    let mut h1 = vec![0.0; nt];
    if verbose != 0 {
        eprint!("memory allocated and initialized/n");
    }
    for i in 0..nt {
        t[i] = i as f64 * dt - t0;
    }
    if verbose != 0 {
        eprint!("t[] initialized/n");
    }
    let sigma = n as f64 / (4.0 * PI * PI * fpeak as f64 * fpeak as f64);
    if verbose != 0 {
        eprint!("sigma={:.6}", sigma);
    }
    for i in 0..nt {
        w[i] = (-t[i] * t[i] / (2.0 * sigma)).exp();
    }
    if verbose != 0 {
        eprint!("Gaussian computed/n");
    }
    for i in 0..nt {
        h0[i] = 1.0;
        h1[i] = t[i] / sigma;
    }
    if n == 1 {
        h.copy_from_slice(&h1);
    }
    if n > 1 {
        oxwave_orthogonal_polynomials::hermite_n_polynomial(&mut h, &mut h0, &mut h1, &t, n, sigma);
    }
    if verbose != 0 {
        eprint!("Hermite polynomial H_{} computed/n", n);
    }
    for i in 0..nt {
        w[i] = h[i] * w[i];
    }
    if verbose != 0 {
        eprint!("waveform computed/n");
    }
    let mut c = w[0].abs();
    for &v in &w[1..] {
        if v.abs() > c {
            c = v.abs();
        }
    }
    if n & 1 != 0 {
        c = -c;
    }
    if verbose != 0 {
        eprint!("C={:.6}/n", c);
    }
    for v in w.iter_mut() {
        *v = sign as f64 * *v / c;
    }
    if verbose != 0 {
        eprint!("waveform normalized/nw[o]={:.12}memory freed/n", w[0]);
    }
}
// SEG-END
