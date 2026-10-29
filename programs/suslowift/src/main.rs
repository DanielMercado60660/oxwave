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

//! SUSLOWIFT - inverse Fourier transform by direct DFT.
//! Port of su/main/transforms/suslowift.c.
//! Credits: CWP John Stockwell, based on suifft by Shuki Ronen,
//! Chris Liner, and Jack K. Cohen. The DFT comments cite Numerical Recipes in C.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUSLOWIFT - Fourier Transforms by (SLOW) DFT algorithm (Not an FFT)",
    "             complex frequency to real time domain traces \t",
    " \t\t\t\t\t\t\t\t",
    " suslowift <stdin >sdout sign=-1 \t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tsign=-1\t\t\tsign in exponent of fft\t\t",
    " \tdt=from header\t\tsampling interval\t\t",
    " \t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns, dt\t\t\t\t",
    " Trace header fields modified: ns, dt, trid\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Notes: To facilitate further processing, the sampling interval",
    "       in frequency and first frequency (0) are set in the\t",
    "\toutput header.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Warning: This program is *not* fft based. Use only for demo \t",
    " \t   purposes, *not* for large data processing.\t\t",
    " \t\t\t\t\t\t\t\t",
    " \tNo check is made that the data are real time traces!\t",
    " \t\t\t\t\t\t\t\t",
];

fn dftcr(sign: i32, nsamp: usize, re: &mut [f32], ct: &[Complex]) {
    let mut tmpre = vec![0.0f32; nsamp];
    let mut tmpim = vec![0.0f32; nsamp];
    for ni in 0..nsamp {
        tmpim[ni] = 0.0;
        tmpre[ni] = 0.0;
        re[ni] = 0.0;
    }
    let twopibynsamp = (sign as f64 * 2.0 * std::f64::consts::PI / nsamp as f64) as f32;
    for ni in 0..nsamp {
        let mut retemp = 0.0f32;
        let mut imtemp = 0.0f32;
        for ki in 0..nsamp {
            let angle = twopibynsamp * ki as f32 * ni as f32;
            let tmpcos = (angle as f64).cos() as f32;
            let tmpsin = (angle as f64).sin() as f32;
            retemp += ct[ki].r * tmpcos - ct[ki].i * tmpsin;
            imtemp += ct[ki].i * tmpcos + ct[ki].r * tmpsin;
        }
        tmpre[ni] = retemp;
        tmpim[ni] = imtemp;
    }
    for ni in 0..nsamp { re[ni] = tmpre[ni]; }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let nfft = nt as f32;
    let twonfft = 2.0f32 / nfft;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if dt == 0.0 {
        dt = 0.004;
        warn!("dt not set, assumed to be .004");
    }
    let mut sign = -1i32;
    if getparint("sign", &mut sign) == 0 { sign = -1; }
    if sign != 1 && sign != -1 { err!("sign = %d must be 1 or -1", sign); }
    // The C implementation reads 2 * ns floats from the fixed trace buffer.
    // Larger inputs access beyond that buffer and have undefined behavior.
    if nt > tr.data.len() / 2 { err!("trace sample count %d exceeds complex input capacity", nt as i32); }
    let mut rt = vec![0.0f32; nt];
    let mut ct = vec![Complex::default(); nt];
    if dt != 0.0 { warn!("d1=%f f1=0.0", 1.0f64 / (nt as f32 * dt) as f64); }
    loop {
        for i in 0..nt {
            ct[i].r = tr.data[2 * i];
            ct[i].i = tr.data[2 * i + 1];
        }
        dftcr(sign, nt, &mut rt, &ct);
        for i in 0..nt { tr.data[i] = rt[i] * twonfft; }
        tr.trid = 1;
        tr.ns = nt as u16;
        tr.d1 = 0.0;
        tr.f1 = 0.0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
