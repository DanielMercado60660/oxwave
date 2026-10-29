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
//! SUSLOWFT - Fourier transforms by a slow DFT.
//! Port of su/main/transforms/suslowft.c.
//! Credits: CWP: Shuki, Chris, Jack. DFT: John Stockwell, 26 Feb 1992;
//! based on the well known algorithm, see Numerical Recipes in C, p. 406.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
use oxwave_cwp::cwp::PI;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUSLOWFT - Fourier Transforms by a (SLOW) DFT algorithm (Not an FFT)",
    " \t\t\t\t\t\t\t\t",
    " suslowft <stdin >sdout sign=1 \t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tsign=1\t\t\tsign in exponent of fft\t\t",
    " \tdt=from header\t\tsampling interval\t\t",
    " \t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns, dt\t\t\t\t",
    " Trace header fields modified: ns, dt, trid\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Notes: To facilitate further processing, the sampling interval",
    "       in frequency and first frequency (0) are set in the\t",
    "\toutput header.\t\t\t\t\t\t",
    " Warning: This program is *not* fft based. Use only for demo \t",
    " \t   purposes, *not* for large data processing.\t\t",
    " \t\t\t\t\t\t\t\t",
    " \tNo check is made that the data are real time traces!\t",
    " suslowft | suslowift is not quite a no-op since the trace     ",
    " length will usually be longer due to fft padding.             ",
    "                                                                       ",
    " Caveats:                                                              ",
    " No check is made that the data IS real time traces!                   ",
    "                                                                       ",
    " Output is type complex. To view amplitude, phase or real, imaginary   ",
    " parts, use    suamp                                                   ",
    "                                                                       ",
    " Examples:                                                             ",
    " suslowft < stdin | suamp mode=amp | ....                                 ",
    " suslowft < stdin | suamp mode=phase | ....                               ",
    " suslowft < stdin | suamp mode=real | ....                                ",
    " suslowft < stdin | suamp mode=imag | ....                                ",
    "                                                                       ",
    " \t\t\t\t\t\t\t\t",
];

fn dftrc(sign: i32, nsamp: usize, re: &[f32], out: &mut [Complex]) {
    let mut tmpre = vec![0.0f32; nsamp];
    let mut tmpim = vec![0.0f32; nsamp];
    for ni in 0..nsamp { out[ni] = Complex::default(); }
    let twopibynsamp = (sign as f64 * 2.0 * PI / nsamp as f64) as f32;
    let mut freq = 0.0f32;
    for ni in 0..nsamp {
        let mut retemp = 0.0f32;
        let mut imtemp = 0.0f32;
        for ki in 0..nsamp {
            let phase = freq * ki as f32;
            let tmpcos = (phase as f64).cos() as f32;
            let tmpsin = (phase as f64).sin() as f32;
            retemp += re[ki] * tmpcos;
            imtemp += re[ki] * tmpsin;
        }
        tmpre[ni] = retemp;
        tmpim[ni] = imtemp;
        freq += twopibynsamp;
    }
    for ni in 0..nsamp {
        out[ni].r = tmpre[ni];
        out[ni].i = tmpim[ni];
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if dt == 0.0 {
        dt = 0.004;
        warn!("dt not set, assumed to be .004");
    }
    let d1 = (1.0f64 / (nt as f32 * dt) as f64) as f32;
    let mut sign = 1;
    if getparint("sign", &mut sign) == 0 { sign = 1; }
    if sign != 1 && sign != -1 { err!("sign = %d must be 1 or -1", sign); }
    let mut rt = vec![0.0f32; nt];
    let mut ct = vec![Complex::default(); nt];
    if dt != 0.0 { warn!("d1=%f f1=0.0", 1.0f64 / (nt as f32 * dt) as f64); }
    loop {
        rt.copy_from_slice(&tr.data[..nt]);
        dftrc(sign, nt, &rt, &mut ct);
        // puttr writes only nt floats; the rest is overwritten by the next input.
        for i in 0..nt {
            tr.data[i] = if i % 2 == 0 { ct[i / 2].r } else { ct[i / 2].i };
        }
        tr.trid = 111;
        tr.ns = nt as u16;
        tr.d1 = d1;
        tr.f1 = 0.0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
