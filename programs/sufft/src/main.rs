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

//! SUFFT - fft real time traces to complex frequency traces.
//! Port of su/main/transforms/sufft.c (Revision 1.31).
//! Credits: CWP: Shuki Ronen, Chris Liner, Jack K. Cohen;
//! CENPET: Werner M. Heigl - added well log support.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
use oxwave_pfa_fft::{npfaro, pfarc};
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUFFT - fft real time traces to complex frequency traces\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suftt <stdin >sdout sign=1 \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " sign=1\t\t\tsign in exponent of fft\t\t\t",
    " dt=from header\t\tsampling interval\t\t\t",
    " verbose=1\t\t=0 to stop advisory messages\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: To facilitate further processing, the sampling interval\t",
    " in frequency and first frequency (0) are set in the\t\t\t",
    " output header.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sufft | suifft is not quite a no-op since the trace\t\t\t",
    " length will usually be longer due to fft padding.\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveats: \t\t\t\t\t\t\t\t",
    " No check is made that the data IS real time traces!\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Output is type complex. To view amplitude, phase or real, imaginary\t",
    " parts, use    suamp \t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Examples: \t\t\t\t\t\t\t\t",
    " sufft < stdin | suamp mode=amp | .... \t\t\t\t",
    " sufft < stdin | suamp mode=phase | .... \t\t\t\t",
    " sufft < stdin | suamp mode=uphase | .... \t\t\t\t",
    " sufft < stdin | suamp mode=real | .... \t\t\t\t",
    " sufft < stdin | suamp mode=imag | .... \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 1;
    getparint("verbose", &mut verbose);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let mut dt: f32;
    if seismic {
        if verbose != 0 { warn!("input is seismic data, trid=%d", tr.trid); }
        dt = (tr.dt as f64 / 1000000.0) as f32;
    } else {
        if verbose != 0 { warn!("input is not seismic data, trid=%d", tr.trid); }
        dt = tr.d1;
    }
    if dt == 0.0 {
        dt = 0.004;
        if verbose != 0 { warn!("dt or d1 not set, assumed to be .004"); }
    }
    let nfft = npfaro(nt, 2 * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    let nf = nfft / 2 + 1;
    let d1 = (1.0 / (nfft as f32 * dt) as f64) as f32;
    let mut sign = 1;
    getparint("sign", &mut sign);
    if sign != 1 && sign != -1 { err!("sign = %d must be 1 or -1", sign); }
    checkpars();
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nf as usize];
    if dt != 0.0 && verbose != 0 { warn!("d1=%f", 1.0 / (nfft as f32 * dt) as f64); }
    loop {
        rt[..nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        rt[nt as usize..].fill(0.0);
        pfarc(sign, &rt, &mut ct);
        for i in 0..nf as usize {
            tr.data[2*i] = ct[i].r;
            tr.data[2*i+1] = ct[i].i;
        }
        tr.ns = (2 * nf) as u16;
        tr.trid = 111;
        tr.d1 = d1;
        tr.f1 = 0.0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
