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


//! SUIFFT - fft complex frequency traces to real time traces.
//! Port of su/main/transforms/suifft.c (Revision 1.15).
//! Credits (from the C source): CWP: Shuki, Chris, Jack.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
use oxwave_pfa_fft::pfacr;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUIFFT - fft complex frequency traces to real time traces\t",
    " \t\t\t\t\t\t\t\t",
    " suiftt <stdin >sdout sign=-1\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    " \tsign=-1\t\tsign in exponent of inverse fft\t\t",
    " \t\t\t\t\t\t\t\t",
    " Output traces are normalized by 1/N where N is the fft size.\t",
    " \t\t\t\t\t\t\t\t",
    " Note: sufft | suifft is not quite a no-op since the trace\t",
    " \tlength will usually be longer due to fft padding.\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get info from first trace */
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    if tr.trid != 111 { err!("input not complex freq data, trid=%d", tr.trid); }
    let nfft = tr.ns as i32 - 2;
    let nf = nfft / 2 + 1;
    let onfft = (1.0f64 / nfft as f64) as f32;
    let newd1 = if tr.d1 != 0.0 {
        onfft / tr.d1
    } else if tr.dt != 0 {
        (tr.dt as f64 / 1000000.0) as f32
    } else {
        0.0f32
    };

    /* Set sign in exponent of transform */
    let mut sign = 0;
    if getparint("sign", &mut sign) == 0 { sign = -1; }
    if sign != 1 && sign != -1 { err!("sign = %d must be 1 or -1", sign); }
    checkpars();

    let mut ct = vec![Complex::default(); nf as usize];
    let mut rt = vec![0.0f32; nfft as usize];
    loop {
        for i in 0..nf as usize {
            ct[i].r = tr.data[2 * i];
            ct[i].i = tr.data[2 * i + 1];
        }
        pfacr(sign, &ct, &mut rt);
        for i in 0..nfft as usize { tr.data[i] = rt[i] * onfft; }
        tr.trid = 1;
        tr.ns = nfft as u16;
        tr.f1 = 0.0;
        tr.d1 = newd1;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
