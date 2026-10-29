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


//! SUSPECFX - Fourier spectrum of traces.
//! Port of su/main/transforms/suspecfx.c (Revision 1.19).
//! Credits: CWP: Dave (algorithm), Jack (reformatting for SU).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::rcabs;
use oxwave_pfa_fft::{npfaro, pfarc, Complex};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUSPECFX - Fourier SPECtrum (T -> F) of traces \t\t",
    " \t\t\t\t\t\t\t\t",
    " suspecfx <infile >outfile \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Note: To facilitate further processing, the sampling interval\t",
    "       in frequency and first frequency (0) are set in the\t",
    "\toutput header.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if dt == 0.0 {
        dt = 0.004;
        warn!("dt not set, assumed to be .004");
    }
    checkpars();
    let nfft = npfaro(nt, 2 * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    let nfby2p1 = nfft / 2 + 1;
    let d1 = (1.0 / ((nfft as f32 * dt) as f64)) as f32;
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nfby2p1 as usize];
    loop {
        rt[..nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        rt[nt as usize..].fill(0.0);
        pfarc(1, &rt, &mut ct);
        tr.data[0] = (rcabs(ct[0]) as f64 / 2.0) as f32;
        for i in 1..nfby2p1 as usize { tr.data[i] = rcabs(ct[i]); }
        tr.ns = nfby2p1 as u16;
        tr.dt = 0;
        tr.trid = 118;
        tr.d1 = d1;
        tr.f1 = 0.0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
