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

//! SUPHASE - phase manipulation by linear transformation.
//! Port of su/main/filters/suphase.c (Revision 1.4).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
use oxwave_pfa_fft::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUPHASE - PHASE manipulation by linear transformation\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  suphase  <stdin >sdout      \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " a=90\t\t\tconstant phase shift              \t\t",
    " b=180/PI              linear phase shift\t\t\t\t",
    " c=0.0\t\t\tphase = a +b*(old_phase)+c*f;\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " A program that allows the user to experiment with changes in the phase",
    " spectrum of a signal.\t\t\t\t\t\t\t",
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
        warn!("dt not set, assumed to be .002");
    }
    let _ = dt;
    let mut a = 0.0f32;
    let mut b = 0.0f32;
    let mut c = 0.0f32;
    if getparfloat("a", &mut a) == 0 { a = 0.0; }
    if getparfloat("b", &mut b) == 0 { b = (180.0 / std::f64::consts::PI) as f32; }
    if getparfloat("c", &mut c) == 0 { c = 0.0; }
    a = (a as f64 * (std::f64::consts::PI / 180.0)) as f32;
    b = (b as f64 * (std::f64::consts::PI / 180.0)) as f32;
    let nfft = npfaro(nt, 2 * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    let nf = nfft as usize / 2 + 1;
    let onfft = (1.0 / nfft as f64) as f32;
    checkpars();
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nf];
    let mut amp = vec![0.0f32; nf];
    let mut ph = vec![0.0f32; nf];
    loop {
        rt[..nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        rt[nt as usize..].fill(0.0);
        pfarc(1, &rt, &mut ct);
        for i in 0..nf {
            amp[i] = ((ct[i].r * ct[i].r + ct[i].i * ct[i].i) as f64).sqrt() as f32;
            ph[i] = (a as f64 + b as f64 * (ct[i].i as f64).atan2(ct[i].r as f64) + (c * i as f32) as f64) as f32;
        }
        for i in 0..nf {
            ct[i].r = (amp[i] as f64 * (ph[i] as f64).cos()) as f32;
            ct[i].i = (amp[i] as f64 * (ph[i] as f64).sin()) as f32;
        }
        pfacr(-1, &ct, &mut rt);
        for i in 0..nt as usize { rt[i] *= onfft; }
        tr.data[..nt as usize].copy_from_slice(&rt[..nt as usize]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
