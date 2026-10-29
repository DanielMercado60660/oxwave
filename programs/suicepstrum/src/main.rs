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

//! SUICEPSTRUM - inverse cepstral transform.
//! Port of su/main/transforms/suicepstrum.c (Revision 1.2).
//! Credits (from the C source):
//! CWP: John Stockwell, Dec 2010 based on suifft.c by:
//! CWP: Shuki Ronen, Chris Liner, Jack K. Cohen, c. 1989

use oxwave_core::prelude::*;
use oxwave_pfa_fft::{npfaro, pfacr, pfarc, Complex};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUICEPSTRUM - fft of complex log frequency traces to real time traces",
    " \t\t\t\t\t\t\t\t",
    "  suicepstrum <stdin >sdout sign2=-1\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    " \tsign1=1\t\tsign in exponent of first fft\t\t",
    " \tsign2=-1\tsign in exponent of inverse fft\t\t",
    "\tsym=0\t\t=1 center  output \t\t\t",
    "\tdt=tr.dt\ttime sampling interval (s) from header\t",
    "\t\t\tif not set assumed to be .004s\t\t",
    " Output traces are normalized by 1/N where N is the fft size.\t",
    " \t\t\t\t\t\t\t\t",
    " Note:\t\t\t\t\t\t\t\t",
    " The forward  cepstral transform is the\t\t\t",
    "   F(t_c) = InvFFT[ln[FFT(F(t))]] \t\t\t\t",
    " The inverse  cepstral transform is the\t\t\t",
    "   F(t) = InvFFT[exp[FFT(F(t_c))]] \t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Here t_c is the cepstral time (quefrency) domain \t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let ns = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    if dt == 0.0 {
        if verbose != 0 {
            warn!("dt header field not set!");
        }
        dt = tr.d1;
    }
    if dt == 0.0 {
        dt = 0.004;
    }
    let nfft = npfaro(ns, 2 * ns);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 {
        err!("Padded ns=%d--too big", nfft);
    }
    let nf = nfft / 2 + 1;
    let mut sym = 0;
    let mut sign1 = 1;
    let mut sign2 = -1;
    getparint("sym", &mut sym);
    getparint("sign1", &mut sign1);
    getparint("sign2", &mut sign2);
    if verbose != 0 {
        warn!("sign1 = %d sign2 = %d", sign1, sign2);
    }
    checkpars();
    if sign1 != 1 && sign1 != -1 {
        err!("sign1 = %d must be 1 or -1", sign1);
    }
    if sign2 != -1 && sign2 != 1 {
        err!("sign2 = %d must be -1 or 1", sign2);
    }
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nf as usize];
    if dt != 0.0 && verbose != 0 {
        warn!("d1=%f", 1.0 / ((nfft as f32 * dt) as f64));
    }
    let newns = 2 * nf;
    let newnfft = newns - 2;
    let newnf = newnfft / 2 + 1;
    let onewnfft = (1.0 / newnfft as f64) as f32;
    let mut newrt = vec![0.0f32; newnfft as usize];
    loop {
        ct.fill(Complex::default());
        rt[..ns as usize].copy_from_slice(&tr.data[..ns as usize]);
        rt[ns as usize..].fill(0.0);
        pfarc(sign1, &rt, &mut ct);
        for i in 0..nf as usize {
            let temp = ct[i].r;
            if temp != 0.0 {
                ct[i].r = ((temp as f64).exp() as f32 as f64 * (ct[i].i as f64).cos()) as f32;
                ct[i].i = ((temp as f64).exp() as f32 as f64 * (ct[i].i as f64).sin()) as f32;
                if sym != 0 && i % 2 != 0 {
                    ct[i].r = -ct[i].r;
                    ct[i].i = -ct[i].i;
                }
            } else {
                ct[i].r = 0.0;
                ct[i].i = 0.0;
            }
        }
        newrt.fill(0.0);
        // C passes newnf, including odd lengths, and leaves the rest zeroed.
        pfacr(sign2, &ct, &mut newrt[..newnf as usize]);
        for i in 0..newnfft as usize {
            tr.data[i] = newrt[i] * onewnfft;
        }
        tr.trid = 1;
        tr.ns = newnfft as u16;
        tr.f1 = 0.0;
        tr.d1 = dt;
        tr.dt = ((dt as f64 * 1000000.0) as i32) as u16;
        if sym != 0 {
            // The C build converts to int, then stores the low 16 bits.
            // A direct Rust float-to-i16 cast would saturate instead.
            tr.delrt =
                (-(tr.delrt as f64 + (newnfft as f32 * tr.d1) as f64 * 1000.0 / 2.0)) as i32 as i16;
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
