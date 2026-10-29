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

//! SUACORFRAC - fractional autocorrelation/convolution.
//! Port of su/main/convolution_correlation/suacorfrac.c (Revision 1.3).
//! Credits: UHouston: Chris Liner, Sept 2009; CWP: based on Hale's crpow.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cmplx, cmul, Complex};
use oxwave_pfa_fft::{npfaro, pfarc, pfacr};
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUACORFRAC -- general FRACtional Auto-CORrelation/convolution\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suacorfrac power= [optional parameters] <indata >outdata \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " a=0\t\t\texponent of complex amplitude\t \t\t",
    " b=0\t\t\tmultiplier of complex phase\t \t\t",
    " dt=(from header)\ttime sample interval (in seconds)\t\t",
    " verbose=0\t\t=1 for advisory messages\t\t\t",
    " ntout=tr.ns\t\tnumber of time samples output\t\t\t",
    " sym=0\t\t\tif non-zero, produce a symmetric output from\t",
    "\t\t\tlag -(ntout-1)/2 to lag +(ntout-1)/2\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The calculation is performed in the frequency domain.\t\t\t",
    " The fractional autocorrelation/convolution is obtained by raising\t",
    " Fourier coefficients to seperate real powers \t\t\t\t",
    "\t\t(a,b) for amp and phase:\t\t\t\t",
    "\t\t     Aout exp[-i Pout] = Ain Ain^a exp[-i (1+b) Pin] \t",
    "\t\twhere A=amplitude  P=phase.\t\t\t\t",
    " Some special cases:\t\t\t\t\t\t\t",
    "\t\t(a,b)=(1,1)\t-->\tauto-correlation\t\t",
    "\t\t(a,b)=(0.5,0.5)\t-->\thalf-auto-correlation\t\t",
    "\t\t(a,b)=(0,0)\t-->\tno change to data\t\t",
    "\t\t(a,b)=(0.5,-0.5)-->\thalf-auto-convolution\t\t",
    "\t\t(a,b)=(1,-1)\t-->\tauto-convolution\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn dopow(u: Complex, a: f32, b: f32) -> Complex {
    if a == 0.0 { return cmplx(1.0, 0.0); }
    if u.r == 0.0 && u.i == 0.0 { return cmplx(0.0, 0.0); }
    let ur = u.r;
    let ui = u.i;
    let amp = (0.5 * a as f64 * ((ur * ur + ui * ui) as f64).ln()).exp() as f32;
    let phs = (-(b as f64) * (ui as f64).atan2(ur as f64)) as f32;
    cmplx((amp as f64 * (phs as f64).cos()) as f32,
          (amp as f64 * (phs as f64).sin()) as f32)
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    let mut a = 0.0f32;
    let mut b = 0.0f32;
    let mut sym = 0;
    getparint("verbose", &mut verbose);
    getparfloat("a", &mut a);
    getparfloat("b", &mut b);
    getparint("sym", &mut sym);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let dt = if seismic {
        if verbose != 0 { warn!("input is seismic data, trid=%d", tr.trid); }
        (tr.dt as f64 / 1000000.0) as f32
    } else {
        if verbose != 0 { warn!("input is not seismic data, trid=%d", tr.trid); }
        tr.d1
    };
    if dt == 0.0 { err!("dt or d1 field is zero and not getparred"); }
    let nt = tr.ns as usize;
    let mut ntout = tr.ns as i32;
    getparint("ntout", &mut ntout);
    checkpars();
    let nfft = npfaro(nt as i32, 2 * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d -- too big", nfft); }
    let nf = nfft as usize / 2 + 1;
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nf];
    loop {
        rt[..nt].copy_from_slice(&tr.data[..nt]);
        rt[nt..].fill(0.0);
        pfarc(1, &rt, &mut ct);
        for i in 0..nf {
            let filt = dopow(ct[i], a, b);
            ct[i] = cmul(ct[i], filt);
            if sym != 0 && i % 2 != 0 {
                ct[i].r = -ct[i].r;
                ct[i].i = -ct[i].i;
            }
        }
        pfacr(-1, &ct, &mut rt);
        tr.data[..nt].copy_from_slice(&rt[..nt]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
