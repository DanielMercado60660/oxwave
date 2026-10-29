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

//! SUMINPHASE - convert input to minimum phase.
//! Port of su/main/filters/suminphase.c (Revision 1.2).
//! Credits: SEAM Project: Bruce VerWest c. 2013.
//! Spectral factorization author: Bruce VerWest 2009.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_pfa_fft::{npfao, pfacc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMINPHASE - convert input to minimum phase\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suminphase <stdin >stdout [optional parameters]\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t \t\t",
    "\tif dt is not set in header, then dt is mandatory\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tsign1=1\t\tsign of first transform\t(1 or -1)\t\t",
    "\tsign2=-1\tsign of second transform (-1 or 1)\t\t",
    "    \tpnoise=1.e-9\t   white noise in spectral routine\t\t",
    "\tverbose=0\t\t=1 for advisory messages\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:  Homomorphic wavelet estimation                              ",
    "\t\t\t\t\t\t\t\t\t",
    " suclogfft < shotgather | suamp mode=real | sustack key=dt > real.su   ",
    " suclogfft < shotgather | suamp mode=imag | sustack key=dt > imag.su   ",
    " suop2 real.su imag.su op=zipper | suiclogfft | suminphase > wavelet.su",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0i32;
    let mut sign1 = 1i32;
    let mut sign2 = -1i32;
    getparint("verbose", &mut verbose);
    getparint("sign1", &mut sign1);
    getparint("sign2", &mut sign2);
    if sign1 == 0 { err!("sign1=%d ", sign1); }
    if sign2 == 0 { err!("sign2=%d ", sign2); }
    sign1 = sign1.checked_div(sign1.wrapping_abs()).unwrap_or(0);
    sign2 = sign2.checked_div(sign2.wrapping_abs()).unwrap_or(0);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let mut dt: f32;
    if seismic {
        if verbose != 0 { warn!("input is seismic data, trid=%d", tr.trid); }
        dt = (tr.dt as f64 / 1000000.0) as f32;
    } else {
        if verbose != 0 { warn!("input is not seismic data, trid=%d", tr.trid); }
        dt = tr.d1;
    }
    if dt == 0.0 && getparfloat("dt", &mut dt) == 0 {
        dt = 0.004;
        if verbose != 0 {
            warn!("neither dt nor d1 are set, nor is dt getparred!");
            warn!("assuming .004 sec sampling!");
        }
    }
    let mut pnoise = 1.0e-9f32;
    getparfloat("pnoise", &mut pnoise);
    let nt = tr.ns as usize;
    let nfft = npfao(nt as i32, 2 * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 {
        err!("Padded nt=%d -- too big", nfft);
    }
    if verbose != 0 { warn!("nfft,nt,dt is =%d %d %f", nfft, nt as i32, dt); }
    let nfft = nfft as usize;
    let mut rt = vec![0.0f32; nfft];
    let mut ct = vec![Complex::default(); nfft];
    loop {
        rt[..nt].copy_from_slice(&tr.data[..nt]);
        rt[nt..].fill(0.0);
        for i in 0..nfft { ct[i] = cmplx(rt[i], 0.0); }
        pfacc(sign1, &mut ct);
        kolmogoroff(pnoise, sign1, sign2, &mut ct);
        pfacc(sign2, &mut ct);
        let scale = (1.0f64 / nfft as f64) as f32;
        for i in 0..nfft { ct[i] = crmul(ct[i], scale); }
        for i in 0..nfft { rt[i] = ct[i].r; }
        tr.data[..nt].copy_from_slice(&rt[..nt]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}

// Spectral factorization, Bruce VerWest 2009.
fn kolmogoroff(pnoise: f32, sign1: i32, sign2: i32, cx: &mut [Complex]) {
    let n = cx.len();
    let mut rmax = 0.0f32;
    for i in 0..n {
        let r = rcabs(cx[i]);
        if r > rmax { rmax = r; }
    }
    for i in 0..n {
        let mut r = rcabs(cx[i]) / rmax;
        if r < pnoise { r = pnoise; }
        cx[i] = cmplx(r * r, 0.0);
        cx[i] = cwp_clog(cx[i]);
    }
    pfacc(sign2, cx);
    let scale = (1.0f64 / n as f64) as f32;
    for i in 0..n { cx[i] = crmul(cx[i], scale); }
    cx[0] = crmul(cx[0], 0.5);
    cx[n / 2] = crmul(cx[n / 2], 0.5);
    for i in n / 2 + 1..n { cx[i] = cmplx(0.0, 0.0); }
    pfacc(sign1, cx);
    for i in 0..n { cx[i] = crmul(cwp_cexp(cx[i]), rmax); }
}
