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

//! SUZEROPHASE - convert input to zero phase equivalent.
//! Port of su/main/transforms/suzerophase.c (Revision 1.1).
//! Credits: c. 2011 Bruce VerWest as part of the SEAM project.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cmplx, conjg, crmul, rcabs};
use oxwave_pfa_fft::{npfao, pfacc, Complex};
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUZEROPHASE - convert input to zero phase equivalent\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suzerophase <stdin >stdout [optional parameters]\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t \t\t",
    "\tif dt is not set in header, then dt is mandatory\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "    \tt0=1.0\t\t\ttime of peak value t0 \t\t\t",
    "\tverbose=0\t\t=1 for advisory messages\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
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
    let mut t0 = 1.0f32;
    getparfloat("t0", &mut t0);
    let nt = tr.ns as usize;
    let nfft = npfao(nt as i32, 2 * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 {
        err!("Padded nt=%d -- too big", nfft);
    }
    let dom = (6.2831853f64 / nfft as f64 / dt as f64) as f32;
    // C omits the final floating-point argument in this advisory.
    if verbose != 0 { warn!("nfft,nt,dt,dom is =%d %d %f %f", nfft, nt as i32, dt, 0.0f64); }
    let nfft = nfft as usize;
    let mut rt = vec![0.0f32; nfft];
    let mut ct = vec![Complex::default(); nfft];
    loop {
        rt[..nt].copy_from_slice(&tr.data[..nt]);
        rt[nt..].fill(0.0);
        for i in 0..nfft { ct[i] = cmplx(rt[i], 0.0); }
        pfacc(1, &mut ct);
        for i in 0..nfft / 2 + 1 {
            let r = rcabs(ct[i]);
            let om = i as f32 * dom;
            let phase = (om * t0) as f64;
            ct[i] = cmplx((r as f64 * phase.cos()) as f32, (r as f64 * phase.sin()) as f32);
        }
        for i in nfft / 2 + 1..nfft { ct[i] = conjg(ct[nfft - i]); }
        pfacc(-1, &mut ct);
        let scale = (1.0f64 / nfft as f64) as f32;
        for i in 0..nfft { ct[i] = crmul(ct[i], scale); }
        for i in 0..nfft { rt[i] = ct[i].r; }
        tr.data[..nt].copy_from_slice(&rt[..nt]);
        tr.delrt = (-t0 * 1000.0) as i16;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
