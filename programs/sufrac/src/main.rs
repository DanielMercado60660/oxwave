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
//! SUFRAC - fractional time derivative, integral, and phase shift.
//! Port of su/main/filters/sufrac.c, revision 1.29.
//! Credits: Chris Liner, Jack K. Cohen, Dave Hale, Zhenyue Liu,
//! John Stockwell, and Werner M. Heigl.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cmplx, cmul, crmul, cwp_cexp, Complex};
use oxwave_pfa_fft::{npfaro, pfarc, pfacr};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUFRAC -- take general (fractional) time derivative or integral of\t",
    "\t    data, plus a phase shift.  Input is CAUSAL time-indexed\t",
    "\t    or depth-indexed data.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sufrac power= [optional parameters] <indata >outdata \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tpower=0\t\texponent of (-i*omega)\t \t\t\t",
    "\t\t\t=0  ==> phase shift only\t\t\t",
    "\t\t\t>0  ==> differentiation\t\t\t\t",
    "\t\t\t<0  ==> integration\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tsign=-1\t\t\tsign in front of i * omega\t\t",
    "\tdt=(from header)\ttime sample interval (in seconds)\t",
    "\tphasefac=0\t\tphase shift by phase=phasefac*PI\t",
    "\tverbose=0\t\t=1 for advisory messages\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "  preprocess to correct 3D data for 2.5D migration\t\t\t",
    "         sufrac < sudata power=.5 sign=1 | ...\t\t\t\t",
    "  preprocess to correct susynlv, susynvxz, etc. (2D data) for 2D migration",
    "         sufrac < sudata phasefac=.25 | ...\t\t\t\t",
    " The filter is applied in frequency domain.\t\t\t\t",
    " if dt is not set in header, then dt is mandatory\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Algorithm:\t\t\t\t\t\t\t\t",
    "\t\tg(t) = Re[INVFTT{ ( (sign) iw)^power FFT(f)}]\t\t",
    " Caveat:\t\t\t\t\t\t\t\t",
    " Large amplitude errors will result if the data set has too few points.",
    "\t\t\t\t\t\t\t\t\t",
    " Good numerical integration routine needed!\t\t\t\t",
    " For example, see Gnu Scientific Library.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut verbose = 0;
    let mut power = 0.0f32;
    let mut sign = -1.0f32;
    let mut phasefac = 0.0f32;
    getparint("verbose", &mut verbose);
    getparfloat("power", &mut power);
    getparfloat("sign", &mut sign);
    getparfloat("phasefac", &mut phasefac);
    let phase = (phasefac as f64 * std::f64::consts::PI) as f32;

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
    let nfft = npfaro(nt as i32, 2 * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 {
        err!("Padded nt=%d -- too big", nfft);
    }
    let nf = (nfft / 2 + 1) as usize;
    let onfft = (1.0f64 / nfft as f64) as f32;
    let domega = (2.0 * std::f64::consts::PI * onfft as f64 / dt as f64) as f32;
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nf];
    let mut filt = vec![Complex::default(); nf];
    let arg = (sign as f64 * (0.5 * std::f64::consts::PI) * power as f64 + phase as f64) as f32;
    let exparg = cwp_cexp(crmul(cmplx(0.0, 1.0), arg));
    for i in 0..nf {
        let mut omega = i as f32 * domega;
        if power < 0.0 && i == 0 { omega = f32::MAX; }
        let amp = ((omega as f64).powf(power as f64) * onfft as f64) as f32;
        filt[i] = crmul(exparg, amp);
    }
    loop {
        rt[..nt].copy_from_slice(&tr.data[..nt]);
        rt[nt..].fill(0.0);
        pfarc(1, &rt, &mut ct);
        for i in 0..nf { ct[i] = cmul(ct[i], filt[i]); }
        pfacr(-1, &ct, &mut rt);
        tr.data[..nt].copy_from_slice(&rt[..nt]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
