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

//! SUPHIDECON - phase inversion deconvolution, port of su/main/decon_shaping/suphidecon.c (Revision 1.9).
//! Credits: Potash Corporation, Saskatchewan, Balasz Nemeth; supplied to CWP in 2008 as supid.c.
//! Reference: Lichman and Northwood, SEG Deconvolution 2, Geophysics reprint series no. 17, pp. 701-718.

use oxwave_array_wrap::wrapArray;
use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::{closeto, nint};
use oxwave_hilbert::hilbert;
use oxwave_pfa_fft::{npfao, pfacc, pfacr, pfarc};
use oxwave_sblas::sscal;
use oxwave_windows::hanningnWindow;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUPHIDECON - PHase Inversion Deconvolution\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "    suphidecon < stdin > stdout\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t  \t",
    "\tnone\t\t\t\t\t\t\t   \t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " ... time range used for wavelet extraction:\t\t\t   \t",
    " tm=-0.1\tPre zero time (maximum phase component )\t\t",
    " tp=+0.4\tPost zero time (minimum phase component + multiples)    ",
    " pad=.1\tpercentage padding for nt prior to cepstrum calculation\t",
    "\t\t\t\t\t\t\t\t\t",
    " pnoise=0.001\tPre-withening (assumed noise to prevent division by zero)",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The wavelet is separated from the reflectivity and noise based on\t",
    " their different 'smoothness' in the pseudo cepstrum domain.\t\t",
    " The extracted wavelet also includes multiples. \t\t\t",
    " The wavelet is reconstructed in frequency domain, end removed\t\t",
    " from the trace. (Method by Lichman and Northwood, 1996.)\t\t",
];

fn compute_pseudo_cepstrum(nt: usize, x: &[f32], c: &mut [Complex], init: bool) -> usize {
    let fftpadd = 50.0f32;
    let ntp = nint(nt as f64 * (1.0 + fftpadd as f64 / 100.0));
    let nfftc = npfao(ntp, 4 * ntp);
    if nfftc >= SU_NFLTS as i32 || nfftc >= 720720 {
        err!("Padded nt=%d--too big", nfftc);
    }
    let n = nfftc as usize;
    let nf = n / 2 + 1;
    let snfftc = (1.0 / nfftc as f64) as f32;
    if init {
        return n;
    }
    let mut w = vec![Complex::default(); n];
    let mut a = vec![0.0f32; n];
    let mut p = vec![0.0f32; n];
    a[..nt].copy_from_slice(&x[..nt]);
    sscal(nt as i32, snfftc, &mut a, 1);
    pfarc(-1, &a, &mut w);
    for iw in 0..nf {
        a[iw] = w[iw].r;
    }
    hilbert(&a[..nf], &mut p[..nf]);
    for iw in 0..nf {
        w[iw].i = ((-1.0 / std::f64::consts::PI) * p[iw] as f64) as f32;
    }
    for iw in 0..nf {
        if !closeto(rcabs(w[iw]), 0.0) {
            w[iw] = cwp_clog(w[iw]);
        } else {
            w[iw] = cmplx(0.0, 0.0);
        }
    }
    for iw in 0..nf - 2 {
        w[nf + iw] = cmplx(w[nf - iw - 2].r, -w[nf - iw - 2].i);
    }
    pfacc(1, &mut w);
    c.copy_from_slice(&w);
    n
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tm = -0.1f32;
    let mut tp = 0.4f32;
    let mut pad = 0.1f32;
    getparfloat("tm", &mut tm);
    getparfloat("tp", &mut tp);
    getparfloat("pad", &mut pad);
    let mut tr = Trace::new();
    let mut trout = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let itmin = (-tp / dt) as i32;
    let itmax = (-tm / dt) as i32;
    // C leaves the window unchecked. Reject overflowing lengths and windows
    // beyond c before allocating or indexing; preserve its malloc diagnostic
    // for negative lengths and its valid zero-length windows.
    let nwind = itmax
        .checked_sub(itmin)
        .and_then(|n| n.checked_add(1))
        .unwrap_or_else(|| err!("time window length overflows"));
    if nwind < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let nwind = nwind as usize;
    let nfft = compute_pseudo_cepstrum(nt, &tr.data, &mut [], true);
    if nwind > nfft {
        err!("time window length=%d exceeds nfft=%d", nwind, nfft);
    }
    // wrapArray reduces any representable shift modulo nfft. Checking the
    // addition also ensures its result can safely be negated for unwrapping.
    let shift = itmax
        .checked_add(1)
        .unwrap_or_else(|| err!("time window wrap offset overflows"));
    let mut window = vec![0.0f32; nwind];
    hanningnWindow(&mut window);
    let mut c = vec![Complex::default(); nfft];
    let mut ct = vec![Complex::default(); nfft];
    loop {
        compute_pseudo_cepstrum(nt, &tr.data, &mut c, false);
        wrapArray(&mut c, shift);
        for it in 0..nwind {
            c[it] = crmul(c[it], window[it]);
        }
        c[nwind..].fill(Complex::default());
        wrapArray(&mut c, -shift);
        pfacc(-1, &mut c);
        let scale = (1.0 / nfft as f64) as f32;
        // Complex has the same repr(C) pair of floats as SU.
        let c_floats =
            unsafe { std::slice::from_raw_parts_mut(c.as_mut_ptr() as *mut f32, nfft * 2) };
        sscal((nfft * 2) as i32, scale, c_floats, 1);
        for v in &mut c {
            let tmpc = cwp_cexp(*v);
            *v = if tmpc.r.is_finite() && tmpc.i.is_finite() {
                tmpc
            } else {
                cmplx(0.0, 0.0)
            };
        }
        tr.data[nt..nfft].fill(0.0);
        pfarc(-1, &tr.data[..nfft], &mut ct);
        let ct_floats =
            unsafe { std::slice::from_raw_parts_mut(ct.as_mut_ptr() as *mut f32, nfft * 2) };
        sscal((nfft * 2) as i32, scale, ct_floats, 1);
        let nf = nfft / 2 + 1;
        let mut spsum = 0.0f32;
        for v in &c[..nf] {
            spsum += rcabs(*v);
        }
        pad = (pad as f64 / 100.0 * (spsum / nf as f32) as f64) as f32;
        for it in 0..nf - 2 {
            ct[nf + it] = cmplx(ct[nf - it - 2].r, -ct[nf - it - 2].i);
        }
        for it in 0..nfft {
            ct[it] = cdiv(ct[it], cmplx(c[it].r + pad, c[it].i));
        }
        pfacr(1, &ct, &mut trout.data[..nfft]);
        trout.hdr = tr.hdr.clone();
        puttr(&mut trout);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
