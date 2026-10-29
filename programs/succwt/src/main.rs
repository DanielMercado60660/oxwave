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

//! SUCCWT - complex continuous wavelet transform, port of su/main/transforms/succwt.c.
//! Copyright (c) University of Tulsa, 2003-4. All rights reserved.
//! Authors: Chris Liner and Bob Clapp.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_pfa_fft::{npfa, pfacc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUCCWT - Complex continuous wavelet transform of seismic traces\t",
    "\t\t\t\t\t\t\t\t\t",
    " succwt < tdata.su > tfdata.su\t[optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " None\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " noct=5\tNumber of octaves (int)\t\t\t\t\t",
    " nv=10\t\tNumber of voices per octave (int)\t\t\t",
    " fmax=Nyq\tHighest frequency in transform\t\t\t\t",
    " p=-0.5\tPower of scale value normalizing CWT\t\t\t",
    "\t\t=0 for amp-preserved spec. decomp.\t\t\t",
    " c=1/(2*fmax)\tTime-domain inverse gaussian damping parameter\t\t",
    "\t\t(bigger c means more wavelet oscillations,\t\t",
    "\t\tdefault gives minimal oscillations)\t\t\t",
    " k=1\t\tUse complex Morlet as wavelet transform kernel\t\t",
    "\t\t=2 use Fourier kernel ... Exp[i 2 pi f t]\t\t",
    " fs=1\t\tUse dyadic freq sampling (CWT standard, honors\t\t",
    "\t\tnoct, nv)\t\t\t\t\t\t",
    "\t\t=2 use linear freq sampling (Fourier standard)\t\t",
    " df=1\t\tFrequency sample interval in Hz (used only for fs=2)\t",
    "\t\t\tNOTE: not yet implimented (hardwired to df=1) \t",
    " dt=(from tr.dt)\tSample interval override (in secs, if time data)",
    " verbose=0\t Run silent, except echo c value. (=1 for more info)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    " This generates amplitude spec of the CWT impulse response (IR).\t",
    "  suspike ntr=1 ix1=1 nt=125 | succwt | suamp | suximage & \t\t",
    " Real part of Fourier IR with linear freq sampling:\t\t\t",
    " suspike ntr=1 ix1=1 nt=125 | succwt k=2 fs=2 | suamp mode=real | suximage &",
    " Real part of Fourier IR with dyadic freq sampling: \t\t\t",
    " suspike ntr=1 ix1=1 nt=125 | succwt k=2 | suamp mode=real | suximage &",
    "\t\t\t\t\t\t\t\t\t",
    " Inverse CWT: (within a constant scale factor)\t\t\t\t",
    "\t... | succwt p=-1 | suamp mode=real | sustack key=cdp > inv.su\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " 1. Total number of scales: nscale = noct*nv\t\t\t\t",
    " 2. Each input trace spawns nscale complex output traces\t\t",
    " 3. Lowest frequency in the transform is fmax/( 2^(noct-1/nv) )\t",
    " 4. Header field (cdp) used as cwt spectrum counter\t\t\t",
    " 5. Header field (cdpt) used as scale counter within cwt spectrum\t",
    " 6. Header field (gut) holds number of cwt scales `na'\t\t\t",
    " 7. Header field (unscale) holds CWT scale `a'\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Header fields set: tracl, cdp, cdpt, unscale, gut\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn wavelet(a: f32, fmax: f32, c: f32, p: f32, dt: f32, nt: usize, w: &mut [Complex], k: i32) {
    let tc = (nt as i32 / 2 - 1) as f32 * dt;
    let cc = (a as f64).powf(p as f64) as f32;
    let fc = fmax / a;
    for it in 1..=nt {
        let t = (it - 1) as f32 * dt;
        let arg = (-2.0 * std::f64::consts::PI * fc as f64 * (t - tc) as f64) as f32;
        let exparg = cwp_cexp(crmul(cmplx(0.0, 1.0), arg));
        let rexp = if k == 2 {
            cc
        } else {
            let mut arg = (t - tc) / (a * c);
            arg = -arg * arg;
            (cc as f64 * (arg as f64).exp()) as f32
        };
        w[it] = crmul(exparg, rexp);
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut noct, mut nv, mut verbose, mut out, mut fs): (i32, i32, i32, i32, i32) =
        (5, 10, 0, 1, 1);
    let mut p = -0.5f32;
    getparint("noct", &mut noct);
    getparint("nv", &mut nv);
    getparint("verbose", &mut verbose);
    getparint("out", &mut out);
    getparint("fs", &mut fs);
    getparfloat("p", &mut p);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let npad = npfa(nt as i32) as usize;
    let mut dt = (tr.dt as f64 / 1000000.0) as f32;
    getparfloat("dt", &mut dt);
    let mut na = noct.wrapping_mul(nv);
    let mut fmax = 1.0f32 / (2.0f32 * dt);
    getparfloat("fmax", &mut fmax);
    let mut c = 1.0f32 / (2.0f32 * fmax);
    getparfloat("c", &mut c);
    warn!("c = %g", c);
    let mut df = 1.0f32;
    getparfloat("df", &mut df);
    if fs == 2 {
        na = ((fmax as i32) as f32 / df) as i32;
    }
    let mut k = 1;
    getparint("k", &mut k);
    if k != 1 && k != 2 {
        k = 1;
        warn!("Bad kernal wavelet choice ... use Morlet");
    }
    checkpars();
    if na < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let na = na as usize;
    if na > 0 && nt * 2 > SU_NFLTS {
        err!("complex output trace exceeds trace buffer");
    }
    let mut a = Vec::new();
    if a.try_reserve_exact(na + 1).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    a.resize(na + 1, 0.0f32);
    let mut ctr = vec![Complex::default(); npad];
    // The original kernel writes w[npad], one past its allocation.
    let mut ctmp = vec![Complex::default(); npad + 1];
    let mut basis = Vec::new();
    if basis.try_reserve_exact(na).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    for _ in 0..na {
        let mut row = Vec::new();
        if row.try_reserve_exact(npad + 1).is_err() {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        row.resize(npad + 1, Complex::default());
        basis.push(row);
    }
    if fs == 2 {
        for ia in 0..na {
            a[ia] = na as f32 / (ia + 1) as f32;
            if verbose == 1 {
                warn!("check a again %d %f", ia as i32, a[ia]);
            }
        }
    } else {
        let mut ia = 0;
        for ioct in 0..noct {
            for iv in 0..nv {
                let e = ioct as f32 + iv as f32 / nv as f32;
                a[ia] = 2.0f64.powf(e as f64) as f32;
                if verbose == 1 {
                    warn!("aval %d %f %f", ia as i32, a[ia], e);
                }
                ia += 1;
            }
        }
    }
    if verbose == 1 {
        warn!("fmax = %g, fmin = %g\n", fmax, fmax / a[na]);
        warn!("scale_index / scale_value / frequency");
        for ia in 0..na {
            warn!("\t%i / %g / %g", ia as i32, a[ia], fmax / a[ia]);
        }
    }
    for ia in 0..na {
        wavelet(a[ia], fmax, c, p, dt, npad, &mut ctmp, k);
        let ishift = npad / 2;
        for i in 0..=ishift {
            basis[ia][i] = ctmp[i + ishift];
            basis[ia][i + ishift] = ctmp[i];
        }
        pfacc(1, &mut basis[ia][..npad]);
    }
    let mut itracl = tr.tracl;
    let mut icdp = 1i32;
    let mut cwt = Trace::new();
    loop {
        let mut rms = 0.0f32;
        for it in 0..nt {
            ctr[it].r = tr.data[it];
            ctr[it].i = 0.0;
            rms = (rms as f64 + ((tr.data[it] * tr.data[it]) as f64).sqrt()) as f32;
        }
        if verbose == 1 {
            warn!("RMS IN %f", (rms as f64 / nt as f64).sqrt());
        }
        for it in nt..npad {
            ctr[it] = Complex::default();
        }
        pfacc(1, &mut ctr);
        for ia in 0..na {
            for i in 0..npad {
                ctmp[i] = cmul(basis[ia][i], ctr[i]);
            }
            pfacc(-1, &mut ctmp[..npad]);
            cwt.hdr = tr.hdr.clone();
            cwt.tracl = itracl;
            cwt.cdp = icdp;
            cwt.cdpt = ia as i32;
            cwt.unscale = a[ia];
            cwt.sut = noct as i16;
            cwt.gut = noct.wrapping_mul(nv) as i16;
            cwt.trid = 111;
            cwt.ns = (nt * 2) as u16;
            for i in 0..nt {
                cwt.data[2 * i] = ctmp[i].r;
                cwt.data[2 * i + 1] = ctmp[i].i;
            }
            puttr(&mut cwt);
            itracl = itracl.wrapping_add(1);
        }
        icdp += 1;
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
