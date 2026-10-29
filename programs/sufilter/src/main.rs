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

//! SUFILTER - zero-phase sine-squared tapered filter.
//! Port of su/main/filters/sufilter.c (Revision 1.23).
//! Credits: CWP: John Stockwell, Jack Cohen; CENPET: Werner M. Heigl.
use oxwave_core::prelude::*;
use oxwave_cwp::complex::crmul;
use oxwave_pfa_fft::{npfaro, pfarc, pfacr, Complex};
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUFILTER - applies a zero-phase, sine-squared tapered filter\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sufilter <stdin >stdout [optional parameters]         \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:                                         \t\t",
    "       if dt is not set in header, then dt is mandatory        \t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "       f=f1,f2,...             array of filter frequencies(HZ) \t",
    "       amps=a1,a2,...          array of filter amplitudes\t\t",
    "       dt = (from header)      time sampling interval (sec)        \t",
    "\tverbose=0\t\t=1 for advisory messages\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Defaults:f=.10*(nyquist),.15*(nyquist),.45*(nyquist),.50*(nyquist)\t",
    "                        (nyquist calculated internally)\t\t",
    "          amps=0.,1.,...,1.,0.  trapezoid-like bandpass filter\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples of filters:\t\t\t\t\t\t\t",
    " Bandpass:   sufilter <data f=10,20,40,50 | ...\t\t\t",
    " Bandreject: sufilter <data f=10,20,30,40 amps=1.,0.,0.,1. | ..\t",
    " Lowpass:    sufilter <data f=10,20,40,50 amps=1.,1.,0.,0. | ...\t",
    " Highpass:   sufilter <data f=10,20,40,50 amps=0.,0.,1.,1. | ...\t",
    " Notch:      sufilter <data f=10,12.5,35,50,60 amps=1.,.5,0.,.5,1. |..\t",
];

fn polygonal_filter(f: &[f32], amps: &[f32], nfft: i32, dt: f32, filter: &mut [f32]) {
    let npoly = f.len();
    let nf = nfft / 2 + 1;
    let onfft = (1.0 / nfft as f64) as f32;
    let df = onfft / dt;
    let mut intfr = vec![0i32; npoly];
    for ifs in 0..npoly {
        let v = f[ifs] / df;
        intfr[ifs] = (v as f64 + if v >= 0.0 { 0.5 } else { -0.5 }) as i32;
        if intfr[ifs] > nf - 1 { intfr[ifs] = nf - 1; }
    }
    let mut taper = 0;
    for i in 0..intfr[0] { filter[i as usize] = amps[0] * onfft; }
    for j in 0..npoly-1 {
        if amps[j] < amps[j+1] {
            taper += 1;
            for i in intfr[j]..=intfr[j+1] {
                let c = (1.57079632679490 / (intfr[j+1] - intfr[j] + 2) as f64) as f32;
                let s = ((c * (i - intfr[j] + 1) as f32) as f64).sin() as f32;
                let adiff = amps[j+1] - amps[j];
                filter[i as usize] = (amps[j] + adiff*s*s) * onfft;
            }
        } else if amps[j] > amps[j+1] {
            taper += 1;
            for i in intfr[j]..=intfr[j+1] {
                let c = (1.57079632679490 / (intfr[j+1] - intfr[j] + 2) as f64) as f32;
                let s = ((c * (intfr[j+1] - i + 1) as f32) as f64).sin() as f32;
                let adiff = amps[j] - amps[j+1];
                filter[i as usize] = (amps[j+1] + adiff*s*s) * onfft;
            }
        } else {
            let start = intfr[j] + if taper == 0 { 0 } else { 1 };
            for i in start..=intfr[j+1] { filter[i as usize] = amps[j] * onfft; }
        }
    }
    for i in intfr[npoly-1]+1..nf { filter[i as usize] = amps[npoly-1] * onfft; }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let mut dt = if seismic {
        if verbose != 0 { warn!("input is seismic data, trid=%d", tr.trid); }
        (tr.dt as f64 / 1000000.0) as f32
    } else {
        if verbose != 0 { warn!("input is not seismic data, trid=%d", tr.trid); }
        tr.d1
    };
    if dt == 0.0 && getparfloat("dt", &mut dt) == 0 {
        dt = 0.004;
        if verbose != 0 {
            warn!("neither dt nor d1 are set, nor is dt getparred!");
            warn!("assuming .004 sec sampling!");
        }
    }
    let nt = tr.ns as i32;
    let nyq = (0.5 / dt as f64) as f32;
    let nfft = npfaro(nt, 2 * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d -- too big", nfft); }
    let nf = (nfft / 2 + 1) as usize;
    let npoly = countparval("f");
    let f: Vec<f32> = if npoly != 0 {
        let mut f = vec![0.0; npoly as usize];
        getparfloat("f", &mut f);
        f
    } else {
        [0.10, 0.15, 0.45, 0.50].iter().map(|x| (x * nyq as f64) as f32).collect()
    };
    let npoly = f.len();
    if npoly < 2 { warn!("Only %d value defining filter", npoly as i32); }
    for i in 0..npoly-1 {
        if f[i] < 0.0 || f[i] > f[i+1] { err!("Bad filter parameters"); }
    }
    let namps = countparval("amps");
    let amps = if namps != 0 {
        let mut amps = vec![0.0; namps as usize];
        getparfloat("amps", &mut amps);
        amps
    } else {
        let mut amps = vec![1.0; npoly];
        amps[0] = 0.0; amps[npoly-1] = 0.0;
        amps
    };
    if amps.len() != npoly { err!("number of f values must = number of amps values"); }
    let mut icount = 0;
    for &amp in &amps {
        if amp > 0.0 { icount += 1; }
        if amp < 0.0 { err!("amp values must be positive"); }
    }
    if icount == 0 { err!("All amps values are zero"); }
    icount = 0;
    for i in 0..npoly-1 { if amps[i] != amps[i+1] { icount += 1; } }
    if icount == 0 { warn!("All amps values are the same"); }
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nf];
    let mut filter = vec![0.0f32; nf];
    polygonal_filter(&f, &amps, nfft, dt, &mut filter);
    loop {
        rt[..nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        rt[nt as usize..].fill(0.0);
        pfarc(1, &rt, &mut ct);
        for i in 0..nf { ct[i] = crmul(ct[i], filter[i]); }
        pfacr(-1, &ct, &mut rt);
        tr.data[..nt as usize].copy_from_slice(&rt[..nt as usize]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
