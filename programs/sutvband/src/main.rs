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
//! SUTVBAND - time-variant bandpass filter. Port of su/main/filters/sutvband.c (Revision 1.18).
//! Credits: CWP: Jack, Ken.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{Complex, crmul};
use oxwave_cwp::cwp::nint;
use oxwave_pfa_fft::{npfaro, pfarc, pfacr};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                               ",
    " SUTVBAND - time-variant bandpass filter (sine-squared taper)  ",
    "                                                               ",
    " sutvband <stdin >stdout tf= f=\t\t\t        ",
    "                                                               ",
    " Required parameters:                                          ",
    "       dt = (from header)      time sampling interval (sec)    ",
    "       tf=             times for which f-vector is specified   ",
    "       f=f1,f2,f3,f4   Corner frequencies corresponding to the ",
    "                       times in tf. Specify as many f= as      ",
    "                       there are entries in tf.                ",
    "                                                               ",
    " The filters are applied in frequency domain.                  ",
    "                                                               ",
    " Example:                                                      ",
    " sutvband <data tf=.2,1.5 f=10,12.5,40,50 f=10,12.5,30,40 | ...",
];

fn makefilter(f: &[f32], nfft: i32, nfreq: i32, dt: f32, filter: &mut [f32]) {
    let onfft = (1.0 / nfft as f64) as f32;
    let df = onfft / dt;
    let if1 = nint((f[0] / df) as f64);
    let if2 = nint((f[1] / df) as f64);
    let if3 = nint((f[2] / df) as f64).min(nfreq - 1);
    let if4 = nint((f[3] / df) as f64).min(nfreq - 1);
    // The C implementation indexes outside filter for these bins.
    if !(0..nfreq).contains(&if1) || !(0..nfreq).contains(&if2)
        || !(0..nfreq).contains(&if3) || !(0..nfreq).contains(&if4) {
        err!("filter frequency bin outside valid range");
    }
    let c = (1.57079632679490 / (if2 - if1 + 2) as f64) as f32;
    for i in if1..=if2 {
        let s = ((c * (i - if1 + 1) as f32) as f64).sin() as f32;
        filter[i as usize] = s * s * onfft;
    }
    let c = (1.57079632679490 / (if4 - if3 + 2) as f64) as f32;
    for i in if3..=if4 {
        let s = ((c * (if4 - i + 1) as f32) as f64).sin() as f32;
        filter[i as usize] = s * s * onfft;
    }
    for i in if2+1..if3 { filter[i as usize] = onfft; }
    for i in 0..if1 { filter[i as usize] = 0.0; }
    for i in if4+1..nfreq { filter[i as usize] = 0.0; }
}

fn bandpass(data: &[f32], nt: usize, nfft: usize, nfreq: usize, filter: &[f32], ftrace: &mut [f32]) {
    let mut rt = vec![0.0; nfft];
    let mut ct = vec![Complex::default(); nfreq];
    rt[..nt].copy_from_slice(&data[..nt]);
    pfarc(1, &rt, &mut ct);
    for i in 0..nfreq { ct[i] = crmul(ct[i], filter[i]); }
    pfacr(-1, &ct, &mut rt);
    ftrace[..nt].copy_from_slice(&rt[..nt]);
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    if !(0..=17).contains(&tr.trid) && tr.trid != 109 { warn!("input is not seismic data, trid=%d", tr.trid); }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if dt == 0.0 { err!("dt field is zero and not getparred"); }
    let tmin = (tr.delrt as f64 / 1000.0) as f32;
    let nfilter = countparval("tf");
    if nfilter == 0 { mustgetparfloat("tf", &mut 0.0); }
    if countparname("f") != nfilter { err!("must give one f 4-tuple for each (%d) tf value", nfilter); }
    let mut tf = vec![0.0f32; nfilter as usize];
    let mut itf = vec![0i32; nfilter as usize + 4];
    getparfloat("tf", &mut tf);
    let mut jmin = 2usize;
    let mut jmax = nfilter as usize + 1;
    for j in jmin..=jmax {
        itf[j] = nint(((tf[j-2] - tmin) / dt) as f64);
        // Outside this range the C implementation indexes outside data.
        if !(0..nt).contains(&itf[j]) { err!("filter time outside trace"); }
    }
    let nfft = npfaro(nt, 2 * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d -- too big", nfft); }
    let nfreq = nfft / 2 + 1;
    let mut filter = vec![vec![0.0f32; nfreq as usize]; nfilter as usize + 4];
    for j in jmin..=jmax {
        let mut f = vec![0.0f32; 4];
        if getnparfloat((j-jmin+1) as i32, "f", &mut f) != 4 { err!("must give 4 corner frequencies in f="); }
        if f[0] < 0.0 || f[0] > f[1] || f[1] >= f[2] || f[2] > f[3] { err!("Filter #%d has bad frequencies", (j-jmin+1) as i32); }
        makefilter(&f, nfft, nfreq, dt, &mut filter[j]);
    }
    if itf[jmin] > 0 { filter[jmin-1] = filter[jmin].clone(); itf[jmin-1] = 0; jmin -= 1; }
    if itf[jmax] < nt-1 { filter[jmax+1] = filter[jmax].clone(); itf[jmax+1] = nt-1; jmax += 1; }
    itf[jmin-1] = 0;
    itf[jmax+1] = nt-1;
    let mut ftrace = vec![vec![0.0f32; nt as usize]; nfilter as usize + 4];
    let mut data = vec![0.0f32; nt as usize];
    loop {
        for j in jmin..=jmax {
            data.fill(0.0);
            for i in itf[j-1]..=itf[j+1] { data[i as usize] = tr.data[i as usize]; }
            bandpass(&data, nt as usize, nfft as usize, nfreq as usize, &filter[j], &mut ftrace[j]);
        }
        for j in jmin..jmax {
            let fitfj = itf[j] as f32;
            for i in itf[j]..=itf[j+1] {
                let a = (i as f32 - fitfj) / (itf[j+1] as f32 - fitfj);
                tr.data[i as usize] = (1.0f32-a)*ftrace[j][i as usize] + a*ftrace[j+1][i as usize];
            }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
