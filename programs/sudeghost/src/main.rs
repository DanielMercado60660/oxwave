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

//! SUDEGHOST - applies a deghosting filter.
//! Port of su/main/filters/sudeghost.c (Revision 1.7).
//! Credits: CWP: John Stockwell (November 2021).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{Complex, cmul};
use oxwave_pfa_fft::{npfaro, pfarc, pfacr};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUDEGHOST - applies a DEGHOSTING filter in (t,x) or (tau,p) data\t",
    "\t\t\t\t\t\t\t\t\t",
    " sudeghost <stdin >stdout [optional parameters]         \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:                                         \t\t",
    " h=\t\t\tsource or receiver depth\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " v=1500.0\t\tspeed of sound in the top layer\t\t\t",
    " r=0.5\t\t\tsurface reflectivity 0 < r < .8\t\t\t",
    " lambert=1 \t\tLambert's cosine law obliquity factor\t\t",
    " dt= (from header)     time sampling interval (sec)        \t\t",
    " verbose=0\t\t=1 for advisory messages, =2 debugging\t\t",
    " deghost=1\t\tdeghosting filter; =0 ghosting filter for modeling",
    " pnoise=1.e-9\t\twhite noise parameter\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The input data are assumed to be shot gathers with no missing traces.\t",
    "\t\t\t\t\t\t\t\t\t",
    " If the input data are in the (x,t) domain, and the header fields\t",
    " f2=0 and d2=0, then the program assumes that the ghost signal is normally",
    " incident and the ghost delay time is 2h/v.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " If the input data are in the (tau,p) domain, then it is assumed that  ",
    " f2='first p value' and d2='increment in p' and the ghost delay time",
    " is assumed to be (2h/v) sqrt{1 - v^2 p^2}, addressing the angular\t",
    " dependence of the ghost delay. \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples: \t\t\t\t\t\t\t\t",
    " (t,x) domain, streamer depth 10m, surface reflectivity r=0.5\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  sudeghost < input h=10 r=0.5 > output\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " (tau,p) domain, streamer depth 10m surface reflectivity r=0.5\t\t",
    " min offset=-636 max offset -3237 \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suradon < input choose=0 igopt=3 interoff=-262 offref=-3237 \t\t",
    "     pmin=-400 pmax=2000 dp=10 f1=90 f2=125 cdpkey=ep anderson=0  |\t",
    "     sudeghost h=10 r=.5 deghost=1 |\t\t\t\t\t",
    " suradon choose=4 igopt=3 interoff=-262 offref=-3237\t\t\t",
    "   pmin=-400 pmax=2000 dp=10  f1=90 f2=125 cdpkey=ep anderson=0  > output",
    " \t\t\t\t\t\t\t\t\t",
    " Caveats:\t\t\t\t\t\t\t\t",
    " The value of r is the reflectivity of the sea surface, which may be   ",
    " both dependent on frequency and on the angle of incidence. The program",
    " may be unstable if r > .8 is chosen.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Smaller values of dp may be needed to retain high frequencies\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " For lambert=1 a cosine obliquity law is chosen. Frequency dependence  ",
    " of reflectivity is not addressed.\t\t \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " If tr.dt is not set in the header, then dt is mandatory\t\t",
];

fn deghosting_filter(lambert: i32, deghost: i32, pnoise: f32, verbose: i32,
    r: f32, h: f32, v: f32, p: f32, nfft: i32, dt: f32, filter: &mut [Complex]) {
    let nf = nfft / 2 + 1;
    let onfft = (1.0 / nfft as f64) as f32;
    let df = onfft / dt;
    let mut amp = 0.0f32;
    for ifreq in 0..nf {
        let mut vvpp = v * v * p * p;
        let f = ifreq as f32 * df;
        if vvpp > 1.0 { vvpp = 1.0; }
        if lambert == 1 {
            amp = (r as f64 * ((1.0f32 - vvpp) as f64).sqrt()) as f32;
            if amp < pnoise { amp = pnoise; }
        } else if lambert == 0 { amp = r; }
        let omega = (2.0 * std::f64::consts::PI * f as f64) as f32;
        let twowaytt = ((2.0f32 * h) as f64 * ((1.0f32 - vvpp) as f64).sqrt() / v as f64) as f32;
        let real = (1.0 - amp as f64 * ((omega * twowaytt) as f64).cos()) as f32;
        let imag = (amp as f64 * ((omega * twowaytt) as f64).sin()) as f32;
        let denom = real * real + imag * imag;
        if verbose == 2 {
            warn!("amp = %f lambert = %d deghost = %d f = %f df = %f omega = %f twowaytt = %f dt = %f ", amp, lambert, deghost, f, df, omega, twowaytt, dt);
        }
        if deghost == 1 {
            filter[ifreq as usize].r = onfft * real / denom;
            filter[ifreq as usize].i = onfft * imag / denom;
        } else if deghost == 0 {
            filter[ifreq as usize].r = onfft * real;
            filter[ifreq as usize].i = -onfft * imag;
        }
    }
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
    let fp = tr.f2;
    let dp = tr.d2;
    let mut h = 0.0f32;
    mustgetparfloat("h", &mut h);
    let mut r = 0.5f32; getparfloat("r", &mut r);
    let mut v = 1500.0f32; getparfloat("v", &mut v);
    let mut deghost = 1; getparint("deghost", &mut deghost);
    let mut lambert = 1; getparint("lambert", &mut lambert);
    let mut pnoise = 1.0e-9f32; getparfloat("pnoise", &mut pnoise);
    let nfft = npfaro(nt, 2 * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d -- too big", nfft); }
    let nf = (nfft / 2 + 1) as usize;
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nf];
    let mut filter = vec![Complex::default(); nf];
    let mut ntr = 0i32;
    loop {
        ntr += 1;
        let p = ntr as f32 * dp + fp;
        deghosting_filter(lambert, deghost, pnoise, verbose, r, h, v, p, nfft, dt, &mut filter);
        rt[..nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        rt[nt as usize..].fill(0.0);
        pfarc(1, &rt, &mut ct);
        for i in 0..nf { ct[i] = cmul(ct[i], filter[i]); }
        pfacr(-1, &ct, &mut rt);
        tr.data[..nt as usize].copy_from_slice(&rt[..nt as usize]);
        if verbose != 0 { warn!("r = %f nyq = %f tr.f2 = %f  tr.d2 = %f", r, nyq, tr.f2, tr.d2); }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
