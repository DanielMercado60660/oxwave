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
//! SULOG - time axis log-stretch of seismic traces.
//! Port of su/main/stretching_moveout_resamp/sulog.c (Revision 1.19).
//! Credits: CWP: Shuki, Chris.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SULOG -- time axis log-stretch of seismic traces\t\t",
    "\t\t\t\t\t\t\t\t",
    " sulog [optional parameters] <stdin >stdout \t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "\tnone\t\t\t\t \t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "\tntmin= .1*nt\t\tminimum time sample of interest\t",
    "\toutpar=/dev/tty\t\toutput parameter file, contains:",
    "\t\t\t\tnumber of samples (nt=)\t\t",
    "\t\t\t\tminimum time sample (ntmin=)\t",
    "\t\t\t\toutput number of samples (ntau=)",
    "\tm=3\t\t\tlength of stretched data\t",
    "\t\t\t\tis set according to\t\t",
    "\t\t\t\t\tntau = nextpow(m*nt)\t",
    "\tntau= pow of 2\t\toverride for length of stretched",
    "\t\t\t\tdata (useful for padding zeros\t",
    "\t\t\t\tto avoid aliasing)\t\t",
    "\t\t\t\t\t\t\t\t",
    " NOTES:\t\t\t\t\t\t\t",
    "\tntmin is required to avoid taking log of zero and to \t",
    "\tkeep number of outsamples (ntau) from becoming enormous.",
    "        Data above ntmin is zeroed out.\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "\tThe output parameters will be needed by suilog to \t",
    "\treconstruct the original data. \t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " EXAMPLE PROCESSING SEQUENCE:\t\t\t\t\t",
    "\t\tsulog outpar=logpar <data1 >data2\t\t",
    "\t\tsuilog par=logpar <data2 >data3\t\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn nextpower(p: i32, n: i32) -> i32 {
    if n == 0 { return 0; }
    // The next power of two cannot fit in an i32 above 2^30.
    // Reject it before the C loop would overflow and cease to terminate.
    if p == 2 && n > (1 << 30) { err!("ntau calculation overflow"); }
    let mut nn: i32 = 1;
    while nn < n {
        nn = nn.checked_mul(p).unwrap_or_else(|| err!("ntau calculation overflow"));
    }
    nn
}

fn lintrp(q: &[f32], w: &mut [f32], it: &mut [i32], lp: i32, lq: i32) {
    for i in 0..lq as usize {
        if q[i] >= 0.0 && q[i] < (lp - 1) as f32 {
            it[i] = q[i] as i32;
            let delta = q[i] - it[i] as f32;
            w[i*2] = (1.0f64 - delta as f64) as f32;
            w[i*2+1] = delta;
        } else {
            it[i] = 0;
            w[i*2] = 0.0;
            w[i*2+1] = 0.0;
        }
    }
}

fn stretch(q: &mut [f32], p: &[f32], w: &[f32], it: &[i32], lq: i32, nw: usize) {
    for i in 0..lq as usize {
        q[i] = 0.0;
        for j in 0..nw { q[i] += w[i*nw+j] * p[it[i] as usize+j]; }
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let mut ntmin = 0;
    if getparint("ntmin", &mut ntmin) == 0 { ntmin = (0.1 * nt as f64) as i32; }
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let outparfp = efopen(&outpar, "w");
    let nw = 2;
    let mut m = 0;
    if getparint("m", &mut m) == 0 { m = 3; }
    let ntautemp = m.checked_mul(nt).unwrap_or_else(|| err!("m*nt overflow"));
    let dtau = (((nt as f32 / ntmin as f32) as f64).ln() / ntautemp as f64) as f32;
    let mut ntau = 0;
    if getparint("ntau", &mut ntau) == 0 { ntau = nextpower(2, ntautemp); }
    checkpars();
    CHECK_NT("ntau", ntau);
    if ntau < 0 { err!("ntau must not be negative"); }
    fprintf!(outparfp, "nt=%d ntmin=%d dt=%g ntau=%d\nt=%d n1=%d d1=%g", nt, ntmin, dt as f64 * 1000000.0, ntau, nt, nt, dt as f64 * 1000000.0);
    let mut t = vec![0.0f32; ntau as usize];
    let mut it = vec![0i32; ntau as usize];
    let mut w = vec![0.0f32; nw * ntau as usize];
    let mut buf = vec![0.0f32; ntau as usize];
    for itau in 0..ntau as usize {
        t[itau] = (ntmin as f64 * ((itau as f32 * dtau) as f64).exp()) as f32;
    }
    lintrp(&t, &mut w, &mut it, nt, ntau);
    loop {
        stretch(&mut buf, &tr.data, &w, &it, ntau, nw);
        tr.data[..ntau as usize].copy_from_slice(&buf);
        tr.ns = ntau as u16;
        tr.dt = (dtau as f64 * 1000000.0) as u32 as u16;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
