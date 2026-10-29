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


//! SUILOG - time axis inverse log-stretch of seismic traces.
//! Port of su/main/stretching_moveout_resamp/suilog.c (Revision 1.17).
//! Credits: CWP: Shuki Ronen, Chris Liner.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUILOG -- time axis inverse log-stretch of seismic traces\t",
    "\t\t\t\t\t\t\t\t",
    " suilog nt= ntmin=  <stdin >stdout \t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnt= \tnt output from sulog prog\t\t\t",
    " \tntmin= \tntmin output from sulog prog\t\t\t",
    " \tdt= \tdt output from sulog prog\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " NOTE:  Parameters needed by suilog to reconstruct the \t",
    "\toriginal data may be input via a parameter file.\t",
    " \t\t\t\t\t\t\t\t",
    " EXAMPLE PROCESSING SEQUENCE:\t\t\t\t\t",
    " \t\tsulog outpar=logpar <data1 >data2\t\t",
    " \t\tsuilog par=logpar <data2 >data3\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " \twhere logpar is the parameter file\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

fn stretch(q: &mut [f32], p: &[f32], w: &[f32], it: &[i32], lq: usize, nw: usize) {
    for i in 0..lq {
        q[i] = 0.0;
        for j in 0..nw {
            q[i] += w[i * nw + j] * p[it[i] as usize + j];
        }
    }
}

fn lintrp(q: &[f32], w: &mut [f32], it: &mut [i32], lp: i32, lq: usize) {
    for i in 0..lq {
        if q[i] >= 0.0 && (q[i] as f64) < (lp - 1) as f64 {
            it[i] = q[i] as i32;
            let delta = q[i] - it[i] as f32;
            w[i * 2] = (1.0 - delta as f64) as f32;
            w[i * 2 + 1] = delta;
        } else {
            it[i] = 0;
            w[i * 2] = 0.0;
            w[i * 2 + 1] = 0.0;
        }
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let ntau = tr.ns as i32;
    let dtau = (tr.dt as f64 / 1000000.0) as f32;
    let mut nt = 0;
    let mut ntmin = 0;
    let mut dt = 0.0f32;
    mustgetparint("nt", &mut nt);
    CHECK_NT("nt", nt);
    mustgetparint("ntmin", &mut ntmin);
    mustgetparfloat("dt", &mut dt);
    checkpars();

    let nw = 2usize;
    if nt < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    // The C loops write outside tau for these values of ntmin.
    if ntmin < 0 || ntmin > nt {
        err!("ntmin must be between 0 and nt");
    }
    let n = nt as usize;
    let mut tau = vec![0.0f32; n];
    let mut itau = vec![0i32; n];
    let mut w = vec![0.0f32; nw * n];
    let mut buf = vec![0.0f32; n];

    for it in 0..ntmin { tau[it as usize] = -1.0; }
    for it in ntmin..nt {
        let ratio = it as f32 / ntmin as f32;
        tau[it as usize] = ((ratio as f64).ln() / dtau as f64) as f32;
    }
    lintrp(&tau, &mut w, &mut itau, ntau, n);

    loop {
        stretch(&mut buf, &tr.data, &w, &itau, n, nw);
        tr.data[..n].copy_from_slice(&buf);
        tr.ns = nt as u16;
        tr.dt = dt as u16;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
