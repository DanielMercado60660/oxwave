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

//! SUTSQ - time axis time-squared stretch of seismic traces.
//! Port of su/main/stretching_moveout_resamp/sutsq.c (Revision 1.8).

use oxwave_core::prelude::*;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUTSQ -- time axis time-squared stretch of seismic traces\t",
    "\t\t\t\t\t\t\t\t",
    " sutsq [optional parameters] <stdin >stdout \t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "\tnone\t\t\t\t \t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "       tmin= .1*nt*dt  minimum time sample of interest\t\t",
    "                       (only needed for forward transform)\t",
    "       dt= .004       output sample rate\t\t\t",
    "                       (only needed for inverse transform)\t",
    "       flag= 1        1=forward transform: time to time squared",
    "                     -1=inverse transform: time squared to time",
    "\t\t\t\t\t\t\t\t",
    " Note: The output of the forward transform always starts with\t",
    " time squared equal to zero.  'tmin' is used to avoid aliasing\t",
    " the early times.",
    " ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let mut dt = (tr.dt as f64 / 1000000.0) as f32;
    let mut tmin = 0.0f32;
    let mut flag = 0i32;
    let mut dtout = 0.0f32;
    if getparfloat("tmin", &mut tmin) == 0 {
        tmin = (0.1 * nt as f64 * dt as f64) as f32;
    }
    if getparint("flag", &mut flag) == 0 { flag = 1; }
    let ntout: i32;
    let mut ttn: Vec<f32>;
    if flag == 1 {
        dtout = (tmin as f64 * 2.0 * dt as f64) as f32;
        let tmax = nt as f32 * dt;
        ntout = (1.0f32 + tmax * tmax / dtout) as i32;
        CHECK_NT("ntout", ntout);
        if ntout < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
        ttn = vec![0.0; ntout as usize];
        for j in 0..ntout { ttn[j as usize] = ((j as f32 * dtout) as f64).sqrt() as f32; }
    } else {
        if getparfloat("dt", &mut dt) == 0 { dtout = 0.004; }
        ntout = (1.0 + ((nt as f32 * dt) as f64).sqrt() / dtout as f64) as i32;
        CHECK_NT("ntout", ntout);
        if ntout < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
        ttn = vec![0.0; ntout as usize];
        for j in 0..ntout { ttn[j as usize] = j.wrapping_mul(j) as f32 * dtout * dtout; }
    }
    checkpars();
    let mut buf = vec![0.0f32; nt as usize];
    fprintf!(stderr(), "sutsq: ntin=%d dtin=%f ntout=%d dtout=%f\n", nt, dt, ntout, dtout);
    loop {
        for j in 0..nt as usize { buf[j] = tr.data[j]; }
        tr.ns = ntout as u16;
        tr.dt = (dtout as f64 * 1000000.0) as i64 as u16;
        ints8r(dt, 0.0, &buf, 0.0, 0.0, &ttn, &mut tr.data[..ntout as usize]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
