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
//! SUADDEVENT - add a linear or hyperbolic moveout event to seismic data.
//! Port of su/main/synthetics_waveforms_testpatterns/suaddevent.c (Revision 1.6).
//! Credits: Gary Billings, Talisman Energy, May 1996, Apr 2000, June 2001.

use oxwave_core::prelude::*;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t       ",
    " SUADDEVENT - add a linear or hyperbolic moveout event to seismic data ",
    "\t\t\t\t\t\t\t\t       ",
    " suaddevent <stdin >stdout [optional parameters]\t\t       ",
    "\t\t\t\t\t\t\t\t       ",
    " Required parameters:\t\t\t\t\t\t  ",
    "       none\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t       ",
    " Optional parameters:\t\t\t\t\t\t  ",
    "     type=nmo    =lmo for linear event \t\t\t\t",
    "     t0=1.0      zero-offset intercept time IN SECONDS\t\t\t",
    "     vel=3000.   moveout velocity in m/s\t\t\t\t",
    "     amp=1.      amplitude\t\t\t\t\t\t",
    "     dt=\t must provide if 0 in headers (seconds)\t\t",
    "   \t\t\t\t\t\t\t\t\t",
    " Typical usage: ",
    "     sunull nt=500 dt=0.004 ntr=100 | sushw key=offset a=-1000 b=20 \\ ",
    "     | suaddevent v=1000 t0=0.05 type=lmo | suaddevent v=1800 t0=0.8 \\",
    "     | sufilter f=8,12,75,90 | suxwigb clip=1 &\t     \t\t",
    "\t\t\t\t\t\t\t\t       ",
    "\t\t\t\t\t\t\t\t       ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut nsegy = gettr(&mut tr);
    if nsegy == 0 { err!("cannot get first trace"); }
    let nt = tr.ns as usize;
    let mut dt = (tr.dt as f64 / 1000000.0) as f32;
    let tmin = (tr.delrt as f64 / 1000.0) as f32;
    let mut addtr = vec![0.0f32; nt];
    let mut times = vec![0.0f32; nt];
    for i in 0..nt { times[i] = tmin + i as f32 * dt; addtr[i] = 0.0; }
    if dt == 0.0 { mustgetparfloat("dt", &mut dt); }
    let mut type_ = String::new();
    let mut typecode = 0;
    if getparstring("type", &mut type_) == 0 { type_ = "nmo".to_string(); }
    if type_ == "lmo" { typecode = 1; }
    let mut t0 = 1.0f32;
    let mut vel = 3000.0f32;
    let mut amp = 1.0f32;
    getparfloat("t0", &mut t0);
    getparfloat("vel", &mut vel);
    getparfloat("amp", &mut amp);
    while nsegy != 0 {
        let tx = if typecode == 1 {
            t0 + tr.offset.wrapping_abs() as f32 / vel
        } else {
            let x = tr.offset as f32 / vel;
            ((t0 * t0 + x * x) as f64).sqrt() as f32
        };
        ints8r(dt, tx, &[amp], 0.0, 0.0, &times, &mut addtr);
        for i in 0..nt { tr.data[i] += addtr[i]; addtr[i] = 0.0; }
        puttr(&mut tr);
        nsegy = gettr(&mut tr);
    }
    exit(CWP_Exit());
}
