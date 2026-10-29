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

//! SURAMP - linearly taper the start and/or end of traces to zero.
//! Port of su/main/tapering/suramp.c (Revision 1.13).
//! Credits: CWP: Jack K. Cohen, Ken Larner.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{max, nint};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SURAMP - Linearly taper the start and/or end of traces to zero.\t",
    " \t\t\t\t\t\t\t\t\t",
    " suramp <stdin >stdout [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " \tif dt is not set in header, then dt is mandatory\t\t",
    " \t\t\t\t\t\t\t       \t\t",
    " Optional parameters\t\t\t\t\t\t\t",
    "\ttmin=tr.delrt/1000\tend of starting ramp (sec)\t\t",
    "\ttmax=(nt-1)*dt\t\tbeginning of ending ramp (sec)\t\t",
    " \tdt = (from header)\tsampling interval (sec)\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The taper is a linear ramp from 0 to tmin and/or tmax to the\t\t",
    " end of the trace.  Default is a no-op!\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if dt == 0.0 { err!("dt field is zero and not getparred"); }
    let mut tmin = 0.0f32;
    let mut tmax = 0.0f32;
    if getparfloat("tmin", &mut tmin) == 0 { tmin = (tr.delrt as f64 / 1000.0) as f32; }
    if getparfloat("tmax", &mut tmax) == 0 { tmax = (nt - 1) as f32 * dt; }
    checkpars();
    let ntaper1 = nint(max(0.0, (tmin as f64 - tr.delrt as f64 / 1000.0) / dt as f64));
    let ntaper2 = nint(max(0.0, (tr.delrt as f64 / 1000.0 + ((nt - 1) as f32 * dt) as f64 - tmax as f64) / dt as f64));
    // C has undefined behavior outside tr.data. Check before allocating weights
    // or converting the ending ramp's potentially negative index to usize.
    // Starting ramps may extend past ns while still fitting in the C buffer.
    if ntaper1 as usize > tr.data.len() {
        err!("starting taper length %d exceeds sample buffer size %d", ntaper1, tr.data.len() as i32);
    }
    if ntaper2 > nt {
        err!("ending taper length %d exceeds trace length %d", ntaper2, nt);
    }
    let mut taper1 = vec![0.0f32; ntaper1 as usize];
    let mut taper2 = vec![0.0f32; ntaper2 as usize];
    for i in 0..ntaper1 { taper1[i as usize] = (i + 1) as f32 / ntaper1 as f32; }
    for i in 0..ntaper2 { taper2[i as usize] = (ntaper2 - i) as f32 / ntaper2 as f32; }
    loop {
        for i in 0..ntaper1 { tr.data[i as usize] *= taper1[i as usize]; }
        for i in 0..ntaper2 { tr.data[(nt - ntaper2 + i) as usize] *= taper2[i as usize]; }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
