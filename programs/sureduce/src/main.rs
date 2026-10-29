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
//! SUREDUCE - convert traces to display in reduced time.
//! Port of su/main/stretching_moveout_resamp/sureduce.c (Revision 1.13).
//! Credits (from the C source): UC Davis: Mike Begnaud, March 1995.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t   \t\t\t\t\t\t\t",
    " SUREDUCE - convert traces to display in reduced time\t\t",
    " \t   \t\t\t\t\t\t\t",
    " sureduce <stdin >stdout rv=\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "\tdt=tr.dt\tif not set in header, dt is mandatory\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "\trv=8.0\t\treducing velocity in km/sec\t\t",
    "\t\t\t\t\t\t\t\t",
    " Note: Useful for plotting refraction seismic data. \t\t",
    " To remove reduction, do:\t\t\t\t\t",
    " suflip < reduceddata.su flip=3 | sureduce rv=RV > flip.su\t",
    " suflip < flip.su flip=3 > unreduceddata.su\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Trace header fields accessed: dt, ns, offset\t\t\t",
    " Trace header fields modified: none\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut rv: f32 = 0.0;
    if getparfloat("rv", &mut rv) == 0 { rv = 8.0; }
    if rv <= 0.0 { err!("rv=%f, must be >= 0.0", rv); }
    if fgettr(&stdin(), &mut tr) == 0 { err!("can't read first trace"); }
    if tr.dt == 0 {
        let mut dt: f32 = 0.0;
        mustgetparfloat("dt", &mut dt);
    }
    let nt = tr.ns as i32;
    let dt = ((tr.dt as f64) / 1000000.0) as f32;
    checkpars();
    loop {
        let off = tr.offset.wrapping_abs() as f32;
        let bt = ((off as f64).abs() / ((rv as f64) * 1000.0)) as f32;
        let rnt = nint((bt / dt) as f64);
        for i in 0..nt {
            let j = i.wrapping_add(rnt);
            tr.data[i as usize] = if i < nt.wrapping_sub(rnt) {
                tr.data[j as usize]
            } else { 0.0 };
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
