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
//! SUGAUSSTAPER - Multiply traces with gaussian taper.
//! Port of su/main/tapering/sugausstaper.c (Revision 1.5).
//!
//! Credits (from the C source): Thomas Bohlen, formerly of TU Bergakademie,
//! Freiberg GDR, most recently of U Karlsruhe, 04.01.2002.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUGAUSSTAPER - Multiply traces with gaussian taper\t\t",
    " \t\t\t\t\t\t\t\t",
    " sugausstaper < stdin > stdout [optional parameters]\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t   \t",
    "   <none>\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t   \t",
    " key=offset    keyword of header field to weight traces by \t",
    " x0=300        key value defining the center of gaussian window",
    " xw=50         width of gaussian window in units of key value \t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " Traces are multiplied with a symmetrical gaussian taper \t",
    "  \tw(t)=exp(-((key-x0)/xw)**2)\t\t\t\t",
    " unlike \"sutaper\" the value of x0 defines center of the taper",
    " rather than the edges of the data.\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut key = String::new();
    let mut x0: f32 = 0.0;
    let mut xw: f32 = 0.0;

    initargs(&argv);
    requestdoc(1, SDOC);

    if getparstring("key", &mut key) == 0 { key = "offset".to_string(); }
    if getparfloat("x0", &mut x0) == 0 { x0 = 300.0; }
    if getparfloat("xw", &mut xw) == 0 { xw = 50.0; }
    checkpars();

    let type_ = hdtype(&key);
    let index = getindex(&key);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let ns = tr.ns as usize;

    loop {
        let mut val = Value::default();
        gethval(&tr, index, &mut val);
        let fval = vtof(type_, val);
        for i in 0..ns {
            let ratio = (fval - x0) / xw;
            tr.data[i] = ((tr.data[i] as f64) * (-(ratio as f64).powf(2.0)).exp()) as f32;
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
