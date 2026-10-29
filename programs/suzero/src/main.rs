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

//! SUZERO - zero-out (or set constant) data within a time window.
//! Port of su/main/amplitudes/suzero.c (Revision 1.14).
//! Credits: CWP: Chris; Geocon: Garry Perratt.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUZERO -- zero-out (or set constant) data within a time window\t",
    " \t\t\t\t\t\t\t\t",
    " suzero itmax= < indata > outdata\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters\t\t\t\t\t\t",
    " \titmax=\t\tlast time sample to zero out\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters\t\t\t\t\t\t",
    " \titmin=0\t\tfirst time sample to zero out\t\t",
    " \tvalue=0\t\tvalue to set\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " See also: sukill, sumute\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut itmin: i32 = 0; /* first sample to zero out */
    let mut itmax: i32 = 0; /* last sample to zero out */
    let mut value: f32 = 0.0; /* value to set within window */
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    mustgetparint("itmax", &mut itmax);
    if getparint("itmin", &mut itmin) == 0 { itmin = 0; }
    if getparfloat("value", &mut value) == 0 { value = 0.0; }
    checkpars();
    if itmax > nt { err!("itmax = %d, must be < nt", itmax); }
    if itmin < 0 { err!("itmin = %d, must not be negative", itmin); }
    if itmax < itmin { err!("itmax < itmin, not allowed"); }
    loop {
        for i in itmin..=itmax {
            // At maximum ns the permitted extra sample lies beyond the data
            // array and is never written to the trace.
            if (i as usize) < tr.data.len() { tr.data[i as usize] = value; }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
