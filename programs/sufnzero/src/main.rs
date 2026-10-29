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
//! SUFNZERO - get time of first non-zero sample by trace.
//! Port of su/main/picking/sufnzero.c (Revision 1.5).
//! Credits: Geocon: Garry Perratt; based on surms by the same,
//! itself based on sugain and sumax by CWP: John Stockwell.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUFNZERO - get Time of First Non-ZERO sample by trace              ",
    "\t\t\t\t\t\t\t\t\t",
    "  sufnzero <stdin >stdout [optional parameters] \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters: \t\t\t\t\t\t\t",
    "\tmode=first   \tOutput time of first non-zero sample\t\t",
    "\t             \t=last for last non-zero sample\t\t\t",
    "\t             \t=both for both first & last non-zero samples    ",
    "\t\t\t\t\t\t\t\t\t",
    "\tmin=0   \tThreshold value for considering as zero         ",
    "\t\t\tAny abs(sample)<min is considered to be zero\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tkey=key1,...\tKeyword(s) to print\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut key = vec![String::new(); SU_NKEYS];
    let mut mode = String::new();
    let mut min = 0.0f32;
    initargs(&argv);
    requestdoc(1, SDOC);
    if getparstring("mode", &mut mode) == 0 { mode = "first".to_string(); }
    if mode != "first" && mode != "last" && mode != "both" {
        err!("%s unknown mode", mode);
    }
    if getparfloat("min", &mut min) == 0 { min = 0.0; }
    let nkeys = countparval("key");
    if nkeys != 0 { getparstringarray("key", &mut key); }
    checkpars();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    loop {
        let mut ftime = -1.0f32;
        let mut ltime = -1.0f32;
        for i in 0..nt {
            if (tr.data[i] as f64).abs() > min as f64 && ftime == -1.0 {
                ftime = ((i as i32).wrapping_mul(tr.dt as i32) / 1000) as f32;
            }
            if (tr.data[i] as f64).abs() > min as f64 {
                ltime = ((i as i32).wrapping_mul(tr.dt as i32) / 1000) as f32;
            }
        }
        for ikey in 0..nkeys as usize {
            let index = getindex(&key[ikey]);
            let mut keyval = Value::default();
            gethval(&tr, index, &mut keyval);
            printfval(hdtype(&key[ikey]), keyval);
            printf!(" ");
        }
        if mode == "first" { printf!("%g\n", ftime); }
        else if mode == "last" { printf!("%g\n", ltime); }
        else if mode == "both" { printf!("%g %g\n", ftime, ltime); }
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
