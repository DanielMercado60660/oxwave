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
 
//! SUKILL - zero out traces. Port of su/main/windowing_sorting_muting/sukill.c (Revision 1.18).
//! Credits (from the C source): CWP: Chris Liner, Jack K. Cohen;
//! header-based trace selection: Florian Bleibinhaus.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUKILL - zero out traces\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " sukill <stdin >stdout [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "\tkey=trid\theader name to select traces to kill\t",
    "\ta=2\t\theader value identifying tracces to kill",
    " or\t\t\t\t\t\t\t\t",
    " \tmin= \t\tfirst trace to kill (one-based)\t\t",
    " \tcount=1\t\tnumber of traces to kill \t\t",
    " \t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    "\tIf min= is set it overrides selecting traces by header.\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut min = 0i32;
    let mut count = 1i32;
    let mut key = String::new();
    let mut a = 2.0f64;

    initargs(&argv);
    requestdoc(1, SDOC);
    if getparint("min", &mut min) == 0 { min = 0; }
    if min < 0 { err!("min = %d, must be >= 1", min); }
    if getparint("count", &mut count) == 0 { count = 1; }
    if getparstring("key", &mut key) == 0 { key = "trid".to_string(); }
    if getpardouble("a", &mut a) == 0 { a = 2.0; }
    checkpars();

    let type_ = hdtype(&key);
    let index = getindex(&key);
    if min > 0 {
        let mut itr = 1i32;
        while itr < min {
            if gettr(&mut tr) != 0 { puttr(&mut tr); }
            else { err!("failed to get requested trace #%ld", itr as i64); }
            itr = itr.wrapping_add(1);
        }
        let nt;
        if gettr(&mut tr) != 0 {
            nt = tr.ns as usize;
            tr.data[..nt].fill(0.0);
            puttr(&mut tr);
            itr = itr.wrapping_add(1);
        } else { err!("failed to get requested trace #%ld", itr as i64); }
        while itr < min.wrapping_add(count) {
            if gettr(&mut tr) != 0 {
                tr.data[..nt].fill(0.0);
                puttr(&mut tr);
            } else { err!("failed to get requested trace #%ld", itr as i64); }
            itr = itr.wrapping_add(1);
        }
        while gettr(&mut tr) != 0 { puttr(&mut tr); }
    } else {
        while gettr(&mut tr) != 0 {
            let nt = tr.ns as usize;
            let mut val = Value::default();
            gethval(&tr, index, &mut val);
            let dval = vtod(type_, val);
            if dval == a { tr.data[..nt].fill(0.0); }
            puttr(&mut tr);
        }
    }
    exit(CWP_Exit());
}
