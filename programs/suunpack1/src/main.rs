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
//! SUUNPACK1 - unpack segy trace data from chars to floats.
//! Port of su/main/data_compression/suunpack1.c (Revision 1.19).
//! Credits: CWP: Jack K. Cohen, Shuki Ronen, Brian Sumner.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, closeto};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    "SUUNPACK1 - unpack segy trace data from chars to floats\t",
    "\t\t\t\t\t\t\t\t",
    "    suunpack1 <packed_file >unpacked_file\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "suunpack1 is the approximate inverse of supack1\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    if tr.trid != CHARPACK { err!("Not char packed traces"); }
    let nt = tr.ns as usize;
    let ungpow = tr.ungpow;
    let isone = closeto(ungpow, 1.0);
    let istwo = closeto(ungpow, 2.0);

    loop {
        // Expand from the back because the input bytes share the float array.
        for i in (0..nt).rev() {
            let mut val = tr.data_bytes()[i] as i8 as f32;
            val *= tr.unscale;
            tr.data[i] = if istwo {
                val * abs(val)
            } else if isone {
                val
            } else if val >= 0.0 {
                (val as f64).powf(ungpow as f64) as f32
            } else {
                (-((-val) as f64).powf(ungpow as f64)) as f32
            };
        }
        tr.trid = 1;
        tr.ungpow = 0.0;
        tr.unscale = 0.0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
