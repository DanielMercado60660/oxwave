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
//! SUPACK1 - pack segy trace data into chars.
//! Port of su/main/data_compression/supack1.c (Revision 1.20).
//! Credits: CWP: Jack K. Cohen, Shuki Ronen, Brian Sumner.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, closeto, max};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUPACK1 - pack segy trace data into chars\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " supack1 <segy_file >packed_file\tgpow=0.5 \t\t",
    "\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t        \t",
    " Optional parameter: \t\t\t\t\t\t",
    "\tgpow=0.5\texponent used to compress the dynamic\t",
    "\t\t\trange of the traces\t\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut gpow: f32 = 0.0;
    if getparfloat("gpow", &mut gpow) == 0 { gpow = 0.5; }
    checkpars();
    if gpow <= 0.0 { err!("gpow = %g must be positive", gpow); }
    let isone = closeto(gpow, 1.0);
    let ishalf = closeto(gpow, 0.5);

    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    loop {
        if !isone {
            for i in 0..nt {
                let val = tr.data[i];
                tr.data[i] = if ishalf {
                    if val >= 0.0 { (val as f64).sqrt() as f32 }
                    else { -((-val as f64).sqrt()) as f32 }
                } else if val >= 0.0 {
                    (val as f64).powf(gpow as f64) as f32
                } else {
                    -((-val as f64).powf(gpow as f64)) as f32
                };
            }
        }
        tr.ungpow = (1.0 / gpow as f64) as f32;
        let mut absmax = abs(tr.data[0]);
        for i in 1..nt { absmax = max(absmax, abs(tr.data[i])); }
        let scale = if absmax != 0.0 { 127.0f32 / absmax } else { 0.0 };
        tr.unscale = if absmax != 0.0 { (1.0 / scale as f64) as f32 } else { 0.0 };
        for i in 0..nt {
            tr.data[i] *= scale;
            let packed = tr.data[i] as i32 as i8 as u8;
            tr.data_bytes_mut()[i] = packed;
        }
        tr.trid = CHARPACK;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
