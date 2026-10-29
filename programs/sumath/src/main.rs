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

//! SUMATH - do math operation on SU data.
//! Port of su/main/operations/sumath.c (Revision 1.1).
//! Credits (from the C source): U Arkansas: Chris Liner Jun 2013.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, sgn};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUMATH - do math operation on su data \t\t",
    " \t\t\t\t\t\t\t\t",
    " suop <stdin >stdout op=mult\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    "\top=mult\t\toperation flag\t\t\t\t",
    "\t\t\t--------- operations -----------------\t",
    "\t\t\tadd   : o = i + a    (o=out; i=in)\t",
    "\t\t\tsub   : o = i - a\t\t\t",
    "\t\t\tmult  : o = i * a  \t\t\t",
    "\t\t\tdiv   : o = i / a  \t\t\t",
    "\t\t\tpow   : o = i ^ a\t\t\t",
    "\t\t\tspow  : o = sgn(i) * abs(i) ^ a  \t",
    "\t\t\t--------- operation parameter --------\t",
    "\ta=1\t\t\t\t\t\t\t",
    "\tcopy=1\t\tn>1 copy each trace n times\t\t",
    "\t\t\t\t\t\t\t\t",
    " Note:\t\t\t\t\t\t\t\t",
    " There is overlap between this program and \"sugain\" and\t",
    " \"suop\"\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut tracl = tr.tracl;
    let mut op = String::from("mult");
    getparstring("op", &mut op);
    let iop = match op.as_str() {
        "add" => 1, "sub" => 2, "mult" => 3, "div" => 4,
        "pow" => 5, "spow" => 6,
        _ => err!("unknown operation=\"%s\", see self-doc", op),
    };
    let mut a = 1.0f32;
    if getparfloat("a", &mut a) == 0 { a = 1.0; }
    let mut copy = 1i32;
    if getparint("copy", &mut copy) == 0 { copy = 1; }
    if copy > 1 { tracl = 1; }
    loop {
        for i in 0..nt {
            match iop {
                1 => tr.data[i] += a,
                2 => tr.data[i] -= a,
                3 => tr.data[i] *= a,
                4 => tr.data[i] /= a,
                5 => tr.data[i] = (tr.data[i] as f64).powf(a as f64) as f32,
                6 => tr.data[i] = (sgn(tr.data[i]) * (abs(tr.data[i]) as f64).powf(a as f64)) as f32,
                _ => err!("mysterious operation=\"%s\"", op),
            }
        }
        if copy == 1 { puttr(&mut tr); }
        else {
            for _j in 1..=copy {
                tr.tracl = tracl;
                puttr(&mut tr);
                tracl = tracl.wrapping_add(1);
            }
        }
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
