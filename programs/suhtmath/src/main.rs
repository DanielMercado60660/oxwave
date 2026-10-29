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

//! SUHTMATH - arithmetic on traces using header values.
//! Port of su/main/headers/suhtmath.c (Revision 1.4).
//! Credits: Matthias Imhof, Virginia Tech, Fri Dec 27 09:17:29 EST 2002.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUHTMATH - do unary arithmetic operation on segy traces with \t",
    "\t     headers values\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " suhtmath <stdin >stdout\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    "\tkey=tracl\theader word to use\t\t\t",
    "\top=nop\t\toperation flag\t\t\t\t",
    "\t\t\tnop   : no operation\t\t\t",
    "\t\t\tadd   : add header to trace\t\t",
    "\t\t\tmult  : multiply trace with header\t",
    "\t\t\tdiv   : divide trace by header\t\t",
    "\t\t\t\t\t\t\t\t",
    "\tscale=1.0\tscalar multiplier for header value\t",
    "\tconst=0.0\tadditive constant for header value\t",
    "\t\t\t\t\t\t\t\t",
    " Operation order:\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " op=add:\tout(t) = in(t) + (scale * key + const)\t\t",
    " op=mult:\tout(t) = in(t) * (scale * key + const)\t\t",
    " op=div:\tout(t) = in(t) / (scale * key + const)\t\t",
];

const ADD: i32 = 1;
const MULT: i32 = 2;
const DIV: i32 = 3;

fn getval(type_: &str, val: &Value) -> f32 {
    match type_.as_bytes()[0] {
        b's' => err!("can't get char header word"),
        b'h' => val.h() as f32,
        b'u' => val.u() as f32,
        b'l' => val.l() as f32,
        b'v' => val.v() as f32,
        b'i' => val.i() as f32,
        b'p' => val.p() as f32,
        b'f' => val.f(),
        b'd' => val.d() as f32,
        _ => err!("unknown type %s", type_),
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut c = 0.0f32;
    let mut s = 1.0f32;
    getparfloat("const", &mut c);
    getparfloat("scale", &mut s);
    let mut key = String::new();
    if getparstring("key", &mut key) == 0 { key = "tracl".to_string(); }
    let type_ = hdtype(&key);
    let indx = getindex(&key);
    let mut op = "nop".to_string();
    getparstring("op", &mut op);
    checkpars();

    let iop = if op == "add" { ADD } else if op == "mult" { MULT }
        else if op == "div" { DIV }
        else { err!("unknown operation=\"%s\", see self-doc", op); };

    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    loop {
        let mut val = Value::default();
        gethval(&tr, indx, &mut val);
        let f = s * getval(type_, &val) + c;
        match iop {
            ADD => for i in 0..nt { tr.data[i] += f; },
            MULT => for i in 0..nt { tr.data[i] *= f; },
            DIV => for i in 0..nt { tr.data[i] /= f; },
            _ => err!("mysterious operation=\"%s\"", op),
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
