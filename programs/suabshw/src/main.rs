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


//! SUABSHW - replace a header key word by its absolute value.
//! Port of su/main/headers/suabshw.c (Revision 1.13).
//! Credits: CWP: Jack K. Cohen

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUABSHW - replace header key word by its absolute value\t",
    " \t\t\t\t\t\t\t\t",
    " suabshw <stdin >stdout key=offset\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    " \tkey=offset\t\theader key word\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut key = String::new();
    let mut val = Value::default();
    initargs(&argv);
    requestdoc(1, SDOC);
    if getparstring("key", &mut key) == 0 {
        key = "offset".to_string();
    }
    checkpars();
    let type_ = hdtype(&key);
    let index = getindex(&key);
    while gettr(&mut tr) != 0 {
        gethval(&tr, index, &mut val);
        absval(type_, &mut val);
        puthval(&mut tr, index, &val);
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

fn absval(type_: &str, val: &mut Value) {
    match type_.as_bytes().first().copied() {
        Some(b's') => err!("can't absval char header word"),
        Some(b'u' | b'v' | b'p') => {},
        Some(b'h') => if val.h() < 0 { val.set_h(val.h().wrapping_neg()); },
        Some(b'l') => if val.l() < 0 { val.set_l(val.l().wrapping_neg()); },
        Some(b'i') => if val.i() < 0 { val.set_i(val.i().wrapping_neg()); },
        Some(b'f') => if val.f() < 0.0 { val.set_f(-val.f()); },
        Some(b'd') => if val.d() < 0.0 { val.set_d(-val.d()); },
        _ => err!("unknown type %s", type_),
    }
}
