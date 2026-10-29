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

// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! SUSEHW - set the trace number within an ensemble.
//! Port of su/main/headers/susehw.c (Revision 1.4).
//! Credits: CWP: John Stockwell (Feb 2008), in answer to Warren Franz,
//! based on codes including susplit, susshw and suchw. The header-setting
//! routine was hacked from sushw by John Stockwell.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSEHW - Set the value the Header Word denoting trace number within\t",
    "\t     an Ensemble defined by the value of another header word\t",
    "\t\t\t\t\t\t\t\t\t",
    "     susehw <stdin >stdout [options]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " key1=cdp\tKey header word defining the ensemble\t\t\t",
    " key2=cdpt\tKey header word defining the count within the ensemble\t",
    " a=1\t\tstarting value of the count in the ensemble\t\t",
    " b=1\t\tincrement or decrement within the ensemble\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This code was written because suresstat requires cdpt to be set.\t",
    " The computation is \t\t\t\t\t\t\t",
    " \tval(key2) = a + b*i\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The input data must first be sorted into constant key1 gathers.\t",
    " Example: setting the cdpt field\t\t\t\t\t",
    "        susehw < cdpgathers.su a=1 b=1 key1=cdp key2=cdpt > new.su\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key1 = String::new();
    let mut key2 = String::new();
    let mut a = 1.0f64;
    let mut b = 1.0f64;
    let mut verbose = 0i32;
    if getparstring("key1", &mut key1) == 0 { key1 = "cdp".to_string(); }
    if getparstring("key2", &mut key2) == 0 { key2 = "cdpt".to_string(); }
    if getpardouble("a", &mut a) == 0 { a = 1.0; }
    if getpardouble("b", &mut b) == 0 { b = 1.0; }
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    checkpars();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let type1 = hdtype(&key1);
    let index1 = getindex(&key1);
    let type2 = hdtype(&key2);
    let index2 = getindex(&key2);
    if verbose != 0 { warn!("key1 = %s", key1); }
    if verbose != 0 { warn!("key2 = %s", key2); }
    let mut val1 = Value::default();
    let mut val2 = Value::default();
    gethval(&tr, index1, &mut val1);
    let mut ival1old = vtoi(type1, val1);
    let mut i = 0i32;
    loop {
        gethval(&tr, index1, &mut val1);
        let ival1 = vtoi(type1, val1);
        if ival1 != ival1old { ival1old = ival1; i = 0; }
        setval(type2, &mut val2, a, b, i as f64);
        puthval(&mut tr, index2, &val2);
        i = i.wrapping_add(1);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}

// The original helper sets the union member corresponding to the header type.
fn setval(type_: &str, valp: &mut Value, a: f64, b: f64, i: f64) {
    let value = a + b * i;
    match type_.as_bytes()[0] {
        b's' => err!("can't set char header word"),
        b'h' => valp.set_h(value as i32 as i16),
        b'u' => valp.set_u(value as u32 as u16),
        b'l' => valp.set_l(value as i64),
        b'v' => valp.set_v(value as u64),
        b'i' => valp.set_i(value as i32),
        b'p' => valp.set_p(value as u32),
        b'f' => valp.set_f(value as f32),
        b'd' => { valp.set_d(value); err!("unknown type %s", type_); }
        _ => err!("unknown type %s", type_),
    }
}
