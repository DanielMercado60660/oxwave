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
//! SUCHW - Change Header Word using one or two header word fields.
//! Port of su/main/headers/suchw.c (Revision 1.18).
//! Credits: SEP: Einar Kjartansson; CWP: Jack K. Cohen and John Stockwell;
//! Delphi: Alexander Koek.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUCHW - Change Header Word using one or two header word fields\t",
    "\t\t\t\t\t\t\t\t\t",
    "  suchw <stdin >stdout [optional parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " key1=cdp,...\toutput key(s) \t\t\t\t\t\t",
    " key2=cdp,...\tinput key(s) \t\t\t\t\t\t",
    " key3=cdp,...\tinput key(s)  \t\t\t\t\t\t",
    " a=0,...\t\toverall shift(s)\t\t\t\t",
    " b=1,...\t\tscale(s) on first input key(s) \t\t\t",
    " c=0,...\t\tscale on second input key(s) \t\t\t",
    " d=1,...\t\toverall scale(s)\t\t\t\t",
    " e=1,...\t\texponent on first input key(s)",
    " f=1,...\t\texponent on second input key(s)",
    "\t\t\t\t\t\t\t\t\t",
    " The value of header word key1 is computed from the values of\t\t",
    " key2 and key3 by:\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tval(key1) = (a + b * val(key2)^e + c * val(key3)^f) / d\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    " Shift cdp numbers by -1:\t\t\t\t\t\t",
    "\tsuchw <data >outdata a=-1\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Add 1000 to tracr value:\t\t\t\t\t\t",
    " \tsuchw key1=tracr key2=tracr a=1000 <infile >outfile\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " We set the receiver point (gx) field by summing the offset\t\t",
    " and shot point (sx) fields and then we set the cdp field by\t\t",
    " averaging the sx and gx fields (we choose to use the actual\t\t",
    " locations for the cdp fields instead of the conventional\t\t",
    " 1, 2, 3, ... enumeration):\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   suchw <indata key1=gx key2=offset key3=sx b=1 c=1 |\t\t\t",
    "   suchw key1=cdp key2=gx key3=sx b=1 c=1 d=2 >outdata\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Do both operations in one call:\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suchw<indata key1=gx,cdp key2=offset,gx key3=sx,sx b=1,1 c=1,1 d=1,2 >outdata",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let n = countparval("key1");
    let nkeys = if n != 0 { n as usize } else { 1 };
    let mut key1 = vec![String::new(); nkeys];
    if n != 0 { getparstringarray("key1", &mut key1); } else { key1[0] = "cdp".into(); }
    let mut key2 = vec![String::new(); nkeys];
    let n = countparval("key2");
    if n != 0 {
        if n as usize != nkeys { err!("number of key2's and key1's must be equal!"); }
        getparstringarray("key2", &mut key2);
    } else {
        if nkeys != 1 { err!("number of key2's and key1's must be equal!"); }
        key2[0] = "cdp".into();
    }
    let mut key3 = vec![String::new(); nkeys];
    let n = countparval("key3");
    if n != 0 {
        if n as usize != nkeys { err!("number of key3's and key1's must be equal!"); }
        getparstringarray("key3", &mut key3);
    } else {
        if nkeys != 1 { err!("number of key3's and key1's must be equal!"); }
        key3[0] = "cdp".into();
    }
    macro_rules! parameter {
        ($name:literal, $default:expr, $message:literal) => {{
            let n = countparval($name);
            let mut values = vec![$default; nkeys];
            if n != 0 {
                if n as usize != nkeys { err!($message); }
                getpardouble($name, &mut values);
            }
            values
        }};
    }
    let a = parameter!("a", 0.0, "number of a's and key1's must be equal!");
    let b = parameter!("b", 1.0, "number of b's and key1's must be equal!");
    let c = parameter!("c", 0.0, "number of c's and key1's must be equal!");
    let d = parameter!("d", 1.0, "number of d's and key1's must be equal!");
    let e = parameter!("e", 1.0, "number of e's and key1's must be equal!");
    let f = parameter!("f", 1.0, "number of f's and key1's must be equal!");
    checkpars();

    let mut type1 = Vec::new();
    let mut type2 = Vec::new();
    let mut type3 = Vec::new();
    let mut index1 = Vec::new();
    let mut index2 = Vec::new();
    let mut index3 = Vec::new();
    for ikey in 0..nkeys {
        type1.push(hdtype(&key1[ikey]));
        type2.push(hdtype(&key2[ikey]));
        type3.push(hdtype(&key3[ikey]));
        index1.push(getindex(&key1[ikey]));
        index2.push(getindex(&key2[ikey]));
        index3.push(getindex(&key3[ikey]));
    }
    let mut val1 = Value::default();
    let mut val2 = Value::default();
    let mut val3 = Value::default();
    while gettr(&mut tr) != 0 {
        for ikey in 0..nkeys {
            gethval(&tr, index2[ikey], &mut val2);
            gethval(&tr, index3[ikey], &mut val3);
            changeval(type1[ikey], &mut val1, type2[ikey], val2, type3[ikey], val3,
                a[ikey], b[ikey], c[ikey], d[ikey], e[ikey], f[ikey]);
            puthval(&mut tr, index1[ikey], &val1);
        }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

fn changeval(type1: &str, val1: &mut Value, type2: &str, val2: Value,
             type3: &str, val3: Value, a: f64, b: f64, c: f64,
             d: f64, e: f64, f: f64) {
    let dval2 = vtod(type2, val2);
    let dval3 = vtod(type3, val3);
    let dval1 = (a + b * dval2.powf(e) + c * dval3.powf(f)) / d;
    match type1.as_bytes()[0] {
        b's' => err!("can't change char header word"),
        b'h' => val1.set_h(dval1 as i32 as i16),
        b'u' => val1.set_u(dval1 as u32 as u16),
        b'l' => val1.set_l(dval1 as i64),
        b'v' => val1.set_v(dval1 as u64),
        b'i' => val1.set_i(dval1 as i32),
        b'p' => val1.set_p(dval1 as u32),
        b'f' => val1.set_f(dval1 as f32),
        b'd' => val1.set_d(dval1),
        _ => err!("unknown type %s", type1),
    }
}
