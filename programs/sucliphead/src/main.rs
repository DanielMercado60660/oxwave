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

//! SUCLIPHEAD - Clip header values. Port of su/main/headers/sucliphead.c.
//! Credits: Geocon: Garry Perratt.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUCLIPHEAD - Clip header values\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sucliphead <stdin >stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tkey=cdp,...\t\t\theader key word(s) to clip\t",
    "\tmin=0,...\t\t\tminimum value to clip\t\t",
    "\tmax=ULONG_MAX,ULONG_MAX,...\tmaximum value to clip\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn changeval(type_: &str, val: &mut Value, fval: f32) {
    match type_.as_bytes()[0] {
        b's' => err!("cannot change char header word"),
        b'h' => val.set_h((fval as i32) as i16),
        b'u' => val.set_u((fval as u32) as u16),
        b'l' => val.set_l(fval as i64),
        b'v' => val.set_v(fval as u64),
        b'i' => val.set_i(fval as i32),
        b'p' => val.set_p(fval as u32),
        b'f' => val.set_f(fval),
        b'd' => val.set_d(fval as f64),
        _ => err!("unknown type %s", type_),
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let nkeys = countparval("key") as usize;
    let mut key = vec![String::new(); SU_NKEYS];
    if nkeys != 0 {
        getparstringarray("key", &mut key);
    } else {
        key[0] = "cdp".to_string();
    }
    let mut types = Vec::with_capacity(nkeys);
    let mut index = Vec::with_capacity(nkeys);
    for k in 0..nkeys {
        types.push(hdtype(&key[k]));
        index.push(getindex(&key[k]));
    }
    let n = countparval("min") as usize;
    let mut min = vec![0.0f32; nkeys];
    if n != 0 {
        if n != nkeys { err!("number of a values not equal to number of keys"); }
        getparfloat("min", &mut min);
    }
    let n = countparval("max") as usize;
    let mut max = vec![u64::MAX as f32; nkeys];
    if n != 0 {
        if n != nkeys { err!("number of a values not equal to number of keys"); }
        getparfloat("max", &mut max);
    }
    checkpars();
    for k in 0..nkeys {
        types[k] = hdtype(&key[k]);
        index[k] = getindex(&key[k]);
    }
    while gettr(&mut tr) != 0 {
        for k in 0..nkeys {
            let mut val = Value::default();
            gethval(&tr, index[k], &mut val);
            let fval = vtof(types[k], val) as f64;
            if fval < min[k] as f64 {
                changeval(types[k], &mut val, min[k]);
                puthval(&mut tr, index[k], &val);
            } else if fval > max[k] as f64 {
                changeval(types[k], &mut val, max[k]);
                puthval(&mut tr, index[k], &val);
            }
        }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
