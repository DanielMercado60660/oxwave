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
// Copyright (c) Colorado School of Mines, 2010. All rights reserved.
//! SUAHW - Assign Header Word using another header word.
//! Port of su/main/headers/suahw.c (Revision 1.3).
//! Credits: Florian Bleibinhaus, U Salzburg, Austria; cloned from suchw of
//! Einar Kajartansson, SEP.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUAHW - Assign Header Word using another header word\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  suahw <stdin >stdout [optional parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "  key1=ep\t\toutput key \t\t\t\t\t",
    "  key2=fldr\t\tinput key \t\t\t\t\t",
    "  a=\t\t\tarray of key1 output values\t\t\t",
    "  b=\t\t\tarray of key2 input values\t\t\t",
    "  mode=extrapolate\thow to assign a key1-value, when the key2-value\t",
    "\t\t\tis not found in b:\t\t\t\t",
    "\t\t\t=interpolate\tinterpolate\t\t\t",
    "\t\t\t=extrapolate\tinterpolate and extrapolate\t",
    "\t\t\t=zero\t\tzero key1-values\t\t",
    "\t\t\t=preserve\tpreserve key1-values\t\t",
    "\t\t\t=transfer\ttransfer key2-values to key1\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "  key3=tracf\t\tinput key \t\t\t\t\t",
    "  c=\t\t\tarray of key3 input values\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The key1-value is assigned based on the key2-value and the arrays a,b.",
    " If the header value of key2= equals the n'th element in b=, then the\t",
    " header value key1= is set to the n'th element in a=.\t\t\t",
    " The arrays a= and b= must have the same size, and the elements of b=\t",
    " must be in ascending order.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The mode-switch decides what to do when a trace header has a key2-value",
    " that is not an element of the b-array:\t\t\t\t",
    "    zero - the key1-value will be set to zero\t\t\t\t",
    "    preserve - the key1-value will not be modified\t\t\t",
    "    transfer - the key2-value will be assigned to key1\t\t\t",
    "    interpolate - if the key2-value is greater than the n'th element\t",
    "\tand less than the (n+1)'th element of b=, then the key1-value\t",
    "\twill be\tinterpolated accordingly from the n'th and (n+1)'th\t",
    "\telement of a=. Otherwise, key1 will not be changed.\t\t",
    "    extrapolate - same as interpolate, plus, if the key2-value is\t",
    "\tsmaller/greater than the first/last element of b=, then the\t",
    "\tkey1-value will be set to the first/last element of a=\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The array c= can be used to prevent the modification of trace headers\t",
    " with certain key3-values. The number of elements in c= is independent\t",
    " of the other arrays.\t\t\t\t\t\t\t",
    " The key1-value will not be modified, if the mode-switch is set to\t",
    "    zero, preserve, transfer - and the key3-value is an element of c=\t",
    "    interpolate, extrapolate - and the key3-value is outside of c=\t",
    "\t\t\t\t(smaller than the first or greater than\t",
    "\t\t\t\tthe last element of c=)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "  Assign shot numbers 1-3 to field file ID 1009,1011,1015 and 0 to the\t",
    "  remaining FFID (fldr):\t\t\t\t\t\t",
    "    suahw <data a=1,2,3 b=1009,1011,1015 mode=zero\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  Use channel numbers (tracf) to assign stations numbers (tracr) for a\t",
    "  split spread with a gap:\t\t\t\t\t\t",
    "    suahw <data key1=tracr a=151,128,124,101 key2=tracf b=1,24,25,48\t",
    "\t\t\t\t\t\t\t\t\t",
    "  Assign shot-statics:\t\t\t\t\t\t\t",
    "    suahw <data key1=sstat key2=ep a=-32,13,-4 b=1,2,3\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  Set trid to 0 for channel 1-24, but only for the record 1016:\t",
    "    suahw <data key1=trid key2=tracf key3=fldr a=0,0 b=1,24 c=1016\t",
];

const ZERO: i32 = 0;
const INTERPOLATE: i32 = 1;
const EXTRAPOLATE: i32 = 2;
const PRESERVE: i32 = 3;
const TRANSFER: i32 = 4;

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut mode = String::new();
    if getparstring("mode", &mut mode) == 0 { mode = "extrapolate".into(); }
    let imode = match mode.as_str() {
        "zero" => ZERO, "interpolate" => INTERPOLATE, "extrapolate" => EXTRAPOLATE,
        "preserve" => PRESERVE, "transfer" => TRANSFER,
        _ => err!("unknown mode=\"%s\", see self-doc", mode),
    };
    let mut key1 = String::new();
    let mut key2 = String::new();
    let mut key3 = String::new();
    if getparstring("key1", &mut key1) == 0 { key1 = "ep".into(); }
    if getparstring("key2", &mut key2) == 0 { key2 = "fldr".into(); }
    if getparstring("key3", &mut key3) == 0 { key3 = "tracf".into(); }
    let nvals = countparval("a") as usize;
    if nvals == 0 { err!("need array of a-values!"); }
    let mut a = vec![0.0f64; nvals];
    getpardouble("a", &mut a);
    if countparval("b") != nvals as i32 { err!("number of b's and a's must be equal!"); }
    let mut b = vec![0.0f64; nvals];
    getpardouble("b", &mut b);
    for i in 1..nvals { if b[i-1] >= b[i] { err!("b's must be in ascending order!"); } }
    let nvalsc = countparval("c") as usize;
    let mut c = vec![0.0f64; nvalsc];
    if nvalsc > 0 { getpardouble("c", &mut c); }
    checkpars();
    let type1 = hdtype(&key1);
    let type2 = hdtype(&key2);
    let type3 = hdtype(&key3);
    let index1 = getindex(&key1);
    let index2 = getindex(&key2);
    let index3 = getindex(&key3);
    while gettr(&mut tr) != 0 {
        let mut val1 = Value::default();
        let mut val2 = Value::default();
        let mut val3 = Value::default();
        gethval(&tr, index1, &mut val1);
        gethval(&tr, index2, &mut val2);
        gethval(&tr, index3, &mut val3);
        let mut dval1 = vtod(type1, val1);
        let dval2 = vtod(type2, val2);
        let dval3 = vtod(type3, val3);
        let mut modval1 = true;
        if nvalsc > 0 {
            match imode {
                ZERO | PRESERVE | TRANSFER => {
                    modval1 = false;
                    for i in 0..nvalsc { if dval3 == c[i] { modval1 = true; } }
                }
                INTERPOLATE | EXTRAPOLATE => {
                    if dval3 < c[0] { modval1 = false; }
                    if dval3 > c[nvalsc-1] { modval1 = false; }
                }
                _ => {}
            }
        }
        let mut ival = 0;
        for i in 0..nvals { if dval2 >= b[i] { ival = i; } }
        if dval2 == b[ival] { dval1 = a[ival]; }
        else {
            match imode {
                ZERO => dval1 = 0.0,
                INTERPOLATE => {
                    if dval2 >= b[0] && dval2 <= b[nvals-1] {
                        if ival + 1 >= nvals {
                            err!("cannot interpolate key2=\"%s\": no upper table entry", key2);
                        }
                        dval1 = a[ival] + (a[ival+1]-a[ival])*(dval2-b[ival])/(b[ival+1]-b[ival]);
                    }
                }
                EXTRAPOLATE => {
                    if dval2 < b[0] { dval1 = a[0]; }
                    else if dval2 > b[nvals-1] { dval1 = a[nvals-1]; }
                    else {
                        // A NaN header or table entry can bypass both endpoint
                        // comparisons. Reject only an out-of-bounds interval: NaN
                        // with two available entries is defined floating-point math.
                        // The C access past the table is undefined; fail cleanly.
                        if ival + 1 >= nvals {
                            err!("cannot interpolate key2=\"%s\": no upper table entry", key2);
                        }
                        dval1 = a[ival] + (a[ival+1]-a[ival])*(dval2-b[ival])/(b[ival+1]-b[ival]);
                    }
                }
                TRANSFER => dval1 = dval2,
                _ => {}
            }
        }
        if modval1 { assgnval(type1, &mut val1, dval1); puthval(&mut tr, index1, &val1); }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

fn assgnval(type1: &str, valp1: &mut Value, dval1: f64) {
    match type1.as_bytes()[0] {
        b's' => err!("can't change char header word"),
        b'h' => valp1.set_h(dval1 as i32 as i16),
        b'u' => valp1.set_u(dval1 as u32 as u16),
        b'l' => valp1.set_l(dval1 as i64),
        b'v' => valp1.set_v(dval1 as u64),
        b'i' => valp1.set_i(dval1 as i32),
        b'p' => valp1.set_p(dval1 as u32),
        b'f' => valp1.set_f(dval1 as f32),
        b'd' => valp1.set_d(dval1),
        _ => err!("unknown type %s", type1),
    }
}
