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
//! SUPAD - Pad zero traces. Port of su/main/windowing_sorting_muting/supad.c (Revision 1.3).
//! Credits: Florian Bleibinhaus, U Salzburg, Austria.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUPAD - Pad zero traces\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  supad <stdin >stdout min= max= [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "  min=\t\t\ttrace key start\t\t\t\t\t",
    "  max=\t\t\ttrace key end\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "  key1=ep\t\tpanel key \t\t\t\t\t",
    "  key2=tracf\t\ttrace key \t\t\t\t\t",
    "  key3=trid\t\tflag key\t\t\t\t\t",
    "  val3=2\t\tvalue assigned to padded traces\t\t\t",
    "  d=1\t\t\ttrace key spacing\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    "  In contrast to most SU codes, supad recognizes panels, or ensembles.\t",
    "  If the input consists of several panels, each panel will be padded\t",
    "  individually.\t\t\t\t\t\t\t",
    "  key1 and key2 are the primary and secondary sort key of the data set.",
    "  The sort order of key1 does not matter at all.\t\t\t",
    "  The sort order of key2 must be monotonous - if key2 is descending,\t",
    "\tsupply a negative value for the spacing d.\t\t\t",
    "  Traces with a key2-value outside the min/max range will be lost. \t",
    "  Traces with a key2-value that is not a multiple of the spacing from\t",
    "\tthe min-value (the max-value, if the spacing is negative) will\t",
    "\tnot be lost. Instead, they will shift the series of key2-values.",
    "  By default the dead trace flag will be raised for the padded traces.\t",
    "  This should make it easy to remove the zero traces later on, if need be.",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "\tsuplane | supad min=1 max=40 key1=offset key2=tracr | ...\t",
    "\t... appends eight empty traces.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tsuplane | supad min=1 max=32 key1=offset key2=tracr d=0.5 | ...\t",
    "\t... inserts a zero trace after each trace (even though the\t",
    "\theader tracr is integer and cannot properly store the floats)\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tsuplane | supad min=1 max=32 | ...\t\t\t\t",
    "\t... produces an error because the panel and trace key are all 0.",
    "\t\t\t\t\t\t\t\t\t",
];

fn assgnval(type1: &str, val: &mut Value, dval: f64) {
    match type1.as_bytes()[0] {
        b's' => err!("can't change char header word"),
        b'h' => val.set_h(dval as i16),
        b'u' => val.set_u(dval as u16),
        b'l' => val.set_l(dval as i64),
        b'v' => val.set_v(dval as u64),
        b'i' => val.set_i(dval as i32),
        b'p' => val.set_p(dval as u32),
        b'f' => val.set_f(dval as f32),
        b'd' => val.set_d(dval),
        _ => err!("unknown type %s", type1),
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut nulltr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key1 = String::new();
    let mut key2 = String::new();
    let mut key3 = String::new();
    let mut dval3 = 0.0f64;
    let mut dx = 0.0f64;
    let mut dmin = 0.0f64;
    let mut dmax = 0.0f64;
    if getparstring("key1", &mut key1) == 0 { key1 = "ep".into(); }
    if getparstring("key2", &mut key2) == 0 { key2 = "tracf".into(); }
    if getparstring("key3", &mut key3) == 0 { key3 = "trid".into(); }
    if getpardouble("val3", &mut dval3) == 0 { dval3 = 2.0; }
    if getpardouble("d", &mut dx) == 0 { dx = 1.0; }
    if getpardouble("min", &mut dmin) == 0 { err!("need lower panel boundary MIN"); }
    if getpardouble("max", &mut dmax) == 0 { err!("need upper panel boundary MAX"); }
    checkpars();
    if dx == 0.0 { err!("trace spacing d cannot be zero"); }
    if dmax < dmin { err!("max needs to be greater than min"); }
    let isgn = if dx < 0.0 { -1.0 } else { 1.0 };
    let type1 = hdtype(&key1);
    let type2 = hdtype(&key2);
    let type3 = hdtype(&key3);
    let index1 = getindex(&key1);
    let index2 = getindex(&key2);
    let index3 = getindex(&key3);
    let mut val1 = Value::default();
    let mut val2 = Value::default();
    let mut val3 = Value::default();
    let mut dval1 = 0.0;
    let mut dval2 = 0.0;
    let mut panelno = 0.0;
    let mut traceno = 0.0;
    let mut iflag = 1;
    while iflag >= 0 {
        if gettr(&mut tr) != 0 {
            gethval(&tr, index1, &mut val1);
            gethval(&tr, index2, &mut val2);
            dval1 = vtod(type1, val1);
            dval2 = vtod(type2, val2);
            let nt = tr.ns as usize;
            nulltr.data[..nt].fill(0.0);
            if iflag == 1 {
                panelno = dval1;
                traceno = dmin - dx;
            }
            iflag = 0;
            if dval2 < dmin || dval2 > dmax { iflag = 2; }
        } else {
            iflag = -1;
        }
        if panelno != dval1 || iflag == -1 {
            let mut c = traceno + dx;
            while isgn * c <= isgn * dmax {
                assgnval(type2, &mut val2, c);
                puthval(&mut nulltr, index2, &val2);
                assgnval(type3, &mut val3, dval3);
                puthval(&mut nulltr, index3, &val3);
                puttr(&mut nulltr);
                c += dx;
            }
            traceno = dmin - dx;
            panelno = dval1;
        }
        if iflag == 0 {
            nulltr.hdr = tr.hdr;
            let mut c = traceno + dx;
            while isgn * c < isgn * dval2 {
                assgnval(type2, &mut val2, c);
                puthval(&mut nulltr, index2, &val2);
                assgnval(type3, &mut val3, dval3);
                puthval(&mut nulltr, index3, &val3);
                puttr(&mut nulltr);
                c += dx;
            }
        }
        if iflag == 0 {
            puttr(&mut tr);
            panelno = dval1;
            traceno = dval2;
        }
    }
    exit(CWP_Exit());
}
