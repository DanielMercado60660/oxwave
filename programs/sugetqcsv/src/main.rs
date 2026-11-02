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


//! SUGETQCSV - Write trace header key values and sample values to a Q-file.
//! Port of su/main/statics/sugetqcsv.c (Revision 1.01).
//! Credits: Andre Latour.

use oxwave_core::prelude::*;
use oxwave_header_case::{GetCase, fromhead};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUGETQCSV - Write Trace Header Key Values and Sample Values to a Q-file    ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  sugetqcsv [parameters].                                                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t         \t\t\t\t\t\t     ",
    "                                                                            ",
    " qout=qhed.csv  Output csv Q-file name (default is qhed.csv)                ",
    "\t\t\t\t\t\t\t\t\t     ",
    " keys=     Key name list. Must be trace header key names, or must be        ",
    "           one of the special names below:                                  ",
    "             mx = midpoint X coordinate (sx+gx)/2.0                         ",
    "             my = midpoint Y coordinate (sy+gy)/2.0                         ",
    "             mg = midpoint station (grnlof+gaps)/2.0                        ",
    "     Note: Names sx,sy,gx,gy,mx,my are scaled by scalco before output and   ",
    "           selev,gelev,sdepth,sdel,gdel,swdep,gwdep are scaled by scalel.   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " formxy=%.20g  The C format code for printing all values to the q-records.  ",
    "              Note that this default format prints up to 20 digits but      ",
    "              not trailing zeroes to the right of the decimal point.        ",
    "              (This format code is the default because some coordinates     ",
    "               scaled by scalco MIGHT need it for full representation).     ",
    "                                                                            ",
    "                                                                            ",
    " fsamp=0   First sample value to output. The first sample is number 0.      ",
    "                                                                            ",
    " lsamp=-1  Last sample value to output.                                     ",
    "           The default (-1) means do not output any sample values.          ",
    "           If this value is greater than the number of samples in the       ",
    "           input traces, the output is restricted to the number of samples. ",
    "     Note: In accordance with q-file standards, sample values are given the ",
    "           names numb0, numb1, numb2... (numb and sample number appended).  ",
    "      ***  To output sample values, all input traces must have the          ",
    "           same amount of samples as the first input trace.                 ",
    "                                                                            ",
    " formsv=%.10g  The C format code for printing sample values to q-records.   ",
    "               Note that this default format prints up to 10 digits but     ",
    "               not trailing zeroes to the right of the decimal point.       ",
    "               (This format code is the default because sample values are   ",
    "               single precision float and need it for full representation). ",
    "                                                                            ",
    "                                                                            ",
];

fn scaledfromhead(tr: &Trace, k: i32) -> f64 {
    let (mut dval, scale) = match k {
        101 => (tr.sx as f64, tr.scalco), 102 => (tr.sy as f64, tr.scalco),
        103 => (tr.gx as f64, tr.scalco), 104 => (tr.gy as f64, tr.scalco),
        105 => (0.5 * (tr.sx as f64 + tr.gx as f64), tr.scalco),
        106 => (0.5 * (tr.sy as f64 + tr.gy as f64), tr.scalco),
        107 => return 0.5 * (tr.grnlof as f64 + tr.gaps as f64),
        111 => (tr.selev as f64, tr.scalel), 112 => (tr.gelev as f64, tr.scalel),
        113 => (tr.sdepth as f64, tr.scalel), 114 => (tr.sdel as f64, tr.scalel),
        115 => (tr.gdel as f64, tr.scalel), 116 => (tr.swdep as f64, tr.scalel),
        117 => (tr.gwdep as f64, tr.scalel),
        1..=80 => return fromhead(tr, k),
        // The SU build returns zero for this uninitialized fromhead result.
        _ => return 0.0,
    };
    if scale > 1 { dval *= scale as f64; }
    else if scale < 0 { dval /= -(scale as i32) as f64; }
    dval
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut qout = String::new();
    if countparval("qout") > 0 { getparstring("qout", &mut qout); }
    else { qout = "qhed.csv".to_string(); }
    let fp_q = match fopen(&qout, "w") {
        Some(fp) => fp,
        None => err!("error: output q-file did not open correctly."),
    };
    let numd = countparval("keys");
    if numd < 1 { err!("**** Error: andre some keys is needed (default?) "); }
    let mut keyn = Vec::new();
    getparstringarray("keys", &mut keyn);
    let mut kase = Vec::new();
    for key in &keyn {
        let k = match key.as_str() {
            "sx" => 101, "sy" => 102, "gx" => 103, "gy" => 104,
            "mx" => 105, "my" => 106, "mg" => 107,
            "selev" => 111, "gelev" => 112, "sdepth" => 113,
            "sdel" => 114, "gdel" => 115, "swdep" => 116, "gwdep" => 117,
            _ => GetCase(key),
        };
        if k < 1 { err!("**** Error: Specified dimension name %s is not recognized.", key); }
        kase.push(k);
    }
    let mut formxy = String::new();
    if getparstring("formxy", &mut formxy) == 0 { formxy = "%.20g".to_string(); }
    let formxylong = format!(",{}", formxy);
    let mut formsv = String::new();
    if getparstring("formsv", &mut formsv) == 0 { formsv = "%.10g".to_string(); }
    let formsvlong = format!(",{}", formsv);
    let mut fsamp = 0;
    let mut lsamp = -1;
    if getparint("fsamp", &mut fsamp) == 0 { fsamp = 0; }
    if fsamp < 0 { err!("**** Error: fsamp cannot be less than 0"); }
    if getparint("lsamp", &mut lsamp) == 0 { lsamp = -1; }
    if lsamp < -1 { err!("**** Error: lsamp cannot be less than -1"); }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("Error: cannot get first trace"); }
    if lsamp > -1 {
        if lsamp > tr.ns as i32 { lsamp = tr.ns as i32 - 1; }
        if fsamp > lsamp { err!("**** Error: fsamp must be less than or equal to lsamp"); }
    }
    fprintf!(fp_q, "C_SU_MATCH,%s\nC_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID", keyn[0]);
    for _ in 0..numd { fprintf!(fp_q, ",%s", formxy); }
    for _ in fsamp..=lsamp { fprintf!(fp_q, ",%s", formsv); }
    fprintf!(fp_q, "\nC_SU_NAMES\nC_SU_ID");
    for key in &keyn { fprintf!(fp_q, ",%s", key); }
    for i in fsamp..=lsamp { fprintf!(fp_q, ",numb%d", i); }
    fprintf!(fp_q, "\n");
    loop {
        fprintf!(fp_q, "Q");
        for &k in &kase { fprintf!(fp_q, &formxylong, scaledfromhead(&tr, k)); }
        for i in fsamp..=lsamp { fprintf!(fp_q, &formsvlong, tr.data[i as usize] as f64); }
        fprintf!(fp_q, "\n");
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
