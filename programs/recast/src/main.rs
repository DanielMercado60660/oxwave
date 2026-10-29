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

//! RECAST - convert native binary data types. Port of par/main/data_conversion/recast.c (Revision 1.25).
//! Credits: CWP: John Stockwell, Jack K. Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " RECAST - RECAST data type (convert from one data type to another)\t",
    " \t\t\t\t\t\t\t\t\t",
    " recast <stdin [optional parameters]  >stdout \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " in=float\tinput type\t(float)\t\t\t\t\t",
    " \t\t=double\t\t(double)\t\t\t\t",
    " \t\t=int\t\t(int)\t\t\t\t\t",
    " \t\t=char\t\t(char)\t\t\t\t\t",
    "\t\t=uchar\t\t(unsigned char)\t\t\t\t",
    " \t\t=short\t\t(short)\t\t\t\t\t",
    " \t\t=long\t\t(long)\t\t\t\t\t",
    " \t\t=ulong\t\t(unsigned long)\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " out=double\toutput type\t(double)\t\t\t\t",
    " \t\t=float\t\t(float)\t\t\t\t\t",
    " \t\t=int\t\t(int)\t\t\t\t\t",
    " \t\t=char\t\t(char)\t\t\t\t\t",
    " \t\t=uchar\t\t(unsigned char)\t\t\t\t",
    " \t\t=short\t\t(short)\t\t\t\t\t",
    " \t\t=long\t\t(long)\t\t\t\t\t",
    " \t\t=ulong\t\t(unsigned long)\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " outpar=/dev/tty\toutput parameter file, contains the\t\t",
    "\t\t\t\tnumber of values (n1=)\t\t\t",
    " \t\t\tother choices for outpar are: /dev/tty,\t\t",
    " \t\t\t/dev/stderr, or a name of a disk file\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: Converting bigger types to smaller is hazardous. For float\t",
    "\t or double conversions to integer types, the results are \t",
    "\t rounded to the nearest integer.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

// The conversion macros follow the three macros in recast.c.
macro_rules! recast {
    ($input:ty, $output:ty, $n:ident) => {{
        let mut xin: [$input; 1] = [0 as $input];
        while efread(&mut xin, &stdin()) != 0 {
            $n = $n.wrapping_add(1);
            let xout: $output = xin[0] as $output;
            efwrite(&[xout], &stdout());
        }
    }};
}
macro_rules! recast_rounded {
    ($input:ty, $output:ty, $intermediate:ty, $n:ident) => {{
        let mut xin: [$input; 1] = [0 as $input];
        while efread(&mut xin, &stdin()) != 0 {
            $n = $n.wrapping_add(1);
            let x = xin[0] as f64;
            let rounded = if x > 0.0 { x + 0.5 } else { x - 0.5 };
            let xout: $output = (rounded as $intermediate) as $output;
            efwrite(&[xout], &stdout());
        }
    }};
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".into(); }
    let outparfp = efopen(&outpar, "w");
    let mut in_type = String::new();
    let mut out = String::new();
    if getparstring("in", &mut in_type) == 0 { in_type = "float".into(); }
    if getparstring("out", &mut out) == 0 { out = "double".into(); }
    checkpars();
    let types = ["float", "double", "int", "char", "uchar", "short", "long", "ulong"];
    if !types.contains(&in_type.as_str()) && !types.contains(&out.as_str()) {
        err!("%s or %s is an unsupported type", in_type, out);
    }
    let mut n: i32 = 0;
    match (in_type.as_str(), out.as_str()) {
        ("float", "double") => recast!(f32, f64, n),
        ("float", "int") => recast_rounded!(f32, i32, i64, n),
        ("float", "char") => recast_rounded!(f32, i8, i64, n),
        ("float", "uchar") => recast_rounded!(f32, u8, i64, n),
        ("float", "short") => recast_rounded!(f32, i16, i64, n),
        ("float", "long") => recast_rounded!(f32, i64, i64, n),
        ("float", "ulong") => recast_rounded!(f32, u64, u64, n),
        ("double", "float") => recast!(f64, f32, n),
        ("double", "int") => recast_rounded!(f64, i32, i64, n),
        ("double", "char") => recast_rounded!(f64, i8, i64, n),
        ("double", "uchar") => recast_rounded!(f64, u8, i64, n),
        ("double", "short") => recast_rounded!(f64, i16, i64, n),
        ("double", "long") => recast_rounded!(f64, i64, i64, n),
        ("double", "ulong") => recast_rounded!(f64, u64, u64, n),
        ("int", "float") => recast!(i32, f32, n),
        ("int", "double") => recast!(i32, f64, n),
        ("int", "char") => recast!(i32, i8, n),
        ("int", "uchar") => recast!(i32, u8, n),
        ("int", "short") => recast!(i32, i16, n),
        ("int", "long") => recast_rounded!(i64, i64, i64, n),
        ("int", "ulong") => recast!(i64, u64, n),
        ("char", "double") => recast!(i8, f64, n),
        ("char", "float") => recast!(i8, f32, n),
        ("char", "int") => recast!(i8, i32, n),
        ("char", "short") => recast!(i8, i16, n),
        ("char", "long") => recast!(i8, i64, n),
        ("char", "ulong") => recast!(i8, i64, n),
        ("uchar", "double") => recast!(u8, f64, n),
        ("uchar", "float") => recast!(u8, f32, n),
        ("uchar", "int") => recast!(u8, i32, n),
        ("uchar", "short") => recast!(u8, i16, n),
        ("uchar", "long") => recast!(u8, i64, n),
        ("uchar", "ulong") => recast!(u8, i64, n),
        ("short", "double") => recast!(i16, f64, n),
        ("short", "float") => recast!(i16, f32, n),
        ("short", "int") => recast!(i16, i32, n),
        ("short", "char") => recast!(i16, i8, n),
        ("short", "long") => recast!(i16, i64, n),
        ("long", "double") => recast!(i64, f64, n),
        ("long", "float") => recast!(i64, f32, n),
        ("long", "int") => recast!(i64, i32, n),
        ("long", "char") => recast!(i64, i8, n),
        ("long", "short") => recast!(i64, i16, n),
        ("ulong", "double") => recast!(u64, f64, n),
        ("ulong", "float") => recast!(u64, f32, n),
        ("ulong", "int") => recast!(u64, i32, n),
        ("ulong", "char") => recast!(u64, i8, n),
        ("ulong", "short") => recast!(u64, i16, n),
        _ if in_type == out => err!("intype and outtype are the same (%s).", in_type),
        _ => err!("Unsupported pair in=%s out=%s", in_type, out),
    }
    fprintf!(outparfp, "n=%d\n", n);
    exit(CWP_Exit());
}
