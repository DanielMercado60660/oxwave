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

//! A2I - convert ASCII to binary integers. Port of par/main/data_conversion/a2i.c (Revision 1.5).
//! Credits: Potash Corporation: Balazs Nemeth; based on a2b.c by CWP: Jack K. Cohen, Dave Hale.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use std::ffi::{c_char, c_int, CString};

unsafe extern "C" {
    fn sscanf(s: *const c_char, format: *const c_char, ...) -> c_int;
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " A2I - convert Ascii to binary Integers\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " a2i <stdin >stdout outpar=/dev/tty \t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tn1=2\t\tintegers per line in input file\t\t",
    " \t\t\t\t\t\t\t\t",
    " \toutpar=/dev/tty\toutput parameter file, contains the\t",
    "\t\t\tnumber of lines (n=)\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    // C stdout exists before outpar is opened; register it in the same order.
    let _ = stdout();
    match filestat(STDOUT) {
        FileType::BadFileType => { warn!("stdout is illegal filetype"); pagedoc(SDOC); }
        FileType::Tty => { warn!("stdout can't be tty"); pagedoc(SDOC); }
        _ => {}
    }
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let outparfp = efopen(&outpar, "w");
    let mut n1 = 0;
    if getparint("n1", &mut n1) == 0 { n1 = 2; }
    if n1 < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut x = vec![0i32; n1 as usize];
    let mut n2 = 0i32;
    checkpars();
    while let Some(buf) = fgets(131072, &stdin()) {
        let bytes = buf.as_bytes();
        let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        let mut p = 0;
        for i1 in 0..n1 {
            while p < end && (bytes[p] == b' ' || bytes[p] == b'\t') { p += 1; }
            let input = CString::new(&bytes[p..end]).unwrap();
            let mut flt = 0.0f32;
            if unsafe { sscanf(input.as_ptr(), c"%f".as_ptr(), &mut flt) } != 1 {
                err!("line #%d: scan failed on number #%d:\n%s", n2 + 1, i1 + 1, &buf[..end]);
            }
            x[i1 as usize] = nint(flt as f64);
            while p < end && bytes[p] != b' ' && bytes[p] != b'\t' { p += 1; }
        }
        efwrite(&x, &stdout());
        n2 += 1;
    }
    fprintf!(&outparfp, "n=%d\n", n2);
    exit(CWP_Exit());
}
