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

//! PSBBOX - change BoundingBOX of existing PostScript file.
//! Port of psplot/main/psbbox.c (Revision 1.7).
//! Credits: CWP: Craig Artley, Colorado School of Mines, 08/30/91.

use oxwave_core::prelude::*;
use oxwave_core::cio::fputs_bytes;
use oxwave_cwp::cwp::nint;
use std::ffi::{c_char, c_int};

unsafe extern "C" {
    fn sscanf(s: *const c_char, format: *const c_char, ...) -> c_int;
}
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " PSBBOX - change BoundingBOX of existing PostScript file\t",
    "\t\t\t\t\t\t\t\t",
    " psbbox < PostScriptfile [optional parameters] > PostScriptfile",
    "\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " llx=\t\tnew llx\t\t\t\t\t\t",
    " lly=\t\tnew lly\t\t\t\t\t\t",
    " urx=\t\tnew urx\t\t\t\t\t\t",
    " ury=\t\tnew ury\t\t\t\t\t\t",
    " verbose=1\t=1 for info printed on stderr (0 for no info)\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 1;
    getparint("verbose", &mut verbose);
    let infp = stdin();
    let outfp = stdout();
    let mut bbox_found = false;
    loop {
        // fgets reads at most MAXLINE-1 bytes, including any embedded NUL.
        let mut line = Vec::new();
        while line.len() < 2047 {
            let c = fgetc(&infp);
            if c == EOF { break; }
            line.push(c as u8);
            if c == b'\n' as i32 { break; }
        }
        if line.is_empty() { break; }
        // C string operations stop at the first embedded NUL.
        let end = line.iter().position(|&c| c == 0).unwrap_or(line.len());
        line.truncate(end);
        if line.first() == Some(&b'%') && !bbox_found
            && line.windows(14).any(|s| s == b"%%BoundingBox:") {
            if line.windows(5).any(|s| s == b"atend") { continue; }
            let mut c_line = line.clone();
            c_line.push(0);
            let (mut llx, mut lly, mut urx, mut ury) = (0i32, 0i32, 0i32, 0i32);
            let (mut fllx, mut flly, mut furx, mut fury) = (0f32, 0f32, 0f32, 0f32);
            unsafe {
                if sscanf(c_line.as_ptr().cast(), c"%*s %d %d %d %d".as_ptr(),
                    &mut llx, &mut lly, &mut urx, &mut ury) == 4 {
                } else if sscanf(c_line.as_ptr().cast(), c"%*s %f %f %f %f".as_ptr(),
                    &mut fllx, &mut flly, &mut furx, &mut fury) == 4 {
                    llx = nint(fllx as f64);
                    lly = nint(flly as f64);
                    urx = nint(furx as f64);
                    ury = nint(fury as f64);
                } else {
                    err!("Error reading BoundingBox!\n");
                }
            }
            if verbose != 0 {
                fputs("Original:  ", &stderr());
                fputs_bytes(&line, &stderr());
            }
            getparint("llx", &mut llx);
            getparint("lly", &mut lly);
            getparint("urx", &mut urx);
            getparint("ury", &mut ury);
            line = sprintf!("%s %d %d %d %d\n", "%%BoundingBox:", llx, lly, urx, ury).into_bytes();
            bbox_found = true;
            if verbose != 0 {
                fputs("Updated:   ", &stderr());
                fputs_bytes(&line, &stderr());
            }
        }
        fputs_bytes(&line, &outfp);
    }
    exit(CWP_Exit());
}
