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

//! MKPARFILE - convert ascii to par file format.
//! Port of par/main/parameter_file_utilities/mkparfile.c (Revision 1.7).
//! Credits: CWP: Jack K. Cohen.

use oxwave_core::prelude::*;
use std::ffi::{c_char, c_int, CString};

unsafe extern "C" {
    fn sscanf(s: *const c_char, format: *const c_char, ...) -> c_int;
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " MKPARFILE - convert ascii to par file format \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " mkparfile <stdin >stdout \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " \tstring1=\"par1\"\tfirst par string\t\t\t",
    " \tstring2=\"par2\"\tsecond par string\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " This is a tool to convert values written line by line to parameter \t",
    " vectors in the form expected by getpar.  For example, if the input\t",
    " file looks like:\t\t\t\t\t\t\t",
    " \tt0 v0\t\t\t\t\t\t\t\t",
    " \tt1 v1\t\t\t\t\t\t\t\t",
    "\t...\t\t\t\t\t\t\t\t",
    " then\t\t\t\t\t\t\t\t\t",
    "\tmkparfile <input >output string1=tnmo string2=vnmo\t\t",
    " yields:\t\t\t\t\t\t\t\t",
    "\ttnmo=t0,t1,...\t\t\t\t\t\t\t",
    "\tvnmo=v0,v1,...\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut n2: i32 = 0;
    let mut x1: f32 = 0.0;
    let mut x2: f32 = 0.0;
    let mut string1 = String::new();
    let mut string2 = String::new();

    initargs(&argv);
    requestdoc(1, SDOC);
    if getparstring("string1", &mut string1) == 0 { string1 = "par1".to_string(); }
    if getparstring("string2", &mut string2) == 0 { string2 = "par2".to_string(); }
    let datafp = etmpfile();
    checkpars();

    // BUFSIZ is 1024 in the original macOS build.
    let mut buf = match fgets(1024, &stdin()) {
        Some(buf) => buf,
        None => err!("no data found"),
    };
    loop {
        let cbuf = CString::new(buf.split('\0').next().unwrap_or("")).unwrap();
        let scanned = unsafe { sscanf(cbuf.as_ptr(), c"%f %f".as_ptr(), &mut x1, &mut x2) };
        if scanned == 2 {
            if n2 == 0 { printf!("%s=%g", string1, x1); }
            else { printf!(",%g", x1); }
            efwrite(&[x2], &datafp);
        } else {
            err!("line #%d: scan failed:\n%s", n2 + 1, buf);
        }
        n2 += 1;
        match fgets(1024, &stdin()) {
            Some(next) => buf = next,
            None => break,
        }
    }
    putchar('\n' as i32);

    erewind(&datafp);
    let mut value = [0.0f32];
    efread(&mut value, &datafp);
    printf!("%s=%g", string2, value[0]);
    for _i2 in 1..n2 {
        efread(&mut value, &datafp);
        printf!(",%g", value[0]);
    }
    putchar('\n' as i32);
    exit(CWP_Exit());
}
