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

//! FTNSTRIP - remove Fortran record delimiters.
//! Port of par/main/data_conversion/ftnstrip.c (Revision 1.6).
//! Credits: CWP: Jack K. Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " FTNSTRIP - convert a file of binary data plus record delimiters created",
    "      via Fortran to a file containing only binary values (as created via C)",
    " \t\t\t\t\t\t\t\t\t",
    " ftnstrip <ftn_data >c_data \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat: this code assumes the conventional Fortran format of header\t",
    "         and trailer integer containing the number of byte in the\t",
    "         record.  This is overwhelmingly common, but not universal.\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let input = stdin();
    let output = stdout();
    let mut n1bytes = [0i32; 1];
    while efread(&mut n1bytes, &input) != 0 {
        let size = n1bytes[0] as usize;
        let mut buf = Vec::<u8>::new();
        if buf.try_reserve_exact(size).is_err() {
            err!("%s: malloc failed (Cannot allocate memory)", "ealloc.c");
        }
        buf.resize(size, 0);
        if size == 0 {
            err!("%s: efread: fread given 0 item size", "subcalls.c");
        }
        efread(&mut buf, &input);
        efwrite(&buf, &output);
        efread(&mut n1bytes, &input);
    }
    exit(CWP_Exit());
}
