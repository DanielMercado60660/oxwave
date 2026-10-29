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

//! SEGYCLEAN - zero out unassigned portion of header.
//! Port of su/main/data_conversion/segyclean.c (Revision 1.13).
//!
//! Credits (from the C source):
//! CWP: Jack Cohen

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SEGYCLEAN - zero out unassigned portion of header\t\t",
    " \t\t\t\t\t\t\t\t",
    " segyclean <stdin >stdout \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Since \"foreign\" SEG-Y tapes may use the unassigned portion\t",
    " of the trace headers and since SU now uses it too, this\t",
    " program zeros out the fields meaningful to SU.\t\t",
    " \t\t\t\t\t\t\t\t",
    "  Example:\t\t\t\t\t\t\t",
    "  \tsegyread trmax=200 | segyclean | suximage\t\t",
    "\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();

    /* Initialize */
    initargs(&argv);
    requestdoc(1, SDOC);

    while gettr(&mut tr) != 0 {
        tr.f1 = 0.0;
        tr.d1 = 0.0;
        tr.f2 = 0.0;
        tr.d2 = 0.0;
        tr.ungpow = 0.0;
        tr.unscale = 0.0;
        tr.ntr = 0;
        tr.mark = 0;
        puttr(&mut tr);
    }

    exit(CWP_Exit());
}
