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


//! SUVLENGTH - Adjust variable length traces to common length.
//! Port of su/main/operations/suvlength.c (Revision 1.8).
//! Credits (from the C source): CWP: Jack Cohen, John Stockwell.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUVLENGTH - Adjust variable length traces to common length   \t",
    " \t\t\t\t\t\t\t\t",
    " suvlength <vdata >stdout\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tns=\toutput number of samples (default: 1st trace ns)",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut ns: i32 = 0; /* samples on output traces */
    initargs(&argv);
    requestdoc(1, SDOC);
    if fvgettr(&stdin(), &mut tr) == 0 { err!("can't get first trace"); }
    if getparint("ns", &mut ns) == 0 { ns = tr.ns as i32; }
    checkpars();
    // Larger positive lengths overrun the fixed trace buffer in C.
    // Negative lengths skip padding and retain C's unsigned header conversion.
    if ns > SU_NFLTS as i32 {
        err!("ns=%d exceeds trace buffer capacity %d", ns, SU_NFLTS as i32);
    }
    loop {
        let nt = tr.ns as i32;
        if nt < ns { /* pad with zeros */
            tr.data[nt as usize..ns as usize].fill(0.0);
        }
        tr.ns = ns as u16;
        puttr(&mut tr);
        if fvgettr(&stdin(), &mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}

