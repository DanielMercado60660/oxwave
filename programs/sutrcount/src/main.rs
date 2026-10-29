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

//! SUTRCOUNT - count traces in an SU input stream.
//! Port of su/main/headers/sutrcount.c (Revision 1.5).
//! Credits: B. Nemeth, Potash Corporation, Saskatchewan;
//! given to CWP in 2008 with permission of Potash Corporation.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUTRCOUNT - SU program to count the TRaces in infile\t\t",
    "       \t\t\t\t\t\t\t",
    "   sutrcount < infile\t\t\t\t\t     \t",
    " Required parameters:\t\t\t\t\t\t",
    "       none\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    "    outpar=stdout\t\t\t\t\t\t",
    " Notes:       \t\t\t\t\t\t\t",
    " Once you have the value of ntr, you may set the ntr header field",
    " via:      \t\t\t\t\t\t\t",
    "       sushw key=ntr a=NTR < datain.su  > dataout.su \t\t",
    " Where NTR is the value of the count obtained with sutrcount \t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut ntr: i32 = 0;
    let mut outpar = String::new();

    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    if getparstring("outpar", &mut outpar) == 0 {
        outpar = "/dev/stdout".to_string();
    }
    let outparfp = efopen(&outpar, "w");
    checkpars();
    loop {
        ntr = ntr.wrapping_add(1);
        if gettr(&mut tr) == 0 { break; }
    }
    fprintf!(outparfp, "%d", ntr);
    exit(CWP_Exit());
}
