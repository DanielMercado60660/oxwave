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

//! BHEDTOPAR - convert a binary tape header to parameter-file format.
//! Port of su/main/headers/bhedtopar.c (Revision 1.11).
//! Credits: CWP: John Stockwell, 11 Nov 1994; revised by Stewart A. Levin.

use oxwave_core::prelude::*;
use oxwave_core::segy::{Bhed, BHDR, BHED_NKEYS, BNYBYTES, getbhval, swapbhval};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " BHEDTOPAR - convert a Binary tape HEaDer file to PAR file format\t",
    " \t\t\t\t\t\t\t\t\t",
    "     bhedtopar < stdin outpar=parfile\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameter:\t\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tswap=0 \t\t\t=1 to swap bytes\t\t\t",
    " \toutpar=/dev/tty\t\t=parfile  name of output param file\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " This program dumps the contents of a SEGY binary tape header file, as ",
    " would be produced by segyread and segyhdrs to a file in \"parfile\" format.",
    " A \"parfile\" is an ASCII file containing entries of the form param=value.",
    " Here \"param\" is the keyword for the binary tape header field and\t",
    " \"value\" is the value of that field. The parfile may be edited as\t",
    " any ASCII file. The edited parfile may then be made into a new binary tape ",
    " header file via the program    setbhed.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " See    sudoc  setbhed   for examples\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let mut swap = 0i32;
    if getparint("swap", &mut swap) == 0 { swap = 0; }
    checkpars();
    let outparfp = efopen(&outpar, "w");

    // XDR reads 4-byte words. On EOF the first three fields remain zero;
    // subsequent pairs reuse the last successfully decoded word. The original
    // build starts that temporary at 1 when no pair can be read.
    let mut raw = [0u8; BNYBYTES];
    let n = fread(&mut raw, &stdin());
    for index in 0..3 {
        if n < (index + 1) * 4 { raw[index * 4..index * 4 + 4].fill(0); }
    }
    let mut last = [0, 0, 0, 1];
    for offset in (12..BNYBYTES).step_by(4) {
        if n >= offset + 4 { last.copy_from_slice(&raw[offset..offset + 4]); }
        else { raw[offset..offset + 4].copy_from_slice(&last); }
    }
    let mut bh = Bhed::from_xdr(&raw);
    for index in 0..BHED_NKEYS {
        let mut val = Value::default();
        if swap != 0 { swapbhval(&mut bh, index as i32); }
        getbhval(&bh, index as i32, &mut val);
        fprintfparval(&outparfp, BHDR[index].key, BHDR[index].type_, val);
    }
    efclose(&outparfp);
    exit(CWP_Exit());
}

fn fprintfparval(stream: &CFile, key: &str, type_: &str, val: Value) {
    match type_ {
        "h" => fprintf!(stream, "%s=%d\n", key, val.h()),
        "i" => fprintf!(stream, "%s=%d\n", key, val.i()),
        _ => err!("fprintfparval: unknown type %s", type_),
    };
}
