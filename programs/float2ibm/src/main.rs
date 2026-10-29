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

//! FLOAT2IBM - convert native binary floats to IBM tape floats.
//! Port of par/main/data_conversion/float2ibm.c (Revision 1.4).
//! Credits: CWP: John Stockwell, based on code by Jack K. Cohen;
//! conversion routine: CWP: Brian Sumner.

use oxwave_core::prelude::*;
use oxwave_core::segy::float_to_ibm;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " FLOAT2IBM - convert native binary FLOATS to IBM tape FLOATS\t",
    " \t\t\t\t\t\t\t\t",
    " float2ibm <stdin >stdout \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " endian=\tbyte order of your system (autodetected)\t",
    " outpar=/dev/tty output parameter file, contains the\t\t",
    "\t\t\tnumber of values (n=)\t\t\t",
    "\t\t       other choices for outpar are: /dev/tty, ",
    "\t\t       /dev/stderr, or a name of a disk file   ",
    " \t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " endian=1 (big endian) endian=0 (little endian) byte order \t",
    " You probably will not have to set this, as the byte order of  ",
    " your system is autodetected by the program. \t\t\t",
    " This program is usable for writing SEG Y traces with the headers",
    " stripped off.\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    match filestat(STDOUT) {
        FileType::BadFileType => {
            warn!("stdout is illegal filetype");
            pagedoc(SDOC);
        }
        FileType::Tty => {
            warn!("stdout can't be tty");
            pagedoc(SDOC);
        }
        _ => {}
    }

    let mut endian = 0;
    if getparint("endian", &mut endian) == 0 {
        endian = if cfg!(target_endian = "big") { 1 } else { 0 };
    }
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 {
        outpar = "/dev/tty".to_string();
    }
    let outparfp = efopen(&outpar, "w");
    checkpars();

    let mut x = [0.0f32; 1];
    let mut n1: i32 = 0;
    while efread(&mut x, &stdin()) != 0 {
        float_to_ibm(&mut x, 1, endian);
        fwrite(&x, &stdout());
        n1 = n1.wrapping_add(1);
    }
    fprintf!(outparfp, "n=%d\n", n1);
    exit(CWP_Exit());
}
