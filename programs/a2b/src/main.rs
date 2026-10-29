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

//! A2B - convert ASCII floats to binary.
//! Port of par/main/data_conversion/a2b.c (Revision 1.16).
//! Credits: CWP: Jack K. Cohen, Dave Hale; Hans Ecke: float-wise reading.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " A2B - convert ascii floats to binary \t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " a2b <stdin >stdout outpar=/dev/null \t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tn1=2\t\tfloats per line in input file\t\t",
    " \t\t\t\t\t\t\t\t",
    " \toutpar=/dev/null output parameter file, contains the\t",
    "\t\t\tnumber of lines (n=)\t\t\t",
    " \t\t\tother choices for outpar are: /dev/tty,\t",
    " \t\t\t/dev/stderr, or a name of a disk file\t",
];

const N1_DEFAULT: i32 = 2;
const OUTPAR_DEFAULT: &str = "/dev/null";

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    match filestat(STDOUT) {
        FileType::BadFileType => { warn!("stdout is illegal filetype"); pagedoc(SDOC); }
        FileType::Tty => { warn!("stdout can't be tty"); pagedoc(SDOC); }
        _ => {}
    }
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = OUTPAR_DEFAULT.to_string(); }
    let outparfp = efopen(&outpar, "w");
    let mut n1 = 0i32;
    if getparint("n1", &mut n1) == 0 { n1 = N1_DEFAULT; }
    let mut x = vec![0.0f32; n1 as usize];
    checkpars();
    let mut vnum = 0usize;
    let mut vec = 0u64;
    let input = stdin();
    loop {
        // Distinguish EOF before a conversion from a failed conversion that
        // consumed a partial token and reached EOF.
        let mut c = fgetc(&input);
        while c != EOF && (c as u8).is_ascii_whitespace() { c = fgetc(&input); }
        if c == EOF {
            if vnum != 0 {
                err!("vector #%lu, float #%lu: we encountered an EOF but the vector has not been read in completely", vec + 1, vnum as u64 + 1);
            }
            break;
        }
        ungetc(c, &input);
        match fscanf_f32(&input) {
            None => err!("could not read float in vector #%lu, float #%lu", vec + 1, vnum as u64 + 1),
            Some(value) => {
                // C writes past its zero-sized allocation for n1=0; retain a slot.
                if vnum >= x.len() { x.push(value); } else { x[vnum] = value; }
            }
        }
        vnum += 1;
        if vnum as u64 == n1 as u64 {
            vnum = 0;
            fwrite(&x, &stdout());
            vec += 1;
        }
    }
    fprintf!(outparfp, "n2=%lu\n", vec);
    exit(CWP_Exit());
}
