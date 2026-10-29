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

//! FTNUNSTRIP - convert C binary floats to Fortran style floats.
//! Port of par/main/data_conversion/ftnunstrip.c (Revision 1.7).
//! Credits: CWP: John Stockwell, Feb 1998, based on ftnstrip by Jack K. Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " FTNUNSTRIP - convert C binary floats to Fortran style floats\t",
    " \t\t\t\t\t\t\t\t",
    " ftnunstrip <stdin >stdout \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tn1=1\t\tfloats per line in output file \t\t",
    " \t\t\t\t\t\t\t\t",
    " \toutpar=/dev/tty output parameter file, contains the\t",
    "\t\t\tnumber of lines (n=)\t\t\t",
    " \t\t\tother choices for outpar are: /dev/tty,\t",
    " \t\t\t/dev/stderr, or a name of a disk file\t",
    " \t\t\t\t\t\t\t\t",
    " Notes: This program assumes that the record length is constant",
    " throughout the input and output files. \t\t\t",
    " In fortran code reading these floats, the following implied\t",
    " do loop syntax would be used: \t\t\t\t",
    "        DO i=1,n2\t\t\t\t\t\t",
    "                 READ (10) (someARRAY(j), j=1,n1) \t\t",
    "        END DO\t\t\t\t\t\t\t",
    " Here n1 is the number of samples per record, n2 is the number ",
    " of records, 10 is some default file (fort.10, for example), and",
    " someArray(j) is an array dimensioned to size n1\t\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 {
        outpar = "/dev/tty".to_string();
    }
    let outparfp = efopen(&outpar, "w");
    let mut n1: i32 = 1;
    if getparint("n1", &mut n1) == 0 { n1 = 1; }
    let n1bytes = 4i32.wrapping_mul(n1);
    if n1 < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut x = vec![0.0f32; n1 as usize];
    let mut n2: i32 = 0;
    checkpars();
    loop {
        let n1read = efread(&mut x, &stdin());
        if n1read == 0 { break; }
        if n1read != n1 as usize {
            err!("out of data in forming line #%d", n2.wrapping_add(1));
        }
        efwrite(&[n1bytes], &stdout());
        efwrite(&x, &stdout());
        efwrite(&[n1bytes], &stdout());
        n2 = n2.wrapping_add(1);
    }
    fprintf!(&outparfp, "n=%d\n", n2);
    exit(CWP_Exit());
}
