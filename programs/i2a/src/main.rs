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

//! I2A - convert binary integers to ascii.
//! Port of par/main/data_conversion/i2a.c (Revision 1.4).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: Potash Corporation: c. 2008, Balazs Nemeth, Saskatoon, Saskatchewan.
//! Based on b2a.c by CWP: Jack K. Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " I2A - convert binary integers to ascii\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " i2a <stdin >stdout \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tn1=2\t\tfloats per line in output file \t\t",
    " \t\t\t\t\t\t\t\t",
    " \toutpar=/dev/tty\toutput parameter file, contains the\t",
    "\t\t\tnumber of lines (n=)\t\t\t",
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
    let mut n1 = 0i32;
    if getparint("n1", &mut n1) == 0 {
        n1 = 2;
    }
    // Negative sizes convert to an impossible size_t allocation in C.
    if n1 < 0 {
        err!("%s: malloc failed (Cannot allocate memory)", "ealloc.c");
    }
    let mut x = vec![0i32; n1 as usize];
    let mut n2 = 0i32;
    checkpars();

    loop {
        let n1read = efread(&mut x, &stdin());
        if n1read == 0 { break; }
        if n1read != n1 as usize {
            err!("out of data in forming line #%d", n2.wrapping_add(1));
        }
        for i1 in 0..n1 as usize {
            printf!("%15d ", x[i1]);
        }
        putchar(b'\n' as i32);
        n2 = n2.wrapping_add(1);
    }
    fprintf!(outparfp, "n=%d\n", n2);
    exit(CWP_Exit());
}
