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
//! B2A - convert binary floats to ASCII. Port of par/main/data_conversion/b2a.c (Revision 1.16).
//! Credits: CWP: Jack K. Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " B2A - convert binary floats to ascii\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " b2a <stdin >stdout \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tn1=2\t\tfloats per line in output file \t\t",
    "       format=0\tscientific notation\t \t\t",
    " \t\t\t=1 long decimal float form\t\t",
    " \t\t\t\t\t\t\t\t",
    " \toutpar=/dev/tty output parameter file, contains the\t",
    "\t\t\tnumber of lines (n=)\t\t\t",
    "                       other choices for outpar are: /dev/tty, ",
    "                       /dev/stderr, or a name of a disk file   ",
    " \t\t\t\t\t\t\t\t",
    " Note: \t\t\t\t\t\t\t",
    " Parameter:\t\t\t\t\t\t\t",
    "  format=0 uses printf(\"%15.10e \", x[i1])\t\t\t",
    "  format=1 uses printf(\"%15.15f \", x[i1])\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let outparfp = efopen(&outpar, "w");
    let mut n1 = 0;
    if getparint("n1", &mut n1) == 0 { n1 = 2; }
    let mut format = 0;
    if getparint("format", &mut format) == 0 { format = 0; }
    let mut x = Vec::<f32>::new();
    if x.try_reserve_exact(n1 as usize).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    x.resize(n1 as usize, 0.0);
    checkpars();

    let mut n2 = 0i32;
    loop {
        let n1read = efread(&mut x, &stdin());
        if n1read == 0 { break; }
        if n1read != n1 as usize {
            err!("out of data in forming line #%d", n2.wrapping_add(1));
        }
        if format == 0 {
            for i1 in 0..n1 as usize { printf!("%15.10e ", x[i1]); }
        } else {
            for i1 in 0..n1 as usize { printf!("%20.20f ", x[i1]); }
        }
        putchar('\n' as i32);
        n2 = n2.wrapping_add(1);
    }
    fprintf!(outparfp, "n=%d\n", n2);
    exit(CWP_Exit());
}
