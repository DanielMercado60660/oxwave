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

//! SUR2AI - reflectivity to acoustic impedance.
//! Port of su/main/amplitudes/sur2ai.c.
//! Credits: Shuki Ronen, 2017. Partial implementation of Aki & Richards,
//! volume II, pages 661-662; transmission effects are ignored.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " sur2ai <stdin >sdout \t\t\t\t\t\t",
    "       Reflectivity to Acoustic Impedance (Inversion)\t\t",
    "\t\t\t\t\t\t\t        ",
    " Required parameters: (none)\t\t\t\t        ",
    "\t\t\t\t\t\t\t        ",
    " Optional parameters: (none)\t\t\t\t        ",
    "\t\t\t\t\t\t\t        ",
    " Notes:\t\t\t\t\t\t        ",
    " Reflectivity traces in acoustic impedence out\t\t        ",
    "\t\t\t\t\t\t\t        ",
];

fn main() {
    let argv = args();
    let mut r = Trace::new();
    let mut ai = Trace::new();
    let mut verbose = 0;

    initargs(&argv);
    requestdoc(1, SDOC);
    getparint("verbose", &mut verbose);

    let mut itr = gettr(&mut r);
    if itr == 0 {
        warn!("can't get first trace");
    }
    assert!(itr != 0);
    let nt = r.ns as usize;
    ai.hdr = r.hdr.clone();

    itr = 1;
    loop {
        ai.data[0] = 1.0;
        for it in 0..nt.saturating_sub(1) {
            if r.data[it] > 1.0 || r.data[it] < -1.0 {
                err!(" Illegal Reflectivity=%g ... abs value must be <1\n", r.data[it]);
            }
            ai.data[it + 1] = ((ai.data[it] as f64) * (1.0 - r.data[it] as f64)
                / (1.0 + r.data[it] as f64)) as f32;
        }
        puttr(&mut ai);
        itr += 1;
        if gettr(&mut r) == 0 {
            break;
        }
    }
    let _ = itr;
    exit(CWP_Exit());
}
