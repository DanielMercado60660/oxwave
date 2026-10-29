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

//! SUAI2R - acoustic impedance to reflectivity.
//! Port of su/main/amplitudes/suai2r.c.
//! Credits: Shuki Ronen, 2017. Partial implementation of Aki & Richards,
//! volume II, pages 661-662; transmission effects are ignored.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " suai2r <stdin >sdout \t\t\t\t\t\t",
    "\tAcoustic Impedance to Reflectivity (Forward Modeling)   ",
    "\t\t\t\t\t\t\t        ",
    "\t\t\t\t\t\t\t        ",
    " Required parameters: (none)\t\t\t\t        ",
    "\t\t\t\t\t\t\t        ",
    " Optional parameters: (none)\t\t\t\t        ",
    "\t\t\t\t\t\t\t        ",
    " Notes:\t\t\t\t\t\t        ",
    "  Acoustic impedence traces in Reflectivity out \t\t",
    "\t\t\t\t\t\t\t        ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut itr = gettr(&mut tr);
    if itr == 0 { warn!("can't get first trace"); }
    assert!(itr != 0);
    let nt = tr.ns as usize;
    itr = 1;
    loop {
        for it in 0..nt - 1 {
            tr.data[it] = -(tr.data[it + 1] - tr.data[it]) / (tr.data[it + 1] + tr.data[it]);
        }
        tr.data[nt - 1] = 0.0;
        puttr(&mut tr);
        itr += 1;
        if gettr(&mut tr) == 0 { break; }
    }
    let _ = itr;
    exit(CWP_Exit());
}
