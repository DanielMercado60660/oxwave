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

//! SUNAN - remove NaNs and Infs from the input stream.
//! Port of su/main/amplitudes/sunan.c (Revision 1.8).
//! Credits: Reginald H. Beardsley (2003); interpolation idea from a
//! version modified by Balasz Nemeth at Potash Corporation in Saskatchewan.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUNAN - remove NaNs & Infs from the input stream\t\t",
    "\t\t\t\t\t\t\t\t",
    "    sunan < in.su >out.su\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " verbose=1\techo locations of NaNs or Infs to stderr\t",
    "\t        =0 silent\t\t\t\t\t",
    " ...user defined ... \t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " value=0.0\tNaNs and Inf replacement value\t\t\t",
    " ... and/or....\t\t\t\t\t\t",
    " interp=0\t=1 replace NaNs and Infs by interpolating\t",
    "                   neighboring finite values\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " A simple program to remove NaNs and Infs from an input stream.",
    " The program sets NaNs and Infs to \"value\" if interp=0. When\t",
    " interp=1 NaNs are replaced with the average of neighboring values",
    " provided that the neighboring values are finite, otherwise\t",
    " NaNs and Infs are replaced by \"value\".\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut itr: i32 = 0;
    let mut verbose = 1;
    let mut interp = 0;
    let mut value = 0.0f32;
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("Can't get first trace \n"); }
    getparint("verbose", &mut verbose);
    getparint("interp", &mut interp);
    getparfloat("value", &mut value);
    checkpars();
    loop {
        itr = itr.wrapping_add(1);
        for i in 0..tr.ns as usize {
            if !tr.data[i].is_finite() {
                if verbose != 0 { warn!("found NaN trace = %d  sample = %d", itr, i as i32); }
                // The C interpolation result is unconditionally overwritten below.
                // Avoid its out-of-bounds neighbor reads for short traces.
                let _ = interp;
                tr.data[i] = value;
            }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
