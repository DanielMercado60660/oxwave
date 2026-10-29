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
//! SUMEAN - get the mean values of data traces.
//! Port of su/main/attributes_parameter_estimation/sumean.c (Revision 1.5).
//! Credits: Bjoern E. Rommel, IKU, Petroleumsforskning / October 1997.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMEAN - get the mean values of data traces\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sumean < stdin > stdout [optional parameters] \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "   power = 2.0\t\tmean to the power\t\t\t\t",
    "\t\t\t(e.g. = 1.0 mean amplitude, = 2.0 mean energy)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters: \t\t\t\t\t\t\t",
    "   verbose = 0\t\twrites mean value of section to outpar\t   \t",
    "\t\t\t= 1 writes mean value of each trace / section to",
    "\t\t\t\toutpar\t\t\t\t\t",
    "   outpar=/dev/tty   output parameter file\t\t\t\t",
    "   abs = 1             average absolute value ",
    "                       = 0 preserve sign if power=1.0",
    "\t\t\t \t\t\t\t\t\t",
    " Notes:\t\t\t \t\t\t\t\t",
    " Each sample is raised to the requested power, and the sum of all those",
    " values is averaged for each trace (verbose=1) and the section.\t",
    " The values power=1.0 and power=2.0 are physical, however other powers\t",
    " represent other mathematical L-p norms and may be of use, as well.\t",
    "\t\t\t \t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut verbose = 0i32;
    let mut outpar = String::new();
    let mut abs = 1i32;
    let mut power = 0.0f32;
    initargs(&argv);
    requestdoc(1, SDOC);
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    if verbose < 0 || verbose > 1 { err!("unknown mode for verbose!"); }
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let outparfp = efopen(&outpar, "w");
    if getparint("abs", &mut abs) == 0 { abs = 0; }
    if getparfloat("power", &mut power) == 0 { power = 2.0; }
    checkpars();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let mut gmean = 0.0f32;
    let mut itr = 0i32;
    loop {
        itr = itr.wrapping_add(1);
        let mut tmean = 0.0f32;
        for it in 0..tr.ns as usize {
            let sample = tr.data[it] as f64;
            let value = if abs != 0 { sample.abs() } else { sample };
            tmean = (tmean as f64 + value.powf(power as f64)) as f32;
        }
        tmean = ((tmean / tr.ns as f32) as f64).powf(1.0 / power as f64) as f32;
        gmean += tmean;
        if verbose == 1 { fprintf!(outparfp, "trace: %i   mean value: %e\n", itr, tmean); }
        if gettr(&mut tr) == 0 { break; }
    }
    gmean /= itr as f32;
    fprintf!(outparfp, "global mean: %e\n", gmean);
    exit(CWP_Exit());
}
