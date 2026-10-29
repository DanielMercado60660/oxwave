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
//! SUNULL - create null (all zeroes) traces.
//! Port of su/main/synthetics_waveforms_testpatterns/sunull.c (Revision 1.16).
//! Credits (from the C source): CWP: Jack K. Cohen.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUNULL - create null (all zeroes) traces\t \t\t",
    " \t\t\t\t\t\t\t\t",
    " sunull nt=   [optional parameters] >outdata\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameter\t\t\t\t\t\t",
    " \tnt=\t\tnumber of samples per trace\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters\t\t\t\t\t\t",
    " \tntr=5\t\tnumber of null traces to create\t\t",
    " \tdt=0.004\ttime sampling interval\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Rationale: It is sometimes useful to insert null traces\t",
    "\t between \"panels\" in a shell loop.\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " See also: sukill, sumute, suzero\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut nt = 0i32;
    let mut ntr = 0i32;
    let mut dt = 0.0f32;

    initargs(&argv);
    requestdoc(0, SDOC);

    mustgetparint("nt", &mut nt);
    CHECK_NT("nt", nt);
    tr.ns = nt as u16;

    if getparint("ntr", &mut ntr) == 0 { ntr = 5; }
    if getparfloat("dt", &mut dt) == 0 { dt = 0.004; }
    let idt = nint(1000000.0 * dt as f64);

    // Trace::new initializes the data to zero, as the original memset does.
    tr.dt = idt as u16;
    tr.ns = nt as u16;
    for itr in 0..ntr {
        tr.tracl = itr + 1;
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
