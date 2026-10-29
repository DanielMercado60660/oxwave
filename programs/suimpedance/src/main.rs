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


//! SUIMPEDANCE - convert reflection coefficients to impedances.
//! Port of su/main/synthetics_waveforms_testpatterns/suimpedance.c (Revision 1.4).
//! Credits: SEP: Stew Levin.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUIMPEDANCE - Convert reflection coefficients to impedances.  ",
    "\t\t\t\t\t\t\t\t",
    " suimpedance <stdin >stdout [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t  \t",
    " v0=1500.\tVelocity at first sample (m/sec)\t\t",
    " rho0=1.0e6\tDensity at first sample  (g/m^3)\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " Implements recursion [1-R(k)]Z(k) = [1+R(k)]Z(k-1).\t\t",
    " The input traces are assumed to be reflectivities, and thus are",
    " expected to have amplitude values between -1.0 and 1.0.\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut v0: f32 = 0.0;
    let mut rho0: f32 = 0.0;

    initargs(&argv);
    requestdoc(0, SDOC);

    if getparfloat("v0", &mut v0) == 0 { v0 = 1500.0; }
    if getparfloat("rho0", &mut rho0) == 0 { rho0 = 1000000.0; }

    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    loop {
        rctoimp(v0, rho0, nt, &mut tr.data);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}

// Convert reflection coefficient trace to impedance trace.
fn rctoimp(v0: f32, rho0: f32, nt: usize, trace: &mut [f32]) {
    let mut zold = (v0 * rho0) as f64;
    for it in 0..nt {
        let mut rc = trace[it] as f64;
        if rc <= -1.0 { rc = -0.9999; }
        if rc >= 1.0 { rc = 0.9999; }
        trace[it] = zold as f32;
        zold *= (1.0 + rc) / (1.0 - rc);
    }
}
