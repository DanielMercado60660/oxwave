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

//! VTLVZ - velocity as a function of time for linear V(Z).
//! Port of par/main/velocity_model_building/vtlvz.c (Revision 1.9).

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " VTLVZ -- Velocity as function of Time for Linear V(Z);\t\t",
    "          writes out a vector of velocity = v0 exp(a t/2)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " vtlvz > velfile nt= dt= v0= a= \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters\t\t\t\t\t\t\t",
    " nt=\tnumber of time samples\t\t\t\t\t\t",
    " dt=\ttime sampling interval\t\t\t\t\t\t",
    " v0=\tvelocity at the surface\t\t\t\t\t\t",
    " a=\tvelocity gradient\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut nt = 0i32;
    let mut dt = 0.0f32;
    let mut v0 = 0.0f32;
    let mut a = 0.0f32;
    initargs(&argv);
    requestdoc(1, SDOC);
    mustgetparint("nt", &mut nt);
    mustgetparfloat("dt", &mut dt);
    mustgetparfloat("v0", &mut v0);
    mustgetparfloat("a", &mut a);
    checkpars();
    if nt < 0 {
        err!("nt must be nonnegative");
    }
    let mut v = vec![0.0f32; nt as usize];
    let mut t = 0.0f32;
    for it in 0..nt {
        v[it as usize] = ((v0 as f64) * (((a * t) as f64) / 2.0).exp()) as f32;
        t += dt;
    }
    fwrite(&v, &stdout());
    exit(CWP_Exit());
}
