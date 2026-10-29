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

//! ENTROPY - compute the entropy of a signal.
//! Port of comp/dct/main/entropy.c (Revision 1.5).
//! Credits (from the C source): CWP: Tong Chen, 1995.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " ENTROPY - compute the ENTROPY of a signal\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "  entropy < stdin n= > stdout\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required Parameter:\t\t\t\t\t\t",
    "  n\t\tnumber of values in data set\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    "  none\t\t\t\t\t\t\t\t",
    "  \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut n = 0i32;
    mustgetparint("n", &mut n);
    // The C loops skip negative n; its failed allocation is never dereferenced
    // with empty input. Retain those signed loop bounds below.
    let mut f = vec![0.0f32; n.max(0) as usize];
    // The C ignores the read count. Its unread samples are uninitialized;
    // this initialized buffer supplies zeros if the input is short.
    fread(&mut f, &stdin());
    let mut hist = [0i32; 1024];
    let mut fhist = [0.0f32; 1024];
    let mut rdev = 0.0f32;
    for i in 0..n {
        let i = i as usize;
        rdev += f[i] * f[i];
    }
    rdev /= n as f32;
    rdev = (rdev as f64).sqrt() as f32;
    let mut dev = 0.0f32;
    if getparfloat("dev", &mut dev) == 0 {
        dev = rdev;
    }
    checkpars();
    fprintf!(stderr(), "dev=%f\n", dev);
    // C's unsuffixed constants make both scale expressions double precision.
    let step = ((dev as f64) * 3.464 * 0.01) as f32;
    let rstep = (1.0f64 / step as f64) as f32;
    let mut error = 0.0f32;
    for i in 0..n {
        let i = i as usize;
        let max = f[i] * rstep;
        // NINT returns int; the products and differences below stay float.
        let diff = (nint(max as f64) as f32) * step - f[i];
        error += diff * diff;
        let index = nint(max as f64).wrapping_add(512) as usize;
        hist[index] = hist[index].wrapping_add(1);
    }
    error /= n as f32;
    error = (error as f64).sqrt() as f32;
    error /= rdev;
    let mut ent = 0.0f32;
    for j in 0..1024 {
        fhist[j] = hist[j] as f32 / n as f32;
        if hist[j] != 0 {
            ent = (ent as f64 + fhist[j] as f64 * (fhist[j] as f64).ln() / 2.0f64.ln()) as f32;
        }
    }
    ent = -ent;
    fprintf!(
        stderr(),
        "entropy of the signal is=%f, average error=%f\n",
        ent,
        error
    );
    fwrite(&fhist, &stdout());
    exit(CWP_Exit());
}
