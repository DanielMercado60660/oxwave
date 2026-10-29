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

//! XY2Z - convert coordinate pairs to spikes on a grid.
//! Port of par/main/resampling/xy2z.c (Revision 1.8).
//! Credits: CWP: John Stockwell, Nov 1995.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t \t\t\t\t\t\t\t\t",
    " XY2Z - converts (X,Y)-pairs to spike Z values on a uniform grid\t",
    " \t      \t\t\t\t\t\t\t\t",
    "    xy2z < stdin npairs= [optional parameters] >stdout \t\t",
    " \t      \t\t\t\t\t\t\t\t",
    " Required parameter:\t\t\t\t\t\t\t",
    " npairs= \tnumber of pairs input\t\t\t\t\t",
    " \t      \t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t\t",
    " scale=1.0\tvalue to scale spikes by\t\t\t\t",
    " nx1=100 \tdimension of first (fast) dimension of output array\t",
    " nx2=100 \tdimension of second (slow) dimension of output array\t",
    " x1pad=2\tzero padding in x1 dimension\t\t\t\t",
    " x2pad=2\tzero padding in x2 dimension\t\t\t\t",
    " yx=0\t\tassume (x,y) pairs \t\t\t\t\t",
    " \t\t=1\tassume (y,x) pairs \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " Converts ordered (x,y) pairs to spike x1values, of height=scale on a \t",
    " uniform grid.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut npairs = 0;
    mustgetparint("npairs", &mut npairs);
    let mut scale = 1.0f32;
    let mut nx1 = 100;
    let mut nx2 = 100;
    let mut x1pad = 2;
    let mut x2pad = 2;
    let mut yx = 0;
    getparfloat("scale", &mut scale);
    getparint("nx1", &mut nx1);
    getparint("nx2", &mut nx2);
    getparint("x1pad", &mut x1pad);
    getparint("x2pad", &mut x2pad);
    getparint("yx", &mut yx);
    checkpars();
    if x1pad >= nx1 { err!("x1pad=%d is too large!", x1pad); }
    if x2pad >= nx2 { err!("x2pad=%d is too large!", x2pad); }
    let mut x1 = vec![0.0f32; npairs as usize];
    let mut x2 = vec![0.0f32; npairs as usize];
    let mut outarray = vec![0.0f32; (nx1 * nx2) as usize];
    let mut z = [0.0f32; 2];
    for i in 0..npairs as usize {
        efread(&mut z, &stdin());
        if yx != 0 {
            x1[i] = z[0]; x2[i] = z[1];
        } else {
            x2[i] = z[0]; x1[i] = z[1];
        }
    }
    let mut x1min = x1[0]; let mut x1max = x1[0];
    let mut x2min = x2[0]; let mut x2max = x2[0];
    for i in 0..npairs as usize {
        x1min = if x1min < x1[i] { x1min } else { x1[i] };
        x1max = if x1max > x1[i] { x1max } else { x1[i] };
        x2min = if x2min < x2[i] { x2min } else { x2[i] };
        x2max = if x2max > x2[i] { x2max } else { x2[i] };
    }
    for i in 0..npairs as usize {
        x1[i] = (nx1 - 2*x1pad) as f32 * (x1[i] - x1min) / (x1max - x1min) + x1pad as f32;
        x2[i] = (nx2 - 2*x2pad) as f32 * (x2[i] - x2min) / (x2max - x2min) + x2pad as f32;
    }
    for i in 0..npairs as usize {
        let j1 = nint(x1[i] as f64);
        let j2 = nint(x2[i] as f64);
        outarray[(j2 * nx1 + j1) as usize] = scale;
    }
    efwrite(&outarray, &stdout());
    exit(CWP_Exit());
}
