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

//! UNISAM2 - uniformly sample a two-dimensional function.
//! Port of par/main/resampling/unisam2.c (Revision 1.12).
//! Author: Dave Hale, Colorado School of Mines, 01/12/91.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_f_cubicspline::{cakima, cmonot, csplin};
use oxwave_linear_interp::intlin;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " UNISAM2 - UNIformly SAMple a 2-D function f(x1,x2)\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " unisam2 [optional parameters] <inputfile >outputfile\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " x1=             array of x1 values at which input f(x1,x2) is sampled\t",
    " ... Or specify a unform linear set of values for x1 via:\t\t",
    " nx1=1           number of input samples in 1st dimension\t\t",
    " dx1=1           input sampling interval in 1st dimension\t\t",
    " fx1=0           first input sample in 1st dimension\t\t\t",
    " ...\t\t\t\t\t\t\t\t\t",
    " n1=1            number of output samples in 1st dimension\t\t",
    " d1=             output sampling interval in 1st dimension\t\t",
    " f1=             first output sample in 1st dimension\t\t\t",
    " x2=             array of x2 values at which input f(x1,x2) is sampled\t",
    " ... Or specify a unform linear set of values for x2 via:\t\t",
    " nx2=1           number of input samples in 2nd dimension\t\t",
    " dx2=1           input sampling interval in 2nd dimension\t\t",
    " fx2=0           first input sample in 2nd dimension\t\t\t",
    " ...\t\t\t\t\t\t\t\t\t",
    " n2=1            number of output samples in 2nd dimension\t\t",
    " d2=             output sampling interval in 2nd dimension\t\t",
    " f2=             first output sample in 2nd dimension\t\t\t",
    " ... \t\t\t\t\t\t\t\t\t",
    " method1=linear  =linear for linear interpolation\t\t\t",
    "                 =mono for monotonic bicubic interpolation\t\t",
    "                 =akima for Akima bicubic interpolation\t\t",
    "                 =spline for bicubic spline interpolation\t\t",
    " method2=linear  =linear for linear interpolation\t\t\t",
    "                 =mono for monotonic bicubic interpolation\t\t",
    "                 =akima for Akima bicubic interpolation\t\t",
    "                 =spline for bicubic spline interpolation\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " NOTES:\t\t\t\t\t\t\t\t",
    " The number of input samples is the number of x1 values times the\t",
    " number of x2 values.  The number of output samples is n1 times n2.\t",
    " The output sampling intervals (d1 and d2) and first samples (f1 and f2)",
    " default to span the range of input x1 and x2 values.  In other words,\t",
    " d1=(x1max-x1min)/(n1-1) and f1=x1min; likewise for d2 and f2.\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Interpolation is first performed along the 2nd dimension for each\t",
    " value of x1 specified.  Interpolation is then performed along the\t",
    " 1st dimension.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut nx1 = countparval("x1");
    let mut x1;
    if nx1 != 0 {
        x1 = vec![0.0; nx1 as usize];
        getparfloat("x1", &mut x1);
    } else {
        nx1 = 1;
        let mut dx1 = 1.0f32;
        let mut fx1 = 0.0f32;
        getparint("nx1", &mut nx1);
        getparfloat("dx1", &mut dx1);
        getparfloat("fx1", &mut fx1);
        if nx1 < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
        x1 = vec![0.0; nx1 as usize];
        for ix1 in 0..nx1 { x1[ix1 as usize] = fx1 + ix1 as f32 * dx1; }
    }
    let mut x1min = x1[0];
    let mut x1max = x1[0];
    for ix1 in 1..nx1 as usize {
        if x1[ix1] < x1min { x1min = x1[ix1]; }
        if x1[ix1] > x1max { x1max = x1[ix1]; }
    }
    let mut n1 = 1;
    getparint("n1", &mut n1);
    let mut d1 = if n1 == 1 { 0.0 } else { (x1max-x1min)/(n1-1) as f32 };
    let mut f1 = x1min;
    getparfloat("d1", &mut d1);
    getparfloat("f1", &mut f1);

    let mut nx2 = countparval("x2");
    let mut x2;
    if nx2 != 0 {
        x2 = vec![0.0; nx2 as usize];
        getparfloat("x2", &mut x2);
    } else {
        nx2 = 1;
        let mut dx2 = 1.0f32;
        let mut fx2 = 0.0f32;
        getparint("nx2", &mut nx2);
        getparfloat("dx2", &mut dx2);
        getparfloat("fx2", &mut fx2);
        if nx2 < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
        x2 = vec![0.0; nx2 as usize];
        for ix2 in 0..nx2 { x2[ix2 as usize] = fx2 + ix2 as f32 * dx2; }
    }
    let mut x2min = x2[0];
    let mut x2max = x2[0];
    for ix2 in 1..nx2 as usize {
        if x2[ix2] < x2min { x2min = x2[ix2]; }
        if x2[ix2] > x2max { x2max = x2[ix2]; }
    }
    let mut n2 = 1;
    getparint("n2", &mut n2);
    let mut d2 = if n2 == 1 { 0.0 } else { (x2max-x2min)/(n2-1) as f32 };
    let mut f2 = x2min;
    getparfloat("d2", &mut d2);
    getparfloat("f2", &mut f2);
    let mut method1 = String::from("linear");
    let mut method2 = String::from("linear");
    getparstring("method1", &mut method1);
    getparstring("method2", &mut method2);
    checkpars();
    let m1 = method1.as_bytes().first().copied().unwrap_or(0);
    let m2 = method2.as_bytes().first().copied().unwrap_or(0);
    if !b"lmas".contains(&m1) { err!("%s is an unknown interpolation method!\n", method1); }
    if !b"lmas".contains(&m2) { err!("%s is an unknown interpolation method!\n", method2); }

    // A negative C int becomes an enormous size_t in ealloc1float, so malloc fails.
    if n1 < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut u1 = vec![0.0; n1 as usize];
    for i1 in 0..n1 { u1[i1 as usize] = f1 + i1 as f32 * d1; }
    if n2 < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut u2 = vec![0.0; n2 as usize];
    for i2 in 0..n2 { u2[i2 as usize] = f2 + i2 as f32 * d2; }
    let (nx1, nx2, n1, n2) = (nx1 as usize, nx2 as usize, n1 as usize, n2 as usize);
    let mut g = vec![0.0f32; nx1 * nx2];
    let mut h = vec![vec![0.0f32; nx1]; n2];
    if fread(&mut g, &stdin()) != nx1 * nx2 { err!("error reading input file!"); }
    let mut g2 = vec![0.0f32; nx2];
    let mut h2 = vec![0.0f32; n2];
    for ix1 in 0..nx1 {
        for ix2 in 0..nx2 { g2[ix2] = g[ix2 * nx1 + ix1]; }
        if m1 == b'l' || nx2 == 1 {
            intlin(&x2, &g2, g2[0], g2[nx2-1], &u2, &mut h2);
        } else {
            let mut g2d = vec![[0.0; 4]; nx2];
            match m1 { b'm' => cmonot(&x2, &g2, &mut g2d), b'a' => cakima(&x2, &g2, &mut g2d), b's' => csplin(&x2, &g2, &mut g2d), _ => {} }
            intcub(0, &x2, &g2d, &u2, &mut h2);
        }
        for i2 in 0..n2 { h[i2][ix1] = h2[i2]; }
    }
    g = vec![0.0; n1 * n2];
    for i2 in 0..n2 {
        let output = &mut g[i2*n1..(i2+1)*n1];
        if m2 == b'l' || nx1 == 1 {
            intlin(&x1, &h[i2], h[i2][0], h[i2][nx1-1], &u1, output);
        } else {
            let mut h1d = vec![[0.0; 4]; nx1];
            match m2 { b'm' => cmonot(&x1, &h[i2], &mut h1d), b'a' => cakima(&x1, &h[i2], &mut h1d), b's' => csplin(&x1, &h[i2], &mut h1d), _ => {} }
            intcub(0, &x1, &h1d, &u1, output);
        }
    }
    if fwrite(&g, &stdout()) != n1*n2 { err!("error writing output file!"); }
    exit(CWP_Exit());
}
