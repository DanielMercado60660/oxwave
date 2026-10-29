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

//! RESAMP - resample the first dimension of a two-dimensional function.
//! Port of par/main/resampling/resamp.c (Revision 1.11).
//! Credits: Dave Hale, Colorado School of Mines, 07/07/89.

use oxwave_core::prelude::*;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " RESAMP - RESAMPle the 1st dimension of a 2-dimensional function f(x1,x2)",
    " \t\t\t\t\t\t\t\t\t",
    " resamp <infile >outfile [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " n1=all                 number of samples in 1st (fast) dimension\t",
    " n2=all                 number of samples in 2nd (slow) dimension\t",
    " d1=1.0                 sampling interval in 1st dimension\t\t",
    " f1=d1                  first sample in 1st dimension\t\t\t",
    " n1r=n1                 number of samples in 1st dimension after resampling",
    " d1r=d1                 sampling interval in 1st dimension after resampling",
    " f1r=f1                 first sample in 1st dimension after resampling\t",
    " \t\t\t\t\t\t\t\t\t",
    " NOTE:  resamp currently performs NO ANTI-ALIAS FILTERING before resampling!",
    " Caveat: this program resamples data that are oscillatory in the fast\t",
    "    dimension only, such as seismic data with no SU headers. To resample",
    "    other 2d data, such as velocity profiles, use \"unisam\" or \"unisam2",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let infp = stdin();
    let outfp = stdout();
    let mut n1 = 0i32;
    let mut n2 = -1i32;
    let mut d1 = 1.0f32;
    let mut f1 = 0.0f32;
    let mut n1r = 0i32;
    let mut d1r = 0.0f32;
    let mut f1r = 0.0f32;
    if getparint("n1", &mut n1) == 0 {
        if efseeko(&infp, 0, SEEK_END) == -1 { err!("must specify n1!\n"); }
        n1 = (eftello(&infp) / 4) as i32;
        efseeko(&infp, 0, SEEK_SET);
    }
    getparint("n2", &mut n2);
    getparfloat("d1", &mut d1);
    if getparfloat("f1", &mut f1) == 0 { f1 = d1; }
    if getparint("n1r", &mut n1r) == 0 { n1r = n1; }
    if getparfloat("d1r", &mut d1r) == 0 { d1r = d1; }
    if getparfloat("f1r", &mut f1r) == 0 { f1r = f1; }
    checkpars();
    let mut p = vec![0.0f32; n1 as usize];
    let mut pr = vec![0.0f32; n1r as usize];
    let mut xr = vec![0.0f32; n1r as usize];
    let mut x1r = f1r;
    for x in &mut xr { *x = x1r; x1r += d1r; }
    let mut i2 = 0i32;
    while i2 < n2 || n2 < 0 {
        if efread(&mut p, &infp) != n1 as usize { break; }
        ints8r(d1, f1, &p, 0.0, 0.0, &xr, &mut pr);
        efwrite(&pr, &outfp);
        i2 = i2.wrapping_add(1);
    }
    exit(CWP_Exit());
}
