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

//! SUSTKVEL - convert interval velocities to stacking velocities.
//! Port of su/main/data_conversion/sustkvel.c (Revision 1.14).
//! Credits (from the C source): CWP: Jack

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSTKVEL - convert constant dip layer interval velocity model to the\t",
    "\t   stacking velocity model required by sunmo\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sustkvel v= h= dip=0.0 outpar=/dev/tty\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t        \t",
    "\tv=\tinterval velocities \t\t\t\t\t",
    "\th=\tlayer thicknesses at the cmp\t \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tdip=0.0\t\t\t(constant) dip of the layers (degrees)\t",
    "\toutpar=/dev/tty\t\toutput parameter file in the form\t",
    "\t\t\t\trequired by sunmo:\t\t\t",
    "\t\t\t\ttv=zero incidence time pick vector\t",
    "\t\t\t\tv=stacking velocities vector\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "    sustkvel v=5000,6000,8000,10000 h=1000,1200,1300,1500 outpar=stkpar",
    "    sunmo <data.cdp par=stkpar >data.nmo\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "    sustkvel par=intpar outpar=stkpar\t\t\t\t\t",
    "    sunmo <data.cdp par=stkpar >data.nmo\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " If the file, intpar, contains:\t\t\t\t\t",
    "    v=5000,6000,8000,10000\t\t\t\t\t\t",
    "    h=1000,1200,1300,1500\t\t\t\t\t\t",
    " then the two examples are equivalent.  The created parameter file,\t",
    " stkpar, is in the form of the velocity model required by sunmo.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: sustkvel does not have standard su syntax since it does not\t",
    "      operate on seismic data.  Hence stdin and stdout are not used.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat: Does not accept a series of interval velocity models to\t",
    "\tproduce a variable velocity file for sunmo.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut dip = 0.0f32;
    if getparfloat("dip", &mut dip) == 0 { dip = 0.0; }
    dip = (dip as f64 * PI / 180.0) as f32;
    let cosdip = (dip as f64).cos() as f32;
    let ocosdip = (1.0 / cosdip as f64) as f32;
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let outparfp = efopen(&outpar, "w");
    let n = countparval("v");
    if n == 0 { err!("no v's specified"); }
    let mut v = vec![0.0f32; n as usize];
    let mut h = vec![0.0f32; n as usize];
    let mut vs = vec![0.0f32; n as usize];
    let mut t0 = vec![0.0f32; n as usize];
    if n != getparfloat("h", &mut h) { err!("expected %d intervals", n); }
    if n != getparfloat("v", &mut v) { err!("expected %d velocities", n); }
    checkpars();
    for i in 0..n as usize {
        if v[i] <= 0.0 { err!("v's must be positive: v[%d] = %f", i as i32, v[i]); }
        if h[i] <= 0.0 { err!("h's must be positive: h[%d] = %f", i as i32, h[i]); }
    }
    let mut t0by2 = 0.0f32;
    let mut vh = 0.0f32;
    for i in 0..n as usize {
        t0by2 += h[i] / v[i];
        vh += v[i] * h[i];
        t0[i] = (2.0 * t0by2 as f64 * cosdip as f64) as f32;
        vs[i] = (((vh / t0by2) as f64).sqrt() * ocosdip as f64) as f32;
    }
    fprintf!(outparfp, "tnmo=");
    for i in 0..(n as usize - 1) { fprintf!(outparfp, "%g,", t0[i]); }
    fprintf!(outparfp, "%g\n", t0[n as usize - 1]);
    fprintf!(outparfp, "vnmo=");
    for i in 0..(n as usize - 1) { fprintf!(outparfp, "%g,", vs[i]); }
    fprintf!(outparfp, "%g\n", vs[n as usize - 1]);
    exit(CWP_Exit());
}
