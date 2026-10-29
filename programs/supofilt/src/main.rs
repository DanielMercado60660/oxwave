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

//! SUPOFILT - polarization filter for three-component data.
//! Port of su/main/multicomponent/supofilt.c (Revision 1.5).
//! Credits: Nils Maercklin, GeoForschungsZentrum Potsdam, 1999-2000.
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUPOFILT - POlarization FILTer for three-component data               ",
    "                                                                       ",
    " supofilt <stdin >stdout [optional parameters]                         ",
    "                                                                       ",
    " Required parameters:                                                  ",
    "    dfile=polar.dir   file containing the 3 components of the          ",
    "                      direction of polarization                        ",
    "    wfile=polar.rl    file name of weighting polarization parameter    ",
    "                                                                       ",
    " Optional parameters:                                                  ",
    "    dt=(from header)  time sampling intervall in seconds               ",
    "    smooth=1          1 = smooth filter operators, 0 do not            ",
    "    sl=0.05           smoothing window length in seconds               ",
    "    wpow=1.0          raise weighting function to power wpow           ",
    "    dpow=1.0          raise directivity functions to power dpow        ",
    "    verbose=0         1 = echo additional information                  ",
    "                                                                       ",
    "                                                                       ",
    " Notes:                                                                ",
    "    Three adjacent traces are considered as one three-component        ",
    "    dataset.                                                           ",
    "                                                                       ",
    "    This program SUPOFILT is an extension to the polarization analysis ",
    "    program supolar. The files wfile and dfile are SU files as written ",
    "    by SUPOLAR.                                                        ",
    "                                                                       ",
];

fn do_smooth(data: &mut [f32], nt: usize, isl: i32) {
    let mut tmpdata = vec![0.0f32; nt];
    for it in 0..nt as i32 {
        let mut sval = 0.0f32;
        if it >= isl / 2 && it < nt as i32 - isl / 2 {
            for jt in it - isl / 2..it + isl / 2 { sval += data[jt as usize]; }
            tmpdata[it as usize] = sval / isl as f32;
        }
    }
    data[..nt].copy_from_slice(&tmpdata);
}
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut dtr = Trace::new();
    let mut wtr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut wfile = String::new();
    let mut dfile = String::new();
    let mut sl = 0.0f32;
    let mut dt = 0.0f32;
    let mut wpow = 0.0f32;
    let mut dpow = 0.0f32;
    let mut smooth = 0i32;
    let mut verbose = 0i32;
    if getparstring("wfile", &mut wfile) == 0 { wfile = "polar.rl".into(); }
    if getparstring("dfile", &mut dfile) == 0 { dfile = "polar.dir".into(); }
    if getparfloat("sl", &mut sl) == 0 { sl = 0.05; }
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if getparfloat("wpow", &mut wpow) == 0 { wpow = 1.0; }
    if getparfloat("dpow", &mut dpow) == 0 { dpow = 1.0; }
    if getparint("smooth", &mut smooth) == 0 { smooth = 1; }
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    checkpars();
    if dt == 0.0 { dt = 0.004; warn!("dt not set, assuming dt=0.004"); }
    let isl = nint((sl / dt) as f64);
    if verbose != 0 && smooth != 0 { warn!("smoothing window = %d samples", isl); }
    let dfp = fopen(&dfile, "r").unwrap_or_else(|| err!("file %s does not exist", dfile));
    let wfp = fopen(&wfile, "r").unwrap_or_else(|| err!("file %s does not exist", wfile));
    if fgettr(&dfp, &mut dtr) == 0 { err!("can't get first trace of %s", dfile); }
    if fgettr(&wfp, &mut wtr) == 0 { err!("can't get first trace of %s", wfile); }
    if smooth != 0 && isl < 1 { err!("sl=%g must be positive", sl); }
    if smooth != 0 && isl > nt as i32 { err!("sl=%g too long for trace", sl); }
    if nt != dtr.ns as usize { err!("stdin and %s have different ns (%d vs %d)", dfile, nt as i32, dtr.ns as i32); }
    if nt != wtr.ns as usize { err!("stdin and %s have different ns (%d vs %d)", wfile, nt as i32, wtr.ns as i32); }
    let mut icomp = 0;
    let mut itr = 1;
    loop {
        for i in 0..nt { dtr.data[i] = (dtr.data[i] as f64).abs().powf(dpow as f64) as f32; }
        if smooth != 0 { do_smooth(&mut dtr.data, nt, isl); }
        if icomp == 0 {
            for i in 0..nt { wtr.data[i] = (wtr.data[i] as f64).powf(wpow as f64) as f32; }
            if smooth != 0 { do_smooth(&mut wtr.data, nt, isl); }
        }
        for i in 0..nt { tr.data[i] *= dtr.data[i] * wtr.data[i]; }
        puttr(&mut tr);
        icomp += 1;
        if icomp == 3 {
            icomp = 0;
            if fgettr(&wfp, &mut wtr) == 0 { break; }
        }
        itr += 1;
        if gettr(&mut tr) == 0 || fgettr(&dfp, &mut dtr) == 0 { break; }
    }
    if verbose != 0 { warn!("processed %d traces", itr); }
    exit(CWP_Exit());
}
