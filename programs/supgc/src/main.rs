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

//! SUPGC - programmed gain control. Port of su/main/amplitudes/supgc.c.
//! Author: John Anderson (visitor to CWP from Mobil).

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUPGC   -   Programmed Gain Control--apply agc like function\t",
    "              but the same function to all traces preserving\t",
    "              relative amplitudes spatially.\t\t\t",
    " Required parameter:\t\t\t\t\t\t",
    " file=             name of input file\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " ntrscan=200       number of traces to scan for gain function\t",
    " lwindow=1.0       length of time window in seconds\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut ntrscan = 200;
    getparint("ntrscan", &mut ntrscan);
    let mut file = String::new();
    getparstring("file", &mut file);
    let mut lwindow = 1.0f32;
    getparfloat("lwindow", &mut lwindow);
    checkpars();

    let fp = match fopen(&file, "r") {
        Some(fp) => fp,
        None => err!("could not open input file"),
    };
    if fgettr(&fp, &mut tr) == 0 { err!("Can't get first trace"); }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let lw = nint(0.5 * lwindow as f64 / dt as f64);
    warn!("window length= %d samples", lw);

    let mut gain = vec![0.0f32; nt as usize];
    let mut g = vec![0.0f32; nt as usize];
    let mut icount = 0i32;
    loop {
        for j in 0..nt as usize {
            gain[j] = (gain[j] as f64 + (tr.data[j] as f64).abs()) as f32;
        }
        icount += 1;
        if fgettr(&fp, &mut tr) == 0 || icount >= ntrscan { break; }
    }
    erewind(&fp);

    for j in 0..nt {
        let mut sum = 0.0f32;
        let mut kk = 0i32;
        let mut k = j.wrapping_sub(lw).max(0);
        while k < j.wrapping_add(lw).min(nt) {
            kk += 1;
            sum += gain[k as usize];
            k += 1;
        }
        if sum == 0.0 { sum = 1.0; }
        g[j as usize] = kk.wrapping_mul(icount) as f32 / sum;
    }
    warn!("scan covered %d traces", icount);
    if fgettr(&fp, &mut tr) == 0 { err!("failure after rewind"); }
    loop {
        for j in 0..nt as usize { tr.data[j] *= g[j]; }
        puttr(&mut tr);
        if fgettr(&fp, &mut tr) == 0 { break; }
    }
    fclose(&fp);
    exit(CWP_Exit());
}
