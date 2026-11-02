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

//! SULPRIME - find an appropriate Backus averaging length.
//! Port of su/main/well_logs/sulprime.c (Revision 1.4).
//! Credits: UHouston: Chris Liner, September 2008. He acknowledged Saudi
//! Aramco for permission to release code developed at EXPEC-ARC.
//! Reference: The Backus Number (Liner and Fei, 2007, TLE).

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SULPRIME - find appropriate Backus average length for  \t",
    " \t\ta given log suite, frequency, and purpose\t\t",
    " \t\t\t\t\t\t\t\t",
    " sulprime < vp_vs_rho.su  [options]\t\tor\t\t",
    " sulprime < vp_vs_rho_gr.su   [options]\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    " b=2.0\t\ttarget value of Backus number\t\t \t",
    "\t\tb=2 is transmission limit (ok for proc, mig, etc.)",
    "\t\tb=0.3 is scattering limit (ok for modeling)\t",
    " dz=1\t\tinput depth sample interval (ft)\t\t ",
    " f=60\t\tfrequency (Hz)... dominant or max (to be safe) \t",
    " nmax=301\tmaximum averaging length (samples)\t\t",
    " verbose=1\tprint intermediate results\t\t\t",
    "\t\t=0 print final result only\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " 1. input is in sync with subackus, but vp and gr not used\t",
    "     (gr= gamma ray log)\t\t\t\t\t",
    " Related codes:  subackus, subackush\t\t\t\t",
];

fn dobackus(navg: i32, nz: usize, p: &[f32], avg: &mut [f32]) {
    let m = (navg - 1) / 2;
    for i in 0..nz {
        let mut val = 0.0f32;
        let ii = i as i32;
        if ii - m > 0 && ii + m < nz as i32 {
            for j in (ii - m)..=(ii + m) { val += p[j as usize]; }
            avg[i] = val / navg as f32;
        }
        if ii <= m { avg[i] = p[i]; }
        if ii >= nz as i32 - m { avg[i] = p[i]; }
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut btarg = 2.0f32;
    let mut dz = 1.0f32;
    let mut freq = 60.0f32;
    let mut nmax = 301i32;
    let mut verbose = 0i32;
    getparfloat("b", &mut btarg);
    getparfloat("dz", &mut dz);
    getparfloat("f", &mut freq);
    getparint("nmax", &mut nmax);
    getparint("verbose", &mut verbose);
    checkpars();
    warn!("Target value: B=%g  Frequency: %g Hz", btarg, freq);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nz = tr.ns as usize;
    let mut vs = vec![0.0f32; nz];
    let mut rhoo = vec![0.0f32; nz];
    gettr(&mut tr);
    for i in 0..nz { vs[i] = tr.data[i]; }
    gettr(&mut tr);
    for i in 0..nz { rhoo[i] = tr.data[i]; }
    let mut lo = vec![0.0f32; nz];
    let mut lb = vec![0.0f32; nz];
    let mut rhob = vec![0.0f32; nz];
    let mut tmp = vec![0.0f32; nz];
    let mut vs0 = vec![0.0f32; nz];
    for i in 0..nz { lo[i] = rhoo[i] * vs[i] * vs[i]; }
    let mut errmin = 999.0f32;
    let mut bmin = 0.0f32;
    let mut nmin = 0i32;
    let mut vmin = 0.0f32;
    let mut navg = 1i32;
    while navg <= nmax {
        for i in 0..nz { tmp[i] = (1.0f64 / lo[i] as f64) as f32; }
        dobackus(navg, nz, &tmp, &mut lb);
        for i in 0..nz { lb[i] = (1.0f64 / lb[i] as f64) as f32; }
        dobackus(navg, nz, &rhoo, &mut rhob);
        let mut minvs0 = 999999.0f32;
        for i in 0..nz {
            vs0[i] = ((lb[i] / rhob[i]) as f64).sqrt() as f32;
            if vs0[i] < minvs0 { minvs0 = vs0[i]; }
        }
        let bb = ((freq * navg as f32 * dz) as f64 * 0.3048 / minvs0 as f64) as f32;
        let errf = ((btarg - bb) as f64).abs() as f32;
        if errf < errmin {
            errmin = errf;
            nmin = navg;
            bmin = bb;
            vmin = minvs0;
        }
        if verbose == 1 {
            warn!("err=%g; B = %g; L=%i (samples); minvs0=%g", errf, bb, navg, minvs0);
        }
        navg = navg.wrapping_add(2);
    }
    warn!("Results:\tB=%g   navg=%i  min[vs0]=%g", bmin, nmin, vmin);
    exit(CWP_Exit());
}
