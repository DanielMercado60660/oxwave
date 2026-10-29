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
//! SU3DCHART - plot x-midpoints vs. y-midpoints for 3-D data.
//! Port of su/main/headers/su3dchart.c (Revision 1.15).
//! Credits: CWP: Shuki Ronen; Toralf Foerster.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SU3DCHART - plot x-midpoints vs. y-midpoints for 3-D data\t",
    " \t\t\t\t\t\t\t\t",
    " su3dchart <stdin >stdout\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " outpar=null\tname of parameter file\t\t\t\t",
    " degree=0\t=1 convert seconds of arc to degrees\t\t",
    " \t\t\t\t\t\t\t\t",
    " The output is the (x, y) pairs of binary floats\t\t",
    " \t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t",
    " su3dchart <segy_data outpar=pfile >plot_data\t\t\t",
    " psgraph <plot_data par=pfile \\\t\t\t\t",
    "\tlinewidth=0 marksize=2 mark=8 | ...\t\t\t",
    " rm plot_data \t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " su3dchart <segy_data | psgraph n=1024 d1=.004 \\\t\t",
    "\tlinewidth=0 marksize=2 mark=8 | ...\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut degree = 0;
    if getparint("degree", &mut degree) == 0 { degree = 0; }
    let mut npairs: i32 = 0;
    while gettr(&mut tr) != 0 {
        let sx = tr.sx as f64;
        let sy = tr.sy as f64;
        let gx = tr.gx as f64;
        let gy = tr.gy as f64;
        let unit = tr.counit;
        let mut factor = if tr.scalco == 0 { 1.0 } else { tr.scalco as f64 };
        if factor < 0.0 { factor = -1.0 / factor; }
        if unit == 2 && degree == 1 { factor /= 3600.0; }
        let mx = (0.5 * (sx + gx) * factor) as f32;
        let my = (0.5 * (sy + gy) * factor) as f32;
        efwrite(&[mx], &stdout());
        efwrite(&[my], &stdout());
        npairs = npairs.wrapping_add(1);
    }
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) != 0 {
        let fp = efopen(&outpar, "w");
        fprintf!(fp, "n=%d\n", npairs);
    }
    checkpars();
    exit(CWP_Exit());
}
