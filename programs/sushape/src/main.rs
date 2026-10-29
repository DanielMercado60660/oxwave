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
//! SUSHAPE - Wiener shaping filter. Port of su/main/decon_shaping/sushape.c (Revision 1.16).
//! Credits: CWP: Jack Cohen; John Stockwell added wfile and dfile options.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::closeto;
use oxwave_correlation::xcor;
use oxwave_toeplitz::stoepf;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUSHAPE - Wiener shaping filter\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  sushape <stdin >stdout  [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " w=\t\tvector of input wavelet to be shaped or ...\t\t",
    " ...or ... \t\t\t\t\t\t\t\t",
    " wfile=        ... file containing input wavelet in SU (SEGY trace) format",
    " d=\t\tvector of desired output wavelet or ...\t\t\t",
    " ...or ... \t\t\t\t\t\t\t\t",
    " dfile=        ... file containing desired output wavelet in SU format\t",
    " dt=tr.dt\t\tif tr.dt is not set in header, then dt is mandatory",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " nshape=trace\t\tlength of shaping filter\t\t\t",
    " pnoise=0.001\t\trelative additive noise level\t\t\t",
    " showshaper=0\t\t=1 to show shaping filter \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " verbose=0\t\tsilent; =1 chatty\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "Notes:\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Example of commandline input wavelets: \t\t\t\t",
    "sushape < indata  w=0,-.1,.1,... d=0,-.1,1,.1,... > shaped_data\t",
    " \t\t\t\t\t\t\t\t\t",
    "sushape < indata  wfile=inputwavelet.su dfile=desire.su > shaped_data\t",
    " \t\t\t\t\t\t\t\t\t",
    " To get the shaping filters into an ascii file:\t\t\t",
    " ... | sushape ... showwshaper=1 2>file | ...   (sh or ksh)\t\t",
    " (... | sushape ... showshaper=1 | ...) >&file  (csh)\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    if gettr(&mut intrace) == 0 { err!("can't get first trace"); }
    let nt = intrace.ns as usize;
    let mut dt = (intrace.dt as f64 / 1000000.0) as f32;
    if dt == 0.0 { mustgetparfloat("dt", &mut dt); }
    let mut showshaper = 0;
    let mut nshape = nt as i32;
    let mut pnoise = 0.001f32;
    let mut verbose = 0;
    getparint("showshaper", &mut showshaper);
    getparint("nshape", &mut nshape);
    getparfloat("pnoise", &mut pnoise);
    getparint("verbose", &mut verbose);
    let mut dfile = String::new();
    let mut wfile = String::new();
    getparstring("dfile", &mut dfile);
    getparstring("wfile", &mut wfile);
    let d = if dfile.is_empty() {
        let nd = countparval("d");
        if nd == 0 { err!("must specify d= desired wavelet"); }
        let mut d = vec![0.0f32; nd as usize];
        getparfloat("d", &mut d);
        d
    } else {
        let dfp = fopen(&dfile, "r").unwrap_or_else(|| err!("cannot open dfile=%s\n", dfile));
        let mut dtr = Trace::new();
        if fgettr(&dfp, &mut dtr) == 0 { err!("can't get input wavelet"); }
        dtr.data[..dtr.ns as usize].to_vec()
    };
    let w = if wfile.is_empty() {
        let nw = countparval("w");
        if nw == 0 { err!("must specify w= desired wavelet"); }
        let mut w = vec![0.0f32; nw as usize];
        getparfloat("w", &mut w);
        w
    } else {
        let wfp = fopen(&wfile, "r").unwrap_or_else(|| err!("cannot open wfile=%s\n", wfile));
        let mut wtr = Trace::new();
        if fgettr(&wfp, &mut wtr) == 0 { err!("can't get desired output wavelet"); }
        wtr.data[..wtr.ns as usize].to_vec()
    };
    checkpars();
    if nshape < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    if nshape == 0 {
        err!("nshape must be positive");
    }
    let mut shaper = vec![0.0f32; nshape as usize];
    let mut spiker = vec![0.0f32; nshape as usize];
    let mut crosscorr = vec![0.0f32; nshape as usize];
    let mut autocorr = vec![0.0f32; nshape as usize];
    xcor(0, &w, 0, &w, 0, &mut autocorr);
    xcor(0, &w, 0, &d, 0, &mut crosscorr);
    if closeto(autocorr[0], 0.0) { err!("can't shape with zero wavelet"); }
    autocorr[0] = (autocorr[0] as f64 * (1.0 + pnoise as f64)) as f32;
    stoepf(&autocorr, &crosscorr, &mut shaper, &mut spiker);
    if showshaper != 0 {
        if verbose != 0 { warn!("Shaping filter:"); }
        for i in 0..nshape as usize {
            fprintf!(stderr(), "%10g%c", shaper[i], if i % 6 == 5 || i == nshape as usize - 1 { '\n' as i32 } else { ' ' as i32 });
        }
    }
    // xcor reverses its first input internally; reversing it here gives convolution.
    let reversed: Vec<f32> = shaper.iter().rev().copied().collect();
    let origin = (w.len() as i32 - d.len() as i32) / 2;
    loop {
        xcor(1 - nshape - origin, &reversed, 0, &intrace.data[..nt], 0, &mut outtrace.data[..nt]);
        outtrace.hdr = intrace.hdr.clone();
        puttr(&mut outtrace);
        if gettr(&mut intrace) == 0 { break; }
    }
    exit(CWP_Exit());
}
