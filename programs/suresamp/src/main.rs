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
//! SURESAMP - Resample in time.
//! Port of su/main/stretching_moveout_resamp/suresamp.c (Revision 1.16).
//! Credits: CWP: Dave and Jack; CENPET: Werner M. Heigl;
//! RISSC: Nils Maercklin.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                   ",
    " SURESAMP - Resample in time                                       ",
    "                                                                   ",
    " suresamp <stdin >stdout  [optional parameters]                    ",
    "                                                                   ",
    " Required parameters:                                              ",
    "     none                                                          ",
    "                                                                   ",
    " Optional Parameters:                                              ",
    "    nt=tr.ns    number of time samples on output                   ",
    "    dt=         time sampling interval on output                   ",
    "                default is:                                        ",
    "                tr.dt/10^6     seismic data                        ",
    "                tr.d1          non-seismic data                    ",
    "    tmin=       time of first sample in output                     ",
    "                default is:                                        ",
    "                tr.delrt/10^3  seismic data                        ",
    "                tr.f1          non-seismic data                    ",
    "    rf=         resampling factor;                                 ",
    "                if defined, set nt=nt_in*rf and dt=dt_in/rf        ",
    "    verbose=0   =1 for advisory messages                           ",
    "                                                                   ",
    "                                                                   ",
    " Example 1: (assume original data had dt=.004 nt=256)              ",
    "    sufilter <data f=40,50 amps=1.,0. |                            ",
    "    suresamp nt=128 dt=.008 | ...                                  ",
    " Using the resampling factor rf, this example translates to:       ",
    "    sufilter <data f=40,50 amps=1.,0. | suresamp rf=0.5 | ...      ",
    "                                                                   ",
    " Note the typical anti-alias filtering before sub-sampling!        ",
    "                                                                   ",
    " Example 2: (assume original data had dt=.004 nt=256)              ",
    "    suresamp <data nt=512 dt=.002 | ...                            ",
    " or use:                                                           ",
    "    suresamp <data rf=2 | ...                                      ",
    "                                                                   ",
    " Example 3: (assume original data had d1=.1524 nt=8192)            ",
    "    sufilter <data f=0,1,3,3.28 amps=1,1,1,0 |                     ",
    "    suresamp <data nt=4096 dt=.3048 | ...                          ",
    "                                                                   ",
    " Example 4: (assume original data had d1=.5 nt=4096)               ",
    "    suresamp <data nt=8192 dt=.25 | ...                            ",
    "                                                                   ",
];

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    if gettr(&mut intrace) == 0 { err!("can't get first trace"); }
    let nt_in = intrace.ns as i32;
    if nt_in == 0 { err!("ns not set in header"); }
    let seismic = (0..=17).contains(&intrace.trid) || intrace.trid == 109;
    let (mut dt_in, tmin_in): (f32, f32);
    if seismic {
        if verbose != 0 { warn!("input is seismic data, trid=%d", intrace.trid); }
        dt_in = (intrace.dt as f64 / 1000000.0) as f32;
        tmin_in = (intrace.delrt as f64 / 1000.0) as f32;
    } else {
        if verbose != 0 { warn!("input is not seismic data, trid=%d", intrace.trid); }
        dt_in = intrace.d1;
        tmin_in = intrace.f1;
    }
    if dt_in == 0.0 { getparfloat("dt_in", &mut dt_in); }
    if dt_in == 0.0 { err!("dt or d1 not set in header or not getparred"); }
    if tmin_in == 0.0 && verbose != 0 { warn!("delrt or f1 not set in header"); }
    let mut rf = 0.0f32;
    getparfloat("rf", &mut rf);
    if rf < 0.0 { err!("factor rf=%g must be positive", rf); }
    let mut nt = 0;
    let mut dt = 0.0f32;
    if rf != 0.0 {
        if getparint("nt", &mut nt) == 0 { nt = nint((nt_in as f32 * rf) as f64); }
        if getparfloat("dt", &mut dt) == 0 { dt = dt_in / rf; }
    } else {
        if getparint("nt", &mut nt) == 0 { nt = nt_in; }
        if getparfloat("dt", &mut dt) == 0 { dt = dt_in; }
    }
    let mut tmin = 0.0f32;
    let tmin_is_set = getparfloat("tmin", &mut tmin) != 0;
    checkpars();
    CHECK_NT("nt", nt);
    let idt = nint(dt as f64 * 1000000.0);
    // A negative int converts to an enormous size_t in C's ealloc1float.
    if nt < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut t = vec![0.0f32; nt as usize];
    loop {
        if !tmin_is_set { tmin = tmin_in; }
        let mut tvalue = tmin;
        for time in &mut t { *time = tvalue; tvalue += dt; }
        outtrace.hdr = intrace.hdr;
        outtrace.ns = nt as u16;
        if seismic { outtrace.dt = idt as u16; } else { outtrace.d1 = dt; }
        if tmin_is_set {
            if seismic { outtrace.delrt = nint(tmin as f64 * 1000.0) as i16; }
            else { outtrace.f1 = tmin; }
        }
        ints8r(dt_in, tmin_in, &intrace.data[..nt_in as usize], 0.0, 0.0, &t, &mut outtrace.data[..nt as usize]);
        puttr(&mut outtrace);
        if gettr(&mut intrace) == 0 { break; }
    }
    exit(CWP_Exit());
}
