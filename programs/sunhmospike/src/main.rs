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
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! SUNHMOSPIKE - non-hyperbolic moveout spike test data.
//! Port of su/main/synthetics_waveforms_testpatterns/sunhmospike.c (Revision 1.4).
//! Credits: CWP: Shuki Ronen, Chris Liner; modified by John Anderson, April 1994.
//! gofx author: John Anderson (visitor to CSM from Mobil), Spring 1993.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUNHMOSPIKE - generates SPIKE test data set with a choice of several   ",
    "   Non-Hyperbolic MOveouts\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   sunhmospike [optional parameters] > out_data_file  \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tnt=300\t  number of time samples\t\t\t\t",
    "\tntr=20\t  number of traces\t\t\t\t\t",
    "\tdt=0.001\ttime sampling interval in seconds\t\t",
    "\toffref=2000\treference offset\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tgopt=\t\t1 = parabolic transform model\t\t\t",
    "\t\t\t2 = Foster/Mosher pseudo hyperbolic option model",
    "\t\t\t3 = linear tau-p model\t\t\t\t",
    "\tdepthref=400\treference depth used when gopt=2\t\t",
    "\toffinc=100\toffset increment\t\t\t\t",
    "\tnspk=4\t  number of events\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tp1 = 0\t\tevent moveout for event #1 in ms on reference offset",
    "\tt1 = 100\tintercept time ms event #1\t\t\t",
    "\ta1 = 1.0\tamplitude for event #1\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tp2 = 200\tevent moveout for event #2 in ms on reference offset",
    "\tt2 = 100\tintercept time ms for spike #2\t\t\t",
    "\ta2 = 1.0\tamplitude for event #2\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tp3 = 0;\t\tevent moveout for event #3 in ms on reference offset",
    "\tt3 = 200\tintercept time for spike #3\t\t\t",
    "\ta3 = 1.0\tamplitude for event #3\t\t\t\t",
    "\tp4 = 120\t event moveout for event #4 in ms on reference offset",
    "\t\t\t\t\t\t\t\t\t",
    "\tt4 = 200\tintercept time s for spike #4\t\t\t",
    "\ta4 = 1.0\tamplitude for event #4\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tcdp = 1\t output cdp number\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Creates a common cdp su data file with up to four spike events\t",
    " for impulse response studies for suradon, and sutifowler\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn gofx(igopt: i32, offset: f32, intercept_off: f32, refdepth: f32) -> f32 {
    let offset = offset - intercept_off;
    if igopt == 1 { return offset * offset; }
    if igopt == 2 { return ((refdepth * refdepth + offset * offset) as f64).sqrt() as f32; }
    if igopt == 3 { return offset; }
    if igopt == 4 { return (offset as f64).abs() as f32; }
    offset
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut igopt, mut ntr, mut nt, mut cdp) = (1, 20, 300, 1);
    getparint("gopt", &mut igopt); getparint("ntr", &mut ntr);
    getparint("nt", &mut nt); getparint("cdp", &mut cdp);
    let (mut dt, mut offref, mut offinc, mut refdepth) = (0.001f32, 2000.0f32, 100.0f32, 400.0f32);
    let mut nspk = 4;
    getparfloat("dt", &mut dt); getparfloat("offref", &mut offref);
    getparfloat("offinc", &mut offinc); getparfloat("refdepth", &mut refdepth);
    getparint("nspk", &mut nspk);
    let (mut p1, mut t1, mut a1) = (0.0f32, 100.0f32, 1.0f32);
    let (mut p2, mut t2, mut a2) = (200.0f32, 100.0f32, 1.0f32);
    let (mut p3, mut t3, mut a3) = (0.0f32, 200.0f32, 1.0f32);
    let (mut p4, mut t4, mut a4) = (100.0f32, 200.0f32, 1.0f32);
    getparfloat("p1", &mut p1); getparfloat("t1", &mut t1); getparfloat("a1", &mut a1);
    getparfloat("p2", &mut p2); getparfloat("t2", &mut t2); getparfloat("a2", &mut a2);
    getparfloat("p3", &mut p3); getparfloat("t3", &mut t3); getparfloat("a3", &mut a3);
    getparfloat("p4", &mut p4); getparfloat("t4", &mut t4); getparfloat("a4", &mut a4);
    t1 = (t1 as f64 / (1000.0 * dt as f64)) as f32;
    t2 = (t2 as f64 / (1000.0 * dt as f64)) as f32;
    t3 = (t3 as f64 / (1000.0 * dt as f64)) as f32;
    t4 = (t4 as f64 / (1000.0 * dt as f64)) as f32;
    p1 /= gofx(igopt, offref, 0.0, refdepth) * 1000.0f32 * dt;
    p2 /= gofx(igopt, offref, 0.0, refdepth) * 1000.0f32 * dt;
    p3 /= gofx(igopt, offref, 0.0, refdepth) * 1000.0f32 * dt;
    p4 /= gofx(igopt, offref, 0.0, refdepth) * 1000.0f32 * dt;
    for itr in 0..ntr {
        for j in 0..nt { tr.data[j as usize] = 0.0; }
        let off = (itr + 1) as f32 * offinc;
        let it1 = (t1 + p1 * gofx(igopt, off, 0.0, refdepth)) as i32;
        let it2 = (t2 + p2 * gofx(igopt, off, 0.0, refdepth)) as i32;
        let it3 = (t3 + p3 * gofx(igopt, off, 0.0, refdepth)) as i32;
        let it4 = (t4 + p4 * gofx(igopt, off, 0.0, refdepth)) as i32;
        if nspk > 0 && it1 > 0 && it1 < nt { tr.data[it1 as usize] += a1; }
        if nspk > 1 && it2 > 0 && it2 < nt { tr.data[it2 as usize] += a2; }
        if nspk > 2 && it3 > 0 && it3 < nt { tr.data[it3 as usize] += a3; }
        if nspk > 3 && it4 > 0 && it4 < nt { tr.data[it4 as usize] += a4; }
        tr.tracl = itr + 1; tr.cdp = cdp; tr.offset = off as i32;
        tr.ns = nt as u16; tr.dt = (dt * 1000000.0f32) as u32 as u16;
        tr.f1 = 0.0; tr.d1 = dt; tr.f2 = offinc; tr.d2 = offinc;
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
