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
//! SURECMO - compensate for a continuously moving marine streamer.
//! Port of su/main/operations/surecmo.c (Revision 1.1).
//! Credits (from the C source): CWP: Taylor Goss, May 2017.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SURECMO - compensate for the continuously moving streamer in marine   ",
    "           seismic acquisition (assume far offset is first channel)    ",
    "           right now valid for data where tmax < 2*dx/vb               ",
    "\t\t\t\t\t\t\t\t\t",
    " surecmo <stdin >sdout [vb= ] [dx= ] [dt= ] [fill= ] [tmax= ]          ",
    "                                                                       ",
    " vb=         boat speed in m/s                                         ",
    " dx=         channel spacing in m                                      ",
    " dt=         sample rate in microseconds                               ",
    " fill=0.0    value to pad                                              ",
    "\t\t\t\t\t\t\t\t\t",
    " Examples: \t\t\t\t\t\t\t\t",
    " surecmo <stdin >sdout vb=2.5 dt=25 dt=4000                            ",
    "                 (valid for tmax < 20s) \t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut dt: i32 = 0;
    if getparint("dt", &mut dt) == 0 { dt = tr.dt as i32; }
    let mut vb: f32 = 0.0;
    let mut dx: f32 = 0.0;
    getparfloat("vb", &mut vb);
    getparfloat("dx", &mut dx);
    let nmix = 3;
    let mut mix = vec![0.0f32; nmix];
    checkpars();
    let mut data = vec![vec![0.0f32; nt]; nmix];
    let mut temp = vec![0.0f32; nt];
    loop {
        temp.fill(0.0);
        data[0].copy_from_slice(&tr.data[..nt]);
        for it in 0..nt {
            mix[1] = (((vb * it as f32 * dt as f32) as f64) / (dx as f64 * 1e6)) as f32;
            mix[0] = (1.0f64 - mix[1] as f64) as f32;
            mix[2] = 0.0;
            if mix[1] > 1.0 {
                mix[2] = mix[1] - 1.0;
                mix[1] = (1.0f64 - mix[2] as f64) as f32;
                mix[0] = 0.0;
            }
            for imix in 0..nmix { temp[it] += data[imix][it] * mix[imix]; }
            tr.data[it] = temp[it];
        }
        for imix in (1..nmix).rev() {
            for it in 0..nt { data[imix][it] = data[imix - 1][it]; }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
