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
//! SUCENTSAMP - centroid sample seismic traces.
//! Port of su/main/amplitudes/sucentsamp.c.
//! Credits: Providence Technologies: Tom Morgan.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUCENTSAMP - CENTRoid SAMPle seismic traces\t\t\t",
    " \t\t\t\t\t\t\t\t",
    "  sucentsamp <stdin [optional parameters] >sdout  \t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tdt=from header\t\tsampling interval\t\t",
    "\tverbose=1\t\t=0 to stop advisory messages\t",
    " \t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t",
    " This program takes seismic traces as input, and returns traces",
    " consisting of spikes of height equal to the area of each lobe ",
    " of each oscillation, located at the centroid of the lobe in\t",
    " question. The height of each spike equal to the area of the   ",
    " corresponding lobe.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Caveat: No check is made that the data are real time traces!\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 1;
    getparint("verbose", &mut verbose);
    let mut nvals_min = 1;
    getparint("nvals_min", &mut nvals_min);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    checkpars();
    if dt == 0.0 {
        dt = 0.004;
        if verbose != 0 { warn!("dt not set, assumed to be .004"); }
    }
    let invdt = (1.0f64 / dt as f64) as f32;
    let hdt = (0.5f64 * dt as f64) as f32;
    // The original can write one sample beyond ct when the centroid rounds to nt.
    let mut rt = vec![0.0f32; nt.max(1)];
    let mut ct = vec![0.0f32; nt + 1];
    let mut mt = vec![0.0f32; nt.max(1)];
    let mut time = vec![0.0f32; nt.max(1)];
    for i in 0..nt { time[i] = (i + 1) as f32 * dt; }
    let mut inflect = false;
    let mut zero_cross = false;
    let mut max_passed = false;
    loop {
        rt[..nt].copy_from_slice(&tr.data[..nt]);
        ct[..nt].fill(0.0);
        let mut first = 0usize;
        mt[0] = rt[0].abs();
        let mut sum_amp = rt[0];
        let mut t_mom = time[0] * mt[0];
        for i in 1..nt {
            mt[i] = rt[i].abs();
            if rt[i] * rt[i-1] > 0.0 {
                if mt[i] > mt[i-1] {
                    if max_passed { inflect = true; }
                } else { max_passed = true; }
            } else { zero_cross = true; }
            if !zero_cross && !inflect {
                sum_amp = sum_amp + rt[i];
                t_mom = t_mom + time[i] * mt[i];
            } else {
                let mut last = i - 1;
                if inflect {
                    last = i;
                    rt[last] = (rt[last] as f64 * 0.5) as f32;
                    mt[last] = rt[last].abs();
                    sum_amp = sum_amp + rt[last];
                    t_mom = t_mom + time[last] * mt[last];
                }
                let nvals = last - first + 1;
                if nvals == 1 {
                    if nvals as i32 >= nvals_min { ct[i] = (rt[i] as f64 * 0.5) as f32; }
                    first = i;
                    inflect = false;
                    max_passed = false;
                    zero_cross = false;
                    sum_amp = rt[i];
                    t_mom = time[i] * mt[i];
                } else {
                    let mut a_height = mt[first];
                    if mt[first] > mt[last] { a_height = mt[last]; }
                    let mut t_width = time[last] - time[first] + dt;
                    let mut a_mom = (a_height as f64 * 0.5) as f32;
                    let mut a_cen = t_width * a_height * a_mom;
                    let mut small = first;
                    if mt[first] > mt[last] { small = last; }
                    for _k in 1..nvals {
                        let prev = small;
                        if prev == first { first += 1; } else { last -= 1; }
                        small = first;
                        if mt[first] > mt[last] { small = last; }
                        a_height = rt[small] - rt[prev];
                        a_mom = (rt[prev] as f64 + a_height as f64 * 0.5) as f32;
                        t_width = t_width - dt;
                        a_cen = a_cen + a_mom * t_width * a_height;
                    }
                    if sum_amp == 0.0 { warn!("i = %d, sum_amp = %f,divide by zero !", i as i32, sum_amp); }
                    a_cen = a_cen / (sum_amp * dt);
                    let t_cen = (t_mom as f64 / (sum_amp as f64).abs()) as f32;
                    sum_amp = rt[i];
                    t_mom = time[i] * mt[i];
                    let isamp = (t_cen * invdt + hdt) as i32;
                    if nvals as i32 >= nvals_min {
                        // Out-of-range C writes have undefined behavior.
                        if let Some(sample) = ct.get_mut(isamp as usize) { *sample = a_cen; }
                    }
                    first = i;
                    inflect = false;
                    max_passed = false;
                    zero_cross = false;
                }
            }
        }
        tr.data[..nt].copy_from_slice(&ct[..nt]);
        tr.ns = nt as u16;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
