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

//! SURANDSPIKE - make a small data set of random spikes.
//! Port of su/main/synthetics_waveforms_testpatterns/surandspike.c (Revision 1.4).
//! Credits: ARAMCO: Chris Liner.

use oxwave_core::prelude::*;
use oxwave_f_franuni::*;
use std::time::{SystemTime, UNIX_EPOCH};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SURANDSPIKE - make a small data set of RANDom SPIKEs \t\t",
    "\t\t\t\t\t\t\t\t",
    "   surandspike [optional parameters] > out_data_file  \t\t",
    "\t\t\t\t\t\t\t\t",
    " Creates a common offset su data file with random spikes\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "\tn1=500 \t\t\tnumber of time samples\t\t",
    "\tn2=200\t\t\tnumber of traces\t\t",
    " \tdt=0.002 \t\ttime sampling interval in seconds\t",
    "\tnspk=20\t\t\tnumber of spikes per trace\t",
    "\tamax=0.2\t\tabs(max) spike value\t\t",
    "\tmode=1\t\t\tdifferent spikes on each trace\t",
    "\t\t\t\t=2 same spikes on each trace\t",
    " \tseed=from_clock    \trandom number seed (integer)    ",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();

    // Initialize; stdin is not used.
    initargs(&argv);
    requestdoc(0, SDOC);

    // Get mode from user.
    let mut mode = 1;
    getparint("mode", &mut mode);
    if mode < 1 || mode > 2 {
        err!("mode must be 1 or 2");
    }
    // Set seed (use the clock if not supplied).
    let mut seed = 0u32;
    if getparuint("seed", &mut seed) == 0 {
        seed = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs() as u32;
        if seed == u32::MAX {
            err!("time() failed to set seed");
        }
    }
    sranuni(seed as i32);

    // Other getpars.
    let mut n1 = 500;
    let mut dt = 0.002f32;
    let mut n2 = 100;
    let mut nspk = 20;
    let mut amax = 0.2f32;
    getparint("n1", &mut n1);
    getparfloat("dt", &mut dt);
    getparint("n2", &mut n2);
    getparint("nspk", &mut nspk);
    getparfloat("amax", &mut amax);
    checkpars();
    if amax < 0.0 {
        amax = -amax;
    }
    let tmax = (n1.wrapping_sub(1) as f32) * dt;

    // Set output header words.
    // The C build converts to an unsigned register before storing its low
    // 16 bits. Casting directly to u16 would saturate intervals above 65535 us.
    tr.dt = ((dt as f64) * 1000000.0) as u32 as u16;
    tr.ns = n1 as u16;

    // Loop over traces.
    for itr in 0..n2 {
        if mode == 1 || (mode == 2 && itr == 0) {
            // Zero trace.
            tr.data[..n1 as usize].fill(0.0);

            // Loop over spikes.
            for _ispk in 0..nspk {
                // Get spike time sample; both operations stay in float.
                let t = franuni() * tmax;
                let it = (t / dt + 1.0f32) as i32;

                // Get spike amplitude; the C literal 0.5 promotes to double.
                let a = (2.0f64 * (0.5f64 - franuni() as f64) * amax as f64) as f32;

                // Put spike in trace (ignore collisions). The full SU sample
                // buffer also retains placements at n1, outside the output.
                tr.data[it as usize] = a;
            }
        }

        // Bump trace counter and write trace.
        tr.tracl = itr + 1;
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
