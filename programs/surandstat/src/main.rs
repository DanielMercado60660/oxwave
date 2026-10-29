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

//! SURANDSTAT - add random static shifts to seismic traces.
//! Port of su/main/statics/surandstat.c (Revision 1.5).
//! Credits (from the C source): U Houston: Chris Liner c. 2009.
//! Original source: Copyright (c) Colorado School of Mines, 2011.
//! All rights reserved.

use oxwave_core::prelude::*;
use oxwave_f_franuni::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SURANDSTAT - Add RANDom time shifts STATIC errors to seismic traces\t",
    "\t\t\t\t\t\t\t\t\t",
    "     surandstat <stdin >stdout  [optional parameters]\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " \tseed=from_clock    \trandom number seed (integer)            ",
    "\tmax=tr.dt \t\tmaximum random time shift (ms)\t\t",
    "\tscale=1.0\t\tscale factor for shifts\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    // The positive-shift loop writes data[nt], including when nt == SU_NFLTS.
    // Keep a scratch cell for that otherwise out-of-bounds store; puttr never
    // emits it. The original binary accepts a full-length trace this way.
    tr.data.push(0.0);
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let mut max = 0i32;
    if getparint("max", &mut max) == 0 {
        max = ((dt as i32) as f64 * 1000.0) as i32;
    }
    let mut scale = 0.0f32;
    if getparfloat("scale", &mut scale) == 0 {
        scale = 1.0;
    }
    let mut seed = 0u32;
    if getparuint("seed", &mut seed) == 0 {
        seed = match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
            Ok(t) => t.as_secs() as u32,
            Err(_) => err!("time() failed to set seed"),
        };
        if seed == u32::MAX {
            err!("time() failed to set seed");
        }
    }
    checkpars();
    let mut itr = 0i32;
    loop {
        // C reseeds with sranuni(seed) on EVERY trace, then calls franuni().
        sranuni(seed as i32);
        // scale*2*max is evaluated in float, then 0.5 promotes to double.
        let its = ((scale * 2.0f32 * max as f32) as f64 * (0.5 - franuni() as f64) + 1.0) as i32;

        if its <= 0 {
            let mut it = 0i32;
            while it < nt.wrapping_sub(its) {
                // Like C, retain samples past ns in the full trace buffer.
                tr.data[it as usize] = tr.data[it.wrapping_sub(its) as usize];
                it = it.wrapping_add(1);
            }
        } else {
            let mut it = nt;
            while it > its {
                // C starts at data[nt] and leaves vacated samples unchanged.
                tr.data[it as usize] = tr.data[(it - its) as usize];
                it -= 1;
            }
        }
        puttr(&mut tr);
        itr = itr.wrapping_add(1);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
