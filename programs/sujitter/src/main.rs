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

//! SUJITTER - Add random time shifts to seismic traces.
//! Port of su/main/noise/sujitter.c (Revision 1.3).
//! Credits (from the C source):
//! U of Houston: Chris Liner
//! UH: Chris added fldr, min, pon options 12/10/08
//
// Original source: Copyright (c) Colorado School of Mines, 2011.
// All rights reserved.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::sgn;
use oxwave_f_franuni::*;
use std::time::{SystemTime, UNIX_EPOCH};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUJITTER - Add random time shifts to seismic traces\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "     sujitter <stdin >stdout  [optional parameters]\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " \tseed=from_clock    \trandom number seed (integer)            ",
    "\tmin=1 \t\t\tminimum random time shift (samples)\t",
    "\tmax=1 \t\t\tmaximum random time shift (samples)\t",
    "\tpon=1 \t\t\tshift can be positive or negative\t",
    "\t\t\t\t=0 shift is positive only\t\t",
    "\tfldr=0 \t\t\teach trace has new shift\t\t",
    "\t\t\t\t=1 new shift when fldr header field changes",
    " Notes:\t\t\t\t\t\t\t\t",
    " Useful for simulating random statics. See also:  suaddstatics\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();

    // Hook up getpar.
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get information from first trace, before parsing parameters.
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let mut ishotold = tr.fldr;
    let _dt = (tr.dt as f64 / 1000000.0) as f32;

    // Get parameters.
    let mut min = 1i32;
    let mut max = 1i32;
    let mut pon = 1i32;
    let mut fldr = 0i32;
    let mut seed = 0u32;
    getparint("min", &mut min);
    getparint("max", &mut max);
    getparint("pon", &mut pon);
    getparint("fldr", &mut fldr);
    if min > max {
        err!("min>max... exit");
    }

    // Set seed; if not supplied, use clock.
    if getparuint("seed", &mut seed) == 0 {
        seed = match SystemTime::now().duration_since(UNIX_EPOCH) {
            Ok(elapsed) => elapsed.as_secs() as u32,
            Err(_) => err!("time() failed to set seed"),
        };
        if seed == u32::MAX {
            err!("time() failed to set seed");
        }
    }
    // The C prototype takes int, converting the unsigned parameter to signed.
    sranuni(seed as i32);

    // Loop on traces. The C program never calls checkpars().
    let mut itr = 1i32;
    let mut its = 0i32;
    loop {
        // Get fldr (shot) number.
        let ishot = tr.fldr;

        if itr == 1 {
            // Initial shift for shot 1 (used if fldr == 1).
            // C evaluates the shift in float, then truncates before signing.
            its = (min as f32 + max.wrapping_sub(min) as f32 * franuni()) as i32;
            if pon == 1 {
                // Include random sign to shift; 0.5 promotes to double in C.
                let sits = sgn(franuni() as f64 - 0.5) as i32;
                its = its.wrapping_mul(sits);
            }
        }

        if fldr == 0 {
            // Each trace gets random shift, including another draw on trace 1.
            its = (min as f32 + max.wrapping_sub(min) as f32 * franuni()) as i32;
            if pon == 1 {
                let sits = sgn(franuni() as f64 - 0.5) as i32;
                its = its.wrapping_mul(sits);
            }
        }

        if fldr == 1 && ishot != ishotold {
            // New shot needs new shift.
            its = (min as f32 + max.wrapping_sub(min) as f32 * franuni()) as i32;
            if pon == 1 {
                let sits = sgn(franuni() as f64 - 0.5) as i32;
                its = its.wrapping_mul(sits);
            }
            ishotold = ishot;
        }

        // Apply shift and output trace.
        if its <= 0 {
            let mut it = 0i32;
            // Preserve the C loop's access past ns, including stale samples
            // in the full SU_NFLTS allocation. The C does not zero the edges.
            // Stop before the C's out-of-allocation read. For the verified
            // nt=65534, its=-1 case, only the un-emitted data[65534] is skipped;
            // the minimal C bounds fix leaves stdout byte-identical.
            while it < nt.wrapping_sub(its) && it.wrapping_sub(its) < SU_NFLTS as i32 {
                tr.data[it as usize] = tr.data[it.wrapping_sub(its) as usize];
                it += 1;
            }
        } else {
            // The C starts at nt, even though puttr emits only indices < nt.
            // At nt=SU_NFLTS its first write is out of allocation; skipping
            // that write, as in the minimal C fix, does not change stdout.
            let mut it = nt.min(SU_NFLTS as i32 - 1);
            while it > its {
                tr.data[it as usize] = tr.data[it.wrapping_sub(its) as usize];
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
