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

//! SUSHIFT - shifted/windowed traces in time.
//! Port of su/main/stretching_moveout_resamp/sushift.c (Revision 1.9).
//! Credits: Toralf Foerster, Institut fuer Ostseeforschung Warnemuende.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSHIFT - shifted/windowed traces in time\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sushift <stdin >stdout [tmin= ] [tmax= ]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmin=\t\t\tmin time to pass\t\t\t\t",
    " tmax=\t\t\tmax time to pass\t\t\t\t",
    " dt=                    sample rate in microseconds \t\t\t",
    " fill=0.0               value to place in padded samples \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " (defaults for tmin and tmax are calculated from the first trace.\t",
    " verbose=\t\t1 echos parameters to stdout\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Background :\t\t\t\t\t\t\t\t",
    " tmin and tmax must be given in seconds\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " In the high resolution single channel seismic profiling the sample \t",
    " interval is short, the shot rate and the number of samples are high.\t",
    " To reduce the file size the delrt time is changed during a profiling\t",
    " trip. To process and display a seismic section a constant delrt is\t",
    " needed. This program does this job.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The SEG-Y header variable delrt (delay in ms) is a short integer.\t",
    " That's why in the example shown below delrt is rounded to 123 !\t",
    "\t\t\t\t\t\t\t\t\t",
    "   ... | sushift tmin=0.1234 tmax=0.2234 | ...\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The dt= and fill= options are intended for manipulating velocity\t",
    " volumes in trace format.  In particular models which were hung\t",
    " from the water bottom when created & which then need to have the\t",
    " water layer added.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let mut dt = 0i32;
    if getparint("dt", &mut dt) == 0 { dt = tr.dt as i32; }
    let mut fill = 0.0f32;
    getparfloat("fill", &mut fill);
    let mut tmin = 0.0f32;
    let utmin = if getparfloat("tmin", &mut tmin) != 0 {
        nint(1e06 * tmin as f64)
    } else { 1000 * tr.delrt as i32 };
    let mut tmax = 0.0f32;
    let utmax = if getparfloat("tmax", &mut tmax) != 0 {
        nint(1e06 * tmax as f64)
    } else { utmin.wrapping_add((tr.ns as i32).wrapping_mul(dt)) };
    let mut verbose = 0i16;
    getparshort("verbose", &mut verbose);
    checkpars();
    if utmin >= utmax { err!("utmin (%i) >= utmax (%i) ?!", utmin, utmax); }
    if nint(1e-03 * utmin as f64) > 32767 || nint(1e-03 * utmin as f64) < -32768 {
        err!(" delrt exceeds short");
    }
    // C's integer division by zero yields zero on the target arm64 platform.
    let ns_out = utmax.wrapping_sub(utmin).checked_div(dt).unwrap_or(0) as u32;
    if verbose != 0 { warn!("IN %u samples  OUT %u samples", tr.ns as u32, ns_out); }
    if ns_out > SU_NFLTS as u32 { err!("ns_out=%u exceeds SU_NFLTS=%u", ns_out, SU_NFLTS as u32); }
    loop {
        let ut1_in = 1000 * tr.delrt as i32;
        let ut2_in = ut1_in.wrapping_add(dt.wrapping_mul(tr.ns as i32));
        if ut1_in != utmin || ut2_in != utmax {
            let (mut ut1, mut ut2, mut ut) = (0i32, 0i32, 0i32);
            if ut1_in < utmin && ut2_in <= utmax {
                ut1 = utmin.wrapping_sub(ut1_in);
                ut = ut2_in.wrapping_sub(utmin);
            } else if ut1_in < utmin && ut2_in > utmax {
                ut1 = utmin.wrapping_sub(ut1_in);
                ut = utmax.wrapping_sub(utmin);
            } else if ut1_in >= utmin && ut2_in <= utmax {
                ut2 = ut1_in.wrapping_sub(utmin);
                ut = ut2_in.wrapping_sub(ut1_in);
            } else if ut1_in >= utmin && ut2_in > utmax {
                ut2 = ut1_in.wrapping_sub(utmin);
                ut = utmax.wrapping_sub(ut1_in);
            } else {
                err!(" internal error : ut1=%i ut2=%i ut=%i ut1_in=%i ut2_in=%i", ut1, ut2, ut, ut1_in, ut2_in);
            }
            let it1 = ut1.checked_div(dt).unwrap_or(0);
            let it2 = ut2.checked_div(dt).unwrap_or(0);
            let it = ut.checked_div(dt).unwrap_or(0);
            let mut temp = vec![fill; ns_out as usize];
            for i in 0..it { temp[(it2 + i) as usize] = tr.data[(it1 + i) as usize]; }
            for i in 0..ns_out as usize { tr.data[i] = temp[i]; }
        }
        tr.delrt = (1e-03 * utmin as f64) as i16;
        tr.ns = ns_out as u16;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
