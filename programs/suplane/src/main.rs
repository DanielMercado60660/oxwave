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

//! SUPLANE - create common offset data with up to three planes.
//! Port of su/main/synthetics_waveforms_testpatterns/suplane.c (Revision 1.20).
//! Credits: CWP: Chris Liner.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    "SUPLANE - create common offset data file with up to 3 planes\t",
    "\t\t\t\t\t\t\t\t",
    "suplane [optional parameters] >stdout\t \t\t\t",
    "\t\t\t\t\t\t\t\t",
    "Optional Parameters:\t\t\t\t\t\t",
    " npl=3\t\t\tnumber of planes\t\t\t",
    " nt=64 \t\tnumber of time samples\t\t\t",
    " ntr=32\t\tnumber of traces\t\t\t",
    " taper=0\t\tno end-of-plane taper\t\t\t",
    "\t\t\t= 1 taper planes to zero at the end\t",
    " offset=400 \t\toffset\t\t\t\t\t",
    " dt=0.004\t \ttime sample interval in seconds\t\t",
    "...plane 1 ...\t\t\t\t\t\t\t",
    "\tdip1=0\t\tdip of plane #1 (ms/trace)\t\t",
    " \tlen1= 3*ntr/4\tHORIZONTAL extent of plane (traces)\t",
    "\tct1= nt/2\ttime sample for center pivot\t \t",
    "\tcx1= ntr/2\ttrace for center pivot\t\t\t",
    "...plane 2 ...\t\t\t\t\t\t\t",
    "\tdip2=4\t\tdip of plane #2 (ms/trace)\t\t",
    "\tlen2= 3*ntr/4\tHORIZONTAL extent of plane (traces)\t",
    "\tct2= nt/2\ttime sample for center pivot \t\t",
    "\tcx2= ntr/2\ttrace for center pivot\t\t\t",
    "...plane 3 ...\t\t\t\t\t\t\t",
    "\tdip3=8\t\tdip of plane #3 (ms/trace)\t\t",
    "\tlen3= 3*ntr/4\tHORIZONTAL extent of plane (traces)\t",
    "\tct3= nt/2\ttime sample for center pivot\t\t",
    "\tcx3= ntr/2\ttrace for center pivot\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " liner=0\tuse parameters\t\t\t\t\t",
    "\t\t\t= 1 parameters set for 64x64 data set   ",
    "\t\t\twith separated dipping planes.\t\t",
];

fn plane(tr: &mut Trace, itr: i32, nt: i32, cx: i32, ct: i32, len: i32, dip: f32, msps: f32, taper: i32, tfe: &mut i32) {
    if itr >= cx - len / 2 && itr <= cx + len / 2 {
        *tfe += 1;
        let fit = ct as f32 - ((cx - itr) as f32 * dip) / msps;
        if fit >= 0.0 && fit <= nt as f32 {
            let itless = fit as usize;
            let eps = fit - itless as f32;
            // C's integer literal 1 converts to float: round the sum before
            // truncating, even when fit is immediately below an integer.
            let itmore = (fit + 1.0f32) as usize;
            // Guarding C's one-past-the-buffer += 0 leaves its output unchanged
            // (ub-2026-11 round 1). Skip that no-op when taper is off. Keep
            // writes beyond nt within the array, and reject other overruns.
            if itless >= tr.data.len()
                || (itmore >= tr.data.len() && (eps != 0.0 || taper == 1))
            {
                err!("plane interpolation exceeds trace data buffer");
            }
            tr.data[itless] = (tr.data[itless] as f64 + (1.0 - eps as f64)) as f32;
            if itmore < tr.data.len() {
                tr.data[itmore] += eps;
            }
            if taper == 1 {
                if *tfe == 1 || *tfe == len + 1 {
                    tr.data[itless] = (tr.data[itless] as f64 / 6.0) as f32;
                    tr.data[itmore] = (tr.data[itmore] as f64 / 6.0) as f32;
                }
                if *tfe == 2 || *tfe == len {
                    tr.data[itless] = (tr.data[itless] as f64 / 3.0) as f32;
                    tr.data[itmore] = (tr.data[itmore] as f64 / 3.0) as f32;
                }
            }
        }
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut tr = Trace::new();
    let mut nt = 64; getparint("nt", &mut nt);
    CHECK_NT("nt", nt); tr.ns = nt as u16;
    let mut ntr = 32; getparint("ntr", &mut ntr);
    let mut dt = 0.004f32; getparfloat("dt", &mut dt); tr.dt = nint(dt as f64 * 1000000.0) as u16;
    let mut offset = 400.0f32; getparfloat("offset", &mut offset); tr.offset = offset as i32;
    let mut npl = 3; getparint("npl", &mut npl);
    let mut taper = 0; getparint("taper", &mut taper);
    let mut dip1 = 0.0f32; getparfloat("dip1", &mut dip1);
    let mut len1 = 3 * ntr / 4; getparint("len1", &mut len1);
    let mut ct1 = nt / 2; getparint("ct1", &mut ct1); ct1 -= 1;
    let mut cx1 = ntr / 2; getparint("cx1", &mut cx1); cx1 -= 1;
    let mut dip2 = 4.0f32; getparfloat("dip2", &mut dip2);
    let mut len2 = 3 * ntr / 4; getparint("len2", &mut len2);
    let mut ct2 = nt / 2; getparint("ct2", &mut ct2); ct2 -= 1;
    let mut cx2 = ntr / 2; getparint("cx2", &mut cx2); cx2 -= 1;
    let mut dip3 = 8.0f32; getparfloat("dip3", &mut dip3);
    let mut len3 = 3 * ntr / 4; getparint("len3", &mut len3);
    let mut ct3 = nt / 2; getparint("ct3", &mut ct3); ct3 -= 1;
    let mut cx3 = ntr / 2; getparint("cx3", &mut cx3); cx3 -= 1;
    let mut liner = 0; getparint("liner", &mut liner);
    if liner == 1 {
        nt = 64; tr.ns = nt as u16; ntr = 64; npl = 3;
        dip1 = 0.0; len1 = ntr / 4; ct1 = nt / 2 - 1; cx1 = 3 * ntr / 4 - 1;
        dip2 = 4.0; len2 = ntr / 4; ct2 = nt / 2 - 1; cx2 = ntr / 2 - 1;
        dip3 = 8.0; len3 = ntr / 4; ct3 = nt / 2 - 1; cx3 = ntr / 4 - 1;
    }
    let msps = (dt as f64 * 1000.0) as f32;
    let (mut tfe1, mut tfe2, mut tfe3) = (0, 0, 0);
    for itr in 0..ntr {
        // A negative nt makes C's memset overrun the buffer. Check here so
        // liner=1 and empty gathers retain their original behavior.
        if nt < 0 {
            err!("nt=%d must not be negative", nt);
        }
        tr.data[..nt as usize].fill(0.0);
        plane(&mut tr, itr, nt, cx1, ct1, len1, dip1, msps, taper, &mut tfe1);
        if npl > 1 { plane(&mut tr, itr, nt, cx2, ct2, len2, dip2, msps, taper, &mut tfe2); }
        if npl > 2 { plane(&mut tr, itr, nt, cx3, ct3, len3, dip3, msps, taper, &mut tfe3); }
        tr.tracl = itr + 1; tr.tracr = itr + 1;
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
