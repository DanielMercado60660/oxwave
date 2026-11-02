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
//! SUSORTY - generate a shot-ordered geometry visualization data set.
//! Port of su/main/windowing_sorting_muting/susorty.c (Revision 1.5).
//! Credits (from the C source): CWP: Chris Liner 10.09.01.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUSORTY - make a small 2-D common shot off-end  \t\t",
    "\t    data set in which the data show geometry \t\t",
    "\t    values to help visualize data sorting.\t\t",
    "\t\t\t\t\t\t\t\t",
    "  susorty [optional parameters] > out_data_file  \t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "\tnt=100 \t\tnumber of time samples\t\t\t",
    "\tnshot=10 \tnumber of shots\t\t\t\t",
    "\tdshot=10 \tshot interval (m)\t\t\t",
    "\tnoff=20 \tnumber of offsets\t\t\t",
    "\tdoff=20 \toffset increment (m)\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " Creates a common shot su data file for sort visualization\t",
    "\t       time samples           quantity\t\t\t",
    "\t       ----------------      ----------\t\t\t",
    "\t       first   25%           shot coord\t\t\t",
    "\t       second  25%           rec coord\t\t\t",
    "\t       third   25%           offset\t\t\t",
    "\t       fourth  25%           cmp coord\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " 1. default is shot ordered (hsv2 cmap looks best to me)\t",
    " susorty | suximage legend=1 units=meters cmap=hsv2\t\t",
    "\t\t\t\t\t\t\t\t",
    " 2. sort on cmp (note random order within a cmp)\t\t",
    " susorty | susort cdp > junk.su \t\t\t\t",
    " suximage < junk.su legend=1 units=meters cmap=hsv2\t\t",
    "\t\t\t\t\t\t\t\t",
    " 3. sort to cmp and subsort on offset \t \t\t\t",
    " susorty | susort cdp offset > junk.su \t\t\t",
    " suximage < junk.su legend=1 units=meters cmap=hsv2\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut nt = 100;
    getparint("nt", &mut nt);
    CHECK_NT("nt", nt);
    if nt < 0 { err!("nt must be nonnegative"); }
    tr.ns = nt as u16;
    let mut nshot = 10;
    getparint("nshot", &mut nshot);
    let mut noff = 24;
    getparint("noff", &mut noff);
    let mut dt: f32 = 0.004;
    getparfloat("dt", &mut dt);
    tr.dt = (dt * 1000000.0) as i32 as u16;
    let mut dshot: f32 = 10.0;
    getparfloat("dshot", &mut dshot);
    let mut doff: f32 = 20.0;
    getparfloat("doff", &mut doff);
    checkpars();

    for ishot in 0..nshot {
        let sx = ishot as f32 * dshot;
        for ioff in 0..noff {
            let offset = (ioff + 1) as f32 * doff;
            let gx = sx + offset;
            let cmp = ((sx + gx) as f64 / 2.0) as f32;
            tr.data[..nt as usize].fill(0.0);
            for it in 0..nt / 4 { tr.data[it as usize] = sx; }
            for it in nt / 4..nt / 2 { tr.data[it as usize] = gx; }
            for it in nt / 2..3 * nt / 4 { tr.data[it as usize] = offset; }
            for it in 3 * nt / 4..nt { tr.data[it as usize] = cmp; }
            tr.sx = sx as i32;
            tr.gx = gx as i32;
            tr.offset = offset as i32;
            tr.cdp = cmp as i32;
            tr.tracl = ishot.wrapping_mul(nshot).wrapping_add(ioff).wrapping_add(1);
            puttr(&mut tr);
        }
    }
    exit(CWP_Exit());
}
