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

//! SUMIX - compute weighted moving average on a panel of seismic data.
//! Port of su/main/operations/sumix.c (Revision 1.11).
//! Credits: CWP: John Stockwell; UArk: Chris Liner.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIX - compute weighted moving average (trace MIX) on a panel\t",
    "\t  of seismic data\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sumix <stdin >sdout \t\t\t\t\t\t\t",
    " mix=.6,1,1,1,.6\tarray of weights for weighted average\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: \t\t\t\t\t\t\t\t",
    " The number of values defined by mix=val1,val2,... determines the number",
    " of traces to be averaged, the values determine the weights.\t\t",
    " Mix applied left-to-right (in 2D data) and only begins after the\t",
    " nmix traces, where nmix is the number of mix weights to avoid \t",
    " left end edge effects. \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Examples: \t\t\t\t\t\t\t\t",
    " sumix <stdin mix=.6,1,1,1,.6 >sdout \t(default) mix over 5 traces weights",
    " sumix <stdin mix=1,1,1 >sdout \tsimple 3 trace moving average\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let mut nmix = countparval("mix") as usize;
    let mut mix;
    if nmix != 0 {
        mix = vec![0.0f32; nmix];
        getparfloat("mix", &mut mix);
    } else {
        nmix = 5;
        mix = vec![0.6, 1.0, 1.0, 1.0, 0.6];
    }
    checkpars();
    for imix in 0..nmix {
        mix[imix] /= nmix as f32;
    }
    let mut data = vec![vec![0.0f32; nt]; nmix];
    let mut temp = vec![0.0f32; nt];
    let mut itr = 0i32;
    loop {
        itr = itr.wrapping_add(1);
        temp.fill(0.0);
        data[0].copy_from_slice(&tr.data[..nt]);
        if itr >= nmix as i32 {
            for it in 0..nt {
                for imix in 0..nmix {
                    temp[it] += data[imix][it] * mix[imix];
                }
                tr.data[it] = temp[it];
            }
        }
        for imix in (1..nmix).rev() {
            for it in 0..nt {
                data[imix][it] = data[imix - 1][it];
            }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
