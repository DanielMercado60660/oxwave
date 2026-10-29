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

//! SUACOR - auto-correlation. Port of su/main/convolution_correlation/suacor.c (Revision 1.16).
//! Credits (from the C source): CWP: Dave Hale.

use oxwave_core::prelude::*;
use oxwave_correlation::xcor;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUACOR - auto-correlation\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suacor <stdin >stdout [optional parms]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " ntout=101\todd number of time samples output\t\t\t",
    " norm=1\tif non-zero, normalize maximum absolute output to 1\t",
    " sym=1\t\tif non-zero, produce a symmetric output from\t\t",
    "\t\t\tlag -(ntout-1)/2 to lag +(ntout-1)/2\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let mut ntout = 101i32;
    let mut norm = 1i32;
    let mut sym = 1i32;
    getparint("ntout", &mut ntout);
    getparint("norm", &mut norm);
    getparint("sym", &mut sym);
    checkpars();
    if ntout < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    // C's memcpy writes past tr.data for larger outputs (undefined behavior).
    // Reject them before allocation, indexing, or narrowing ntout into tr.ns.
    if ntout as usize > tr.data.len() {
        err!("ntout=%d exceeds trace capacity (%d)", ntout, tr.data.len() as i32);
    }
    // C's zero-byte allocation is still accessed at zero lag; no samples are output.
    let mut temp = vec![0.0f32; (ntout as usize).max(1)];
    let istart = if sym == 0 { 0 } else { -(ntout - 1) / 2 };
    let izero = -istart;
    loop {
        xcor(0, &tr.data[..nt], 0, &tr.data[..nt], istart, &mut temp[..ntout as usize]);
        if norm != 0 {
            let scale = (1.0f64 / if temp[izero as usize] == 0.0 { 1.0 } else { temp[izero as usize] as f64 }) as f32;
            for it in 0..ntout as usize { temp[it] *= scale; }
        }
        tr.data[..ntout as usize].copy_from_slice(&temp[..ntout as usize]);
        tr.ns = ntout as u16;
        tr.f1 = ((-dt * ntout as f32) as f64 / 2.0) as f32;
        tr.delrt = 0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
