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
//! SUDIVSTACK - diversity stacking. Port of su/main/stacking/sudivstack.c (Revision 1.12).
//! Credits: Mary Palen-Murphy; key option: Nils Maercklin.

use oxwave_core::prelude::*;
use oxwave_linear_interp::intlin;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUDIVSTACK -  Diversity Stacking using either average power or peak   ",
    "               power within windows                                    ",
    "                                                                       ",
    "                                                                       ",
    " Required parameters:                                                  ",
    "    none                                                               ",
    "                                                                       ",
    " Optional parameters:                                                  ",
    "    key=tracf        key header word to stack on                       ",
    "    winlen=0.064     window length in seconds.                         ",
    "                     typical choices: 0.064, 0.128, 0.256,             ",
    "                                      0.512, 1.024, 2.048, 4.096       ",
    "                     if not specified the entire trace is used         ",
    "                                                                       ",
    "    peak=1           peak power option default is average power        ",
    "                                                                       ",
    " Notes:                                                                ",
    "    Diversity stacking is a noise reduction technique used in the      ",
    "    summation of duplicate data. Each trace is scaled by the inverse   ",
    "    of its average power prior to stacking.  The composite trace is    ",
    "    then renormalized by dividing by the sum of the scalers used.      ",
    "                                                                       ",
    "    This program stacks adjacent traces having the same key header     ",
    "    word, which can be specified by the key parameter. The default     ",
    "    is \"tracf\" (trace number within field record).                   ",
    "    For more information on key header words, type \"sukeyword -o\".   ",
    "                                                                       ",
    " Examples:                                                             ",
    "    For duplicate field shot records:                                  ",
    "        susort < field.data tracf | sudivstack > stack.data            ",
    "    For CDP ordered data:                                              ",
    "        sudivstack < cdp.data key=cdp > stack.data                     ",
    "                                                                       ",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut intrace = Trace::new();
    if gettr(&mut intrace) == 0 { err!("can not get first trace"); }
    let nt = intrace.ns as usize;
    let dt = (intrace.dt as f64 / 1000000.0) as f32;
    let entire = nt as f32 * dt - dt;
    let mut key = String::new();
    if getparstring("key", &mut key) == 0 { key = "tracf".to_string(); }
    let mut winlen = 0.0f32;
    if getparfloat("winlen", &mut winlen) == 0 { winlen = entire; }
    let mut peak = 0i32;
    getparint("peak", &mut peak);
    checkpars();
    let mut ntwin = ((winlen / dt) as i32).wrapping_add(1);
    if ntwin > nt as i32 { ntwin = nt as i32; }
    let type_ = hdtype(&key);
    let index = getindex(&key);
    let mut intpapw = vec![0.0f32; nt];
    let mut sumscale = vec![0.0f32; nt];
    let mut sumdata = vec![0.0f32; nt];
    if ntwin < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    // A zero-sized window makes the original's window loop nonterminating.
    if ntwin == 0 { err!("window length gives zero samples"); }
    let mut xout = vec![0.0f32; ntwin as usize];
    let mut data = true;
    let mut newtracl = 1i32;
    while data {
        sumscale.fill(0.0);
        sumdata.fill(0.0);
        let mut outtrace = intrace.clone();
        let mut val = Value::default();
        gethval(&intrace, index, &mut val);
        let trnum = vtoi(type_, val);
        loop {
            let nw = ntwin as usize;
            let mut avepwa;
            if peak != 0 {
                let mut maxval = 0.0f32;
                for i in 0..nw { let x = intrace.data[i]; if x*x > maxval { maxval = x*x; } }
                avepwa = maxval;
            } else {
                let mut totalpwa = 0.0f32;
                for i in 0..nw { let x = intrace.data[i]; totalpwa += x*x; }
                avepwa = totalpwa / ntwin as f32;
            }
            for i in 0..nw { intpapw[i] = avepwa; }
            let mut j = 1usize;
            while nw * (j + 1) < nt {
                let avepwb;
                if peak != 0 {
                    let mut maxval = 0.0f32;
                    for i in nw*j..nw*(j+1) { let x = intrace.data[i]; if x*x > maxval { maxval = x*x; } }
                    avepwb = maxval;
                } else {
                    let mut totalpwb = 0.0f32;
                    for i in nw*j..nw*(j+1) { let x = intrace.data[i]; totalpwb += x*x; }
                    avepwb = totalpwb / ntwin as f32;
                }
                let xin = [(nw*j-1) as f32, (nw*(j+1)-1) as f32];
                let yin = [avepwa, avepwb];
                for i in nw*j..nw*(j+1) { xout[i-nw*j] = i as f32; }
                intlin(&xin, &yin, avepwa, avepwb, &xout[..nw], &mut intpapw[nw*j..nw*(j+1)]);
                avepwa = avepwb;
                j += 1;
            }
            if nw < nt {
                let avepwb;
                if peak != 0 {
                    let mut maxval = 0.0f32;
                    for i in nw*j..nt { let x = intrace.data[i]; if x*x > maxval { maxval = x*x; } }
                    avepwb = maxval;
                } else {
                    let mut count = 0i32;
                    let mut totalpwb = 0.0f32;
                    for i in nw*j..nt { let x = intrace.data[i]; totalpwb += x*x; count += 1; }
                    avepwb = totalpwb / count as f32;
                }
                let xin = [(nw*j-1) as f32, (nt-1) as f32];
                let yin = [avepwa, avepwb];
                let count = nt - nw*j;
                for i in nw*j..nt { xout[i-nw*j] = i as f32; }
                intlin(&xin, &yin, avepwa, avepwb, &xout[..count], &mut intpapw[nw*j..nt]);
            }
            for i in 0..nt {
                if intpapw[i] != 0.0 {
                    sumdata[i] += intrace.data[i] / intpapw[i];
                    sumscale[i] += 1.0f32 / intpapw[i];
                }
            }
            if gettr(&mut intrace) == 0 { data = false; break; }
            gethval(&intrace, index, &mut val);
            if vtoi(type_, val) != trnum { break; }
        }
        for i in 0..nt {
            if sumscale[i] == 0.0 { sumscale[i] = 1.0; }
            outtrace.data[i] = sumdata[i] / sumscale[i];
        }
        outtrace.tracl = newtracl;
        newtracl = newtracl.wrapping_add(1);
        puttr(&mut outtrace);
    }
    exit(CWP_Exit());
}
