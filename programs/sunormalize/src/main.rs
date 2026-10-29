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

//! SUNORMALIZE - Trace normalization. Port of su/main/amplitudes/sunormalize.c (Revision 1.8).
//! Credits: Ramone Carbonell; Nils Maercklin.

/* this was originally the source code for div
 * Copyright (c) 1990, 1993
 *    The Regents of the University of California.  All rights reserved.
 *
 * This code is derived from software contributed to Berkeley by
 * Chris Torek.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the University of
 *    California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

use oxwave_core::prelude::*;
use oxwave_sort::qksort;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                ",
    " SUNORMALIZE - Trace NORMALIZation by rms, max, or median       ",
    "               or median balancing                              ",
    "                                                                ",
    "   sunormalize <stdin >stdout t0=0 t1=TMAX norm=rms             ",
    "                                                                ",
    " Required parameters:                                           ",
    "    dt=tr.dt    if not set in header, dt is mandatory           ",
    "    ns=tr.ns    if not set in header, ns is mandatory           ",
    "                                                                ",
    " Optional parameters:                                           ",
    "    norm=rms    type of norm rms, max, med , balmed             ",
    "    t0=0.0      startimg time for window                        ",
    "    t1=TMAX     ending time for window                          ",
    "                                                                ",
    " Notes:                                                         ",
    " Traces are divided by either the root mean squared amplitude,  ",
    " trace maximum, or the median value. The option \"balmed\" is   ",
    " median balancing which is a shift of the amplitudes by the\t ",
    " median value of the amplitudes.\t\t\t\t ",
    "                                                                ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut rms: f32 = 0.0;
    initargs(&argv);
    requestdoc(1, SDOC);
    if fgettr(&stdin(), &mut tr) == 0 { err!("can't read first trace"); }
    let mut dt: f32 = 0.0;
    if tr.dt == 0 { mustgetparfloat("dt", &mut dt); }
    else { dt = (tr.dt as f64 / 1000000.0) as f32; }
    let mut ns: i32 = 0;
    if tr.ns == 0 { mustgetparint("ns", &mut ns); }
    else { ns = tr.ns as i32; }
    let mut t0: f32 = 0.0;
    let mut t1: f32 = 0.0;
    if getparfloat("t0", &mut t0) == 0 { t0 = 0.0; }
    if getparfloat("t1", &mut t1) == 0 { t1 = ns as f32 * dt; }
    let it0 = (t0 / dt) as i32;
    let it1 = (t1 / dt) as i32;
    let n = it1.wrapping_sub(it0);
    // The C even median can read one float past the allocated array.
    if n < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut z = vec![0.0f32; n as usize + 2];
    let mut norm = String::from("rms");
    getparstring("norm", &mut norm);
    checkpars();
    let inorm = match norm.as_str() {
        "rms" => 0, "max" => 1, "med" | "median" => 2, "balmed" => 3,
        _ => err!("unknown operation=\"%s\", see self-doc", norm),
    };
    loop {
        match inorm {
            0 => {
                for j in 0..n { let x = tr.data[(it0 + j) as usize]; rms += x * x; }
                rms = ((rms / n as f32) as f64).sqrt() as f32;
                if rms == 0.0 { rms = 1.0; }
                for i in 0..ns { tr.data[i as usize] /= rms; }
            }
            1 => {
                rms = tr.data[it0 as usize].abs();
                for j in 0..n {
                    let value = tr.data[(it0 + j) as usize].abs();
                    if value > rms { rms = value; }
                }
                if rms == 0.0 { rms = 1.0; }
                for i in 0..ns { tr.data[i as usize] /= rms; }
            }
            2 | 3 => {
                for j in 0..n { z[j as usize] = tr.data[(it0 + j) as usize]; }
                qksort(&mut z[..n as usize]);
                let quot = n / 2;
                if n % 2 != 0 { rms = z[quot as usize]; }
                else { rms = ((z[quot as usize] + z[quot as usize + 1]) as f64 / 2.0) as f32; }
                if rms == 0.0 { rms = 1.0; }
                for i in 0..ns {
                    if inorm == 2 { tr.data[i as usize] /= rms; }
                    else { tr.data[i as usize] -= rms; }
                }
            }
            _ => err!("mysterious operation=\"%s\"", norm),
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
