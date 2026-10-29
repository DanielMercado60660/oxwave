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
//! SUTAPER - Taper the edge traces of a data panel to zero.
//! Port of su/main/tapering/sutaper.c (Revision 1.18).
//! Credits: CWP: Chris Liner, Jack K. Cohen; rewrite: Tagir Galikeev.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUTAPER - Taper the edge traces of a data panel to zero.\t",
    "\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " sutaper <stdin >stdout [optional parameters]\t\t  ",
    "\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t  ",
    " ntr=tr.ntr\tnumber of traces. If tr.ntr is not set, then\t",
    " \t\tntr is mandatory\t\t\t\t",
    " tr1=0\t number of traces to be tapered at beginning\t",
    " tr2=tr1\tnumber of traces to be tapered at end\t\t",
    " min=0.\t\tminimum amplitude factor of taper\t\t",
    " tbeg=0\t\tlength of taper (ms) at trace start\t\t",
    " tend=0\t\tlength of taper (ms) at trace end\t\t",
    " taper=1\ttaper type\t\t\t\t\t",
    "\t\t =1 linear (default)\t\t\t   ",
    "\t\t =2 sine\t\t\t\t\t",
    "\t\t =3 cosine\t\t\t\t\t",
    "\t\t =4 gaussian (+/-3.8)\t\t\t  ",
    "\t\t =5 gaussian (+/-2.0)\t\t\t  ",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    "   To eliminate the taper, choose tbeg=0. and tend=0. and tr1=0",
    "\t\t\t\t\t\t\t\t",
];

fn envelope(f: f32, kind: i32, min: f32, linear_weights: bool) -> f32 {
    match kind as i8 {
        1 => if linear_weights { min + (1.0f32 - min) * f } else { f },
        2 => (std::f64::consts::PI * f as f64 / 2.0).sin() as f32,
        3 => (0.5 * (1.0 - (std::f64::consts::PI * f as f64).cos())) as f32,
        4 | 5 => {
            let scale = if kind as i8 == 4 { 3.8090232 } else { 2.0 };
            let x = (scale * (1.0f32 - f) as f64) as f32;
            (-(x * x) as f64).exp() as f32
        }
        _ => err!(" taper ?!"),
    }
}

fn weights(tr1: i32, tr2: i32, min: f32, kind: i32, w: &mut [f32]) {
    if tr2 != 0 && tr1 != tr2 {
        for i in 0..=tr2 {
            let f = i as f32 / tr2 as f32;
            w[(1 + tr1 + i) as usize] = envelope(f, kind, min, true);
        }
    }
    if tr1 != 0 {
        for i in 0..=tr1 {
            let f = i as f32 / tr1 as f32;
            w[i as usize] = envelope(f, kind, min, true);
        }
    }
}

fn taper(t1: f32, t2: f32, kind: i32, t: f32, dt: f32, trace: &mut [f32]) {
    let nt = ((t / dt) as f64 + 1.0) as i32;
    let nt1 = ((t1 / dt) as f64 + 1.0) as i32;
    let nt2 = ((t2 / dt) as f64 + 1.0) as i32;
    if nt1 > 1 {
        for i in 0..nt1 {
            let f = i as f32 / nt1 as f32;
            trace[i as usize] *= envelope(f, kind, 0.0, false);
        }
    }
    if nt2 > 1 {
        for i in 0..nt2 {
            let f = i as f32 / nt2 as f32;
            trace[(nt-i) as usize] *= envelope(f, kind, 0.0, false);
        }
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut verbose: i16 = 0;
    getparshort("verbose", &mut verbose);
    let mut min = 0.0f32;
    getparfloat("min", &mut min);
    if min > 1.0 { err!("min must be less than 1"); }
    let mut t1 = 0.0f32;
    let mut t2 = 0.0f32;
    let mut ttaper = 1;
    getparfloat("tbeg", &mut t1);
    getparfloat("tend", &mut t2);
    getparint("taper", &mut ttaper);
    let mut tr1 = 0;
    getparint("tr1", &mut tr1);
    let mut tr2 = tr1;
    getparint("tr2", &mut tr2);
    let ntap = if tr1 != tr2 { tr1 + tr2 + 2 } else { tr1 + 1 };
    let mut taperv = vec![0.0f32; ntap as usize];
    weights(tr1, tr2, min, ttaper, &mut taperv);
    if verbose != 0 {
        fprintf!(stderr(), "print %i taper weights: \n", ntap);
        for i in 0..ntap { fprintf!(stderr(), "Taper %i = %g \n", i, taperv[i as usize]); }
    }
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000.0) as f32;
    let mut ntr = tr.ntr;
    let mut have_ntr = ntr != 0;
    if getparint("n2", &mut ntr) != 0 || getparint("ntr", &mut ntr) != 0 { have_ntr = true; }
    if !have_ntr { err!("ntr neither set nor getparred"); }
    checkpars();
    let tlen = (nt - 1) as f32 * dt;
    if t1 + t2 > tlen {
        err!("sum of tapers tbeg=%f, tend=%f exceeds trace length(%f ms)", t1, t2, tlen);
    }
    let mut itr = 0;
    let mut pos: i32 = 0;
    loop {
        let mut fac = 1.0f32;
        if itr < tr1 { fac *= taperv[pos as usize]; pos += 1; }
        if itr == tr1 { pos += ntap - tr1 - 2; }
        itr += 1;
        if itr > ntr - tr2 { fac *= taperv[pos as usize]; pos -= 1; }
        if verbose != 0 && tr1 + tr2 != 0 && fac < 1.0 {
            fprintf!(stderr(), "trace %i factor %g \n", itr, fac);
        }
        if tr1 != 0 || tr2 != 0 {
            for i in 0..nt as usize { tr.data[i] *= fac; }
        }
        if t1 != 0.0 || t2 != 0.0 { taper(t1, t2, ttaper, tlen, dt, &mut tr.data); }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
