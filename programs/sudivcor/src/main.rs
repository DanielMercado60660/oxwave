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

//! SUDIVCOR - divergence (spreading) correction.
//! Port of su/main/amplitudes/sudivcor.c (Revision 1.16).
//! Credits: CWP: Jack K. Cohen, Francesca Fazarri.

use oxwave_core::prelude::*;
use oxwave_linear_interp::intlin;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUDIVCOR - Divergence (spreading) correction\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sudivcor <stdin >stdout  [optional parms]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " trms=0.0\ttimes corresponding to rms velocities in vrms\t\t",
    " vrms=1500.0\trms velocities corresponding to times in trms\t",
    " vfile=\tbinary (non-ascii) file containing velocities vrms(t)\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The trms, vrms arrays specify an rms velocity function of time.\t",
    " Linear interpolation and constant extrapolation is used to determine\t",
    " rms velocities at times not specified.  Values specified in trms\t",
    " must increase monotonically.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Alternatively, rms velocities may be stored in a binary file\t\t",
    " containing one velocity for every time sample.  If vfile is specified,",
    " then the trms and vrms arrays are ignored.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The time of the first sample is assumed to be constant, and is taken\t",
    " as the value of the first trace header field delrt. \t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let tmin = (tr.delrt as f64 / 1000.0) as f32;
    // C reads vt[1] even for a one-sample velocity file.
    let mut vt = vec![0.0f32; nt.max(2)];
    let mut vfile = String::new();
    let denom: f32;
    if getparstring("vfile", &mut vfile) == 0 {
        let mut ntrms = countparval("trms") as usize;
        let mut nvrms = countparval("vrms") as usize;
        if nvrms != ntrms { err!("number of trms and vrms must be equal"); }
        if ntrms == 0 { ntrms = 2; }
        // The C source writes element 1 even when only one pick was allocated.
        let mut trms = vec![0.0f32; ntrms.max(2)];
        let npar = getparfloat("trms", &mut trms);
        if npar == 0 { trms[0] = 0.0; trms[1] = dt; }
        if npar == 1 { trms[1] = trms[0] + dt; }
        if nvrms == 0 { nvrms = 2; }
        let mut vrms = vec![0.0f32; nvrms.max(2)];
        let npar = getparfloat("vrms", &mut vrms);
        if npar == 0 { vrms[0] = 1500.0; vrms[1] = 1500.0; }
        if npar == 1 { vrms[1] = vrms[0]; }
        for itrms in 1..ntrms {
            if trms[itrms] <= trms[itrms - 1] { err!("trms must increase monotonically"); }
        }
        for it in 0..nt {
            let t = tmin + it as f32 * dt;
            intlin(&trms[..ntrms], &vrms[..ntrms], vrms[0], vrms[ntrms - 1], &[t], &mut vt[it..it+1]);
        }
        denom = trms[1] * vrms[1] * vrms[1];
    } else {
        let fp = efopen(&vfile, "r");
        if efread(&mut vt[..nt], &fp) != nt {
            err!("cannot read %d velocities from file %s", nt, vfile);
        }
        denom = (tmin + dt) * vt[1] * vt[1];
    }
    checkpars();
    let mut divcor = vec![0.0f32; nt];
    for it in 0..nt {
        let t = tmin + it as f32 * dt;
        divcor[it] = t * vt[it] * vt[it] / denom;
    }
    loop {
        for it in 0..nt { tr.data[it] *= divcor[it]; }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
