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

//! SUINTVEL - convert stacking velocity model to interval velocity model.
//! Port of su/main/data_conversion/suintvel.c (Revision 1.14).
//! Credits: CWP: Jack.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUINTVEL - convert stacking velocity model to interval velocity model\t",
    "\t\t\t\t\t\t\t\t\t",
    " suintvel vs= t0= outpar=/dev/tty\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t        \t",
    "\tvs=\tstacking velocities \t\t\t\t\t",
    "\tt0=\tnormal incidence times\t\t \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tmode=0\t\t\toutput h= v= ; =1 output v=  t= \t",
    "\toutpar=/dev/tty\t\toutput parameter file in the form:\t",
    "\t\t\t\th=layer thicknesses vector\t\t",
    "\t\t\t\tv=interval velocities vector\t\t",
    "\t\t\t\t....or ...\t\t\t\t",
    "\t\t\t\tt=vector of times from t0\t\t",
    "\t\t\t\tv=interval velocities vector\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "    suintvel vs=5000,5523,6339,7264 t0=.4,.8,1.125,1.425 outpar=intpar\t",
    "\t\t\t\t\t\t\t\t\t",
    "    suintvel par=stkpar outpar=intpar\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " If the file, stkpar, contains:\t\t\t\t\t",
    "    vs=5000,5523,6339,7264\t\t\t\t\t\t",
    "    t0=.4,.8,1.125,1.425\t\t\t\t\t\t",
    " then the two examples are equivalent.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: suintvel does not have standard su syntax since it does not\t",
    "      operate on seismic data.  Hence stdin and stdout are not used.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: may go away in favor of par program, velconv, by Dave\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut outpar = "/dev/tty".to_string();
    getparstring("outpar", &mut outpar);
    let outparfp = efopen(&outpar, "w");

    let n = countparval("t0");
    if n == 0 { err!("no t0's specified"); }
    let n = n as usize;
    let mut v = vec![0.0f32; n];
    let mut h = vec![0.0f32; n];
    let mut vs = vec![0.0f32; n];
    let mut t0 = vec![0.0f32; n];
    if getparfloat("t0", &mut t0) != n as i32 { err!("expected %d intervals", n as i32); }
    if getparfloat("vs", &mut vs) != n as i32 { err!("expected %d velocities", n as i32); }
    let mut mode = 0i32;
    getparint("mode", &mut mode);
    checkpars();

    for i in 0..n {
        if vs[i] <= 0.0 { err!("vs's must be positive: vs[%d] = %f", i as i32, vs[i]); }
        if t0[i] < 0.0 { err!("t0's must be positive: t0[%d] = %f", i as i32, t0[i]); }
    }
    h[0] = (0.5 * vs[0] as f64 * t0[0] as f64) as f32;
    v[0] = vs[0];
    for i in 1..n {
        let t2 = (0.5 * t0[i] as f64) as f32;
        let t1 = (0.5 * t0[i - 1] as f64) as f32;
        let v2 = vs[i];
        let v1 = vs[i - 1];
        let dt = t2 - t1;
        let dvh = v2 * v2 * t2 - v1 * v1 * t1;
        h[i] = ((dvh * dt) as f64).sqrt() as f32;
        v[i] = ((dvh / dt) as f64).sqrt() as f32;
    }
    if mode == 0 {
        fprintf!(outparfp, "h=");
        for i in 0..n - 1 { fprintf!(outparfp, "%g,", h[i]); }
        fprintf!(outparfp, "%g\n", h[n - 1]);
    } else {
        fprintf!(outparfp, "t=");
        for i in 0..n - 1 { fprintf!(outparfp, "%g,", t0[i]); }
        fprintf!(outparfp, "%g\n", t0[n - 1]);
    }
    fprintf!(outparfp, "v=");
    for i in 0..n - 1 { fprintf!(outparfp, "%g,", v[i]); }
    fprintf!(outparfp, "%g\n", v[n - 1]);
    exit(CWP_Exit());
}
