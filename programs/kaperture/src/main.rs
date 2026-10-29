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


//! KAPERTURE - generate the k domain of a line scatterer for a seismic array.
//! Port of par/main/apertures/kaperture.c (Revision 1.10).
//! Credits: CWP: Jack K. Cohen, Sebastien Geoltrain, Norm Bleistein.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t \t\t\t\t\t\t\t\t",
    " KAPERTURE - generate the k domain of a line scatterer for a seismic array",
    " \t      \t\t\t\t\t\t\t\t",
    " kaperture [optional parameters] >stdout \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters\t\t\t\t\t\t\t",
    " \tx0=1000\t\tpoint scatterer location\t\t\t",
    " \tz0=1000\t\tpoint scatterer location\t\t\t",
    " \tnshot=1\t\tnumber of shots\t\t\t\t\t",
    " \tsxmin=0\t\tfirst shot location\t\t\t\t",
    " \tszmin=0\t\tfirst shot location\t\t\t\t",
    " \tdsx=100\t\tx-steps in shot location\t\t\t",
    " \tdsz=0\t\tz-steps in shot location\t\t\t",
    " \tngeo=1\t\tnumber of receivers\t\t\t\t",
    " \tgxmin=0\t\tfirst receiver location\t\t\t\t",
    " \tgzmin=0\t\tfirst receiver location\t\t\t\t",
    " \tdgx=100\t\tx-steps in receiver location\t\t\t",
    " \tdgz=0\t\tz-steps in receiver location\t\t\t",
    " \tfnyq=125\tNyquist frequency  (Hz)\t\t\t\t",
    " \tfmax=125\tmaximum frequency  (Hz)\t\t\t\t",
    " \tfmin=5\t\tminimum frequency  (Hz)\t\t\t\t",
    " \tnfreq=2\t\tnumber of frequencies   \t\t\t",
    " \tboth=0\t\t= 1 gives negative freqs too\t\t\t",
    " \tnstep=60\tpoints on Nyquist circle\t\t\t",
    " \tc=5000\t\tspeed\t\t\t\t\t\t",
    " \toutpar=/dev/tty output parameter file, contains:\t\t",
    " \t\t\t\txmin, xmax, ymin, ymax \t\t\t",
    " \t\t\t\tand npairs (needed for psgraph or xgraph)",
    " \t\t\tother choices for outpar are: /dev/tty,\t\t",
    " \t\t\t/dev/stderr, or a name of a disk file\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    "       nfreq=1 produces fmin\t\t\t\t\t\t",
    "       nstep=0 suppresses the Nyquist circle\t\t\t\t",
    " \t\t\t\tand npairs\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Default case: both=0 nfreq=2\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \tkaperture nshot=NSHOT ngeo=NGEO nstep=NSTEP |\t\t\t",
    " \tpsgraph\tn=NPAIRS,NSTEP mark=1,0 marksize=1,0 linewidth=0,1 |...",
    " \t\tWHERE: NPAIRS=NSHOT*NGEO\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Other cases: \t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " both=0 nfreq=NFREQ > 2\t\t\t\t\t\t",
    " \tkaperture both=0 nfreq=NFREQ nshot=NSHOT ngeo=NGEO nstep=NSTEP |",
    " \tpsgraph\tn=NPAIRS,NSTEP mark=1,0 marksize=1,0 linewidth=0,1 |...\t",
    " \t\tWHERE: NPAIRS=NFREQ*NSHOT*NGEO\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " both=1 nfreq=NFREQ > 2\t\t\t\t\t\t",
    " \tkaperture both=1 nfreq=NFREQ nshot=NSHOT ngeo=NGEO nstep=NSTEP |",
    " \tpsgraph\tn=NPAIRS,NSTEP mark=1,0 marksize=1,0 linewidth=0,1 |...\t",
    " \t\t WHERE: NPAIRS=NFREQ*NSHOT*NGEO*2\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " When in doubt to the size of NPAIRS, redirect output of kaperture to\t",
    " /dev/tty the first time to get npairs=:\t\t\t\t",
    "\t\t kaperture [optional parameters] > /dev/tty\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    let (mut nfreq, mut nshot, mut ngeo, mut nstep, mut both) = (1, 1, 1, 60, 0);
    let (mut c, mut x0, mut z0) = (5000.0f32, 1000.0f32, 1000.0f32);
    let (mut sxmin, mut szmin, mut gxmin, mut gzmin) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let (mut dsx, mut dsz, mut dgx, mut dgz) = (100.0f32, 0.0f32, 100.0f32, 0.0f32);
    let (mut fnyq, mut fmax, mut fmin) = (125.0f32, 125.0f32, 125.0f32);
    getparint("nfreq", &mut nfreq);
    getparint("nshot", &mut nshot);
    getparint("ngeo", &mut ngeo);
    getparfloat("c", &mut c);
    getparfloat("x0", &mut x0);
    getparfloat("z0", &mut z0);
    getparfloat("sxmin", &mut sxmin);
    getparfloat("szmin", &mut szmin);
    getparfloat("gxmin", &mut gxmin);
    getparfloat("gzmin", &mut gzmin);
    getparfloat("dsx", &mut dsx);
    getparfloat("dsz", &mut dsz);
    getparfloat("dgx", &mut dgx);
    getparfloat("dgz", &mut dgz);
    getparfloat("fnyq", &mut fnyq);
    getparfloat("fmax", &mut fmax);
    getparfloat("fmin", &mut fmin);
    getparint("nstep", &mut nstep);
    getparint("both", &mut both);
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".into(); }
    checkpars();
    let outparfp = efopen(&outpar, "w");

    let npaths = nshot.wrapping_mul(ngeo);
    let npoints = if both != 0 { 2i32.wrapping_mul(npaths) } else { npaths };
    let mut x = Vec::<f32>::new();
    let mut y = Vec::<f32>::new();
    if x.try_reserve_exact(npoints as usize).is_err() || y.try_reserve_exact(npoints as usize).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    x.resize(npoints as usize, 0.0);
    y.resize(npoints as usize, 0.0);
    for ipoint in 0..npaths {
        let s = (ipoint as f32 / ngeo as f32) as i32;
        let g = (ipoint as f32 - s.wrapping_mul(ngeo) as f32) as i32;
        let sx = x0 - (sxmin + s as f32 * dsx);
        let sz = z0 - (szmin + s as f32 * dsz);
        let gx = x0 - (gxmin + g as f32 * dgx);
        let gz = z0 - (gzmin + g as f32 * dgz);
        let mut rs = ((sx * sx + sz * sz) as f64).sqrt() as f32;
        if rs as f64 <= 1.0e-20 { rs = 1.0e-20; }
        let mut rg = ((gx * gx + gz * gz) as f64).sqrt() as f32;
        if rg as f64 <= 1.0e-20 { rg = 1.0e-20; }
        x[ipoint as usize] = sx / rs + gx / rg;
        y[ipoint as usize] = -(sz / rs + gz / rg);
        if both != 0 {
            x[(npaths + ipoint) as usize] = -x[ipoint as usize];
            y[(npaths + ipoint) as usize] = -y[ipoint as usize];
        }
    }
    let df = if nfreq == 1 { 0.0 } else { ((fmax - fmin) as f64 / (nfreq as f32 as f64 - 1.0)) as f32 };
    let mut npairs = 0i32;
    for f in 0..nfreq {
        let freq = fmin + f as f32 * df;
        let kscale = (6.28318530717959f64 * freq as f64 / c as f64) as f32;
        for ipoint in 0..npoints {
            npairs = npairs.wrapping_add(1);
            let tmpx = kscale * x[ipoint as usize];
            let tmpy = kscale * y[ipoint as usize];
            efwrite(&[tmpx], &stdout());
            efwrite(&[tmpy], &stdout());
        }
    }
    let knyqscale = (12.5663706143592f64 * fnyq as f64 / c as f64) as f32;
    let xmin = (-1.2f64 * knyqscale as f64) as f32;
    let xmax = (1.2f64 * knyqscale as f64) as f32;
    let ymin = xmin;
    let ymax = xmax;
    if nstep != 0 {
        for iphi in 0..nstep {
            let phi = (iphi as f64 * 6.28318530717959f64 / nstep as f64) as f32;
            let tmpx = ((phi as f64).cos() * knyqscale as f64) as f32;
            let tmpy = ((phi as f64).sin() * knyqscale as f64) as f32;
            efwrite(&[tmpx], &stdout());
            efwrite(&[tmpy], &stdout());
        }
    }
    fprintf!(outparfp, "xmin=%f xmax=%f ymin=%f ymax=%f\n", xmin, xmax, ymin, ymax);
    fprintf!(outparfp, "npairs=%d\n", npairs);
    exit(CWP_Exit());
}
