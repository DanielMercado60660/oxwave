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

//! VERHULST - solve the Verhulst logistic equation.
//! Port of par/main/ordinary_differential_equations/verhulst.c.
//! Credits (from the C source): CWP: Aug 2009: John Stockwell.

use oxwave_core::prelude::*;
use oxwave_f_rke::*;
use std::cell::Cell;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " VERHULST - solve the VERHULST logistic equation\t\t",
    "\t\t\t\t\t\t\t\t",
    "  verhulst > [stdout]\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required Parameters: none\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " a1=1.0\t\tparameter for verhulst equation\t\t",
    " a2=2000\t\tparameter for verhulst equation\t\t",
    " y0=10\t\t\tinitial value of y[0]\t\t\t",
    " h=.01\t\t\tincrement in time\t\t\t",
    " tol=1.e-08\t\terror tolerance\t\t\t\t",
    " stepmax=2000\t\tmaximum number of steps to compute\t",
    " mode=x\t\txy-pairs, =yz yz-pairs, =xz xz-pairs,\t",
    "\t\t\t=xyz xyz-triplet, =x only, =y only, =z only",
    " toutfile=\t\tfile of time values\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " This program is really just a demo showing how to use the \t",
    " differential equation solver rke_solve written by Francois \t",
    " Pinard, based on a modified form of the 4th order Runge-Kutta ",
    " method, which employs the error checking method of R. England ",
    " 1969.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " The output consists of unformated C-style binary floats, of\t",
    " either pairs or triplets as specified by the \"mode\" paramerter.",
    "\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t",
    " x is the population:\t\t\t\t\t\t",
    "    verhulst stepmax=2000 mode=x |\t\t\t\t",
    "\txgraph n=2000 nplot=1 d1=0.01 style=normal &\t\t",
    " y is dx/dt, the rate of growth of the population:\t\t",
    "    verhulst stepmax=2000 mode=y |\t\t\t\t",
    "\txgraph n=2000 nplot=1 d1=0.01 style=normal &\t\t",
    " show both:\t\t\t\t\t\t\t",
    "    verhulst stepmax=2000 mode=x > population.bin\t\t",
    "    verhulst stepmax=2000 mode=y > rate.bin\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " cat population.bin rate.bin |\t\t\t\t\t",
    "     xplot n=2000 nplot=2 d1=.01 style=normal &\t\t",
    " \t\t\t\t\t\t\t\t",
    " In the Verhulst equation, a1 is the reproduction rate and\t",
    " a2 is the carrying capacity\t\t\t\t\t",
    " \tx'(t) = a1 * x * ( 1 - x/a2 )\t\t\t \t",
    " Related program: logisticfit \t\t\t\t\t",
];

const SXY: i32 = 0;
const SXZ: i32 = 1;
const SYZ: i32 = 2;
const SXYZ: i32 = 3;
const SX: i32 = 4;
const SY: i32 = 5;
const SZ: i32 = 6;

fn main() {
    let argv = args();
    let mut verbose: i32 = 0;
    let mut stepmax: i32 = 0;
    let mut h: f64 = 0.001;
    let mut tol: f64 = 0.0;
    let mut y = [0.0f64; 3];
    let out_file = stdout();
    let mut mode = String::from("x");
    let mut imode = SXY;
    let mut toutfile = String::new();
    let mut toutfp = None;
    let mut toutt = [0.0f32; 1];

    initargs(&argv);
    requestdoc(0, SDOC);

    // Prevent floats from dumping on screen.
    match filestat(STDOUT) {
        FileType::BadFileType => {
            warn!("stdout is illegal filetype");
            pagedoc(SDOC);
        }
        FileType::Tty => {
            warn!("stdout can't be tty");
            pagedoc(SDOC);
        }
        _ => {}
    }

    if getparint("stepmax", &mut stepmax) == 0 {
        stepmax = 2000;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    if getpardouble("y0", &mut y[..]) == 0 {
        y[0] = 10.0;
    }
    if getpardouble("h", &mut h) == 0 {
        h = 0.01;
    }
    if getpardouble("tol", &mut tol) == 0 {
        tol = 0.00000001; // RKE_ERR_BIAS_INIT
    }

    getparstring("mode", &mut mode);
    if mode == "yz" {
        imode = SYZ;
    } else if mode == "xz" {
        imode = SXZ;
    } else if mode == "x" {
        imode = SX;
    } else if mode == "y" {
        imode = SY;
    } else if mode == "z" {
        imode = SZ;
    } else if mode != "xy" {
        err!("unknown operation=\"%s\", see self-doc", mode);
    }

    getparstring("toutfile", &mut toutfile);
    if !toutfile.is_empty() {
        toutfp = Some(match fopen(&toutfile, "w") {
            Some(fp) => fp,
            None => err!("cannot open toutfile=%s\n", toutfile),
        });
    }

    // C's callback writes y[1] despite integrating just one equation. Only
    // calls on the main state affect the output; writes to the oversized C
    // solver work arrays are never read. The Rust callback takes a mutable
    // one-element slice, so record that side effect separately and copy it
    // back after each solve. Pointer equality distinguishes main/work states
    // without dereferencing raw pointers or modifying the solver library.
    let main_y = y.as_ptr();
    let growth_rate = Cell::new(y[1]);
    let mut p = rke_init(1, |t, state: &mut [f64], yprime: &mut [f64]| {
        let result = verhulst_equation(t, state, yprime);
        if state.as_ptr() == main_y {
            growth_rate.set(yprime[0]);
        }
        result
    });
    p.error_bias = tol;

    for i in 0..stepmax {
        let number = 3;
        let mut yout = [0.0f32; 3];
        let mut t = i as f64 * h;

        if let Some(fp) = &toutfp {
            toutt[0] = t as f32;
            efwrite(&toutt, fp);
        }

        let aimed_t = t + h;
        if verbose != 0 {
            warn!(
                "using %3d accepted and %3d rejected steps",
                p.accepted_steps, p.rejected_steps
            );
            if verbose != 0 {
                warn!("error tolerance = %10.24f", p.error_bias);
            }
        }

        for j in 0..number {
            yout[j] = y[j] as f32;
        }

        let mut tmpout = [0.0f32; 2];
        match imode {
            SXY => {
                tmpout[0] = yout[0];
                tmpout[1] = yout[1];
                efwrite(&tmpout, &out_file);
            }
            SYZ => {
                tmpout[0] = yout[1];
                tmpout[1] = yout[2];
                efwrite(&tmpout, &out_file);
            }
            SXZ => {
                tmpout[0] = yout[0];
                tmpout[1] = yout[2];
                efwrite(&tmpout, &out_file);
            }
            SXYZ => {
                efwrite(&yout, &out_file);
            }
            SX => {
                tmpout[0] = yout[0];
                efwrite(&tmpout[..1], &out_file);
            }
            SY => {
                tmpout[0] = yout[1];
                efwrite(&tmpout[..1], &out_file);
            }
            SZ => {
                tmpout[0] = yout[2];
                efwrite(&tmpout[..1], &out_file);
            }
            _ => err!("mysterious operation=\"%s\"", mode),
        }

        rke_solve(&mut p, &mut t, &mut y, aimed_t);
        y[1] = growth_rate.get();
    }

    rke_term(p);
    if let Some(fp) = &toutfp {
        efclose(fp);
    }
    exit(CWP_Exit());
}

// The system of ODEs describing population growth by the logistic equation.
// This is autonomous: time does not enter the derivative.
fn verhulst_equation(_t: f64, y: &[f64], yprime: &mut [f64]) -> i32 {
    let mut a1: f64 = 0.0;
    let mut a2: f64 = 0.0;
    if getpardouble("a1", &mut a1) == 0 {
        a1 = 1.0;
    }
    if getpardouble("a2", &mut a2) == 0 {
        a2 = 2000.0;
    }
    yprime[0] = a1 * y[0] * (1.0 - y[0] / a2);
    // y[1] = yprime[0] is preserved by the callback bridge in main.
    1
}
