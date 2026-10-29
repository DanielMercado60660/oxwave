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

//! ROSSLER - compute the Rossler attractor.
//! Port of par/main/ordinary_differential_equations/rossler.c.
//! Credits (from the C source): CWP: Aug 2013: John Stockwell.

use oxwave_core::prelude::*;
use oxwave_f_rke::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " ROSSLER - compute the ROSSLER attractor\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "  rossler > [stdout]\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required Parameters: none\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " a=.2\t\t\tparameter for rossler equations\t\t",
    " b=.2\t\t\tparameter for rossler equations\t\t",
    " c=5.7\t\t\tparameter for rossler equations\t\t",
    " y0=1.0\t\tinitial value of y[0]\t\t\t",
    " y1=-1.0\t\tinitial value of y[1]\t\t\t",
    " y2=1.0\t\tinitial value of y[2]\t\t\t",
    " h=.01\t\t\tincrement in time\t\t\t",
    " tol=1.e-08\t\terror tolerance\t\t\t\t",
    " stepmax=500\t\tmaximum number of steps to compute\t",
    " mode=xy\t\txy-pairs, =yz yz-pairs, =xz xz-pairs,\t",
    "\t\t\t=xyz xyz-triplet, =x only, =y only, =z only",
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
    " rossler stepmax=1000 mode=xy | xgraph n=1000\t&\t\t",
    " rossler stepmax=1000 mode=yz | xgraph n=1000\t&\t\t",
    " rossler stepmax=1000 mode=xz | xgraph n=1000\t&\t\t",
    "\t\t\t\t\t\t\t\t",
    " rossler stepmax=1000 mode=x | suaddhead ns=1000 | suxwigb &\t",
    " rossler stepmax=1000 mode=y | suaddhead ns=1000 | suxwigb &\t",
    " rossler stepmax=1000 mode=z | suaddhead ns=1000 | suxwigb &\t",
    "\t\t\t\t\t\t\t\t",
    " GNUPLOT 3D plot example:                                      ",
    " rossler stepmax=2000 mode=xyz > rossler.bin                   ",
    " ...when you run gnuplot type the following command ...        ",
    " splot \"rossler.bin\" binary record=2000:2000:2000 with points pointsize .1 ",
    "                                                               ",
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
    let mut verbose: i32 = 0; /* verbose flag =1 chatty, =0 silent */
    let mut stepmax: i32 = 0; /* maximum number of steps */
    let mut h: f64 = 0.001; /* time increment */
    let mut tol: f64 = 0.0;
    let mut y = [0.0f64; 3]; /* dependent variable of ODE system */
    let out_file = stdout();
    let mut mode = String::from("xy");
    let mut imode = SXY;

    /* Hook up getpar */
    initargs(&argv);
    requestdoc(0, SDOC);

    /* Prevent floats from dumping on screen */
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

    /* Get parameters */
    if getparint("stepmax", &mut stepmax) == 0 {
        stepmax = 500;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    if getpardouble("y0", &mut y[0]) == 0 {
        y[0] = 1.0;
    }
    if getpardouble("y1", &mut y[1]) == 0 {
        y[1] = -1.0;
    }
    if getpardouble("y2", &mut y[2]) == 0 {
        y[2] = 1.0;
    }
    if getpardouble("h", &mut h) == 0 {
        h = 0.01;
    }
    if getpardouble("tol", &mut tol) == 0 {
        // RKE_ERR_BIAS_INIT from rke.h.
        tol = 0.00000001;
    }

    /* Get output mode */
    getparstring("mode", &mut mode);
    if mode == "yz" {
        imode = SYZ;
    } else if mode == "xz" {
        imode = SXZ;
    } else if mode == "xyz" {
        imode = SXYZ;
    } else if mode == "x" {
        imode = SX;
    } else if mode == "y" {
        imode = SY;
    } else if mode == "z" {
        imode = SZ;
    } else if mode != "xy" {
        err!("unknown operation=\"%s\", see self-doc", mode);
    }

    /* Initialize Runge-Kutta-England routines and set tolerance */
    let mut p = rke_init(3, rossler_equations);
    p.error_bias = tol;

    for i in 0..stepmax {
        let number = 3;
        let mut yout = [0.0f32; 3];
        let mut t = i as f64 * h;
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

        /* Convert doubles in y to floats in yout and write out */
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

        /* Run the Runge-Kutta-England solver */
        rke_solve(&mut p, &mut t, &mut y, aimed_t);
    }

    /* End the session with rke */
    rke_term(p);
    exit(CWP_Exit());
}

/* The autonomous system of ODEs describing the Rossler attractor. */
fn rossler_equations(_t: f64, y: &mut [f64], yprime: &mut [f64]) -> i32 {
    let mut a: f64 = 0.0;
    let mut b: f64 = 0.0;
    let mut c: f64 = 0.0;

    if getpardouble("a", &mut a) == 0 {
        a = 0.2;
    }
    if getpardouble("b", &mut b) == 0 {
        b = 0.2;
    }
    if getpardouble("c", &mut c) == 0 {
        c = 5.7;
    }

    yprime[0] = -y[1] - y[2];
    yprime[1] = y[0] + a * y[1];
    yprime[2] = b + y[2] * (y[0] - c);
    1
}
