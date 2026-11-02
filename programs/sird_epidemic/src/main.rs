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

//! SIRD_EPIDEMIC - the SIRD and SIRDS EPIDEMIC models.
//! Port of par/main/ordinary_differential_equations/sird_epidemic.c (Revision 1.4).
//! Credits (from the C source): Author: April 2020: John Stockwell.

use oxwave_core::prelude::*;
use oxwave_f_rke::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SIRD_EPIDEMIC - the SIRD and SIRDS EPIDEMIC models \t\t",
    "\t\t\t\t\t\t\t\t",
    "  sirdepidemic > [stdout]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required Parameters: none\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " normalize=1\t\tnormalize S,I by N; =0 don't normalize\t",
    " scale=0\t\tdon't scale output\t\t\t",
    "\t\t\t=1 scale output S,I,R by N\t\t",
    "\t\t\t=2 scale output S,I,R by s0\t\t",
    " N=1000\t\ttotal population\t\t\t",
    " s0=N\t\t\tinitial number of susceptibles\t\t",
    " i0=1\t\t\tinitial number of infectives\t\t",
    " r0=0\t\t\tinitial number of removed (usually 0)\t",
    "\t \t\t(not the basic reproducion rate R0)\t",
    " d0=0\t\t\tinitial number of dead (usually 0)\t",
    " \t\t\t\t\t\t\t\t",
    " k=.5\t\t\ttransmission rate\t\t\t",
    " b=.3333\t\tremoval rate = death + recovery rates\t",
    " nu=0.01\t\tmortality rate \t\t\t\t",
    " \t\t\t(0.01 = 1% mortality)\t\t\t",
    " gamma=0\t\tno social distancing; gamma > 0 social\t",
    " \t\t\tdistancing (Baker 2020 model)\t\t",
    " h=1\t\t\tincrement in time\t\t\t",
    " tol=1.e-08\t\terror tolerance\t\t\t\t",
    " stepmax=40\t\tmaximum number of steps to compute\t",
    " mode=SIRD\t\tS, I, R, D successively output\t\t",
    "\t\t\t=S only, =I only, =R only, =D only\t",
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
    " About the SIRD epidemic model:\t\t\t\t",
    " R0 = number of new infections per single infected host  \t",
    "  1 < R0 < 1.5 for influenza (2.2 to 2.7 for Covid-19)\t\t",
    "  b, k, s0, and R0 are related via\t\t\t\t",
    "  k = b*s0/R0. When s is normalized by N, k=b/R0\t\t",
    "  It is often easier to determine the recovery rate k (in units",
    "  of h and to determine reasonable estimate of s0 and of R0 \t",
    "  and to calculate the infection rate b = k*R0/s0, which becomes",
    "  b = k*R0 for s normalized by N.\t\t\t\t",
    "  nu= mortality rate\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " S = total number susceptible to the infection\t\t\t",
    " I = total number of those capable of passing on the infection\t",
    " R = total number recovered\t\t\t\t\t",
    " D = total number dead\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t",
    " Default:\t\t\t\t\t\t\t",
    " sird_epidemic | xgraph n=40 nplot=3 d1=1 style=normal &\t",
    " \t\t\t\t\t\t\t\t",
    " Hypothetical Influenza in a boarding school:\t\t\t",
    " N=762 i0=1,  2 students infected per day, 1/2 of the infected\t",
    " population removed per day. Take b=2 k=0.5\t\t\t",
    " assume 10% mortality nu=0.1\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Normalized by N:\t\t\t\t\t\t",
    " sird_epidemic h=0.1 stepmax=200 i0=1 b=2 k=.5 nu=.1 N=762 mode=SIRD |",
    "  xgraph n=200 nplot=4 d1=.1 style=normal label1=\"days\"  &\t",
    " \t\t\t\t\t\t\t\t",
    " Normalized by N, output scaled by N:\t\t\t\t",
    " sird_epidemic h=0.1 stepmax=200 i0=1 b=2 k=.5 nu=.1 N=762 mode=SIRD scale=1 |",
    "  xgraph n=200 nplot=4 d1=.1 style=normal label1=\"days\" &\t",
    " \t\t\t\t\t\t\t\t",
    " Kong Flu 1968-1969:\t\t\t\t\t\t",
    " https://services.math.duke.edu/education/ccp/materials/diffcalc/sir/sir1.html",
    " Population is N=s0=7.9 million, r0=1.5, the average period of\t",
    " infectiveness is  3 days so k=1/3, b=R0*k=(3/2)(1/3)=0.5, and initial",
    " infected is i0=10.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "  Normalized by N\t\t\t\t\t\t",
    "  sird_epidemic h=1 i0=10 stepmax=200 k=.3333 b=.5 N=7.9e6 mode=SIR |\t",
    "      xgraph n=200 nplot=3 d1=1 style=normal &\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "  Normalized by N, with scaling of the output by N:\t\t",
    "  sird_epidemic h=1 scale=1 i0=10 stepmax=200 k=.3333 b=.5 N=7.9e6 mode=SIR |",
    "      xgraph n=200 nplot=3 d1=1 style=normal &\t\t\t",
    "     \t\t\t\t\t\t\t\t",
    " Related programs: sir_epidemic, seir_epidemic\t\t\t",
];

const SIRD_MODE: i32 = 0;
const S_MODE: i32 = 1;
const I_MODE: i32 = 2;
const R_MODE: i32 = 3;
const D_MODE: i32 = 4;

fn main() {
    let argv = args();
    let number = 4;
    let mut verbose = 0i32;
    let mut stepmax = 0i32;

    /* Initial values of S, I, R, D. */
    let mut normalize = 1i32;
    let mut scale = 0i32;
    let mut scalar = 0.0f32;
    let mut n = 0.0f64;
    let mut s0 = 0.0f64;
    let mut i0 = 0.0f64;
    let mut r0 = 0.0f64;
    let mut d0 = 0.0f64;
    let mut h = 0.001f64;
    let mut tol = 0.0f64;
    let mut y = [0.0f64; 4];
    let out_file = stdout();
    let mut mode = String::from("SIRD");
    let mut imode = SIRD_MODE;

    /* Hook up getpar. */
    initargs(&argv);
    requestdoc(0, SDOC);

    /* Prevent floats from dumping on screen. */
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

    /* Get parameters. */
    if getparint("stepmax", &mut stepmax) == 0 {
        stepmax = 40;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* Initial conditions y[0] = S, y[1] = I, y[2] = R, y[3] = D. */
    if getparint("normalize", &mut normalize) == 0 {
        normalize = 1;
    }
    if getpardouble("N", &mut n) == 0 {
        n = 1000.0;
    }
    if getpardouble("s0", &mut s0) == 0 {
        s0 = n;
    }
    y[0] = if normalize != 0 { s0 / n } else { s0 };
    if getpardouble("i0", &mut i0) == 0 {
        i0 = 1.0;
    }
    y[1] = if normalize != 0 { i0 / n } else { i0 };
    if getpardouble("r0", &mut r0) == 0 {
        r0 = 0.0;
    }
    y[2] = if normalize != 0 { r0 / n } else { r0 };
    if getpardouble("d0", &mut d0) == 0 {
        d0 = 0.01;
    }
    y[3] = if normalize != 0 { d0 / n } else { d0 };

    /* Scale output factor. */
    if getparint("scale", &mut scale) == 0 {
        scale = 0;
    }
    if getpardouble("h", &mut h) == 0 {
        h = 1.0;
    }
    if getpardouble("tol", &mut tol) == 0 {
        tol = 0.00000001; /* RKE_ERR_BIAS_INIT */
    }

    /* Get output mode. The second read of h is also in the C source. */
    if getpardouble("h", &mut h) == 0 {
        h = 1.0;
    }
    getparstring("mode", &mut mode);
    if mode == "S" {
        imode = S_MODE;
    } else if mode == "I" {
        imode = I_MODE;
    } else if mode == "R" {
        imode = R_MODE;
    } else if mode == "D" {
        imode = D_MODE;
    } else if mode != "SIRD" {
        err!("unknown operation=\"%s\", see self-doc", mode);
    }
    // The C program never calls checkpars(), including in its ODE callback.

    // Match the fixed C's bounds check before its signed 4*stepmax products.
    // Without this check INT_MIN wraps to zero and the output loops panic.
    if stepmax.checked_mul(4).is_none() {
        err!("stepmax*4 overflows int");
    }

    /* Allocate space in the output array, preserving C's excess rows. */
    let nrows = (stepmax * 4) as usize;
    let mut yout: Vec<[f32; 4]> = Vec::new();
    if yout.try_reserve_exact(nrows).is_err() {
        // Negative stepmax becomes a huge size_t in ealloc2float. Rust's
        // capacity check does not set errno; preserve C's ENOMEM diagnostic.
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    yout.resize(nrows, [0.0; 4]);

    /* Initialize Runge-Kutta-England routines and set tolerance. */
    let mut p = rke_init(4, sird_epidemic_equations);
    p.error_bias = tol;

    for i in 0..stepmax {
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

        /* Convert doubles in y to floats in yout before solving. */
        for j in 0..number {
            yout[i as usize][j] = y[j] as f32;
        }
        rke_solve(&mut p, &mut t, &mut y, aimed_t);
    }

    /* Write out according to the mode. */
    let mut tempout = vec![0.0f32; nrows];
    if scale == 1 {
        scalar = n as f32;
    } else if scale == 2 {
        scalar = s0 as f32;
    }

    if imode == S_MODE {
        for i in 0..stepmax as usize {
            tempout[i] = if scale != 0 {
                scalar * yout[i][0]
            } else {
                yout[i][0]
            };
        }
    } else if imode == I_MODE {
        for i in 0..stepmax as usize {
            tempout[i] = if scale != 0 {
                scalar * yout[i][1]
            } else {
                yout[i][1]
            };
        }
    } else if imode == R_MODE {
        for i in 0..stepmax as usize {
            tempout[i] = if scale != 0 {
                scalar * yout[i][2]
            } else {
                yout[i][2]
            };
        }
    } else if imode == D_MODE {
        for i in 0..stepmax as usize {
            tempout[i] = if scale != 0 {
                scalar * yout[i][3]
            } else {
                yout[i][3]
            };
        }
    } else if imode == SIRD_MODE {
        for i in 0..stepmax as usize {
            tempout[i] = if scale != 0 {
                scalar * yout[i][0]
            } else {
                yout[i][0]
            };
        }
        for i in 0..stepmax as usize {
            tempout[i + stepmax as usize] = if scale != 0 {
                scalar * yout[i][1]
            } else {
                yout[i][1]
            };
        }
        for i in 0..stepmax as usize {
            tempout[i + 2 * stepmax as usize] = if scale != 0 {
                scalar * yout[i][2]
            } else {
                yout[i][2]
            };
        }
        for i in 0..stepmax as usize {
            tempout[i + 3 * stepmax as usize] = if scale != 0 {
                scalar * yout[i][3]
            } else {
                yout[i][3]
            };
        }
    }

    if imode == SIRD_MODE {
        efwrite(&tempout[..nrows], &out_file);
    } else {
        efwrite(&tempout[..stepmax as usize], &out_file);
    }
    rke_term(p);
    exit(CWP_Exit());
}

/// The autonomous system of ODEs describing the SIRD epidemic model.
fn sird_epidemic_equations(_t: f64, y: &mut [f64], yprime: &mut [f64]) -> i32 {
    let mut b = 0.0f64; /* Infection rate. */
    let mut k = 0.0f64; /* Recovery rate. */
    let mut nu = 0.0f64; /* Mortality rate. */
    let mut gamma = 0.0f64; /* Social distancing parameter. */

    /* Parameters are read anew on each evaluation, as in C. */
    if getpardouble("b", &mut b) == 0 {
        b = 0.5;
    }
    if getpardouble("k", &mut k) == 0 {
        k = 0.333;
    }
    if getpardouble("nu", &mut nu) == 0 {
        nu = 0.0;
    }
    if getpardouble("gamma", &mut gamma) == 0 {
        gamma = 0.0;
    }

    yprime[0] = -b * y[0] * y[1] / (1.0 + gamma * y[1]);
    yprime[1] = b * y[0] * y[1] / (1.0 + gamma * y[1]) - k * y[1] - nu * y[1];
    yprime[2] = k * y[1];
    yprime[3] = nu * y[1];
    1
}
