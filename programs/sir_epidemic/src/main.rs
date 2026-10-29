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

//! SIR_EPIDEMIC - the SIR and SIRS epidemic models with and without vital dynamics.
//! Port of par/main/ordinary_differential_equations/sir_epidemic.c (Revision 1.4).
//! Credits (from the C source): Author: April 2020: John Stockwell.

use oxwave_core::prelude::*;
use oxwave_f_rke::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SIR_EPIDEMIC - the SIR and SIRS EPIDEMIC models with and without",
    "\t\t vital dynamics\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "  sir_epidemic > [stdout]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required Parameters: none\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " normalize=1\t\tNormalize s, i by N; =0 don't normalize\t",
    " N=1000\t\ttotal population size\t\t\t",
    " s0=N\t\t\tinitial number of susceptibles\t\t",
    " i0=1\t\t\tinitial number of infectives\t\t",
    " r0=0.0\t\tinitial number of removed (should be 0)\t",
    "\t \t\t(not the basic reproducion rate R0)\t",
    " scale=0\t\tdon't scale output\t\t\t",
    "\t\t\t=1 scale output s,i,r by N\t\t",
    "\t\t\t=2 scale output s,i,r by s0\t\t",
    " k=.5\t\t\ttransmission rate\t\t\t",
    " b=.3333\t\tremoval rate = death + recovery rates\t",
    "  ... with vital dynamics\t\t\t\t\t",
    " mu=0.0\t\tbirth rate\t\t\t\t",
    " nu=0.0\t\tdeath rate\t\t\t\t",
    "  ... SIRS ... with reinfection\t\t\t\t",
    " xi=0.0\t\tre-infection parameter\t\t\t",
    " ... with social distancing (Baker 2020)\t\t\t",
    " gamma=0\t\tno social distancing; gamma > 0 distancing",
    " \t\t\t(Baker 2020 uses value of 20)\t\t",
    " h=1\t\t\tincrement in time\t\t\t",
    " tol=1.e-08\t\terror tolerance\t\t\t\t",
    " stepmax=40\t\tmaximum number of steps to compute\t",
    " mode=SIR\t\tS followed by I, followed by R\t\t",
    "\t\t\t=S only, =I only, =R only\t\t",
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
    " About compartmentalized models: The population is assumed to  ",
    " move from being Susceptible, to Infective, and finally to the ",
    " Removed, who are dead and the recovered.\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Important quantities:\t\t\t\t\t\t",
    " R0 = number of new infections per single infected host  \t",
    "  1 < R0 < 1.5 for influenza, (2.2 to 2.7 for Covid-19), 12 to\t",
    " 18 for measles.\t\t\t\t\t\t",
    "  b, k, s0, and R0 are related via\t\t\t\t",
    "  k = b*s0/R0 = b/R0 when  s0=s0/N and s0=N \t\t\t",
    "  \t\t\t\t\t\t\t\t",
    "  It is often easier to determine the recovery rate k (in units",
    "  of h and to determine reasonable estimate of s0 and of R0 \t",
    "  and to calculate the infection rate b = k*R0/s0 or b=k*R0\t",
    "  when s0=N and is normalized by N.\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " S = total number susceptible to the infection\t\t\t",
    " I = total number of those capable of passing on the infection\t",
    " R = total number removed = dead + recovered\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " When xi is nonzero, then there is a potential that fraction of",
    " the removed population can be reinfected.\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t",
    " Default:\t\t\t\t\t\t\t",
    " sir_epidemic | xgraph n=40 nplot=3 d1=1 style=normal &\t",
    " \t\t\t\t\t\t\t\t",
    " Influenza in an English boarding school, 1978:\t\t",
    " s0=N=762 i0=1,  2 students infected per day, 1/2 of the infected\t",
    " population removed per day. Take b=2 k=0.5 \t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Normalized by N:\t\t\t\t\t\t",
    " sir_epidemic h=0.1 stepmax=200 i0=1 b=2 k=.5 N=762 mode=SIR |\t",
    "  xgraph n=200 nplot=3 d1=.1 style=normal label1=\"days\"  &\t",
    " \t\t\t\t\t\t\t\t",
    " Normalized by N, output scaled by N:\t\t\t\t",
    " sir_epidemic h=0.1 stepmax=200 i0=1 b=2 k=.5 N=762 mode=SIR scale=1 |",
    "  xgraph n=200 nplot=3 d1=.1 style=normal label1=\"days\" &\t",
    " \t\t\t\t\t\t\t\t",
    " Kong Flu 1968-1969:\t\t\t\t\t\t",
    " https://services.math.duke.edu/education/ccp/materials/diffcalc/sir/sir1.html",
    " Population is N=s0=7.9 million, R0=1.5, the average period of\t",
    " infectiveness is  3 days so k=1/3, b=R0*k=(3/2)(1/3)=0.5, and initial",
    " infected is i0=10.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "  Normalized by N\t\t\t\t\t\t",
    "  sir_epidemic h=1 stepmax=200 k=.3333 b=.5 N=7.9e6 mode=SIR |\t",
    "      xgraph n=200 nplot=3 d1=1 style=normal &\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "  Normalized by N, with scaling of the output by N:\t\t",
    "  sir_epidemic h=1 scale=1 stepmax=200 k=.3333 b=.5 N=7.9e6 mode=SIR |",
    "      xgraph n=200 nplot=3 d1=1 style=normal &\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Related programs: sird_epidemic, seir_epidemic \t\t",
];

const SIR_MODE: i32 = 0;
const S_MODE: i32 = 1;
const I_MODE: i32 = 2;
const R_MODE: i32 = 3;

fn main() {
    let argv = args();
    let number = 3;
    let mut verbose = 0;
    let mut stepmax: i32 = 0;
    let mut normalize = 1;
    let mut scale = 0;
    let mut scalar = 0.0f32;
    let mut n = 0.0f64;
    let mut s0 = 0.0f64;
    let mut i0 = 0.0f64;
    let mut r0 = 0.0f64;
    let mut h = 0.001f64;
    let mut tol = 0.0f64;
    let mut y = [0.0f64; 3];
    let out_file = stdout();
    let mut mode = String::from("SIR");
    let mut imode = SIR_MODE;

    initargs(&argv);
    requestdoc(0, SDOC);
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
        stepmax = 40;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    // Initial conditions y[0] = s, y[1] = I, y[2] = R.
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

    if getparint("scale", &mut scale) == 0 {
        scale = 0;
    }
    if getpardouble("h", &mut h) == 0 {
        h = 1.0;
    }
    if getpardouble("tol", &mut tol) == 0 {
        tol = 0.00000001;
    }

    if getpardouble("h", &mut h) == 0 {
        h = 1.0;
    }
    getparstring("mode", &mut mode);
    if mode == "S" {
        imode = S_MODE;
    } else if mode == "S" {
        // The repeated S test is present in the C source.
        imode = I_MODE;
    } else if mode == "R" {
        imode = R_MODE;
    } else if mode != "SIR" {
        err!("unknown operation=\"%s\", see self-doc", mode);
    }

    // ealloc2float(3, 3*stepmax), including its allocation failure for
    // negative dimensions. Vec capacity errors do not set libc errno.
    let rows = stepmax.wrapping_mul(3) as usize;
    let len = rows.checked_mul(3).unwrap_or(usize::MAX);
    let mut yout = Vec::<f32>::new();
    if yout.try_reserve_exact(len).is_err() {
        err!("%s: malloc failed (Cannot allocate memory)", "ealloc.c");
    }
    yout.resize(len, 0.0);

    let mut p = rke_init(3, sir_epidemic_equations);
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
        for j in 0..number {
            yout[i as usize * 3 + j] = y[j] as f32;
        }
        rke_solve(&mut p, &mut t, &mut y, aimed_t);
    }

    let mut tempout = vec![0.0f32; rows];
    if scale == 1 {
        scalar = n as f32;
    } else if scale == 2 {
        scalar = s0 as f32;
    }
    let steps = stepmax as usize;
    if imode == S_MODE {
        for i in 0..steps {
            tempout[i] = if scale != 0 {
                scalar * yout[i * 3]
            } else {
                yout[i * 3]
            };
        }
    } else if imode == I_MODE {
        for i in 0..steps {
            tempout[i] = if scale != 0 {
                scalar * yout[i * 3 + 1]
            } else {
                yout[i * 3 + 1]
            };
        }
    } else if imode == R_MODE {
        for i in 0..steps {
            tempout[i] = if scale != 0 {
                scalar * yout[i * 3 + 2]
            } else {
                yout[i * 3 + 2]
            };
        }
    } else if imode == SIR_MODE {
        for i in 0..steps {
            tempout[i] = if scale != 0 {
                scalar * yout[i * 3]
            } else {
                yout[i * 3]
            };
        }
        for i in 0..steps {
            tempout[i + steps] = if scale != 0 {
                scalar * yout[i * 3 + 1]
            } else {
                yout[i * 3 + 1]
            };
        }
        for i in 0..steps {
            tempout[i + 2 * steps] = if scale != 0 {
                scalar * yout[i * 3 + 2]
            } else {
                yout[i * 3 + 2]
            };
        }
    }
    if imode == SIR_MODE {
        efwrite(&tempout[..rows], &out_file);
    } else {
        efwrite(&tempout[..steps], &out_file);
    }
    rke_term(p);
    exit(CWP_Exit());
}

/// The autonomous system of ODEs describing the SIR epidemic model.
fn sir_epidemic_equations(_t: f64, y: &mut [f64], yprime: &mut [f64]) -> i32 {
    let mut b = 0.0f64;
    let mut k = 0.0f64;
    let mut xi = 0.0f64;
    let mut mu = 0.0f64;
    let nu = 0.0f64; // The C program never reads the documented nu parameter.
    let mut gamma = 0.0f64;
    if getpardouble("b", &mut b) == 0 {
        b = 0.5;
    }
    if getpardouble("k", &mut k) == 0 {
        k = 0.333;
    }
    if getpardouble("xi", &mut xi) == 0 {
        xi = 0.0;
    }
    if getpardouble("mu", &mut mu) == 0 {
        mu = 0.0;
    }
    if getpardouble("gamma", &mut gamma) == 0 {
        gamma = 0.0;
    }

    yprime[0] = mu - b * y[0] * y[1] / (1.0 + gamma * y[1]) + xi * y[2] - nu * y[0];
    yprime[1] = b * y[0] * y[1] / (1.0 + gamma * y[1]) - k * y[1] - nu * y[1];
    yprime[2] = k * y[1] - xi * y[2] - nu * y[2];
    1
}
