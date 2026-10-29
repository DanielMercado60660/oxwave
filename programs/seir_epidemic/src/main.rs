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

//! SEIREPIDEMIC - the SEIR and SEIRS epidemic models.
//! Port of par/main/ordinary_differential_equations/seir_epidemic.c (Revision 1.4).
//! Credits (from the C source): April 2020: John Stockwell.

use oxwave_core::prelude::*;
use oxwave_f_rke::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SEIREPIDEMIC - the SEIR and SEIRS EPIDEMIC models\t\t",
    "\t\t\t\t\t\t\t\t",
    "  seirepidemic > [stdout]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required Parameters: none\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " normalize=1\t\tnormalize S,E,I by N; =0 don't normalize",
    " scale=0\t\tdon't scale\t\t\t\t",
    "\t\t\t=1 scale output S,E,I by N\t\t",
    "\t\t\t=2 scale output S,E,I by s0\t\t",
    " N=1000\t\ttotal population\t\t\t",
    " s0=N\t\t\tinitial number of susceptibles\t\t",
    " e0=1\t\t\tinitial number of exposed\t\t",
    " i0=1\t\t\tinitial number of infectives\t\t",
    " r0=0.0\t\tinitial number of removed (should be 0)\t",
    " \t\t\t\t\t\t\t\t",
    " k=.5\t\t\ttransmission rate\t\t\t",
    " eti=.3333\t\texposure to infective rate\t\t",
    " b=.3333\t\tremoval rate = death + recovery rates\t",
    " \t\t\t\t\t\t\t\t",
    " .... with vital dynamics (i.e. birth and death rate) \t\t",
    " mu=0.0\t\tbirth rate \t\t\t\t",
    " nu=0.0\t\tdeath rate \t\t\t\t",
    "  ... SEIRS ... with reinfection\t\t\t\t",
    " xi=0.0\t\tre-infection parameter\t\t  ",
    " \t\t\t\t\t\t\t\t",
    " h=1\t\t\tincrement in time\t\t\t",
    " tol=1.e-08\t\terror tolerance\t\t\t\t",
    " stepmax=40\t\tmaximum number of steps to compute\t",
    " mode=SEIR\t\tS -> E -> I -> R\t\t\t",
    "\t\t\t=S only, =E only, =I only, =R only\t",
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
    " r0 = number of new infections per single infected host\t",
    "  1 < r0 < 1.5 for influenza, (2.2 to 2.7 for Covid-19), 12 to ",
    " 18 for measles.\t\t\t\t\t\t",
    "  b, k, s0, and r0 are related via\t\t\t\t",
    "  k = b*s0/r0 = b/r0 when s0/N and s0=N\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " S = number susceptible to the infection\t\t\t",
    " E = number of exposed\t\t\t\t\t\t",
    " I = number of those capable of passing on the infection\t",
    " R = number removed = dead + recovered\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t",
    " Default:\t\t\t\t\t\t\t",
    " seirpidemic | xgraph n=40 nplot=4 d1=1 style=normal &\t\t",
    "\t\t\t\t\t\t\t\t",
    " Hypothetical flu in a school models:\t\t\t\t",
    " here there is a 0.3 exposure to infective rate\t\t",
    " seirepidemic h=0.1 stepmax=100 s=.3 b=.3 k=.9 N=100 mode=SIR |",
    "\txgraph n=100 nplot=4 d1=1 style=normal &\t\t",
    "\t\t\t\t\t\t\t\t",
    " Hong Kong Flu 1968-1969:\t\t\t\t\t",
    " https://services.math.duke.edu/education/ccp/materials/diffcalc/sir/sir1.html",
    " Population is N=s0=7.9 million, r0=1.5, the average period of \t",
    " infectiveness is  3 days so k=1/3, b=r0*k=(3/2)(1/3)=0.54, and initial",
    " infected is i0=10. An exposed to infective rate s=1/3 is assumed\t",
    "\t\t\t\t\t\t\t\t\t",
    "  seirepidemic h=1 stepmax=200 s=.3333 k=.3333 b=.5 N=7.9e6 mode=SIR |\t",
    "\txgraph n=200 nplot=4 d1=1 style=normal &\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Related programs: sir_epidemic, sird_epidemic \t\t\t",
];

const SEIR_MODE: i32 = 0;
const S_MODE: i32 = 1;
const E_MODE: i32 = 2;
const I_MODE: i32 = 3;
const R_MODE: i32 = 4;

#[allow(non_snake_case, unused_assignments)]
fn main() {
    let argv = args();
    let number = 4;
    let mut verbose = 0;
    let mut stepmax: i32 = 0;
    let mut normalize = 1;
    let mut scale = 0;
    let mut scalar = 0.0f32;
    let mut N = 0.0f64;
    let mut s0 = 0.0f64;
    let mut e0 = 0.0f64;
    let mut i0 = 0.0f64;
    let mut r0 = 0.0f64;
    let mut t = 0.0f64;
    let mut h = 0.001f64;
    let mut tol = 0.0f64;
    let mut y = [0.0f64; 4];
    let out_file = stdout();
    let mut mode = String::from("SEIR");
    let mut imode = SEIR_MODE;

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
        stepmax = 40;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    // Initial conditions: preserve the C assignments, including the unused e0
    // and y[2] remaining zero. i0 is assigned to the exposed state y[1].
    if getparint("normalize", &mut normalize) == 0 {
        normalize = 1;
    }
    if getpardouble("N", &mut N) == 0 {
        N = 1000.0;
    }
    if getpardouble("s0", &mut s0) == 0 {
        s0 = N;
    }
    y[0] = if normalize != 0 { s0 / N } else { s0 };
    if getpardouble("e0", &mut e0) == 0 {
        e0 = 1.0;
    }
    y[1] = if normalize != 0 { i0 / N } else { i0 };
    if getpardouble("i0", &mut i0) == 0 {
        i0 = 1.0;
    }
    y[1] = if normalize != 0 { i0 / N } else { i0 };
    if getpardouble("r0", &mut r0) == 0 {
        r0 = 0.0;
    }
    y[3] = if normalize != 0 { r0 / N } else { r0 };

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
    } else if mode == "E" {
        imode = E_MODE;
    } else if mode == "I" {
        imode = I_MODE;
    } else if mode == "R" {
        imode = R_MODE;
    } else if mode != "SEIR" {
        err!("unknown operation=\"%s\", see self-doc", mode);
    }

    // ealloc2float(4, 4*stepmax): a negative dimension becomes an
    // unallocatable size_t in C. Preserve its allocation failure message.
    let rows = stepmax.wrapping_mul(4) as usize;
    let mut yout = Vec::<[f32; 4]>::new();
    if yout.try_reserve_exact(rows).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    yout.resize(rows, [0.0; 4]);

    let mut p = rke_init(4, sir_epidemic_equations);
    p.error_bias = tol;

    for i in 0..stepmax {
        t = i as f64 * h;
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
            yout[i as usize][j] = y[j] as f32;
        }
        rke_solve(&mut p, &mut t, &mut y, aimed_t);
    }

    let n = stepmax as usize;
    let mut tempout = vec![0.0f32; stepmax.wrapping_mul(4) as usize];
    if scale == 1 {
        scalar = N as f32;
    } else if scale == 2 {
        scalar = s0 as f32;
    }

    if imode == S_MODE {
        for i in 0..n {
            tempout[i] = if scale != 0 {
                scalar * yout[i][0]
            } else {
                yout[i][0]
            };
        }
    } else if imode == E_MODE {
        for i in 0..n {
            tempout[i] = if scale != 0 {
                scalar * yout[i][1]
            } else {
                yout[i][1]
            };
        }
    } else if imode == I_MODE {
        for i in 0..n {
            tempout[i] = if scale != 0 {
                scalar * yout[i][2]
            } else {
                yout[i][2]
            };
        }
    } else if imode == R_MODE {
        for i in 0..n {
            tempout[i] = if scale != 0 {
                scalar * yout[i][3]
            } else {
                yout[i][3]
            };
        }
    } else if imode == SEIR_MODE {
        for i in 0..n {
            tempout[i] = if scale != 0 {
                scalar * yout[i][0]
            } else {
                yout[i][0]
            };
        }
        for i in 0..n {
            tempout[i + n] = if scale != 0 {
                scalar * yout[i][1]
            } else {
                yout[i][1]
            };
        }
        for i in 0..n {
            tempout[i + 2 * n] = if scale != 0 {
                scalar * yout[i][2]
            } else {
                yout[i][2]
            };
        }
        for i in 0..n {
            tempout[i + 3 * n] = if scale != 0 {
                scalar * yout[i][3]
            } else {
                yout[i][3]
            };
        }
    }

    if imode == SEIR_MODE {
        efwrite(&tempout[..4 * n], &out_file);
    } else {
        efwrite(&tempout[..n], &out_file);
    }
    rke_term(p);
    exit(CWP_Exit());
}

// Autonomous system of ODEs describing the SEIR epidemic model.
fn sir_epidemic_equations(_t: f64, y: &mut [f64], yprime: &mut [f64]) -> i32 {
    let mut b = 0.0f64;
    let mut eti = 0.0f64;
    let mut k = 0.0f64;
    let mut mu = 0.0f64;
    let mut nu = 0.0f64;
    let mut xi = 0.0f64;

    if getpardouble("b", &mut b) == 0 {
        b = 0.5;
    }
    // C reads eti into k, leaving eti at zero when the parameter is supplied.
    if getpardouble("eti", &mut k) == 0 {
        eti = 0.333;
    }
    if getpardouble("k", &mut k) == 0 {
        k = 0.333;
    }
    if getpardouble("mu", &mut mu) == 0 {
        mu = 0.0;
    }
    if getpardouble("nu", &mut nu) == 0 {
        nu = 0.0;
    }
    if getpardouble("xi", &mut xi) == 0 {
        xi = 0.0;
    }

    yprime[0] = mu - nu * y[0] - b * y[0] * y[1] + xi * y[3];
    yprime[1] = b * y[0] * y[1] - nu * y[1] - eti * y[1];
    yprime[2] = eti * y[1] - k * y[2] - nu * y[2];
    yprime[3] = k * y[2] - nu * y[3] - xi * y[3];
    1
}
