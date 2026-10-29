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

//! VOLTLOTKA - modified Volterra-Lotka predator-prey model.
//! Port of par/main/ordinary_differential_equations/voltlotka.c (Revision 1.3).
//! Credits (from the C source): May 2020: John Stockwell.

use oxwave_core::prelude::*;
use oxwave_f_rke::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " VOLTLOTKA -  VOLTerra-LOTKA predator-prey model, modified \t",
    "\t        slightly to include a carrying capacity of the  ",
    "\t\tprey species\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "  voltlotka  > [stdout]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required Parameters: none\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " a=0.5\t\t\tprey increase parameter\t\t\t",
    " b=0.02\t\tprey reduction by predation parameter\t",
    " c=0.75\t\tpredator reduction rate\t\t\t",
    " d=0.005\t\tpredator increase by prey capture\t",
    "\t\t\t\t\t\t\t\t",
    " K=1000000\t\tcarrying capacity of prey environment \t",
    "\t\t\treduces to classical Volterra-Lotka when",
    "\t\t\tK >> x and K >> y\t",
    "\t\t\t\t\t\t\t\t",
    " x0=500\t\tinitial number of prey\t\t\t",
    " y0=25\t\t\tinitial number of predator\t\t",
    "\t\t\t\t\t\t\t\t",
    " thresh=1.0\t\textinction if x <= thresh or y <= thresh",
    "\t\t\t\t\t\t\t\t",
    " h=1.0\t\t\tstep size\t\t\t\t",
    " stepmax=100\t\tmaximum number of steps to compute\t",
    " mode=PP\t\tnumber of Pred. followed by number of Prey",
    "\t\t\t=R predator only, =Y prey only\t\t",
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
    " Modified Volterra-Lotka equations: \t\t\t\t",
    " dx/dt = a*x*(1 - x/K) - b*x*y\t\t\t\t\t",
    " dy/dt = -c*y + p*x*y\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " x = number of members of the prey species\t\t\t",
    " y = number of members of the predator species \t\t",
    " bxy = prey population reduction due to predation\t\t",
    " pxy = predator population increase due due to predation\t",
    " a = prey increase parameter (prey have unlimited food)\t",
    " c = predator reduction (predators starve without prey)\t",
    " K = carrying capacity of for the prey specie, when K large\t",
    "     (1 - x/K)-> 1 and the system reduces to original VL system",
    "\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t",
    " Default: classic Snowshoe hare versus Canadian lynx\t\t",
    " voltlotka | xgraph n=100 nplot=2 d1=1 style=normal &\t\t",
    " \t\t\t\t\t\t\t\t",
    " Caveat: if there is weird behavior, try reducing the h= value\t",
];

const PP_MODE: i32 = 0;
const R_MODE: i32 = 1;
const Y_MODE: i32 = 2;

fn main() {
    let argv = args();
    let number = 2;
    let mut verbose = 0;
    let mut stepmax: i32 = 0;
    let mut x0 = 0.0;
    let mut y0 = 0.0;
    let mut h = 0.001;
    let mut tol = 0.0;
    // C supplies two doubles to a three-equation solver. Reserve the extra
    // slot explicitly. Zero padding matches the checked reference outputs;
    // the original third state and derivative are undefined.
    let mut y = [0.0; 3];
    let out_file = stdout();
    let mut mode = String::from("PP");
    let mut imode = PP_MODE;

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
        stepmax = 100;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    if getpardouble("x0", &mut x0) == 0 {
        x0 = 500.0;
    }
    y[0] = x0;
    if getpardouble("y0", &mut y0) == 0 {
        y0 = 25.0;
    }
    y[1] = y0;
    if getpardouble("h", &mut h) == 0 {
        h = 1.0;
    }
    if getpardouble("tol", &mut tol) == 0 {
        tol = 0.00000001; // RKE_ERR_BIAS_INIT
    }
    if getpardouble("h", &mut h) == 0 {
        h = 1.0;
    }
    getparstring("mode", &mut mode);
    if mode == "R" {
        imode = R_MODE;
    } else if mode == "Y" {
        imode = Y_MODE;
    } else if mode != "PP" {
        err!("unknown operation=\"%s\", see self-doc", mode);
    }

    // Negative dimensions are converted to size_t by ealloc2float and fail
    // allocation in the reference, rather than entering the integration loop.
    let rows = stepmax.wrapping_mul(2) as usize;
    let mut yout = Vec::<[f32; 2]>::new();
    if yout.try_reserve_exact(rows).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    yout.resize(rows, [0.0; 2]);

    // The callback zeroes extinct populations in the live state and solver
    // stages; subsequent calculations must use those mutated values.
    let mut p = rke_init(3, modified_volterra_lotka_equations);
    p.error_bias = tol;

    for i in 0..stepmax as usize {
        let mut t = i as f64 * h;
        let aimed_t = t + h;
        if verbose != 0 {
            warn!(
                "using %3d accepted and %3d rejected steps",
                p.accepted_steps, p.rejected_steps
            );
            warn!("error tolerance = %10.24f", p.error_bias);
        }
        for j in 0..number {
            yout[i][j] = y[j] as f32;
        }

        rke_solve(&mut p, &mut t, &mut y, aimed_t);
    }

    let mut tempout = vec![0.0f32; rows];
    if imode == R_MODE {
        for i in 0..stepmax as usize {
            tempout[i] = yout[i][0];
        }
    } else if imode == Y_MODE {
        for i in 0..stepmax as usize {
            tempout[i] = yout[i][1];
        }
    } else if imode == PP_MODE {
        for i in 0..stepmax as usize {
            tempout[i] = yout[i][0];
        }
        for i in 0..stepmax as usize {
            tempout[i + stepmax as usize] = yout[i][1];
        }
    }
    if imode == PP_MODE {
        efwrite(&tempout, &out_file);
    } else {
        efwrite(&tempout[..stepmax as usize], &out_file);
    }
    rke_term(p);
    exit(CWP_Exit());
}

/// Autonomous modified Volterra-Lotka system from the C program.
#[allow(non_snake_case)]
fn modified_volterra_lotka_equations(_t: f64, y: &mut [f64], yprime: &mut [f64]) -> i32 {
    let mut a = 0.0;
    let mut b = 0.0;
    let mut c = 0.0;
    let mut p = 0.0;
    let mut thresh = 1.0;
    let mut K = 0.0;

    if getpardouble("a", &mut a) == 0 {
        a = 0.5;
    }
    if getpardouble("b", &mut b) == 0 {
        b = 0.02;
    }
    if getpardouble("c", &mut c) == 0 {
        c = 0.75;
    }
    if getpardouble("p", &mut p) == 0 {
        p = 0.005;
    }
    if getpardouble("K", &mut K) == 0 {
        K = 1000000.0;
    }
    if getpardouble("thresh", &mut thresh) == 0 {
        thresh = 1.0;
    }

    if y[0] <= thresh {
        y[0] = 0.0;
    }
    if y[1] <= thresh {
        y[1] = 0.0;
    }
    yprime[0] = a * y[0] * (1.0 - y[0] / K) - b * y[0] * y[1];
    yprime[1] = -c * y[1] + p * y[0] * y[1];
    1
}
