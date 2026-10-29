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

//! LOGISTICFIT - least squares fit of a logistic function to input data.
//! Port of par/main/ordinary_differential_equations/logisticfit.c (Revision 1.3).
//! Credits: John Stockwell, Adjunct Faculty, Colorado School of Mines, April 2020.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_linear_regression::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " LOGISTICFIT - least squares fit of a LOGISTIC function to input data \t",
    " \t\t\t\t\t\t\t\t\t",
    "  logisticfit <infile >outfile [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameter:\t\t\t\t\t\t\t",
    " n=\t\t\tthe number of values of t and P(t)\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nstart=0\t\tstart of input data window\t\t\t",
    " nend=n\t\tend of input data window\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " outpar=/dev/tty \toutput parameter file\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " If the data are a population P(t) as a function of time t, then the \t",
    " input consists of unformatted C-style floats in the form of a vector of",
    " length 2n consisting of the vector of the values of \t\t\t",
    " [t0,t1,t2,...tn-1)] followed by [P(t0),P(t1),P(t2),...,P(tn-1)].\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Notes about logistic curve fitting:\t\t\t\t\t",
    " The logistic function is a nontrivial solution to the Verhulst equation",
    " \t\t\t\t\t\t\t\t\t",
    "\t\t dP/dt = r P(1-P/K)\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " which describes the growth of a population P(t), with growth rate r   ",
    " and carrying capacity K. This equation has analytic solution:\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "\t       P(t)= K/(1 + (K-P0)/P0 exp(-r(t-t0))\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Here P0=P(t0). \t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The goal of this program is, that given population data that is assumed",
    " to be governed by the logistic equation, to estimate K and r. These   ",
    " may then be used as input to the program 'verhulst'. In that program  ",
    " a1 corresponds to r and a2 corresponds to K here. \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat: this code will give results that get you close, but you may   ",
    " have to hand adjust r and K to suit your data. \t\t\t",
    " Related program: verhulst \t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    /* Hook up getpar (the original does this twice). */
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut outpar = String::new();
    let mut n = 0i32;
    let mut nstart = 0i32;
    let mut nend = 0i32;
    let mut ncount = 0i32;
    let mut tp = [0.0f32];

    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get parameters and do set up. */
    if getparstring("outpar", &mut outpar) == 0 {
        outpar = "/dev/tty".to_string();
    }
    let outparfp = efopen(&outpar, "w");

    mustgetparint("n", &mut n);
    if getparint("nstart", &mut nstart) == 0 {
        nstart = 0;
    }
    if getparint("nend", &mut nend) == 0 {
        nend = n;
    }
    if nend < nstart {
        err!("nend cannot be greater than nstart!");
    }
    if nend > n {
        err!("nend cannot be greater than n!");
    }
    if nstart < 0 {
        err!("nstart cannot be less than 0!");
    }
    if nend < 0 {
        err!("nstart cannot be less than 0!");
    }
    let nwindow = nend.wrapping_sub(nstart) as usize;

    // C reads p[1], t[1], p[-1], and t[-1] for a one-point window.
    // Explicit zero padding reproduces the fresh allocation contents observed
    // in the reference without an out-of-bounds Rust access. Index 1 is C's 0.
    let mut t = vec![0.0f32; nwindow.max(2) + 1];
    let mut p = vec![0.0f32; nwindow.max(2) + 1];

    /* Loop over data reading time and p data. */
    while efread(&mut tp, &stdin()) != 0 {
        if ncount >= nstart && ncount < nend {
            t[ncount.wrapping_sub(nstart) as usize + 1] = tp[0];
        }
        if ncount >= n.wrapping_add(nstart) {
            let i = ncount.wrapping_sub(n.wrapping_add(nstart)) as usize + 1;
            // The C population write has no upper window bound. Keep surplus
            // values in backing storage; the regressions still use nwindow.
            if i >= p.len() {
                p.resize(i + 1, 0.0);
            }
            p[i] = tp[0];
        }
        ncount = ncount.wrapping_add(1);
    }

    /* Part I - linear fit of (P,d(ln(P))/dt). */
    let (k, kfit, kerr) = {
        let mut coeff = [0.0f32; 4];
        let mut logpprime = vec![0.0f32; nwindow.max(1)];

        if nwindow != 0 {
            /* Differentiate p[i] and divide by p[i]. */
            /* First as a leading difference. */
            logpprime[0] = ((p[2] - p[1]) / (t[2] - t[1])) / p[1];
            /* Middle values as a centered difference. */
            for i in 1..nwindow - 1 {
                logpprime[i] = ((p[i + 2] - p[i]) / (2.0f32 * (t[i + 2] - t[i]))) / p[i + 1];
            }
            /* Last as a lagging difference. */
            logpprime[nwindow - 1] =
                ((p[nwindow] - p[nwindow - 1]) / (t[nwindow] - t[nwindow - 1])) / p[nwindow];

            linear_regression(&logpprime[..nwindow], &p[1..nwindow + 1], &mut coeff);
        } else {
            // A zero window has undefined accesses in C, and the regression
            // crate requires a nonempty slice. Preserve the reference's
            // observed degenerate result, rather than calling it illegally.
            coeff = [f32::NAN, f32::NAN, f32::NAN, -0.0];
        }

        let k = nint((-coeff[1] / coeff[0]) as f64) as f32;
        if k < p[nwindow] {
            warn!(
                "K = %f is less than largest value p[%d]=%f!",
                k,
                nwindow as i32 - 1,
                p[nwindow]
            );
        }
        (k, coeff[2], coeff[3])
    };

    /* Part II - linear fit of (log(p[i]/(K-p[i])),t). */
    let (r, rfit, rerr) = {
        let mut coeff = [0.0f32; 4];
        let mut logpp = vec![0.0f32; nwindow];
        for i in 0..nwindow {
            logpp[i] = ((p[i + 1] / (k - p[i + 1])) as f64).ln() as f32;
        }
        if nwindow != 0 {
            linear_regression(&logpp, &t[1..nwindow + 1], &mut coeff);
        } else {
            // Same zero-window C boundary behavior as in Part I.
            coeff = [f32::NAN, f32::NAN, 999.0, -0.0];
        }
        (coeff[0], coeff[2], coeff[3])
    };

    /* Make par file. */
    fprintf!(outparfp, "\n");
    fprintf!(outparfp, "Carrying capacity estimate: K = %0.0f\n", k);
    fprintf!(outparfp, "Growth rate estimate: r = %f\n", r);
    fprintf!(outparfp, "\n");
    fprintf!(
        outparfp,
        "Goodness of fit: K-fit = %f, percent error: K-err = %f\n",
        kfit,
        kerr
    );
    fprintf!(
        outparfp,
        "Goodness of fit: r-fit = %f, percent error: r-err = %f\n",
        rfit,
        rerr
    );

    exit(CWP_Exit());
}
