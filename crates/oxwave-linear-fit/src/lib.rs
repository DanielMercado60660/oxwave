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

//! Port of `cwp/lib/linfit.c` (Colorado School of Mines, 2011).
//! Author: Balasz Nemeth (inspired by similar Numerical Recipes code,
//! including `fit.c` for `linfit`).

use oxwave_core::prelude::{exit, fprintf, stderr, EXIT_FAILURE};

/// Natural logarithm of the gamma function, using the original coefficients.
///
/// C prototype: `float gammln(float xx);`
pub fn gammln(xx: f32) -> f32 {
    let cof = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ];
    let x = xx as f64;
    let mut y = x;
    let mut tmp = x + 5.5;
    tmp -= (x + 0.5) * tmp.ln();
    let mut ser = 1.000000000190015;
    for j in 0..=5 {
        y += 1.0;
        ser += cof[j] / y;
    }
    (-tmp + (2.5066282746310005 * ser / x).ln()) as f32
}

/// Upper incomplete gamma probability, selecting the C series/fraction branch.
///
/// C prototype: `float gammq(float a, float x);`
///
/// Like C, invalid arguments (`x < 0` or `a <= 0`) terminate the process
/// with status 1 after printing the original diagnostic.
pub fn gammq(a: f32, x: f32) -> f32 {
    if x < 0.0 || a <= 0.0 {
        fprintf!(stderr(), "Invalid arguments in routine gammq\n");
        exit(EXIT_FAILURE);
    }
    if (x as f64) < a as f64 + 1.0 {
        // C leaves gamser uninitialized if gser fails to converge. Initialize
        // to NaN here so that this otherwise undefined C case is deterministic.
        let mut gamser = f32::NAN;
        let mut gln = 0.0;
        gser(&mut gamser, a, x, &mut gln);
        (1.0 - gamser as f64) as f32
    } else {
        let mut gammcf = 0.0;
        let mut gln = 0.0;
        gcf(&mut gammcf, a, x, &mut gln);
        gammcf
    }
}

/// Upper incomplete gamma probability from a continued fraction.
///
/// C prototype: `void gcf(float *gammcf, float a, float x, float *gln);`
///
/// Writes `gammln(a)` to `gln`. Failure to converge within 100 iterations
/// prints the C diagnostic and terminates the process with status 1.
pub fn gcf(gammcf: &mut f32, a: f32, x: f32, gln: &mut f32) {
    const ITMAX: i32 = 100;
    const EPS: f64 = 3.0e-7;
    const FPMIN: f64 = 1.0e-30;

    *gln = gammln(a);
    let mut b = (x as f64 + 1.0 - a as f64) as f32;
    let mut c = (1.0 / FPMIN) as f32;
    let mut d = (1.0 / b as f64) as f32;
    let mut h = d;
    let mut i = 1;
    while i <= ITMAX {
        let an = (-i) as f32 * (i as f32 - a);
        b = (b as f64 + 2.0) as f32;
        d = an * d + b;
        if (d.abs() as f64) < FPMIN {
            d = FPMIN as f32;
        }
        c = b + an / c;
        if (c.abs() as f64) < FPMIN {
            c = FPMIN as f32;
        }
        d = (1.0 / d as f64) as f32;
        let del = d * c;
        h *= del;
        // fabsf narrows the double expression del-1.0 back to float.
        if (((del as f64 - 1.0) as f32).abs() as f64) < EPS {
            break;
        }
        i += 1;
    }
    if i > ITMAX {
        fprintf!(stderr(), "a too large, ITMAX too small in gcf\n");
        exit(EXIT_FAILURE);
    }
    *gammcf = ((-(x as f64) + a as f64 * (x as f64).ln() - *gln as f64).exp() * h as f64) as f32;
}

/// Lower incomplete gamma probability from the original series approximation.
///
/// C prototype: `void gser(float *gamser, float a, float x, float *gln);`
///
/// Writes `gammln(a)` to `gln`. Negative `x` terminates with status 1;
/// failure to converge prints a diagnostic and leaves `gamser` unchanged.
pub fn gser(gamser: &mut f32, a: f32, x: f32, gln: &mut f32) {
    const ITMAX: i32 = 100;
    const EPS: f64 = 3.0e-7;

    *gln = gammln(a);
    if x <= 0.0 {
        if x < 0.0 {
            fprintf!(stderr(), "x less than 0 in routine gser\n");
            exit(EXIT_FAILURE);
        }
        *gamser = 0.0;
    } else {
        let mut ap = a;
        let mut sum = (1.0 / a as f64) as f32;
        let mut del = sum;
        for _n in 1..=ITMAX {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (del.abs() as f64) < sum.abs() as f64 * EPS {
                *gamser = (sum as f64
                    * (-(x as f64) + a as f64 * (x as f64).ln() - *gln as f64).exp())
                    as f32;
                return;
            }
        }
        fprintf!(stderr(), "a too large, ITMAX too small in routine gser\n");
    }
}

/// Fit `y = a + b*x`, returning the same six outputs as the C routine.
///
/// C prototype:
/// `void linfit(float *x, float *y, int ndata, float *sig, int mwt,
/// float *a, float *b, float *siga, float *sigb, float *chi2, float *q);`
///
/// `x.len()` supplies `ndata`; `y` must have the same length. For nonzero
/// `mwt`, `sig` must also have that length; otherwise it is unused and may
/// be empty. The original expression `1.0/sig[i]*sig[i]` is intentional.
/// Degenerate inputs retain C's floating-point NaNs/infinities. Weighted
/// fits with at most two samples terminate through `gammq`, as in C.
#[allow(clippy::too_many_arguments)]
pub fn linfit(
    x: &[f32],
    y: &[f32],
    sig: &[f32],
    mwt: i32,
    a: &mut f32,
    b: &mut f32,
    siga: &mut f32,
    sigb: &mut f32,
    chi2: &mut f32,
    q: &mut f32,
) {
    assert_eq!(x.len(), y.len());
    if mwt != 0 {
        assert_eq!(x.len(), sig.len());
    }
    let ndata = i32::try_from(x.len()).expect("ndata exceeds C int range");
    let mut sx = 0.0_f32;
    let mut sy = 0.0_f32;
    let mut st2 = 0.0_f32;
    let mut ss;

    *b = 0.0;
    if mwt != 0 {
        ss = 0.0;
        for i in 0..x.len() {
            let wt = (1.0 / sig[i] as f64 * sig[i] as f64) as f32;
            ss += wt;
            sx += x[i] * wt;
            sy += y[i] * wt;
        }
    } else {
        for i in 0..x.len() {
            sx += x[i];
            sy += y[i];
        }
        ss = ndata as f32;
    }
    let sxoss = sx / ss;
    if mwt != 0 {
        for i in 0..x.len() {
            let t = (x[i] - sxoss) / sig[i];
            st2 += t * t;
            *b += t * y[i] / sig[i];
        }
    } else {
        for i in 0..x.len() {
            let t = x[i] - sxoss;
            st2 += t * t;
            *b += t * y[i];
        }
    }
    *b /= st2;
    *a = (sy - sx * *b) / ss;
    *siga = (((1.0 + (sx * sx / (ss * st2)) as f64) / ss as f64).sqrt()) as f32;
    *sigb = (1.0 / st2 as f64).sqrt() as f32;
    *chi2 = 0.0;
    if mwt == 0 {
        for i in 0..x.len() {
            *chi2 += (y[i] - *a - *b * x[i]) * (y[i] - *a - *b * x[i]);
        }
        *q = 1.0;
        let sigdat = ((*chi2 / (ndata - 2) as f32) as f64).sqrt() as f32;
        *siga *= sigdat;
        *sigb *= sigdat;
    } else {
        for i in 0..x.len() {
            *chi2 += ((y[i] - *a - *b * x[i]) / sig[i]) * ((y[i] - *a - *b * x[i]) / sig[i]);
        }
        *q = gammq(
            (0.5 * (ndata - 2) as f64) as f32,
            (0.5 * *chi2 as f64) as f32,
        );
    }
}
