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

//! Port of cwp/lib/printerplot.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89 and 12/26/89.

use oxwave_core::{cio::CFile, err, fprintf};

/// Print a one-dimensional array. The array must not be empty.
///
/// C prototype: `void pp1d (FILE *fp, char *title, int lx, int ifx, float x[]);`
pub fn pp1d(fp: &CFile, title: &str, ifx: i32, x: &[f32]) {
    plot(fp, title, ifx, x);
}

/// Print a one-dimensional array, indexed from zero. The array must not be empty.
///
/// C prototype: `void pplot1 (FILE *fp, char *title, int nx, float ax[]);`
pub fn pplot1(fp: &CFile, title: &str, ax: &[f32]) {
    plot(fp, title, 0, ax);
}

fn plot(fp: &CFile, title: &str, ifx: i32, x: &[f32]) {
    if x.is_empty() || x.iter().any(|v| !v.is_finite()) {
        err!("printerplot.c: array must be nonempty and finite");
    }
    let title = title.split('\0').next().unwrap();
    fprintf!(fp, "\n");
    fprintf!(fp, "%s\n", title);
    let mut xmin = x[0];
    let mut xmax = x[0];
    for &value in &x[1..] {
        xmin = if xmin < value { xmin } else { value };
        xmax = if xmax > value { xmax } else { value };
    }
    fprintf!(fp, "minimum = %g\n", xmin);
    fprintf!(fp, "maximum = %g\n", xmax);
    let xscale: f32 = if xmin == xmax { 1.0 } else { 60.0f32 / (xmax - xmin) };
    // C's float-to-int conversion is undefined if this scale overflows.
    // Retain its float arithmetic everywhere else, and use a finite scale
    // for tiny ranges so each row remains within the 60-column plot.
    let wide_scale = if xscale.is_infinite() {
        Some(60.0 / (xmax as f64 - xmin as f64))
    } else {
        None
    };
    let (ibase, xbase) = if xmin < 0.0 && xmax < 0.0 {
        (60i32, xmax)
    } else if xmin < 0.0 && xmax >= 0.0 {
        ((((0.0f64 - xmin as f64) * wide_scale.unwrap_or(xscale as f64)) as i32).clamp(0, 60), 0.0f32)
    } else {
        (0, xmin)
    };
    for (ix, &value) in x.iter().enumerate() {
        let v = match wide_scale {
            Some(scale) => (value as f64 - xbase as f64) * scale,
            None => ((value - xbase) * xscale) as f64,
        };
        let rounded = if v > 0.0 { (v + 0.5) as i32 } else { (v - 0.5) as i32 };
        let icx = ibase.wrapping_add(rounded).clamp(0, 60);
        fprintf!(fp, "%4d %13.6e ", ifx.wrapping_add(ix as i32), value);
        for _ in 0..ibase.min(icx) { fprintf!(fp, " "); }
        for _ in ibase.min(icx)..=ibase.max(icx) { fprintf!(fp, "*"); }
        fprintf!(fp, "\n");
    }
}
