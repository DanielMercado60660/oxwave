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

//! Port of psplot/lib/psWiggle.c.
//! Credits: Dave Hale, Craig Artley, and Paul Michaels.
#![allow(non_snake_case)]
use oxwave_curve_clip::yclip;
use oxwave_f_basic::*;

fn path(xc: &[f32], zc: &[f32], base: f32, xb: f32, xs: f32, zs: f32) {
    let mut xl = xc[0];
    let mut zl = base;
    for (&x, &z) in xc.iter().zip(zc) {
        if z != base {
            if zl == base { moveto(xb + xl * xs, 0.0 + zl * zs); }
            lineto(xb + x * xs, 0.0 + z * zs);
        } else if zl != base { lineto(xb + x * xs, 0.0 + z * zs); }
        xl = x; zl = z;
    }
    if zl != base { lineto(xb + xl * xs, 0.0 + base * zs); }
    fill();
}

/// Draw a clipped wiggle trace with optional filled areas.
///
/// C prototype: `void psWiggle(int n, float z[], float zmin, float zmax,
/// float zbase, float yzmin, float yzmax, float xfirst, float xlast,
/// int ifill, const char *tracecolor);`
///
/// The clip limits must be ordered, and for filled traces the base must lie
/// between them, as required by `yclip`.
#[allow(clippy::too_many_arguments)]
pub fn psWiggle(z: &[f32], zmin: f32, zmax: f32, zbase: f32,
    yzmin: f32, yzmax: f32, xfirst: f32, xlast: f32, ifill: i32, tracecolor: &str) {
    let n = z.len();
    let xs = if n > 1 { (xlast - xfirst) / (n - 1) as f32 } else { 1.0 };
    let zs = if zmax != zmin { (yzmax - yzmin) / (zmax - zmin) } else { 1.0 };
    let zb = if zmax != zmin { yzmin - zmin * zs } else { (0.5 * (yzmin + yzmax) as f64) as f32 };
    let zeps = (0.001 * (if zmax > zmin { zmax - zmin } else { zmin - zmax }) as f64) as f32;
    gsave(); translate(0.0, zb); newpath(); setlinejoin(2);
    let mut xc = [0.0; 400]; let mut zc = [0.0; 400];
    let mut k1 = 0;
    let mut k2 = 0;
    while k2 < n.saturating_sub(1) {
        k2 = (k1 + 199).min(n - 1);
        let segment = &z[k1..=k2];
        if ifill != 0 {
            if ifill.unsigned_abs() >= 2 {
                let (lo, hi) = if ifill < 0 { (zbase, zmax) } else { (zmin, zbase) };
                let nc = yclip(1.0, k1 as f32, segment, lo, hi, &mut xc, &mut zc) as usize;
                let mut shade = (1.0 - 0.2 * ifill.unsigned_abs() as f32 as f64) as f32;
                if shade < 0.0 { shade = 0.0; }
                setrgbcolor(shade, shade, shade);
                path(&xc[..nc], &zc[..nc], zbase, xfirst, xs, zs);
                setcolor(tracecolor);
            }
            let (lo, hi) = if ifill > 0 { (zbase, zmax) } else { (zmin, zbase) };
            let nc = yclip(1.0, k1 as f32, segment, lo, hi, &mut xc, &mut zc) as usize;
            path(&xc[..nc], &zc[..nc], zbase, xfirst, xs, zs);
        }
        let nc = yclip(1.0, k1 as f32, segment, zmin, zmax, &mut xc, &mut zc) as usize;
        moveto(xfirst + xc[0] * xs, 0.0 + zc[0] * zs);
        let mut il = 0; let mut zl = zc[0];
        for ic in 1..nc {
            let dz = zc[ic] - zl;
            if dz < -zeps || dz > zeps {
                if il != ic - 1 { lineto(xfirst + xc[ic-1] * xs, 0.0 + zc[ic-1] * zs); }
                lineto(xfirst + xc[ic] * xs, 0.0 + zc[ic] * zs);
                il = ic; zl = zc[ic];
            }
        }
        if il != nc - 1 { lineto(xfirst + xc[nc-1] * xs, 0.0 + zc[nc-1] * zs); }
        stroke(); k1 += 198;
    }
    grestore();
}
