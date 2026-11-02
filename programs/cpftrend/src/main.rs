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

//! CPFTREND - cumulative probability picks, par/main/statistics/cpftrend.c (Revision 1.4).
//! Credits: Reginald H. Beardsley, rhb@acm.org.
//! Copyright 2006 Exploration Software Consultants Inc.

use oxwave_core::prelude::*;
use std::ffi::{c_char, c_int, CString};

unsafe extern "C" {
    fn sscanf(s: *const c_char, format: *const c_char, ...) -> c_int;
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    "   CPFTREND - generate picks of the Cumulate Probability Function \t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "  ix=      - column containing X variable\t\t\t\t",
    "  iy=      - column containing Y variable\t\t\t\t",
    "  min_x=   - minimum X bin\t\t\t\t\t\t",
    "  max_x=   - maximum X bin\t\t\t\t\t\t",
    "  min_y=   - minimum Y bin \t\t\t\t\t\t",
    "  max_y=   - maximum Y bin\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "  nx=100    - number of X bins \t\t\t\t\t",
    "  ny=100    - number of Y bins \t\t\t\t\t",
    "  logx=0   - =1 use logarithmic scale for X axis\t\t\t",
    "  logy=0   - =1 use logarithmic scale for Y axis\t\t\t",
    "  ir=       - column containing reject variable \t\t\t",
    "  rmin=     - reject values below rmin \t\t\t\t",
    "  rmax=     - reject values above rmax \t\t\t\t",
    "              NOTE: only one, rmin or rmax, may be used\t\t",
    " NOTES:\t\t\t\t\t\t\t\t",
    "  cpftrend makes picks on the 2D cumulate representing the\t\t",
    "  probability density function of the input data.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   Commandline options allow selecting any of several normalizations\t",
    "   to apply to the distributions.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut ix, mut iy, mut ir) = (0, 0, -1);
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let (mut rmin, mut rmax) = (0.0f32, 0.0f32);
    let (mut nx, mut ny, mut logx, mut logy) = (100, 100, 0, 0);
    if getparint("ix", &mut ix) == 0 { err!("Missing X column index, ix="); }
    if getparint("iy", &mut iy) == 0 { err!("Missing Y column index, ix="); }
    getparint("ir", &mut ir);
    if ir > 0 && getparfloat("rmax", &mut rmax) == 0 && getparfloat("rmin", &mut rmin) == 0 {
        err!("Failed to specify rejection threshold");
    }
    if getparfloat("min_x", &mut min_x) == 0 { err!("min_x= not specified"); }
    if getparfloat("max_x", &mut max_x) == 0 { err!("max_x= not specified"); }
    if getparfloat("min_y", &mut min_y) == 0 { err!("min_y= not specified"); }
    if getparfloat("max_y", &mut max_y) == 0 { err!("max_y= not specified"); }
    getparint("nx", &mut nx);
    getparint("ny", &mut ny);
    getparint("logx", &mut logx);
    getparint("logy", &mut logy);
    checkpars();
    if logx != 0 {
        if min_x < 0.0 || max_x < 0.0 { err!("Cannot use logscale for X"); }
        min_x = (min_x as f64).ln() as f32;
        max_x = (max_x as f64).ln() as f32;
    }
    if logy != 0 {
        if min_y < 0.0 || max_y < 0.0 { err!("Cannot use logscale for Y"); }
        min_y = (min_y as f64).ln() as f32;
        max_y = (max_y as f64).ln() as f32;
    }
    let dx = (max_x - min_x) / nx as f32;
    let dy = (max_y - min_y) / ny as f32;
    let mut bin = vec![vec![0.0f32; ny as usize]; nx as usize];
    let mut xcnt = vec![0i32; nx as usize];
    let mut a = vec![0.0f32; nx as usize];
    let mut b = vec![0.0f32; nx as usize];
    let mut total = 0i32;
    let mut v = [0.0f32; 128];
    while let Some(buf) = fgets(1024, &stdin()) {
        // strtok sees the first NUL as the end of the C buffer.
        let buf = buf.split('\0').next().unwrap_or("");
        for (k, token) in buf.split(|c| c == ' ' || c == '\t' || c == '\n').filter(|s| !s.is_empty()).enumerate() {
            if k >= v.len() { break; }
            let token = CString::new(token).unwrap();
            // C scans directly into float; parsing as double first can round twice.
            unsafe { sscanf(token.as_ptr(), c"%f".as_ptr(), &mut v[k] as *mut f32); }
        }
        let mut x = v[(ix - 1) as usize];
        let mut y = v[(iy - 1) as usize];
        if logx != 0 && x > 0.0 { x = (x as f64).ln() as f32; }
        if logy != 0 && y > 0.0 { y = (y as f64).ln() as f32; }
        let i = ((x - min_x) / dx) as i32;
        let j = ((y - min_y) / dy) as i32;
        let accept = if ir > 0 {
            let r = v[(ir - 1) as usize];
            (rmin != 0.0 && r > rmin || rmax != 0.0 && r < rmax) && i > 0 && i < nx && j > 0 && j < ny
        } else { i > 0 && i < nx && j > 0 && j < ny };
        if accept {
            let (i, j) = (i as usize, j as usize);
            bin[i][j] += 1.0;
            xcnt[i] += 1;
            a[i] += y * y;
            b[i] += y;
            total += 1;
        }
    }
    for i in 0..nx as usize {
        if (xcnt[i] as f64) > 0.0001 * total as f64 {
            let mut sum = 0.0f32;
            let xpos = min_x + i as f32 * dx;
            if logx != 0 { printf!("%f ", (xpos as f64).exp()); }
            else { printf!("%f ", xpos as f64); }
            printf!("%d ", xcnt[i]);
            let mut t = [false; 5];
            for j in 0..ny as usize {
                sum += bin[i][j];
                let m = (100.0 * sum as f64 / xcnt[i] as f64) as f32;
                let thresholds = [2.28, 15.87, 50.00, 84.13, 97.72];
                for n in 0..5 {
                    if !t[n] && m as f64 >= thresholds[n] {
                        t[n] = true;
                        let ypos = min_y + j as f32 * dy;
                        if logy != 0 { printf!("%f ", (ypos as f64).exp()); }
                        else { printf!("%f ", ypos as f64); }
                        break;
                    }
                }
            }
            let denom = xcnt[i].wrapping_mul(xcnt[i]).wrapping_sub(1);
            let sigma = (((xcnt[i] as f32 * a[i] - b[i] * b[i]) / denom as f32) as f64).sqrt() as f32;
            if logy != 0 && sigma > 0.0 { printf!("%f ", (sigma as f64).ln()); }
            else { printf!("%f ", sigma as f64); }
            printf!("\n");
        }
    }
    exit(CWP_Exit());
}
