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

//! PDFHISTOGRAM - histogram of a probability density function.
//! Port of par/main/statistics/pdfhistogram.c (Revision 1.3).
//! Credits: Reginald H. Beardsley; Copyright 2006 Exploration Software Consultants Inc.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " PDFHISTOGRAM - generate a HISTOGRAM of the Probability Density function",
    " \t\t\t\t\t\t\t\t\t",
    "  pdfhistogram < stdin > sdtout [Required params] (Optional params)\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters: \t\t\t\t\t\t\t",
    " ix=\t\tcolumn containing X variable\t\t\t\t",
    " iy=\t\tcolumn containing Y variable \t\t\t\t",
    " min_x=\tminimum X bin\t\t\t\t\t\t",
    " max_x=\tmaximum X bin\t\t\t\t\t\t",
    " min_y=\tminimum Y bin\t\t\t\t\t\t",
    " max_y=\tmaximum Y bin\t\t\t\t\t\t",
    " logx=0\t=1 use logarithmic scale for X axis\t\t\t",
    " logy=0\t=1 use logarithmic scale for Y axis\t\t\t",
    " norm=\tselected normalization type \t\t\t\t\t",
    "\t\tsqrt\t- bin / sqrt( xnct*ycnt) \t\t\t",
    "\t\tavg_cnt\t- 0.5* bin / (xcnt + ycnt) \t\t\t",
    "\t\tavg_sum\t- (bin / xcnt + bin / ycnt ) / 2 \t\t",
    "\t\txcnt\t- bin / xcnt \t\t\t\t\t",
    "\t\tycnt\t- bin / ycnt \t\t\t\t\t",
    "\t\tlog\t- log(bin) \t\t\t\t\t",
    "\t\ttotal\t- bin / total \t\t\t\t\t",
    " Optional parameters: \t\t\t\t\t\t\t",
    "  nx=100\t- number of X bins\t\t\t\t\t",
    "  ny=100\t- number of Y bins \t\t\t\t\t",
    "  ir=\t- column containing reject variable \t\t\t\t",
    "  rmin=\t- reject values below rmin \t\t\t\t",
    "  rmax=\t- reject values above rmax \t\t\t\t",
    "\t\tNOTE: only one, rmin or rmax, may be used \t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " PDFHISTOGRAM creates a 2D histogram representing the probability density",
    " function of the input data. The output is in the form of a binary array",
    " that can then be plotted via ximage.\t\t\t\t\t",
    " Commandline options allow selecting any of several normalizations\t",
    " to apply to the distributions.\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut ix, mut iy, mut ir) = (0, 0, -1);
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let (mut rmin, mut rmax) = (0.0f32, 0.0f32);
    let (mut logx, mut logy, mut nx, mut ny) = (0, 0, 100, 100);
    if getparint("ix", &mut ix) == 0 { err!("Missing X column index, ix="); }
    if getparint("iy", &mut iy) == 0 { err!("Missing Y column index, ix="); }
    getparint("ir", &mut ir);
    if ir > 0 && getparfloat("rmax", &mut rmax) == 0 && getparfloat("rmin", &mut rmin) == 0 { err!("Failed to specify rejection threshold"); }
    if getparfloat("min_x", &mut min_x) == 0 { err!("min_x= not specified"); }
    if getparfloat("max_x", &mut max_x) == 0 { err!("max_x= not specified"); }
    if getparfloat("min_y", &mut min_y) == 0 { err!("min_y= not specified"); }
    if getparfloat("max_y", &mut max_y) == 0 { err!("max_y= not specified"); }
    getparint("logx", &mut logx); getparint("logy", &mut logy);
    getparint("nx", &mut nx); getparint("ny", &mut ny);
    let mut norm = String::new();
    if getparstring("norm", &mut norm) == 0 { norm = "total".to_owned(); }
    if !["sqrt", "avg_cnt", "avg_sum", "xcnt", "ycnt", "log", "total"].contains(&norm.as_str()) { err!("Unknown normalization"); }
    checkpars();
    if logx != 0 {
        if min_x < 0.0 || max_x < 0.0 { err!("Cannot use logscale for X"); }
        min_x = (min_x as f64).ln() as f32; max_x = (max_x as f64).ln() as f32;
    }
    if logy != 0 {
        if min_y < 0.0 || max_y < 0.0 { err!("Cannot use logscale for Y"); }
        min_y = (min_y as f64).ln() as f32; max_y = (max_y as f64).ln() as f32;
    }
    let dx = (max_x - min_x) / nx as f32;
    let dy = (max_y - min_y) / ny as f32;
    let (nxu, nyu) = (nx as usize, ny as usize);
    let mut bin = vec![vec![0.0f32; nxu]; nyu];
    let mut xcnt = vec![0i32; nxu];
    // The C code reads ycnt[i] in YCNT mode, even when i exceeds ny.
    let mut ycnt = vec![0i32; nyu.max(nxu)];
    let mut total = 0i32;
    let mut v = [0.0f32; 128];
    while let Some(buf) = fgets(1024, &stdin()) {
        for (k, token) in buf.split(|c| c == ' ' || c == '\t' || c == '\n').filter(|s| !s.is_empty()).enumerate() {
            if k >= v.len() { break; }
            let (number, used) = strtod(token);
            if used != 0 {
                // sscanf("%f") rounds decimal input directly to float, not via double.
                v[k] = token[..used].parse::<f32>().unwrap_or(number as f32);
            }
        }
        let mut x = v[(ix - 1) as usize]; let mut y = v[(iy - 1) as usize];
        if logx != 0 && x > 0.0 { x = (x as f64).ln() as f32; }
        if logy != 0 && y > 0.0 { y = (y as f64).ln() as f32; }
        let i = ((x - min_x) / dx) as i32; let j = ((y - min_y) / dy) as i32;
        let accepted = if ir > 0 {
            let r = v[(ir - 1) as usize];
            (rmin != 0.0 && r > rmin) || (rmax != 0.0 && r < rmax)
        } else { true };
        if accepted && i > 0 && i < nx && j > 0 && j < ny {
            let (i, j) = (i as usize, j as usize);
            bin[j][i] += 1.0; xcnt[i] += 1; ycnt[j] += 1; total += 1;
        }
    }
    let mut clip = 0.0f32;
    for i in 0..nxu { for j in 0..nyu {
        let b = bin[j][i];
        if b > 0.0 {
            let xok = (xcnt[i] as f64) > 1e-4 * total as f64;
            let yok = (ycnt[j] as f64) > 1e-4 * total as f64;
            bin[j][i] = if norm == "sqrt" && xok && yok {
                (b as f64 / ((xcnt[i] as f32 * ycnt[j] as f32) as f64).sqrt()) as f32
            } else if norm == "avg_cnt" && xok && yok {
                ((b / xcnt[i] as f32 + b / ycnt[j] as f32) as f64 * 0.5) as f32
            } else if norm == "avg_sum" && xok && yok {
                (2.0 * b as f64 / (xcnt[i] as f32 + ycnt[j] as f32) as f64) as f32
            } else if norm == "xcnt" && xok { b / xcnt[i] as f32
            } else if norm == "ycnt" && (ycnt[i] as f64) > 1e-4 * total as f64 { b / ycnt[i] as f32
            } else if norm == "log" { (b as f64).ln() as f32
            } else if norm == "total" { b / total as f32
            } else { 0.0 };
            if bin[j][i].is_nan() {
                fprintf!(stderr(), "Nan at ");
                let xx = min_x + i as f32 * dx; let yy = min_y + j as f32 * dy;
                if logx != 0 { fprintf!(stderr(), "%f ", (xx as f64).exp()); } else { fprintf!(stderr(), "%f ", xx); }
                if logy != 0 { fprintf!(stderr(), "%f ", (yy as f64).exp()); } else { fprintf!(stderr(), "%f ", yy); }
                fprintf!(stderr(), "\n");
            }
            if bin[j][i] > clip { clip = bin[j][i]; }
        } else { bin[j][i] = 0.0; }
    }}
    clip = (clip as f64 * 0.9) as f32;
    let out = stdout();
    fwrite(&[nx], &out); fwrite(&[min_x, dx], &out);
    fwrite(&[ny], &out); fwrite(&[min_y, dy, 0.0f32, clip], &out);
    for row in &bin { fwrite(row, &out); }
    exit(CWP_Exit());
}
