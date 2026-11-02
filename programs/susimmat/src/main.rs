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

//! SUSIMMAT - port of su/main/operations/susimmat.c, revision 1.3.
//! Copyright (c) Colorado School of Mines, 2015. All rights reserved.
//! Credits: Chris Liner, U Arkansas; John Stockwell, CWP.
//! References: Liner and Clapp, The Leading Edge 23.11 (2004), 1146-1150;
//! Liner and Clapp, GEOPHYSICS 69.6 (2004), 1552-1559, doi:10.1190/1.1836828.

use oxwave_core::prelude::*;
use oxwave_linear_regression::linear_regression;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUSIMAT - Correlation similarity matrix for two traces.\t\t",
    " \t\t Output is zero lag of cross-correlation of traces,\t",
    " \t\t or linear regression correlation coefficient.\t\t",
    " \t\t Horizontal axis is time of trace 1, vertical is trace 2.",
    " \t\t\t\t\t\t\t\t\t",
    " susimmat <data12 sufile=data2 >dataout \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " sufile=\t\tfile containing SU traces to use as filter\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " panel=0\t\tuse only the first trace of sufile as filter\t",
    " \t\t      \t=1 compute sim matrix trace by trace an entire  ",
    "                        gather\t\t\t\t\t\t",
    " mt=21\t\t\toperator window length (odd integer)\t\t",
    " eps=1e-3\t\tstability parameter\t\t\t\t",
    " taper=0\t\tno taper to data fragments\t\t\t",
    " \t\t\t=1 apply exponential taper (1/e at ends)\t",
    " method=1\t\tuse xcorrelation as similarity meausure\t\t",
    " \t\t\t=2 same but normalized by (rms+eps)\t\t",
    " \t\t\t=3 use linear regression CC\t\t\t",
    " \t\t\t=4 use mt-dimensional vector angle\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " EXAMPLE PROCESSING SEQUENCES:\t\t\t\t\t\t",
    "   1. Look for all possible alignments of OBC P and Z data  \t\t",
    "   \tsusimmat < P_Z.su sufile=OBC_P.su  mt=71 > P_Zxcor.su\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note:  xcor window is collapsed as needed to compute edge values\t",
    " It is quietly assumed that the time sampling interval on the  single  ",
    " trace\t\t\t\t\t\t\t\t\t",
    " and the output traces is the same as that on the traces in the input  ",
    " file.  \t\t\t\t\t\t\t\t",
    " The sufile may actually have more than one trace, but only the first  ",
    " trace is used when panel=0. When panel=1 the number of traces in the ",
    " sufile MUST be  the same as the number of traces in the input.\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut mt = 21;
    getparint("mt", &mut mt);
    if mt != 21 && mt != 1 {
        mt = mt / 2 * 2 + 1;
        if verbose != 0 { warn!("Using mt=%i", mt); }
    }
    let mtorig = mt;
    let mut method = 1;
    getparint("method", &mut method);
    if !(1..=4).contains(&method) { err!("Bad method"); }
    if verbose != 0 { warn!("method = %i", method); }
    if mt < 5 && method == 2 { err!("Need more than 5 points for linear regression"); }
    let mut taper = 0;
    getparint("taper", &mut taper);
    let mut eps = 1e-3f32;
    getparfloat("eps", &mut eps);
    let mut intrace = Trace::new();
    let mut sutrace = Trace::new();
    let mut outtrace = Trace::new();
    if gettr(&mut intrace) == 0 { err!("can't get first trace"); }
    let nt = intrace.ns as usize;
    let dt = intrace.dt;
    // C's high-edge window can read beyond the allocated sample vectors.
    let pad = (mtorig.max(1) as usize) + 1;
    let mut trace = vec![0.0f32; nt + pad];
    let mut data = vec![0.0f32; nt + pad];
    let mut mat = vec![0.0f32; nt * nt];
    // ealloc1float converts a negative mt to size_t and its malloc fails.
    if mt < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut x = vec![0.0f32; mt.max(1) as usize];
    let mut y = vec![0.0f32; mt.max(1) as usize];
    let mut sufile = String::new();
    mustgetparstring("sufile", &mut sufile);
    let sufilefp = efopen(&sufile, "r");
    fgettr(&sufilefp, &mut sutrace);
    data[..nt].copy_from_slice(&sutrace.data[..nt]);
    loop {
        mat.fill(0.0);
        trace[..nt].copy_from_slice(&intrace.data[..nt]);
        for icol in 0..nt {
            for irow in 0..nt {
                let mut hw = (mtorig - 1) / 2;
                let (row, col, n) = (irow as i32, icol as i32, nt as i32);
                let (mut x0, mut rms1, mut rms2) = (0.0f32, 0.0f32, 0.0f32);
                if row <= hw || col <= hw {
                    if row <= col { hw = row; } else { hw = col; }
                }
                if row >= n - hw || col >= n - hw {
                    if row <= hw || col <= hw {
                        if row <= col { hw = n - col; } else { hw = n - row; }
                    }
                }
                mt = 2 * hw + 1;
                for j in 0..mt {
                    // Integer division by zero on arm64 produces zero.
                    let q = (hw - j + 1).checked_div(hw + 1).unwrap_or(0);
                    let alpha = (4 * q * q) as f32;
                    let a = trace[(col - hw + j) as usize];
                    let b = data[(row - hw + j) as usize];
                    if taper == 0 { x0 += a * b; }
                    if taper == 1 { x0 = (x0 as f64 + (a * b) as f64 * (-(alpha as f64)).exp()) as f32; }
                    rms1 += a * a;
                    rms2 += b * b;
                    x[j as usize] = a;
                    y[j as usize] = b;
                }
                mat[icol * nt + irow] = match method {
                    1 => x0,
                    2 | 4 => ((x0 as f64 / (((rms1 * rms2) as f64).sqrt() + eps as f64)).acos() * 180.0 / std::f64::consts::PI) as f32,
                    3 => { let mut coeff = [0.0; 4]; linear_regression(&y[..mt as usize], &x[..mt as usize], &mut coeff); coeff[2] },
                    _ => unreachable!(),
                };
            }
        }
        outtrace.trid = 1;
        outtrace.ns = nt as u16;
        outtrace.dt = dt;
        for icol in 0..nt {
            outtrace.data[..nt].copy_from_slice(&mat[icol * nt..(icol + 1) * nt]);
            outtrace.tracl = (icol + 1) as i32;
            puttr(&mut outtrace);
        }
        if gettr(&mut intrace) == 0 { break; }
    }
    exit(CWP_Exit());
}
