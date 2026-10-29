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
//! SUHISTOGRAM - create histogram of input amplitudes.
//! Port of su/main/attributes_parameter_estimation/suhistogram.c (Revision 1.10).
//! Credits: Reginald H. Beardsley, 2006, rhb@acm.org.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUHISTOGRAM - create histogram of input amplitudes\t\t",
    "\t\t\t\t\t\t\t\t",
    "    suhistogram <in.su >out.dat\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " min=\t\tminimum bin \t\t\t\t\t",
    " max=\t\tmaximum bin \t\t\t\t\t",
    " bins=\t\tnumber of bins\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters\t\t\t\t\t\t",
    " trend=0\t=0 1-D histogram\t\t\t\t",
    "\t   =1 2-D histogram picks on cumulate\t\t\t",
    "\t   =2 2-D histogram in trace format\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " clip=     threshold value to drop outliers\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " dt=\tsample rate in feet or milliseconds.  Defaults  to\t",
    "    \ttr.dt*1e-3\t\t\t\t\t  \t",
    " datum=  header key to get datum shift if desired (e.g. to\t",
    "\t hang from water bottom)\t\t\t    \t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " trend=0 produces a two column ASCII output for use w/ gnuplot.",
    " Extreme values are counted in the end bins.\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " trend=1 produces a 6 column ASCII output for use w/ gnuplot   ",
    " The columns are time/depth and picks on the cumulate\t\t",
    " at 2.28%, 15.87%, 50%, 84.13% & 97.72% of the total points    ",
    " corresponding to the median and +- 1 or 2 standard deviations ",
    " for a Gaussian distribution.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " trend=2 produces an SU trace panel w/ one trace per bin that  ",
    " can be displayed w/ suximage, etc.\t\t\t\t",
    "\t\t\t\t\t\t\t   \t",
    " Example for plotting with xgraph:\t\t\t\t",
    " suhistogram < data.su min=MIN max=MAX bins=BINS |\t\t",
    " a2b n1=2 | xgraph n=BINS nplot=1\t\t\t \t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("Can't get first trace \n"); }
    let mut min = 0.0f32;
    let mut max = 0.0f32;
    let mut bins = 0i32;
    let mut trend = 0i32;
    let mut datum = String::new();
    let mut clip = 1e38f32;
    let mut dt = 0.0f32;
    if getparfloat("min", &mut min) == 0 { err!("min must be specified\n"); }
    if getparfloat("max", &mut max) == 0 { err!("max must be specified\n"); }
    if getparint("bins", &mut bins) == 0 { err!("bins must be specified\n"); }
    getparint("trend", &mut trend);
    let has_datum = getparstring("datum", &mut datum) != 0;
    getparfloat("clip", &mut clip);
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 * 1e-3) as f32; }
    checkpars();
    if bins < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    // The original indexes outside its allocation for a zero-bin 1-D histogram.
    if bins == 0 && trend == 0 { err!("bins must be positive"); }
    let bin = (max - min) / bins as f32;
    let ns = tr.ns as usize;
    let mut cnt = vec![0.0f32; ns];
    let mut a = vec![0.0f32; ns];
    let mut b = vec![0.0f32; ns];
    if trend != 0 {
        let mut histogram = vec![0i32; ns * bins as usize];
        loop {
            let fval = if has_datum {
                let mut val = Value::default();
                gethdval(&tr, &datum, &mut val);
                vtof(hdtype(&datum), val)
            } else { 0.0 };
            for i in 0..tr.ns as usize {
                if tr.data[i] < clip {
                    let j = ((tr.data[i] - min) / bin) as i32;
                    let k = ((i as f32 * dt - fval) / dt) as i32;
                    if j >= 0 && j < bins && k >= 0 && k < tr.ns as i32 {
                        let k = k as usize;
                        let idx = k * bins as usize + j as usize;
                        histogram[idx] = histogram[idx].wrapping_add(1);
                        a[k] += tr.data[i] * tr.data[i];
                        b[k] += tr.data[i];
                        cnt[k] += 1.0;
                    }
                }
            }
            if gettr(&mut tr) == 0 { break; }
        }
        if trend == 1 {
            for i in 0..tr.ns as usize {
                if cnt[i] > 0.0 {
                    let (mut p1, mut p2, mut p3, mut p4, mut p5) = (false, false, false, false, false);
                    let mut sum = 0i32;
                    printf!("%f ", (i as i32 * tr.dt as i32) as f64 / 1000.0);
                    for j in 0..bins as usize {
                        sum = sum.wrapping_add(histogram[i * bins as usize + j]);
                        let m = (100.0 * sum as f64 / cnt[i] as f64) as f32;
                        let x = min + j as f32 * bin;
                        if !p1 && m as f64 >= 2.28 { p1 = true; printf!("%f ", x); }
                        else if !p2 && m as f64 >= 15.87 { p2 = true; printf!("%f ", x); }
                        else if !p3 && m as f64 >= 50.00 { p3 = true; printf!("%f ", x); }
                        else if !p4 && m as f64 >= 84.13 { p4 = true; printf!("%f ", x); }
                        else if !p5 && m as f64 >= 97.72 { p5 = true; printf!("%f ", x); }
                    }
                    let c = (((cnt[i] * a[i] - b[i] * b[i]) / (cnt[i] * (cnt[i] - 1.0))) as f64).sqrt() as f32;
                    printf!("%f ", c);
                    printf!("\n");
                }
            }
        } else {
            for j in 0..bins as usize {
                for i in 0..tr.ns as usize {
                    tr.data[i] = if cnt[i] > 0.0 { histogram[i * bins as usize + j] as f32 / cnt[i] } else { histogram[i * bins as usize + j] as f32 };
                }
                tr.tracl = j as i32;
                fputtr(&stdout(), &mut tr);
            }
        }
    } else {
        let mut histogram = vec![0i32; bins as usize];
        loop {
            for i in 0..tr.ns as usize {
                let mut j = ((tr.data[i] - min) / bin) as i32;
                j = if j < 0 { 0 } else { j };
                j = if j > bins - 1 { bins - 1 } else { j };
                histogram[j as usize] = histogram[j as usize].wrapping_add(1);
                cnt[0] += 1.0;
            }
            if gettr(&mut tr) == 0 { break; }
        }
        for i in 0..bins as usize {
            printf!("%15f ", min + i as f32 * bin);
            printf!("%15f ", histogram[i] as f32 / cnt[0]);
            printf!("\n");
        }
    }
    exit(CWP_Exit());
}
