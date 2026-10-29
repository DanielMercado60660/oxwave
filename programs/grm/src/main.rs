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

// Original source: Copyright (c) Colorado School of Mines, 2011.
// All rights reserved.

//! GRM - Generalized Reciprocal refraction analysis for a single layer.
//! Port of par/main/refraction/grm.c (Revision 1.4).
//!
//! Credits (from the C source):
//! CWP: Steven D. Sheaffer
//! D. Palmer, "The Generalized Reciprocal Method of Seismic
//! Refraction Interpretation", SEG, 1982.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " GRM - Generalized Reciprocal refraction analysis for a single layer\t",
    " \t\t\t\t\t\t\t\t\t",
    "     grm <stdin >stdout  [parameters]    \t\t \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " nt=\t\tNumber of arrival time pairs\t\t\t\t",
    " dx=\t\tGeophone spacing (m)\t\t\t\t\t",
    " v0=\t\tVelocity in weathering layer (m/s)\t\t\t",
    " abtime=\tIf set to 0, use last time as a-b, else give time (ms)  ",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " XY=      Value of XY if you want to override the optimum XY\t\t",
    "\t  algorithm in the program. If it is not an integer multiple of ",
    "\t dx, then it will be converted to the closest\t\t\t",
    "\t\t one.\t\t\t\t\t\t\t",
    "\tXYmax   Maximum offset distance allowed when searching for      ",
    "\t\toptimum XY (m)  (Default is 2*dx*10)\t\t\t",
    "\tdepthres  Size of increment in x during verical depth search(m) ",
    "\t\t  (Default is 0.5m)\t\t\t\t\t",
    " Input file:\t\t\t\t\t\t\t\t",
    "\t4 column ASCII - x,y, forward time, reverse time \t\t",
    " Output file:\t\t\t\t\t\t\t\t",
    "\t1) XYoptimum  \t\t\t\t\t\t\t",
    "\t2) apparent refractor velcocity\t\t\t\t\t",
    "\t3) x, y, z(x,y), y-z(x,y)\t\t\t\t\t",
    "\t\tz(x,y) = calculated (GRM) depth below (x y) \t\t",
    "\t\ty-z(x,y) = GRM depth subtracted from y - absolute depth ",
    "      .............\t\t\t\t\t\t\t",
    "      4) x, y, d(x,y), y-d(x,y), (error)  \t\t\t\t",
    "\t\td(x,y) = dip corrected depth estimate below (x,y)       ",
    "\t\ty-d(x,y) = dip corrected absolute depth \t\t",
    "\t\terror = estimated error in depth due only to the inexact",
    "\t\t      matching of tangents to arcs in dip estimate.\t",
    "\t\t\t\t\t\t\t\t\t",
    "      If the XY calculation is bypassed and XY specified, the values\t",
    "      used will precede 1) above.  XYoptimum will still be calculated\t",
    "      and displayed for reference.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t       \t",
    "      Uses average refactor velocity along interface.\t\t\t",
];

fn main() {
    let argv = args();
    let mut offset: i32 = 0;
    let mut nt: i32 = 0;
    let mut dx = 0.0f32;
    let mut v0 = 0.0f32;
    let mut abtime = 0.0f32;
    let mut xymax = 0.0f32;
    let mut xy = 0.0f32;
    let mut depthres = 0.0f32;
    let mut errd = 0.0f32;

    /* Initialize */
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get parameters */
    if getparint("nt", &mut nt) == 0 {
        err!("must specify nt");
    }
    if getparfloat("dx", &mut dx) == 0 {
        err!("must specify dx");
    }
    if getparfloat("v0", &mut v0) == 0 {
        err!("must specify v0");
    }
    if getparfloat("abtime", &mut abtime) == 0 {
        err!("must specify abtime");
    }
    if getparfloat("xy", &mut xy) == 0 {
        xy = 99999.9f32;
    }
    if getparfloat("xymax", &mut xymax) == 0 {
        xymax = 2.0f32 * dx * 10.0f32;
    }
    if getparfloat("depthres", &mut depthres) == 0 {
        depthres = 0.5f32;
    }

    checkpars();

    /* Turn XYmax into an integer # of samples on each side. */
    let jmax = ((xymax as f64 / 2.0) / dx as f64) as i32;

    /* Allocate memory to arrays (the original uses one-based indices). */
    let mut x = vec![0.0f32; nt.wrapping_add(1) as usize];
    let mut y = vec![0.0f32; nt.wrapping_add(1) as usize];
    let mut tf = vec![0.0f32; nt.wrapping_add(1) as usize];
    let mut tr = vec![0.0f32; nt.wrapping_add(1) as usize];
    let mut lap = vec![0.0f32; jmax.wrapping_add(1) as usize];
    let mut tv = vec![0.0f32; nt.wrapping_add(1) as usize];
    let mut zz = vec![0.0f32; nt.wrapping_add(1) as usize];
    let mut d = vec![0.0f32; nt.wrapping_add(1) as usize];

    /* Read arrival time data. scanf leaves previously read values intact
     * on a short row/EOF and stops the four conversions on a mismatch.
     * Values before the first successful scan are undefined in C; zero
     * initializes those otherwise indeterminate values here. */
    let input = stdin();
    let (mut r1, mut r2, mut r3, mut r4) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    for i in 1..=nt {
        if let Some(value) = fscanf_f32(&input) {
            r1 = value;
            if let Some(value) = fscanf_f32(&input) {
                r2 = value;
                if let Some(value) = fscanf_f32(&input) {
                    r3 = value;
                    if let Some(value) = fscanf_f32(&input) {
                        r4 = value;
                    }
                }
            }
        }
        let i = i as usize;
        x[i] = r1;
        y[i] = r2;
        tf[i] = (r3 as f64 * 0.001) as f32;
        tr[i] = (r4 as f64 * 0.001) as f32;
    }

    /* Set A-B time. Use the given abtime, or the last forward time. */
    let tab = if abtime != 0.0f32 {
        (abtime as f64 * 0.001) as f32
    } else {
        tf[nt as usize]
    };

    /* Generate the velocity analysis function for each half-offset j.
     * Use the same midpoint range at all offsets for equal lengths. */
    for j in 0..=jmax {
        for i in jmax + 1..=nt - jmax {
            let tay = tf[(i + j) as usize];
            let tbx = tr[(i - j) as usize];
            tv[i as usize] = (tay - tbx + tab) / 2.0f32;
        }

        /* Sum the absolute finite-difference Laplacian of tv. */
        lap[j as usize] = 0.0f32;
        for i in jmax + 2..=nt - jmax - 1 {
            let i = i as usize;
            let dlap = (1000.0 * (tv[i - 1] - 2.0f32 * tv[i] + tv[i + 1]) as f64 / (dx * dx) as f64)
                as f32;
            lap[j as usize] = (lap[j as usize] as f64 + (dlap as f64).abs()) as f32;
        }
    }

    /* The smoothest tv function has the smallest lap. Scanning backwards
     * means equal minima keep the larger offset. */
    let mut hold = 99999.0f32;
    for j in (0..=jmax).rev() {
        if lap[j as usize] < hold {
            offset = j;
            hold = lap[j as usize];
        }
    }
    printf!(
        "%s %5.1f \n",
        "Optimum XY(m): ",
        offset as f32 * dx * 2.0f32
    );

    if (xy as f64) < 99999.7 {
        offset = (xy / (2.0f32 * dx)) as i32;
        printf!("%s %5.1f \n", "XY Specified(m): ", xy);
        printf!(
            "%s %5.1f \n",
            "XY Used(m): ",
            offset.wrapping_mul(2) as f32 * dx
        );
    }

    /* Average the slopes before taking the reciprocal for Vapp. */
    for i in offset + 1..=nt - offset {
        let tay = tf[(i + offset) as usize];
        let tbx = tr[(i - offset) as usize];
        tv[i as usize] = (tay - tbx + tab) / 2.0f32;
    }

    let mut slope = 0.0f32;
    for i in offset + 1..=nt - offset - 1 {
        let dslope = (tv[(i + 1) as usize] - tv[i as usize]) / dx;
        slope += dslope;
    }
    let v1 = (nt - 2 * offset - 1) as f32 / slope;
    printf!("%s %5.1f \n", "Avg. Refractor Vapp: ", v1);

    /* Calculate time depth and depth at each usable geophone. */
    xy = offset.wrapping_mul(2) as f32 * dx;
    printf!("%s \n", "GRM RESULTS - X,Y,Z,(Y-Z):");
    for i in offset + 1..=nt - offset {
        let tay = tf[(i + offset) as usize];
        let tbx = tr[(i - offset) as usize];
        let tg = (tay + tbx - tab - (xy / v1)) / 2.0f32;
        let i = i as usize;
        zz[i] = (tg as f64 * ((v0 * v1) as f64 / (((v1 * v1) - (v0 * v0)) as f64).sqrt())) as f32;
        printf!("%f %f %f %f \n", x[i], y[i], zz[i], y[i] - zz[i]);
    }

    /* Estimate vertical depth using a line tangent to adjacent arcs. */
    printf!("%s \n", "DIP CORRECTED - X,Y,D,Y-D,error:");
    for i in offset + 1..=nt - 1 - offset {
        let i = i as usize;
        d[i] = 0.0f32;
        let mut resid = 99999.0f32;
        let jdx = (2.0f32 * zz[i] / depthres) as i32;

        /* Move x1 along its semicircle; x2 has the same polar angle.
         * Relative elevations suffice, so treat yi as zero. */
        for j in 0..=jdx {
            let x1 = x[i] - zz[i] + j as f32 * depthres;
            let y1 = (((zz[i] * zz[i]) - ((x1 - x[i]) * (x1 - x[i]))) as f64)
                .abs()
                .sqrt() as f32;
            let x2 = ((x1 - x[i]) * zz[i + 1] / zz[i]) + x[i + 1];
            let mut y2 = (((zz[i + 1] * zz[i + 1]) - ((x2 - x[i + 1]) * (x2 - x[i + 1]))) as f64)
                .abs()
                .sqrt() as f32;
            y2 += y[i] - y[i + 1];

            /* Compare the arc tangent slope and the connecting slope. */
            let mc = if y1 != 0.0f32 {
                -(x1 - x[i]) * 1.0f32 / y1
            } else {
                0.00001f32
            };
            let mt = (y2 - y1) / (x2 - x1);
            let newresid = ((mc - mt) as f64).abs() as f32;

            /* Retain the smallest residual and its depth error estimate. */
            if newresid < resid {
                resid = newresid;
                hold = (mc as f64).atan() as f32;
                d[i] = (zz[i] as f64 / (hold as f64).cos()) as f32;
                hold = ((mc + resid) as f64).atan() as f32;
                errd = (d[i] as f64 - zz[i] as f64 / (hold as f64).cos()).abs() as f32;
            }
        }
        printf!("%f %f %f %f (%f) \n", x[i], y[i], d[i], y[i] - d[i], errd);
    }

    exit(CWP_Exit());
}
