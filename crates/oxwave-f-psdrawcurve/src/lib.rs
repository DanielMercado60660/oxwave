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

//! Port of psplot/lib/psDrawCurve.c.
//! Credits: Brian Macy, Phillips Petroleum Co.; adapted from Dave Hale and others' psAxesBox.

use oxwave_f_basic::*;

/// Draw a clipped curve in PostScript.
///
/// C: `void psDrawCurve(float x, float y, float width, float height,
/// float x1Beg, float x1End, float p1Beg, float p1End,
/// float x2Beg, float x2End, float p2Beg, float p2End,
/// float *x1curve, float *x2curve, int ncurve, char *curveColor,
/// float curvewidth, int curvedash, int style);`
///
/// Both coordinate slices must have the same length; their length is `ncurve`.
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn psDrawCurve(
    x: f32,
    y: f32,
    width: f32,
    height: f32,
    x1Beg: f32,
    x1End: f32,
    p1Beg: f32,
    p1End: f32,
    x2Beg: f32,
    x2End: f32,
    p2Beg: f32,
    p2End: f32,
    x1curve: &[f32],
    x2curve: &[f32],
    curveColor: &str,
    curvewidth: f32,
    curvedash: i32,
    style: i32,
) {
    let mut dashes: [&mut [f32]; 12] = [
        &mut [],
        &mut [10., 3.],
        &mut [10., 5.],
        &mut [5.],
        &mut [5., 10.],
        &mut [1., 3.],
        &mut [1., 4.],
        &mut [1., 7.],
        &mut [5., 3., 1., 3.],
        &mut [5., 3., 1., 3., 1., 3.],
        &mut [5., 3., 1., 3., 1., 3., 1., 3.],
        &mut [1., 3., 1., 3., 5., 3., 5., 3.],
    ];
    if curvedash < 0 || curvedash as usize >= dashes.len() {
        eprintln!("No such dash index: {}", curvedash);
        std::process::exit(-1);
    }
    assert_eq!(x1curve.len(), x2curve.len());
    for dash in &mut dashes {
        for value in dash.iter_mut() {
            *value *= curvewidth;
        }
    }
    let (size1, size2) = if style == 0 {
        (width, height)
    } else {
        (height, width)
    };
    gsave();
    translate(x, y);
    if style != 0 {
        rotate(-90.0);
        translate(-height, 0.0);
    }
    newpath();
    setcolor(curveColor);
    setlinewidth(curvewidth);
    setdash(dashes[curvedash as usize], 0.0);
    let xscale = size1 / (x1End + p1End - x1Beg - p1Beg);
    let xbase = -xscale * (x1Beg + p1Beg);
    let yscale = size2 / (x2End + p2End - x2Beg - p2Beg);
    let ybase = -yscale * (x2Beg + p2Beg);
    let mut lpoints = Vec::with_capacity(x1curve.len());
    for i in 0..x1curve.len() {
        lpoints.push((xbase + xscale * x1curve[i], ybase + yscale * x2curve[i]));
    }
    rectclip(0.0, 0.0, size1, size2);
    if lpoints.len() > 1 {
        moveto(lpoints[0].0, lpoints[0].1);
        for p in &lpoints[1..] {
            lineto(p.0, p.1);
        }
        stroke();
    }
    grestore();
}
