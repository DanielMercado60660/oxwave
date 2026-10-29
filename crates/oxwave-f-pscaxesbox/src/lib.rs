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

//! Port of psplot/lib/psCAxesBox.c.
//! Authors: Zhiming Li and Dave Hale; modified by Craig Artley.
#![allow(non_snake_case, clippy::too_many_arguments)]
use oxwave_axis_scale::scaxis;
use oxwave_core::sprintf;
use oxwave_f_basic::*;

/// Draw a cube axes box as PostScript.
/// C prototype: `void psCubeAxesBox(float x, float y, float size1, float size2, float size3, float angle, float x1Beg, float x1End, float p1Beg, float p1End, float d1Num, float f1Num, int n1Tic, int grid1, char *label1, float x2Beg, float x2End, float p2Beg, float p2End, float d2Num, float f2Num, int n2Tic, int grid2, char *label2, float x3Beg, float x3End, float p3Beg, float p3End, float d3Num, float f3Num, int n3Tic, int grid3, char *label3, char *labelFont, float labelSize, char *title, char *titleFont, float titleSize, char *titleColor, char *axesColor, char *gridColor);`
pub fn psCubeAxesBox(
    x: f32,
    y: f32,
    size1: f32,
    size2: f32,
    size3: f32,
    angle: f32,
    x1Beg: f32,
    x1End: f32,
    p1Beg: f32,
    p1End: f32,
    mut d1Num: f32,
    mut f1Num: f32,
    n1Tic: i32,
    grid1: i32,
    label1: &str,
    x2Beg: f32,
    x2End: f32,
    p2Beg: f32,
    p2End: f32,
    mut d2Num: f32,
    mut f2Num: f32,
    n2Tic: i32,
    grid2: i32,
    label2: &str,
    x3Beg: f32,
    x3End: f32,
    p3Beg: f32,
    p3End: f32,
    mut d3Num: f32,
    mut f3Num: f32,
    n3Tic: i32,
    grid3: i32,
    label3: &str,
    labelFont: &str,
    labelSize: f32,
    title: &str,
    titleFont: &str,
    titleSize: f32,
    titleColor: &str,
    axesColor: &str,
    gridColor: &str,
) {
    let labelCH = fontheight(labelFont, labelSize);
    let labelCW = fontwidth(labelFont, labelSize);
    let labelCD = fontdescender(labelFont, labelSize);
    let titleCH = fontheight(titleFont, titleSize);
    let titleCA = fontascender(titleFont, titleSize);
    if d1Num == 0.0 {
        let mut n = (size1 / (4.0 * labelCW)) as i32;
        scaxis(x1Beg, x1End, &mut n, &mut d1Num, &mut f1Num);
    }
    if d2Num == 0.0 {
        let mut n = (size2 / (4.0 * labelCW)) as i32;
        scaxis(x2Beg, x2End, &mut n, &mut d2Num, &mut f2Num);
    }
    if d3Num == 0.0 {
        let mut n = (size3 / (4.0 * labelCW)) as i32;
        scaxis(x3Beg, x3End, &mut n, &mut d3Num, &mut f3Num);
    }
    gsave();
    translate(x, y);
    setfont(labelFont, labelSize);
    let ticsize = (0.3 * labelCH as f64) as f32;
    let s = (angle as f64).sin();
    let c = (angle as f64).cos();
    let height = (size1 as f64 + size3 as f64 * s) as f32;
    let width = (size2 as f64 + size3 as f64 * c) as f32;
    let yyy = (angle as f64 * 180.0 / 3.141592654) as f32;
    let axes = [
        (
            size1, x1Beg, x1End, p1Beg, p1End, d1Num, f1Num, n1Tic, grid1, label1,
        ),
        (
            size2, x2Beg, x2End, p2Beg, p2End, d2Num, f2Num, n2Tic, grid2, label2,
        ),
        (
            size3, x3Beg, x3End, p3Beg, p3End, d3Num, f3Num, n3Tic, grid3, label3,
        ),
    ];
    for (i, &(size, beg, end, pbeg, pend, dnum, fnum, ntic, grid, label)) in axes.iter().enumerate()
    {
        gsave();
        setcolor(axesColor);
        match i {
            0 => {
                translate(0.0, size1);
                rotate(-90.0);
            }
            1 => {
                translate(width - size2, height);
                rotate(-90.0);
            }
            _ => {
                translate(0.0, size1);
                rotate(yyy);
            }
        }
        let amin = if beg < end { beg } else { end };
        let amax = if beg > end { beg } else { end };
        let azero = (0.0001 * (amax - amin) as f64) as f32;
        let scale = size / (end + pend - beg - pbeg);
        let base = -scale * (beg + pbeg);
        let ticb = if i == 2 { ticsize } else { -ticsize };
        let numb = if i == 2 {
            ticb - labelCD
        } else {
            ticb + labelCD
        };
        let labelb = match i {
            0 => (numb as f64 - 2.0 * labelCW as f64) as f32,
            1 => numb - labelCW,
            _ => numb + labelCW,
        };
        let mut anum = fnum;
        while anum <= amax {
            if anum >= amin {
                let a = base + scale * anum;
                match i {
                    0 => {
                        moveto(a, 0.0);
                        lineto(a, ticb);
                    }
                    1 => {
                        moveto(0.0, a);
                        lineto(ticb, a);
                    }
                    _ => {
                        moveto(a, 0.0);
                        lineto(
                            (a as f64 - ticb as f64 * c) as f32,
                            (ticb as f64 * s) as f32,
                        );
                    }
                }
                let text = sprintf!(
                    "%.4g",
                    if anum > -azero && anum < azero {
                        0.0
                    } else {
                        anum
                    }
                );
                match i {
                    0 => {
                        moveto((a as f64 + 0.3 * labelCH as f64) as f32, numb);
                        rotate(90.0);
                        justshow(-1.0, &text);
                        rotate(-90.0);
                    }
                    1 => {
                        moveto(numb, a);
                        rotate(90.0);
                        justshow(-0.5, &text);
                        rotate(-90.0);
                    }
                    _ => {
                        moveto(
                            (a as f64 - numb as f64 * c + labelCD as f64 * s) as f32,
                            (numb as f64 * s + labelCD as f64 * c) as f32,
                        );
                        rotate(-yyy);
                        justshow(-1.0, &text);
                        rotate(yyy);
                    }
                }
            }
            anum += dnum;
        }
        stroke();
        let dtic = dnum / (ntic as f32);
        let mut atic = fnum - (ntic as f32) * dtic - dtic;
        while atic <= amax {
            if atic >= amin {
                let a = base + scale * atic;
                match i {
                    0 => {
                        moveto(a, 0.0);
                        lineto(a, ticb / 2.0);
                    }
                    1 => {
                        moveto(0.0, a);
                        lineto(ticb / 2.0, a);
                    }
                    _ => {
                        moveto(a, 0.0);
                        lineto(
                            (a as f64 - ticb as f64 / 2.0 * c) as f32,
                            (ticb as f64 / 2.0 * s) as f32,
                        );
                    }
                }
            }
            atic += dtic;
        }
        stroke();
        match i {
            0 => {
                moveto((size1 as f64 / 2.0) as f32, labelb);
                rotate(180.0);
                justshow(-0.5, label);
                rotate(-180.0);
            }
            1 => {
                moveto(labelb, (size2 as f64 / 2.0) as f32);
                rotate(90.0);
                justshow(-0.5, label);
                rotate(-90.0);
            }
            _ => {
                moveto((size3 as f64 / 2.0) as f32, labelb);
                justshow(-0.5, label);
            }
        }
        setcolor(gridColor);
        if (1..=3).contains(&grid) {
            let mut a = fnum;
            while a <= amax {
                if a >= amin {
                    let mut xa = base + scale * a;
                    match i {
                        0 => {
                            moveto(xa, 0.0);
                            lineto(xa, size2);
                            xa = (xa as f64 - size3 as f64 * s) as f32;
                            let ya = (size2 as f64 + size3 as f64 * c) as f32;
                            lineto(xa, ya);
                        }
                        1 => {
                            let mut ya = xa;
                            moveto(0.0, ya);
                            ya = (ya as f64 - size3 as f64 * c) as f32;
                            xa = (size3 as f64 * s) as f32;
                            lineto(xa, ya);
                            xa += size1;
                            lineto(xa, ya);
                        }
                        _ => {
                            moveto(xa, 0.0);
                            xa = (xa as f64 + size2 as f64 * c) as f32;
                            let mut ya = (-size2 as f64 * s) as f32;
                            lineto(xa, ya);
                            xa = (xa as f64 - size1 as f64 * s) as f32;
                            ya = (ya as f64 - size1 as f64 * c) as f32;
                            lineto(xa, ya);
                        }
                    }
                }
                a += dnum;
            }
            let dash: &[f32] = match grid {
                1 => &[1.0, 5.0],
                2 => &[10.0],
                _ => &[],
            };
            setdash(dash, 0.0);
            stroke();
            setdash(&[], 0.0);
        }
        grestore();
    }
    setcolor(titleColor);
    setfont(titleFont, titleSize);
    moveto(
        ((size2 as f64 + size3 as f64 * c) / 2.0) as f32,
        (-0.5 * titleCH as f64 - titleCA as f64) as f32,
    );
    justshow(-0.5, title);
    setcolor(axesColor);
    moveto(0.0, 0.0);
    rlineto(size2, 0.0);
    rlineto(0.0, size1);
    rlineto(-size2, 0.0);
    rlineto(0.0, -size1);
    moveto(0.0, size1);
    rlineto((size3 as f64 * c) as f32, (size3 as f64 * s) as f32);
    rlineto(size2, 0.0);
    rlineto((-size3 as f64 * c) as f32, (-size3 as f64 * s) as f32);
    moveto(size2, 0.0);
    rlineto((size3 as f64 * c) as f32, (size3 as f64 * s) as f32);
    rlineto(0.0, size1);
    stroke();
    grestore();
}
