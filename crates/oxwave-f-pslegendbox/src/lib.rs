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

//! Port of psplot/lib/psLegendBox.c.
//! Credits: Dave Hale, Ken Larner, Torsten Schoenfelder.
#![allow(non_snake_case)]
use oxwave_axis_scale::scaxis;
use oxwave_f_basic::{fontascender, fontdescender, fontheight, fontwidth};
use oxwave_f_basic::{
    grestore, gsave, justshow, lineto, moveto, newpath, rotate, setcolor, setdash, setfont, show,
    stroke, translate,
};

fn iround(x: f64) -> i32 {
    (if x > 0.0 { x + 0.5 } else { x - 0.5 }) as i32
}

/// Draw a legend as PostScript.
/// C: `void psLegendBox(float x, float y, float width, float height, float x1Beg, float x1End, float p1Beg, float p1End, float d1Num, float f1Num, int n1Tic, int grid1, char *label1, char *labelFont, float labelSize, char *axesColor, char *gridColor, int style);`
#[allow(clippy::too_many_arguments)]
pub fn psLegendBox(
    x: f32,
    y: f32,
    width: f32,
    height: f32,
    x1Beg: f32,
    x1End: f32,
    p1Beg: f32,
    p1End: f32,
    mut d1Num: f32,
    mut f1Num: f32,
    n1Tic: i32,
    grid1: i32,
    label1: &str,
    labelFont: &str,
    labelSize: f32,
    axesColor: &str,
    gridColor: &str,
    style: i32,
) {
    let vert = style == 0 || style == 1;
    let left = style != 1;
    let ch = fontheight(labelFont, labelSize);
    let cw = fontwidth(labelFont, labelSize);
    let ca = fontascender(labelFont, labelSize);
    let cd = fontdescender(labelFont, labelSize);
    if d1Num == 0.0 {
        let mut n = if vert {
            iround(height as f64 / (2.0 * ch as f64))
        } else {
            iround(width as f64 / (3.5 * cw as f64))
        };
        scaxis(x1Beg, x1End, &mut n, &mut d1Num, &mut f1Num);
    }
    gsave();
    translate(x, y);
    newpath();
    setfont(labelFont, labelSize);
    let ticsize = (0.3 * (if vert { cw } else { ch }) as f64) as f32;
    setcolor(axesColor);
    let amin = if x1Beg < x1End { x1Beg } else { x1End };
    let amax = if x1Beg > x1End { x1Beg } else { x1End };
    let azero = (0.0001 * (amax - amin) as f64) as f32;
    let dnum = d1Num;
    let fnum = f1Num;
    let scale = (if vert { height } else { width }) / (x1End + p1End - x1Beg - p1Beg);
    let base = -scale * (x1Beg + p1Beg);
    let ticb = -ticsize;
    let mut numb = (1.2 * ticb as f64 - ca as f64) as f32;
    if vert {
        numb = (1.2 * ticb as f64) as f32;
    }
    let mut labelb = numb - ch;
    let mut ndigits = 0;
    let mut ndig = 0;
    let mut nexpmax: i32 = 0;
    let mut fdexp = (1.001 * dnum.abs() as f64).log10() as f32;
    if fdexp < 0.0 {
        fdexp = (fdexp as f64 - 1.0) as f32;
    }
    let mut anum = fnum;
    while anum <= amax {
        if anum == fnum && anum == 0.0 {
            nexpmax = 0;
        }
        if anum != 0.0 {
            let e = (1.001 * anum.abs() as f64).log10() as f32;
            let nexp = if e > 0.0 { e as i32 } else { e as i32 - 1 };
            if anum == fnum {
                nexpmax = nexp;
            }
            if nexpmax < nexp || nexpmax == 0 {
                nexpmax = nexp;
            }
            ndig = 1 + nexp - fdexp as i32;
        }
        if ndigits < ndig {
            ndigits = ndig;
        }
        anum += dnum;
    }
    let mut nexplot = 0;
    if nexpmax.abs() > 3 {
        nexplot = nexpmax;
    }
    if nexpmax < 0 && ndigits + nexpmax.abs() > 4 {
        nexplot = nexpmax;
    }
    let pnorm = 10.0f64.powf(nexplot as f64) as f32;
    let azeronorm = azero / pnorm;
    fdexp = (1.001 * (dnum / pnorm).abs() as f64).log10() as f32;
    if fdexp < 0.0 {
        fdexp = (fdexp as f64 - 1.0) as f32;
    }
    let mut nplacesmax = 0;
    anum = fnum;
    while anum <= amax {
        if anum >= amin {
            let xa = base + scale * anum;
            let norm = anum / pnorm;
            let mut places = 1;
            let mut nexp = 0;
            if norm < -azeronorm || norm > azeronorm {
                let e = (1.001 * norm.abs() as f64).log10() as f32;
                if e > 0.0 {
                    nexp = e as i32;
                    places = nexp + 1;
                } else {
                    nexp = e as i32 - 1;
                    places = -nexp + 2;
                }
            }
            ndigits = 1 + nexp - fdexp as i32;
            if ndigits > 4 {
                ndigits = 4;
            }
            if nexp >= 0 && ndigits + 1 > places {
                places = ndigits + 1;
            }
            if anum < 0.0 {
                places += 1;
            }
            if nexp < 0 {
                places += ndigits - 1;
            }
            if nplacesmax < places {
                nplacesmax = places;
            }
            let nf = (ndigits - (nexp + 1)).max(0);
            let s = if norm > -azeronorm && norm < azeronorm {
                oxwave_core::sprintf!("%.4g", 0.0f64)
            } else {
                oxwave_core::sprintf!("%.*f", nf, norm)
            };
            if vert {
                if left {
                    moveto(numb, (xa as f64 - ca as f64 / 2.0) as f32);
                    justshow(-1.0, &s);
                } else {
                    moveto(width - numb, (xa as f64 - ca as f64 / 2.0) as f32);
                    justshow(0.0, &s);
                }
            } else {
                moveto(xa, numb);
                justshow(-0.5, &s);
            }
        }
        anum += dnum;
    }
    if nexplot != 0 {
        if !vert {
            moveto(width - numb, 0.0);
            show("x10");
            moveto(
                ((width - numb) as f64 + 1.5 * cw as f64) as f32,
                (0.5 * ca as f64) as f32,
            );
        } else if left {
            moveto((-1.5 * cw as f64) as f32, labelb);
            show("x10");
            moveto(0.0, (labelb as f64 + 0.5 * ca as f64) as f32);
        } else {
            moveto(width, labelb);
            show("x10");
            moveto(
                (width as f64 + 1.5 * cw as f64) as f32,
                (labelb as f64 + 0.5 * ca as f64) as f32,
            );
        }
        show(&nexplot.to_string());
    }
    stroke();
    if vert {
        labelb = (numb as f64 - 0.5 * nplacesmax as f64 * cw as f64) as f32;
    }
    let dtic = dnum / n1Tic as f32;
    let mut atic = fnum - n1Tic as f32 * dtic - dtic;
    while atic <= amax {
        if atic >= amin {
            let xa = base + scale * atic;
            if vert {
                if left {
                    moveto(0.0, xa);
                    lineto(ticb / 2.0, xa);
                } else {
                    moveto(width, xa);
                    lineto(width - ticb / 2.0, xa);
                }
            } else {
                moveto(xa, 0.0);
                lineto(xa, ticb / 2.0);
            }
        }
        atic += dtic;
    }
    stroke();
    if !vert {
        moveto((width as f64 / 2.0) as f32, labelb);
        justshow(-0.5, label1);
    } else {
        if left {
            moveto(labelb + cd, (height as f64 / 2.0) as f32);
        } else {
            moveto(width - labelb + ch, (height as f64 / 2.0) as f32);
        }
        rotate(90.0);
        justshow(-0.5, label1);
        rotate(-90.0);
    }
    setcolor(gridColor);
    if (1..=3).contains(&grid1) {
        anum = fnum;
        while anum <= amax {
            if anum >= amin {
                let xa = base + scale * anum;
                if vert {
                    moveto(0.0, xa);
                    lineto(width, xa);
                } else {
                    moveto(xa, 0.0);
                    lineto(xa, height);
                }
            }
            anum += dnum;
        }
        let dash: &[f32] = match grid1 {
            1 => &[1.0, 5.0],
            2 => &[10.0],
            _ => &[],
        };
        setdash(dash, 0.0);
        stroke();
        setdash(&[], 0.0);
    }
    setcolor(axesColor);
    moveto(0.0, 0.0);
    lineto(width, 0.0);
    lineto(width, height);
    lineto(0.0, height);
    lineto(0.0, 0.0);
    stroke();
    grestore();
}

/// Estimate a legend's bounding box.
/// C: `void psLegendBBox(float x, float y, float width, float height, char *labelFont, float labelSize, int style, int bbox[]);`
pub fn psLegendBBox(
    x: f32,
    y: f32,
    width: f32,
    height: f32,
    labelFont: &str,
    labelSize: f32,
    style: i32,
    bbox: &mut [i32],
) {
    let ch = fontheight(labelFont, labelSize);
    let cw = fontwidth(labelFont, labelSize);
    let ca = fontascender(labelFont, labelSize);
    let cd = fontdescender(labelFont, labelSize);
    if style == 0 || style == 1 {
        bbox[1] = iround(y as f64 - 3.0 * ch as f64) - 1;
        bbox[3] = iround((y + height + ch) as f64) + 1;
        if style == 1 {
            bbox[0] = (x - 2.0) as i32;
            bbox[2] = iround((x + width) as f64 + 4.0 * cw as f64 - cd as f64) + 1;
        } else {
            bbox[0] = iround(x as f64 - 3.36 * cw as f64 + 2.0 * cd as f64 - ch as f64) - 1;
            bbox[2] = (x + width + 3.0) as i32;
        }
    } else {
        bbox[0] = iround(x as f64 - 3.5 * cw as f64 - ch as f64) - 1;
        bbox[1] = iround(y as f64 - 1.36 * ch as f64 - ca as f64 + cd as f64) - 1;
        bbox[2] = iround((x + width) as f64 + 4.0 * cw as f64) + 1;
        bbox[3] = iround((y + height) as f64) + 2;
    }
}
