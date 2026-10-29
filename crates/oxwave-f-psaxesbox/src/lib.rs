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

//! Port of psplot/lib/psAxesBox.c.
//! Author: Dave Hale, Colorado School of Mines, 06/27/89.
//! Modified: Ken Larner, Colorado School of Mines, 08/30/90;
//! Dave Hale, Advance Geophysical, 10/18/92 (color parameters);
//! Morten Wendell Pedersen, Aarhus University, 23/3-97 (line widths).

#![allow(non_snake_case, clippy::too_many_arguments)]

use oxwave_axis_scale::scaxis;
use oxwave_core::sprintf;
use oxwave_f_basic::{
    fontascender, fontdescender, fontheight, fontwidth, grestore, gsave, justshow, lineto, moveto,
    newpath, rotate, setcolor, setdash, setfont, setlinewidth, show, stroke, translate,
};

/// C: `void psAxesBBox(float x, float y, float width, float height, char *labelFont, float labelSize, char *titleFont, float titleSize, int style, int bbox[]);`
pub fn psAxesBBox(
    x: f32,
    y: f32,
    width: f32,
    height: f32,
    labelFont: &str,
    labelSize: f32,
    titleFont: &str,
    titleSize: f32,
    style: i32,
    bbox: &mut [i32],
) {
    let labelCH = fontheight(labelFont, labelSize);
    let labelCW = fontwidth(labelFont, labelSize);
    let titleCH = fontheight(titleFont, titleSize);
    bbox[0] = (x as f64 - 3.5 * labelCW as f64 - labelCH as f64) as i32;
    if style == 0 {
        bbox[1] = (y as f64 - 3.0 * labelCH as f64) as i32;
        bbox[2] = ((x + width) as f64 + 3.0 * labelCW as f64) as i32;
        bbox[3] = ((y + height) as f64 + 2.0 * titleCH as f64 + labelCH as f64) as i32;
    } else {
        bbox[1] = (y as f64 - 2.0 * titleCH as f64 - labelCH as f64) as i32;
        bbox[2] = ((x + width) as f64 + 3.0 * labelCW as f64) as i32;
        bbox[3] = ((y + height) as f64 + 3.0 * labelCH as f64) as i32;
    }
}

/// Draw axes, numbered ticks, grids, labels, and a title to standard output.
/// Grid codes are 0 (none), 1 (dots), 2 (dashes), and 3 (solid).
/// Style 0 is normal; every other value uses seismic orientation.
/// As in C, tick increments must advance and loops must terminate.
///
/// C prototype:
/// ```c
/// void psAxesBox(
///     float x, float y, float width, float height,
///     float x1Beg, float x1End, float p1Beg, float p1End,
///     float d1Num, float f1Num, int n1Tic, int grid1, char *label1,
///     float x2Beg, float x2End, float p2Beg, float p2End,
///     float d2Num, float f2Num, int n2Tic, int grid2, char *label2,
///     char *labelFont, float labelSize,
///     char *title, char *titleFont, float titleSize,
///     char *titleColor, char *axesColor, char *gridColor,
///     float ticwidth, float axeswidth, float gridwidth, int style);
/// ```
pub fn psAxesBox(
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
    x2Beg: f32,
    x2End: f32,
    p2Beg: f32,
    p2End: f32,
    mut d2Num: f32,
    mut f2Num: f32,
    n2Tic: i32,
    grid2: i32,
    label2: &str,
    labelFont: &str,
    labelSize: f32,
    title: &str,
    titleFont: &str,
    titleSize: f32,
    titleColor: &str,
    axesColor: &str,
    gridColor: &str,
    ticwidth: f32,
    axeswidth: f32,
    gridwidth: f32,
    style: i32,
) {
    // Determine font dimensions.
    let labelCH = fontheight(labelFont, labelSize);
    let labelCW = fontwidth(labelFont, labelSize);
    let labelCA = fontascender(labelFont, labelSize);
    let labelCD = fontdescender(labelFont, labelSize);
    let titleCH = fontheight(titleFont, titleSize);
    let titleCA = fontascender(titleFont, titleSize);
    let titleCD = fontdescender(titleFont, titleSize);

    let (size1, size2) = if style == 0 {
        (width, height)
    } else {
        (height, width)
    };
    if d1Num == 0.0 {
        let mut n1num = (size1 / (4.0f32 * labelCW)) as i32;
        scaxis(x1Beg, x1End, &mut n1num, &mut d1Num, &mut f1Num);
    }
    if d2Num == 0.0 {
        let mut n2num = (size2 / (4.0f32 * labelCW)) as i32;
        scaxis(x2Beg, x2End, &mut n2num, &mut d2Num, &mut f2Num);
    }
    gsave();
    translate(x, y);
    if style != 0 {
        rotate(-90.0);
        translate(-height, 0.0);
    }
    newpath();
    setfont(labelFont, labelSize);
    let ticsize = (0.3 * labelCH as f64) as f32;

    // Draw axis 1.
    setcolor(axesColor);
    setlinewidth(ticwidth);
    let mut amin = if x1Beg < x1End { x1Beg } else { x1End };
    let mut amax = if x1Beg > x1End { x1Beg } else { x1End };
    let mut azero = (0.0001 * (amax - amin) as f64) as f32;
    let mut dnum = d1Num;
    let mut fnum = f1Num;
    let mut ntic = n1Tic;
    let mut scale = size1 / (x1End + p1End - x1Beg - p1Beg);
    let mut base = -scale * (x1Beg + p1Beg);
    let mut ticb = -ticsize;
    let mut numb = (1.2 * ticb as f64 - labelCA as f64) as f32;
    if style != 0 {
        numb = (1.2 * ticb as f64) as f32;
    }
    let mut labelb = numb - labelCH;
    let mut nexpmax = 0i32;

    // Determine axis 1 exponent for scientific notation.
    let mut ndigits = 0i32;
    let mut ndig = 0i32;
    let mut fdexp = (1.001 * c_abs(dnum) as f64).log10() as f32;
    if fdexp < 0.0 {
        fdexp = (fdexp as f64 - 1.0) as f32;
    }
    let mut anum = fnum;
    while anum <= amax {
        if anum == fnum && anum == 0.0 {
            nexpmax = 0;
        }
        if anum != 0.0 {
            let fnexp = (1.001 * c_abs(anum) as f64).log10() as f32;
            let nexp = if fnexp > 0.0 {
                fnexp as i32
            } else {
                (fnexp as i32).wrapping_sub(1)
            };
            if anum == fnum {
                nexpmax = nexp;
            }
            if nexpmax < nexp || nexpmax == 0 {
                nexpmax = nexp;
            }
            ndig = 1i32.wrapping_add(nexp).wrapping_sub(fdexp as i32);
        }
        if ndigits < ndig {
            ndigits = ndig;
        }
        anum += dnum;
    }
    let mut nexplot1 = 0;
    if nexpmax.wrapping_abs() > 3 {
        nexplot1 = nexpmax;
    }
    if nexpmax < 0 && ndigits.wrapping_add(nexpmax.wrapping_abs()) > 4 {
        nexplot1 = nexpmax;
    }

    // Loop for axis 1 numbering.
    let mut nplacesmax = 0i32;
    let pnorm = 10.0f64.powf(nexplot1 as f64) as f32;
    let azeronorm = azero / pnorm;
    fdexp = (1.001 * c_abs(dnum / pnorm) as f64).log10() as f32;
    if fdexp < 0.0 {
        fdexp = (fdexp as f64 - 1.0) as f32;
    }
    anum = fnum;
    while anum <= amax {
        if anum < amin {
            anum += dnum;
            continue;
        }
        let xa = base + scale * anum;
        let xas = (xa as f64 + 0.25 * labelCH as f64) as f32;
        moveto(xa, 0.0);
        lineto(xa, ticb);
        let anumnorm = anum / pnorm;

        // Find number of places; numbers are limited to four significant digits.
        let mut nplaces = 1i32;
        let mut nexp = 0i32;
        if anumnorm < -azeronorm || anumnorm > azeronorm {
            let fnexp = (1.001 * c_abs(anumnorm) as f64).log10() as f32;
            if fnexp > 0.0 {
                nexp = fnexp as i32;
                nplaces = nexp.wrapping_add(2);
            } else {
                nexp = (fnexp as i32).wrapping_sub(1);
                nplaces = nexp.wrapping_neg().wrapping_add(2);
            }
        }
        ndigits = 1i32.wrapping_add(nexp).wrapping_sub(fdexp as i32);
        if ndigits > 4 {
            ndigits = 4;
        }
        if nexp >= 0 && ndigits.wrapping_add(1) > nplaces {
            nplaces = ndigits.wrapping_add(1);
        }
        if anum < 0.0 {
            nplaces = nplaces.wrapping_add(1);
        }
        if nexp < 0 {
            nplaces = nplaces.wrapping_add(ndigits).wrapping_sub(1);
        }
        if nplacesmax < nplaces {
            nplacesmax = nplaces;
        }
        let mut nformat = ndigits.wrapping_sub(nexp.wrapping_add(1));
        if nformat < 0 {
            nformat = 0;
        }
        let text = if anumnorm > -azeronorm && anumnorm < azeronorm {
            sprintf!("%.4g", 0.0f64)
        } else {
            let sformat = sprintf!("%%.%df", nformat);
            sprintf!(&sformat, anumnorm)
        };
        if style != 0 {
            moveto(xas, numb);
            rotate(90.0);
            justshow(-1.0, &text);
            rotate(-90.0);
        } else {
            moveto(xa, numb);
            justshow(-0.5, &text);
        }
        anum += dnum;
    }
    // Draw exponential multiplier for axis 1.
    if nexplot1 != 0 {
        if style == 0 {
            moveto(size1, labelb);
            show("x10");
            moveto(
                (size1 as f64 + 1.5 * labelCW as f64) as f32,
                (labelb as f64 + 0.5 * labelCA as f64) as f32,
            );
            show(&sprintf!("%d", nexplot1));
        } else {
            moveto(
                (size1 as f64 + 1.4 * labelCH as f64) as f32,
                (-2.0 * labelCW as f64) as f32,
            );
            rotate(90.0);
            show("x10");
            rotate(-90.0);
            moveto(
                (size1 as f64 + 1.05 * labelCH as f64) as f32,
                (-0.5 * labelCW as f64) as f32,
            );
            let text = sprintf!("%d", nexplot1);
            rotate(90.0);
            show(&text);
            rotate(-90.0);
        }
    }
    stroke();

    // Draw axis 1 tick marks.
    if style != 0 {
        labelb = (numb as f64 - 0.5 * nplacesmax as f64 * labelCW as f64) as f32;
    }
    let dtic = dnum / ntic as f32;
    let mut atic = fnum - ntic as f32 * dtic - dtic;
    while atic <= amax {
        if atic < amin {
            atic += dtic;
            continue;
        }
        let xa = base + scale * atic;
        moveto(xa, 0.0);
        lineto(xa, ticb / 2.0f32);
        atic += dtic;
    }
    stroke();
    moveto((size1 as f64 / 2.0) as f32, labelb);
    if style != 0 {
        rotate(180.0);
        justshow(-0.5, label1);
        rotate(-180.0);
    } else {
        justshow(-0.5, label1);
    }

    // Draw axis 1 grid.
    setcolor(gridColor);
    setlinewidth(gridwidth);
    let mut dash = [0.0f32; 2];
    let ndash = if grid1 == 3 {
        0
    } else if grid1 == 2 {
        dash[0] = 10.0;
        1
    } else if grid1 == 1 {
        dash[0] = 1.0;
        dash[1] = 5.0;
        2
    } else {
        0
    };
    if grid1 == 3 || grid1 == 2 || grid1 == 1 {
        anum = fnum;
        while anum <= amax {
            if anum < amin {
                anum += dnum;
                continue;
            }
            let xa = base + scale * anum;
            moveto(xa, 0.0);
            lineto(xa, size2);
            anum += dnum;
        }
        setdash(&dash[..ndash], 0.0);
        stroke();
        setdash(&[], 0.0);
    }

    // Draw axis 2.
    setcolor(axesColor);
    setlinewidth(ticwidth);
    amin = if x2Beg < x2End { x2Beg } else { x2End };
    amax = if x2Beg > x2End { x2Beg } else { x2End };
    azero = (0.0001 * (amax - amin) as f64) as f32;
    dnum = d2Num;
    fnum = f2Num;
    ntic = n2Tic;
    scale = size2 / (x2End + p2End - x2Beg - p2Beg);
    base = -scale * (x2Beg + p2Beg);
    ticb = -ticsize;
    numb = ticb + labelCD;
    labelb = numb - labelCH;

    // Determine axis 2 exponent for scientific notation.
    let mut ndigits = 0i32;
    let mut ndig = 0i32;
    let mut fdexp = (1.001 * c_abs(dnum) as f64).log10() as f32;
    if fdexp < 0.0 {
        fdexp = (fdexp as f64 - 1.0) as f32;
    }
    let mut anum = fnum;
    while anum <= amax {
        if anum == fnum && anum == 0.0 {
            nexpmax = 0;
        }
        if anum != 0.0 {
            let fnexp = (1.001 * c_abs(anum) as f64).log10() as f32;
            let nexp = if fnexp > 0.0 {
                fnexp as i32
            } else {
                (fnexp as i32).wrapping_sub(1)
            };
            if anum == fnum {
                nexpmax = nexp;
            }
            if nexpmax < nexp || nexpmax == 0 {
                nexpmax = nexp;
            }
            ndig = 1i32.wrapping_add(nexp).wrapping_sub(fdexp as i32);
        }
        if ndigits < ndig {
            ndigits = ndig;
        }
        anum += dnum;
    }
    let mut nexplot2 = 0;
    if nexpmax.wrapping_abs() > 3 {
        nexplot2 = nexpmax;
    }
    if nexpmax < 0 && ndigits.wrapping_add(nexpmax.wrapping_abs()) > 4 {
        nexplot2 = nexpmax;
    }

    // Loop for axis 2 numbering.
    let mut nplacesmax = 0i32;
    let pnorm = 10.0f64.powf(nexplot2 as f64) as f32;
    let azeronorm = azero / pnorm;
    fdexp = (1.001 * c_abs(dnum / pnorm) as f64).log10() as f32;
    if fdexp < 0.0 {
        fdexp = (fdexp as f64 - 1.0) as f32;
    }
    anum = fnum;
    while anum <= amax {
        if anum < amin {
            anum += dnum;
            continue;
        }
        let ya = base + scale * anum;
        moveto(0.0, ya);
        lineto(ticb, ya);
        let anumnorm = anum / pnorm;
        let yas = (ya as f64 - 0.25 * labelCH as f64) as f32;

        // Find number of places; numbers are limited to four significant digits.
        let mut nplaces = 1i32;
        let mut nexp = 1i32;
        if anumnorm < -azeronorm || anumnorm > azeronorm {
            let fnexp = (1.001 * c_abs(anumnorm) as f64).log10() as f32;
            if fnexp > 0.0 {
                nexp = fnexp as i32;
                nplaces = nexp.wrapping_add(2);
            } else {
                nexp = (fnexp as i32).wrapping_sub(1);
                nplaces = nexp.wrapping_neg().wrapping_add(2);
            }
        }
        ndigits = 1i32.wrapping_add(nexp).wrapping_sub(fdexp as i32);
        if ndigits > 4 {
            ndigits = 4;
        }
        if nexp >= 0 && ndigits.wrapping_add(1) > nplaces {
            nplaces = ndigits.wrapping_add(1);
        }
        if anum < 0.0 {
            nplaces = nplaces.wrapping_add(1);
        }
        if nexp < 0 {
            nplaces = nplaces.wrapping_add(ndigits).wrapping_sub(1);
        }
        if nplacesmax < nplaces {
            nplacesmax = nplaces;
        }
        let mut nformat = ndigits.wrapping_sub(nexp.wrapping_add(1));
        if nformat < 0 {
            nformat = 0;
        }
        let text = if anumnorm > -azeronorm && anumnorm < azeronorm {
            sprintf!("%.4g", 0.0f64)
        } else {
            let sformat = sprintf!("%%.%df", nformat);
            sprintf!(&sformat, anumnorm)
        };
        if style != 0 {
            moveto(numb, ya);
            rotate(90.0);
            justshow(-0.5, &text);
            rotate(-90.0);
        } else {
            moveto(numb, yas);
            justshow(-1.0, &text);
        }
        anum += dnum;
    }
    // Draw exponential multiplier for axis 2.
    if nexplot2 != 0 {
        if style == 0 {
            let mut numb2 = (size2 as f64 + 0.5 * labelCH as f64) as f32;
            moveto(labelb, numb2);
            show("x10");
            let xb = (labelb as f64 + 1.5 * labelCW as f64) as f32;
            numb2 = (numb2 as f64 + 0.5 * labelCA as f64) as f32;
            moveto(xb, numb2);
            show(&sprintf!("%d", nexplot2));
        } else {
            moveto(labelb, (size2 as f64 - 0.5 * labelCW as f64) as f32);
            rotate(90.0);
            show("x10");
            rotate(-90.0);
            moveto(
                (labelb as f64 - 0.35 * labelCH as f64) as f32,
                (size2 as f64 + 1.0 * labelCW as f64) as f32,
            );
            let text = sprintf!("%d", nexplot2);
            rotate(90.0);
            show(&text);
            rotate(-90.0);
        }
    }
    stroke();

    // Draw axis 2 tick marks.
    if style == 0 {
        labelb = (numb as f64 - 0.5 * nplacesmax as f64 * labelCW as f64) as f32;
    }
    let dtic = dnum / ntic as f32;
    let mut atic = fnum - ntic as f32 * dtic - dtic;
    while atic <= amax {
        if atic < amin {
            atic += dtic;
            continue;
        }
        let ya = base + scale * atic;
        moveto(0.0, ya);
        lineto(ticb / 2.0f32, ya);
        atic += dtic;
    }
    stroke();
    moveto(labelb, (size2 as f64 / 2.0) as f32);
    rotate(90.0);
    justshow(-0.5, label2);
    rotate(-90.0);

    // Draw axis 2 grid.
    setcolor(gridColor);
    setlinewidth(gridwidth);
    let ndash = if grid2 == 3 {
        0
    } else if grid2 == 2 {
        dash[0] = 10.0;
        1
    } else if grid2 == 1 {
        dash[0] = 1.0;
        dash[1] = 5.0;
        2
    } else {
        0
    };
    if grid2 == 3 || grid2 == 2 || grid2 == 1 {
        anum = fnum;
        while anum <= amax {
            if anum < amin {
                anum += dnum;
                continue;
            }
            let ya = base + scale * anum;
            moveto(0.0, ya);
            lineto(size1, ya);
            anum += dnum;
        }
        setdash(&dash[..ndash], 0.0);
        stroke();
        setdash(&[], 0.0);
    }

    // Draw axes box.
    setcolor(axesColor);
    setlinewidth(axeswidth);
    moveto(0.0, 0.0);
    lineto(size1, 0.0);
    lineto(size1, size2);
    lineto(0.0, size2);
    lineto(0.0, 0.0);
    stroke();

    // Draw title.
    setcolor(titleColor);
    setfont(titleFont, titleSize);
    if style == 0 {
        let mut numb2 = (size2 as f64 + 0.5 * titleCH as f64 - titleCD as f64) as f32;
        if nexplot2 != 0 {
            numb2 = (numb2 as f64 + 0.5 * (labelCH + labelCA) as f64) as f32;
        }
        moveto((size1 as f64 / 2.0) as f32, numb2);
        justshow(-0.5, title);
    } else {
        let mut numb2 = (size1 as f64 + 0.5 * titleCH as f64 + titleCA as f64) as f32;
        if nexplot1 != 0 {
            numb2 += labelCH;
        }
        moveto(numb2, (size2 as f64 / 2.0) as f32);
        rotate(90.0);
        justshow(-0.5, title);
        rotate(-90.0);
    }
    grestore();
}

// C's ABS macro retains negative zero (unlike f32::abs).
fn c_abs(x: f32) -> f32 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}
