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

//! Draw contours via PostScript. Port of `psplot/lib/psContour.c`.
//! Author: Dave Hale, Colorado School of Mines, 06/28/89.
//! Contour labeling added by: Zhenyue Liu, August 1993.
#![allow(non_snake_case)]
use oxwave_core::sprintf;
use oxwave_f_basic::*;
fn bit(z: f32, b: u32) -> bool {
    z.to_bits() & b != 0
}
fn set(z: &mut f32, b: u32) {
    *z = f32::from_bits(z.to_bits() | b)
}
fn clear(z: &mut f32, b: u32) {
    *z = f32::from_bits(z.to_bits() & !b)
}
fn min(a: f32, b: f32) -> f32 {
    if a < b {
        a
    } else {
        b
    }
}
fn max(a: f32, b: f32) -> f32 {
    if a > b {
        a
    } else {
        b
    }
}
fn delta(c: f32, a: f32, b: f32) -> f32 {
    if a != b {
        (c - a) / (b - a)
    } else {
        1.0
    }
}
fn interp(d: f32, a: f32, b: f32) -> f32 {
    // C promotes (1.0-d)*a to double, but evaluates d*b in float first.
    ((1.0 - d as f64) * a as f64 + (d * b) as f64) as f32
}
// PATHMAX = PATHLIMIT/2 in psplot.h: 1500/2.
struct Draw(i32);
impl Draw {
    fn mov(&mut self, x: f32, y: f32) {
        self.0 += 1;
        if self.0 > 750 {
            stroke();
            self.0 = 0
        }
        moveto(x, y)
    }
    fn line(&mut self, x: f32, y: f32) {
        lineto(x, y);
        self.0 += 1;
        if self.0 > 750 {
            stroke();
            moveto(x, y);
            self.0 = 0
        }
    }
}
fn connect(
    c: f32,
    x: &[f32],
    y: &[f32],
    z: &mut [f32],
    cell: &mut usize,
    draw: &mut Draw,
) -> Option<(f32, f32, bool)> {
    let nx = x.len();
    let ny = y.len();
    let mut k = *cell;
    let mut ix = k % nx;
    let mut iy = k / nx;
    let (xd, yd, more);
    if bit(z[k + nx], 1) {
        k += nx;
        iy += 1;
        let d = delta(c, z[k], z[k + 1]);
        xd = interp(d, x[ix], x[ix + 1]);
        yd = y[iy];
        clear(&mut z[k], 1);
        *cell = k;
        more = iy < ny - 1
    } else if bit(z[k + 1], 2) {
        k += 1;
        ix += 1;
        let d = delta(c, z[k], z[k + nx]);
        xd = x[ix];
        yd = interp(d, y[iy], y[iy + 1]);
        clear(&mut z[k], 2);
        *cell = k;
        more = ix < nx - 1
    } else if bit(z[k], 1) {
        let d = delta(c, z[k], z[k + 1]);
        xd = interp(d, x[ix], x[ix + 1]);
        yd = y[iy];
        clear(&mut z[k], 1);
        *cell = k.wrapping_sub(nx);
        more = iy > 0
    } else if bit(z[k], 2) {
        let d = delta(c, z[k], z[k + nx]);
        xd = x[ix];
        yd = interp(d, y[iy], y[iy + 1]);
        clear(&mut z[k], 2);
        *cell = k.wrapping_sub(1);
        more = ix > 0
    } else {
        return None;
    }
    draw.line(xd, yd);
    Some((xd, yd, more))
}
fn wcell(x: &[f32], y: &[f32], w: &mut [f32], cell: usize, xc: f32, yc: f32, xw: f32, yw: f32) {
    let nx = x.len();
    let (mut xl, mut xh, mut yl, mut yh) = (cell % nx, cell % nx, cell / nx, cell / nx);
    while x[xl] as f64 >= xc as f64 - 2.0 * xw as f64
        && x[xl] as f64 <= xc as f64 + 2.0 * xw as f64
        && xl > 0
    {
        xl -= 1
    }
    while x[xh] as f64 >= xc as f64 - 2.0 * xw as f64
        && x[xh] as f64 <= xc as f64 + 2.0 * xw as f64
        && xh < nx - 1
    {
        xh += 1
    }
    while y[yl] as f64 >= yc as f64 - 1.5 * yw as f64
        && y[yl] as f64 <= yc as f64 + 1.5 * yw as f64
        && yl > 0
    {
        yl -= 1
    }
    while y[yh] as f64 >= yc as f64 - 1.5 * yw as f64
        && y[yh] as f64 <= yc as f64 + 1.5 * yw as f64
        && yh < y.len() - 1
    {
        yh += 1
    }
    for iy in yl..=yh {
        for ix in xl..=xh {
            w[iy * nx + ix] += 1.0
        }
    }
}
struct Label<'a> {
    xw: f32,
    yw: f32,
    size: f32,
    text: &'a str,
    font: &'a str,
    color: &'a str,
    bounds: (f32, f32, f32, f32),
}
fn follow(
    c: f32,
    x: &[f32],
    y: &[f32],
    z: &mut [f32],
    w: &mut [f32],
    mut cell: usize,
    non: &mut i32,
    draw: &mut Draw,
    label: &Label<'_>,
) {
    let (xmin, xmax, ymin, ymax) = label.bounds;
    let (mut xdmin, mut xdmax, mut ydmin, mut ydmax) = (xmax, xmin, ymax, ymin);
    let (mut xc, mut yc, mut cells, mut id) = (0.0, 0.0, 0, true);
    while let Some((xd, yd, more)) = connect(c, x, y, z, &mut cell, draw) {
        if !more {
            break;
        }
        *non -= 1;
        if w[cell] < 0.5 && id {
            xc = xd;
            yc = yd;
            cells = cell;
            id = false
        }
        xdmin = min(xdmin, xd);
        xdmax = max(xdmax, xd);
        ydmin = min(ydmin, yd);
        ydmax = max(ydmax, yd)
    }
    if label.size > 1.0 && !id && xdmax + ydmax - xdmin - ydmin > label.xw + label.yw {
        wcell(x, y, w, cells, xc, yc, label.xw, label.yw);
        let a = xc - label.xw / 2.0;
        let b = yc - label.yw / 2.0;
        stroke();
        gsave();
        newpath();
        setcolor("white");
        rectfill(a, b, label.xw, label.yw);
        setcolor(label.color);
        setfont(label.font, label.size);
        moveto(a + label.xw, b);
        rotate(90.0);
        show(label.text);
        rotate(-90.0);
        stroke();
        grestore()
    }
}
/// Draw a contour, modifying the flag bits of `z` and the label mask `w`.
///
/// C: `void psContour(float c, int nx, float x[], int ny, float y[], float z[], float lcs, char *lcf, char *lcc, float *w, int nplaces);`
///
/// Coordinates may be nonuniform. Values use row-major order with x varying
/// fastest. The two least significant bits of each z value mark intersections;
/// the mask w accumulates regions excluded from contour labeling. PostScript is
/// written to standard output through `oxwave_f_basic`.
///
/// Both coordinates must be nonempty; grids must contain `x.len()*y.len()` values.
/// A grid with one row or column must have no contour intersections: the C
/// implementation accesses outside the grid when tracing such an intersection.
pub fn psContour(
    c: f32,
    x: &[f32],
    y: &[f32],
    z: &mut [f32],
    lcs: f32,
    lcf: &str,
    lcc: &str,
    w: &mut [f32],
    nplaces: i32,
) {
    let nx = x.len();
    let ny = y.len();
    assert!(nx > 0 && ny > 0 && z.len() >= nx * ny && w.len() >= nx * ny);
    let (xmin, xmax, ymin, ymax) = (
        min(x[0], x[nx - 1]),
        max(x[0], x[nx - 1]),
        min(y[0], y[ny - 1]),
        max(y[0], y[ny - 1]),
    );
    let text = sprintf!("%.*g", nplaces, c);
    let yw = (lcs as f64 * 0.55 * text.len() as f64) as f32;
    let xw = (lcs as f64 * 0.8) as f32;
    for iy in 0..ny - 1 {
        for ix in 0..nx - 1 {
            if (x[ix] as f64) < xmin as f64 + 2.0 * xw as f64
                || (x[ix] as f64) > xmax as f64 - 2.0 * xw as f64
                || y[iy] < ymin + yw
                || y[iy] > ymax - yw
            {
                w[iy * nx + ix] += 1.0
            }
        }
    }
    let mut non = 0;
    for iy in 0..ny {
        for ix in 0..nx {
            let k = iy * nx + ix;
            if iy < ny - 1 && min(z[k], z[k + nx]) <= c && c < max(z[k], z[k + nx]) {
                set(&mut z[k], 2);
                non += 1
            } else {
                clear(&mut z[k], 2)
            }
            if ix < nx - 1 && min(z[k], z[k + 1]) <= c && c < max(z[k], z[k + 1]) {
                set(&mut z[k], 1);
                non += 1
            } else {
                clear(&mut z[k], 1)
            }
        }
    }
    gsave();
    newpath();
    let mut draw = Draw(0);
    let label = Label {
        xw,
        yw,
        size: lcs,
        text: &text,
        font: lcf,
        color: lcc,
        bounds: (xmin, xmax, ymin, ymax),
    };
    for side in 0..5 {
        let count = match side {
            0 | 2 => nx - 1,
            1 | 3 => ny - 1,
            _ => ny.saturating_sub(2) * (nx - 1),
        };
        for i in 0..count {
            if non <= 0 {
                break;
            }
            let (k, b, cell, px, py) = match side {
                0 => {
                    let k = (ny - 1) * nx + i;
                    (k, 1, k.wrapping_sub(nx), i, ny - 1)
                }
                1 => {
                    let k = i * nx + nx - 1;
                    (k, 2, k - 1, nx - 1, i)
                }
                2 => (i, 1, i, i, 0),
                3 => {
                    let k = i * nx;
                    (k, 2, k, 0, i)
                }
                _ => {
                    let iy = i / (nx - 1) + 1;
                    let ix = i % (nx - 1);
                    let k = iy * nx + ix;
                    (k, 1, k, ix, iy)
                }
            };
            if !bit(z[k], b) {
                continue;
            }
            let d = if b == 1 {
                delta(c, z[k], z[k + 1])
            } else {
                delta(c, z[k], z[k + nx])
            };
            if b == 1 {
                draw.mov(interp(d, x[px], x[px + 1]), y[py])
            } else {
                draw.mov(x[px], interp(d, y[py], y[py + 1]))
            }
            clear(&mut z[k], b);
            non -= 1;
            let mut cell = cell;
            if side == 4 {
                if let Some((_, _, more)) = connect(c, x, y, z, &mut cell, &mut draw) {
                    if more {
                        set(&mut z[k], 1);
                        non += 1
                    } else {
                        continue;
                    }
                } else {
                    continue;
                }
            }
            follow(c, x, y, z, w, cell, &mut non, &mut draw, &label);
        }
    }
    stroke();
    grestore();
}
