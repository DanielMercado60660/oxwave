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

//! PSMERGE - merge PostScript files. Port of psplot/main/psmerge.c (Revision 1.10).
//! Credits: Dave Hale and Craig Artley, Colorado School of Mines.

use oxwave_core::prelude::*;
use oxwave_f_basic::*;
use std::ffi::CString;

unsafe extern "C" {
    fn sscanf(s: *const i8, format: *const i8, ...) -> i32;
}
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " PSMERGE - MERGE PostScript files\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " psmerge in= [optional parameters] >postscriptfile\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " in=                    postscript file to merge\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " origin=0.0,0.0         x,y origin in inches\t\t\t\t",
    " scale=1.0,1.0          x,y scale factors\t\t\t\t",
    " rotate=0.0             rotation angle in degrees\t\t\t",
    " translate=0.0,0.0      x,y translation in inches\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " More than one set of in, origin, scale, rotate, and translate\t\t",
    " parameters may be specified.  Output x and y coordinates are\t\t",
    " determined by:\t\t\t\t\t\t\t",
    "          x = tx + (x-ox)*sx*cos(d) - (y-oy)*sy*sin(d)\t\t\t",
    "          y = ty + (x-ox)*sx*sin(d) + (y-oy)*sy*cos(d)\t\t\t",
    " where tx,ty are translate coordinates, ox,oy are origin coordinates,\t",
    " sx,sy are scale factors, and d is the rotation angle.  Note that the\t",
    " order of operations is shift (origin), scale, rotate, and translate.\t",
    " \t\t\t\t\t\t\t\t\t",
    " If the number of occurrences of a given parameter is less than the number",
    " of input files, then the last occurrence of that parameter will apply to",
    " all subsequent files.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

#[derive(Clone, Copy)]
struct Transform {
    sx: f32,
    sy: f32,
    tx: f32,
    ty: f32,
    ox: f32,
    oy: f32,
    d: f32,
}

fn updatebbox(t: Transform, b: &mut [i32; 4]) {
    let s = (std::f64::consts::PI * t.d as f64 / 180.0).sin() as f32;
    let c = (std::f64::consts::PI * t.d as f64 / 180.0).cos() as f32;
    let [llx, lly, urx, ury] = *b;
    let x1 = c * (llx as f32 - t.ox) * t.sx - s * (lly as f32 - t.oy) * t.sy + t.tx;
    let y1 = s * (llx as f32 - t.ox) * t.sx + c * (lly as f32 - t.oy) * t.sy + t.ty;
    let x2 = c * (urx as f32 - t.ox) * t.sx - s * (lly as f32 - t.oy) * t.sy + t.tx;
    let y2 = s * (urx as f32 - t.ox) * t.sx + c * (lly as f32 - t.oy) * t.sy + t.ty;
    let x3 = c * (urx as f32 - t.ox) * t.sx - s * (ury as f32 - t.oy) * t.sy + t.tx;
    let y3 = s * (urx as f32 - t.ox) * t.sx + c * (ury as f32 - t.oy) * t.sy + t.ty;
    let x4 = c * (llx as f32 - t.ox) * t.sx - s * (ury as f32 - t.oy) * t.sy + t.tx;
    let y4 = s * (llx as f32 - t.ox) * t.sx + c * (ury as f32 - t.oy) * t.sy + t.ty;
    fn min(a: f32, b: f32) -> f32 {
        if a < b { a } else { b }
    }
    fn max(a: f32, b: f32) -> f32 {
        if a > b { a } else { b }
    }
    *b = [
        min(x1, min(x2, min(x3, x4))) as i32,
        min(y1, min(y2, min(y3, y4))) as i32,
        max(x1, max(x2, max(x3, x4))) as i32,
        max(y1, max(y2, max(y3, y4))) as i32,
    ];
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let nin = countparname("in");
    if nin <= 0 {
        err!("must specify at least one input file");
    }
    let mut input = Vec::new();
    let mut transforms = Vec::new();
    let mut limits = [i32::MAX, i32::MAX, i32::MIN, i32::MIN];
    for iin in 0..nin {
        let mut name = String::new();
        getnparstring(iin + 1, "in", &mut name);
        let fp = efopen(&name, "r");
        let mut t = transforms.last().copied().unwrap_or(Transform {
            sx: 1.0,
            sy: 1.0,
            tx: 0.0,
            ty: 0.0,
            ox: 0.0,
            oy: 0.0,
            d: 0.0,
        });
        let mut fval = vec![0.0f32; 100];
        if getnparfloat(iin + 1, "rotate", &mut fval) == 1 {
            t.d = fval[0];
        }
        if getnparfloat(iin + 1, "scale", &mut fval) == 2 {
            t.sx = fval[0];
            t.sy = fval[1];
        }
        if getnparfloat(iin + 1, "translate", &mut fval) == 2 {
            t.tx = (fval[0] as f64 * 72.0) as f32;
            t.ty = (fval[1] as f64 * 72.0) as f32;
        }
        if getnparfloat(iin + 1, "origin", &mut fval) == 2 {
            t.ox = (fval[0] as f64 * 72.0) as f32;
            t.oy = (fval[1] as f64 * 72.0) as f32;
        }
        while let Some(line) = fgets(2048, &fp) {
            let line = line.split('\0').next().unwrap();
            if line.starts_with('%') && line.contains("%%BoundingBox:") {
                if line.contains("atend") {
                    continue;
                }
                let c = CString::new(line.split('\0').next().unwrap()).unwrap();
                let mut b = [0i32; 4];
                let n = unsafe {
                    sscanf(
                        c.as_ptr(),
                        c"%*s %d %d %d %d".as_ptr(),
                        &mut b[0] as *mut i32,
                        &mut b[1] as *mut i32,
                        &mut b[2] as *mut i32,
                        &mut b[3] as *mut i32,
                    )
                };
                if n != 4 {
                    let mut f = [0f32; 4];
                    let n = unsafe {
                        sscanf(
                            c.as_ptr(),
                            c"%*s %f %f %f %f".as_ptr(),
                            &mut f[0] as *mut f32,
                            &mut f[1] as *mut f32,
                            &mut f[2] as *mut f32,
                            &mut f[3] as *mut f32,
                        )
                    };
                    if n != 4 {
                        err!("Error reading BoundingBox!\n");
                    }
                    for j in 0..4 {
                        b[j] = if f[j] > 0.0 {
                            (f[j] as f64 + 0.5) as i32
                        } else {
                            (f[j] as f64 - 0.5) as i32
                        };
                    }
                }
                updatebbox(t, &mut b);
                limits[0] = limits[0].min(b[0]);
                limits[1] = limits[1].min(b[1]);
                limits[2] = limits[2].max(b[2]);
                limits[3] = limits[3].max(b[3]);
                break;
            }
        }
        fclose(&fp);
        input.push(name);
        transforms.push(t);
    }
    boundingbox(limits[0], limits[1], limits[2], limits[3]);
    begineps();
    for (name, t) in input.iter().zip(transforms.iter()) {
        let fp = efopen(name, "r");
        gsave();
        translate(t.tx, t.ty);
        rotate(t.d);
        scale(t.sx, t.sy);
        translate(-t.ox, -t.oy);
        while let Some(line) = fgets(2048, &fp) {
            let line = line.split('\0').next().unwrap();
            if line.starts_with('%') || line.starts_with("showpage") {
                continue;
            }
            fputs(line, &stdout());
        }
        grestore();
        fclose(&fp);
    }
    showpage();
    endeps();
    exit(CWP_Exit());
}
