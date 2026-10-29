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

//! SUWINDPOLY, port of su/main/windowing_sorting_muting/suwindpoly.c (Revision 1.7).
//! Credits: Reginald H. Beardsley; originally suxarb.c adapted from the SLT/SU package.
//! The polygon algorithm is based on Joseph O'Rourke's Computational Geometry in C.
use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUWINDPOLY - WINDow data to extract traces on or within a respective\t",
    "\tPOLYgonal line or POLYgon with coordinates specified by header\t",
    "\tkeyword values \t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  suwindpoly <stdin [Required parameters] [Optional params] file=outfile",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " x=x1,x2,...\tlist of X coordinates for vertices\t\t\t",
    " y=y1,y2,...\tlist of Y coordinates for vertices\t\t\t",
    " file=file1,file2,..\toutput filename(s)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters\t\t\t\t\t\t\t",
    " xkey=fldr\tX coordinate header key\t\t\t\t\t",
    " ykey=ep\tY coordinate header key\t\t\t\t\t",
    " pass=0 \tpolyline mode: pass traces near the polygonal line\t",
    "\t\t=1 pass all traces interior to polygon\t\t\t",
    "\t\t=2 pass all traces exterior to polygon\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters used in polyline pass=0 mode only:\t\t",
    " The following need to be given if the unit increments in the X & Y\t",
    " directions are not equal.  For example, if fldr increments by 1 and\t",
    " ep increments by 4 to form 25 x 25 m bins specify dx=25.0 & dy=6.25.\t",
    " The output binning key will be converted to integers by the scaling\t",
    " with the smaller of the two values.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " dx=1.0\tunit increment distance in X direction\t\t\t",
    " dy=1.0\tunit increment distance in Y direction\t\t\t",
    " ilkey=tracl\tkey for resulting inline index in polyline mode\t\t",
    " xlkey=tracr\tkey for resulting xline index in polyline mode\t\t",
    " dw=1.0\tdistance in X-Y coordinate units of extracted line\t",
    "\t\tto pass points to output.  Width of resulting line is\t",
    "\t\t2*dw.  Ignored if polygon mode is specified.\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " In polyline mode (pass=0), a single trace may be output multiple times",
    " if it meets the acceptance criteria (distance from line segment < dw)\t",
    " for multiple line segments. However, the headers will be distinct\t",
    " and will associate the output trace with a line segment. This\t\t",
    " behavior facilitates creation of 3D supergathers from polyline\t",
    " output. Use susort after running in polyline mode.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " x=& y=lists should be repeated for as many polygons as needed when  ",
    " pass=1 or pass=2. \t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " In polygon mode, the polygon closes itself from the last vertex to\t",
    " the first.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    "  suwindpoly <input.su x=10,20,50 y=0,30,60 dw=10 pass=0 file=out.su\t",
    "\t\t\t\t\t\t\t\t\t",
];
fn in_polygon(x0: f32, y0: f32, x: &[f32], y: &[f32]) -> i32 {
    let n = x.len();
    let (mut r, mut l) = (0, 0);
    for i in 0..n {
        if x[i] == x0 && y[i] == y0 {
            return 1;
        }
        let j = (i + n - 1) % n;
        if (y[i] > y0 && y[j] <= y0) || (y[j] > y0 && y[i] <= y0) {
            let v = ((x[i] - x0) * (y[j] - y0) - (x[j] - x0) * (y[i] - y0))
                / ((y[j] - y0) - (y[i] - y0));
            if v > 0.0 {
                r += 1
            }
            if v < 0.0 {
                l += 1
            }
        }
    }
    if r % 2 != l % 2 {
        2
    } else if r % 2 == 1 {
        3
    } else {
        0
    }
}
fn polyline(x0: f32, y0: f32, x: &[f32], y: &[f32], s: &[f32]) -> Vec<(f32, f32)> {
    let mut out = Vec::new();
    for i in 1..x.len() {
        if x[i] != x[i - 1] {
            let b = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            let a = y[i - 1] - x[i - 1] * b;
            let xx = ((x0 + (y0 - a) * b) as f64 / (1.0 + (b * b) as f64)) as f32;
            let yy = a + b * xx;
            let m1 = (yy - y[i - 1]) / (xx - x[i - 1]);
            let m2 = (y0 - y[i - 1]) / (x0 - x[i - 1]);
            let m = ((m1 - m2) as f64 / (1.0 + (m1 * m2) as f64)) as f32;
            if (x[i] >= xx && xx >= x[i - 1]) || (x[i] <= xx && xx <= x[i - 1]) {
                let il = (s[i - 1] as f64
                    + (((yy - y[i - 1]) * (yy - y[i - 1]) + (xx - x[i - 1]) * (xx - x[i - 1]))
                        as f64)
                        .sqrt()) as f32;
                let d = (((x0 - xx) * (x0 - xx) + (y0 - yy) * (y0 - yy)) as f64).sqrt() as f32;
                out.push((il, if m >= 0.0 { d } else { -d }));
            }
        } else if (y[i] >= y0 && y0 >= y[i - 1]) || (y[i] <= y0 && y0 <= y[i - 1]) {
            let il = (s[i - 1] as f64 + ((y0 - y[i - 1]) as f64).abs()) as f32;
            let d = ((x0 - x[i]) as f64).abs() as f32;
            out.push((il, if y0 >= y[i - 1] && x0 > x[i] { -d } else { d }));
        }
    }
    out
}
fn changeval(t: &str, f: i32) -> Value {
    let mut v = Value::default();
    match t {
        "s" => err!("can't change char header word"),
        "h" => v.set_h(f as i16),
        "u" => v.set_u(f as u16),
        "l" => v.set_l(f as i64),
        "v" => v.set_v(f as u64),
        "i" => v.set_i(f),
        "p" => v.set_p(f as u32),
        "f" => v.set_f(f as f32),
        "d" => v.set_d(f as f64),
        _ => err!("unknown type %s", t),
    }
    v
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut xkey, mut ykey, mut ilkey, mut xlkey) = (
        "fldr".to_string(),
        "ep".to_string(),
        "tracl".to_string(),
        "tracr".to_string(),
    );
    getparstring("xkey", &mut xkey);
    getparstring("ykey", &mut ykey);
    getparstring("ilkey", &mut ilkey);
    getparstring("xlkey", &mut xlkey);
    let (mut dw, mut dx, mut dy) = (1.0f32, 1.0f32, 1.0f32);
    getparfloat("dw", &mut dw);
    getparfloat("dx", &mut dx);
    getparfloat("dy", &mut dy);
    let mut pass = 0;
    getparint("pass", &mut pass);
    checkpars();
    let nx = countparname("x");
    let ny = countparname("y");
    let np = if pass == 0 {
        if nx > 1 || ny > 1 {
            err!("only one set of x=& y=permitted if pass=0")
        }
        if nx != 1 {
            err!("x values missing for polyline")
        }
        if ny != 1 {
            err!("y values missing for polyline")
        }
        1
    } else {
        if nx > ny {
            err!("y values missing for polygon")
        }
        if ny > nx {
            err!("x values missing for polygon")
        }
        nx
    };
    let nf = countparname("file");
    if pass == 1 && nf != np {
        err!("One output filename must be given for each polygon")
    }
    if pass == 2 && nf != 1 {
        err!("Only one output filename may be given when pass=2")
    }
    if pass > 2 || pass < 0 {
        err!("Invalid pass=option")
    }
    let mut names = Vec::new();
    for i in 0..nf {
        let mut name = String::new();
        getnparstring(i + 1, "file", &mut name);
        names.push(name);
    }
    let mut polys = Vec::new();
    for i in 0..np {
        let nx = countnparval(i + 1, "x");
        let ny = countnparval(i + 1, "y");
        if nx > ny {
            err!("Y values missing from polygon %d", i + 1)
        }
        if ny > nx {
            err!("X values missing from polygon %d", i + 1)
        }
        let mut x = vec![0.0f32; nx as usize];
        let mut y = vec![0.0f32; nx as usize];
        getnparfloat(i + 1, "x", &mut x);
        getnparfloat(i + 1, "y", &mut y);
        if pass != 0 && !x.is_empty() && (x[x.len() - 1] != x[0] || y[y.len() - 1] != y[0]) {
            x.push(x[0]);
            y.push(y[0]);
        } else {
            for j in 0..x.len() {
                x[j] *= dx;
                y[j] *= dy;
            }
        }
        let mut s = vec![0.0f32; x.len()];
        for j in 1..x.len() {
            s[j] = (s[j - 1] as f64
                + (((x[j] - x[j - 1]) as f64).powf(2.0) + ((y[j] - y[j - 1]) as f64).powf(2.0))
                    .sqrt()) as f32;
        }
        polys.push((x, y, s));
    }
    for (i, (x, y, s)) in polys.iter().enumerate() {
        if pass == 1 {
            fprintf!(stderr(), "Passing interior of ");
        } else if pass == 2 {
            fprintf!(stderr(), "Passing exterior of ");
        }
        fprintf!(stderr(), "Polygon: %4d ", (i + 1) as i32);
        if names.len() > i {
            fprintf!(stderr(), "Filename: %s ", names[i].as_str());
        }
        warn!("");
        for j in 0..x.len() {
            fprintf!(stderr(), "%14s", "");
            fprintf!(stderr(), "Vertex: %d ", (j + 1) as i32);
            fprintf!(stderr(), "X: %f ", x[j]);
            fprintf!(stderr(), "Y: %f ", y[j]);
            fprintf!(stderr(), "S: %f ", s[j]);
            fprintf!(stderr(), "\n");
        }
        fprintf!(stderr(), "\n");
    }
    let mut files = Vec::new();
    for name in &names {
        match fopen(name, "w") {
            Some(fp) => files.push(fp),
            None => err!("fopen(3c) failed %s %d", "suwindpoly.c", 387),
        }
    }
    let mut tr = Trace::new();
    if fgettr(&stdin(), &mut tr) == 0 {
        err!("can't get first trace")
    }
    let (xt, yt, ilt, xlt) = (hdtype(&xkey), hdtype(&ykey), hdtype(&ilkey), hdtype(&xlkey));
    let (xi, yi, ili, xli) = (
        getindex(&xkey),
        getindex(&ykey),
        getindex(&ilkey),
        getindex(&xlkey),
    );
    loop {
        let (mut xv, mut yv) = (Value::default(), Value::default());
        gethval(&tr, xi, &mut xv);
        gethval(&tr, yi, &mut yv);
        let (mut x0, mut y0) = (vtof(xt, xv), vtof(yt, yv));
        if pass == 0 {
            x0 *= dx;
            y0 *= dy;
            let (x, y, s) = &polys[0];
            for (il, xl) in polyline(x0, y0, x, y, s) {
                if (xl as f64).abs() <= dw as f64 {
                    // The original dereferences a missing output file here (undefined behavior).
                    if files.is_empty() {
                        err!("output filename required for selected polyline trace")
                    }
                    puthval(&mut tr, ili, &changeval(ilt, il as i32));
                    puthval(&mut tr, xli, &changeval(xlt, xl as i32));
                    fputtr(&files[0], &mut tr);
                }
            }
        } else if pass == 1 {
            for (i, (x, y, _)) in polys.iter().enumerate() {
                if in_polygon(x0, y0, x, y) != 0 {
                    fputtr(&files[i], &mut tr);
                }
            }
        } else if pass == 2 {
            let mut inside = false;
            for (x, y, _) in &polys {
                if in_polygon(x0, y0, x, y) != 0 {
                    inside = true;
                    break;
                }
            }
            if !inside {
                fputtr(&files[0], &mut tr);
            }
        }
        if fgettr(&stdin(), &mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
