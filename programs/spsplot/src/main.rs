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
//! SPSPLOT, port of tri/graphics/psplot/spsplot.c.
//! Credits: Dave Hale; Craig Artley; Morten Wendell Pedersen.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SPSPLOT - plot a triangulated sloth function s(x,z) via PostScript\t",
    "\t\t\t\t\t\t\t\t\t",
    " spsplot <modelfile >postscriptfile [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " gedge=0.0             gray to draw fixed edges (in interval [0.0,1.0])",
    " gtri=1.0              gray to draw non-fixed edges of triangles \t",
    " gmin=0.0              min gray to shade triangles (in interval [0.0,1.0])",
    " gmax=1.0              max gray to shade triangles (in interval [0.0,1.0])",
    " sgmin=minimum s(x,z)  s(x,y) corresponding to gmin \t\t\t",
    " sgmax=maximum s(x,z)  s(x,y) corresponding to gmax \t\t\t",
    " xbox=1.5              offset in inches of left side of axes box \t",
    " ybox=1.5              offset in inches of bottom side of axes box\t",
    " wbox=6.0              width in inches of axes box\t\t\t",
    " hbox=8.0              height in inches of axes box\t\t\t",
    " xbeg=xmin             value at which x axis begins\t\t\t",
    " xend=xmax             value at which x axis ends\t\t\t",
    " dxnum=0.0             numbered tic interval on x axis (0.0 for automatic)",
    " fxnum=xmin            first numbered tic on x axis (used if dxnum not 0.0)",
    " nxtic=1               number of tics per numbered tic on x axis\t",
    " gridx=none            grid lines on x axis - none, dot, dash, or solid",
    " labelx=               label on x axis\t\t\t\t\t",
    " zbeg=zmin             value at which z axis begins\t\t\t",
    " zend=zmax             value at which z axis ends\t\t\t",
    " dznum=0.0             numbered tic interval on z axis (0.0 for automatic)",
    " fznum=zmin            first numbered tic on z axis (used if dynum not 0.0)",
    " nztic=1               number of tics per numbered tic on z axis\t",
    " gridz=none            grid lines on z axis - none, dot, dash, or solid",
    " labelz=               label on z axis\t\t\t\t\t",
    " labelfont=Helvetica   font name for axes labels\t\t\t",
    " labelsize=12          font size for axes labels\t\t\t",
    " title=                title of plot\t\t\t\t\t",
    " titlefont=Helvetica-Bold  font name for title\t\t\t\t",
    " titlesize=24          font size for title\t\t\t\t",
    " titlecolor=black      color of title\t\t\t\t\t",
    " axescolor=black       color of axes\t\t\t\t\t",
    " gridcolor=black       color of grid\t\t\t\t\t",
    " style=seismic         normal (z axis horizontal, x axis vertical) or\t",
    "                       seismic (z axis vertical, x axis horizontal)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note:  A value of gedge or gtri outside the interval [0.0,1.0]\t",
    " results in that class of edge not being drawn.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

use oxwave_core::prelude::*;
use oxwave_f_basic::*;
use oxwave_f_psaxesbox::*;
use oxwave_lib_tri::*;
use std::ffi::{c_char, c_void};
unsafe extern "C" {
    fn fdopen(fd: i32, mode: *const c_char) -> *mut c_void;
}
#[repr(C)]
struct FaceAttributes {
    s00: f32,
    dsdx: f32,
    dsdz: f32,
    dens: f32,
    qfac: f32,
}
unsafe fn vertices(f: *mut Face) -> [*mut EdgeUse; 3] {
    unsafe {
        let e = (*f).eu;
        [e, (*e).euCW, (*e).euCCW]
    }
}
unsafe fn coords(e: *mut EdgeUse) -> (f32, f32) {
    unsafe {
        let v = (*(*e).vu).v;
        ((*v).y, (*v).x)
    }
}
unsafe fn minmax(m: *mut Model) -> (f32, f32) {
    unsafe {
        let (mut lo, mut hi) = (f32::MAX, -f32::MAX);
        let mut f = (*m).f;
        loop {
            if !(*f).fa.is_null() {
                let fa = &*((*f).fa as *const FaceAttributes);
                let mut e = (*f).eu;
                loop {
                    let (x, z) = coords(e);
                    let s = fa.s00 + x * fa.dsdx + z * fa.dsdz;
                    if s < lo {
                        lo = s;
                    }
                    if s > hi {
                        hi = s;
                    }
                    e = (*e).euCW;
                    if e == (*f).eu {
                        break;
                    }
                }
            }
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
        (lo, hi)
    }
}
fn shaded(
    x1: f32,
    y1: f32,
    g1: f32,
    mut x2: f32,
    mut y2: f32,
    g2: f32,
    mut x3: f32,
    mut y3: f32,
    g3: f32,
) {
    gsave();
    newpath();
    moveto(x1, y1);
    lineto(x2, y2);
    lineto(x3, y3);
    clip();
    translate(x1, y1);
    x2 -= x1;
    y2 -= y1;
    x3 -= x1;
    y3 -= y1;
    let a = y2 * (g3 - g1) - y3 * (g2 - g1);
    let b = x3 * (g2 - g1) - x2 * (g3 - g1);
    let c = x2 * y3 - y2 * x3;
    let dx = if c != 0.0 { -a / c } else { 0.0 };
    let dy = if c != 0.0 { -b / c } else { 0.0 };
    let angle = if dx != 0.0 || dy != 0.0 {
        (dy as f64).atan2(dx as f64) as f32
    } else {
        0.0
    };
    let ca = (angle as f64).cos() as f32;
    let sa = (angle as f64).sin() as f32;
    rotate((angle as f64 * 180.0 / 3.141592654) as f32);
    let (x1r, y1r) = (0.0f32, 0.0f32);
    let x2r = x2 * ca + y2 * sa;
    let y2r = -x2 * sa + y2 * ca;
    let x3r = x3 * ca + y3 * sa;
    let y3r = -x3 * sa + y3 * ca;
    let xmin = if x1r < x2r { x1r } else { x2r };
    let xmin = if xmin < x3r { xmin } else { x3r };
    let ymin = if y1r < y2r { y1r } else { y2r };
    let ymin = if ymin < y3r { ymin } else { y3r };
    let gmin = if xmin == x1r {
        g1
    } else if xmin == x2r {
        g2
    } else {
        g3
    };
    let xmax = if x1r > x2r { x1r } else { x2r };
    let xmax = if xmax > x3r { xmax } else { x3r };
    let ymax = if y1r > y2r { y1r } else { y2r };
    let ymax = if ymax > y3r { ymax } else { y3r };
    let gmax = if xmax == x1r {
        g1
    } else if xmax == x2r {
        g2
    } else {
        g3
    };
    let eps = 0.0166666f32;
    let nlines = ((1.0f32 + (gmax - gmin).abs() / eps) as i32).max(2);
    let gs = (gmax - gmin) / (nlines - 1) as f32;
    let xs = (xmax - xmin) / (nlines - 1) as f32;
    setlinewidth((1.1f64 * xs.abs() as f64) as f32);
    newpath();
    let (mut gray, mut x) = (gmin, xmin);
    for _ in 0..nlines {
        setgray(gray);
        moveto(x, ymin);
        lineto(x, ymax);
        stroke();
        gray += gs;
        x += xs;
    }
    grestore();
}
unsafe fn interior(f: *mut Face, xs: f32, zs: f32, base: f32, scale: f32) {
    unsafe {
        if (*f).fa.is_null() {
            return;
        }
        let fa = &*((*f).fa as *const FaceAttributes);
        let es = vertices(f);
        let mut p = [(0.0, 0.0, 0.0); 3];
        for i in 0..3 {
            let (x, z) = coords(es[i]);
            let s = fa.s00 + x * fa.dsdx + z * fa.dsdz;
            p[i] = (z * zs, x * xs, base + s * scale);
        }
        shaded(
            p[0].0, p[0].1, p[0].2, p[1].0, p[1].1, p[1].2, p[2].0, p[2].1, p[2].2,
        );
    }
}
unsafe fn edges(f: *mut Face, xs: f32, zs: f32, fixed: bool) {
    unsafe {
        let es = vertices(f);
        let mut p = [(0.0, 0.0); 3];
        for i in 0..3 {
            let (x, z) = coords(es[i]);
            p[i] = (z * zs, x * xs);
        }
        for i in 0..3 {
            if ((*(*es[i]).e).fixed != 0) == fixed {
                moveto(p[i].0, p[i].1);
                lineto(p[(i + 1) % 3].0, p[(i + 1) % 3].1);
            }
        }
    }
}
fn grid(s: &str) -> i32 {
    match s {
        "dot" => 1,
        "dash" => 2,
        "solid" => 3,
        _ => 0,
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    unsafe {
        let m = readModel(fdopen(0, c"r".as_ptr()));
        let (smin, smax) = minmax(m);
        macro_rules! fp {
            ($n:ident,$v:expr) => {
                let mut $n = $v;
                getparfloat(stringify!($n), &mut $n);
            };
        }
        macro_rules! ip {
            ($n:ident,$v:expr) => {
                let mut $n = $v;
                getparint(stringify!($n), &mut $n);
            };
        }
        macro_rules! sp {
            ($n:ident,$v:expr) => {
                let mut $n = String::from($v);
                getparstring(stringify!($n), &mut $n);
            };
        }
        fp!(gedge, 0.0);
        fp!(gtri, 1.0);
        fp!(gmin, 0.0);
        fp!(gmax, 1.0);
        fp!(sgmin, smin);
        fp!(sgmax, smax);
        sp!(labelx, "");
        sp!(labelz, "");
        sp!(title, "");
        sp!(style, "seismic");
        let style = if style == "normal" { 0 } else { 1 };
        sp!(labelfont, "Helvetica");
        sp!(titlefont, "Helvetica-Bold");
        fp!(labelsize, 18.0);
        fp!(titlesize, 24.0);
        sp!(titlecolor, "black");
        sp!(axescolor, "black");
        sp!(gridcolor, "black");
        fp!(xbox, 1.5);
        fp!(ybox, 1.5);
        fp!(wbox, 6.0);
        fp!(hbox, 8.0);
        fp!(xbeg, (*m).ymin);
        fp!(xend, (*m).ymax);
        fp!(zbeg, (*m).xmin);
        fp!(zend, (*m).xmax);
        fp!(dxnum, 0.0);
        fp!(fxnum, (*m).ymin);
        ip!(nxtic, 1);
        sp!(gridx, "none");
        let gridx = grid(&gridx);
        fp!(dznum, 0.0);
        fp!(fznum, (*m).xmin);
        ip!(nztic, 1);
        sp!(gridz, "none");
        let gridz = grid(&gridz);
        fp!(axeswidth, 1.0);
        fp!(ticwidth, axeswidth);
        fp!(gridwidth, axeswidth);
        xbox *= 72.0;
        ybox *= 72.0;
        wbox *= 72.0;
        hbox *= 72.0;
        let mut bbox = [0; 4];
        psAxesBBox(
            xbox, ybox, wbox, hbox, &labelfont, labelsize, &titlefont, titlesize, style, &mut bbox,
        );
        boundingbox(bbox[0], bbox[1], bbox[2], bbox[3]);
        begineps();
        gsave();
        rectclip(xbox, ybox, wbox, hbox);
        let xsize = if style == 0 { hbox } else { wbox };
        let zsize = if style == 0 { wbox } else { hbox };
        translate(xbox, ybox);
        if style != 0 {
            rotate(-90.0);
            translate(-hbox, 0.0);
        }
        let xscale = xsize / (xend - xbeg);
        let zscale = zsize / (zend - zbeg);
        translate(-zbeg * zscale, -xbeg * xscale);
        if sgmax == sgmin {
            sgmax += 1.0;
        }
        let gscale = (gmax - gmin) / (sgmax - sgmin);
        let gbase = gmin - sgmin * gscale;
        gsave();
        let mut f = (*m).f;
        loop {
            interior(f, xscale, zscale, gbase, gscale);
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
        grestore();
        if gtri >= 0.0 && gtri <= 1.0 {
            gsave();
            setgray(gtri);
            let mut f = (*m).f;
            loop {
                edges(f, xscale, zscale, false);
                f = (*f).fNext;
                if f == (*m).f {
                    break;
                }
            }
            stroke();
            grestore();
        }
        if gedge >= 0.0 && gedge <= 1.0 {
            gsave();
            setgray(gedge);
            let mut f = (*m).f;
            loop {
                edges(f, xscale, zscale, true);
                f = (*f).fNext;
                if f == (*m).f {
                    break;
                }
            }
            stroke();
            grestore();
        }
        grestore();
        psAxesBox(
            xbox,
            ybox,
            wbox,
            hbox,
            zbeg,
            zend,
            0.0,
            0.0,
            dznum,
            fznum,
            nztic,
            gridz,
            &labelz,
            xbeg,
            xend,
            0.0,
            0.0,
            dxnum,
            fxnum,
            nxtic,
            gridx,
            &labelx,
            &labelfont,
            labelsize,
            &title,
            &titlefont,
            titlesize,
            &titlecolor,
            &axescolor,
            &gridcolor,
            ticwidth,
            axeswidth,
            gridwidth,
            style,
        );
        showpage();
        endeps();
    }
    exit(CWP_Exit());
}
