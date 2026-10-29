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

//! ELAPS - plot a triangulated elastic model. Port of Trielas/graphics/psplot/elaps.c.
//! Credits: Dave Hale and Andreas Rueger, Colorado School of Mines.

use oxwave_core::prelude::*;
use oxwave_f_basic::*;
use oxwave_f_psaxesbox::*;
use oxwave_lib_tri::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " ELAPS - plot a triangulated function p(x,z) via PostScript \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " elaps <modelfile >postscriptfile [optional parameters] \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters: \t\t\t\t\t\t\t",
    " p=0\t\t\tplot sqrt(a3333) (vertical P-wave velocity)\t",
    " \t\t     =1 plot sqrt(a1313)(vertical S-wave velocity)\t",
    " \t\t     =2 plot (v_Ph-v_Pv)/v_Pv   \t\t\t",
    "                    =3 plot (a1212-a1313)/(2*a1313) = gamma            ",
    " gedge=0.0\t\tgray to draw fixed edges (in interval [0.0,1.0])",
    " gtri=1.0\t\tgray to draw non-fixed edges of triangles \t",
    "\t\t\t   put negative number for not drawing  \t",
    " gmin=0.0\t\tmin gray to shade triangles (in interval [0.0,1.0])",
    " gmax=1.0\t\tmax gray to shade triangles (in interval [0.0,1.0])",
    " pgmin=minimum p(x,z)\tp(x,y)-value corresponding to gmin \t\t",
    " pgmax=maximum p(x,z)\tp(x,y)-value corresponding to gmax \t\t",
    " xbox=1.5\t\toffset in inches of left side of axes box \t",
    " ybox=1.5\t\toffset in inches of bottom side of axes box\t",
    " wbox=6.0\t\twidth in inches of axes box\t\t\t",
    " hbox=8.0\t\theight in inches of axes box\t\t\t",
    " xbeg=xmin\t\tvalue at which x axis begins\t\t\t",
    " xend=xmax\t\tvalue at which x axis ends\t\t\t",
    " dxnum=0.0\t\tnumbered tic interval on x axis (0.0 for automatic)",
    " fxnum=xmin\t\tfirst numbered tic on x axis (used if dxnum not 0.0)",
    " nxtic=1\t\tnumber of tics per numbered tic on x axis\t",
    " gridx=none\t\tgrid lines on x axis - none, dot, dash, or solid",
    " labelx=\t\tlabel on x axis\t\t\t\t\t",
    " zbeg=zmin\t\tvalue at which z axis begins\t\t\t",
    " zend=zmax\t\tvalue at which z axis ends\t\t\t",
    " dznum=0.0\t\tnumbered tic interval on z axis (0.0 for automatic)",
    " fznum=zmin\t\tfirst numbered tic on z axis (used if dynum not 0.0)",
    " nztic=1\t\tnumber of tics per numbered tic on z axis\t",
    " gridz=none\t\tgrid lines on z axis - none, dot, dash, or solid",
    " labelz=\t\tlabel on z axis\t\t\t\t\t",
    " labelfont=Helvetica\tfont name for axes labels\t\t\t",
    " labelsize=12\t\tfont size for axes labels\t\t\t",
    " title=\t\ttitle of plot\t\t\t\t\t",
    " titlecolor=black      color of title                                  ",
    " axescolor=black       color of axes                                   ",
    " gridcolor=black       color of grid                                   ",
    " titlefont=Helvetica-Bold font name for title\t\t\t\t",
    " titlesize=24\t\tfont size for title\t\t\t\t",
    " style=seismic\t  normal (z axis horizontal, x axis vertical) or\t",
    "\t\t\tseismic (z axis vertical, x axis horizontal)\t",
    "\t\t\t\t\t\t\t\t\t",
];

extern "C" {
    fn fdopen(fd: i32, mode: *const i8) -> *mut std::ffi::c_void;
    fn dup(fd: i32) -> i32;
}
#[repr(C)]
struct FaceAttributes {
    a1111: f32,
    a3333: f32,
    a1133: f32,
    a1313: f32,
    a1113: f32,
    a3313: f32,
    a1212: f32,
    a2323: f32,
    a1223: f32,
    rho: f32,
    mindex: i32,
}
fn floatpar(name: &str, default: f32) -> f32 {
    let mut v = default;
    getparfloat(name, &mut v);
    v
}
fn intpar(name: &str, default: i32) -> i32 {
    let mut v = default;
    getparint(name, &mut v);
    v
}
fn stringpar(name: &str, default: &str) -> String {
    let mut v = default.to_string();
    getparstring(name, &mut v);
    v
}
fn grid(s: &str) -> i32 {
    match s {
        "dot" => 1,
        "dash" => 2,
        "solid" => 3,
        _ => 0,
    }
}
unsafe fn value(fa: *mut std::ffi::c_void, p: i32, minmax: bool) -> f32 {
    unsafe {
        let a = &*(fa as *const FaceAttributes);
        match p {
            0 => (a.a3333 as f64).sqrt() as f32,
            1 => (a.a1313 as f64).sqrt() as f32,
            2 => (((a.a1111 as f64).sqrt() - (a.a3333 as f64).sqrt()) / a.a3333 as f64) as f32,
            _ => {
                if minmax {
                    (0.5 * (a.a1212 - a.a1313) as f64 / a.a1313 as f64) as f32
                } else {
                    (((a.a1212 - a.a1313) / a.a1313) as f64 * 0.5) as f32
                }
            }
        }
    }
}
unsafe fn drawtri(
    t: *mut Face,
    xscale: f32,
    zscale: f32,
    gedge: f32,
    gtri: f32,
    gbase: f32,
    gscale: f32,
    p: i32,
) {
    unsafe {
        let eu = (*t).eu;
        let e2 = (*eu).euCW;
        let e3 = (*eu).euCCW;
        let v1 = &*(*(*eu).vu).v;
        let v2 = &*(*(*e2).vu).v;
        let v3 = &*(*(*e3).vu).v;
        let (x1, z1, x2, z2, x3, z3) = (
            v1.y * xscale,
            v1.x * zscale,
            v2.y * xscale,
            v2.x * zscale,
            v3.y * xscale,
            v3.x * zscale,
        );
        let fa = (*t).fa;
        let mut g = 0.0;
        if !fa.is_null() {
            g = gbase + value(fa, p, false) * gscale;
        }
        gsave();
        if !fa.is_null() {
            gsave();
            newpath();
            moveto(z1, x1);
            lineto(z2, x2);
            lineto(z3, x3);
            closepath();
            setgray(g);
            fill();
            grestore();
        }
        for (e, za, xa, zb, xb) in [
            (eu, z1, x1, z2, x2),
            (e2, z2, x2, z3, x3),
            (e3, z3, x3, z1, x1),
        ] {
            let fixed = (*(*e).e).fixed != 0;
            if fixed && gedge >= 0.0 && gedge <= 1.0 {
                setgray(gedge);
                moveto(za, xa);
                lineto(zb, xb);
                stroke();
            } else if !fixed && gtri >= 0.0 && gtri <= 1.0 {
                setgray(gtri);
                moveto(za, xa);
                lineto(zb, xb);
                stroke();
            }
        }
        grestore();
    }
}
fn main() {
    unsafe {
        let argv = args();
        initargs(&argv);
        requestdoc(0, SDOC);
        let m = readModel(fdopen(dup(0), c"r".as_ptr()));
        let p = intpar("p", 0);
        match p {
            0 => fprintf!(stderr(), "\n **** v_P(x,z) field will be plotted ****\n\n"),
            1 => fprintf!(stderr(), "\n **** v_S(x,z) field will be plotted ****\n\n"),
            2 => fprintf!(stderr(), "\n **** epsilon field will be plotted ****\n\n"),
            3 => fprintf!(stderr(), "\n **** gamma field will be plotted ****\n\n"),
            _ => err!("ERROR: parameter <p> must be either 0,1,2 or 3 "),
        };
        let mut pmin = f32::MAX;
        let mut pmax = -f32::MAX;
        let mut f = (*m).f;
        loop {
            let fa = (*f).fa;
            if !fa.is_null() {
                let v = value(fa, p, true);
                if v < pmin {
                    pmin = v;
                }
                if v > pmax {
                    pmax = v;
                }
                f = (*f).fNext;
            }
            if f == (*m).f {
                break;
            }
        }
        fprintf!(stderr(), "\n minimum value = %g \n\n", pmin);
        fprintf!(stderr(), "\n maximum value = %g \n\n", pmax);
        let gedge = floatpar("gedge", 0.0);
        let gtri = floatpar("gtri", 1.0);
        let gmin = floatpar("gmin", 0.0);
        let gmax = floatpar("gmax", 1.0);
        let pgmin = floatpar("pgmin", pmin);
        let mut pgmax = floatpar("pgmax", pmax);
        let labelx = stringpar("labelx", "");
        let labelz = stringpar("labelz", "");
        let title = stringpar("title", "");
        let styles = stringpar("style", "seismic");
        let style = if styles == "normal" { 0 } else { 1 };
        let labelfont = stringpar("labelfont", "Helvetica");
        let titlefont = stringpar("titlefont", "Helvetica-Bold");
        let titlecolor = stringpar("titlecolor", "black");
        let axescolor = stringpar("axescolor", "black");
        let gridcolor = stringpar("gridcolor", "black");
        let labelsize = floatpar("labelsize", 18.0);
        let titlesize = floatpar("titlesize", 24.0);
        let mut xbox = floatpar("xbox", 1.5);
        let mut ybox = floatpar("ybox", 1.5);
        let mut wbox = floatpar("wbox", 6.0);
        let mut hbox = floatpar("hbox", 8.0);
        let xbeg = floatpar("xbeg", (*m).ymin);
        let xend = floatpar("xend", (*m).ymax);
        let zbeg = floatpar("zbeg", (*m).xmin);
        let zend = floatpar("zend", (*m).xmax);
        let dxnum = floatpar("dxnum", 0.0);
        let fxnum = floatpar("fxnum", (*m).ymin);
        let nxtic = intpar("nxtic", 1);
        let gridx = grid(&stringpar("gridx", "none"));
        let dznum = floatpar("dznum", 0.0);
        let fznum = floatpar("fznum", (*m).xmin);
        let nztic = intpar("nztic", 1);
        let gridz = grid(&stringpar("gridz", "none"));
        let axeswidth = floatpar("axeswidth", 1.0);
        let ticwidth = floatpar("ticwidth", axeswidth);
        let gridwidth = floatpar("gridwidth", axeswidth);
        beginps();
        newpage("1", 1);
        xbox *= 72.0;
        ybox *= 72.0;
        wbox *= 72.0;
        hbox *= 72.0;
        let mut bbox = [0; 4];
        psAxesBBox(
            xbox, ybox, wbox, hbox, &labelfont, labelsize, &titlefont, titlesize, style, &mut bbox,
        );
        boundingbox(bbox[0], bbox[1], bbox[2], bbox[3]);
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
        if pgmax == pgmin {
            pgmax += 1.0;
        }
        let gscale = (gmin - gmax) / (pgmax - pgmin);
        let gbase = gmax - pgmin * gscale;
        f = (*m).f;
        loop {
            drawtri(f, xscale, zscale, gedge, gtri, gbase, gscale, p);
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
        showpage();
        endps();
        exit(1);
    }
}
