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
//! SULCTHW - Linear Coordinate Transformation of Header Words.
//! Port of su/main/headers/sulcthw.c (Revision 1.4).
//! Credits: Potash Corporation of Saskatchewan: Balasz Nemeth, c. 2008.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SULCTHW - Linear Coordinate Transformation of Header Words\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   sulcthw <infile >outfile\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " xt=0.0\tTranslation of X\t\t\t\t\t",
    " yt=0.0\tTranslation of Y\t\t\t\t\t",
    " zt=0.0\tTranslation of Z\t\t\t\t\t",
    " xr=0.0\tRotation around X in degrees\t \t\t\t",
    " yr=0.0\tRotation aroun Y  in degrees\t \t\t\t",
    " zr=0.0\tRotation around Z in degrees \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Translation:\t\t\t\t\t\t\t\t",
    " x = x'+ xt;y = y'+ yt;z = z' + zt;\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Rotations:\t\t\t\t\t  \t\t\t",
    " Around Z axis\t\t\t\t\t\t\t\t",
    " X = x*cos(zr)+y*sin(zr);\t\t\t  \t\t\t",
    " Y = y*cos(zr)-x*sin(zr);\t\t\t  \t\t\t",
    " Around Y axis\t\t\t\t\t\t\t\t",
    " Z = z*cos(yr)+x*sin(yr);\t\t\t  \t\t\t",
    " X = x*cos(yr)-z*sin(yr);\t\t\t  \t\t\t",
    " Around X axis\t\t\t\t\t\t\t\t",
    " Y = y*cos(xr)+z*sin(xr);\t\t\t  \t\t\t",
    " Z = Z*cos(xr)-y*sin(xr);\t\t\t  \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Header words triplets that are transformed\t\t\t\t",
    " sx,sy,selev\t\t\t\t\t\t\t\t",
    " gx,gy,gelev\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The header words restored as 32 bit integers using SEG-Y\t\t",
    " convention (with coordinate scalers scalco and scalel).\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " After transformation they are converted back to integers and stored.\t",
    "                       \t\t\t \t\t\t",
];

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_segy_scalar::*;

#[derive(Clone, Copy)]
struct Ctrp { x: f64, y: f64, z: f64 }

fn transf(mut p: Ctrp, t: Ctrp, r: Ctrp) -> Ctrp {
    p.x += t.x; p.y += t.y; p.z += t.z;
    let mut pt = Ctrp {
        x: p.x * r.z.cos() + p.y * r.z.sin(),
        y: p.y * r.z.cos() - p.x * r.z.sin(),
        z: p.z,
    };
    pt.z = pt.z * r.y.cos() + pt.x * r.y.sin();
    pt.x = pt.x * r.y.cos() - pt.z * r.y.sin();
    pt.y = pt.y * r.x.cos() + pt.z * r.x.sin();
    pt.z = pt.z * r.x.cos() - pt.y * r.x.sin();
    pt
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tmp = 0.0f32;
    if getparfloat("xt", &mut tmp) == 0 { tmp = 0.0; }
    let tx = tmp as f64;
    if getparfloat("yt", &mut tmp) == 0 { tmp = 0.0; }
    let ty = tmp as f64;
    if getparfloat("zt", &mut tmp) == 0 { tmp = 0.0; }
    let t = Ctrp { x: tx, y: ty, z: tmp as f64 };
    if getparfloat("xr", &mut tmp) == 0 { tmp = 0.0; }
    let rx = tmp as f64;
    if getparfloat("yr", &mut tmp) == 0 { tmp = 0.0; }
    let ry = tmp as f64;
    if getparfloat("zr", &mut tmp) == 0 { tmp = 0.0; }
    let mut r = Ctrp { x: rx, y: ry, z: tmp as f64 };
    checkpars();
    r.x *= 0.017453292; r.y *= 0.017453292; r.z *= 0.017453292;
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    loop {
        let mut cotmp = from_segy_elco_multiplier(tr.scalco);
        let mut eltmp = from_segy_elco_multiplier(tr.scalel);
        let p = Ctrp { x: tr.sx as f64 * cotmp, y: tr.sy as f64 * cotmp, z: tr.selev as f64 * eltmp };
        let pt = transf(p, t, r);
        let q = Ctrp { x: tr.gx as f64 * cotmp, y: tr.gy as f64 * cotmp, z: tr.gelev as f64 * eltmp };
        let qt = transf(q, t, r);
        let mut c = [p.x, p.y, q.x, q.y, 0.0, 0.0, 0.0];
        tr.scalco = elco_scalar(&c[..4]);
        cotmp = to_segy_elco_multiplier(tr.scalco);
        tr.sx = nint(pt.x * cotmp); tr.sy = nint(pt.y * cotmp);
        tr.gx = nint(qt.x * cotmp); tr.gy = nint(qt.y * cotmp);
        c[0] = pt.z; c[1] = qt.z; c[2] = tr.sdepth as f64 * eltmp;
        c[3] = tr.gdel as f64 * eltmp; c[4] = tr.sdel as f64 * eltmp;
        c[5] = tr.swdep as f64 * eltmp; c[6] = tr.gwdep as f64 * eltmp;
        tr.scalel = elco_scalar(&c);
        eltmp = to_segy_elco_multiplier(tr.scalco);
        tr.selev = nint(c[0] * eltmp); tr.gelev = nint(c[1] * eltmp);
        tr.sdepth = nint(c[2] * eltmp); tr.gdel = nint(c[3] * eltmp);
        tr.sdel = nint(c[4] * eltmp); tr.swdep = nint(c[5] * eltmp);
        tr.gwdep = nint(c[6] * eltmp);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
