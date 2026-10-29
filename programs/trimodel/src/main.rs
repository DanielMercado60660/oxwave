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

//! TRIMODEL - make a triangulated sloth (1/velocity^2) model.
//! Port of tri/main/trimodel.c (Revision 1.7).
//! Author: Dave Hale, Colorado School of Mines, 02/12/91.
//! Modified: Andreas Rueger, 01/18/93 (density and Q); Craig Artley,
//! 03/27/94 (s00 correction); Boyi Ou, 04/14/95 (normal-ray interfaces).

#![allow(non_snake_case)]

use oxwave_core::prelude::*;
use oxwave_cubic_interp::*;
use oxwave_f_cubicspline::*;
use oxwave_lib_tri::*;
use std::f64::consts::PI;
use std::{ffi::c_void, ptr::null_mut};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " TRIMODEL - make a triangulated sloth (1/velocity^2) model                  \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " trimodel >modelfile [optional parameters] \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " xmin=0.0               minimum horizontal coordinate (x)\t\t",
    " xmax=1.0               maximum horizontal coordinate (x)\t\t",
    " zmin=0.0               minimum vertical coordinate (z)\t\t",
    " zmax=1.0               maximum vertical coordinate (z)\t\t",
    " xedge=                 x coordinates of an edge\t\t\t",
    " zedge=                 z coordinates of an edge\t\t\t",
    " sedge=                 sloth along an edge\t\t\t\t",
    " kedge=                 array of indices used to identify edges\t",
    " normray               0:do not generate parameters 1: does it   \t",
    " normface              specify which interface to shoot rays   \t",
    " nrays                 number of locations to shoot rays      \t        ",
    " sfill=                 x, z, x0, z0, s00, dsdx, dsdz to fill a region\t",
    " densfill=              x, z, dens to fill a region \t\t\t",
    " qfill=                 x, z, Q-factor to fill a region \t\t",
    " maxangle=5.0           maximum angle (deg) between adjacent edge segments",
    "   \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " More than set of xedge, zedge, and sedge parameters may be \t\t",
    " specified, but the numbers of these parameters must be equal. \t",
    " \t\t\t\t\t\t\t\t\t",
    " Within each set, vertices will be connected by fixed edges. \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Edge indices in the k array are used to identify edges \t\t",
    " specified by the x and z parameters.  The first k index \t\t",
    " corresponds to the first set of x and z parameters, the \t\t",
    " second k index corresponds to the second set, and so on. \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " After all vertices and their corresponding sloth values have \t\t",
    " been inserted into the model, the optional sfill parameters \t\t",
    " are used to fill closed regions bounded by fixed edges. \t\t",
    " Let (x,z) denote any point inside a closed region.  Sloth inside \t",
    " this region is determined by s(x,z) = s00+(x-x0)*dsdx+(z-z0)*dsdz.  \t",
    " The (x,z) component of the sfill parameter is used to identify a \t",
    " closed region. \t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

// Attribute layouts from tri/include/sloth.h, passed opaquely to lib-tri.
#[repr(C)]
struct VertexAttributes {
    s: f32,
}
#[repr(C)]
struct EdgeAttributes {
    k: i32,
}
#[repr(C)]
struct EdgeUseAttributes {
    tx: f32,
    tz: f32,
    c: f32,
}
#[repr(C)]
struct FaceAttributes {
    s00: f32,
    dsdx: f32,
    dsdz: f32,
    dens: f32,
    qfac: f32,
}

extern "C" {
    fn fdopen(fd: i32, mode: *const i8) -> *mut c_void;
    fn fflush(fp: *mut c_void) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
    fn raise(sig: i32) -> i32;
}

fn segfault() -> ! {
    // Missing vertex attributes are dereferenced by the original, even
    // before region fills. Preserve its signal instead of invoking Rust UB.
    unsafe {
        signal(11, 0);
        raise(11);
    }
    exit(139)
}

struct Interfaces {
    xfp: Option<CFile>,
    zfp: Option<CFile>,
    xzfp: Option<CFile>,
    afp: Option<CFile>,
}
struct Segments {
    x: Vec<f32>,
    y: Vec<f32>,
    s: Vec<f32>,
    tx: Vec<f32>,
    ty: Vec<f32>,
    c: Vec<f32>,
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut xmin, mut xmax, mut zmin, mut zmax, mut maxangle) =
        (0.0f32, 1.0f32, 0.0f32, 1.0f32, 1.0f32);
    getparfloat("xmin", &mut xmin);
    getparfloat("xmax", &mut xmax);
    getparfloat("zmin", &mut zmin);
    getparfloat("zmax", &mut zmax);
    getparfloat("maxangle", &mut maxangle);
    let (mut normray, mut normface, mut nrays) = (0, 0, 0);
    getparint("normray", &mut normray);
    getparint("normface", &mut normface);
    getparint("nrays", &mut nrays);
    maxangle = (maxangle as f64 * (PI / 180.0)) as f32;
    let nxe = countparname("xedge");
    let nze = countparname("zedge");
    let nse = countparname("sedge");
    if nxe != nze {
        err!("number of xedge must equal number of zedge");
    }
    if nxe != nse {
        err!("number of xedge must equal number of sedge");
    }
    let ne = nxe;
    let (mut xe, mut ze, mut se) = (Vec::new(), Vec::new(), Vec::new());
    for ie in 0..ne {
        let nxe = countnparval(ie + 1, "xedge");
        let nze = countnparval(ie + 1, "zedge");
        let nse = countnparval(ie + 1, "sedge");
        if nxe != nze {
            err!("number of xedge must equal number of zedge");
        }
        if nxe != nse {
            err!("number of xedge must equal number of sedge");
        }
        let (mut x, mut z, mut s) = (
            vec![0.0f32; nxe as usize],
            vec![0.0f32; nxe as usize],
            vec![0.0f32; nxe as usize],
        );
        getnparfloat(ie + 1, "xedge", &mut x);
        getnparfloat(ie + 1, "zedge", &mut z);
        getnparfloat(ie + 1, "sedge", &mut s);
        xe.push(x);
        ze.push(z);
        se.push(s);
    }
    let nke = countparval("kedge");
    if nke > ne {
        err!("more kedge values than edges specified");
    }
    let mut ke = vec![0i32; ne as usize];
    getparint("kedge", &mut ke);
    ke.resize(ne as usize, 0);
    let interfaces = if normray == 1 {
        Some(Interfaces {
            xfp: fopen("XInterface", "w"),
            zfp: fopen("ZInterface", "w"),
            xzfp: fopen("XZInterface", "w"),
            afp: fopen("AInterface", "w"),
        })
    } else {
        None
    };
    let mut sr = Vec::new();
    for ir in 0..countparname("sfill") {
        let mut row = Vec::<f32>::new();
        if getnparfloat(ir + 1, "sfill", &mut row) != 7 {
            err!("7 values must be specified in sfill parameter");
        }
        sr.push(row);
    }
    let mut dptr = Vec::new();
    for id in 0..countparname("densfill") {
        let mut row = Vec::<f32>::new();
        if getnparfloat(id + 1, "densfill", &mut row) != 3 {
            err!("3 values must be specified in densfill parameter");
        }
        dptr.push(row);
    }
    let mut qptr = Vec::new();
    for iq in 0..countparname("qfill") {
        let mut row = Vec::<f32>::new();
        if getnparfloat(iq + 1, "qfill", &mut row) != 3 {
            err!("3 values must be specified in qfill parameter");
        }
        qptr.push(row);
    }
    let mut segments = Vec::new();
    for ie in 0..ne as usize {
        segments.push(smoothLineSegments(
            maxangle,
            &xe[ie],
            &ze[ie],
            &se[ie],
            if ie as i32 == normface {
                interfaces.as_ref()
            } else {
                None
            },
            nrays,
        ));
    }
    unsafe {
        let m = makemod(xmin, zmin, xmax, zmax, &ke, &segments);
        fillsloth(m, &sr);
        if !dptr.is_empty() {
            filldens(m, &dptr);
        }
        if !qptr.is_empty() {
            fillq(m, &qptr);
        }
        // lib-tri uses native libc FILE*, whereas core owns Rust streams.
        let outfp = fdopen(1, c"w".as_ptr());
        writeModel(m, outfp);
        fflush(outfp);
    }
    exit(CWP_Exit());
}

fn smoothLineSegments(
    maxangle: f32,
    xin: &[f32],
    yin: &[f32],
    sin: &[f32],
    interfaces: Option<&Interfaces>,
    nrays: i32,
) -> Segments {
    let nin = xin.len();
    if nin == 1 {
        return Segments {
            x: xin.to_vec(),
            y: yin.to_vec(),
            s: sin.to_vec(),
            tx: Vec::new(),
            ty: Vec::new(),
            c: Vec::new(),
        };
    }
    if nin == 0 {
        segfault();
    }
    if interfaces.is_some() && nrays < 0 {
        // alloc1float receives a negative count in C and returns NULL;
        // the subsequent uou[0] store raises SIGSEGV.
        segfault();
    }
    let mut uin = vec![0.0f32; nin];
    let mut u = vec![0.0f32; nin];
    let mut x = xin.to_vec();
    let mut y = yin.to_vec();
    for i in 1..nin {
        u[i] = (u[i - 1] as f64
            + (((x[i] - x[i - 1]) as f64).powf(2.0) + ((y[i] - y[i - 1]) as f64).powf(2.0)).sqrt())
            as f32;
        uin[i] = u[i];
    }
    let mut uou = vec![
        0.0f32;
        if interfaces.is_some() {
            nrays.max(1) as usize
        } else {
            0
        }
    ];
    // C stores uou[0] even when nrays==0. Reserve that unused slot to
    // reproduce its zero-ray output without an out-of-bounds Rust store.
    if interfaces.is_some() {
        let maxchord = u[nin - 1];
        let dchord = maxchord / nrays as f32;
        for i in 1..nrays as usize {
            uou[i] = uou[i - 1] + dchord;
        }
    }
    let (mut xind, mut yind, mut sind) = (
        vec![[0.0; 4]; nin],
        vec![[0.0; 4]; nin],
        vec![[0.0; 4]; nin],
    );
    csplin(&uin, xin, &mut xind);
    csplin(&uin, yin, &mut yind);
    csplin(&uin, sin, &mut sind);
    let mut i = 1;
    while i + 1 < u.len() {
        smoothTwoSegments(maxangle, i, &uin, &xind, &yind, &mut u, &mut x, &mut y);
        i += 1;
    }
    let n = u.len();
    let (mut s, mut tx, mut ty, mut c) = (
        vec![0.0f32; n],
        vec![0.0f32; n],
        vec![0.0f32; n],
        vec![0.0f32; n],
    );
    for i in 0..n {
        intcub(0, &uin, &sind, &u[i..i + 1], &mut s[i..i + 1]);
        let (mut dx, mut dy, mut ddx, mut ddy) = ([0.0], [0.0], [0.0], [0.0]);
        intcub(1, &uin, &xind, &u[i..i + 1], &mut dx);
        intcub(1, &uin, &yind, &u[i..i + 1], &mut dy);
        intcub(2, &uin, &xind, &u[i..i + 1], &mut ddx);
        intcub(2, &uin, &yind, &u[i..i + 1], &mut ddy);
        let (dx, dy, ddx, ddy) = (dx[0], dy[0], ddx[0], ddy[0]);
        let ds = ((dx * dx + dy * dy) as f64).sqrt() as f32;
        tx[i] = dx / ds;
        ty[i] = dy / ds;
        c[i] = (dx * ddy - dy * ddx) / (ds * ds * ds);
    }
    if let Some(fp) = interfaces {
        let (mut xou, mut you, mut agou) = (
            vec![0.0f32; nrays as usize],
            vec![0.0f32; nrays as usize],
            vec![0.0f32; nrays as usize],
        );
        for i in 0..nrays as usize {
            let (mut dx, mut dy) = ([0.0], [0.0]);
            intcub(1, &uin, &xind, &uou[i..i + 1], &mut dx);
            intcub(1, &uin, &yind, &uou[i..i + 1], &mut dy);
            intcub(0, &uin, &xind, &uou[i..i + 1], &mut xou[i..i + 1]);
            intcub(0, &uin, &yind, &uou[i..i + 1], &mut you[i..i + 1]);
            agou[i] = (90.0 - (-dx[0] as f64).atan2(dy[0] as f64) * 180.0 / PI) as f32;
            fprintf!(
                fp.xzfp.as_ref().unwrap_or_else(|| segfault()),
                "%f %f\n",
                xou[i],
                you[i]
            );
        }
        // fopen failures are unchecked in C. A nonempty write through the
        // resulting null FILE* signals instead of reporting an efopen error.
        if !xou.is_empty() {
            fwrite(&xou, fp.xfp.as_ref().unwrap_or_else(|| segfault()));
            fwrite(&you, fp.zfp.as_ref().unwrap_or_else(|| segfault()));
            fwrite(&agou, fp.afp.as_ref().unwrap_or_else(|| segfault()));
        }
    }
    Segments { x, y, s, tx, ty, c }
}

fn smoothTwoSegments(
    maxangle: f32,
    i: usize,
    ud: &[f32],
    xd: &[[f32; 4]],
    yd: &[[f32; 4]],
    uo: &mut Vec<f32>,
    xo: &mut Vec<f32>,
    yo: &mut Vec<f32>,
) {
    let no = uo.len();
    if i == 0 || i == no - 1 {
        return;
    }
    let (xm, ym) = (xo[i] - xo[i - 1], yo[i] - yo[i - 1]);
    let (xp, yp) = (xo[i + 1] - xo[i], yo[i + 1] - yo[i]);
    let dot = xm * xp + ym * yp;
    let ams = xm * xm + ym * ym;
    let aps = xp * xp + yp * yp;
    let cosa = (dot as f64 / ((ams * aps) as f64).sqrt()) as f32;
    let cosa = if 1.0 < cosa { 1.0 } else { cosa };
    let cosa = if -1.0 > cosa { -1.0 } else { cosa };
    let angle = (cosa as f64).acos() as f32;
    if angle <= maxangle {
        return;
    }
    let inew = if ams > aps { i } else { i + 1 };
    uo.insert(inew, 0.0);
    xo.insert(inew, 0.0);
    yo.insert(inew, 0.0);
    uo[inew] = (0.5 * (uo[inew - 1] + uo[inew + 1]) as f64) as f32;
    intcub(0, ud, xd, &uo[inew..inew + 1], &mut xo[inew..inew + 1]);
    intcub(0, ud, yd, &uo[inew..inew + 1], &mut yo[inew..inew + 1]);
    smoothTwoSegments(maxangle, inew, ud, xd, yd, uo, xo, yo);
    smoothTwoSegments(maxangle, inew - 1, ud, xd, yd, uo, xo, yo);
    smoothTwoSegments(maxangle, inew + 1, ud, xd, yd, uo, xo, yo);
}

unsafe fn makemod(
    xmin: f32,
    zmin: f32,
    xmax: f32,
    zmax: f32,
    ke: &[i32],
    segments: &[Segments],
) -> *mut Model {
    let m = makeModel(zmin, xmin, zmax, xmax);
    (*m).eps = (1.0e-5
        * (((xmax - xmin) as f64).powf(2.0) + ((zmax - zmin) as f64).powf(2.0)).sqrt())
        as f32;
    (*m).sfa = std::mem::size_of::<FaceAttributes>();
    (*m).seua = std::mem::size_of::<EdgeUseAttributes>();
    (*m).sea = std::mem::size_of::<EdgeAttributes>();
    for (ie, a) in segments.iter().enumerate() {
        let mut vlast = newVertex(m, null_mut(), a.x[0], a.y[0], a.s[0]);
        let mut islast = 0;
        for is in 1..a.x.len() {
            let v = newVertex(m, vlast, a.x[is], a.y[is], a.s[is]);
            if v == vlast {
                continue;
            }
            fixEdgeBetweenVertices(vlast, v);
            setEdgeAttributes(vlast, v, ke[ie]);
            setEdgeUseAttributes(
                vlast,
                v,
                a.tx[islast],
                a.ty[islast],
                a.c[islast],
                -a.tx[is],
                -a.ty[is],
                -a.c[is],
            );
            vlast = v;
            islast = is;
        }
    }
    let mut f = (*m).f;
    loop {
        makeSlothForTri(f);
        f = (*f).fNext;
        if f == (*m).f {
            break;
        }
    }
    m
}

unsafe fn newVertex(m: *mut Model, vlast: *mut Vertex, x: f32, z: f32, s: f32) -> *mut Vertex {
    let mut v = addVertexToModel(m, z, x);
    if v.is_null() {
        v = nearestVertexInModel(m, vlast, z, x);
    }
    if (*v).va.is_null() {
        (*v).va = Box::into_raw(Box::new(VertexAttributes { s })).cast();
    } else {
        (*(*v).va.cast::<VertexAttributes>()).s = s;
    }
    v
}

unsafe fn setEdgeAttributes(v1: *mut Vertex, v2: *mut Vertex, k: i32) {
    let mut vu = (*v1).vu;
    let mut eu;
    loop {
        eu = (*vu).eu;
        if (*(*(*eu).euMate).vu).v == v2 {
            break;
        }
        vu = (*vu).vuNext;
        if vu == (*v1).vu {
            break;
        }
    }
    if (*(*(*eu).euMate).vu).v != v2 {
        return;
    }
    (*(*eu).e).ea = Box::into_raw(Box::new(EdgeAttributes { k })).cast();
}

unsafe fn setEdgeUseAttributes(
    v1: *mut Vertex,
    v2: *mut Vertex,
    tx1: f32,
    tz1: f32,
    c1: f32,
    tx2: f32,
    tz2: f32,
    c2: f32,
) {
    let mut vu = (*v1).vu;
    let mut eu;
    loop {
        eu = (*vu).eu;
        if (*(*(*eu).euMate).vu).v == v2 {
            break;
        }
        vu = (*vu).vuNext;
        if vu == (*v1).vu {
            break;
        }
    }
    if (*(*(*eu).euMate).vu).v != v2 {
        return;
    }
    (*eu).eua = Box::into_raw(Box::new(EdgeUseAttributes {
        tx: tx1,
        tz: tz1,
        c: c1,
    }))
    .cast();
    eu = (*eu).euMate;
    (*eu).eua = Box::into_raw(Box::new(EdgeUseAttributes {
        tx: tx2,
        tz: tz2,
        c: c2,
    }))
    .cast();
}

unsafe fn makeSlothForTri(t: *mut Tri) {
    let eu = (*t).eu;
    let v1 = (*(*eu).vu).v;
    let v2 = (*(*(*eu).euCW).vu).v;
    let v3 = (*(*(*eu).euCCW).vu).v;
    if (*v1).va.is_null() || (*v2).va.is_null() || (*v3).va.is_null() {
        segfault();
    }
    let (x1, z1, s1) = (
        (*v1).y as f64,
        (*v1).x as f64,
        (*(*v1).va.cast::<VertexAttributes>()).s as f64,
    );
    let (x2, z2, s2) = (
        (*v2).y as f64,
        (*v2).x as f64,
        (*(*v2).va.cast::<VertexAttributes>()).s as f64,
    );
    let (x3, z3, s3) = (
        (*v3).y as f64,
        (*v3).x as f64,
        (*(*v3).va.cast::<VertexAttributes>()).s as f64,
    );
    let (x2mx1, z2mz1, s2ms1) = (x2 - x1, z2 - z1, s2 - s1);
    let (x3mx1, z3mz1, s3ms1) = (x3 - x1, z3 - z1, s3 - s1);
    let det = z3mz1 * x2mx1 - z2mz1 * x3mx1;
    let a = z3mz1 * s2ms1 - z2mz1 * s3ms1;
    let b = s3ms1 * x2mx1 - s2ms1 * x3mx1;
    let dsdx = a / det;
    let dsdz = b / det;
    let s00 = s2 - dsdx * x2 - dsdz * z2;
    (*t).fa = Box::into_raw(Box::new(FaceAttributes {
        s00: s00 as f32,
        dsdx: dsdx as f32,
        dsdz: dsdz as f32,
        dens: f32::MAX,
        qfac: f32::MAX,
    }))
    .cast();
}

unsafe fn fillsloth(m: *mut Model, sr: &[Vec<f32>]) {
    for r in sr {
        let (x, z, x0, z0, mut s00, dsdx, dsdz) = (r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
        s00 -= x0 * dsdx + z0 * dsdz;
        let t = insideTriInModel(m, null_mut(), z, x);
        setsloth(t, s00, dsdx, dsdz);
    }
}

unsafe fn setsloth(t: *mut Tri, s00: f32, dsdx: f32, dsdz: f32) {
    let fa = (*t).fa.cast::<FaceAttributes>();
    if !fa.is_null() && (*fa).s00 == s00 && (*fa).dsdx == dsdx && (*fa).dsdz == dsdz {
        return;
    }
    if fa.is_null() {
        (*t).fa = Box::into_raw(Box::new(FaceAttributes {
            s00,
            dsdx,
            dsdz,
            dens: 0.0,
            qfac: 0.0,
        }))
        .cast();
    } else {
        (*fa).s00 = s00;
        (*fa).dsdx = dsdx;
        (*fa).dsdz = dsdz;
    }
    let mut eu = (*t).eu;
    loop {
        if (*(*eu).e).fixed == 0 {
            setsloth((*(*eu).euMate).f, s00, dsdx, dsdz);
        }
        eu = (*eu).euCW;
        if eu == (*t).eu {
            break;
        }
    }
}

unsafe fn filldens(m: *mut Model, dptr: &[Vec<f32>]) {
    for r in dptr {
        let t = insideTriInModel(m, null_mut(), r[1], r[0]);
        setdens(t, r[2]);
    }
    let mut f = (*m).f;
    loop {
        let dens = (*(*f).fa.cast::<FaceAttributes>()).dens;
        let mut filled = false;
        for r in dptr {
            if dens == r[2] {
                filled = true;
                break;
            }
        }
        if !filled {
            err!("WARNING!!\nCHECK THAT ALL REGIONS ARE FILLED WITH A DENSITY VALUE\n");
        }
        f = (*f).fNext;
        if f == (*m).f {
            break;
        }
    }
}

unsafe fn setdens(t: *mut Tri, dens: f32) {
    let fa = (*t).fa.cast::<FaceAttributes>();
    if !fa.is_null() && (*fa).dens == dens {
        return;
    }
    if fa.is_null() {
        (*t).fa = Box::into_raw(Box::new(FaceAttributes {
            s00: 0.0,
            dsdx: 0.0,
            dsdz: 0.0,
            dens,
            qfac: 0.0,
        }))
        .cast();
    } else {
        (*fa).dens = dens;
    }
    let mut eu = (*t).eu;
    loop {
        if (*(*eu).e).fixed == 0 {
            setdens((*(*eu).euMate).f, dens);
        }
        eu = (*eu).euCW;
        if eu == (*t).eu {
            break;
        }
    }
}

unsafe fn fillq(m: *mut Model, qptr: &[Vec<f32>]) {
    for r in qptr {
        let t = insideTriInModel(m, null_mut(), r[1], r[0]);
        setq(t, r[2]);
    }
    let mut f = (*m).f;
    loop {
        let qfac = (*(*f).fa.cast::<FaceAttributes>()).qfac;
        let mut filled = false;
        for r in qptr {
            if qfac == r[2] {
                filled = true;
                break;
            }
        }
        if !filled {
            err!("WARNING!!\nCHECK THAT ALL REGIONS ARE FILLED WITH A Q-FACTOR\n");
        }
        f = (*f).fNext;
        if f == (*m).f {
            break;
        }
    }
}

unsafe fn setq(t: *mut Tri, qfac: f32) {
    let fa = (*t).fa.cast::<FaceAttributes>();
    if !fa.is_null() && (*fa).qfac == qfac {
        return;
    }
    if fa.is_null() {
        (*t).fa = Box::into_raw(Box::new(FaceAttributes {
            s00: 0.0,
            dsdx: 0.0,
            dsdz: 0.0,
            dens: 0.0,
            qfac,
        }))
        .cast();
    } else {
        (*fa).qfac = qfac;
    }
    let mut eu = (*t).eu;
    loop {
        if (*(*eu).e).fixed == 0 {
            setq((*(*eu).euMate).f, qfac);
        }
        eu = (*eu).euCW;
        if eu == (*t).eu {
            break;
        }
    }
}
