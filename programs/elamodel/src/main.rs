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

//! ELAMODEL - make piecewise homogeneous anisotropic model.
//! Port of Trielas/main/elamodel.c (Test Release 1.1).
//! AUTHOR: Dave Hale, Colorado School of Mines, 02/12/91.
//! Modified: Andreas Rueger, Colorado School of Mines, 01/18/94;
//! built anisotropic models.

#![allow(non_snake_case)]

use oxwave_core::prelude::*;
use oxwave_cubic_interp::*;
use oxwave_f_cubicspline::*;
use oxwave_lib_tri::*;
use std::ffi::{c_char, c_void};
use std::ptr::null_mut;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " ELAMODEL - make piecewise homogeneous anisotropic model    \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " elamodel >modelfile fill= [optional parameters]   \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Input Parameters: \t\t\t\t\t\t\t",
    " xmin=0.0               minimum horizontal coordinate (x) \t\t",
    " xmax=1.0               maximum horizontal coordinate (x) \t\t",
    " zmin=0.0               minimum vertical coordinate (z) \t\t",
    " zmax=1.0               maximum vertical coordinate (z) \t\t",
    " xedge=                 x coordinates of an edge \t\t\t",
    " zedge=                 z coordinates of an edge \t\t\t",
    " kedge=                 array of indices used to identify edges \t",
    " fill=    iso      \t x,z,v_p,v_s,rho   \t\t\t\t",
    "          tiso      \t x,z,v_p,v_s,epsilon,delta,gamma,phi,rho\t",
    "          ani           x,z,a1111,a3333,a1133,a1313,a1113,a3313        ",
    "                            a1212,a2323,a1223,rho                      ",
    " maxangle=5.0           maximum angle (deg) between adjacent edge segments ",
    "   \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " More than set of xedge and zedge parameters may be \t\t        ",
    " specified, but the numbers of these parameters must be equal. \t",
    " \t\t\t\t\t\t\t\t\t",
    " Within each set, vertices will be connected by fixed edges. \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Edge indices in the k array are used to identify edges \t\t",
    " specified by the x and z parameters.  The first k index \t\t",
    " corresponds to the first set of x and z parameters, the \t\t",
    " second k index corresponds to the second set, and so on. \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " After all vertices have been inserted into the model,\tthe fill        ",
    " parameters is used to fill closed regions bounded by fixed edges.     ",
    " Three input modes are available:                                      ",
    " Isotropic blocks:     x,z,v_p,v_s,rho                                 ",
    " Transversely iso:     x,z,v_p,v_s,epsilon,delta,gamma,phi,rho         ",
    " General 2D aniso:     x,z,a1111,a3333,a1133,a1313,a1113,a3313         ",
    "                       a1212,a2323,a1223,rho                           ",
    " Hereby:  ",
    " x,z\t\t\t   coordinates of one point in a bounded region ",
    " v_p,v_s\t\t   P, S-wave velocity along symmetry axis       ",
    " epsilon, delta, gammma   Thomsen's parameters              ",
    " rho \t\t\t   density \t\t\t     ",
    " phi\t\t\t   angle of symmetry axes with vertical ",
    " aijkl\t\t\t   density normalized stiffness coefficients ",
    "  ",
    " Each block can be defined by different input modes. The number of     ",
    " input parameters defines the input type. Incorrect number of input    ",
    " parameters result in an Error-message\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

// The tri library uses native C FILE streams for its native model format.
extern "C" {
    fn fdopen(fd: i32, mode: *const c_char) -> *mut c_void;
    fn fflush(fp: *mut c_void) -> i32;
}

#[repr(C)]
#[derive(Clone, Copy)]
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

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    let (mut xmin, mut xmax, mut zmin, mut zmax, mut maxangle) = (0.0, 1.0, 0.0, 1.0, 5.0f32);
    if getparfloat("xmin", &mut xmin) == 0 {
        xmin = 0.0;
    }
    if getparfloat("xmax", &mut xmax) == 0 {
        xmax = 1.0;
    }
    if getparfloat("zmin", &mut zmin) == 0 {
        zmin = 0.0;
    }
    if getparfloat("zmax", &mut zmax) == 0 {
        zmax = 1.0;
    }
    if getparfloat("maxangle", &mut maxangle) == 0 {
        maxangle = 5.0;
    }
    maxangle = (maxangle as f64 * (std::f64::consts::PI / 180.0)) as f32;

    let nxe = countparname("xedge");
    let nze = countparname("zedge");
    if nxe != nze {
        err!("\n ERROR number of xedge must equal number of zedge\n");
    }
    let ne = nxe as usize;
    let mut xe = vec![Vec::<f32>::new(); ne];
    let mut ze = vec![Vec::<f32>::new(); ne];
    for ie in 0..ne {
        let nxe = countnparval(ie as i32 + 1, "xedge");
        let nze = countnparval(ie as i32 + 1, "zedge");
        if nxe != nze {
            err!("\n ERROR number of xedge must equal number of zedge\n");
        }
        getnparfloat(ie as i32 + 1, "xedge", &mut xe[ie]);
        getnparfloat(ie as i32 + 1, "zedge", &mut ze[ie]);
    }
    let nke = countparval("kedge") as usize;
    let mut ke = vec![0i32; ne];
    if nke > ne {
        err!("\n ERROR more <kedge> values than edges specified\n");
    }
    if nke < ne {
        err!("\n ERROR less <kedge> values than edges specified\n");
    }
    getparint("kedge", &mut ke);
    checkpars();

    let nr = countparname("fill");
    if nr == 0 {
        err!("\n ERROR required parameters <fill> missing \n");
    }
    let mut sr = vec![0.0f32; 12];
    let mut sra = vec![[0.0f32; 13]; nr as usize];
    for ir in 0..nr as usize {
        let num = getnparfloat(ir as i32 + 1, "fill", &mut sr);
        if num != 5 && num != 9 && num != 12 {
            err!(
                "\nERROR wrong number of input parameters in <fill> #%i\n",
                ir as i32 + 1
            );
        }
        let s = &mut sra[ir];
        if num == 5 {
            s[0] = sr[0];
            s[1] = sr[1];
            s[3] = sr[2] * sr[2];
            s[2] = s[3];
            s[5] = sr[3] * sr[3];
            s[4] = s[2] - 2.0f32 * s[5];
            s[6] = 0.0;
            s[7] = 0.0;
            s[9] = s[5];
            s[8] = s[9];
            s[10] = 0.0;
            s[11] = sr[4];
            s[12] = 0.0;
        } else if num == 9 {
            s[0] = sr[0];
            s[1] = sr[1];
            s[3] = sr[2] * sr[2];
            s[5] = sr[3] * sr[3];
            s[2] = s[3] + 2.0f32 * sr[4] * s[3];
            let temp = 2.0f32 * sr[5] * s[3] * (s[3] - s[5]) + (s[3] - s[5]) * (s[3] - s[5]);
            if temp < 0.0 {
                err!("ERROR: a1133 not reel \n");
            }
            s[4] = ((temp as f64).sqrt() - s[5] as f64) as f32;
            s[6] = 0.0;
            s[7] = 0.0;
            s[10] = 0.0;
            s[8] = 2.0f32 * s[5] * sr[6] + s[5];
            s[9] = s[5];
            s[11] = sr[8];
            s[12] = 1.0;
            let angle = (sr[7] as f64 * std::f64::consts::PI / 180.0) as f32;
            if angle != 0.0 && rottensors(s, angle) != 1 {
                err!(" Error \n rotation failed \n");
            }
        } else {
            for i in 0..12 {
                s[i] = sr[i];
            }
            s[12] = 2.0;
        }
    }

    let mut edges = Vec::new();
    for ie in 0..ne {
        edges.push(smoothLineSegments(maxangle, &xe[ie], &ze[ie]));
    }
    unsafe {
        let m = makemod(xmin, zmin, xmax, zmax, &ke, &edges);
        fillcoeff(m, &sra);
        checkfill(m, &sra);
        let fp = fdopen(1, c"wb".as_ptr());
        writeModel(m, fp);
        fflush(fp);
    }
    // The original returns 1 even after successfully writing its native model.
    exit(1);
}

struct Segments {
    x: Vec<f32>,
    y: Vec<f32>,
    tx: Vec<f32>,
    ty: Vec<f32>,
    c: Vec<f32>,
}

fn smoothLineSegments(maxangle: f32, xin: &[f32], yin: &[f32]) -> Segments {
    let nin = xin.len();
    if nin == 1 {
        return Segments {
            x: xin.to_vec(),
            y: yin.to_vec(),
            tx: vec![],
            ty: vec![],
            c: vec![],
        };
    }
    let mut uin = vec![0.0f32; nin];
    let mut u = vec![0.0f32; nin];
    let mut x = xin.to_vec();
    let mut y = yin.to_vec();
    uin[0] = 0.0;
    u[0] = 0.0;
    for i in 1..nin {
        let v = (u[i - 1] as f64
            + (((x[i] - x[i - 1]) as f64).powf(2.0) + ((y[i] - y[i - 1]) as f64).powf(2.0)).sqrt())
            as f32;
        uin[i] = v;
        u[i] = v;
    }
    let mut xind = vec![[0.0f32; 4]; nin];
    let mut yind = vec![[0.0f32; 4]; nin];
    csplin(&uin, xin, &mut xind);
    csplin(&uin, yin, &mut yind);
    let mut i = 1;
    while i + 1 < x.len() {
        smoothTwoSegments(maxangle, i, &uin, &xind, &yind, &mut u, &mut x, &mut y);
        i += 1;
    }
    let n = x.len();
    let mut tx = vec![0.0; n];
    let mut ty = vec![0.0; n];
    let mut c = vec![0.0; n];
    for i in 0..n {
        let (mut dx, mut dy, mut ddx, mut ddy) = ([0.0f32], [0.0f32], [0.0f32], [0.0f32]);
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
    Segments { x, y, tx, ty, c }
}

#[allow(clippy::too_many_arguments)]
fn smoothTwoSegments(
    maxangle: f32,
    i: usize,
    ud: &[f32],
    xd: &[[f32; 4]],
    yd: &[[f32; 4]],
    u: &mut Vec<f32>,
    x: &mut Vec<f32>,
    y: &mut Vec<f32>,
) {
    let no = x.len();
    if i == 0 || i == no - 1 {
        return;
    }
    let xm = x[i] - x[i - 1];
    let ym = y[i] - y[i - 1];
    let xp = x[i + 1] - x[i];
    let yp = y[i + 1] - y[i];
    let dot = xm * xp + ym * yp;
    let ams = xm * xm + ym * ym;
    let aps = xp * xp + yp * yp;
    let cosa = (dot as f64 / ((ams * aps) as f64).sqrt()) as f32;
    let cosa = if 1.0f64 < cosa as f64 {
        1.0
    } else {
        cosa as f64
    };
    let cosa = if -1.0f64 > cosa { -1.0 } else { cosa };
    let angle = cosa.acos() as f32;
    if angle <= maxangle {
        return;
    }
    let inew = if ams > aps { i } else { i + 1 };
    u.insert(inew, 0.0);
    x.insert(inew, 0.0);
    y.insert(inew, 0.0);
    u[inew] = (0.5 * (u[inew - 1] + u[inew + 1]) as f64) as f32;
    intcub(0, ud, xd, &u[inew..inew + 1], &mut x[inew..inew + 1]);
    intcub(0, ud, yd, &u[inew..inew + 1], &mut y[inew..inew + 1]);
    smoothTwoSegments(maxangle, inew, ud, xd, yd, u, x, y);
    smoothTwoSegments(maxangle, inew - 1, ud, xd, yd, u, x, y);
    smoothTwoSegments(maxangle, inew + 1, ud, xd, yd, u, x, y);
}

unsafe fn makemod(
    xmin: f32,
    zmin: f32,
    xmax: f32,
    zmax: f32,
    ke: &[i32],
    edges: &[Segments],
) -> *mut Model {
    unsafe {
        let m = makeModel(zmin, xmin, zmax, xmax);
        (*m).eps = (1.0e-5
            * (((xmax - xmin) as f64).powf(2.0) + ((zmax - zmin) as f64).powf(2.0)).sqrt())
            as f32;
        (*m).sfa = std::mem::size_of::<FaceAttributes>();
        (*m).seua = std::mem::size_of::<EdgeUseAttributes>();
        (*m).sea = std::mem::size_of::<EdgeAttributes>();
        for ie in 0..edges.len() {
            let e = &edges[ie];
            let mut vlast = newVertex(m, null_mut(), e.x[0], e.y[0]);
            let mut islast = 0;
            for is in 1..e.x.len() {
                let v = newVertex(m, vlast, e.x[is], e.y[is]);
                if v == vlast {
                    continue;
                }
                fixEdgeBetweenVertices(vlast, v);
                setEdgeAttributes(vlast, v, ke[ie]);
                setEdgeUseAttributes(
                    vlast,
                    v,
                    e.tx[islast],
                    e.ty[islast],
                    e.c[islast],
                    -e.tx[is],
                    -e.ty[is],
                    -e.c[is],
                );
                vlast = v;
                islast = is;
            }
        }
        m
    }
}

unsafe fn newVertex(m: *mut Model, vlast: *mut Vertex, x: f32, z: f32) -> *mut Vertex {
    unsafe {
        let mut v = addVertexToModel(m, z, x);
        if v.is_null() {
            v = nearestVertexInModel(m, vlast, z, x);
        }
        v
    }
}

unsafe fn setEdgeAttributes(v1: *mut Vertex, v2: *mut Vertex, k: i32) {
    unsafe {
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
}

#[allow(clippy::too_many_arguments)]
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
    unsafe {
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
}

unsafe fn fillcoeff(m: *mut Model, sr: &[[f32; 13]]) {
    unsafe {
        for s in sr {
            let t = insideTriInModel(m, null_mut(), s[1], s[0]);
            setcoeff(
                t,
                FaceAttributes {
                    a1111: s[2],
                    a3333: s[3],
                    a1133: s[4],
                    a1313: s[5],
                    a1113: s[6],
                    a3313: s[7],
                    a1212: s[8],
                    a2323: s[9],
                    a1223: s[10],
                    rho: s[11],
                    mindex: s[12] as i32,
                },
            );
        }
    }
}

unsafe fn setcoeff(t: *mut Tri, a: FaceAttributes) {
    unsafe {
        let fa = (*t).fa.cast::<FaceAttributes>();
        if !fa.is_null() {
            let f = &*fa;
            if f.a1111 == a.a1111
                && f.a3333 == a.a3333
                && f.a1133 == a.a1133
                && f.a1313 == a.a1313
                && f.a3313 == a.a3313
                && f.a1113 == a.a1113
                && f.a1212 == a.a1212
                && f.a2323 == a.a2323
                && f.a1223 == a.a1223
                && f.rho == a.rho
            {
                return;
            }
        }
        if fa.is_null() {
            (*t).fa = Box::into_raw(Box::new(a)).cast();
        } else {
            *fa = a;
        }
        let mut eu = (*t).eu;
        loop {
            if (*(*eu).e).fixed == 0 {
                setcoeff((*(*eu).euMate).f, a);
            }
            eu = (*eu).euCW;
            if eu == (*t).eu {
                break;
            }
        }
    }
}

unsafe fn checkfill(m: *mut Model, sra: &[[f32; 13]]) {
    unsafe {
        let mut f = (*m).f;
        loop {
            if (*f).fa.is_null() {
                err!("\n ERROR: You didn't fill all regions \n");
            }
            let fac = (*(*f).fa.cast::<FaceAttributes>()).a1111;
            let mut icheck = false;
            for s in sra {
                if fac == s[2] {
                    icheck = true;
                    break;
                }
            }
            if !icheck {
                err!("\n ERROR: You didn't fill all regions \n");
            }
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
    }
}

fn rottensors(aijkl: &mut [f32; 13], angle: f32) -> i32 {
    let si = (angle as f64).sin() as f32;
    let co = (angle as f64).cos() as f32;
    let ss = si * si;
    let cc = co * co;
    let sc = si * co;
    let a = aijkl[2];
    let c = aijkl[3];
    let f = aijkl[4];
    let l = aijkl[5];
    let m = aijkl[6];
    let n = aijkl[7];
    let o = aijkl[8];
    let p = aijkl[9];
    let q = aijkl[10];
    aijkl[2] = a * cc * cc + c * ss * ss + 2.0 * f * cc * ss + 4.0 * l * cc * ss
        - 4.0 * m * cc * sc
        - 4.0 * n * ss * sc;
    aijkl[3] = a * ss * ss
        + c * cc * cc
        + 2.0 * f * cc * ss
        + 4.0 * l * cc * ss
        + 4.0 * m * ss * sc
        + 4.0 * n * cc * sc;
    aijkl[4] = a * ss * cc + c * ss * cc + f * (ss * ss + cc * cc)
        - 4.0 * l * ss * cc
        - 2.0 * m * (ss * sc - cc * sc)
        - 2.0 * n * (cc * sc - ss * sc);
    aijkl[5] = a * ss * cc + c * ss * cc - 2.0 * f * ss * cc
        + l * (cc * cc + ss * ss - 2.0 * ss * cc)
        - 2.0 * m * (ss * sc - cc * sc)
        - 2.0 * n * (cc * sc - ss * sc);
    aijkl[6] = a * cc * sc - c * ss * sc
        + f * (ss * sc - cc * sc)
        + 2.0 * l * (ss * sc - cc * sc)
        + m * (cc * cc - 3.0 * cc * ss)
        + (3.0 * cc * ss - ss * ss) * n;
    aijkl[7] = a * ss * sc - c * cc * sc
        + f * (cc * sc - ss * sc)
        + 2.0 * l * (cc * sc - ss * sc)
        + m * (3.0 * cc * ss - ss * ss)
        + n * (cc * cc - 3.0 * ss * cc);
    aijkl[8] = cc * o - 2.0 * sc * q + ss * p;
    aijkl[9] = ss * o + 2.0 * sc * q + cc * p;
    aijkl[10] = sc * o + (cc - ss) * q - sc * p;
    1
}
