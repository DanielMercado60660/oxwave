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

//! Port of tri/lib/{check,circum,colinear,create,delete,distance,fixedges,
//! inside,nearest,project,readwrite}.c.
//! Credits: Dave Hale, Colorado School of Mines, 1990–1991;
//! Chris Elmer and Dave Hale (unfixEdge); Jack K. Cohen (model I/O),
//! with modifications by Dave Hale and Craig Artley, 1990–1994.

#![allow(non_snake_case)]

use std::ffi::c_void;

mod check;
mod create;
mod delete;
mod fixedges;
mod readwrite;
pub use check::*;
pub use create::*;
pub use delete::*;
pub use fixedges::*;
pub use readwrite::*;

extern "C" {
    fn calloc(count: usize, size: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    fn exit(status: i32) -> !;
}

// Use the C allocator so raw topology has the same ownership convention as SU.
// Zero padding and fields left unspecified by C (notably Model.ma).
unsafe fn allocate<T>() -> *mut T {
    let p: *mut T = unsafe { calloc(1, std::mem::size_of::<T>()) }.cast();
    assert!(!p.is_null(), "tri: allocation failed");
    p
}

/// C `Vertex` layout from tri/include/tri.h.
#[repr(C)]
pub struct Vertex {
    pub vu: *mut VertexUse,
    pub fixed: i32,
    pub x: f32,
    pub y: f32,
    pub va: *mut c_void,
}

/// C: `void circum(float x1, float y1, float x2, float y2, float x3, float y3, float *xc, float *yc, float *rs);`
/// Port of tri/lib/circum.c; author Dave Hale, Colorado School of Mines, 1990.
pub fn circum(
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    x3: f32,
    y3: f32,
    xc: &mut f32,
    yc: &mut f32,
    rs: &mut f32,
) {
    let x21 = x2 - x1;
    let y21 = y2 - y1;
    let x31 = x3 - x1;
    let y31 = y3 - y1;
    let r21 = x21 * x21 + y21 * y21;
    let r31 = x31 * x31 + y31 * y31;
    let mut det = x21 * y31 - x31 * y21;
    let tiny = (1e-10_f64 * (r21 + r31) as f64) as f32;
    if det < 0.0 && -det < tiny {
        det = -tiny;
    } else if det > 0.0 && det < tiny {
        det = tiny;
    }
    if det != 0.0 {
        let scl = (0.5_f64 / det as f64) as f32;
        let xcr = scl * (r21 * y31 - r31 * y21);
        let ycr = scl * (x21 * r31 - x31 * r21);
        *rs = xcr * xcr + ycr * ycr;
        *xc = x1 + xcr;
        *yc = y1 + ycr;
    } else {
        *rs = f32::MAX;
        *xc = 0.0;
        *yc = 0.0;
    }
}

/// C: `int inCircum(float x, float y, float xc, float yc, float rs);`
/// Port of tri/lib/inside.c; author Dave Hale, Colorado School of Mines, 1991.
pub fn inCircum(x: f32, y: f32, xc: f32, yc: f32, rs: f32) -> i32 {
    if rs == f32::MAX {
        return 1;
    }
    let mut d = xc - x;
    let mut ds = d * d;
    if ds >= rs {
        return 0;
    }
    d = yc - y;
    ds += d * d;
    if ds >= rs {
        return 0;
    }
    1
}

/// C: `int in3Vertices(float x, float y, Vertex *v1, Vertex *v2, Vertex *v3);`
/// Port of tri/lib/inside.c; author Dave Hale, Colorado School of Mines, 1991.
pub fn in3Vertices(x: f32, y: f32, v1: &Vertex, v2: &Vertex, v3: &Vertex) -> i32 {
    let (x1, y1, x2, y2, x3, y3) = (v1.x, v1.y, v2.x, v2.y, v3.x, v3.y);
    let s1 = (x - x1) * (y2 - y1) - (x2 - x1) * (y - y1);
    let s2 = (x - x2) * (y3 - y2) - (x3 - x2) * (y - y2);
    if s1 * s2 <= 0.0 {
        return 0;
    }
    let s3 = (x - x3) * (y1 - y3) - (x1 - x3) * (y - y3);
    if s2 * s3 <= 0.0 || s1 * s3 <= 0.0 {
        return 0;
    }
    1
}

/// C: `int vertexBetweenVertices(Vertex *v, Vertex *v1, Vertex *v2);`
/// Port of tri/lib/colinear.c; author Dave Hale, Center for Wave Phenomena, 1990–1991.
pub fn vertexBetweenVertices(v: &Vertex, v1: &Vertex, v2: &Vertex) -> i32 {
    let x1 = v.x - v1.x;
    let y1 = v.y - v1.y;
    let x2 = v.x - v2.x;
    let y2 = v.y - v2.y;
    let x12 = v2.x - v1.x;
    let y12 = v2.y - v1.y;
    let a1 = x1 * x1 + y1 * y1;
    let a2 = x2 * x2 + y2 * y2;
    let a12 = x12 * x12 + y12 * y12;
    if a1 == 0.0 || a2 == 0.0 {
        return 1;
    }
    if a12 == 0.0 {
        return 0;
    }
    let d1 = x1 * x12 + y1 * y12;
    let d2 = x2 * x12 + y2 * y12;
    if d1 * d2 > 0.0 {
        return 0;
    }
    let c = x1 * y12 - x12 * y1;
    let s1 = (c * c) / (a1 * a12);
    if s1 as f64 > 0.0001 {
        return 0;
    }
    let c = x2 * y12 - x12 * y2;
    let s2 = (c * c) / (a2 * a12);
    if s2 as f64 > 0.0001 {
        return 0;
    }
    1
}

/// C `VertexUse` from tri/include/tri.h.
#[repr(C)]
pub struct VertexUse {
    pub eu: *mut EdgeUse,
    pub vuPrev: *mut VertexUse,
    pub vuNext: *mut VertexUse,
    pub v: *mut Vertex,
    pub vua: *mut c_void,
}
/// C `Edge` from tri/include/tri.h.
#[repr(C)]
pub struct Edge {
    pub eu: *mut EdgeUse,
    pub fixed: i32,
    pub ea: *mut c_void,
}
/// C `EdgeUse` from tri/include/tri.h.
#[repr(C)]
pub struct EdgeUse {
    pub f: *mut Face,
    pub vu: *mut VertexUse,
    pub euMate: *mut EdgeUse,
    pub euCW: *mut EdgeUse,
    pub euCCW: *mut EdgeUse,
    pub e: *mut Edge,
    pub eua: *mut c_void,
}
/// C `Face` (also called `Tri`) from tri/include/tri.h.
#[repr(C)]
pub struct Face {
    pub m: *mut Model,
    pub fPrev: *mut Face,
    pub fNext: *mut Face,
    pub eu: *mut EdgeUse,
    pub xc: f32,
    pub yc: f32,
    pub rs: f32,
    pub fa: *mut c_void,
}
/// C `Tri` is an alias of `Face`.
pub type Tri = Face;
/// C `Model` from tri/include/tri.h.
#[repr(C)]
pub struct Model {
    pub f: *mut Face,
    pub xmin: f32,
    pub ymin: f32,
    pub xmax: f32,
    pub ymax: f32,
    pub eps: f32,
    pub sma: usize,
    pub sfa: usize,
    pub seua: usize,
    pub sea: usize,
    pub svua: usize,
    pub sva: usize,
    pub vAdd: Option<unsafe extern "C" fn(*mut Model, *mut Vertex)>,
    pub vDel: Option<unsafe extern "C" fn(*mut Model, *mut Vertex)>,
    pub tAdd: Option<unsafe extern "C" fn(*mut Model, *mut Tri)>,
    pub tDel: Option<unsafe extern "C" fn(*mut Model, *mut Tri)>,
    pub ma: *mut c_void,
}

unsafe fn endpoints(e: &Edge) -> (&Vertex, &Vertex) {
    unsafe { (&*(*(*e.eu).vu).v, &*(*(*(*e.eu).euCW).vu).v) }
}
unsafe fn triangle_vertices(t: &Tri) -> (&Vertex, &Vertex, &Vertex) {
    unsafe {
        let eu = &*t.eu;
        (&*(*eu.vu).v, &*(*(*eu.euCW).vu).v, &*(*(*eu.euCCW).vu).v)
    }
}

/// C: `void circumTri (Tri *t);` (tri/lib/circum.c, Dave Hale).
/// # Safety
/// The face and its three edge and vertex links must be valid.
pub unsafe fn circumTri(t: &mut Tri) {
    let (a, b, c) = unsafe { triangle_vertices(t) };
    let (mut xc, mut yc, mut rs) = (0.0, 0.0, 0.0);
    circum(a.x, a.y, b.x, b.y, c.x, c.y, &mut xc, &mut yc, &mut rs);
    t.rs = rs;
    t.xc = xc;
    t.yc = yc;
}
/// C: `int inCircumTri (float x, float y, Tri *t);` (tri/lib/inside.c, Dave Hale).
pub fn inCircumTri(x: f32, y: f32, t: &Tri) -> i32 {
    inCircum(x, y, t.xc, t.yc, t.rs)
}
/// C: `int inTri (float x, float y, Tri *t);` (tri/lib/inside.c, Dave Hale).
/// # Safety
/// The face and its three edge and vertex links must be valid.
pub unsafe fn inTri(x: f32, y: f32, t: &Tri) -> i32 {
    let (a, b, c) = unsafe { triangle_vertices(t) };
    in3Vertices(x, y, a, b, c)
}
/// C: `int edgesColinear (Edge *e1, Edge *e2);` (tri/lib/colinear.c, Dave Hale).
/// # Safety
/// Both edges must have valid endpoint links.
pub unsafe fn edgesColinear(e1: &Edge, e2: &Edge) -> i32 {
    let (a, b) = unsafe { endpoints(e1) };
    let (c, d) = unsafe { endpoints(e2) };
    let (x1, y1, x2, y2) = (b.x - a.x, b.y - a.y, d.x - c.x, d.y - c.y);
    let cross = x1 * y2 - x2 * y1;
    let a1 = x1 * x1 + y1 * y1;
    let a2 = x2 * x2 + y2 * y2;
    let s = (cross * cross) / (a1 * a2);
    if (s as f64) < 0.0001 {
        1
    } else {
        0
    }
}
/// C: `float distanceToEdge (Edge *e, float x, float y);` (tri/lib/distance.c, Dave Hale).
/// # Safety
/// The edge must have valid endpoint links.
pub unsafe fn distanceToEdge(e: &Edge, x: f32, y: f32) -> f32 {
    let (a, b) = unsafe { endpoints(e) };
    let (x1mx, y1my) = (a.x - x, a.y - y);
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let mut t = -(x1mx * dx + y1my * dy) / (dx * dx + dy * dy);
    if t < 0.0 {
        t = 0.0
    }
    if t > 1.0 {
        t = 1.0
    }
    let xf = x1mx + t * dx;
    let yf = y1my + t * dy;
    ((xf * xf + yf * yf) as f64).sqrt() as f32
}
/// C: `void projectToEdge (Edge *e, float *x, float *y);` (tri/lib/project.c, Dave Hale).
/// # Safety
/// The edge must have valid endpoint links.
pub unsafe fn projectToEdge(e: &Edge, x: &mut f32, y: &mut f32) {
    let (a, b) = unsafe { endpoints(e) };
    let (x1mx, y1my) = (a.x - *x, a.y - *y);
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let mut t = -(x1mx * dx + y1my * dy) / (dx * dx + dy * dy);
    if t < 0.0 {
        t = 0.0
    }
    if t > 1.0 {
        t = 1.0
    }
    *x = a.x + t * dx;
    *y = a.y + t * dy;
}

/// C: `Tri *insideTriInModel (Model *m, Tri *start, float x, float y);` (tri/lib/create.c, Dave Hale).
/// # Safety
/// The model and all traversed topology links must be valid; `start` may be null.
pub unsafe fn insideTriInModel(m: *mut Model, start: *mut Tri, x: f32, y: f32) -> *mut Tri {
    unsafe {
        let mut t = if start.is_null() { (*m).f } else { start };
        loop {
            let mut eu = (*t).eu;
            loop {
                let v1 = &*(*(*eu).vu).v;
                let v2 = &*(*(*(*eu).euCW).vu).v;
                let v3 = &*(*(*(*eu).euCCW).vu).v;
                let s1 = (v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y);
                let s2 = (v2.x - v1.x) * (y - v1.y) - (x - v1.x) * (v2.y - v1.y);
                if s1 * s2 < 0.0 && !(*(*eu).euMate).f.is_null() {
                    t = (*(*eu).euMate).f;
                    break;
                }
                eu = (*eu).euCW;
                if eu == (*t).eu {
                    break;
                }
            }
            if eu == (*t).eu {
                return t;
            }
        }
    }
}
/// C: `Edge *nearestEdgeInModel (Model *m, Edge *start, float x, float y);` (tri/lib/nearest.c, Dave Hale).
/// # Safety
/// The model and all traversed topology links must be valid; `start` may be null.
pub unsafe fn nearestEdgeInModel(m: *mut Model, start: *mut Edge, x: f32, y: f32) -> *mut Edge {
    unsafe {
        let t = if start.is_null() {
            insideTriInModel(m, std::ptr::null_mut(), x, y)
        } else {
            let eu = (*start).eu;
            let f = if !(*eu).f.is_null() {
                (*eu).f
            } else {
                (*(*eu).euMate).f
            };
            insideTriInModel(m, f, x, y)
        };
        let mut eu = (*t).eu;
        let mut dmin = 0.0;
        let mut emin = std::ptr::null_mut();
        loop {
            let d = distanceToEdge(&*(*eu).e, x, y);
            if eu == (*t).eu || d < dmin {
                dmin = d;
                emin = (*eu).e;
            }
            eu = (*eu).euCW;
            if eu == (*t).eu {
                break;
            }
        }
        emin
    }
}
/// C: `Vertex *nearestVertexInModel (Model *m, Vertex *start, float x, float y);` (tri/lib/nearest.c, Dave Hale).
/// # Safety
/// The model and all traversed topology links must be valid; `start` may be null.
pub unsafe fn nearestVertexInModel(
    m: *mut Model,
    start: *mut Vertex,
    x: f32,
    y: f32,
) -> *mut Vertex {
    unsafe {
        let mut vmin = if start.is_null() {
            (*(*(*(*m).f).eu).vu).v
        } else {
            start
        };
        let dx = (*vmin).x - x;
        let dy = (*vmin).y - y;
        let mut dsmin = dx * dx + dy * dy;
        loop {
            let v = vmin;
            let mut vu = (*v).vu as *mut VertexUse;
            loop {
                let vn = (*(*(*(*vu).eu).euCW).vu).v;
                let dx = (*vn).x - x;
                let dy = (*vn).y - y;
                let ds = dx * dx + dy * dy;
                if dsmin > ds {
                    dsmin = ds;
                    vmin = vn;
                }
                vu = (*vu).vuNext;
                if vu == (*v).vu as *mut VertexUse {
                    break;
                }
            }
            if v == vmin {
                return v;
            }
        }
    }
}
