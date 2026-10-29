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

//! Port of tri/lib/create.c.
//! Author: Dave Hale, Colorado School of Mines, Fall 1990.

#![allow(unused_assignments, unused_parens, unused_mut)]
use super::*;

/// C: `Model *makeModel (float xmin, float ymin, float xmax, float ymax);`
/// Returns a caller-owned model containing two triangles in a rectangle.
pub fn makeModel(mut xmin: f32, mut ymin: f32, mut xmax: f32, mut ymax: f32) -> *mut Model {
    unsafe {
        let mut v1: *mut Vertex = std::ptr::null_mut();
        let mut v2: *mut Vertex = std::ptr::null_mut();
        let mut v3: *mut Vertex = std::ptr::null_mut();
        let mut v4: *mut Vertex = std::ptr::null_mut();
        let mut vu12: *mut VertexUse = std::ptr::null_mut();
        let mut vu14: *mut VertexUse = std::ptr::null_mut();
        let mut vu21: *mut VertexUse = std::ptr::null_mut();
        let mut vu23: *mut VertexUse = std::ptr::null_mut();
        let mut vu24: *mut VertexUse = std::ptr::null_mut();
        let mut vu32: *mut VertexUse = std::ptr::null_mut();
        let mut vu34: *mut VertexUse = std::ptr::null_mut();
        let mut vu41: *mut VertexUse = std::ptr::null_mut();
        let mut vu42: *mut VertexUse = std::ptr::null_mut();
        let mut vu43: *mut VertexUse = std::ptr::null_mut();
        let mut e1: *mut Edge = std::ptr::null_mut();
        let mut e2: *mut Edge = std::ptr::null_mut();
        let mut e3: *mut Edge = std::ptr::null_mut();
        let mut e4: *mut Edge = std::ptr::null_mut();
        let mut e5: *mut Edge = std::ptr::null_mut();
        let mut eu12: *mut EdgeUse = std::ptr::null_mut();
        let mut eu14: *mut EdgeUse = std::ptr::null_mut();
        let mut eu21: *mut EdgeUse = std::ptr::null_mut();
        let mut eu23: *mut EdgeUse = std::ptr::null_mut();
        let mut eu24: *mut EdgeUse = std::ptr::null_mut();
        let mut eu32: *mut EdgeUse = std::ptr::null_mut();
        let mut eu34: *mut EdgeUse = std::ptr::null_mut();
        let mut eu41: *mut EdgeUse = std::ptr::null_mut();
        let mut eu42: *mut EdgeUse = std::ptr::null_mut();
        let mut eu43: *mut EdgeUse = std::ptr::null_mut();
        let mut f124: *mut Face = std::ptr::null_mut();
        let mut f234: *mut Face = std::ptr::null_mut();
        let mut m: *mut Model = std::ptr::null_mut();
        v1 = allocate::<Vertex>();
        v2 = allocate::<Vertex>();
        v3 = allocate::<Vertex>();
        v4 = allocate::<Vertex>();
        vu12 = allocate::<VertexUse>();
        vu14 = allocate::<VertexUse>();
        vu21 = allocate::<VertexUse>();
        vu23 = allocate::<VertexUse>();
        vu24 = allocate::<VertexUse>();
        vu32 = allocate::<VertexUse>();
        vu34 = allocate::<VertexUse>();
        vu41 = allocate::<VertexUse>();
        vu42 = allocate::<VertexUse>();
        vu43 = allocate::<VertexUse>();
        e1 = allocate::<Edge>();
        e2 = allocate::<Edge>();
        e3 = allocate::<Edge>();
        e4 = allocate::<Edge>();
        e5 = allocate::<Edge>();
        eu12 = allocate::<EdgeUse>();
        eu21 = allocate::<EdgeUse>();
        eu23 = allocate::<EdgeUse>();
        eu32 = allocate::<EdgeUse>();
        eu34 = allocate::<EdgeUse>();
        eu43 = allocate::<EdgeUse>();
        eu41 = allocate::<EdgeUse>();
        eu14 = allocate::<EdgeUse>();
        eu24 = allocate::<EdgeUse>();
        eu42 = allocate::<EdgeUse>();
        f124 = allocate::<Face>();
        f234 = allocate::<Face>();
        m = allocate::<Model>();
        (*v1).x = xmin;
        (*v1).y = ymax;
        (*v1).vu = vu12;
        (*v1).fixed = 1;
        (*v1).va = std::ptr::null_mut();
        (*v2).x = xmax;
        (*v2).y = ymax;
        (*v2).vu = vu23;
        (*v2).fixed = 1;
        (*v2).va = std::ptr::null_mut();
        (*v3).x = xmax;
        (*v3).y = ymin;
        (*v3).vu = vu34;
        (*v3).fixed = 1;
        (*v3).va = std::ptr::null_mut();
        (*v4).x = xmin;
        (*v4).y = ymin;
        (*v4).vu = vu41;
        (*v4).fixed = 1;
        (*v4).va = std::ptr::null_mut();
        (*vu12).v = v1;
        (*vu12).eu = eu12;
        (*vu12).vuPrev = vu14;
        (*vu12).vuNext = vu14;
        (*vu12).vua = std::ptr::null_mut();
        (*vu14).v = v1;
        (*vu14).eu = eu14;
        (*vu14).vuPrev = vu12;
        (*vu14).vuNext = vu12;
        (*vu14).vua = std::ptr::null_mut();
        (*vu21).v = v2;
        (*vu21).eu = eu21;
        (*vu21).vuPrev = vu24;
        (*vu21).vuNext = vu23;
        (*vu21).vua = std::ptr::null_mut();
        (*vu23).v = v2;
        (*vu23).eu = eu23;
        (*vu23).vuPrev = vu21;
        (*vu23).vuNext = vu24;
        (*vu23).vua = std::ptr::null_mut();
        (*vu24).v = v2;
        (*vu24).eu = eu24;
        (*vu24).vuPrev = vu23;
        (*vu24).vuNext = vu21;
        (*vu24).vua = std::ptr::null_mut();
        (*vu32).v = v3;
        (*vu32).eu = eu32;
        (*vu32).vuPrev = vu34;
        (*vu32).vuNext = vu34;
        (*vu32).vua = std::ptr::null_mut();
        (*vu34).v = v3;
        (*vu34).eu = eu34;
        (*vu34).vuPrev = vu32;
        (*vu34).vuNext = vu32;
        (*vu34).vua = std::ptr::null_mut();
        (*vu41).v = v4;
        (*vu41).eu = eu41;
        (*vu41).vuPrev = vu43;
        (*vu41).vuNext = vu42;
        (*vu41).vua = std::ptr::null_mut();
        (*vu42).v = v4;
        (*vu42).eu = eu42;
        (*vu42).vuPrev = vu41;
        (*vu42).vuNext = vu43;
        (*vu42).vua = std::ptr::null_mut();
        (*vu43).v = v4;
        (*vu43).eu = eu43;
        (*vu43).vuPrev = vu42;
        (*vu43).vuNext = vu41;
        (*vu43).vua = std::ptr::null_mut();
        (*e1).eu = eu12;
        (*e1).fixed = 1;
        (*e1).ea = std::ptr::null_mut();
        (*e2).eu = eu23;
        (*e2).fixed = 1;
        (*e2).ea = std::ptr::null_mut();
        (*e3).eu = eu34;
        (*e3).fixed = 1;
        (*e3).ea = std::ptr::null_mut();
        (*e4).eu = eu41;
        (*e4).fixed = 1;
        (*e4).ea = std::ptr::null_mut();
        (*e5).eu = eu24;
        (*e5).fixed = 0;
        (*e5).ea = std::ptr::null_mut();
        (*eu12).e = e1;
        (*eu12).euMate = eu21;
        (*eu12).vu = vu12;
        (*eu12).f = f124;
        (*eu21).e = e1;
        (*eu21).euMate = eu12;
        (*eu21).vu = vu21;
        (*eu21).f = std::ptr::null_mut();
        (*eu23).e = e2;
        (*eu23).euMate = eu32;
        (*eu23).vu = vu23;
        (*eu23).f = f234;
        (*eu32).e = e2;
        (*eu32).euMate = eu23;
        (*eu32).vu = vu32;
        (*eu32).f = std::ptr::null_mut();
        (*eu34).e = e3;
        (*eu34).euMate = eu43;
        (*eu34).vu = vu34;
        (*eu34).f = f234;
        (*eu43).e = e3;
        (*eu43).euMate = eu34;
        (*eu43).vu = vu43;
        (*eu43).f = std::ptr::null_mut();
        (*eu41).e = e4;
        (*eu41).euMate = eu14;
        (*eu41).vu = vu41;
        (*eu41).f = f124;
        (*eu14).e = e4;
        (*eu14).euMate = eu41;
        (*eu14).vu = vu14;
        (*eu14).f = std::ptr::null_mut();
        (*eu24).e = e5;
        (*eu24).euMate = eu42;
        (*eu24).vu = vu24;
        (*eu24).f = f124;
        (*eu42).e = e5;
        (*eu42).euMate = eu24;
        (*eu42).vu = vu42;
        (*eu42).f = f234;
        (*eu12).euCW = eu24;
        (*eu12).euCCW = eu41;
        (*eu12).eua = std::ptr::null_mut();
        (*eu21).euCW = eu14;
        (*eu21).euCCW = eu32;
        (*eu21).eua = std::ptr::null_mut();
        (*eu23).euCW = eu34;
        (*eu23).euCCW = eu42;
        (*eu23).eua = std::ptr::null_mut();
        (*eu32).euCW = eu21;
        (*eu32).euCCW = eu43;
        (*eu32).eua = std::ptr::null_mut();
        (*eu34).euCW = eu42;
        (*eu34).euCCW = eu23;
        (*eu34).eua = std::ptr::null_mut();
        (*eu43).euCW = eu32;
        (*eu43).euCCW = eu14;
        (*eu43).eua = std::ptr::null_mut();
        (*eu41).euCW = eu12;
        (*eu41).euCCW = eu24;
        (*eu41).eua = std::ptr::null_mut();
        (*eu14).euCW = eu43;
        (*eu14).euCCW = eu21;
        (*eu14).eua = std::ptr::null_mut();
        (*eu24).euCW = eu41;
        (*eu24).euCCW = eu12;
        (*eu24).eua = std::ptr::null_mut();
        (*eu42).euCW = eu23;
        (*eu42).euCCW = eu34;
        (*eu42).eua = std::ptr::null_mut();
        (*f124).m = {
            (*f234).m = m;
            (*f234).m
        };
        (*f124).fPrev = {
            (*f124).fNext = f234;
            (*f124).fNext
        };
        (*f234).fPrev = {
            (*f234).fNext = f124;
            (*f234).fNext
        };
        (*f124).eu = eu12;
        (*f234).eu = eu23;
        (*f124).fa = std::ptr::null_mut();
        (*f234).fa = std::ptr::null_mut();
        circumTri(&mut *f124);
        circumTri(&mut *f234);
        (*m).f = f124;
        (*m).xmin = xmin;
        (*m).ymin = ymin;
        (*m).xmax = xmax;
        (*m).ymax = ymax;
        (*m).eps = (0.005
            * ((((xmax - xmin) * (xmax - xmin)) + ((ymax - ymin) * (ymax - ymin))) as f64).sqrt())
            as f32;
        (*m).vAdd = None;
        (*m).vDel = None;
        (*m).tAdd = None;
        (*m).tDel = None;
        (*m).sma = {
            (*m).sfa = {
                (*m).seua = {
                    (*m).sea = {
                        (*m).svua = {
                            (*m).sva = (0) as usize;
                            (*m).sva
                        };
                        (*m).svua
                    };
                    (*m).sea
                };
                (*m).seua
            };
            (*m).sfa
        };
        return m;
    }
}

/// C: `void makeBoundaryEdgeTri (Vertex *v, Edge **enew, Tri **tnew);`
/// # Safety
/// v must be on the model boundary, collinear with its boundary neighbors.
/// All topology links and the two output pointers must be valid and uniquely
/// mutable. The resulting zero-area triangle permits boundary vertex deletion.
pub unsafe fn makeBoundaryEdgeTri(
    mut v: *mut Vertex,
    mut enew: *mut *mut Edge,
    mut tnew: *mut *mut Tri,
) {
    unsafe {
        let mut vi: *mut Vertex = std::ptr::null_mut();
        let mut vo: *mut Vertex = std::ptr::null_mut();
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut vui: *mut VertexUse = std::ptr::null_mut();
        let mut vuo: *mut VertexUse = std::ptr::null_mut();
        let mut e: *mut Edge = std::ptr::null_mut();
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        let mut eui: *mut EdgeUse = std::ptr::null_mut();
        let mut euo: *mut EdgeUse = std::ptr::null_mut();
        let mut f: *mut Face = std::ptr::null_mut();
        vu = (*v).vu;
        loop {
            {
                eu = (*vu).eu;
                vu = (*vu).vuNext;
            }
            if !((*eu).f != std::ptr::null_mut()) {
                break;
            }
        }
        e = allocate::<Edge>();
        (*e).ea = std::ptr::null_mut();
        f = allocate::<Face>();
        (*f).fa = std::ptr::null_mut();
        eui = allocate::<EdgeUse>();
        (*eui).eua = std::ptr::null_mut();
        euo = allocate::<EdgeUse>();
        (*euo).eua = std::ptr::null_mut();
        vui = allocate::<VertexUse>();
        (*vui).vua = std::ptr::null_mut();
        vuo = allocate::<VertexUse>();
        (*vuo).vua = std::ptr::null_mut();
        vi = (*(*(*eu).euCW).vu).v;
        vo = (*(*(*eu).euCCW).vu).v;
        (*vui).eu = eui;
        (*vui).vuPrev = (*vi).vu;
        (*vui).vuNext = (*(*vi).vu).vuNext;
        (*(*vui).vuPrev).vuNext = vui;
        (*(*vui).vuNext).vuPrev = vui;
        (*vui).v = vi;
        (*vuo).eu = euo;
        (*vuo).vuPrev = (*vo).vu;
        (*vuo).vuNext = (*(*vo).vu).vuNext;
        (*(*vuo).vuPrev).vuNext = vuo;
        (*(*vuo).vuNext).vuPrev = vuo;
        (*vuo).v = vo;
        (*e).eu = eui;
        (*e).fixed = 1;
        (*eu).f = {
            (*(*eu).euCCW).f = f;
            (*(*eu).euCCW).f
        };
        (*eui).f = f;
        (*eui).vu = vui;
        (*eui).euMate = euo;
        (*eui).euCW = (*eu).euCCW;
        (*eui).euCCW = eu;
        (*eui).e = e;
        (*euo).f = std::ptr::null_mut();
        (*euo).vu = vuo;
        (*euo).euMate = eui;
        (*euo).euCW = (*eu).euCW;
        (*euo).euCCW = (*(*eu).euCCW).euCCW;
        (*euo).e = e;
        (*(*eui).euCW).euCCW = eui;
        (*(*eui).euCCW).euCW = eui;
        (*(*euo).euCW).euCCW = euo;
        (*(*euo).euCCW).euCW = euo;
        (*f).m = (*(*(*eu).euMate).f).m;
        (*f).fPrev = (*(*f).m).f;
        (*f).fNext = (*(*(*f).m).f).fNext;
        (*(*f).fPrev).fNext = f;
        (*(*f).fNext).fPrev = f;
        (*f).eu = eu;
        circumTri(&mut *f);
        (*enew) = e;
        (*tnew) = (f) as *mut Tri;
    }
}

/// C: `void makeEdgeFace (Vertex *v1, Vertex *v2, Edge **enew, Face **fnew);`
/// # Safety
/// The two vertices must belong to one common interior face. All topology
/// links and output pointers must be valid and uniquely mutable. This routine
/// splits a polygon; the caller computes circumcircles after triangulation.
pub unsafe fn makeEdgeFace(
    mut v1: *mut Vertex,
    mut v2: *mut Vertex,
    mut enew: *mut *mut Edge,
    mut fnew: *mut *mut Face,
) {
    unsafe {
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut vum: *mut VertexUse = std::ptr::null_mut();
        let mut vu1: *mut VertexUse = std::ptr::null_mut();
        let mut vu2: *mut VertexUse = std::ptr::null_mut();
        let mut e: *mut Edge = std::ptr::null_mut();
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        let mut eum: *mut EdgeUse = std::ptr::null_mut();
        let mut eu1: *mut EdgeUse = std::ptr::null_mut();
        let mut eu2: *mut EdgeUse = std::ptr::null_mut();
        let mut eut: *mut EdgeUse = std::ptr::null_mut();
        let mut f: *mut Face = std::ptr::null_mut();
        let mut fm: *mut Face = std::ptr::null_mut();
        vu2 = std::ptr::null_mut();
        vu1 = (*v1).vu;
        while (vu2 == std::ptr::null_mut()) {
            {
                vu = (*(*(*vu1).eu).euCW).vu;
                vu
            };
            while (((vu != vu1) && ((*(*vu).eu).f != std::ptr::null_mut())) && ((*vu).v != v2)) {
                {}
                {
                    vu = (*(*(*vu).eu).euCW).vu;
                    vu
                };
            }
            if ((*vu).v == v2) {
                vu2 = vu;
            } else {
                vu1 = (*vu1).vuNext;
            }
        }
        vu = allocate::<VertexUse>();
        (*vu).vua = std::ptr::null_mut();
        (*vu).v = v2;
        (*vu).vuPrev = vu2;
        (*vu).vuNext = (*vu2).vuNext;
        (*(*vu2).vuNext).vuPrev = vu;
        (*vu2).vuNext = vu;
        vum = allocate::<VertexUse>();
        (*vum).vua = std::ptr::null_mut();
        (*vum).v = v1;
        (*vum).vuPrev = vu1;
        (*vum).vuNext = (*vu1).vuNext;
        (*(*vu1).vuNext).vuPrev = vum;
        (*vu1).vuNext = vum;
        eu1 = (*vu1).eu;
        eu2 = (*vu2).eu;
        e = allocate::<Edge>();
        (*e).ea = std::ptr::null_mut();
        (*e).eu = {
            eu = allocate::<EdgeUse>();
            eu
        };
        (*eu).eua = std::ptr::null_mut();
        (*e).fixed = 0;
        eum = allocate::<EdgeUse>();
        (*eum).eua = std::ptr::null_mut();
        (*eu).e = {
            (*eum).e = e;
            (*eum).e
        };
        (*eu).vu = vu;
        (*vu).eu = eu;
        (*eum).vu = vum;
        (*vum).eu = eum;
        (*eu).euMate = eum;
        (*eum).euMate = eu;
        (*eu).euCW = eu1;
        (*eu).euCCW = (*eu2).euCCW;
        (*eum).euCW = eu2;
        (*eum).euCCW = (*eu1).euCCW;
        (*(*eu1).euCCW).euCW = eum;
        (*(*eu2).euCCW).euCW = eu;
        (*eu1).euCCW = eu;
        (*eu2).euCCW = eum;
        f = (*eu1).f;
        fm = allocate::<Face>();
        (*fm).fa = std::ptr::null_mut();
        (*fm).eu = eum;
        (*fm).fPrev = f;
        (*fm).fNext = (*f).fNext;
        (*(*f).fNext).fPrev = fm;
        (*f).fNext = fm;
        (*fm).m = (*f).m;
        (*f).eu = eu1;
        (*eu).f = f;
        (*eum).f = fm;
        {
            eut = eu2;
            eut
        };
        while ((*eut).f != fm) {
            {
                (*eut).f = fm;
            }
            {
                eut = (*eut).euCW;
                eut
            };
        }
        (*enew) = e;
        (*fnew) = fm;
    }
}

/// C: `void makeEdgeVertex (Vertex *v1, float x, float y, Face *f, Edge **enew, Vertex **vnew);`
/// # Safety
/// v1 must be adjacent to f. All topology links and output pointers must be
/// valid and uniquely mutable. The new edge initially has the same face on
/// both sides; makeEdgeFace completes the subdivision.
pub unsafe fn makeEdgeVertex(
    mut v1: *mut Vertex,
    mut x: f32,
    mut y: f32,
    mut f: *mut Face,
    mut enew: *mut *mut Edge,
    mut vnew: *mut *mut Vertex,
) {
    unsafe {
        let mut v: *mut Vertex = std::ptr::null_mut();
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut vum: *mut VertexUse = std::ptr::null_mut();
        let mut vu1: *mut VertexUse = std::ptr::null_mut();
        let mut e: *mut Edge = std::ptr::null_mut();
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        let mut eum: *mut EdgeUse = std::ptr::null_mut();
        let mut eu1: *mut EdgeUse = std::ptr::null_mut();
        {
            vu1 = (*v1).vu;
            vu1
        };
        while ((*(*vu1).eu).f != f) {
            {}
            {
                vu1 = (*vu1).vuNext;
                vu1
            };
        }
        v = allocate::<Vertex>();
        (*v).va = std::ptr::null_mut();
        vu = allocate::<VertexUse>();
        (*vu).vua = std::ptr::null_mut();
        (*v).x = x;
        (*v).y = y;
        (*v).fixed = 0;
        (*v).vu = vu;
        (*vu).v = v;
        (*vu).vuNext = {
            (*vu).vuPrev = vu;
            (*vu).vuPrev
        };
        vum = allocate::<VertexUse>();
        (*vum).vua = std::ptr::null_mut();
        (*vum).v = v1;
        (*vum).vuPrev = vu1;
        (*vum).vuNext = (*vu1).vuNext;
        (*(*vu1).vuNext).vuPrev = vum;
        (*vu1).vuNext = vum;
        eu1 = (*vu1).eu;
        e = allocate::<Edge>();
        (*e).ea = std::ptr::null_mut();
        (*e).eu = {
            eu = allocate::<EdgeUse>();
            eu
        };
        (*eu).eua = std::ptr::null_mut();
        (*e).fixed = 0;
        eum = allocate::<EdgeUse>();
        (*eum).eua = std::ptr::null_mut();
        (*eu).e = {
            (*eum).e = e;
            (*eum).e
        };
        (*eu).f = {
            (*eum).f = f;
            (*eum).f
        };
        (*eu).vu = vu;
        (*vu).eu = eu;
        (*eum).vu = vum;
        (*vum).eu = eum;
        (*eu).euMate = eum;
        (*eum).euMate = eu;
        (*eu).euCW = eu1;
        (*eu).euCCW = eum;
        (*eum).euCW = eu;
        (*eum).euCCW = (*eu1).euCCW;
        (*(*eu1).euCCW).euCW = eum;
        (*eu1).euCCW = eu;
        (*enew) = e;
        (*vnew) = v;
    }
}

/// C: `Vertex *addVertexToModel (Model *m, float x, float y);`
/// Returns null if the point is within m.eps of an existing vertex. Projects
/// nearby points onto fixed edges and preserves both halves of a fixed edge.
/// # Safety
/// m must be a valid triangulated model, uniquely mutable for the call.
/// Callbacks must preserve the topology links needed by this operation.
pub unsafe fn addVertexToModel(mut m: *mut Model, mut x: f32, mut y: f32) -> *mut Vertex {
    unsafe {
        let mut first: i32 = 0;
        let mut fixed: i32 = 0;
        let mut v: *mut Vertex = std::ptr::null_mut();
        let mut vo: *mut Vertex = std::ptr::null_mut();
        let mut v1: *mut Vertex = std::ptr::null_mut();
        let mut v2: *mut Vertex = std::ptr::null_mut();
        let mut v3: *mut Vertex = std::ptr::null_mut();
        let mut v1fixed: *mut Vertex = std::ptr::null_mut();
        let mut v2fixed: *mut Vertex = std::ptr::null_mut();
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut e: *mut Edge = std::ptr::null_mut();
        let mut e1: *mut Edge = std::ptr::null_mut();
        let mut e2: *mut Edge = std::ptr::null_mut();
        let mut e3: *mut Edge = std::ptr::null_mut();
        let mut eb: *mut Edge = std::ptr::null_mut();
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        let mut euo: *mut EdgeUse = std::ptr::null_mut();
        let mut f: *mut Face = std::ptr::null_mut();
        let mut fs: *mut Face = std::ptr::null_mut();
        let mut f1: *mut Face = std::ptr::null_mut();
        let mut f2: *mut Face = std::ptr::null_mut();
        let mut f3: *mut Face = std::ptr::null_mut();
        // Reject vertices closer than the model tolerance.
        v = nearestVertexInModel(m, std::ptr::null_mut(), x, y);
        if ((((((*v).x - x) * ((*v).x - x)) + (((*v).y - y) * ((*v).y - y))) as f64).sqrt()
            < ((*m).eps) as f64)
        {
            return std::ptr::null_mut();
        }
        f = (if ((*(*(*v).vu).eu).f != std::ptr::null_mut()) {
            (*(*(*v).vu).eu).f
        } else {
            (*(*(*(*v).vu).eu).euMate).f
        });
        f = insideTriInModel(m, f, x, y);
        e = nearestEdgeInModel(m, (*(*f).eu).e, x, y);
        if ((*e).fixed != 0 && (distanceToEdge(&*e, x, y) < (*m).eps)) {
            v1fixed = (*(*(*e).eu).vu).v;
            v2fixed = (*(*(*(*e).eu).euCW).vu).v;
            (*e).fixed = 0;
            projectToEdge(&*e, &mut *(&raw mut x), &mut *(&raw mut y));
            fixed = 1;
            if (((*(*e).eu).f == std::ptr::null_mut())
                || ((*(*(*e).eu).euMate).f == std::ptr::null_mut()))
            {
                eb = e;
            }
        }
        if (*m).tDel.is_some() {
            (*m).tDel.unwrap()(m, f);
        }
        // Split the containing triangle into three triangles.
        vu = (*(*f).eu).vu;
        v1 = (*vu).v;
        v2 = (*(*(*(*vu).eu).euCW).vu).v;
        v3 = (*(*(*(*vu).eu).euCCW).vu).v;
        f1 = f;
        makeEdgeVertex(v1, x, y, f1, &raw mut e1, &raw mut v);
        makeEdgeFace(v2, v, &raw mut e2, &raw mut f2);
        makeEdgeFace(v3, v, &raw mut e3, &raw mut f3);
        if (*m).tAdd.is_some() {
            (*m).tAdd.unwrap()(m, f1);
            (*m).tAdd.unwrap()(m, f2);
            (*m).tAdd.unwrap()(m, f3);
        }
        vu = (*v).vu;
        loop {
            {
                circumTri(&mut *(*(*vu).eu).f);
                vu = (*vu).vuNext;
            }
            if !(vu != (*v).vu) {
                break;
            }
        }
        vu = (*v).vu;
        // Visit opposite edges clockwise and swap those failing the Delaunay
        // circumcircle test. Stay on the same vertex-use after a swap.
        first = 1;
        loop {
            {
                eu = (*(*vu).eu).euCW;
                euo = (*eu).euMate;
                e = (*eu).e;
                if ((((*euo).f != std::ptr::null_mut()) && !((*(*euo).e).fixed != 0))
                    && inCircumTri(x, y, &*(*euo).f) != 0)
                {
                    vo = (*(*(*euo).euCCW).vu).v;
                    if (*m).tDel.is_some() {
                        (*m).tDel.unwrap()(m, (*(*e).eu).f);
                        (*m).tDel.unwrap()(m, (*(*(*e).eu).euMate).f);
                    }
                    killEdge(e, &raw mut fs);
                    makeEdgeFace(v, vo, &raw mut e, &raw mut f);
                    if (*m).tAdd.is_some() {
                        (*m).tAdd.unwrap()(m, (*(*e).eu).f);
                        (*m).tAdd.unwrap()(m, (*(*(*e).eu).euMate).f);
                    }
                    circumTri(&mut *(*(*e).eu).f);
                    circumTri(&mut *(*(*(*e).eu).euMate).f);
                } else {
                    vu = (*(*(*(*vu).eu).euCCW).euMate).vu;
                    first = 0;
                }
            }
            if !((vu != (*v).vu) || first != 0) {
                break;
            }
        }
        if fixed != 0 {
            vu = (*v).vu;
            loop {
                {
                    eu = (*vu).eu;
                    if ((*(*(*eu).euCCW).vu).v == v1fixed) {
                        (*(*(*eu).euCCW).e).fixed = 1;
                    }
                    if ((*(*(*eu).euCW).vu).v == v2fixed) {
                        (*(*eu).e).fixed = 1;
                    }
                    vu = (*vu).vuNext;
                }
                if !(vu != (*v).vu) {
                    break;
                }
            }
        }
        if (eb != std::ptr::null_mut()) {
            killBoundaryEdge(eb);
        }
        if (*m).vAdd.is_some() {
            (*m).vAdd.unwrap()(m, v);
        }
        return v;
    }
}
