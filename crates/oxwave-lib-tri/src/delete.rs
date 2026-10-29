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

//! Port of tri/lib/delete.c.
//! Author: Dave Hale, Colorado School of Mines, Fall 1990.

#![allow(unused_assignments, unused_parens, unused_mut)]
use super::*;

/// C: `void deleteVertexFromModel (Model *m, Vertex *v);`
/// Fixed vertices are retained. A boundary vertex must be collinear with its
/// boundary neighbors, as required by makeBoundaryEdgeTri.
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn deleteVertexFromModel(mut m: *mut Model, mut v: *mut Vertex) {
    unsafe {
        let mut containsb: i32 = 0;
        let mut edgedeleted: i32 = 0;
        let mut x: f32 = (0) as f32;
        let mut y: f32 = (0) as f32;
        let mut xb: f32 = (0) as f32;
        let mut yb: f32 = (0) as f32;
        let mut xc: f32 = (0) as f32;
        let mut yc: f32 = (0) as f32;
        let mut rs: f32 = (0) as f32;
        let mut v1: *mut Vertex = std::ptr::null_mut();
        let mut v2: *mut Vertex = std::ptr::null_mut();
        let mut v3: *mut Vertex = std::ptr::null_mut();
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut vub: *mut VertexUse = std::ptr::null_mut();
        let mut e: *mut Edge = std::ptr::null_mut();
        let mut f: *mut Face = std::ptr::null_mut();
        let mut fs: *mut Face = std::ptr::null_mut();
        if (*v).fixed != 0 {
            return;
        }
        vu = (*v).vu;
        loop {
            {
                if ((*(*vu).eu).f == std::ptr::null_mut()) {
                    makeBoundaryEdgeTri(v, &raw mut e, &raw mut f);
                    break;
                }
                vu = (*vu).vuNext;
            }
            if !(vu != (*v).vu) {
                break;
            }
        }
        if (*m).vDel.is_some() {
            (*m).vDel.unwrap()(m, v);
        }
        x = (*v).x;
        y = (*v).y;
        vu = (*v).vu;
        while ((*(*(*vu).vuNext).vuNext).vuNext != vu) {
            edgedeleted = 0;
            v1 = (*(*(*(*(*vu).eu).euMate).euCCW).vu).v;
            v2 = (*(*(*(*vu).eu).euMate).vu).v;
            v3 = (*(*(*(*vu).eu).euCCW).vu).v;
            circum(
                (*v1).x,
                (*v1).y,
                (*v2).x,
                (*v2).y,
                (*v3).x,
                (*v3).y,
                &mut *(&raw mut xc),
                &mut *(&raw mut yc),
                &mut *(&raw mut rs),
            );
            if (inCircum(x, y, xc, yc, rs) != 0 && !(in3Vertices(x, y, &*v1, &*v2, &*v3) != 0)) {
                vub = (*(*(*(*(*vu).eu).euCCW).euMate).euCCW).vu;
                loop {
                    {
                        xb = (*(*vub).v).x;
                        yb = (*(*vub).v).y;
                        containsb = inCircum(xb, yb, xc, yc, rs);
                        vub = (*(*(*(*vub).eu).euMate).euCCW).vu;
                    }
                    if !(((*vub).v != v1) && !(containsb != 0)) {
                        break;
                    }
                }
                if !(containsb != 0) {
                    if (*m).tDel.is_some() {
                        (*m).tDel.unwrap()(m, (*(*vu).eu).f);
                        (*m).tDel.unwrap()(m, (*(*(*vu).eu).euMate).f);
                    }
                    killEdge((*(*vu).eu).e, &raw mut fs);
                    makeEdgeFace(v1, v3, &raw mut e, &raw mut f);
                    if (*m).tAdd.is_some() {
                        (*m).tAdd.unwrap()(m, (*(*e).eu).f);
                        (*m).tAdd.unwrap()(m, (*(*(*e).eu).euMate).f);
                    }
                    circumTri(&mut *(*(*e).eu).f);
                    circumTri(&mut *(*(*(*e).eu).euMate).f);
                    edgedeleted = 1;
                }
            }
            if edgedeleted != 0 {
                vu = (*v).vu;
            } else {
                vu = (*vu).vuNext;
            }
        }
        if (*m).tDel.is_some() {
            vu = (*v).vu;
            loop {
                {
                    (*m).tDel.unwrap()(m, (*(*vu).eu).f);
                    vu = (*vu).vuNext;
                }
                if !(vu != (*v).vu) {
                    break;
                }
            }
        }
        killEdge((*(*(*v).vu).eu).e, &raw mut fs);
        killEdge((*(*(*v).vu).eu).e, &raw mut fs);
        killEdge((*(*(*v).vu).eu).e, &raw mut fs);
        if (*m).tAdd.is_some() {
            (*m).tAdd.unwrap()(m, fs);
        }
        circumTri(&mut *fs);
    }
}

/// C: `void killModel (Model *m);`
/// Preserves C's ownership behavior: attribute buffers and exterior topology
/// allocations are not freed by this routine.
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn killModel(mut m: *mut Model) {
    unsafe {
        let mut f: *mut Face = std::ptr::null_mut();
        let mut fnext: *mut Face = std::ptr::null_mut();
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        let mut eunext: *mut EdgeUse = std::ptr::null_mut();
        f = (*m).f;
        loop {
            {
                fnext = (*f).fNext;
                eu = (*f).eu;
                loop {
                    {
                        eunext = (*eu).euCW;
                        if ((*(*eu).vu).vuNext == (*eu).vu) {
                            free(((*(*eu).vu).v) as *mut c_void);
                        } else {
                            (*(*(*eu).vu).vuPrev).vuNext = (*(*eu).vu).vuNext;
                            (*(*(*eu).vu).vuNext).vuPrev = (*(*eu).vu).vuPrev;
                        }
                        free(((*eu).vu) as *mut c_void);
                        if ((*eu).euMate == std::ptr::null_mut()) {
                            free(((*eu).e) as *mut c_void);
                        } else {
                            (*(*eu).euMate).euMate = std::ptr::null_mut();
                        }
                        free((eu) as *mut c_void);
                        eu = eunext;
                    }
                    if !(eu != (*f).eu) {
                        break;
                    }
                }
                free((f) as *mut c_void);
                f = fnext;
            }
            if !(f != (*m).f) {
                break;
            }
        }
        free((m) as *mut c_void);
    }
}

/// C: `void killEdge (Edge *e, Face **fs);`
/// Merges the two adjacent interior faces. If both uses already belong to the
/// same face, fs is left untouched. As in C, the Edge allocation itself remains
/// allocated; callers owning it may free it after the topology is detached.
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn killEdge(mut e: *mut Edge, mut fs: *mut *mut Face) {
    unsafe {
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut vum: *mut VertexUse = std::ptr::null_mut();
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        let mut eum: *mut EdgeUse = std::ptr::null_mut();
        let mut eut: *mut EdgeUse = std::ptr::null_mut();
        let mut f: *mut Face = std::ptr::null_mut();
        let mut fm: *mut Face = std::ptr::null_mut();
        eu = (*e).eu;
        eum = (*eu).euMate;
        vu = (*eu).vu;
        vum = (*eum).vu;
        f = (*eu).f;
        fm = (*eum).f;
        if (f != fm) {
            (*(*fm).m).f = f;
            (*(*fm).fPrev).fNext = (*fm).fNext;
            (*(*fm).fNext).fPrev = (*fm).fPrev;
            free((fm) as *mut c_void);
            {
                eut = (*eum).euCW;
                eut
            };
            while ((*eut).f != f) {
                {
                    (*eut).f = f;
                }
                {
                    eut = (*eut).euCW;
                    eut
                };
            }
            (*fs) = f;
        }
        (*f).eu = (*(*eu).euCW).euCW;
        (*(*eu).euCCW).euCW = (*eum).euCW;
        (*(*eum).euCCW).euCW = (*eu).euCW;
        (*(*eu).euCW).euCCW = (*eum).euCCW;
        (*(*eum).euCW).euCCW = (*eu).euCCW;
        free((eu) as *mut c_void);
        free((eum) as *mut c_void);
        if ((*vu).vuNext == vu) {
            free(((*vu).v) as *mut c_void);
        } else {
            (*(*vu).v).vu = (*vu).vuNext;
            (*(*vu).vuPrev).vuNext = (*vu).vuNext;
            (*(*vu).vuNext).vuPrev = (*vu).vuPrev;
            free((vu) as *mut c_void);
        }
        if ((*vum).vuNext == vum) {
            free(((*vum).v) as *mut c_void);
        } else {
            (*(*vum).v).vu = (*vum).vuNext;
            (*(*vum).vuPrev).vuNext = (*vum).vuNext;
            (*(*vum).vuNext).vuPrev = (*vum).vuPrev;
            free((vum) as *mut c_void);
        }
    }
}

/// C: `void killBoundaryEdge (Edge *e);`
/// # Safety
/// e must separate the exterior from a face to remove, with other interior
/// faces remaining in the model. All topology links must be valid and uniquely
/// mutable. This normally removes a degenerate triangle after edge splitting.
pub unsafe fn killBoundaryEdge(mut e: *mut Edge) {
    unsafe {
        let mut f: *mut Face = std::ptr::null_mut();
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        let mut eui: *mut EdgeUse = std::ptr::null_mut();
        let mut euo: *mut EdgeUse = std::ptr::null_mut();
        let mut vui: *mut VertexUse = std::ptr::null_mut();
        let mut vuo: *mut VertexUse = std::ptr::null_mut();
        if ((*(*e).eu).f == std::ptr::null_mut()) {
            euo = (*e).eu;
            eui = (*euo).euMate;
        } else {
            eui = (*e).eu;
            euo = (*eui).euMate;
        }
        vui = (*eui).vu;
        vuo = (*euo).vu;
        f = (*eui).f;
        (*(*f).fNext).fPrev = (*f).fPrev;
        (*(*f).fPrev).fNext = (*f).fNext;
        (*(*f).m).f = (*f).fNext;
        free((f) as *mut c_void);
        eu = (*eui).euCW;
        loop {
            {
                (*eu).f = std::ptr::null_mut();
                (*(*eu).e).fixed = 1;
                eu = (*eu).euCW;
            }
            if !(eu != eui) {
                break;
            }
        }
        (*(*vui).vuNext).vuPrev = (*vui).vuPrev;
        (*(*vui).vuPrev).vuNext = (*vui).vuNext;
        (*(*vui).v).vu = (*vui).vuNext;
        free((vui) as *mut c_void);
        (*(*vuo).vuNext).vuPrev = (*vuo).vuPrev;
        (*(*vuo).vuPrev).vuNext = (*vuo).vuNext;
        (*(*vuo).v).vu = (*vuo).vuNext;
        free((vuo) as *mut c_void);
        (*(*eui).euCW).euCCW = (*euo).euCCW;
        (*(*eui).euCCW).euCW = (*euo).euCW;
        (*(*euo).euCW).euCCW = (*eui).euCCW;
        (*(*euo).euCCW).euCW = (*eui).euCW;
        free((eui) as *mut c_void);
        free((euo) as *mut c_void);
        free((e) as *mut c_void);
    }
}
