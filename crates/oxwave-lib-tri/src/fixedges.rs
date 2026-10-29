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

//! Port of tri/lib/fixedges.c.
//! Author: Dave Hale, Colorado School of Mines, 06/04/91.
//! unfixEdge and correctTrisAroundVertex: Chris Elmer & Dave Hale,
//! Colorado School of Mines, 06/04/91.

#![allow(unused_assignments, unused_parens, unused_mut)]
use super::*;

/// C: `int fixEdgeBetweenVertices (Vertex *v1, Vertex *v2);`
/// Fixes existing collinear edge chains, replacing intersecting edges and
/// retriangulating on both sides when necessary.
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn fixEdgeBetweenVertices(mut v1: *mut Vertex, mut v2: *mut Vertex) -> i32 {
    unsafe {
        let mut e: *mut Edge = std::ptr::null_mut();
        let mut f: *mut Face = std::ptr::null_mut();
        fixExistingEdgesBetween(&raw mut v1, &raw mut v2);
        if (v1 == v2) {
            return 1;
        }
        killEdgesBetween(v1, v2);
        makeEdgeFace(v1, v2, &raw mut e, &raw mut f);
        (*e).fixed = 1;
        triangFaceNextToEdge((*e).eu);
        triangFaceNextToEdge((*(*e).eu).euMate);
        return 1;
    }
}
unsafe fn fixExistingEdgesBetween(mut v1: *mut *mut Vertex, mut v2: *mut *mut Vertex) {
    unsafe {
        let mut moved: i32 = 0;
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut v: *mut Vertex = std::ptr::null_mut();
        let mut v1t: *mut Vertex = (*v1);
        let mut v2t: *mut Vertex = (*v2);
        moved = 1;
        while moved != 0 {
            moved = 0;
            vu = (*v1t).vu;
            loop {
                {
                    v = (*(*(*(*vu).eu).euCW).vu).v;
                    if vertexBetweenVertices(&*v, &*v1t, &*v2t) != 0 {
                        (*(*(*vu).eu).e).fixed = 1;
                        v1t = v;
                        moved = 1;
                    }
                    vu = (*vu).vuNext;
                }
                if !(((vu != (*v1t).vu) && !(moved != 0)) && (v1t != v2t)) {
                    break;
                }
            }
            if (v1t == v2t) {
                (*v1) = v1t;
                (*v2) = v2t;
                return;
            }
        }
        moved = 1;
        while moved != 0 {
            moved = 0;
            vu = (*v2t).vu;
            loop {
                {
                    v = (*(*(*(*vu).eu).euCW).vu).v;
                    if vertexBetweenVertices(&*v, &*v1t, &*v2t) != 0 {
                        (*(*(*vu).eu).e).fixed = 1;
                        v2t = v;
                        moved = 1;
                    }
                    vu = (*vu).vuNext;
                }
                if !(((vu != (*v2t).vu) && !(moved != 0)) && (v1t != v2t)) {
                    break;
                }
            }
            if (v1t == v2t) {
                (*v1) = v1t;
                (*v2) = v2t;
                return;
            }
        }
        (*v1) = v1t;
        (*v2) = v2t;
    }
}
unsafe fn killEdgesBetween(mut v1: *mut Vertex, mut v2: *mut Vertex) {
    unsafe {
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        let mut fs: *mut Face = std::ptr::null_mut();
        let mut m: *mut Model = std::ptr::null_mut();
        m = (if ((*(*(*v1).vu).eu).f != std::ptr::null_mut()) {
            (*(*(*(*v1).vu).eu).f).m
        } else {
            (*(*(*(*(*v1).vu).eu).euMate).f).m
        });
        {
            vu = (*v1).vu;
            vu
        };
        while !(edgeBetweenVertices((*(*(*vu).eu).euCW).e, v1, v2) != 0) {
            {}
            {
                vu = (*vu).vuNext;
                vu
            };
        }
        if (*m).tDel.is_some() {
            (*m).tDel.unwrap()(m, (*(*vu).eu).f);
        }
        eu = (*vu).eu;
        loop {
            {
                if (*m).tDel.is_some() {
                    (*m).tDel.unwrap()(m, (*(*(*eu).euCW).euMate).f);
                }
                killEdge((*(*eu).euCW).e, &raw mut fs);
                if !(edgeBetweenVertices((*(*eu).euCW).e, v1, v2) != 0) {
                    eu = (*eu).euCW;
                }
            }
            if !((*(*(*eu).euCW).vu).v != v2) {
                break;
            }
        }
    }
}
unsafe fn edgeBetweenVertices(mut e: *mut Edge, mut v1: *mut Vertex, mut v2: *mut Vertex) -> i32 {
    unsafe {
        let mut e1: *mut Vertex = std::ptr::null_mut();
        let mut e2: *mut Vertex = std::ptr::null_mut();
        let mut dxv2v1: f32 = (0) as f32;
        let mut dyv2v1: f32 = (0) as f32;
        let mut dxe2e1: f32 = (0) as f32;
        let mut dye2e1: f32 = (0) as f32;
        let mut dxe1v1: f32 = (0) as f32;
        let mut dye1v1: f32 = (0) as f32;
        let mut det: f32 = (0) as f32;
        let mut s: f32 = (0) as f32;
        let mut t: f32 = (0) as f32;
        e1 = (*(*(*e).eu).vu).v;
        e2 = (*(*(*(*e).eu).euCW).vu).v;
        if ((((e1 == v1) || (e1 == v2)) || (e2 == v1)) || (e2 == v2)) {
            return 0;
        }
        dxv2v1 = ((*v2).x - (*v1).x);
        dyv2v1 = ((*v2).y - (*v1).y);
        dxe2e1 = ((*e2).x - (*e1).x);
        dye2e1 = ((*e2).y - (*e1).y);
        dxe1v1 = ((*e1).x - (*v1).x);
        dye1v1 = ((*e1).y - (*v1).y);
        det = ((dxe2e1 * dyv2v1) - (dye2e1 * dxv2v1));
        if ((det) as f64 == 0.0) {
            return 0;
        }
        s = (((dxe2e1 * dye1v1) - (dye2e1 * dxe1v1)) / det);
        if (((s) as f64 <= 0.0) || ((s) as f64 >= 1.0)) {
            return 0;
        }
        t = (((dxv2v1 * dye1v1) - (dyv2v1 * dxe1v1)) / det);
        if (((t) as f64 <= 0.0) || ((t) as f64 >= 1.0)) {
            return 0;
        }
        return 1;
    }
}
unsafe fn triangFaceNextToEdge(mut eu: *mut EdgeUse) {
    unsafe {
        let mut xa: f32 = (0) as f32;
        let mut ya: f32 = (0) as f32;
        let mut xb: f32 = (0) as f32;
        let mut yb: f32 = (0) as f32;
        let mut x: f32 = (0) as f32;
        let mut y: f32 = (0) as f32;
        let mut dxa: f32 = (0) as f32;
        let mut dya: f32 = (0) as f32;
        let mut dxb: f32 = (0) as f32;
        let mut dyb: f32 = (0) as f32;
        let mut cota: f32 = (0) as f32;
        let mut cotamin: f32 = (0) as f32;
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut e: *mut Edge = std::ptr::null_mut();
        let mut eua: *mut EdgeUse = std::ptr::null_mut();
        let mut eub: *mut EdgeUse = std::ptr::null_mut();
        let mut f: *mut Face = std::ptr::null_mut();
        let mut m: *mut Model = std::ptr::null_mut();
        m = (*(*eu).f).m;
        if ((*(*(*eu).euCW).euCW).euCW == eu) {
            circumTri(&mut *(*eu).f);
            if (*m).tAdd.is_some() {
                (*m).tAdd.unwrap()(m, (*eu).f);
            }
            return;
        }
        eua = eu;
        eub = (*eu).euCW;
        xa = (*(*(*eua).vu).v).x;
        ya = (*(*(*eua).vu).v).y;
        xb = (*(*(*eub).vu).v).x;
        yb = (*(*(*eub).vu).v).y;
        cotamin = f32::MAX;
        eu = (*eub).euCW;
        loop {
            {
                x = (*(*(*eu).vu).v).x;
                y = (*(*(*eu).vu).v).y;
                dxa = (xa - x);
                dya = (ya - y);
                dxb = (xb - x);
                dyb = (yb - y);
                cota = (((dxa * dxb) + (dya * dyb)) / ((dxb * dya) - (dxa * dyb)));
                if (cota < cotamin) {
                    vu = (*eu).vu;
                    cotamin = cota;
                }
                eu = (*eu).euCW;
            }
            if !(eu != eua) {
                break;
            }
        }
        if ((*(*vu).eu).euCCW != eub) {
            makeEdgeFace((*vu).v, (*(*eub).vu).v, &raw mut e, &raw mut f);
            eu = (if ((*(*e).eu).f == (*eua).f) {
                (*(*e).eu).euMate
            } else {
                (*e).eu
            });
            triangFaceNextToEdge(eu);
        }
        if ((*(*vu).eu).euCW != eua) {
            makeEdgeFace((*vu).v, (*(*eua).vu).v, &raw mut e, &raw mut f);
            eu = (if ((*(*e).eu).f == (*eua).f) {
                (*(*e).eu).euMate
            } else {
                (*e).eu
            });
            triangFaceNextToEdge(eu);
        }
        circumTri(&mut *(*eua).f);
        if (*m).tAdd.is_some() {
            (*m).tAdd.unwrap()(m, (*eua).f);
        }
    }
}

/// C: `int unfixEdge(Edge *e);`
/// Returns zero for boundary edges, otherwise restores the Delaunay criterion.
/// Saves both endpoint pointers before correction to avoid the original C
/// use-after-free when the first correction swaps this edge.
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn unfixEdge(mut e: *mut Edge) -> i32 {
    unsafe {
        if (((*(*e).eu).f == std::ptr::null_mut())
            || ((*(*(*e).eu).euMate).f == std::ptr::null_mut()))
        {
            return 0;
        }
        (*e).fixed = 0;
        // C rereads e->eu after the first correction may have freed both
        // edge uses. Save the endpoints before modifying the topology.
        let v1 = (*(*(*e).eu).vu).v;
        let v2 = (*(*(*(*e).eu).euMate).vu).v;
        correctTrisAroundVertex(v1);
        correctTrisAroundVertex(v2);
        return 1;
    }
}
unsafe fn correctTrisAroundVertex(mut v: *mut Vertex) {
    unsafe {
        let mut v2: *mut Vertex = std::ptr::null_mut();
        let mut vtemp1: *mut Vertex = std::ptr::null_mut();
        let mut vtemp2: *mut Vertex = std::ptr::null_mut();
        let mut vu: *mut VertexUse = std::ptr::null_mut();
        let mut enew: *mut Edge = std::ptr::null_mut();
        let mut fnew: *mut Face = std::ptr::null_mut();
        let mut fs: *mut Face = std::ptr::null_mut();
        let mut t: *mut Face = std::ptr::null_mut();
        let mut m: *mut Model = std::ptr::null_mut();
        m = (if ((*(*(*v).vu).eu).f != std::ptr::null_mut()) {
            (*(*(*(*v).vu).eu).f).m
        } else {
            (*(*(*(*(*v).vu).eu).euMate).f).m
        });
        vu = (*v).vu;
        loop {
            {
                if !((*(*(*vu).eu).e).fixed != 0) {
                    t = (*(*vu).eu).f;
                    vtemp1 = (*(*(*(*(*vu).eu).euMate).euCCW).vu).v;
                    if inCircumTri((*vtemp1).x, (*vtemp1).y, &*t) != 0 {
                        vtemp2 = (*(*(*(*vu).eu).euCCW).vu).v;
                        v2 = (*(*(*(*vu).eu).euMate).vu).v;
                        if (*m).tDel.is_some() {
                            (*m).tDel.unwrap()(m, t);
                            (*m).tDel.unwrap()(m, (*(*(*vu).eu).euMate).f);
                        }
                        killEdge((*(*vu).eu).e, &raw mut fs);
                        makeEdgeFace(vtemp1, vtemp2, &raw mut enew, &raw mut fnew);
                        circumTri(&mut *fs);
                        circumTri(&mut *fnew);
                        if (*m).tAdd.is_some() {
                            (*m).tAdd.unwrap()(m, fs);
                            (*m).tAdd.unwrap()(m, fnew);
                        }
                        correctTrisAroundVertex(v);
                        correctTrisAroundVertex(v2);
                        break;
                    }
                }
                vu = (*vu).vuNext;
            }
            if !(vu != (*v).vu) {
                break;
            }
        }
    }
}
