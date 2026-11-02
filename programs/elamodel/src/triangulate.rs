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

//! ELAMODEL's vertex insertion, from tri/lib/create.c.
//! Author: Dave Hale, Colorado School of Mines, 07/09/90.

use oxwave_lib_tri::*;
use std::ptr::null_mut;

// Keep this insertion adapter local to elamodel: unlike ordinary triangular
// inputs, its finely smoothed curves can create two-edge faces. The shared
// insertion routine leaves a freed edge-use in such a face after a swap.
// All other topology operations still use the shared tri library.
//
// Safety: m and all its topology links must be valid and uniquely mutable.
pub(super) unsafe fn addVertexToModel(m: *mut Model, mut x: f32, mut y: f32) -> *mut Vertex {
    unsafe {
        // Find existing vertex nearest to new vertex.
        let mut v = nearestVertexInModel(m, null_mut(), x, y);
        if ((((*v).x - x) * ((*v).x - x) + ((*v).y - y) * ((*v).y - y)) as f64).sqrt()
            < (*m).eps as f64
        {
            return null_mut();
        }

        // Find triangle that contains new vertex.
        let mut f = if !(*(*(*v).vu).eu).f.is_null() {
            (*(*(*v).vu).eu).f
        } else {
            (*(*(*(*v).vu).eu).euMate).f
        };
        f = insideTriInModel(m, f, x, y);

        // If the edge nearest new vertex is a fixed edge and very nearby.
        let mut e = nearestEdgeInModel(m, (*(*f).eu).e, x, y);
        let mut fixed = false;
        let mut v1fixed = null_mut();
        let mut v2fixed = null_mut();
        let mut eb = null_mut();
        if (*e).fixed != 0 && distanceToEdge(&*e, x, y) < (*m).eps {
            v1fixed = (*(*(*e).eu).vu).v;
            v2fixed = (*(*(*(*e).eu).euCW).vu).v;
            (*e).fixed = 0;
            projectToEdge(&*e, &mut x, &mut y);
            fixed = true;
            if (*(*e).eu).f.is_null() || (*(*(*e).eu).euMate).f.is_null() {
                eb = e;
            }
        }

        if let Some(tDel) = (*m).tDel {
            tDel(m, f);
        }

        // Make edges from new vertex to vertices of this triangle.
        let mut vu = (*(*f).eu).vu;
        let v1 = (*vu).v;
        let v2 = (*(*(*(*vu).eu).euCW).vu).v;
        let v3 = (*(*(*(*vu).eu).euCCW).vu).v;
        let f1 = f;
        let (mut e1, mut e2, mut e3) = (null_mut(), null_mut(), null_mut());
        let (mut f2, mut f3) = (null_mut(), null_mut());
        makeEdgeVertex(v1, x, y, f1, &raw mut e1, &raw mut v);
        makeEdgeFace(v2, v, &raw mut e2, &raw mut f2);
        makeEdgeFace(v3, v, &raw mut e3, &raw mut f3);
        if let Some(tAdd) = (*m).tAdd {
            tAdd(m, f1);
            tAdd(m, f2);
            tAdd(m, f3);
        }

        // Compute circumcircles of the three new triangles.
        vu = (*v).vu;
        loop {
            circumTri(&mut *(*(*vu).eu).f);
            vu = (*vu).vuNext;
            if vu == (*v).vu {
                break;
            }
        }

        // Loop over vertex-uses (triangles) in clockwise direction.
        vu = (*v).vu;
        let mut first = true;
        loop {
            let eu = (*(*vu).eu).euCW;
            let euo = (*eu).euMate;
            e = (*eu).e;
            if !(*euo).f.is_null() && (*(*euo).e).fixed == 0 && inCircumTri(x, y, &*(*euo).f) != 0 {
                let vo = (*(*(*euo).euCCW).vu).v;
                if let Some(tDel) = (*m).tDel {
                    tDel(m, (*(*e).eu).f);
                    tDel(m, (*(*(*e).eu).euMate).f);
                }

                // A two-edge face makes killEdge's eu->euCW->euCW equal
                // to eu itself. Save identities before deletion; never read
                // the retired edge-uses. Explicitly select the corresponding
                // replacement edge-use, preserving the native model produced
                // when SU's allocator reuses the retired allocations.
                let old_eu = (*e).eu;
                let old_eum = (*old_eu).euMate;
                let survivor = (*old_eu).f;
                let mut fs = null_mut();
                killEdge(e, &raw mut fs);
                makeEdgeFace(v, vo, &raw mut e, &raw mut f);
                if (*survivor).eu == old_eu {
                    (*survivor).eu = (*e).eu;
                } else if (*survivor).eu == old_eum {
                    (*survivor).eu = (*(*e).eu).euMate;
                }

                if let Some(tAdd) = (*m).tAdd {
                    tAdd(m, (*(*e).eu).f);
                    tAdd(m, (*(*(*e).eu).euMate).f);
                }
                circumTri(&mut *(*(*e).eu).f);
                circumTri(&mut *(*(*(*e).eu).euMate).f);
            } else {
                vu = (*(*(*(*vu).eu).euCCW).euMate).vu;
                first = false;
            }
            if vu == (*v).vu && !first {
                break;
            }
        }

        // If new vertex was added on an edge that was fixed.
        if fixed {
            vu = (*v).vu;
            loop {
                let eu = (*vu).eu;
                if (*(*(*eu).euCCW).vu).v == v1fixed {
                    (*(*(*eu).euCCW).e).fixed = 1;
                }
                if (*(*(*eu).euCW).vu).v == v2fixed {
                    (*(*eu).e).fixed = 1;
                }
                vu = (*vu).vuNext;
                if vu == (*v).vu {
                    break;
                }
            }
        }
        if !eb.is_null() {
            killBoundaryEdge(eb);
        }
        if let Some(vAdd) = (*m).vAdd {
            vAdd(m, v);
        }
        v
    }
}
