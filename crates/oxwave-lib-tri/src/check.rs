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

//! Port of tri/lib/check.c.
//! Author: Dave Hale, Colorado School of Mines, 07/09/90.

#![allow(unused_assignments, unused_parens, unused_mut)]
use super::*;

/// C: `void badModel(void);`
/// Terminates through C exit(-1), including flushing C streams.
pub fn badModel() -> ! {
    unsafe {
        exit(-1);
    }
}

/// C: `void checkVertexUse (VertexUse *vu);`
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn checkVertexUse(mut vu: *mut VertexUse) {
    unsafe {
        let mut vut: *mut VertexUse = std::ptr::null_mut();
        if ((*vu).v != (*(*vu).vuPrev).v) {
            badModel();
        }
        if ((*vu).v != (*(*vu).vuNext).v) {
            badModel();
        }
        if ((*(*vu).vuNext).vuPrev != vu) {
            badModel();
        }
        if ((*(*vu).vuPrev).vuNext != vu) {
            badModel();
        }
        if ((*(*vu).eu).vu != vu) {
            badModel();
        }
        vut = vu;
        loop {
            {
                if ((*vut).v != (*vu).v) {
                    badModel();
                }
                vut = (*vut).vuNext;
            }
            if !(vut != vu) {
                break;
            }
        }
    }
}

/// C: `void checkEdgeUse (EdgeUse *eu);`
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn checkEdgeUse(mut eu: *mut EdgeUse) {
    unsafe {
        if ((*(*eu).vu).eu != eu) {
            badModel();
        }
        if ((*(*eu).euMate).euMate != eu) {
            badModel();
        }
        if ((*(*eu).euCW).euCCW != eu) {
            badModel();
        }
        if ((*(*eu).euCCW).euCW != eu) {
            badModel();
        }
        if ((*eu).euCW == (*eu).euCCW) {
            badModel();
        }
        if (((*(*eu).e).eu != eu) && ((*(*eu).e).eu != (*eu).euMate)) {
            badModel();
        }
        if (((*(*(*eu).euCW).euCW).euCW != eu) && ((*eu).f != std::ptr::null_mut())) {
            badModel();
        }
        if (((*(*(*eu).euCCW).euCCW).euCCW != eu) && ((*eu).f != std::ptr::null_mut())) {
            badModel();
        }
        checkVertexUse((*eu).vu);
        checkVertexUse((*(*eu).euMate).vu);
    }
}

/// C: `void checkFace (Face *f);`
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn checkFace(mut f: *mut Face) {
    unsafe {
        let mut eu: *mut EdgeUse = std::ptr::null_mut();
        if ((*(*f).fPrev).fNext != f) {
            badModel();
        }
        if ((*(*f).fNext).fPrev != f) {
            badModel();
        }
        if ((*(*f).eu).f != f) {
            badModel();
        }
        eu = (*f).eu;
        loop {
            {
                checkEdgeUse(eu);
                checkEdgeUse((*eu).euMate);
                eu = (*eu).euCW;
            }
            if !(eu != (*f).eu) {
                break;
            }
        }
    }
}

/// C: `void checkModel (Model *m);`
/// # Safety
/// All topology pointers must be valid and uniquely mutable for the call.
pub unsafe fn checkModel(mut m: *mut Model) {
    unsafe {
        let mut f: *mut Face = std::ptr::null_mut();
        if ((*(*m).f).m != m) {
            badModel();
        }
        f = (*m).f;
        loop {
            {
                if ((*f).m != m) {
                    badModel();
                }
                checkFace(f);
                f = (*f).fNext;
            }
            if !(f != (*m).f) {
                break;
            }
        }
    }
}
