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

//! Checked workspace addressing for the pointer aliases in Refl/lib/psvreflect.c
//! and shreflect.c. Credits: Copyright (c) Colorado School of Mines, 2011.
//! C reads outside allocations have no portable result; Rust panics at that access.
use std::cell::RefCell;
use std::collections::HashMap;
thread_local! {static BOUNDS: RefCell<HashMap<usize,(isize,usize)>> = RefCell::new(HashMap::new());}

#[derive(Default)]
pub(super) struct BoundScope {
    previous: Vec<(usize, Option<(isize, usize)>)>,
}
impl BoundScope {
    pub(super) fn track_slice<T>(&mut self, data: &[T]) {
        let p = data.as_ptr() as usize;
        let old = BOUNDS.with(|b| b.borrow_mut().insert(p, (0, data.len())));
        self.previous.push((p, old));
    }
    pub(super) fn track_workspace<T>(&mut self, p: *mut T, n: usize) {
        let p = p as usize;
        let old = BOUNDS.with(|b| b.borrow_mut().insert(p, (-1, n)));
        self.previous.push((p, old));
    }
}
impl Drop for BoundScope {
    fn drop(&mut self) {
        BOUNDS.with(|b| {
            let mut b = b.borrow_mut();
            for (p, old) in self.previous.drain(..).rev() {
                if let Some(n) = old {
                    b.insert(p, n);
                } else {
                    b.remove(&p);
                }
            }
        });
    }
}

/// Resolve an array access only after checking its actual allocation length.
/// The scope keeps registration alive for every buffer and pointer table.
pub(super) unsafe fn checked_at<T>(p: *mut T, index: isize) -> *mut T {
    let (lower, n) = BOUNDS.with(|b| {
        *b.borrow()
            .get(&(p as usize))
            .expect("unregistered reflectivity workspace")
    });
    assert!(
        index >= lower && index < n as isize,
        "C reflectivity array index {index} outside allocation of {n}"
    );
    p.offset(index)
}
