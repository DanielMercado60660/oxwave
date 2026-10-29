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

//! Port of `cwp/lib/wrap.c`.
//! Author: Potash Corporation: Balasz Nemeth, given to CWP 2008.

/// Wrap an array by `f` elements, moving positive shifts to the right.
///
/// C prototype: `void wrapArray(void *base, size_t nmemb, size_t size, int f)`.
///
/// The array must be nonempty and its length must fit in a C `int`, as
/// required by the original function's modulo and indexing arithmetic.
#[allow(non_snake_case)]
pub fn wrapArray<T: Copy>(base: &mut [T], f: i32) {
    assert!(!base.is_empty() && base.len() <= i32::MAX as usize);
    let nmemb = base.len() as i32;
    let mut shift = f % nmemb;
    if shift == 0 {
        return;
    }
    if shift > 0 {
        shift = -nmemb + shift;
    }
    let split = (-shift) as usize;
    let mut base2 = Vec::with_capacity(base.len());
    base2.extend_from_slice(&base[split..]);
    base2.extend_from_slice(&base[..split]);
    base.copy_from_slice(&base2);
}
