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

//! Port of cwp/lib/windows.c.
//! Author: Potash Corporation, Sascatchewan: Balasz Nemeth, given to CWP 2008.

/// Overwrites `w` with an n-element Hanning window, where n is its length.
///
/// C prototype: `void hanningnWindow(int n, float *w);`
///
/// An empty slice performs no writes, like a nonpositive C length.
/// The slice length must be less than `i32::MAX`, so that C's `n + 1`
/// is representable.
#[allow(non_snake_case)]
pub fn hanningnWindow(w: &mut [f32]) {
    let n = i32::try_from(w.len()).expect("window length exceeds C int");
    assert!(n < i32::MAX, "window length makes n + 1 overflow");
    // C's PI is double, but PI2 is stored as float before the loop.
    let pi2 = (2.0 * 3.141592653589793_f64) as f32;
    for i in 0..n {
        // Multiply and divide in float, then promote for C's double cos.
        let angle = (pi2 * (i + 1) as f32) / (n + 1) as f32;
        w[i as usize] = (0.5 * (1.0 - (angle as f64).cos())) as f32;
    }
}
