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

//! Port of cwp/lib/maxmin.c.
//! Author: Balasz Nemeth: Potash Corporation, given to CWP in 2008.

/// Return the index of the maximum among indices 0, 1, 1 + `inc`, ... .
/// Ties retain the first index visited. An empty slice returns 0, as C does.
/// Comparisons use C's strict ordering: NaNs never replace the current value,
/// and a NaN at index 0 keeps the result at 0.
///
/// C prototype: `int max_index(int n, float *a, int inc);`
///
/// # Panics
/// Panics if the length exceeds `i32::MAX`, or `inc` is not positive when
/// the loop would execute.
pub fn max_index(a: &[f32], inc: i32) -> i32 {
    let n = i32::try_from(a.len()).expect("array length exceeds C int");
    assert!(n <= 1 || inc > 0, "increment must be positive");
    let mut im = 0;
    let mut i = 1;
    while i < n {
        if a[im as usize] < a[i as usize] {
            im = i;
        }
        // C signed overflow is undefined; a wider index ends the loop safely.
        let next = i as i64 + inc as i64;
        if next >= n as i64 {
            break;
        }
        i = next as i32;
    }
    im
}

/// Return the index of the minimum among indices 0, 1, 1 + `inc`, ... .
/// Ties retain the first index visited. An empty slice returns 0, as C does.
/// Comparisons use C's strict ordering: NaNs never replace the current value,
/// and a NaN at index 0 keeps the result at 0.
///
/// C prototype: `int min_index(int n, float *a, int inc);`
///
/// # Panics
/// Panics if the length exceeds `i32::MAX`, or `inc` is not positive when
/// the loop would execute.
pub fn min_index(a: &[f32], inc: i32) -> i32 {
    let n = i32::try_from(a.len()).expect("array length exceeds C int");
    assert!(n <= 1 || inc > 0, "increment must be positive");
    let mut im = 0;
    let mut i = 1;
    while i < n {
        if a[im as usize] > a[i as usize] {
            im = i;
        }
        let next = i as i64 + inc as i64;
        if next >= n as i64 {
            break;
        }
        i = next as i32;
    }
    im
}
