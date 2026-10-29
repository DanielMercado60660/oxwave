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

//! Port of cwp/lib/xindex.c.
//! Author: Dave Hale, Colorado School of Mines, 12/25/89.

/// Find the index of `x` in a nonempty monotonic array using a previous index.
///
/// C prototype: `void xindex(int nx, float ax[], float x, int *index)`.
/// For increasing arrays, returns the last index with value at most `x`;
/// otherwise, returns the last index with value at least `x`. If no value
/// qualifies, returns zero. Constant arrays use the decreasing-array branch.
/// The input index is a search hint and is clamped to the array bounds.
/// Comparisons follow C semantics, including for NaN queries.
///
/// # Panics
/// Panics if the slice is empty or its length exceeds `i32::MAX`. The C
/// routine requires a nonempty array whose length fits in a C `int`.
pub fn xindex(ax: &[f32], x: f32, index: &mut i32) {
    assert!(!ax.is_empty() && ax.len() <= i32::MAX as usize);
    let nx = ax.len() as i32;
    // Initialize lower and upper indices and step.
    let mut lower = *index;
    if lower < 0 {
        lower = 0;
    }
    if lower >= nx {
        lower = nx - 1;
    }
    let mut upper = lower + 1;
    let mut step: i32 = 1;
    // If x values are increasing.
    if ax[(nx - 1) as usize] > ax[0] {
        // Find indices such that ax[lower] <= x < ax[upper].
        while lower > 0 && ax[lower as usize] > x {
            upper = lower;
            lower = lower.wrapping_sub(step);
            step = step.wrapping_add(step);
        }
        if lower < 0 {
            lower = 0;
        }
        while upper < nx && ax[upper as usize] <= x {
            lower = upper;
            upper = upper.wrapping_add(step);
            step = step.wrapping_add(step);
        }
        if upper > nx {
            upper = nx;
        }
        // Find index via bisection.
        loop {
            let middle = lower.wrapping_add(upper) >> 1;
            if middle == lower {
                break;
            }
            if x >= ax[middle as usize] {
                lower = middle;
            } else {
                upper = middle;
            }
        }
    } else {
        // Find indices such that ax[lower] >= x > ax[upper].
        while lower > 0 && ax[lower as usize] < x {
            upper = lower;
            lower = lower.wrapping_sub(step);
            step = step.wrapping_add(step);
        }
        if lower < 0 {
            lower = 0;
        }
        while upper < nx && ax[upper as usize] >= x {
            lower = upper;
            upper = upper.wrapping_add(step);
            step = step.wrapping_add(step);
        }
        if upper > nx {
            upper = nx;
        }
        // Find index via bisection.
        loop {
            let middle = lower.wrapping_add(upper) >> 1;
            if middle == lower {
                break;
            }
            if x <= ax[middle as usize] {
                lower = middle;
            } else {
                upper = middle;
            }
        }
    }
    // Return lower index.
    *index = lower;
}
