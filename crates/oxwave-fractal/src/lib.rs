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

//! Port of `par/lib/fractal.c`.
//! Author: Balazs Nemeth, Potash Corporation of Saskatoon Saskatchewan;
//! code donated to CWP c. 2008.

use oxwave_linear_fit::linfit;

/// Return the distance between two points.
///
/// C prototype: `float distance(float x0, float y0, float x1, float y1)`.
pub fn distance(x0: f32, y0: f32, x1: f32, y1: f32) -> f32 {
    (((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)) as f64).sqrt() as f32
}

/// Compute the fractal dimension of an array via the divider method.
///
/// C prototype: `float hausdorff_dimension(float *ar, int n, int minl, int maxl, int dl)`.
/// The array length supplies `n`. Divider lengths and the step must be positive,
/// and `maxl` must be at least `minl`, as required by the C algorithm.
pub fn hausdorff_dimension(ar: &[f32], minl: i32, maxl: i32, dl: i32) -> f32 {
    assert!(minl > 0 && dl > 0 && maxl >= minl);
    let n = i32::try_from(ar.len()).expect("array length exceeds C int range");
    let nm = (maxl - minl) / dl;
    let mut ll = vec![0.0_f32; nm as usize];
    let mut ls = vec![0.0_f32; nm as usize];
    let mut l = minl;
    for i in 0..nm as usize {
        let mut sum = 0.0_f64;
        let mut i1 = 0_i32;
        let mut i2 = l;
        while i2 < n {
            sum += distance(i1 as f32, ar[i1 as usize], i2 as f32, ar[i2 as usize]) as f64;
            i1 = i2;
            i2 = i2.wrapping_add(l);
        }
        ll[i] = (l as f32 as f64).ln() as f32;
        ls[i] = sum.ln() as f32;
        l = l.wrapping_add(dl);
    }
    let (mut a, mut b, mut siga, mut sigb, mut chi2, mut q) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    linfit(
        &ll,
        &ls,
        &[],
        0,
        &mut a,
        &mut b,
        &mut siga,
        &mut sigb,
        &mut chi2,
        &mut q,
    );
    (1.0_f64 - b as f64) as f32
}
