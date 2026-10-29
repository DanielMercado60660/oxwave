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

//! Ports of cwp/lib/holbergd1.c, cwp/lib/mkdiff.c, and cwp/lib/mkhdiff.c.
//! Author: Dave Hale, Colorado School of Mines, 06/06/91 and 06/02/89.

use oxwave_core::err;
use oxwave_vandermonde::vanded;

/// Compute Holberg's first-derivative convolution filter.
/// C prototype: `void holbergd1 (float e, int n, float d[]);`
/// The output length supplies `n`, which must be 2, 4, 6, or 8.
pub fn holbergd1(e: f32, d: &mut [f32]) {
    const B: [[[f32; 6]; 4]; 4] = [
        [[1.0, 1.0, 0.0, 0.0, 0.0, 0.0], [0.0; 6], [0.0; 6], [0.0; 6]],
        [
            [1.125, 0.4330, -0.4583, 0.2566, 0.0, 0.0],
            [-0.04167, -0.1443, -0.1806, -0.0855, 0.0, 0.0],
            [0.0; 6],
            [0.0; 6],
        ],
        [
            [1.172, 0.2742, -0.3006, 0.2637, -0.2391, 0.0],
            [-0.06510, -0.1371, -0.0100, 0.1077, -0.0086, 0.0],
            [0.004688, 0.0274, 0.0661, 0.0826, 0.0530, 0.0],
            [0.0; 6],
        ],
        [
            [1.1963, 0.1987, -0.2195, 0.2096, -0.2038, 0.1779],
            [-0.07975, -0.1192, 0.03929, 0.04049, -0.05841, 0.04389],
            [0.009570, 0.03975, 0.04853, -0.00967, -0.05507, -0.01117],
            [-0.0006975, -0.00568, -0.02014, -0.04039, -0.04937, -0.03625],
        ],
    ];
    let n = d.len();
    if !matches!(n, 2 | 4 | 6 | 8) {
        err!("holbergd1: invalid output length");
    }
    // The C computes values that it subsequently overwrites; for n=8 those
    // unused values involve an out-of-bounds read of its table.
    for i in 0..n / 2 {
        d[i] = 0.0;
        for j in 0..6 {
            d[i] = (d[i] as f64
                + B[n / 2 - 1][i][j] as f64 * (e as f64).powf(2.0 * j as f64 / n as f64))
                as f32;
        }
    }
    for i in n / 2..n {
        d[i] = -d[i - n / 2];
    }
    for i in 0..n / 2 {
        d[i] = -d[n - i - 1];
    }
}

/// Make a Taylor-series approximation to a derivative.
/// C prototype: `void mkdiff (int n, float a, float h, int l, int m, float d[]);`
/// Requires `0 <= n <= m-l`, at most 100 coefficients, and sufficient output space.
pub fn mkdiff(n: i32, a: f32, h: f32, l: i32, m: i32, d: &mut [f32]) {
    let count = i64::from(m) - i64::from(l) + 1;
    if !(1..=100).contains(&count)
        || n < 0
        || i64::from(n) >= count
        || d.len() < count as usize
        || l == i32::MIN
    {
        err!("mkdiff: invalid arguments");
    }
    let nc = count as usize;
    let mut v = vec![0.0; nc];
    let mut b = vec![0.0; nc];
    for i in 0..nc {
        v[i] = (-(l + i as i32) as f32 - a) as f64;
        b[i] = 0.0;
    }
    let mut fn_ = 1.0;
    let mut hn = 1.0;
    for i in 1..=n {
        fn_ *= i as f64;
        hn *= h as f64;
    }
    b[n as usize] = fn_ / hn;
    let mut x = vec![0.0; nc];
    vanded(&v, &b, &mut x);
    for i in 0..nc {
        d[i] = x[i] as f32;
    }
}

/// Compute the tapered bandlimited half-differentiator.
/// C prototype: `void mkhdiff (float h, int l, float d[]);`
/// `l` must be nonnegative and `d` must hold at least `1+2*l` values.
/// As in C, `l = 0` yields a NaN coefficient because the integration and
/// taper divide by zero; use a positive half-length for a usable filter.
pub fn mkhdiff(h: f32, l: i32, d: &mut [f32]) {
    if l < 0 || l > i32::MAX / 8 || d.len() < 1 + 2 * l as usize {
        err!("mkhdiff: invalid arguments");
    }
    const PI: f64 = 3.141592653589793;
    let nw = 8 * l;
    let dw = (PI / nw as f64) as f32;
    for j in -l..=l {
        let t = j as f32;
        let mut sum = (0.0_f64.sqrt() * (0.0_f64.cos() - 0.0_f64.sin())
            - PI.sqrt() * ((PI * t as f64).cos() - (PI * t as f64).sin()))
            as f32;
        let mut iw = 1;
        let mut w = dw;
        while iw <= nw - 1 {
            let wt = (w * t) as f64;
            let wdt = ((w + dw) * t) as f64;
            sum = (sum as f64
                + (4.0 * (w as f64).sqrt() * (wt.cos() - wt.sin())
                    + 2.0 * ((w + dw) as f64).sqrt() * (wdt.cos() - wdt.sin())))
                as f32;
            iw += 2;
            w = (w as f64 + 2.0 * dw as f64) as f32;
        }
        d[(l + j) as usize] =
            ((2.0 / h as f64).sqrt() / (2.0 * PI) * PI / (3.0 * nw as f64) * sum as f64) as f32;
    }
    for j in -l..=l {
        let i = (l + j) as usize;
        d[i] = (d[i] as f64 * (0.54 + 0.46 * (PI * j as f64 / l as f64).cos())) as f32;
    }
}
