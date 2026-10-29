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

//! Port of cwp/lib/abel.c.
//! Credits: Dave Hale and Lydia Deng, Colorado School of Mines, 06/01/90.

const NSE: usize = 9;
/// Precomputed coefficients for an Abel transform.
pub struct AbelTransformer {
    n: usize,
    a: Vec<[f32; NSE]>,
    b0: Vec<[f32; NSE]>,
    b1: Vec<[f32; NSE]>,
}

/// Allocate and precompute an Abel transformer.
/// C prototype: `void *abelalloc (int n);`
/// Transforming a single sample is undefined in C because coefficients are uninitialized.
pub fn abelalloc(n: i32) -> AbelTransformer {
    assert!(n >= 1);
    let n = n as usize;
    let h: [f32; NSE] = [
        1.0,
        0.610926299405048390,
        0.895089852938535935,
        1.34082948787002865,
        2.02532848558443890,
        3.18110895533701843,
        5.90898360396353794,
        77.6000213494180286,
        528.221800846070892,
    ];
    let lambda: [f32; NSE] = [
        0.0,
        -2.08424632126539366,
        -5.78928630565552371,
        -14.6268676854951032,
        -35.0617158334443104,
        -83.3258406398958158,
        -210.358805421311445,
        -6673.64911325382036,
        -34897.7050244132261,
    ];
    let mut a = vec![[0.0; NSE]; n];
    let mut b0 = vec![[0.0; NSE]; n];
    let mut b1 = vec![[0.0; NSE]; n];
    for i in 1..n {
        let fi = (i as f32 as f64 + 1.0) as f32;
        for j in 0..NSE {
            let hj = h[j] as f64;
            let lambdaj = lambda[j] as f64;
            let fi = fi as f64;
            let mut temp = (fi / (fi - 1.0)).powf(lambdaj) as f32;
            a[i][j] = temp;
            temp = (temp as f64 * (fi / (fi - 1.0))) as f32;
            let scale = (2.0 * hj * (fi - 1.0) / ((lambdaj + 1.0) * (lambdaj + 2.0))) as f32;
            b0[i][j] = (scale as f64 * (fi - 1.0 + (lambdaj + 2.0 - fi) * temp as f64)) as f32;
            b1[i][j] =
                (-(scale as f64) * (lambdaj + 1.0 + fi - ((fi as f32) * temp) as f64)) as f32;
        }
    }
    AbelTransformer { n, a, b0, b1 }
}
/// Release an Abel transformer by consuming it.
/// C prototype: `void abelfree (void *at);`
pub fn abelfree(at: AbelTransformer) {
    drop(at);
}
/// Compute the transform into `g`. Both slices must contain at least `n` samples.
/// To transform in place, copy the input before calling this function.
/// C prototype: `void abel (void *at, float f[], float g[]);`
pub fn abel(at: &AbelTransformer, f: &[f32], g: &mut [f32]) {
    let n = at.n;
    assert!(
        n >= 2,
        "the C transform uses uninitialized coefficients for n=1"
    );
    assert!(f.len() >= n && g.len() >= n);
    let mut xi = [0.0f32; NSE];
    let mut fi = f[n - 1];
    g[0] = (0.5 * f[0] as f64 + fi as f64) as f32;
    let mut sum = 0.0f32;
    for j in 0..NSE {
        xi[j] = at.b1[n - 1][j] * fi;
        sum += xi[j];
    }
    g[n - 1] = sum;
    for i in (1..n - 1).rev() {
        let fip1 = fi;
        fi = f[i];
        g[0] += fi;
        sum = 0.0;
        for j in 0..NSE {
            xi[j] = at.a[i][j] * xi[j] + at.b0[i][j] * fip1 + at.b1[i][j] * fi;
            sum += xi[j];
        }
        g[i] = sum;
    }
    g[0] = (g[0] as f64 * 2.0) as f32;
}
