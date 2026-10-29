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

//! Port of `cwp/lib/conjugate_gradient.c`.
//! Credits: CWP: John Stockwell, May 2014.
//! Based on the subroutine cg() which appeared in (Claerbout, 1985).

/// Compute the scalar product of two equally sized vectors.
///
/// C prototype: `float inner_product(float *x, float *y, int n);`
///
/// Accumulates in `f32`, in index order. Empty slices return positive zero.
/// Panics if the slice lengths differ.
pub fn inner_product(x: &[f32], y: &[f32]) -> f32 {
    assert_eq!(x.len(), y.len());
    let mut prod = 0.0_f32;
    for i in 0..x.len() {
        let tmpx = x[i];
        let tmpy = y[i];
        prod += tmpx * tmpy;
    }
    prod
}

/// Solve the least-squares system with the original conjugate gradient routine.
///
/// C prototype: `void simple_conj_gradient(int n, float *x, int m, float *b, float **a, int niter);`
///
/// `a` contains `b.len()` rows of `x.len()` elements. The initial contents
/// of `x` are discarded, as in C. The C routine assumes `b.len() >= x.len()`.
/// Nonpositive `niter` returns an all-zero solution. Degenerate systems can
/// yield NaNs because the C routine does not stop on a zero gradient.
/// Panics if the matrix dimensions do not match the vector lengths.
pub fn simple_conj_gradient(x: &mut [f32], b: &[f32], a: &[&[f32]], niter: i32) {
    let n = x.len();
    let m = b.len();
    assert_eq!(a.len(), m);
    assert!(a.iter().all(|row| row.len() == n));
    let mut s = vec![0.0_f32; n];
    let mut r = b.to_vec();
    let mut g = vec![0.0_f32; n];
    let mut gamma = vec![0.0_f32; m];
    let mut sigma = vec![0.0_f32; m];
    x.fill(0.0);
    for iter in 0..niter {
        for i in 0..n {
            g[i] = 0.0;
            for j in 0..m {
                g[i] += r[j] * a[j][i];
            }
        }
        for j in 0..m {
            gamma[j] = 0.0;
            for i in 0..n {
                gamma[j] += a[j][i] * g[i];
            }
        }
        let (alpha, beta) = if iter == 0 {
            let gambyr = inner_product(&gamma, &r);
            let gamma2 = inner_product(&gamma, &gamma);
            (gambyr / gamma2, 0.0)
        } else {
            let gamma2 = inner_product(&gamma, &gamma);
            let sigma2 = inner_product(&sigma, &sigma);
            let gambysig = inner_product(&gamma, &sigma);
            // CONJ_EPS is an unsuffixed C double: subtract in float, add
            // the epsilon in double, then round back to float on assignment.
            let determ = ((gamma2 * sigma2 - gambysig * gambysig) as f64 + 1e-30) as f32;
            let gambyr = inner_product(&gamma, &r);
            let sigbyr = inner_product(&sigma, &r);
            (
                (sigma2 * gambyr - gambysig * sigbyr) / determ,
                (gamma2 * sigbyr - gambysig * gambyr) / determ,
            )
        };
        for i in 0..n {
            s[i] = alpha * g[i] + beta * s[i];
        }
        for i in 0..m {
            sigma[i] = alpha * gamma[i] + beta * sigma[i];
        }
        for i in 0..n {
            x[i] += s[i];
        }
        for i in 0..m {
            r[i] -= sigma[i];
        }
    }
}
