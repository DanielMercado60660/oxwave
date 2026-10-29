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

//! Symmetric eigensolvers ported from `cwp/lib/symmeigen.c`.
//! Credits from the C source: "(inspired by Press et. al., 1996)".

/// Find eigenvalues and eigenvectors of a symmetric matrix by Jacobi sweeps.
///
/// C prototype: `void eig_jacobi(float **a, float d[], float **v, int n);`
///
/// Here `n = d.len()` and all indices are zero-based. Both matrices must have
/// at least `n` rows with at least `n` entries each. Only the diagonal and upper
/// triangle of `a` are read; rotations overwrite the upper triangle, leaving
/// its diagonal and lower triangle unchanged. Eigenvectors are the columns of
/// `v`, which is initialized to the identity. Eigenvalues are not sorted.
/// As in C, after 50 unsuccessful sweeps a diagnostic is written to stderr
/// and the current estimates are returned.
pub fn eig_jacobi(a: &mut [Vec<f32>], d: &mut [f32], v: &mut [Vec<f32>]) {
    let n = d.len();
    assert!(a.len() >= n && v.len() >= n);
    for ip in 0..n {
        assert!(a[ip].len() >= n && v[ip].len() >= n);
    }
    let mut b = vec![0.0f32; n];
    let mut z = vec![0.0f32; n];

    // Initialize v to the identity matrix.
    for ip in 0..n {
        for iq in 0..n {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
    }
    for ip in 0..n {
        d[ip] = a[ip][ip];
        b[ip] = d[ip];
        z[ip] = 0.0;
    }

    for i in 1..=50 {
        let mut sm = 0.0f32;
        for ip in 0..n.saturating_sub(1) {
            for iq in ip + 1..n {
                // fabs returns double, including for the compound assignment.
                sm = (sm as f64 + (a[ip][iq] as f64).abs()) as f32;
            }
        }
        if sm == 0.0 {
            return;
        }
        let tresh = if i < 4 {
            (0.2 * sm as f64 / (n as i32).wrapping_mul(n as i32) as f64) as f32
        } else {
            0.0
        };
        for ip in 0..n.saturating_sub(1) {
            for iq in ip + 1..n {
                let g = (100.0 * (a[ip][iq] as f64).abs()) as f32;
                if i > 4
                    && ((d[ip] as f64).abs() + g as f64) as f32 == (d[ip] as f64).abs() as f32
                    && ((d[iq] as f64).abs() + g as f64) as f32 == (d[iq] as f64).abs() as f32
                {
                    a[ip][iq] = 0.0;
                } else if (a[ip][iq] as f64).abs() > tresh as f64 {
                    let mut h = d[iq] - d[ip];
                    let mut t;
                    if ((h as f64).abs() + g as f64) as f32 == (h as f64).abs() as f32 {
                        t = a[ip][iq] / h;
                    } else {
                        let theta = (0.5 * h as f64 / a[ip][iq] as f64) as f32;
                        // theta*theta is float before adding the double 1.0.
                        t = (1.0 / ((theta as f64).abs() + (1.0 + (theta * theta) as f64).sqrt()))
                            as f32;
                        if theta < 0.0 {
                            t = -t;
                        }
                    }
                    // C's integer literal in 1+t*t makes this sum float.
                    let c = (1.0 / ((1.0f32 + t * t) as f64).sqrt()) as f32;
                    let s = t * c;
                    let tau = (s as f64 / (1.0 + c as f64)) as f32;
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;

                    for j in 0..ip {
                        rotate(a, j, ip, j, iq, s, tau);
                    }
                    for j in ip + 1..iq {
                        rotate(a, ip, j, j, iq, s, tau);
                    }
                    for j in iq + 1..n {
                        rotate(a, ip, j, iq, j, s, tau);
                    }
                    for j in 0..n {
                        rotate(v, j, ip, j, iq, s, tau);
                    }
                }
            }
        }
        for ip in 0..n {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    // Ignore an I/O error, as C does for fprintf here.
    use std::io::Write;
    let _ = std::io::stderr().write_all(b"jacobi iteration does not converge\n");
}

// C's ROTATE macro; every intermediate operation has type float.
fn rotate(a: &mut [Vec<f32>], i: usize, j: usize, k: usize, l: usize, s: f32, tau: f32) {
    let g = a[i][j];
    let h = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

/// Sort eigenvalues in descending order, permuting eigenvector columns with them.
///
/// C prototype: `void sort_eigenvalues(float d[], float **v, int n);`
///
/// Here `n = d.len()`, indices are zero-based, and `v` must have at least `n`
/// rows of at least `n` entries. The C selection sort's `>=` comparison chooses
/// the last matching maximum on ties, so equal eigenvalues can change order.
pub fn sort_eigenvalues(d: &mut [f32], v: &mut [Vec<f32>]) {
    let n = d.len();
    assert!(v.len() >= n);
    for row in &v[..n] {
        assert!(row.len() >= n);
    }
    for i in 0..n.saturating_sub(1) {
        let mut k = i;
        let mut p = d[i];
        for j in i + 1..n {
            if d[j] >= p {
                k = j;
                p = d[j];
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for j in 0..n {
                p = v[j][i];
                v[j][i] = v[j][k];
                v[j][k] = p;
            }
        }
    }
}
