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

//! Clean-room implementation of cwp/lib/svd.c using Golub–Reinsch SVD.
//! Original C credits: compute_svd was adapted from XLISP-STAT 2.1 svdecomp.c,
//! adapted by Luke Tierney and David Betz from Numerical Recipes. Other
//! contributor: Ian Kay, Canadian Geological Survey, Ottawa, Ontario, 1999.
//! The published algorithm is Golub and Reinsch (1971), Handbook for Automatic
//! Computation II, Linear Algebra, pp. 134–151 (also published in EISPACK).
//! svd_sort: Nils Maercklin, GeoForschungsZentrum (GFZ) Potsdam, Germany, 2001.
//!
//! No restricted C function bodies were viewed or used. This follows
//! the published Householder bidiagonalization and implicit-shift QR algorithm.
//! Algorithm reference: <https://www.netlib.org/eispack/svd.f>.
//! Matrices are row-major: `a[row][column]`. Singular vectors occupy columns.
//! The decomposition uses double working scalars and bidiagonal workspace,
//! with single-precision matrices and singular values. Back substitution
//! accumulates in single precision. Differential tests check these boundaries.
//! The QR cancellation step preserves SU's behavior, including its
//! nonstandard results on rank-deficient input.

fn signed_magnitude(magnitude: f64, sign: f64) -> f64 {
    if sign >= 0.0 {
        magnitude.abs()
    } else {
        -magnitude.abs()
    }
}

fn norm(a: f64, b: f64) -> f64 {
    let a = a.abs();
    let b = b.abs();
    if a > b {
        let r = b / a;
        a * (1.0 + r * r).sqrt()
    } else if b != 0.0 {
        let r = a / b;
        b * (1.0 + r * r).sqrt()
    } else {
        0.0
    }
}

fn check_matrix(a: &[Vec<f32>], rows: usize, columns: usize) {
    assert!(
        a.len() >= rows && a[..rows].iter().all(|row| row.len() >= columns),
        "matrix dimensions do not match"
    );
}

/// Decompose `a = U diag(w) Vᵀ`, replacing `a` by the thin left factor U.
/// Singular values are nonnegative and unordered; V contains right vectors.
/// `a` has m rows and n columns, `w` has n elements, and `v` is n by n.
///
/// C prototype: `int compute_svd(float **a, int m, int n, float w[], float **v);`
///
/// Returns 1 on success and 0 if QR iteration fails to converge.
/// Returns 0 and prints the original diagnostic when m < n.
/// Panics on inconsistent array dimensions.
/// Preserves SU's QR cancellation behavior, which can differ from the
/// mathematical SVD for rank-deficient matrices.
pub fn compute_svd(a: &mut [Vec<f32>], w: &mut [f32], v: &mut [Vec<f32>]) -> i32 {
    let m = a.len();
    let n = w.len();
    if m < n {
        oxwave_core::fprintf!(oxwave_core::cio::stderr(), "#rows must be > #cols \n");
        return 0;
    }
    check_matrix(a, m, n);
    check_matrix(v, n, n);
    // SU keeps the off-diagonal workspace (rv1) in double precision.
    let mut off = vec![0.0f64; n];
    let mut g = 0.0f64;
    let mut scale = 0.0f64;
    let mut anorm = 0.0f64;

    // Reduce the matrix to upper bidiagonal form using two Householder
    // reflectors per column. Store the reflectors in a for accumulation below.
    // Black-box comparisons require each matrix update to round to f32 before
    // addition, both here and when accumulating the singular-vector factors.
    for i in 0..n {
        let l = i + 1;
        off[i] = scale * g;
        g = 0.0;
        let mut s = 0.0f64;
        scale = 0.0;
        for k in i..m {
            scale += a[k][i].abs() as f64;
        }
        if scale != 0.0 {
            for k in i..m {
                a[k][i] = (a[k][i] as f64 / scale) as f32;
                s += (a[k][i] as f64) * (a[k][i] as f64);
            }
            let f = a[i][i] as f64;
            g = -signed_magnitude(s.sqrt(), f);
            let h = f * g - s;
            a[i][i] = (f - g) as f32;
            for j in l..n {
                s = 0.0;
                for k in i..m {
                    s += a[k][i] as f64 * a[k][j] as f64;
                }
                let f = s / h;
                for k in i..m {
                    a[k][j] += (f * a[k][i] as f64) as f32;
                }
            }
            for k in i..m {
                a[k][i] = (a[k][i] as f64 * scale) as f32;
            }
        }
        w[i] = (scale * g) as f32;
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if l < n {
            for k in l..n {
                scale += a[i][k].abs() as f64;
            }
            if scale != 0.0 {
                for k in l..n {
                    a[i][k] = (a[i][k] as f64 / scale) as f32;
                    s += (a[i][k] as f64) * (a[i][k] as f64);
                }
                let f = a[i][l] as f64;
                g = -signed_magnitude(s.sqrt(), f);
                let h = f * g - s;
                a[i][l] = (f - g) as f32;
                for k in l..n {
                    off[k] = a[i][k] as f64 / h;
                }
                for j in l..m {
                    s = 0.0;
                    for k in l..n {
                        s += a[j][k] as f64 * a[i][k] as f64;
                    }
                    for k in l..n {
                        a[j][k] += (s * off[k]) as f32;
                    }
                }
                for k in l..n {
                    a[i][k] = (a[i][k] as f64 * scale) as f32;
                }
            }
        }
        let bound = w[i].abs() as f64 + off[i].abs();
        if bound > anorm {
            anorm = bound;
        }
    }

    // Form the right orthogonal factor from the stored row reflectors.
    for i in (0..n).rev() {
        let l = i + 1;
        if l < n {
            if g != 0.0 {
                for j in l..n {
                    v[j][i] = ((a[i][j] as f64 / a[i][l] as f64) / g) as f32;
                }
                for j in l..n {
                    let mut s = 0.0f64;
                    for k in l..n {
                        s += a[i][k] as f64 * v[k][j] as f64;
                    }
                    for k in l..n {
                        v[k][j] += (s * v[k][i] as f64) as f32;
                    }
                }
            }
            for j in l..n {
                v[i][j] = 0.0;
                v[j][i] = 0.0;
            }
        }
        v[i][i] = 1.0;
        g = off[i];
    }

    // Form the thin left orthogonal factor from the column reflectors.
    for i in (0..n).rev() {
        let l = i + 1;
        g = w[i] as f64;
        for j in l..n {
            a[i][j] = 0.0;
        }
        if g != 0.0 {
            g = 1.0 / g;
            for j in l..n {
                let mut s = 0.0f64;
                for k in l..m {
                    s += a[k][i] as f64 * a[k][j] as f64;
                }
                let f = (s / a[i][i] as f64) * g;
                for k in i..m {
                    a[k][j] += (f * a[k][i] as f64) as f32;
                }
            }
            for j in i..m {
                a[j][i] = (a[j][i] as f64 * g) as f32;
            }
        } else {
            for j in i..m {
                a[j][i] = 0.0;
            }
        }
        a[i][i] += 1.0;
    }

    // Diagonalize by implicit-shift QR, chasing the bulge with Givens
    // rotations and applying the same rotations to both singular factors.
    for k in (0..n).rev() {
        for iteration in 0..=30 {
            let mut l = k;
            let mut cancel = false;
            loop {
                if off[l].abs() + anorm == anorm {
                    break;
                }
                if l == 0 {
                    break;
                }
                if w[l - 1].abs() as f64 + anorm == anorm {
                    cancel = true;
                    break;
                }
                l -= 1;
            }
            if cancel {
                let mut s = 1.0f64;
                for i in l..=k {
                    let f = s * off[i];
                    // SU leaves off[i] unchanged here and keeps scanning
                    // after a negligible f, unlike the textbook cancellation.
                    if f.abs() + anorm == anorm {
                        continue;
                    }
                    let g = w[i] as f64;
                    let h = norm(f, g);
                    w[i] = h as f32;
                    let reciprocal = 1.0 / h;
                    let c = g * reciprocal;
                    s = -f * reciprocal;
                    for j in 0..m {
                        let y = a[j][l - 1] as f64;
                        let z = a[j][i] as f64;
                        a[j][l - 1] = (y * c + z * s) as f32;
                        a[j][i] = (z * c - y * s) as f32;
                    }
                }
            }
            let z = w[k] as f64;
            if l == k {
                if z < 0.0 {
                    w[k] = -z as f32;
                    for j in 0..n {
                        v[j][k] = -v[j][k];
                    }
                }
                break;
            }
            if iteration == 30 {
                return 0;
            }
            let mut x = w[l] as f64;
            let y = w[k - 1] as f64;
            let g = off[k - 1];
            let h = off[k];
            // Keep the factored differences of squares: rearranging the shift
            // changes QR deflation and the singular-vector basis returned by C.
            let mut f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            let g = norm(f, 1.0);
            f = ((x - z) * (x + z) + h * (y / (f + signed_magnitude(g, f)) - h)) / x;
            let mut c = 1.0f64;
            let mut s = 1.0f64;
            for j in l..k {
                let i = j + 1;
                let mut g = off[i];
                let mut y = w[i] as f64;
                let h = s * g;
                g *= c;
                let mut z = norm(f, h);
                off[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                let h = y * s;
                y *= c;
                for jj in 0..n {
                    x = v[jj][j] as f64;
                    z = v[jj][i] as f64;
                    v[jj][j] = (x * c + z * s) as f32;
                    v[jj][i] = (z * c - x * s) as f32;
                }
                z = norm(f, h);
                w[j] = z as f32;
                if z != 0.0 {
                    let reciprocal = 1.0 / z;
                    c = f * reciprocal;
                    s = h * reciprocal;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for jj in 0..m {
                    y = a[jj][j] as f64;
                    z = a[jj][i] as f64;
                    a[jj][j] = (y * c + z * s) as f32;
                    a[jj][i] = (z * c - y * s) as f32;
                }
            }
            off[l] = 0.0;
            off[k] = f;
            w[k] = x as f32;
        }
    }
    1
}

/// Solve `U diag(w) Vᵀ x = b` in the least-squares sense.
/// A zero singular value contributes zero to the solution; the caller can
/// truncate small singular values by setting them to zero before this call.
///
/// C prototype: `void svd_backsubstitute(float **u, float w[], float **v,
/// int m, int n, float b[], float x[]);`
///
/// Dimensions: U is m by n, V is n by n, b has m elements, w and x have n.
pub fn svd_backsubstitute(u: &[Vec<f32>], w: &[f32], v: &[Vec<f32>], b: &[f32], x: &mut [f32]) {
    let m = b.len();
    let n = w.len();
    check_matrix(u, m, n);
    check_matrix(v, n, n);
    assert_eq!(x.len(), n, "solution dimension does not match");
    let mut projected = vec![0.0f32; n];
    for j in 0..n {
        let mut s = 0.0f32;
        if w[j] != 0.0 {
            for i in 0..m {
                s += u[i][j] * b[i];
            }
            s /= w[j];
        }
        projected[j] = s;
    }
    for j in 0..n {
        let mut s = 0.0f32;
        for jj in 0..n {
            s += v[j][jj] * projected[jj];
        }
        x[j] = s;
    }
}

/// Permute singular values and columns of U and V using the original sort.
///
/// C prototype: `void svd_sort(float **u, float *w, float **v, int n, int m);`
///
/// **C compatibility:** black-box differential probes show that the original
/// compares through index n (inclusive), while taking only n-1 selection
/// steps. Therefore w and every row of U and V must include an explicit extra
/// column at index n. `n` is the logical number of singular values; m is the
/// number of U rows. V has at least n rows. The padding value participates in
/// the permutation and can change. This makes the original out-of-bounds
/// behavior deterministic and safe, rather than inventing a value for it.
/// For n <= 1 the original does nothing and no extra column is required.
/// To obtain ordinary descending sorting of finite singular values, initialize
/// `w[n]` to `f32::NEG_INFINITY`. `compute_svd` can fill `&mut w[..n]` while
/// leaving that padding value and the extra matrix columns intact.
pub fn svd_sort(u: &mut [Vec<f32>], w: &mut [f32], v: &mut [Vec<f32>], n: usize) {
    if n <= 1 {
        return;
    }
    let m = u.len();
    assert!(
        w.len() > n,
        "svd_sort requires an explicit extra singular value at index n"
    );
    check_matrix(u, m, n + 1);
    check_matrix(v, n, n + 1);
    for i in 0..n - 1 {
        let mut k = i;
        for j in i + 1..=n {
            if w[j] >= w[k] {
                k = j;
            }
        }
        if k != i {
            w.swap(i, k);
            for j in 0..m {
                u[j].swap(i, k);
            }
            for j in 0..n {
                v[j].swap(i, k);
            }
        }
    }
}
