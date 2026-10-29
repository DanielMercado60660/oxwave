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

//! Differential tests against original `cwp/lib/symmeigen.c`.
//! Credits from the C source: "(inspired by Press et. al., 1996)".

use oxwave_cref::XorShift;
use oxwave_symmetric_eigen as rust;

extern "C" {
    fn eig_jacobi(a: *mut *mut f32, d: *mut f32, v: *mut *mut f32, n: i32);
    fn sort_eigenvalues(d: *mut f32, v: *mut *mut f32, n: i32);
}

const GUARD: f32 = -987654.0;

// Allocate real padding, rather than forming C's traditional pointer-before-
// allocation. C addresses only indices 1..=n; Rust addresses 0..n.
struct CMatrix {
    rows: Vec<Vec<f32>>,
    pointers: Vec<*mut f32>,
}

impl CMatrix {
    fn new(a: &[Vec<f32>]) -> Self {
        let n = a.len();
        let mut rows = vec![vec![GUARD; n + 2]; n + 2];
        for i in 0..n {
            rows[i + 1][1..=n].copy_from_slice(&a[i]);
        }
        let pointers = rows.iter_mut().map(|row| row.as_mut_ptr()).collect();
        Self { rows, pointers }
    }

    fn as_mut_ptr(&mut self) -> *mut *mut f32 {
        self.pointers.as_mut_ptr()
    }

    fn check(&self, a: &[Vec<f32>], label: &str) {
        let n = a.len();
        for i in 0..n + 2 {
            for j in 0..n + 2 {
                let expected = if i > 0 && i <= n && j > 0 && j <= n {
                    a[i - 1][j - 1]
                } else {
                    GUARD
                };
                assert_eq!(
                    self.rows[i][j].to_bits(),
                    expected.to_bits(),
                    "{label}: n={n}, padded index ({i},{j}), C={}, Rust={expected}",
                    self.rows[i][j]
                );
            }
        }
    }
}

fn check_eig(input: Vec<Vec<f32>>) {
    let n = input.len();
    let mut a = input.clone();
    let mut d = vec![GUARD; n];
    let mut v = vec![vec![GUARD; n]; n];
    let mut ca = CMatrix::new(&a);
    let mut cv = CMatrix::new(&v);
    let mut cd = vec![GUARD; n + 2];
    unsafe { eig_jacobi(ca.as_mut_ptr(), cd.as_mut_ptr(), cv.as_mut_ptr(), n as i32) };
    rust::eig_jacobi(&mut a, &mut d, &mut v);

    // Bitwise comparisons are stricter than the requested 1e-5 tolerance.
    for i in 0..n {
        assert_eq!(
            d[i].to_bits(),
            cd[i + 1].to_bits(),
            "eigenvalue n={n}, i={i}"
        );
    }
    assert_eq!(cd[0].to_bits(), GUARD.to_bits());
    assert_eq!(cd[n + 1].to_bits(), GUARD.to_bits());
    ca.check(&a, "Jacobi matrix");
    cv.check(&v, "Jacobi eigenvectors");

    // The algorithm must not touch the diagonal or lower triangle.
    for i in 0..n {
        for j in 0..=i {
            assert_eq!(a[i][j].to_bits(), input[i][j].to_bits());
        }
    }

    // Exercise sorting of actual eigenpairs, too.
    unsafe { sort_eigenvalues(cd.as_mut_ptr(), cv.as_mut_ptr(), n as i32) };
    rust::sort_eigenvalues(&mut d, &mut v);
    for i in 0..n {
        assert_eq!(d[i].to_bits(), cd[i + 1].to_bits());
        if i > 0 {
            assert!(d[i - 1] >= d[i]);
        }
    }
    cv.check(&v, "sorted Jacobi eigenvectors");
}

fn symmetric(rng: &mut XorShift, n: usize, scale: f32, mode: usize) -> Vec<Vec<f32>> {
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i..n {
            let x = match mode {
                0 => 0.0,
                1 if i != j => 0.0,
                2 if i == j => -scale,
                3 if i != j && j != i + 1 => 0.0,
                4 => rng.range_i32(-3, 3) as f32 * scale,
                5 => -rng.uniform_f32(0.0, 1.0) * scale,
                6 if i != j => rng.uniform_f32(-1.0, 1.0) * scale * 1e-7,
                _ => rng.uniform_f32(-1.0, 1.0) * scale,
            };
            a[i][j] = x;
            a[j][i] = x;
        }
    }
    a
}

#[test]
fn eig_jacobi_matches_c() {
    let mut rng = XorShift::new(0x5e19_2026);
    let sizes = [1, 2, 3, 5, 7, 8, 11, 16, 17, 32];
    let scales = [1e-30, 1e-15, 1e-4, 1.0, 100.0, 1e15, 1e30];
    for case in 0..4096 {
        let n = sizes[rng.next_u32() as usize % sizes.len()];
        let scale = scales[rng.next_u32() as usize % scales.len()];
        check_eig(symmetric(&mut rng, n, scale, case % 8));
    }
    // Large primes and powers of two, with dense and structured matrices.
    for n in [64, 97, 128] {
        for mode in [0, 1, 3, 7] {
            check_eig(symmetric(&mut rng, n, 1.0, mode));
        }
    }
    check_eig(vec![]);
    check_eig(vec![vec![-0.0]]);
    check_eig(vec![vec![f32::MAX, 0.0], vec![0.0, -f32::MAX]]);
    for scale in [f32::from_bits(1), 1e-40, f32::MIN_POSITIVE] {
        for n in [1, 2, 3, 5, 8] {
            check_eig(symmetric(&mut rng, n, scale, 7));
        }
    }
    // Diagonal separation makes the small-off-diagonal shortcut applicable.
    check_eig(vec![vec![1e20, 1e-20], vec![1e-20, -1e20]]);
    // Lower triangle is ignored, even if it differs from the upper triangle.
    check_eig(vec![vec![1.0, 0.25], vec![12345.0, -2.0]]);
}

fn check_sort(mut d: Vec<f32>, mut v: Vec<Vec<f32>>) {
    let n = d.len();
    let mut cd = vec![GUARD; n + 2];
    cd[1..n + 1].copy_from_slice(&d);
    let mut cv = CMatrix::new(&v);
    unsafe { sort_eigenvalues(cd.as_mut_ptr(), cv.as_mut_ptr(), n as i32) };
    rust::sort_eigenvalues(&mut d, &mut v);
    for i in 0..n {
        assert_eq!(d[i].to_bits(), cd[i + 1].to_bits(), "sort n={n}, i={i}");
    }
    assert_eq!(cd[0].to_bits(), GUARD.to_bits());
    assert_eq!(cd[n + 1].to_bits(), GUARD.to_bits());
    cv.check(&v, "sort eigenvectors");
}

#[test]
fn sort_eigenvalues_matches_c() {
    let mut rng = XorShift::new(0x5077_e19e);
    let sizes = [0, 1, 2, 3, 5, 7, 8, 11, 16, 31, 32, 64, 127, 128];
    for case in 0..4096 {
        let n = sizes[rng.next_u32() as usize % sizes.len()];
        let mut d = Vec::with_capacity(n);
        for i in 0..n {
            d.push(match case % 8 {
                0 => 0.0,
                1 => rng.range_i32(-4, 4) as f32, // many ties
                2 => rng.uniform_f32(-1e30, 1e30),
                3 => rng.uniform_f32(-1e-30, 1e-30),
                4 => -(i as f32),
                5 => i as f32,
                6 => rng.any_f32(), // infinities, NaNs, signed zero, subnormals
                _ => rng.uniform_f32(-1000.0, 1000.0),
            });
        }
        // Arbitrary column payloads make every permutation observable.
        let v = (0..n)
            .map(|_| (0..n).map(|_| rng.any_f32()).collect())
            .collect();
        check_sort(d, v);
    }
    for n in [257, 512] {
        let d = (0..n).map(|_| rng.uniform_f32(-100.0, 100.0)).collect();
        let v = (0..n)
            .map(|i| (0..n).map(|j| (i * n + j) as f32).collect())
            .collect();
        check_sort(d, v);
    }
    let d = vec![
        f32::NAN,
        f32::INFINITY,
        -0.0,
        0.0,
        -f32::INFINITY,
        f32::from_bits(1),
        -f32::from_bits(1),
        2.0,
        2.0,
    ];
    let n = d.len();
    check_sort(
        d,
        (0..n)
            .map(|i| (0..n).map(|j| (i * n + j) as f32).collect())
            .collect(),
    );
}

// A mathematical check independent of the C oracle guards the row/column
// convention used by the differential adapter.
#[test]
fn eigenpairs_have_small_residuals() {
    let mut rng = XorShift::new(0xace1);
    for n in [1, 2, 3, 7, 16, 31] {
        let input = symmetric(&mut rng, n, 1.0, 7);
        let mut a = input.clone();
        let mut d = vec![0.0; n];
        let mut v = vec![vec![0.0; n]; n];
        rust::eig_jacobi(&mut a, &mut d, &mut v);
        for k in 0..n {
            for i in 0..n {
                let mut av = 0.0f64;
                for j in 0..n {
                    av += input[i][j] as f64 * v[j][k] as f64;
                }
                let vd = v[i][k] as f64 * d[k] as f64;
                assert!((av - vd).abs() <= 1e-5 * n as f64);
                let mut dot = 0.0f64;
                for j in 0..n {
                    dot += v[j][i] as f64 * v[j][k] as f64;
                }
                assert!((dot - if i == k { 1.0 } else { 0.0 }).abs() <= 1e-5);
            }
        }
    }
}
