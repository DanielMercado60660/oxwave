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

//! Differential tests against `cwp/lib/conjugate_gradient.c`.
//! Credits: CWP: John Stockwell, May 2014.
//! Based on the subroutine cg() which appeared in (Claerbout, 1985).

use oxwave_cref::XorShift;
use oxwave_f_conjugate_gradient as rust;

extern "C" {
    // float inner_product(float *x, float *y, int n);
    fn inner_product(x: *mut f32, y: *mut f32, n: i32) -> f32;
    // void simple_conj_gradient(int n, float *x, int m, float *b, float **a, int niter);
    fn simple_conj_gradient(n: i32, x: *mut f32, m: i32, b: *mut f32, a: *mut *mut f32, niter: i32);
}

fn sample(rng: &mut XorShift) -> f32 {
    ((rng.next_f64() * 2.0 - 1.0) * 2.0) as f32
}

// This build of C has FP contraction disabled. Require identical finite
// results (stronger than 1e-5 relative); NaN payloads/signs are not portable.
fn compare(actual: f32, expected: f32, context: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{context}: expected NaN, got {actual}");
    } else {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{context}: actual {actual:?}, expected {expected:?}"
        );
    }
}

fn check_inner_product(x: &[f32], y: &[f32], context: &str) {
    let mut c_x = x.to_vec();
    let mut c_y = y.to_vec();
    // C reads exactly x.len() initialized elements from each vector.
    let c = unsafe { inner_product(c_x.as_mut_ptr(), c_y.as_mut_ptr(), x.len() as i32) };
    compare(rust::inner_product(x, y), c, context);
}

fn check_solver(
    matrix: &[Vec<f32>],
    b: &[f32],
    initial: &[f32],
    niter: i32,
    context: &str,
) -> Vec<f32> {
    let mut c_matrix = matrix.to_vec();
    let mut c_b = b.to_vec();
    let mut ptrs: Vec<*mut f32> = c_matrix.iter_mut().map(|row| row.as_mut_ptr()).collect();
    let mut c = initial.to_vec();
    let mut actual = initial.to_vec();
    // All m row pointers refer to initialized n-element buffers that remain
    // alive throughout the call, as do the m-element b and n-element x.
    unsafe {
        simple_conj_gradient(
            c.len() as i32,
            c.as_mut_ptr(),
            c_b.len() as i32,
            c_b.as_mut_ptr(),
            ptrs.as_mut_ptr(),
            niter,
        );
    }
    let rows: Vec<&[f32]> = matrix.iter().map(Vec::as_slice).collect();
    rust::simple_conj_gradient(&mut actual, b, &rows, niter);
    for (i, (&got, &expected)) in actual.iter().zip(&c).enumerate() {
        compare(got, expected, &format!("{context}, x[{i}]"));
    }
    if niter <= 0 {
        assert!(actual.iter().all(|v| v.to_bits() == 0), "{context}");
    }
    actual
}

#[test]
fn inner_product_random_matches_c() {
    let mut rng = XorShift::new(0x12345678);
    let sizes = [0, 1, 2, 3, 7, 16, 31, 64, 127, 256, 1024, 4096];
    let scales = [0.0, 1e-20, 1e-15, 1e-4, 1.0, 1e4, 1e15, 1e20];
    for case in 0..6600 {
        let n = sizes[case % sizes.len()];
        let sx = scales[(case / sizes.len()) % scales.len()];
        let sy = scales[(case / 7) % scales.len()];
        let x: Vec<_> = (0..n).map(|_| sample(&mut rng) * sx).collect();
        let y: Vec<_> = (0..n).map(|_| sample(&mut rng) * sy).collect();
        check_inner_product(&x, &y, &format!("random case {case}, n={n}"));
    }
}

#[test]
fn inner_product_cancellation_and_ieee_edges_match_c() {
    let cases = [
        (vec![], vec![]),
        (vec![-0.0, 0.0], vec![1.0, -1.0]),
        (vec![1e20, 1.0, -1e20], vec![1.0; 3]),
        (vec![1e20, -1e20, 1.0], vec![1.0; 3]),
        (vec![f32::from_bits(1), -f32::from_bits(1)], vec![1.0; 2]),
        (vec![f32::MIN_POSITIVE], vec![0.25]),
        (vec![f32::MAX], vec![2.0]),
        (vec![-f32::MAX], vec![2.0]),
        (vec![f32::INFINITY, f32::NEG_INFINITY], vec![1.0; 2]),
        (vec![f32::INFINITY], vec![0.0]),
        (vec![f32::NAN], vec![1.0]),
    ];
    for (case, (x, y)) in cases.iter().enumerate() {
        check_inner_product(x, y, &format!("IEEE case {case}"));
    }
    let mut rng = XorShift::new(0x41cafebabe);
    // Include arbitrary signs, exponents, subnormals and NaN payloads.
    for case in 0..2048 {
        let n = [1, 2, 3, 7, 32, 127][case % 6];
        let x: Vec<_> = (0..n).map(|_| f32::from_bits(rng.next_u32())).collect();
        let y: Vec<_> = (0..n).map(|_| f32::from_bits(rng.next_u32())).collect();
        check_inner_product(&x, &y, &format!("raw bits case {case}"));
    }
}

#[test]
fn simple_conj_gradient_random_matches_c() {
    let mut rng = XorShift::new(0x87654321);
    let sizes = [1, 2, 3, 5, 8, 13, 16, 31, 64];
    let scales = [0.0, 1e-12, 1e-5, 0.1, 1.0, 10.0, 1e8];
    let b_scales = [0.0, 1e-12, 0.5, 1.0, 1e6];
    for case in 0..3780 {
        let n = sizes[case % sizes.len()];
        let m = n + [0, 1, 3, n][(case / sizes.len()) % 4];
        let niter = [-3, 0, 1, 2, 3, 5, 8, 17][(case / 3) % 8];
        let scale = scales[(case / 11) % scales.len()];
        let b_scale = b_scales[(case / 13) % b_scales.len()];
        let matrix: Vec<Vec<f32>> = (0..m)
            .map(|_| (0..n).map(|_| sample(&mut rng) * scale).collect())
            .collect();
        let b: Vec<f32> = (0..m).map(|_| sample(&mut rng) * b_scale).collect();
        // Test that previous contents, including NaNs, are discarded.
        let initial: Vec<_> = (0..n)
            .map(|_| {
                if case % 7 == 0 {
                    f32::NAN
                } else {
                    sample(&mut rng)
                }
            })
            .collect();
        check_solver(
            &matrix,
            &b,
            &initial,
            niter,
            &format!("random case {case}, n={n}, m={m}, niter={niter}"),
        );
    }
}

#[test]
fn simple_conj_gradient_structured_matches_c() {
    let mut rng = XorShift::new(0x71c0ffee);
    let sizes = [1, 2, 3, 7, 16, 31];
    for case in 0..1920 {
        let n = sizes[case % sizes.len()];
        let m = n + (case / 7) % 4;
        let mode = (case / sizes.len()) % 8;
        let niter = [-1, 0, 1, 2, 7, 20][(case / 11) % 6];
        let mut b: Vec<f32> = (0..m).map(|_| sample(&mut rng)).collect();
        let mut matrix = vec![vec![0.0; n]; m];
        let common: Vec<_> = (0..n).map(|_| sample(&mut rng)).collect();
        for j in 0..m {
            for i in 0..n {
                matrix[j][i] = match mode {
                    0 => 0.0,              // Zero matrix, independently nonzero data.
                    1 => sample(&mut rng), // Zero data below.
                    2 => {
                        if j == i {
                            1.0
                        } else {
                            0.0
                        }
                    } // Tall identity.
                    3 => common[i],        // Rank one: identical rows.
                    4 => common[0],        // Rank one: identical columns.
                    5 => {
                        if j == i {
                            [1e-4, 1.0, 2.0, 4.0][i % 4]
                        } else {
                            0.0
                        }
                    }
                    6 => {
                        if j == i && i % 2 == 0 {
                            -2.0
                        } else {
                            0.0
                        }
                    }
                    _ => {
                        if j == i {
                            -1.0
                        } else {
                            0.0
                        }
                    }
                };
            }
        }
        if mode == 1 {
            b.fill(0.0);
        }
        let result = check_solver(
            &matrix,
            &b,
            &vec![42.0; n],
            niter,
            &format!("structured case {case}, mode={mode}, niter={niter}"),
        );
        if mode == 2 && niter == 1 {
            for i in 0..n {
                compare(result[i], b[i], "identity solution after one step");
            }
        }
        if (mode == 0 || mode == 1) && niter > 0 {
            assert!(result.iter().all(|v| v.is_nan()));
        }
    }
}

#[test]
fn simple_conj_gradient_large_matches_c() {
    let mut rng = XorShift::new(0x81923745);
    for case in 0..24 {
        let n = [127, 128, 256, 257][case % 4];
        let m = n + 17;
        let niter = [1, 2, 5, 12][(case / 4) % 4];
        let matrix: Vec<Vec<f32>> = (0..m)
            .map(|_| (0..n).map(|_| sample(&mut rng) * 0.01).collect())
            .collect();
        let b: Vec<_> = (0..m).map(|_| sample(&mut rng)).collect();
        let result = check_solver(
            &matrix,
            &b,
            &vec![-99.0; n],
            niter,
            &format!("large case {case}, n={n}, m={m}, niter={niter}"),
        );
        assert!(result.iter().all(|v| v.is_finite()));
    }
}
