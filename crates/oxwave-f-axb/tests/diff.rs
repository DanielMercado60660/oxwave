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

//! Differential tests for cwp/lib/axb.c.
//! Credits: Adapted from discussions in Numerical Recipes, by Gabriel Alvarez (1995).
//! Only the C declarations, header documentation, and linked library are used.
use oxwave_cref::XorShift;
use oxwave_f_axb as rust;

extern "C" {
    fn LU_decomposition(n: i32, a: *mut *mut f32, idx: *mut i32, d: *mut f32);
    fn backward_substitution(n: i32, a: *mut *mut f32, idx: *mut i32, b: *mut f32);
    fn inverse_matrix(n: i32, a: *mut *mut f32);
    fn inverse_matrix_multiply(
        n: i32,
        a: *mut *mut f32,
        cols: i32,
        rows: i32,
        b: *mut *mut f32,
        out: *mut *mut f32,
    );
}

fn ptrs(a: &mut [Vec<f32>]) -> Vec<*mut f32> {
    a.iter_mut().map(|r| r.as_mut_ptr()).collect()
}
fn matrix(rng: &mut XorShift, n: usize, trial: usize) -> Vec<Vec<f32>> {
    let scale = [1.0, 1e-4, 1e4, 1e-12, 1e12][trial % 5];
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = if (i + j + trial) % 7 == 0 {
                0.0
            } else {
                rng.uniform_f32(-1.0, 1.0) * scale
            };
        }
        if trial % 5 != 0 || n == 1 {
            a[i][i] += (n as f32 + 1.0) * scale;
        }
    }
    if trial % 11 == 0 {
        // Signed diagonal matrices include exact off-diagonal zeros and cover
        // negative scalar systems. Every diagonal is bounded away from zero.
        for i in 0..n {
            a[i].fill(0.0);
            let sign = if (i + trial) % 2 == 0 { -1.0 } else { 1.0 };
            a[i][i] = sign * rng.uniform_f32(0.5, 2.0) * scale;
        }
    }
    // Force permutations, while keeping the matrix well conditioned.
    if n > 1 && trial % 2 == 0 {
        a.swap(0, n - 1);
    }
    a
}
fn compare(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    assert!(a.iter().chain(b).all(|x| x.is_finite()));
    let max = b.iter().fold(0.0f32, |m, x| m.max(x.abs()));
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() <= max * 1e-5, "{x} != {y}, scale {max}");
    }
}
fn compare_matrix(a: &[Vec<f32>], b: &[Vec<f32>]) {
    compare(
        &a.iter().flatten().copied().collect::<Vec<_>>(),
        &b.iter().flatten().copied().collect::<Vec<_>>(),
    );
}

// Nonfinite outputs are part of the reference behavior for singular systems.
// Compare their classification separately so NaNs/infinities cannot mask a
// mismatch or loosen the tolerance for finite entries.
fn compare_nonfinite(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    let max = b
        .iter()
        .filter(|x| x.is_finite())
        .fold(0.0f64, |m, x| m.max(f64::from(x.abs())));
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        if y.is_nan() {
            assert!(x.is_nan(), "entry {i}: {x} != NaN");
        } else if y.is_infinite() {
            assert_eq!(x.to_bits(), y.to_bits(), "entry {i}: infinity sign");
        } else {
            assert!(x.is_finite(), "entry {i}: {x} != {y}");
            assert!(
                (f64::from(x) - f64::from(y)).abs() <= max * 1e-5,
                "entry {i}: {x} != {y}, scale {max}"
            );
        }
    }
}

fn compare_nonfinite_matrix(a: &[Vec<f32>], b: &[Vec<f32>]) {
    assert_eq!(a.len(), b.len());
    for (row, reference) in a.iter().zip(b) {
        assert_eq!(row.len(), reference.len());
    }
    compare_nonfinite(
        &a.iter().flatten().copied().collect::<Vec<_>>(),
        &b.iter().flatten().copied().collect::<Vec<_>>(),
    );
}

fn check_lu(mut a: Vec<Vec<f32>>) -> Vec<i32> {
    let n = a.len();
    let mut c = a.clone();
    let mut index = vec![-1; n];
    let mut c_index = index.clone();
    let mut d = 0.0;
    let mut c_d = 0.0;
    unsafe {
        LU_decomposition(
            n as i32,
            ptrs(&mut c).as_mut_ptr(),
            c_index.as_mut_ptr(),
            &mut c_d,
        );
    }
    rust::LU_decomposition(&mut a, &mut index, &mut d);
    assert_eq!(index, c_index);
    assert_eq!(d.to_bits(), c_d.to_bits());
    compare_nonfinite_matrix(&a, &c);
    index
}
const SIZES: [usize; 12] = [1, 2, 3, 5, 7, 8, 16, 31, 32, 64, 127, 128];

fn trials(n: usize) -> usize {
    if n > 64 {
        20
    } else {
        300
    }
}

#[test]
fn lu_zero_pivot_matches_c() {
    // Each row is nonzero, but the rows are identical: the last pivot is zero.
    let mut a = vec![vec![1.0, -2.0], vec![1.0, -2.0]];
    let mut c = a.clone();
    let mut index = vec![0; 2];
    let mut c_index = index.clone();
    let mut d = 0.0;
    let mut c_d = 0.0;
    unsafe {
        LU_decomposition(2, ptrs(&mut c).as_mut_ptr(), c_index.as_mut_ptr(), &mut c_d);
    }
    rust::LU_decomposition(&mut a, &mut index, &mut d);
    assert_eq!(index, c_index);
    assert_eq!(d, c_d);
    compare_matrix(&a, &c);
    assert_eq!(a[1][1], c[1][1]);
}

#[test]
fn lu_zero_rows_match_c() {
    for a in [
        vec![vec![0.0]],
        vec![vec![-0.0]],
        vec![vec![0.0; 2]; 2],
        vec![vec![0.0, 0.0], vec![1.0, -2.0]],
        vec![vec![1.0, -2.0], vec![0.0, 0.0]],
        vec![vec![2.0, 1.0, 0.0], vec![0.0; 3], vec![-1.0, 3.0, 4.0]],
        // A zero row and NaN scores occur together in sufxdecon inputs.
        vec![vec![2.0, f32::NAN], vec![0.0, 0.0]],
    ] {
        check_lu(a);
    }
}

#[test]
fn lu_nan_scores_keep_previous_pivot_matches_c() {
    // Column zero chooses row zero. Every score in column one is NaN,
    // so the remembered row must remain zero and be exchanged with row one.
    assert_eq!(
        check_lu(vec![vec![2.0, f32::NAN], vec![1.0, 2.0]]),
        vec![0, 0]
    );
    // The same remembered row can survive several columns.
    assert_eq!(
        check_lu(vec![
            vec![2.0, f32::NAN, 1.0],
            vec![1.0, 3.0, 1.0],
            vec![1.0, 1.0, 4.0],
        ]),
        vec![0, 0, 0]
    );
}

#[test]
fn lu_singular_diagnostic_matches_c() {
    const CHILD: &str = "OXWAVE_AXB_DIAGNOSTIC_CHILD";
    if let Some(backend) = std::env::var_os(CHILD) {
        let mut a = vec![vec![0.0; 3]; 3];
        let mut index = vec![-1; 3];
        let mut d = 0.0;
        match backend.to_str().unwrap() {
            "c" => unsafe {
                LU_decomposition(3, ptrs(&mut a).as_mut_ptr(), index.as_mut_ptr(), &mut d);
            },
            "rust" => rust::LU_decomposition(&mut a, &mut index, &mut d),
            _ => panic!("unknown diagnostic backend"),
        }
        // A successful test result also proves the call returned normally.
        assert!(index.iter().all(|&pivot| (0..3).contains(&pivot)));
        return;
    }

    // Separate processes capture actual stderr without redirecting global
    // descriptors or interfering with concurrently running differential tests.
    let run = |backend| {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "lu_singular_diagnostic_matches_c",
                "--nocapture",
                "--test-threads=1",
            ])
            .env(CHILD, backend)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{backend}: {}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        output.stderr
    };
    let c_stderr = run("c");
    assert_eq!(
        c_stderr,
        b"error, singular matrix in LU decomposition\n".repeat(3)
    );
    assert_eq!(run("rust"), c_stderr);
}

#[test]
fn singular_and_nonfinite_systems_match_c() {
    const CHILD: &str = "OXWAVE_AXB_NONFINITE_CHILD";
    if std::env::var_os(CHILD).is_none() {
        // Singular rows deliberately emit many diagnostics. Capture them in a
        // child so ordinary test output remains readable; the separate stderr
        // regression checks the exact diagnostic and successful continuation.
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "singular_and_nonfinite_systems_match_c",
                "--nocapture",
                "--test-threads=1",
            ])
            .env(CHILD, "1")
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert!(
            output.status.success() && stdout.contains("1 passed; 0 failed"),
            "{stdout}\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }

    let mut rng = XorShift::new(7);
    for &n in &SIZES {
        for trial in 0..trials(n) {
            let mut a = matrix(&mut rng, n, trial);
            let row = rng.next_u32() as usize % n;
            let col = rng.next_u32() as usize % n;
            match trial % 6 {
                0 => a[row].fill(0.0),
                1 => a.iter_mut().for_each(|r| r.fill(-0.0)),
                2 => {
                    for r in &mut a {
                        r[col] = f32::NAN;
                    }
                }
                3 => {
                    // Zero rows and NaN scores must coexist, as in
                    // the caller that exposed the two LU defects.
                    a[row].fill(0.0);
                    a[(row + 1) % n][col] = f32::NAN;
                }
                4 => {
                    a[row][col] = f32::INFINITY;
                    a[(row + 1) % n][(col + 1) % n] = f32::NEG_INFINITY;
                }
                _ => {
                    a[row].fill(-0.0);
                    a[(row + 1) % n].fill(f32::NAN);
                }
            }
            check_lu(a.clone());

            let mut factors = a.clone();
            let mut index = vec![-1; n];
            let mut d = 0.0;
            unsafe {
                LU_decomposition(
                    n as i32,
                    ptrs(&mut factors).as_mut_ptr(),
                    index.as_mut_ptr(),
                    &mut d,
                );
            }
            let mut b: Vec<_> = (0..n)
                .map(|i| {
                    if (i + trial) % 3 == 0 {
                        0.0
                    } else {
                        rng.uniform_f32(-10.0, 10.0)
                    }
                })
                .collect();
            let mut c_b = b.clone();
            rust::backward_substitution(&factors, &index, &mut b);
            unsafe {
                backward_substitution(
                    n as i32,
                    ptrs(&mut factors).as_mut_ptr(),
                    index.as_mut_ptr(),
                    c_b.as_mut_ptr(),
                );
            }
            compare_nonfinite(&b, &c_b);

            let mut inverse = a.clone();
            let mut c_inverse = a.clone();
            rust::inverse_matrix(&mut inverse);
            unsafe {
                inverse_matrix(n as i32, ptrs(&mut c_inverse).as_mut_ptr());
            }
            compare_nonfinite_matrix(&inverse, &c_inverse);

            let cols = 1 + trial % 7;
            let rhs: Vec<Vec<_>> = (0..n)
                .map(|i| {
                    (0..cols)
                        .map(|j| {
                            if (i + j + trial) % 3 == 0 {
                                0.0
                            } else {
                                rng.uniform_f32(-10.0, 10.0)
                            }
                        })
                        .collect()
                })
                .collect();
            let mut out = vec![vec![f32::NAN; cols]; n];
            let mut c_out = out.clone();
            let mut c_a = a.clone();
            let mut c_rhs = rhs.clone();
            rust::inverse_matrix_multiply(&a, &rhs, &mut out);
            unsafe {
                inverse_matrix_multiply(
                    n as i32,
                    ptrs(&mut c_a).as_mut_ptr(),
                    cols as i32,
                    n as i32,
                    ptrs(&mut c_rhs).as_mut_ptr(),
                    ptrs(&mut c_out).as_mut_ptr(),
                );
            }
            compare_nonfinite_matrix(&out, &c_out);
            assert_eq!(c_rhs, rhs);
            // NaNs do not compare equal, but preserved input bits must match.
            for (row, reference) in c_a.iter().zip(&a) {
                for (x, y) in row.iter().zip(reference) {
                    assert_eq!(x.to_bits(), y.to_bits());
                }
            }
        }
    }
}

#[test]
fn lu_pivot_ties_and_row_padding_match_c() {
    // Equal scaled pivot scores exercise the tie rule. The extra column is
    // caller storage outside the active square matrix and must stay in its row.
    for mut a in [
        vec![vec![1.0, 1.0, 100.0], vec![-1.0, 1.0, 200.0]],
        vec![vec![0.0, -2.0, 300.0], vec![3.0, 0.0, 400.0]],
    ] {
        let mut c = a.clone();
        let mut index = vec![0; 2];
        let mut c_index = index.clone();
        let mut d = 0.0;
        let mut c_d = 0.0;
        unsafe {
            LU_decomposition(2, ptrs(&mut c).as_mut_ptr(), c_index.as_mut_ptr(), &mut c_d);
        }
        rust::LU_decomposition(&mut a, &mut index, &mut d);
        assert_eq!(index, c_index);
        assert_eq!(d, c_d);
        for i in 0..2 {
            assert_eq!(a[i][2].to_bits(), c[i][2].to_bits());
        }
        compare_matrix(&a, &c);
    }
}

#[test]
fn lu_matches_c() {
    let mut rng = XorShift::new(1);
    for &n in &SIZES {
        for trial in 0..trials(n) {
            let mut a = matrix(&mut rng, n, trial);
            let mut c = a.clone();
            let mut idx = vec![0; n];
            let mut ci = idx.clone();
            let mut d = 0.0;
            let mut cd = 0.0;
            unsafe {
                LU_decomposition(
                    n as i32,
                    ptrs(&mut c).as_mut_ptr(),
                    ci.as_mut_ptr(),
                    &mut cd,
                );
            }
            rust::LU_decomposition(&mut a, &mut idx, &mut d);
            assert_eq!(idx, ci);
            assert_eq!(d, cd);
            compare_matrix(&a, &c);
        }
    }
}

#[test]
fn substitution_matches_c() {
    let mut rng = XorShift::new(2);
    for &n in &SIZES {
        for trial in 0..trials(n) {
            let mut a = matrix(&mut rng, n, trial);
            let mut idx = vec![0; n];
            let mut d = 0.0;
            // Feed the very same independently C-factored matrix to both
            // solvers so decomposition errors cannot mask solve errors.
            unsafe {
                LU_decomposition(
                    n as i32,
                    ptrs(&mut a).as_mut_ptr(),
                    idx.as_mut_ptr(),
                    &mut d,
                );
            }
            let mut b: Vec<_> = (0..n)
                .map(|i| {
                    if trial % 13 == 0 || i < trial % n || (i + trial) % 4 == 0 {
                        0.0
                    } else {
                        rng.uniform_f32(-10.0, 10.0)
                    }
                })
                .collect();
            let mut cb = b.clone();
            let mut ca = a.clone();
            let mut ci = idx.clone();
            unsafe {
                backward_substitution(
                    n as i32,
                    ptrs(&mut ca).as_mut_ptr(),
                    ci.as_mut_ptr(),
                    cb.as_mut_ptr(),
                );
            }
            rust::backward_substitution(&a, &idx, &mut b);
            compare(&b, &cb);
            assert_eq!(a, ca);
            assert_eq!(idx, ci);
        }
    }
}

#[test]
fn inverse_matches_c() {
    let mut rng = XorShift::new(3);
    for &n in &SIZES {
        for trial in 0..trials(n) {
            let mut a = matrix(&mut rng, n, trial);
            let mut c = a.clone();
            unsafe {
                inverse_matrix(n as i32, ptrs(&mut c).as_mut_ptr());
            }
            rust::inverse_matrix(&mut a);
            compare_matrix(&a, &c);
        }
    }
}

#[test]
fn inverse_multiply_matches_c() {
    let mut rng = XorShift::new(4);
    for &n in &SIZES {
        for trial in 0..trials(n) {
            let a = matrix(&mut rng, n, trial);
            let cols = 1 + trial % 7;
            let b: Vec<Vec<f32>> = (0..n)
                .map(|i| {
                    (0..cols)
                        .map(|j| {
                            if trial % 13 == 0 || (i + j + trial) % 5 == 0 {
                                0.0
                            } else {
                                rng.uniform_f32(-10.0, 10.0)
                            }
                        })
                        .collect()
                })
                .collect();
            let mut ca = a.clone();
            let mut cb = b.clone();
            let mut out = vec![vec![0.0; cols]; n];
            let mut co = out.clone();
            unsafe {
                inverse_matrix_multiply(
                    n as i32,
                    ptrs(&mut ca).as_mut_ptr(),
                    cols as i32,
                    n as i32,
                    ptrs(&mut cb).as_mut_ptr(),
                    ptrs(&mut co).as_mut_ptr(),
                );
            }
            rust::inverse_matrix_multiply(&a, &b, &mut out);
            compare_matrix(&out, &co);
            assert_eq!(a, ca);
            assert_eq!(b, cb);
        }
    }
}

#[test]
fn padded_storage_matches_c() {
    let mut rng = XorShift::new(5);
    for n in [1, 2, 3, 5, 8, 31, 32] {
        for trial in 0..32 {
            let original = matrix(&mut rng, n, trial);
            // Distinct suffixes detect row swaps or writes beyond the active
            // square. Keep them out of the floating-point comparison scale.
            let mut a = original.clone();
            for (i, row) in a.iter_mut().enumerate() {
                row.extend_from_slice(&[-0.0, 1000.0 + i as f32]);
            }
            let padded = a.clone();
            let mut c = a.clone();
            let mut index = vec![-99; n + 2];
            let mut c_index = index.clone();
            let mut d = 0.0;
            let mut c_d = 0.0;
            unsafe {
                LU_decomposition(
                    n as i32,
                    ptrs(&mut c).as_mut_ptr(),
                    c_index.as_mut_ptr(),
                    &mut c_d,
                );
            }
            rust::LU_decomposition(&mut a, &mut index, &mut d);
            assert_eq!(index, c_index);
            assert_eq!(&index[n..], &[-99; 2]);
            assert_eq!(d.to_bits(), c_d.to_bits());
            for i in 0..n {
                compare(&a[i][..n], &c[i][..n]);
                for j in n..n + 2 {
                    assert_eq!(a[i][j].to_bits(), padded[i][j].to_bits());
                    assert_eq!(c[i][j].to_bits(), padded[i][j].to_bits());
                }
            }

            let mut b: Vec<_> = (0..n).map(|_| rng.uniform_f32(-1.0, 1.0)).collect();
            b.extend_from_slice(&[-0.0, 1234.0]);
            let before_b = b.clone();
            let mut c_b = b.clone();
            // Use identical C factors to isolate substitution from LU error.
            let factors = c.clone();
            let pivots = c_index.clone();
            rust::backward_substitution(&factors, &pivots, &mut b);
            unsafe {
                backward_substitution(
                    n as i32,
                    ptrs(&mut c).as_mut_ptr(),
                    c_index.as_mut_ptr(),
                    c_b.as_mut_ptr(),
                );
            }
            compare(&b[..n], &c_b[..n]);
            assert_eq!(c_index, pivots);
            for (row, expected) in c.iter().zip(&factors) {
                for (x, y) in row.iter().zip(expected) {
                    assert_eq!(x.to_bits(), y.to_bits());
                }
            }
            for j in n..n + 2 {
                assert_eq!(b[j].to_bits(), before_b[j].to_bits());
                assert_eq!(c_b[j].to_bits(), before_b[j].to_bits());
            }

            let mut a = padded.clone();
            let mut c = padded.clone();
            rust::inverse_matrix(&mut a);
            unsafe {
                inverse_matrix(n as i32, ptrs(&mut c).as_mut_ptr());
            }
            for i in 0..n {
                compare(&a[i][..n], &c[i][..n]);
                for j in n..n + 2 {
                    assert_eq!(a[i][j].to_bits(), padded[i][j].to_bits());
                    assert_eq!(c[i][j].to_bits(), padded[i][j].to_bits());
                }
            }
        }
    }
}

#[test]
fn inverse_multiply_wide_and_empty_rhs_match_c() {
    let mut rng = XorShift::new(6);
    for n in [1, 2, 3, 5, 8, 31, 32] {
        for cols in [0, 1, 2, 7, 32, 67, 128] {
            for trial in 0..8 {
                let a = matrix(&mut rng, n, trial);
                let b: Vec<Vec<f32>> = (0..n)
                    .map(|_| (0..cols).map(|_| rng.uniform_f32(-10.0, 10.0)).collect())
                    .collect();
                let mut c_a = a.clone();
                let mut c_b = b.clone();
                // NaNs in active output storage detect unwritten entries;
                // signed-zero suffixes must be preserved bit for bit.
                let mut out = vec![vec![f32::NAN; cols + 2]; n];
                for row in &mut out {
                    row[cols] = -0.0;
                    row[cols + 1] = 4321.0;
                }
                let mut c_out = out.clone();
                rust::inverse_matrix_multiply(&a, &b, &mut out);
                unsafe {
                    inverse_matrix_multiply(
                        n as i32,
                        ptrs(&mut c_a).as_mut_ptr(),
                        cols as i32,
                        n as i32,
                        ptrs(&mut c_b).as_mut_ptr(),
                        ptrs(&mut c_out).as_mut_ptr(),
                    );
                }
                assert_eq!(c_a, a);
                assert_eq!(c_b, b);
                for i in 0..n {
                    compare(&out[i][..cols], &c_out[i][..cols]);
                    for j in cols..cols + 2 {
                        let expected: f32 = if j == cols { -0.0 } else { 4321.0 };
                        assert_eq!(out[i][j].to_bits(), expected.to_bits());
                        assert_eq!(c_out[i][j].to_bits(), expected.to_bits());
                    }
                }
            }
        }
    }
}
