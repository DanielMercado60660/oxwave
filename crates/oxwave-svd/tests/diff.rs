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

//! Differential tests for cwp/lib/svd.c. Original credits: Luke Tierney,
//! David Betz, Ian Kay (1999), and Nils Maercklin (2001); published algorithm:
//! Golub and Reinsch (1971). No restricted C function bodies were viewed or used.

use oxwave_cref::XorShift;
use oxwave_svd as rust;
use std::sync::Mutex;

static C_LOCK: Mutex<()> = Mutex::new(());

// Keep the original C names in the declarations; Rust calls use the `rust` namespace.
extern "C" {
    fn compute_svd(a: *mut *mut f32, m: i32, n: i32, w: *mut f32, v: *mut *mut f32) -> i32;
    fn svd_backsubstitute(
        u: *mut *mut f32,
        w: *mut f32,
        v: *mut *mut f32,
        m: i32,
        n: i32,
        b: *mut f32,
        x: *mut f32,
    );
    fn svd_sort(u: *mut *mut f32, w: *mut f32, v: *mut *mut f32, n: i32, m: i32);
}

fn pointers(a: &mut [Vec<f32>]) -> Vec<*mut f32> {
    a.iter_mut().map(|row| row.as_mut_ptr()).collect()
}

fn close(a: &[f64], b: &[f64], label: &str) {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().fold(0.0f64, |s, &x| s.max(x.abs()));
    for (i, (&a, &b)) in a.iter().zip(b).enumerate() {
        assert!(
            (a - b).abs() <= 1e-5 * scale.max(1e-30),
            "{label}[{i}]: Rust={a}, C/reference={b}, scale={scale}"
        );
    }
}

fn reconstruct(u: &[Vec<f32>], w: &[f32], v: &[Vec<f32>]) -> Vec<f64> {
    let mut result = vec![0.0; u.len() * w.len()];
    for i in 0..u.len() {
        for j in 0..w.len() {
            for k in 0..w.len() {
                result[i * w.len() + j] += u[i][k] as f64 * w[k] as f64 * v[j][k] as f64;
            }
        }
    }
    result
}

fn gram(a: &[Vec<f32>], n: usize) -> Vec<f64> {
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            for row in a {
                result[i * n + j] += row[i] as f64 * row[j] as f64;
            }
        }
    }
    result
}

// Compare every output directly against C, as well as the defining equations
// on matrices where SU produces a mathematical SVD. Retain the independent
// reconstruction, orthogonality, and separated-spectrum alignment checks.
fn check_svd(input: Vec<Vec<f32>>, compare_vectors: bool, label: &str) {
    let m = input.len();
    let n = input.first().map_or(0, Vec::len);
    let mut u = input.clone();
    let mut cu = input.clone();
    let mut w = vec![-999.0; n];
    let mut cw = w.clone();
    let mut v = vec![vec![-999.0; n]; n];
    let mut cv = v.clone();
    let status = rust::compute_svd(&mut u, &mut w, &mut v);
    let c_status = unsafe {
        compute_svd(
            pointers(&mut cu).as_mut_ptr(),
            m as i32,
            n as i32,
            cw.as_mut_ptr(),
            pointers(&mut cv).as_mut_ptr(),
        )
    };
    assert_eq!(status, c_status, "{label}: return status");
    assert_eq!(status, 1, "{label}: convergence");
    assert!(
        w.iter().all(|x| x.is_finite() && *x >= 0.0),
        "{label}: singular values"
    );
    close_outputs(&u, &w, &v, &cu, &cw, &cv, label);

    let mut order: Vec<usize> = (0..n).collect();
    let mut corder = order.clone();
    order.sort_by(|&i, &j| w[j].total_cmp(&w[i]));
    corder.sort_by(|&i, &j| cw[j].total_cmp(&cw[i]));
    let ordered: Vec<f64> = order.iter().map(|&i| w[i] as f64).collect();
    let cordered: Vec<f64> = corder.iter().map(|&i| cw[i] as f64).collect();
    close(&ordered, &cordered, &format!("{label}: singular values"));

    let original: Vec<f64> = input.iter().flatten().map(|&x| x as f64).collect();
    let reconstructed = reconstruct(&u, &w, &v);
    let c_reconstructed = reconstruct(&cu, &cw, &cv);
    close(
        &reconstructed,
        &c_reconstructed,
        &format!("{label}: Rust/C reconstruction"),
    );
    close(
        &reconstructed,
        &original,
        &format!("{label}: reconstruction"),
    );
    close(
        &c_reconstructed,
        &original,
        &format!("{label}: C reconstruction"),
    );
    let mut identity = vec![0.0; n * n];
    for i in 0..n {
        identity[i * n + i] = 1.0;
    }
    let ug = gram(&u, n);
    let vg = gram(&v, n);
    close(&ug, &identity, &format!("{label}: U orthogonality"));
    close(&vg, &identity, &format!("{label}: V orthogonality"));
    close(
        &ug,
        &gram(&cu, n),
        &format!("{label}: Rust/C U Gram matrix"),
    );
    close(
        &vg,
        &gram(&cv, n),
        &format!("{label}: Rust/C V Gram matrix"),
    );

    if compare_vectors {
        let mut aligned_u = Vec::with_capacity(m * n);
        let mut ordered_cu = Vec::with_capacity(m * n);
        let mut aligned_v = Vec::with_capacity(n * n);
        let mut ordered_cv = Vec::with_capacity(n * n);
        for j in 0..n {
            let col = order[j];
            let ccol = corder[j];
            let mut dot = 0.0f64;
            for i in 0..m {
                dot += u[i][col] as f64 * cu[i][ccol] as f64;
            }
            let sign = if dot < 0.0 { -1.0 } else { 1.0 };
            for i in 0..m {
                aligned_u.push(sign * u[i][col] as f64);
                ordered_cu.push(cu[i][ccol] as f64);
            }
            for i in 0..n {
                aligned_v.push(sign * v[i][col] as f64);
                ordered_cv.push(cv[i][ccol] as f64);
            }
        }
        close(
            &aligned_u,
            &ordered_cu,
            &format!("{label}: individual U vectors"),
        );
        close(
            &aligned_v,
            &ordered_cv,
            &format!("{label}: individual V vectors"),
        );
    }
}

#[test]
fn compute_svd_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0x5ad);
    for case in 0..3000 {
        let n = [1, 2, 3, 5, 7, 11, 16, 17, 31, 32, 64][case % 11];
        let m = n + [0, 1, 2, 7, 32, 97][case % 6];
        let scale = [1e-12, 1e-6, 1.0, 1e6, 1e12][case % 5];
        let input = (0..m).map(|_| rng.vec_f32(n, -scale, scale)).collect();
        check_svd(input, false, &format!("dense {case}, {m}x{n}"));
    }
    for &(m, n) in &[(257, 3), (257, 31), (1024, 7), (128, 128)] {
        let input = (0..m).map(|_| rng.vec_f32(n, -1.0, 1.0)).collect();
        check_svd(input, false, &format!("large {m}x{n}"));
    }
}

#[test]
fn compute_svd_single_columns_match_c_bits() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0xc010_5ca1);
    for case in 0..4000 {
        let m = [1, 2, 3, 5, 7, 16, 31, 64, 127, 257, 1024][case % 11];
        let scale = [1e-30, 1e-15, 1e-6, 1.0, 1e6, 1e15, 1e30][case % 7];
        let mut input = rng.vec_f32(m, -scale, scale);
        match case % 17 {
            0 => input.fill(0.0),
            1 => input.fill(-0.0),
            2 => input.fill(-scale),
            3 => {
                input.fill(0.0);
                input[m - 1] = -scale;
            }
            _ => (),
        }
        let mut u: Vec<Vec<f32>> = input.into_iter().map(|x| vec![x]).collect();
        let mut cu = u.clone();
        let mut w = [0.0];
        let mut cw = w;
        let mut v = vec![vec![0.0]];
        let mut cv = v.clone();
        assert_eq!(rust::compute_svd(&mut u, &mut w, &mut v), 1);
        assert_eq!(
            unsafe {
                compute_svd(
                    pointers(&mut cu).as_mut_ptr(),
                    m as i32,
                    1,
                    cw.as_mut_ptr(),
                    pointers(&mut cv).as_mut_ptr(),
                )
            },
            1
        );
        assert_eq!(w[0].to_bits(), cw[0].to_bits(), "w case {case}");
        assert_eq!(v[0][0].to_bits(), cv[0][0].to_bits(), "v case {case}");
        for i in 0..m {
            assert_eq!(u[i][0].to_bits(), cu[i][0].to_bits(), "u[{i}] case {case}");
        }
    }
}

#[test]
fn compute_svd_individual_vectors_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0x5ed_a7a);
    for case in 0..3000 {
        let n = [1, 2, 3, 5, 7][case % 5];
        let m = n + [0, 1, 7, 32][case % 4];
        let scale = [1e-12, 1.0, 1e12][case % 3];
        let mut input = vec![vec![0.0; n]; m];
        for j in 0..n {
            input[j][j] = (10.0 * (n - j) as f32 + rng.uniform_f32(-0.5, 0.5)) * scale;
        }
        // Random orthogonal row and column rotations make the factors dense,
        // preserving a deliberately separated, nonzero singular spectrum.
        for _ in 0..8 * (m + n) {
            let row1 = rng.range_i32(0, m as i32 - 1) as usize;
            let row2 = rng.range_i32(0, m as i32 - 1) as usize;
            let (s, c) = rng.uniform(-3.0, 3.0).sin_cos();
            if row1 != row2 {
                for j in 0..n {
                    let x = input[row1][j] as f64;
                    let y = input[row2][j] as f64;
                    input[row1][j] = (c * x + s * y) as f32;
                    input[row2][j] = (c * y - s * x) as f32;
                }
            }
            let col1 = rng.range_i32(0, n as i32 - 1) as usize;
            let col2 = rng.range_i32(0, n as i32 - 1) as usize;
            if col1 != col2 {
                for row in &mut input {
                    let x = row[col1] as f64;
                    let y = row[col2] as f64;
                    row[col1] = (c * x + s * y) as f32;
                    row[col2] = (c * y - s * x) as f32;
                }
            }
        }
        check_svd(input, true, &format!("separated spectrum {case}, {m}x{n}"));
    }
}

#[test]
fn compute_svd_degenerate_and_sparse_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    check_svd(Vec::new(), true, "empty");
    for n in [1, 2, 3, 5, 7, 16, 31, 32, 64] {
        for m in [n, n + 5] {
            for kind in 0..8 {
                let mut a = vec![vec![0.0; n]; m];
                for i in 0..m {
                    for j in 0..n {
                        a[i][j] = match kind {
                            0 => 0.0,
                            1 => {
                                if i == j {
                                    1.0
                                } else {
                                    0.0
                                }
                            }
                            2 => {
                                if i == j {
                                    -(j as f32 + 1.0)
                                } else {
                                    0.0
                                }
                            }
                            3 => 1.0, // rank one and repeated zero singular values
                            4 => {
                                if j < 2 {
                                    i as f32 - 2.0
                                } else {
                                    0.0
                                }
                            }
                            5 => {
                                if i == j {
                                    if j == 0 {
                                        1e-6
                                    } else {
                                        1.0
                                    }
                                } else {
                                    0.0
                                }
                            }
                            6 => {
                                if (i + j) % 3 == 0 {
                                    ((i * 7 + j * 3) % 13) as f32 - 6.0
                                } else {
                                    0.0
                                }
                            }
                            _ => {
                                if i == j {
                                    -0.0
                                } else {
                                    0.0
                                }
                            }
                        };
                    }
                }
                check_svd(a, kind == 2, &format!("edge {kind}, {m}x{n}"));
            }
        }
    }
}

// SU's cancellation step can produce a factorization different from the
// mathematical SVD. Compare its actual outputs, without reconstruction against
// the input or any sign/permutation alignment that could hide a difference.
fn check_svd_outputs(input: Vec<Vec<f32>>, label: &str) -> Vec<f32> {
    let m = input.len();
    let n = input[0].len();
    let mut u = input.clone();
    let mut cu = input;
    let mut w = vec![0.0; n];
    let mut cw = w.clone();
    let mut v = vec![vec![0.0; n]; n];
    let mut cv = v.clone();
    let status = rust::compute_svd(&mut u, &mut w, &mut v);
    let c_status = unsafe {
        compute_svd(
            pointers(&mut cu).as_mut_ptr(),
            m as i32,
            n as i32,
            cw.as_mut_ptr(),
            pointers(&mut cv).as_mut_ptr(),
        )
    };
    assert_eq!(status, c_status, "{label}: return status");
    assert_eq!(status, 1, "{label}: convergence");
    close_outputs(&u, &w, &v, &cu, &cw, &cv, label);
    w
}

fn close_outputs(
    u: &[Vec<f32>],
    w: &[f32],
    v: &[Vec<f32>],
    cu: &[Vec<f32>],
    cw: &[f32],
    cv: &[Vec<f32>],
    label: &str,
) {
    close(
        &w.iter().copied().map(f64::from).collect::<Vec<_>>(),
        &cw.iter().copied().map(f64::from).collect::<Vec<_>>(),
        &format!("{label}: unpermuted singular values"),
    );
    for (a, b, factor) in [(u, cu, "U"), (v, cv, "V")] {
        close(
            &a.iter()
                .flatten()
                .copied()
                .map(f64::from)
                .collect::<Vec<_>>(),
            &b.iter()
                .flatten()
                .copied()
                .map(f64::from)
                .collect::<Vec<_>>(),
            &format!("{label}: unaligned {factor}"),
        );
    }
}

#[test]
fn compute_svd_cancellation_regression_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut w = check_svd_outputs(
        vec![vec![0.0, 1.0, 0.0], vec![0.0, 1.0, 1.0], vec![0.0; 3]],
        "review cancellation regression",
    );
    w.sort_by(|a, b| b.total_cmp(a));
    close(
        &w.into_iter().map(f64::from).collect::<Vec<_>>(),
        &[1.78498, 0.560231, 0.0],
        "review singular values",
    );
}

#[test]
fn compute_svd_negligible_cancellation_continues_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    // After cancelling the first off-diagonal, s is about coupling. The next
    // f = s * gap can be negligible even though the following f is significant.
    // SU must continue scanning and rotate the final pair. An exact-zero gap
    // additionally exercises splitting before cancellation instead.
    for scale in [1e-12, 1e-6, 1.0, 1e6, 1e12] {
        for coupling in [-1e-2, -1e-3, 1e-3, 1e-2] {
            for gap in [0.0, 1e-20, 1e-16, 1e-14, 1e-12] {
                let mut a = vec![
                    vec![0.0, coupling, 0.0, 0.0],
                    vec![0.0, 1.0, gap, 0.0],
                    vec![0.0, 0.0, 2.0, 1.0],
                    vec![0.0, 0.0, 0.0, 3.0],
                ];
                for x in a.iter_mut().flatten() {
                    *x *= scale;
                }
                check_svd_outputs(
                    a,
                    &format!("cancellation scale={scale}, coupling={coupling}, gap={gap}"),
                );
            }
        }
    }
}

#[test]
fn compute_svd_rank_deficient_outputs_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0xca11_ce11);
    for case in 0..3000 {
        let n = [2, 3, 5, 7, 16, 31, 32, 64][case % 8];
        let m = n + [0, 1, 7, 31][case % 4];
        let scale = [1e-12, 1e-6, 1.0, 1e6, 1e12][case % 5];
        let mut a: Vec<Vec<f32>> = (0..m).map(|_| rng.vec_f32(n, -scale, scale)).collect();
        match case % 3 {
            0 => {
                // A zero leading column with nonzero entries to its right
                // forces cancellation of a nonzero bidiagonal off-diagonal.
                for row in &mut a {
                    row[0] = 0.0;
                }
            }
            1 => {
                let column = rng.range_i32(0, n as i32 - 1) as usize;
                for row in &mut a {
                    row[column] = 0.0;
                }
            }
            _ => {
                for row in &mut a {
                    row[n - 1] = row[0];
                }
            }
        }
        let label = format!("rank deficient {case}, {m}x{n}");
        check_svd_outputs(a, &label);
    }
    for (m, n) in [(257, 31), (257, 127), (1024, 7), (128, 128)] {
        let mut a: Vec<Vec<f32>> = (0..m).map(|_| rng.vec_f32(n, -1.0, 1.0)).collect();
        for row in &mut a {
            row[0] = 0.0;
        }
        check_svd_outputs(a, &format!("large rank deficient {m}x{n}"));
    }
}

#[test]
fn compute_svd_rejects_wide_matrices_like_c() {
    // Re-execute this test to capture C and Rust diagnostics independently,
    // without redirecting process-global file descriptors in a test thread.
    let mode = std::env::var("OXWAVE_SVD_WIDE_DIAGNOSTIC");
    if mode.is_err() {
        let run = |implementation| {
            std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "compute_svd_rejects_wide_matrices_like_c",
                    "--nocapture",
                ])
                .env("OXWAVE_SVD_WIDE_DIAGNOSTIC", implementation)
                .output()
                .unwrap()
        };
        let c = run("C");
        let rust = run("Rust");
        assert!(c.status.success() && rust.status.success());
        assert_eq!(rust.stderr, c.stderr);
        assert_eq!(c.stderr, b"#rows must be > #cols \n".repeat(3));
        assert!(!String::from_utf8_lossy(&rust.stdout).contains("#rows"));
        assert!(!String::from_utf8_lossy(&c.stdout).contains("#rows"));
        return;
    }
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    for (m, n) in [(0, 1), (1, 2), (2, 3)] {
        let mut u = vec![vec![7.0; n]; m];
        let mut w = vec![-999.0; n];
        let mut v = vec![vec![-999.0; n]; n];
        let original_u = u.clone();
        let original_w = w.clone();
        let original_v = v.clone();
        let status = if mode.as_deref() == Ok("C") {
            unsafe {
                compute_svd(
                    pointers(&mut u).as_mut_ptr(),
                    m as i32,
                    n as i32,
                    w.as_mut_ptr(),
                    pointers(&mut v).as_mut_ptr(),
                )
            }
        } else {
            rust::compute_svd(&mut u, &mut w, &mut v)
        };
        assert_eq!(status, 0);
        assert_eq!(u, original_u);
        assert_eq!(w, original_w);
        assert_eq!(v, original_v);
    }
}

#[test]
fn svd_backsubstitute_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0xbacc);
    for case in 0..6000 {
        let m = [0, 1, 2, 3, 7, 16, 31, 64, 127, 257][case % 10];
        let n = [0, 1, 2, 3, 5, 11, 16, 31, 32, 64, 128][case % 11];
        let bscale = [1e-15, 1e-6, 1.0, 1e6, 1e15][case % 5];
        let wscale = [1e-12, 1.0, 1e12][case % 3];
        let mut u: Vec<Vec<f32>> = (0..m).map(|_| rng.vec_f32(n, -1.0, 1.0)).collect();
        let mut v: Vec<Vec<f32>> = (0..n).map(|_| rng.vec_f32(n, -1.0, 1.0)).collect();
        let mut w = rng.vec_f32(n, -wscale, wscale);
        let mut b = rng.vec_f32(m, -bscale, bscale);
        for j in 0..n {
            if rng.next_u32() % 5 == 0 {
                w[j] = 0.0;
            }
        }
        match case % 13 {
            0 => w.fill(0.0),
            1 => b.fill(0.0),
            2 => {
                for row in &mut u {
                    row.fill(0.0);
                }
            }
            3 => {
                for row in &mut v {
                    row.fill(0.0);
                }
            }
            _ => (),
        }
        let mut x = vec![-999.0; n];
        let mut cx = x.clone();
        rust::svd_backsubstitute(&u, &w, &v, &b, &mut x);
        unsafe {
            svd_backsubstitute(
                pointers(&mut u).as_mut_ptr(),
                w.as_mut_ptr(),
                pointers(&mut v).as_mut_ptr(),
                m as i32,
                n as i32,
                b.as_mut_ptr(),
                cx.as_mut_ptr(),
            );
        }
        let bits: Vec<u32> = x.iter().map(|x| x.to_bits()).collect();
        let cbits: Vec<u32> = cx.iter().map(|x| x.to_bits()).collect();
        assert_eq!(bits, cbits, "back substitution case {case}");
    }
}

#[test]
fn svd_sort_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0x5077);
    for case in 0..4000 {
        let n = [0, 1, 2, 3, 5, 7, 16, 31, 32, 64, 127][case % 11];
        let m = [0, 1, 2, 7, 16, 31, 64, 131][case % 8];
        // Guarded black-box probes establish that C reads/writes column n.
        // Include it explicitly, plus untouched sentinel columns and V rows.
        let mut u: Vec<Vec<f32>> = (0..m).map(|_| rng.vec_f32(n + 3, -10.0, 10.0)).collect();
        let mut v: Vec<Vec<f32>> = (0..n + 2)
            .map(|_| rng.vec_f32(n + 3, -10.0, 10.0))
            .collect();
        let mut w: Vec<f32> = (0..n + 3).map(|_| rng.range_i32(-3, 3) as f32).collect();
        if case % 6 == 0 {
            w.fill(0.0);
        }
        if case % 6 == 1 {
            w.fill(2.0);
        }
        if case % 6 == 2 {
            w.sort_by(|a, b| a.total_cmp(b));
        }
        if case % 6 == 3 {
            w.sort_by(|a, b| b.total_cmp(a));
        }
        if case % 6 == 4 {
            w = rng.vec_f32(n + 3, -1e12, 1e12);
        }
        if case % 6 == 5 {
            for x in &mut w {
                if *x == 0.0 {
                    *x = -0.0;
                }
            }
        }
        let mut cu = u.clone();
        let mut cw = w.clone();
        let mut cv = v.clone();
        rust::svd_sort(&mut u, &mut w, &mut v, n);
        unsafe {
            svd_sort(
                pointers(&mut cu).as_mut_ptr(),
                cw.as_mut_ptr(),
                pointers(&mut cv).as_mut_ptr(),
                n as i32,
                m as i32,
            );
        }
        assert_eq!(
            w.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            cw.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            "w case {case}"
        );
        assert_eq!(u, cu, "u case {case}");
        assert_eq!(v, cv, "v case {case}");
    }
}

#[test]
fn svd_complete_least_squares_pipeline_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0x5017_bacc);
    for case in 0..2000 {
        let n = [1, 2, 3, 5, 7, 11, 16, 31, 32, 64][case % 10];
        let m = n + case % 17;
        let scale = [1e-12, 1.0, 1e12][case % 3];
        let mut a: Vec<Vec<f32>> = (0..m)
            .map(|_| {
                let mut row = rng.vec_f32(n, -scale, scale);
                row.push(0.0);
                row
            })
            .collect();
        for j in 0..n {
            a[j][j] += 3.0 * n as f32 * scale;
        }
        let mut expected = rng.vec_f32(n, -10.0, 10.0);
        if case % 11 == 0 {
            expected.fill(0.0);
        }
        let mut b = vec![0.0; m];
        for i in 0..m {
            for j in 0..n {
                b[i] += a[i][j] * expected[j];
            }
        }
        let mut ca = a.clone();
        // C's inclusive sort needs explicit padding. A negative-infinity
        // sentinel keeps the extra column outside the finite spectrum.
        let mut w = vec![0.0; n + 1];
        w[n] = f32::NEG_INFINITY;
        let mut cw = w.clone();
        let mut v = vec![vec![0.0; n + 1]; n];
        let mut cv = v.clone();
        assert_eq!(rust::compute_svd(&mut a, &mut w[..n], &mut v), 1);
        assert_eq!(
            unsafe {
                compute_svd(
                    pointers(&mut ca).as_mut_ptr(),
                    m as i32,
                    n as i32,
                    cw.as_mut_ptr(),
                    pointers(&mut cv).as_mut_ptr(),
                )
            },
            1
        );
        rust::svd_sort(&mut a, &mut w, &mut v, n);
        unsafe {
            svd_sort(
                pointers(&mut ca).as_mut_ptr(),
                cw.as_mut_ptr(),
                pointers(&mut cv).as_mut_ptr(),
                n as i32,
                m as i32,
            );
        }
        assert_eq!(w[n], f32::NEG_INFINITY);
        assert_eq!(cw[n], f32::NEG_INFINITY);
        let mut x = vec![0.0; n];
        let mut cx = x.clone();
        rust::svd_backsubstitute(&a, &w[..n], &v, &b, &mut x);
        unsafe {
            svd_backsubstitute(
                pointers(&mut ca).as_mut_ptr(),
                cw.as_mut_ptr(),
                pointers(&mut cv).as_mut_ptr(),
                m as i32,
                n as i32,
                b.as_mut_ptr(),
                cx.as_mut_ptr(),
            );
        }
        let x: Vec<f64> = x.into_iter().map(f64::from).collect();
        let cx: Vec<f64> = cx.into_iter().map(f64::from).collect();
        let expected: Vec<f64> = expected.into_iter().map(f64::from).collect();
        close(&x, &cx, &format!("pipeline {case}: Rust/C solution"));
        close(&x, &expected, &format!("pipeline {case}: known solution"));
        close(
            &cx,
            &expected,
            &format!("pipeline {case}: C known solution"),
        );
    }
}
