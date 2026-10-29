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

//! Differential tests against `cwp/lib/mnewt.c`.
//! Original author: Dave Hale, Colorado School of Mines, 06/06/91.

use oxwave_cref::XorShift;
use oxwave_newton_solver::mnewt as rust_mnewt;
use std::ffi::c_void;

extern "C" {
    fn mnewt(
        maxiter: i32,
        ftol: f32,
        dxtol: f32,
        n: i32,
        x: *mut f32,
        aux: *mut c_void,
        fdfdx: extern "C" fn(i32, *mut f32, *mut f32, *mut *mut f32, *mut c_void),
    ) -> i32;
}

#[derive(Clone)]
enum Mode {
    Linear,
    Cubic(f32),
    Constant,
    // The two-equation example supplied in the original mnewt.c.
    Trig,
}

#[derive(Clone)]
struct Problem {
    mode: Mode,
    a: Vec<Vec<f32>>,
    root: Vec<f32>,
    // Record every callback's input, including the final convergence check.
    trace: Vec<Vec<f32>>,
}

impl Problem {
    fn evaluate(&mut self, x: &mut [f32], f: &mut [f32], dfdx: &mut [&mut [f32]]) {
        self.trace.push(x.to_vec());
        let n = x.len();
        if let Mode::Trig = self.mode {
            // Preserve C's double math calls and promotions in the example.
            f[0] = ((x[0] as f64).sin() - (x[1] as f64).cos()) as f32;
            f[1] = ((x[0] + x[1]) as f64 - 1.0) as f32;
            dfdx[0][0] = (x[0] as f64).cos() as f32;
            dfdx[0][1] = 1.0;
            dfdx[1][0] = (x[1] as f64).sin() as f32;
            dfdx[1][1] = 1.0;
            return;
        }
        for i in 0..n {
            if let Mode::Constant = self.mode {
                f[i] = self.root[i];
                for j in 0..n {
                    dfdx[j][i] = self.a[j][i];
                }
            } else {
                f[i] = 0.0;
                for j in 0..n {
                    let y = x[j] - self.root[j];
                    let (value, derivative) = match self.mode {
                        Mode::Cubic(c) => (y + c * y * y * y, 1.0f32 + 3.0f32 * c * y * y),
                        _ => (y, 1.0f32),
                    };
                    f[i] += self.a[j][i] * value;
                    dfdx[j][i] = self.a[j][i] * derivative;
                }
            }
        }
    }
}

extern "C" fn c_fdfdx(n: i32, x: *mut f32, f: *mut f32, dfdx: *mut *mut f32, aux: *mut c_void) {
    // mnewt allocates n disjoint columns and residuals and passes back the
    // live Problem supplied by check_case. No references escape this call.
    unsafe {
        let n = n as usize;
        let x = std::slice::from_raw_parts_mut(x, n);
        let f = std::slice::from_raw_parts_mut(f, n);
        let columns = std::slice::from_raw_parts(dfdx, n);
        let mut columns: Vec<&mut [f32]> = columns
            .iter()
            .map(|&p| std::slice::from_raw_parts_mut(p, n))
            .collect();
        let problem = &mut *aux.cast::<Problem>();
        problem.evaluate(x, f, &mut columns);
    }
}

fn assert_values(rust: &[f32], c: &[f32], context: &str) {
    assert_eq!(rust.len(), c.len());
    let mut magnitude = 0.0f64;
    for (&a, &b) in rust.iter().zip(c) {
        if a.is_finite() && b.is_finite() {
            magnitude = magnitude.max((a as f64).abs()).max((b as f64).abs());
        }
    }
    for (i, (&a, &b)) in rust.iter().zip(c).enumerate() {
        if a.is_nan() || b.is_nan() {
            assert!(
                a.is_nan() && b.is_nan(),
                "{context}, element {i}: {a} vs {b}"
            );
        } else {
            if a.is_finite() && b.is_finite() {
                assert!(
                    (a as f64 - b as f64).abs() <= 1e-5 * magnitude,
                    "{context}, element {i}: Rust {a}, C {b}, scale {magnitude}"
                );
            }
            // Stronger than the required tolerance: every non-NaN value,
            // including signed zero and infinity, must match bit for bit.
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "{context}, element {i}: {a} vs {b}"
            );
        }
    }
}

fn check_case(
    problem: Problem,
    initial: &[f32],
    maxiter: i32,
    ftol: f32,
    dxtol: f32,
) -> (i32, Vec<f32>, usize) {
    let mut c_problem = problem.clone();
    let mut rust_problem = problem;
    let mut cx = initial.to_vec();
    let mut rx = initial.to_vec();
    let c_status = unsafe {
        mnewt(
            maxiter,
            ftol,
            dxtol,
            cx.len() as i32,
            cx.as_mut_ptr(),
            (&mut c_problem as *mut Problem).cast(),
            c_fdfdx,
        )
    };
    let rust_status = rust_mnewt(
        maxiter,
        ftol,
        dxtol,
        &mut rx,
        &mut rust_problem,
        |x, f, j, p| {
            p.evaluate(x, f, j);
        },
    );
    assert_eq!(
        rust_status,
        c_status,
        "n={}, limit={maxiter}, ftol={ftol}, dxtol={dxtol}, x={initial:?}",
        initial.len()
    );
    assert_eq!(
        rust_problem.trace.len(),
        c_problem.trace.len(),
        "callback count"
    );
    assert_values(&rx, &cx, "final x");
    for (r, c) in rust_problem.trace.iter().zip(&c_problem.trace) {
        assert_values(r, c, "iteration x");
    }
    (rust_status, rx, rust_problem.trace.len())
}

fn diagonal(diag: &[f32], residual: &[f32]) -> Problem {
    let n = diag.len();
    let mut a = vec![vec![0.0; n]; n];
    for j in 0..n {
        a[j][j] = diag[j];
    }
    Problem {
        mode: Mode::Constant,
        a,
        root: residual.to_vec(),
        trace: Vec::new(),
    }
}

fn random_problem(rng: &mut XorShift, n: usize, mode: Mode) -> Problem {
    let mut a = vec![vec![0.0f32; n]; n];
    // Rotate rows of a diagonally dominant, nonsymmetric matrix to force
    // pivoting and expose a transposed Jacobian in either implementation.
    let shift = rng.range_i32(0, n as i32 - 1) as usize;
    let scale = [1e-12f32, 1e-4, 1.0, 1e4, 1e12][rng.next_u32() as usize % 5];
    for j in 0..n {
        for i in 0..n {
            a[j][(i + shift) % n] = scale
                * if i == j {
                    let sign = if rng.next_u32() & 1 == 0 { -1.0 } else { 1.0 };
                    sign * (n as f32 + rng.uniform_f32(1.0, 3.0))
                } else {
                    rng.uniform_f32(-0.5, 0.5)
                };
        }
    }
    Problem {
        mode,
        a,
        root: rng.vec_f32(n, -10.0, 10.0),
        trace: Vec::new(),
    }
}

#[test]
fn mnewt_dense_linear_matches_c() {
    let mut rng = XorShift::new(0x6d6e65777401);
    let sizes = [1, 2, 3, 5, 7, 8, 11, 16];
    for case in 0..3200 {
        let n = sizes[case % sizes.len()];
        let problem = random_problem(&mut rng, n, Mode::Linear);
        let mut initial = rng.vec_f32(n, -100.0, 100.0);
        if case % 11 == 0 {
            initial.fill(0.0);
        }
        if case % 13 == 0 {
            initial.clone_from(&problem.root);
        }
        let maxiter = [0, 1, 2, 5, 10][case % 5];
        let ftol = [0.0, 1e-20, 1e-6, 1e-3, -1.0][case % 5];
        let dxtol = [0.0, 1e-6, 1e-4, 1000.0, -1.0][(case / 5) % 5];
        check_case(problem, &initial, maxiter, ftol, dxtol);
    }
    // Larger primes, powers of two, and adjacent sizes exercise LU's long
    // columns and BLAS unrolling without making thousands of large solves.
    for n in [17, 31, 32, 63, 64, 97, 127, 128] {
        for _ in 0..8 {
            let p = random_problem(&mut rng, n, Mode::Linear);
            let initial = rng.vec_f32(n, -100.0, 100.0);
            check_case(p, &initial, 6, 1e-6, 1e-5);
        }
    }
}

#[test]
fn mnewt_coupled_nonlinear_matches_c() {
    let mut rng = XorShift::new(0x6d6e65777402);
    let sizes = [1, 2, 3, 5, 7, 8, 13, 16];
    let mut converged = 0;
    let mut exhausted = 0;
    for case in 0..3000 {
        let n = sizes[case % sizes.len()];
        let curvature = [0.01, 0.25, 1.0, 10.0][case % 4];
        let p = random_problem(&mut rng, n, Mode::Cubic(curvature));
        let mut initial = rng.vec_f32(n, -25.0, 25.0);
        if case % 17 == 0 {
            initial.fill(-0.0);
        }
        let maxiter = [1, 2, 8, 20, 40][case % 5];
        let ftol = [0.0, 1e-8, 1e-4][case % 3];
        let dxtol = [0.0, 1e-5, 1e-3][(case / 3) % 3];
        let (status, _, _) = check_case(p, &initial, maxiter, ftol, dxtol);
        if status >= 0 {
            converged += 1;
        } else {
            exhausted += 1;
        }
    }
    assert!(converged > 100 && exhausted > 100);
}

#[test]
fn mnewt_original_trigonometric_example_matches_c() {
    let mut rng = XorShift::new(0x6d6e65777403);
    for case in 0..1200 {
        let p = Problem {
            mode: Mode::Trig,
            a: Vec::new(),
            root: Vec::new(),
            trace: Vec::new(),
        };
        let initial = rng.vec_f32(2, -10.0, 10.0);
        check_case(p, &initial, [1, 10, 50][case % 3], 1e-6, 1e-6);
    }
}

#[test]
fn mnewt_singular_jacobians_match_c() {
    let mut rng = XorShift::new(0x6d6e65777404);
    for case in 0..1200 {
        let n = [1, 2, 3, 5, 8, 11, 16][case % 7];
        let mut diag = rng.vec_f32(n, 0.5, 10.0);
        // Hit first, middle and last singular pivots, including signed zero.
        let pivot = rng.range_i32(0, n as i32 - 1) as usize;
        diag[pivot] = if case % 2 == 0 { 0.0 } else { -0.0 };
        let residual = rng.vec_f32(n, -5.0, 5.0);
        let p = diagonal(&diag, &residual);
        let initial = rng.vec_f32(n, -5.0, 5.0);
        let (status, x, calls) = check_case(p, &initial, 10, -1.0, f32::INFINITY);
        assert_eq!(status, -1);
        assert_eq!(calls, 1);
        assert_values(&x, &initial, "singular leaves x unchanged");
    }
    // Residual convergence precedes factorization even for an all-zero matrix.
    let (status, _, calls) = check_case(diagonal(&[0.0; 5], &[0.0; 5]), &[3.0; 5], 10, 0.0, -1.0);
    assert_eq!((status, calls), (0, 1));
}

#[test]
fn mnewt_stopping_order_and_boundaries_match_c() {
    let p = diagonal(&[2.0, -4.0], &[0.5, 2.0]);
    let initial = [7.0, -9.0];
    // sum(abs(f)) == 2.5, while sum(abs(dx)) == 0.75.
    let (status, x, calls) = check_case(p.clone(), &initial, 5, 2.5, -1.0);
    assert_eq!((status, calls), (0, 1));
    assert_values(&x, &initial, "ftol equality");
    let (status, x, calls) = check_case(p.clone(), &initial, 5, -1.0, 0.75);
    assert_eq!((status, calls), (0, 1));
    assert_values(&x, &initial, "dxtol equality before update");
    let below_step = f32::from_bits(0.75f32.to_bits() - 1);
    let (status, x, calls) = check_case(p.clone(), &initial, 1, -1.0, below_step);
    assert_eq!((status, calls), (-1, 1));
    assert_values(&x, &[6.75, -8.5], "exhaustion retains last update");
    let below_residual = f32::from_bits(2.5f32.to_bits() - 1);
    let (status, _, _) = check_case(p.clone(), &initial, 1, below_residual, below_step);
    assert_eq!(status, -1);
    for maxiter in [i32::MIN, -10, -1, 0] {
        let (status, x, calls) = check_case(p.clone(), &initial, maxiter, 0.0, 0.0);
        assert_eq!(status, if maxiter == 0 { -1 } else { 0 });
        assert_eq!(calls, 0);
        assert_values(&x, &initial, "no iteration");
    }
    for (ftol, dxtol, limit, expected, calls) in [
        (0.0, -1.0, 5, 0, 1),
        (-1.0, 0.0, 5, 0, 1),
        (-1.0, -1.0, 5, -1, 5),
    ] {
        let (status, _, actual_calls) = check_case(
            diagonal(&[1.0, 1.0], &[0.0, -0.0]),
            &[0.0, -0.0],
            limit,
            ftol,
            dxtol,
        );
        assert_eq!((status, actual_calls), (expected, calls));
    }
    // Convergence on the second callback returns 1, not a callback count of 2.
    let p = Problem {
        mode: Mode::Linear,
        a: vec![vec![2.0]],
        root: vec![3.0],
        trace: Vec::new(),
    };
    let (status, x, calls) = check_case(p, &[1.0], 2, 0.0, 0.0);
    assert_eq!((status, calls), (1, 2));
    assert_values(&x, &[3.0], "one Newton update");
}

#[test]
fn mnewt_extreme_values_match_c() {
    let mut rng = XorShift::new(0x6d6e65777405);
    let values = [
        0.0,
        -0.0,
        f32::from_bits(1),
        -f32::from_bits(1),
        f32::MIN_POSITIVE,
        -f32::MIN_POSITIVE,
        1e-30,
        -1e-30,
        1e-10,
        -1e-10,
        1.0,
        -1.0,
        1e10,
        -1e10,
        1e30,
        -1e30,
        f32::MAX,
        -f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ];
    for case in 0..2000 {
        let n = [1, 2, 3, 5, 8][case % 5];
        let mut diag = Vec::new();
        let mut residual = Vec::new();
        let mut initial = Vec::new();
        for _ in 0..n {
            diag.push(values[rng.next_u32() as usize % values.len()]);
            residual.push(values[rng.next_u32() as usize % values.len()]);
            initial.push(values[rng.next_u32() as usize % values.len()]);
        }
        let ftol = values[case % values.len()];
        let dxtol = values[(case / values.len()) % values.len()];
        check_case(diagonal(&diag, &residual), &initial, 3, ftol, dxtol);
    }
}
