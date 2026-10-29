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

//! Clean-room differential tests for par/lib/minfunc.c.
//! Credits: Colorado School of Mines (2011); the original notes adaptation
//! from Numerical Recipes. No original implementation bodies are used here.

use oxwave_cref::XorShift;
use oxwave_f_minfunc as rust;
use std::cell::RefCell;

// The original line minimizer has process-global callback state.
static C_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
type Scalar = extern "C" fn(f32) -> f32;
type Vector = extern "C" fn(*const f32) -> f32;
extern "C" {
    fn bracket_minimum(
        a: *mut f32,
        b: *mut f32,
        c: *mut f32,
        fa: *mut f32,
        fb: *mut f32,
        fc: *mut f32,
        f: Scalar,
    );
    fn golden_bracket(a: f32, b: f32, c: f32, f: Scalar, tol: f32, xmin: *mut f32) -> f32;
    fn brent_bracket(a: f32, b: f32, c: f32, f: Scalar, tol: f32, xmin: *mut f32) -> f32;
    fn linmin(p: *mut f32, xi: *mut f32, n: i32, fret: *mut f32, f: Vector);
    fn oneDfunction(x: f32) -> f32;
    fn powell_minimization(
        p: *mut f32,
        xi: *mut *mut f32,
        n: i32,
        ftol: f32,
        iter: *mut i32,
        fret: *mut f32,
        f: Vector,
    );
}

#[derive(Clone, Debug, Default)]
struct Objective {
    center: Vec<f32>,
    weight: Vec<f32>,
    kind: usize,
    offset: f32,
    coupling: f32,
    // Probe oneDfunction only while linmin has initialized its live globals.
    probes: Vec<f32>,
    probe_results: Vec<f32>,
    probing: bool,
}
thread_local! { static OBJECTIVE: RefCell<Objective> = RefCell::new(Objective::default()); }

fn shape(x: f32, kind: usize) -> f32 {
    match kind {
        0 => x * x,
        1 => x * x + 0.125f32 * x * x * x * x,
        2 => x.abs() + 0.25f32 * x * x,
        _ => {
            if x < 0.0 {
                0.5f32 * x * x
            } else {
                2.0f32 * x * x
            }
        }
    }
}
fn evaluate(p: &[f32], o: &Objective) -> f32 {
    let mut v = 0.0;
    let mut total = 0.0;
    for j in 0..p.len() {
        let d = p[j] - o.center[j];
        v += o.weight[j] * shape(d, o.kind);
        total += d;
    }
    v + o.coupling * total * total + o.offset
}
fn objective(p: &[f32]) -> f32 {
    OBJECTIVE.with(|o| evaluate(p, &o.borrow()))
}
extern "C" fn scalar(x: f32) -> f32 {
    objective(&[x])
}
extern "C" fn vector(p: *const f32) -> f32 {
    let (value, probes) = OBJECTIVE.with(|o| {
        let mut o = o.borrow_mut();
        let value = evaluate(unsafe { std::slice::from_raw_parts(p, o.center.len()) }, &o);
        let probes = if !o.probing && o.probe_results.is_empty() {
            o.probes.clone()
        } else {
            vec![]
        };
        if !probes.is_empty() {
            o.probing = true;
        }
        (value, probes)
    });
    if !probes.is_empty() {
        let results = probes.iter().map(|&x| unsafe { oneDfunction(x) }).collect();
        OBJECTIVE.with(|o| {
            let mut o = o.borrow_mut();
            o.probe_results = results;
            o.probing = false;
        });
    }
    value
}
fn setup(r: &mut XorShift, case: usize, n: usize) -> Objective {
    let scale = [0.01, 1.0, 100.0][case % 3];
    let o = Objective {
        center: (0..n)
            .map(|j| {
                if case % 19 == 0 {
                    0.0
                } else {
                    r.uniform_f32(-5.0, 5.0) * scale + (j % 3) as f32 * scale
                }
            })
            .collect(),
        weight: (0..n)
            .map(|_| r.uniform_f32(0.25, 4.0) / (scale * scale))
            .collect(),
        kind: case % 4,
        // A cusp at exactly zero with zero baseline cannot satisfy the C
        // golden search's relative-only stopping rule at subnormal spacing.
        offset: if case % 19 == 0 && case % 4 == 2 {
            1.0
        } else {
            [0.0, 1.0, -2.0][case % 3]
        },
        coupling: if case % 5 == 0 {
            0.125 / (scale * scale)
        } else {
            0.0
        },
        ..Objective::default()
    };
    OBJECTIVE.with(|c| *c.borrow_mut() = o.clone());
    o
}
fn compare(a: &[f32], b: &[f32], context: &str) {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().chain(b).fold(1e-30f32, |s, v| s.max(v.abs()));
    for (j, (&a, &b)) in a.iter().zip(b).enumerate() {
        assert!(
            a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-5 * scale,
            "{context} element {j}: C {a:?}, Rust {b:?}, max magnitude {scale}"
        );
    }
}
fn bracket_c(a: f32, b: f32) -> [f32; 6] {
    let [mut a, mut b, mut c, mut fa, mut fb, mut fc] = [a, b, 0.0, 0.0, 0.0, 0.0];
    unsafe {
        bracket_minimum(&mut a, &mut b, &mut c, &mut fa, &mut fb, &mut fc, scalar);
    }
    [a, b, c, fa, fb, fc]
}
fn bracket_r(a: f32, b: f32) -> [f32; 6] {
    let [mut a, mut b, mut c, mut fa, mut fb, mut fc] = [a, b, 0.0, 0.0, 0.0, 0.0];
    rust::bracket_minimum(&mut a, &mut b, &mut c, &mut fa, &mut fb, &mut fc, |x| {
        scalar(x)
    });
    [a, b, c, fa, fb, fc]
}

#[test]
fn bracket_minimum_matches_c() {
    let mut r = XorShift::new(11);
    for case in 0..6000 {
        let o = setup(&mut r, case, 1);
        let a = r.uniform_f32(-8.0, 8.0);
        let b = if case % 31 == 0 {
            a
        } else {
            a + r.uniform_f32(0.01, 4.0) * if case % 2 == 0 { -1.0 } else { 1.0 }
        };
        let c = bracket_c(a, b);
        let d = bracket_r(a, b);
        let context = format!("bracket case {case}, {a}, {b}, {o:?}");
        compare(&c[..3], &d[..3], &context);
        compare(&c[3..], &d[3..], &context);
    }
}

fn scalar_minima(golden: bool) {
    let mut r = XorShift::new(if golden { 12 } else { 13 });
    for case in 0..6000 {
        let o = setup(&mut r, case, 1);
        let width = r.uniform_f32(0.1, 20.0);
        let a = o.center[0] - width;
        let c = o.center[0] + 1.75 * width;
        let b = o.center[0] + 0.1 * width;
        let (a, c) = if case % 2 == 0 { (a, c) } else { (c, a) };
        let tol = [1e-2, 1e-3, 1e-4][case % 3];
        let (mut x, mut y) = (0.0, 0.0);
        let (v, w) = if golden {
            (
                unsafe { golden_bracket(a, b, c, scalar, tol, &mut x) },
                rust::golden_bracket(a, b, c, |x| scalar(x), tol, &mut y),
            )
        } else {
            (
                unsafe { brent_bracket(a, b, c, scalar, tol, &mut x) },
                rust::brent_bracket(a, b, c, |x| scalar(x), tol, &mut y),
            )
        };
        let context = format!("golden={golden}, case {case}, {o:?}, bracket {a} {b} {c}");
        compare(&[x], &[y], &context);
        compare(&[v], &[w], &context);
    }
}
#[test]
fn golden_bracket_matches_c() {
    scalar_minima(true);
}
#[test]
fn brent_bracket_matches_c() {
    scalar_minima(false);
}

const SIZES: &[usize] = &[1, 2, 3, 5, 7, 8, 16, 31, 64];
// Line evaluation is linear in the dimension; also exercise large prime and
// power-of-two dimensions without making Powell's dense sweeps cubic in 1024.
const LINE_SIZES: &[usize] = &[1, 2, 3, 5, 7, 8, 16, 31, 64, 127, 256, 1024];
fn line_input(r: &mut XorShift, case: usize, n: usize) -> (Vec<f32>, Vec<f32>) {
    let p = (0..n)
        .map(|_| {
            if case % 23 == 0 {
                0.0
            } else {
                r.uniform_f32(-10.0, 10.0)
            }
        })
        .collect();
    let xi = (0..n)
        .map(|_| {
            if case % 29 == 0 {
                0.0
            } else {
                r.uniform_f32(-2.0, 2.0)
            }
        })
        .collect();
    (p, xi)
}
#[test]
fn linmin_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(14);
    for case in 0..3000 {
        let n = LINE_SIZES[(case / 4) % LINE_SIZES.len()];
        let o = setup(&mut r, case, n);
        let (mut p, mut xi) = line_input(&mut r, case, n);
        let (mut q, mut yi) = (p.clone(), xi.clone());
        let (mut a, mut b) = (0.0, 0.0);
        unsafe {
            linmin(p.as_mut_ptr(), xi.as_mut_ptr(), n as i32, &mut a, vector);
        }
        rust::linmin(&mut q, &mut yi, &mut b, objective);
        let context = format!("linmin case {case}, {o:?}");
        compare(&p, &q, &context);
        compare(&xi, &yi, &context);
        compare(&[a], &[b], &context);
    }
}

#[test]
fn one_d_function_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(15);
    for case in 0..2000 {
        let n = LINE_SIZES[(case / 4) % LINE_SIZES.len()];
        setup(&mut r, case, n);
        let (mut p, mut xi) = line_input(&mut r, case, n);
        let (base, direction) = (p.clone(), xi.clone());
        let probes = vec![
            0.0,
            -0.0,
            1.0,
            -1.0,
            r.uniform_f32(-100.0, 100.0),
            r.uniform_f32(-0.001, 0.001),
        ];
        OBJECTIVE.with(|o| o.borrow_mut().probes = probes.clone());
        let mut value = 0.0;
        unsafe {
            linmin(
                p.as_mut_ptr(),
                xi.as_mut_ptr(),
                n as i32,
                &mut value,
                vector,
            );
        }
        let c = OBJECTIVE.with(|o| o.borrow().probe_results.clone());
        assert_eq!(c.len(), probes.len());
        for (&x, &expected) in probes.iter().zip(&c) {
            let actual = rust::oneDfunction(x, &base, &direction, objective);
            // Both paths must form the identical single-precision point.
            assert_eq!(
                expected.to_bits(),
                actual.to_bits(),
                "oneD case {case}, x={x}"
            );
        }
    }
}

#[test]
fn powell_minimization_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(16);
    for case in 0..2200 {
        let n = SIZES[case % SIZES.len()];
        let o = setup(&mut r, case, n);
        let (mut p, _) = line_input(&mut r, case, n);
        let mut q = p.clone();
        let mut xi = vec![vec![0.0; n]; n];
        for i in 0..n {
            xi[i][i] = r.uniform_f32(0.5, 2.0);
            if case % 7 == 0 && i + 1 < n {
                xi[i][i + 1] = 0.1;
            }
        }
        let mut yi = xi.clone();
        let mut rows: Vec<_> = xi.iter_mut().map(|row| row.as_mut_ptr()).collect();
        let (mut a, mut b, mut it, mut jt) = (0.0, 0.0, 0, 0);
        let tol = [1e-2, 1e-3, 1e-4][case % 3];
        unsafe {
            powell_minimization(
                p.as_mut_ptr(),
                rows.as_mut_ptr(),
                n as i32,
                tol,
                &mut it,
                &mut a,
                vector,
            );
        }
        rust::powell_minimization(&mut q, &mut yi, tol, &mut jt, &mut b, objective);
        let context = format!("powell case {case}, {o:?}, C it={it}, Rust it={jt}");
        compare(&p, &q, &context);
        compare(&[a], &[b], &context);
        for i in 0..n {
            compare(&xi[i], &yi[i], &context);
        }
        assert_eq!(it, jt, "{context}");
    }
}

#[test]
fn powell_round_two_regression() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    OBJECTIVE.with(|o| {
        *o.borrow_mut() = Objective {
            center: vec![3.25, -2.0],
            weight: vec![1.0, 1.0],
            offset: 1.0,
            ..Objective::default()
        }
    });
    let (mut p, mut q) = ([0.0, 0.0], [0.0, 0.0]);
    let mut xi = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let mut yi = xi.clone();
    let mut rows: Vec<_> = xi.iter_mut().map(|row| row.as_mut_ptr()).collect();
    let (mut a, mut b, mut it, mut jt) = (0.0, 0.0, 0, 0);
    unsafe {
        powell_minimization(
            p.as_mut_ptr(),
            rows.as_mut_ptr(),
            2,
            1e-4,
            &mut it,
            &mut a,
            vector,
        );
    }
    rust::powell_minimization(&mut q, &mut yi, 1e-4, &mut jt, &mut b, objective);
    compare(&p, &q, "round two Powell regression");
    compare(&[a], &[b], "round two Powell regression objective");
    assert_eq!(it, jt);
    for i in 0..2 {
        compare(&xi[i], &yi[i], "round two Powell regression directions");
    }
}

#[test]
fn constant_objectives_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(17);
    for case in 0..2000 {
        // Constant objectives exercise equality branches and immediate Powell
        // convergence, including an objective that is identically zero.
        let value = match case % 3 {
            0 => 0.0,
            1 => r.uniform_f32(-100.0, -0.01),
            _ => r.uniform_f32(0.01, 100.0),
        };
        OBJECTIVE.with(|o| {
            *o.borrow_mut() = Objective {
                center: vec![0.0],
                weight: vec![0.0],
                offset: value,
                ..Objective::default()
            };
        });
        let a = r.uniform_f32(-20.0, -1.0);
        let c = r.uniform_f32(1.0, 20.0);
        let (a, c) = if case % 2 == 0 { (a, c) } else { (c, a) };
        let context = format!("constant objective case {case}, value {value}");
        let expected = bracket_c(a, c);
        let actual = bracket_r(a, c);
        compare(&expected[..3], &actual[..3], &context);
        compare(&expected[3..], &actual[3..], &context);
        let tol = [1e-2, 1e-3, 1e-4][case % 3];
        for golden in [false, true] {
            let (mut x, mut y) = (0.0, 0.0);
            let (expected, actual) = if golden {
                (
                    unsafe { golden_bracket(a, 0.0, c, scalar, tol, &mut x) },
                    rust::golden_bracket(a, 0.0, c, |x| scalar(x), tol, &mut y),
                )
            } else {
                (
                    unsafe { brent_bracket(a, 0.0, c, scalar, tol, &mut x) },
                    rust::brent_bracket(a, 0.0, c, |x| scalar(x), tol, &mut y),
                )
            };
            compare(&[x], &[y], &context);
            compare(&[expected], &[actual], &context);
        }

        let n = SIZES[case % SIZES.len()];
        OBJECTIVE.with(|o| {
            let mut o = o.borrow_mut();
            o.center = vec![0.0; n];
            o.weight = vec![0.0; n];
        });
        let (base, direction) = line_input(&mut r, case, n);
        let (mut p, mut q) = (base.clone(), base.clone());
        let (mut xi, mut yi) = (direction.clone(), direction);
        let (mut expected, mut actual) = (0.0, 0.0);
        unsafe {
            linmin(
                p.as_mut_ptr(),
                xi.as_mut_ptr(),
                n as i32,
                &mut expected,
                vector,
            );
        }
        rust::linmin(&mut q, &mut yi, &mut actual, objective);
        compare(&p, &q, &context);
        compare(&xi, &yi, &context);
        compare(&[expected], &[actual], &context);

        let (mut p, mut q) = (base.clone(), base);
        let mut xi = vec![vec![0.0; n]; n];
        for i in 0..n {
            xi[i][i] = r.uniform_f32(-2.0, -0.5);
        }
        let mut yi = xi.clone();
        let mut rows: Vec<_> = xi.iter_mut().map(|row| row.as_mut_ptr()).collect();
        let (mut it, mut jt) = (0, 0);
        unsafe {
            powell_minimization(
                p.as_mut_ptr(),
                rows.as_mut_ptr(),
                n as i32,
                tol,
                &mut it,
                &mut expected,
                vector,
            );
        }
        rust::powell_minimization(&mut q, &mut yi, tol, &mut jt, &mut actual, objective);
        compare(&p, &q, &context);
        compare(&[expected], &[actual], &context);
        assert_eq!(it, jt, "{context}");
        assert_eq!(it, 1, "{context}");
        for i in 0..n {
            compare(&xi[i], &yi[i], &context);
        }
    }
}
