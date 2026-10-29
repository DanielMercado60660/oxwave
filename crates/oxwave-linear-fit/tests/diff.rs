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

//! Differential tests against `cwp/lib/linfit.c` (CSM, 2011).
//! Author of the C routines: Balasz Nemeth (inspired by similar Numerical
//! Recipes code, including `fit.c`).

use oxwave_cref::XorShift;
use oxwave_linear_fit as rust;
use std::process::Command;

extern "C" {
    fn gammln(xx: f32) -> f32;
    fn gammq(a: f32, x: f32) -> f32;
    fn gcf(gammcf: *mut f32, a: f32, x: f32, gln: *mut f32);
    fn gser(gamser: *mut f32, a: f32, x: f32, gln: *mut f32);
    fn linfit(
        x: *mut f32,
        y: *mut f32,
        ndata: i32,
        sig: *mut f32,
        mwt: i32,
        a: *mut f32,
        b: *mut f32,
        siga: *mut f32,
        sigb: *mut f32,
        chi2: *mut f32,
        q: *mut f32,
    );
}

// Use a scalar-relative bound for each output, including tiny probabilities.
// NaNs and infinities must match classification; no nonfinite case is skipped.
fn same(actual: f32, expected: f32, context: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{context}: Rust {actual:?}, C NaN");
    } else if expected.is_infinite() || expected == 0.0 {
        assert_eq!(actual.to_bits(), expected.to_bits(), "{context}");
    } else {
        let error = (actual as f64 - expected as f64).abs();
        let scale = (actual as f64).abs().max((expected as f64).abs());
        assert!(
            error <= 1e-5 * scale,
            "{context}: Rust {actual:?}, C {expected:?}, relative error {}",
            error / scale
        );
    }
}

fn positive_log(rng: &mut XorShift, lo: f64, hi: f64) -> f32 {
    10.0_f64.powf(rng.uniform(lo, hi)) as f32
}

#[test]
fn gammln_matches_c() {
    let mut rng = XorShift::new(0x10_600d);
    for i in 0..24_000 {
        let x = match i % 4 {
            0 => rng.uniform_f32(0.001, 200.0),
            1 => positive_log(&mut rng, -44.0, 38.0),
            2 => rng.uniform_f32(-100.0, 0.0),
            _ => rng.any_f32(),
        };
        same(
            rust::gammln(x),
            unsafe { gammln(x) },
            &format!("gammln({x:?})"),
        );
    }
    for x in [
        0.0,
        -0.0,
        1.0,
        2.0,
        0.5,
        -0.5,
        -1.0,
        -2.0,
        -5.5,
        -6.0,
        f32::from_bits(1),
        f32::MIN_POSITIVE,
        f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ] {
        same(
            rust::gammln(x),
            unsafe { gammln(x) },
            &format!("gammln({x:?}) edge"),
        );
    }
}

fn check_gammq(a: f32, x: f32) {
    same(
        rust::gammq(a, x),
        unsafe { gammq(a, x) },
        &format!("gammq({a:?}, {x:?})"),
    );
}

#[test]
fn gammq_matches_c() {
    let mut rng = XorShift::new(0x20_600d);
    for i in 0..16_000 {
        // a <= 96 keeps even the series at the branch boundary within ITMAX.
        let a = match i % 3 {
            0 => positive_log(&mut rng, -8.0, 1.9),
            _ => rng.uniform_f32(0.01, 96.0),
        };
        let boundary = (a as f64 + 1.0) as f32;
        let x = match i % 8 {
            0 => 0.0,
            1 => -0.0,
            2 => rng.uniform_f32(0.0, a),
            3 => boundary,
            4 => f32::from_bits(boundary.to_bits() - 1),
            5 => f32::from_bits(boundary.to_bits() + 1),
            6 => positive_log(&mut rng, -44.0, 38.0),
            _ => rng.uniform_f32(a + 1.0, 1000.0),
        };
        check_gammq(a, x);
    }
    for a in [0.001, 0.5, 1.0, 2.0, 10.0, 1000.0, 10_000.0, 1e6] {
        for x in [0.0, f32::from_bits(1), a * 0.1, a * 2.0, 1e20, f32::MAX] {
            check_gammq(a, x);
        }
    }
}

fn check_gcf(a: f32, x: f32) {
    let (mut c, mut cgln) = (-123.0, -456.0);
    let (mut r, mut rgln) = (c, cgln);
    unsafe { gcf(&mut c, a, x, &mut cgln) };
    rust::gcf(&mut r, a, x, &mut rgln);
    let context = format!("gcf({a:?}, {x:?})");
    same(r, c, &context);
    same(rgln, cgln, &format!("{context} gln"));
}

#[test]
fn gcf_matches_c() {
    let mut rng = XorShift::new(0x30_600d);
    for i in 0..16_000 {
        let a = match i % 4 {
            0 => positive_log(&mut rng, -8.0, 5.0),
            1 => rng.uniform_f32(-30.0, -0.01),
            _ => rng.uniform_f32(0.01, 100.0),
        };
        let x = match i % 5 {
            0 => positive_log(&mut rng, 1.0, 38.0) + a.abs(),
            1 => (a as f64 + 1.0) as f32,
            _ => a.abs() + rng.uniform_f32(1.0, 300.0),
        };
        // For negative a the x=a+1 cases stay positive by construction here.
        let x = if a < 0.0 {
            x.abs() + 1.0
        } else if a > 100.0 {
            x + a
        } else {
            x
        };
        check_gcf(a, x);
    }
    // Integer a terminates the fraction exactly, even outside gammq's branch.
    // Tiny b exercises the FPMIN guard on c and d at the second iteration.
    for a in [2.0, 3.0, 10.0] {
        for x in [f32::from_bits(1), 1e-35, 1e-20, 0.25, 1.5, 20.0, 1e30] {
            check_gcf(a, x);
        }
    }
    for x in [0.25, 1.5, 20.0, 1e30] {
        check_gcf(1.0, x);
    }
    check_gcf(0.0, 1.0);
    check_gcf(-1.0, 2.0);
    check_gcf(1.0, -2.0);
}

fn check_gser(a: f32, x: f32) {
    let (mut c, mut cgln) = (-123.0, -456.0);
    let (mut r, mut rgln) = (c, cgln);
    unsafe { gser(&mut c, a, x, &mut cgln) };
    rust::gser(&mut r, a, x, &mut rgln);
    let context = format!("gser({a:?}, {x:?})");
    same(r, c, &context);
    same(rgln, cgln, &format!("{context} gln"));
}

#[test]
fn gser_matches_c() {
    let mut rng = XorShift::new(0x40_600d);
    for i in 0..16_000 {
        let a = match i % 4 {
            0 => positive_log(&mut rng, -12.0, 1.9),
            1 => rng.uniform_f32(-0.99, -0.01),
            _ => rng.uniform_f32(0.001, 96.0),
        };
        let x = match i % 6 {
            0 => 0.0,
            1 => -0.0,
            2 => positive_log(&mut rng, -44.0, -1.0),
            3 => ((a as f64 + 1.0) as f32).max(0.001),
            _ => rng.uniform_f32(0.0, a.abs() + 0.5),
        };
        check_gser(a, x);
    }
    for a in [0.0, -1.0, -2.0, 1e10, f32::INFINITY, f32::NAN] {
        check_gser(a, 0.0);
    }
    for a in [1.0, 0.5, 1000.0, 1e20] {
        check_gser(a, f32::from_bits(1));
        check_gser(a, 0.01);
    }
}

fn fit_c(x: &[f32], y: &[f32], sig: &[f32], mwt: i32) -> [f32; 6] {
    // Mutable copies honor the original prototypes; C does not change inputs.
    let mut cx = x.to_vec();
    let mut cy = y.to_vec();
    let mut csig = sig.to_vec();
    let (mut a, mut b, mut siga, mut sigb, mut chi2, mut q) = (99.0, 98.0, 97.0, 96.0, 95.0, 94.0);
    unsafe {
        linfit(
            cx.as_mut_ptr(),
            cy.as_mut_ptr(),
            x.len() as i32,
            csig.as_mut_ptr(),
            mwt,
            &mut a,
            &mut b,
            &mut siga,
            &mut sigb,
            &mut chi2,
            &mut q,
        );
    }
    for (actual, expected) in cx
        .iter()
        .chain(&cy)
        .chain(&csig)
        .zip(x.iter().chain(y).chain(sig))
    {
        assert_eq!(actual.to_bits(), expected.to_bits());
    }
    [a, b, siga, sigb, chi2, q]
}

fn fit_rust(x: &[f32], y: &[f32], sig: &[f32], mwt: i32) -> [f32; 6] {
    let (mut a, mut b, mut siga, mut sigb, mut chi2, mut q) = (99.0, 98.0, 97.0, 96.0, 95.0, 94.0);
    rust::linfit(
        x, y, sig, mwt, &mut a, &mut b, &mut siga, &mut sigb, &mut chi2, &mut q,
    );
    [a, b, siga, sigb, chi2, q]
}

fn check_fit(x: &[f32], y: &[f32], sig: &[f32], mwt: i32) {
    let c = fit_c(x, y, sig, mwt);
    let r = fit_rust(x, y, sig, mwt);
    for (i, name) in ["a", "b", "siga", "sigb", "chi2", "q"].iter().enumerate() {
        same(
            r[i],
            c[i],
            &format!("linfit n={} mwt={mwt} {name}", x.len()),
        );
        // Fit arithmetic without transcendentals should be bit-identical.
        if i <= 1 || i == 4 {
            assert!(
                r[i].to_bits() == c[i].to_bits() || (r[i].is_nan() && c[i].is_nan()),
                "linfit n={} mwt={mwt} {name}: Rust {:?}, C {:?}",
                x.len(),
                r[i],
                c[i]
            );
        }
    }
}

#[test]
fn linfit_matches_c() {
    let mut rng = XorShift::new(0x50_600d);
    let sizes = [1, 2, 3, 5, 7, 16, 31, 64, 127, 256, 1024, 4093, 8192];
    for case in 0..4096 {
        let n = sizes[case % sizes.len()];
        let scale = positive_log(&mut rng, -10.0, 10.0);
        let offset = rng.uniform_f32(-3.0, 3.0) * scale;
        let slope = rng.uniform_f32(-3.0, 3.0);
        let noise = match case % 4 {
            0 => 0.0,
            1 => 0.001,
            2 => 0.25,
            _ => 4.0,
        };
        let mut x = Vec::with_capacity(n);
        let mut y = Vec::with_capacity(n);
        let mut sig = Vec::with_capacity(n);
        for _ in 0..n {
            let xi = rng.uniform_f32(-2.0, 2.0) * scale;
            x.push(xi);
            y.push(offset + slope * xi + rng.uniform_f32(-noise, noise) * scale);
            // Keep large fits well away from the nonconvergent series region.
            sig.push(rng.uniform_f32(0.5, 2.0) * scale * if n > 256 { 100.0 } else { 1.0 });
        }
        check_fit(&x, &y, &[], 0);
        if n >= 3 {
            let mwt = [1, -1, 2, i32::MIN, i32::MAX][case % 5];
            check_fit(&x, &y, &sig, mwt);
            for value in &mut sig {
                *value = -*value;
            }
            check_fit(&x, &y, &sig, mwt);
        }
    }
    // Empty, singular, zero and negative coordinates, exact fits, cancellation,
    // underflow/overflow, NaNs and infinities all have defined IEEE outcomes.
    check_fit(&[], &[], &[], 0);
    for n in sizes {
        let x: Vec<_> = (0..n).map(|i| i as f32 - 2.0).collect();
        let y: Vec<_> = x.iter().map(|v| -3.0 + 2.0 * v).collect();
        check_fit(&x, &y, &[], 0);
        check_fit(&x, &vec![0.0; n], &[], 0);
        check_fit(&vec![0.0; n], &y, &[], 0);
        check_fit(&vec![-2.0; n], &y, &[], 0);
        if n >= 3 {
            check_fit(&x, &y, &vec![1.0; n], 1);
            check_fit(&x, &vec![0.0; n], &vec![-1.0; n], -7);
        }
    }
    for value in [
        f32::from_bits(1),
        f32::MIN_POSITIVE,
        1e-30,
        1e30,
        f32::MAX,
        f32::NAN,
        f32::INFINITY,
        f32::NEG_INFINITY,
    ] {
        check_fit(&[-value, 0.0, value], &[1.0, -2.0, 3.0], &[], 0);
        check_fit(&[-1.0, 0.0, 1.0], &[-value, 0.0, value], &[], 0);
    }
}

// Run the same test executable in isolation so C exit() cannot kill the suite.
// This test also acts as its own subprocess entry point; there are no ignored
// tests or filtered-out input cases in the random differential tests above.
#[test]
fn diagnostics_match_c() {
    const CASE_ENV: &str = "OXWAVE_LINEAR_FIT_DIAGNOSTIC_CASE";
    const SIDE_ENV: &str = "OXWAVE_LINEAR_FIT_DIAGNOSTIC_SIDE";
    if let Ok(case) = std::env::var(CASE_ENV) {
        let c = std::env::var(SIDE_ENV).unwrap() == "c";
        let mut output = -123.0;
        let mut gln = -456.0;
        match case.as_str() {
            "gammq_x" | "gammq_a_zero" | "gammq_a_negative" => {
                let (a, x) = match case.as_str() {
                    "gammq_x" => (1.0, -1.0),
                    "gammq_a_zero" => (0.0, 1.0),
                    _ => (-1.0, 0.0),
                };
                if c {
                    unsafe { gammq(a, x) };
                } else {
                    rust::gammq(a, x);
                }
            }
            "gser_negative" => {
                if c {
                    unsafe { gser(&mut output, 1.0, -1.0, &mut gln) };
                } else {
                    rust::gser(&mut output, 1.0, -1.0, &mut gln);
                }
            }
            "gcf_singular" | "gcf_tiny" | "gcf_nan" | "gcf_limit" => {
                let (a, x) = match case.as_str() {
                    "gcf_singular" => (1.0, 0.0),
                    "gcf_tiny" => (1.0, 1e-35),
                    "gcf_nan" => (f32::NAN, 1.0),
                    _ => (1e6, 1e6),
                };
                if c {
                    unsafe { gcf(&mut output, a, x, &mut gln) };
                } else {
                    rust::gcf(&mut output, a, x, &mut gln);
                }
            }
            "gser_limit" | "gser_zero_a" => {
                let (a, x) = if case == "gser_limit" {
                    (1e4, 1e4)
                } else {
                    (0.0, 1.0)
                };
                if c {
                    unsafe { gser(&mut output, a, x, &mut gln) };
                } else {
                    rust::gser(&mut output, a, x, &mut gln);
                }
                assert_eq!(output.to_bits(), (-123.0_f32).to_bits());
                same(gln, unsafe { gammln(a) }, "gser failure gln");
                return;
            }
            "linfit_one" | "linfit_two" | "linfit_zero_sig" => {
                let n = match case.as_str() {
                    "linfit_one" => 1,
                    "linfit_two" => 2,
                    _ => 3,
                };
                let x = [0.0, 1.0, 2.0];
                let y = [1.0, 2.0, 3.0];
                let sig = if case == "linfit_zero_sig" {
                    [0.0; 3]
                } else {
                    [1.0; 3]
                };
                if c {
                    fit_c(&x[..n], &y[..n], &sig[..n], 1);
                } else {
                    fit_rust(&x[..n], &y[..n], &sig[..n], 1);
                }
            }
            _ => panic!("unknown diagnostic case {case}"),
        }
        panic!("{case} should have exited");
    }

    let cases = [
        ("gammq_x", 1, "Invalid arguments in routine gammq\n"),
        ("gammq_a_zero", 1, "Invalid arguments in routine gammq\n"),
        (
            "gammq_a_negative",
            1,
            "Invalid arguments in routine gammq\n",
        ),
        ("gser_negative", 1, "x less than 0 in routine gser\n"),
        ("gcf_singular", 1, "a too large, ITMAX too small in gcf\n"),
        ("gcf_tiny", 1, "a too large, ITMAX too small in gcf\n"),
        ("gcf_nan", 1, "a too large, ITMAX too small in gcf\n"),
        ("gcf_limit", 1, "a too large, ITMAX too small in gcf\n"),
        (
            "gser_limit",
            0,
            "a too large, ITMAX too small in routine gser\n",
        ),
        (
            "gser_zero_a",
            0,
            "a too large, ITMAX too small in routine gser\n",
        ),
        ("linfit_one", 1, "Invalid arguments in routine gammq\n"),
        ("linfit_two", 1, "Invalid arguments in routine gammq\n"),
        (
            "linfit_zero_sig",
            1,
            "a too large, ITMAX too small in gcf\n",
        ),
    ];
    for (case, code, message) in cases {
        for side in ["c", "rust"] {
            let output = Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "diagnostics_match_c", "--nocapture"])
                .env(CASE_ENV, case)
                .env(SIDE_ENV, side)
                .output()
                .unwrap();
            assert_eq!(
                output.status.code(),
                Some(code),
                "{side} {case}: {output:?}"
            );
            assert_eq!(output.stderr, message.as_bytes(), "{side} {case}");
        }
    }
}
