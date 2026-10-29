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

//! Differential tests for cwp/lib/unwrapphase.c.
//! Author: John Stockwell, CWP, 2010.

use oxwave_cref::XorShift;
use std::process::Command;

extern "C" {
    fn oppenheim_unwrap_phase(
        n: i32,
        trend: i32,
        zeromean: i32,
        df: f32,
        xr: *mut f32,
        xi: *mut f32,
        phase: *mut f32,
    );
    fn simple_unwrap_phase(n: i32, trend: i32, zeromean: i32, w: f32, phase: *mut f32);
}

const SIZES: &[usize] = &[
    1, 2, 3, 5, 7, 13, 16, 31, 32, 61, 64, 97, 127, 128, 251, 256, 509, 512, 1021, 1024, 2048,
    4093, 4096, 8192,
];

fn compare(actual: &[f32], expected: &[f32], context: &str) {
    assert_eq!(actual.len(), expected.len());
    for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
        // NaN payloads/signs are not specified by C. All other results,
        // including infinities and signed zeros, must be bit-exact: this
        // is stronger than 1e-5 relative to the maximum output magnitude.
        assert!(
            (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits(),
            "{context}, index {i}: Rust {a:?} ({:08x}), C {b:?} ({:08x})",
            a.to_bits(),
            b.to_bits(),
        );
    }
}

fn check_oppenheim(xr: &[f32], xi: &[f32], df: f32, trend: i32, zeromean: i32) {
    let n = xr.len();
    // For n=1, C's differentiate reads both f[-1] and f[1]. Place the
    // logical signal between initialized guard samples so these reads are
    // defined. Neither read affects phase[0], which is never integrated.
    let mut cr = vec![13.0f32; n + 2];
    let mut ci = vec![-7.0f32; n + 2];
    cr[1..n + 1].copy_from_slice(xr);
    ci[1..n + 1].copy_from_slice(xi);
    let mut expected = vec![12345.0f32; n + 2];
    let mut actual = vec![-12345.0f32; n];
    unsafe {
        oppenheim_unwrap_phase(
            n as i32,
            trend,
            zeromean,
            df,
            cr.as_mut_ptr().add(1),
            ci.as_mut_ptr().add(1),
            expected.as_mut_ptr().add(1),
        );
    }
    oxwave_phase_unwrap::oppenheim_unwrap_phase(trend, zeromean, df, xr, xi, &mut actual);
    let context = format!("oppenheim n={n}, df={df}, trend={trend}, zeromean={zeromean}");
    compare(&actual, &expected[1..n + 1], &context);
    assert_eq!(expected[0], 12345.0);
    assert_eq!(expected[n + 1], 12345.0);
    compare(&cr[1..n + 1], xr, "C real input unchanged");
    compare(&ci[1..n + 1], xi, "C imaginary input unchanged");
}

fn check_simple(input: &[f32], w: f32, trend: i32, zeromean: i32) {
    let n = input.len();
    let mut expected = vec![12345.0f32; n + 2];
    expected[1..n + 1].copy_from_slice(input);
    let mut actual = input.to_vec();
    unsafe {
        simple_unwrap_phase(n as i32, trend, zeromean, w, expected.as_mut_ptr().add(1));
    }
    oxwave_phase_unwrap::simple_unwrap_phase(trend, zeromean, w, &mut actual);
    let context = format!("simple n={n}, w={w}, trend={trend}, zeromean={zeromean}");
    compare(&actual, &expected[1..n + 1], &context);
    assert_eq!(expected[0], 12345.0);
    assert_eq!(expected[n + 1], 12345.0);
}

#[test]
fn oppenheim_unwrap_phase_matches_c() {
    let mut rng = XorShift::new(0x0bbe_1234_5678);
    let scales = [0.0, 1e-20, 1e-10, 1e-3, 1.0, 1e3, 1e10, 1e18];
    let intervals = [0.0, -0.0, 1e-20, 1e-5, 0.004, 0.1, 1.0, 10.0, -0.125, 1e20];
    for case in 0..4096 {
        let n = SIZES[case % SIZES.len()];
        let scale = scales[(case / SIZES.len()) % scales.len()];
        let df = if case % 3 == 0 {
            rng.uniform_f32(-10.0, 10.0)
        } else {
            intervals[(case / 3) % intervals.len()]
        };
        let trend = rng.range_i32(-1, 1);
        let zeromean = rng.range_i32(-1, 1);
        let mut xr = rng.vec_f32(n, -scale, scale);
        let mut xi = rng.vec_f32(n, -scale, scale);
        match case % 7 {
            0 => {
                xr.fill(0.0);
                xi.fill(0.0);
            }
            1 => {
                xr.fill(-scale);
                xi.fill(scale);
            }
            2 => {
                // Smooth complex spectra with a frequency-dependent phase.
                let slope = rng.uniform(-3.0, 3.0);
                let offset = rng.uniform(-10.0, 10.0);
                for i in 0..n {
                    let angle = slope * i as f64 + offset;
                    xr[i] = scale * angle.cos() as f32;
                    xi[i] = scale * angle.sin() as f32;
                }
            }
            3 => {
                xr.fill(0.0);
                xi[n / 2] = 0.0;
            }
            4 => {
                // Isolated spectral nulls exercise the exact-zero denominator.
                for i in (0..n).step_by(3) {
                    xr[i] = 0.0;
                    xi[i] = -0.0;
                }
            }
            _ => {}
        }
        check_oppenheim(&xr, &xi, df, trend, zeromean);
    }

    // Nonfinite intermediates, underflow, overflow, signed zeros, and n=1.
    let edges = [
        0.0,
        -0.0,
        f32::from_bits(1),
        f32::MIN_POSITIVE,
        1.0,
        -1.0,
        f32::MAX,
        -f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ];
    for &n in &[1, 2, 3, 17, 64] {
        for &value in &edges {
            let xr = vec![value; n];
            let xi: Vec<_> = (0..n).map(|i| edges[i % edges.len()]).collect();
            for trend in [0, 1, -7] {
                for zeromean in [0, 1] {
                    for df in [0.0, -0.0, 0.25, -2.0, f32::INFINITY, f32::NAN] {
                        check_oppenheim(&xr, &xi, df, trend, zeromean);
                    }
                }
            }
        }
    }
    for case in 0..256 {
        let n = SIZES[case % SIZES.len()];
        let xr: Vec<_> = (0..n).map(|_| rng.any_f32()).collect();
        let xi: Vec<_> = (0..n).map(|_| rng.any_f32()).collect();
        let df = rng.any_f32();
        check_oppenheim(&xr, &xi, df, (case % 2) as i32, ((case / 2) % 2) as i32);
    }
}

#[test]
fn simple_unwrap_phase_matches_c() {
    let mut rng = XorShift::new(0x51ab_1234_5678);
    let scales = [
        0.0,
        1e-30,
        1e-10,
        1e-3,
        std::f32::consts::PI,
        1e3,
        1e18,
        1e37,
    ];
    let widths = [
        0.1,
        1.0,
        2.0,
        -1.0,
        1e-30,
        1e30,
        f32::MAX,
        f32::MIN_POSITIVE,
        f32::from_bits(1),
    ];
    for case in 0..4096 {
        let n = SIZES[case % SIZES.len()];
        let scale = scales[(case / SIZES.len()) % scales.len()];
        let w = if case % 3 == 0 {
            let magnitude = rng.uniform_f32(0.001, 20.0);
            if case % 2 == 0 {
                magnitude
            } else {
                -magnitude
            }
        } else {
            widths[(case / 3) % widths.len()]
        };
        let trend = rng.range_i32(-1, 1);
        let zeromean = rng.range_i32(-1, 1);
        let mut input = rng.vec_f32(n, -scale, scale);
        match case % 6 {
            0 => input.fill(0.0),
            1 => input.fill(-scale),
            2 => {
                let slope = rng.uniform(-2.0, 2.0);
                let offset = rng.uniform(-10.0, 10.0);
                for i in 0..n {
                    let angle = slope * i as f64 + offset;
                    input[i] = (angle.sin().atan2(angle.cos())) as f32;
                }
            }
            3 => {
                input.fill(0.0);
                input[n / 2] = scale;
            }
            _ => {}
        }
        check_simple(&input, w, trend, zeromean);
    }

    // At, just below, and just above the >= PI/w decision boundary.
    for w in [0.25f32, 1.0, 2.0, 100.0] {
        let threshold = (std::f64::consts::PI / w as f64) as f32;
        for jump in [
            f32::from_bits(threshold.to_bits() - 1),
            threshold,
            f32::from_bits(threshold.to_bits() + 1),
        ] {
            for sign in [-1.0, 1.0] {
                for n in [2, 3, 7, 32] {
                    let mut input = vec![0.0f32; n];
                    input[1] = sign * jump;
                    for trend in [0, 1, -7] {
                        for zeromean in [0, 1] {
                            check_simple(&input, w, trend, zeromean);
                        }
                    }
                }
            }
        }
    }
    let edges = [
        0.0,
        -0.0,
        f32::from_bits(1),
        f32::MIN_POSITIVE,
        1.0,
        -1.0,
        f32::MAX,
        -f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ];
    for &n in &[1, 2, 3, 17, 64] {
        for &value in &edges {
            let mut input: Vec<_> = (0..n).map(|i| edges[i % edges.len()]).collect();
            input[0] = value;
            for trend in [0, 1, -7] {
                for zeromean in [0, 1] {
                    for w in [1.0, -1.0, f32::INFINITY, f32::NEG_INFINITY, f32::NAN] {
                        check_simple(&input, w, trend, zeromean);
                    }
                }
            }
        }
    }
    for case in 0..256 {
        let n = SIZES[case % SIZES.len()];
        let input: Vec<_> = (0..n).map(|_| rng.any_f32()).collect();
        check_simple(&input, 1.0, (case % 2) as i32, ((case / 2) % 2) as i32);
    }
}

#[test]
fn simple_zero_w_matches_c() {
    // Isolate C's exit(1) and the Rust equivalent in child test processes.
    if let Ok(backend) = std::env::var("OXWAVE_PHASE_ZERO_CHILD") {
        let w = if std::env::var("OXWAVE_PHASE_NEGATIVE_ZERO").unwrap() == "yes" {
            -0.0
        } else {
            0.0
        };
        let mut phase = [1.0f32, -1.0];
        if backend == "c" {
            unsafe {
                simple_unwrap_phase(2, 0, 0, w, phase.as_mut_ptr());
            }
        } else {
            oxwave_phase_unwrap::simple_unwrap_phase(0, 0, w, &mut phase);
        }
        panic!("zero wrapping parameter must exit");
    }
    for negative_zero in ["no", "yes"] {
        let run = |backend| {
            Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "simple_zero_w_matches_c", "--nocapture"])
                .env("OXWAVE_PHASE_ZERO_CHILD", backend)
                .env("OXWAVE_PHASE_NEGATIVE_ZERO", negative_zero)
                .output()
                .unwrap()
        };
        let c = run("c");
        let rust = run("rust");
        assert_eq!(c.status.code(), Some(1));
        assert_eq!(rust.status.code(), c.status.code());
        assert_eq!(c.stderr, b"wrapping parameter is zero");
        assert_eq!(rust.stderr, c.stderr);
        assert_eq!(rust.stdout, c.stdout);
    }
}
