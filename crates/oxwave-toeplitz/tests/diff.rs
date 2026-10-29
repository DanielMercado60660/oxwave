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

//! Differential tests for cwp/lib/stoep.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

use oxwave_cref::XorShift;

extern "C" {
    fn stoepd(n: i32, r: *mut f64, g: *mut f64, f: *mut f64, a: *mut f64);
    fn stoepf(n: i32, r: *mut f32, g: *mut f32, f: *mut f32, a: *mut f32);
}

#[test]
fn stoepd_matches_c() {
    let mut rng = XorShift::new(0x1234_5678);
    let sizes = [
        1, 2, 3, 5, 7, 11, 16, 31, 32, 61, 64, 127, 128, 257, 512, 1024,
    ];
    let scales: [f64; 5] = [1e-100, 1e-8, 1.0, 1e8, 1e100];
    for case in 0..6000 {
        let n = sizes[case % sizes.len()];
        let block = case / sizes.len();
        let scale = scales[(block / 12) % scales.len()];
        let rhs_scale = scales[(block / 12 + 2) % scales.len()];
        let mut r: Vec<f64> = (0..n)
            .map(|_| (rng.uniform(-1.0, 1.0) as f64) * scale / n as f64)
            .collect();
        let mut g: Vec<f64> = (0..n)
            .map(|_| (rng.uniform(-1.0, 1.0) as f64) * rhs_scale)
            .collect();
        r[0] = 2.0 * scale;
        match block % 12 {
            0 => r.fill(0.0),
            1 => r[0] = 0.0,
            2 => r[0] = -0.0,
            3 => {
                r.fill(0.0);
                r[0] = scale;
            }
            4 => r[0] = -2.0 * scale,
            5 | 6 => {
                // Autocorrelations with either sign, including nearly singular systems.
                let rho = if block % 12 == 5 {
                    rng.uniform(0.0, 0.9999) as f64
                } else {
                    rng.uniform(-0.9999, 0.0) as f64
                };
                r[0] = scale;
                for i in 1..n {
                    r[i] = r[i - 1] * rho;
                }
            }
            7 => r.fill(scale), // Singular for n > 1.
            8 => {
                // Singular leading 2x2 system, even if the full system is invertible.
                r.fill(0.0);
                r[0] = scale;
                if n > 1 {
                    r[1] = scale;
                }
            }
            9 => g.fill(0.0),
            10 => {
                // General indefinite systems, with positive and negative lags.
                for x in &mut r {
                    *x = (rng.uniform(-1.0, 1.0) as f64) * scale;
                }
            }
            _ => {} // Strictly diagonally dominant systems.
        }
        compare_stoepd(&r, &g, &format!("random case {case}, n={n}"));
    }

    // No cases are filtered out when the recurrence becomes singular.
    for r in [
        vec![0.0],
        vec![-0.0],
        vec![f64::from_bits(1)],
        vec![-f64::from_bits(1)],
        vec![f64::MIN_POSITIVE, f64::from_bits(1)],
        vec![f64::MAX, 0.0],
        vec![f64::INFINITY, 0.0],
        vec![f64::NAN, 1.0],
        vec![1.0, f64::INFINITY],
        vec![1.0, f64::NAN],
    ] {
        for rhs in [0.0, -0.0, 1.0, -1.0, f64::MAX, f64::INFINITY, f64::NAN] {
            let g = vec![rhs; r.len()];
            compare_stoepd(&r, &g, &format!("special r={r:?}, rhs={rhs:?}"));
        }
    }
}

fn compare_stoepd(r: &[f64], g: &[f64], context: &str) {
    let n = r.len();
    let mut cr = r.to_vec();
    let mut cg = g.to_vec();
    let mut cf = vec![17.0; n];
    let mut ca = vec![-13.0; n];
    // Check output tails independently; padding never contributes to a tolerance.
    let padding = [-0.0, f64::NAN, f64::INFINITY, -f64::INFINITY];
    cf.extend_from_slice(&padding);
    ca.extend_from_slice(&padding);
    let mut rf = cf.clone();
    let mut ra = ca.clone();
    unsafe {
        stoepd(
            n as i32,
            cr.as_mut_ptr(),
            cg.as_mut_ptr(),
            cf.as_mut_ptr(),
            ca.as_mut_ptr(),
        );
    }
    oxwave_toeplitz::stoepd(r, g, &mut rf, &mut ra);
    for (label, actual, expected) in [("f", &rf, &cf), ("a", &ra, &ca)] {
        for (i, (&x, &y)) in actual[..n].iter().zip(&expected[..n]).enumerate() {
            // Same operations and precision permit a stricter check than 1e-5.
            // NaN payloads are not part of the numerical contract.
            assert!(
                x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()),
                "{context}, {label}[{i}]: Rust {x:?}, C {y:?}"
            );
            if r[0] == 0.0 {
                let initial: f64 = if label == "f" { 17.0 } else { -13.0 };
                assert_eq!(x.to_bits(), initial.to_bits(), "{context}: output changed");
            }
        }
        for (i, &sentinel) in padding.iter().enumerate() {
            assert_eq!(
                actual[n + i].to_bits(),
                sentinel.to_bits(),
                "{context}: Rust tail"
            );
            assert_eq!(
                expected[n + i].to_bits(),
                sentinel.to_bits(),
                "{context}: C tail"
            );
        }
    }
    // C's non-const pointer prototypes do not imply mutation of its inputs.
    for (label, actual, expected) in [("r", &cr, r), ("g", &cg, g)] {
        for (&x, &y) in actual.iter().zip(expected) {
            assert_eq!(x.to_bits(), y.to_bits(), "{context}: C changed {label}");
        }
    }
}

#[test]
fn stoepf_matches_c() {
    let mut rng = XorShift::new(0x8765_4321);
    let sizes = [
        1, 2, 3, 5, 7, 11, 16, 31, 32, 61, 64, 127, 128, 257, 512, 1024,
    ];
    let scales: [f32; 5] = [1e-12, 1e-6, 1.0, 1e6, 1e12];
    for case in 0..6000 {
        let n = sizes[case % sizes.len()];
        let block = case / sizes.len();
        let scale = scales[(block / 12) % scales.len()];
        let rhs_scale = scales[(block / 12 + 2) % scales.len()];
        let mut r: Vec<f32> = (0..n)
            .map(|_| (rng.uniform(-1.0, 1.0) as f32) * scale / n as f32)
            .collect();
        let mut g: Vec<f32> = (0..n)
            .map(|_| (rng.uniform(-1.0, 1.0) as f32) * rhs_scale)
            .collect();
        r[0] = 2.0 * scale;
        match block % 12 {
            0 => r.fill(0.0),
            1 => r[0] = 0.0,
            2 => r[0] = -0.0,
            3 => {
                r.fill(0.0);
                r[0] = scale;
            }
            4 => r[0] = -2.0 * scale,
            5 | 6 => {
                // Autocorrelations with either sign, including nearly singular systems.
                let rho = if block % 12 == 5 {
                    rng.uniform(0.0, 0.9999) as f32
                } else {
                    rng.uniform(-0.9999, 0.0) as f32
                };
                r[0] = scale;
                for i in 1..n {
                    r[i] = r[i - 1] * rho;
                }
            }
            7 => r.fill(scale), // Singular for n > 1.
            8 => {
                // Singular leading 2x2 system, even if the full system is invertible.
                r.fill(0.0);
                r[0] = scale;
                if n > 1 {
                    r[1] = scale;
                }
            }
            9 => g.fill(0.0),
            10 => {
                // General indefinite systems, with positive and negative lags.
                for x in &mut r {
                    *x = (rng.uniform(-1.0, 1.0) as f32) * scale;
                }
            }
            _ => {} // Strictly diagonally dominant systems.
        }
        compare_stoepf(&r, &g, &format!("random case {case}, n={n}"));
    }

    // No cases are filtered out when the recurrence becomes singular.
    for r in [
        vec![0.0],
        vec![-0.0],
        vec![f32::from_bits(1)],
        vec![-f32::from_bits(1)],
        vec![f32::MIN_POSITIVE, f32::from_bits(1)],
        vec![f32::MAX, 0.0],
        vec![f32::INFINITY, 0.0],
        vec![f32::NAN, 1.0],
        vec![1.0, f32::INFINITY],
        vec![1.0, f32::NAN],
    ] {
        for rhs in [0.0, -0.0, 1.0, -1.0, f32::MAX, f32::INFINITY, f32::NAN] {
            let g = vec![rhs; r.len()];
            compare_stoepf(&r, &g, &format!("special r={r:?}, rhs={rhs:?}"));
        }
    }
}

fn compare_stoepf(r: &[f32], g: &[f32], context: &str) {
    let n = r.len();
    let mut cr = r.to_vec();
    let mut cg = g.to_vec();
    let mut cf = vec![17.0; n];
    let mut ca = vec![-13.0; n];
    // Check output tails independently; padding never contributes to a tolerance.
    let padding = [-0.0, f32::NAN, f32::INFINITY, -f32::INFINITY];
    cf.extend_from_slice(&padding);
    ca.extend_from_slice(&padding);
    let mut rf = cf.clone();
    let mut ra = ca.clone();
    unsafe {
        stoepf(
            n as i32,
            cr.as_mut_ptr(),
            cg.as_mut_ptr(),
            cf.as_mut_ptr(),
            ca.as_mut_ptr(),
        );
    }
    oxwave_toeplitz::stoepf(r, g, &mut rf, &mut ra);
    for (label, actual, expected) in [("f", &rf, &cf), ("a", &ra, &ca)] {
        for (i, (&x, &y)) in actual[..n].iter().zip(&expected[..n]).enumerate() {
            // Same operations and precision permit a stricter check than 1e-5.
            // NaN payloads are not part of the numerical contract.
            assert!(
                x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()),
                "{context}, {label}[{i}]: Rust {x:?}, C {y:?}"
            );
            if r[0] == 0.0 {
                let initial: f32 = if label == "f" { 17.0 } else { -13.0 };
                assert_eq!(x.to_bits(), initial.to_bits(), "{context}: output changed");
            }
        }
        for (i, &sentinel) in padding.iter().enumerate() {
            assert_eq!(
                actual[n + i].to_bits(),
                sentinel.to_bits(),
                "{context}: Rust tail"
            );
            assert_eq!(
                expected[n + i].to_bits(),
                sentinel.to_bits(),
                "{context}: C tail"
            );
        }
    }
    // C's non-const pointer prototypes do not imply mutation of its inputs.
    for (label, actual, expected) in [("r", &cr, r), ("g", &cg, g)] {
        for (&x, &y) in actual.iter().zip(expected) {
            assert_eq!(x.to_bits(), y.to_bits(), "{context}: C changed {label}");
        }
    }
}
