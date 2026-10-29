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

//! Differential tests for cwp/lib/linregress.c.
//! Author: Chris Liner, UTulsa, 11/16/03.

use oxwave_cref::XorShift;

extern "C" {
    /// C prototype: `void linear_regression(float *y, float *x, int n, float coeff[4]);`
    fn linear_regression(y: *mut f32, x: *mut f32, n: i32, coeff: *mut f32);
}

fn compare(y: &[f32], x: &[f32]) -> [f32; 4] {
    assert_eq!(y.len(), x.len());
    assert!(!x.is_empty());
    let mut cy = y.to_vec();
    let mut cx = x.to_vec();
    let mut c = [f32::NAN; 4];
    let mut rust = [f32::NAN; 4];
    // The C routine accesses n initialized elements and writes four coefficients.
    unsafe {
        linear_regression(
            cy.as_mut_ptr(),
            cx.as_mut_ptr(),
            x.len() as i32,
            c.as_mut_ptr(),
        )
    };
    oxwave_linear_regression::linear_regression(y, x, &mut rust);
    for i in 0..4 {
        let (a, b) = (rust[i], c[i]);
        if b.is_nan() {
            // NaN sign and payload are not part of the numerical result.
            assert!(
                a.is_nan(),
                "n={} coefficient={i} Rust={a:?} C={b:?}",
                x.len()
            );
        } else {
            // Stronger than 1e-5 relative: also checks signed zero and infinity.
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "n={} coefficient={i} Rust={a:?} C={b:?}",
                x.len()
            );
        }
    }
    rust
}

#[test]
fn linear_regression_matches_c() {
    let mut rng = XorShift::new(0x42a32);
    let sizes = [1, 2, 3, 5, 7, 16, 31, 64, 127, 256, 1024, 4093, 8192];
    let scales = [1e-12, 0.0001, 1.0, 1000.0, 1e10];
    for case in 0..7800 {
        let n = sizes[case % sizes.len()];
        let group = case / sizes.len();
        let xscale = scales[group % scales.len()];
        let yscale = scales[(group / scales.len()) % scales.len()];
        let mut x = rng.vec_f32(n, -xscale, xscale);
        let mut y = rng.vec_f32(n, -yscale, yscale);
        match (group / (scales.len() * scales.len())) % 12 {
            0 => x.fill(0.0),
            1 => y.fill(0.0),
            2 => {
                x.fill(-3.0);
                y.fill(5.0);
            }
            3 => {
                for i in 0..n {
                    y[i] = x[i];
                }
            }
            4 => {
                for i in 0..n {
                    y[i] = -x[i];
                }
            }
            5 => {
                for i in 0..n {
                    if i % 3 == 0 {
                        x[i] = 0.0;
                        y[i] = 0.0;
                    }
                }
            }
            6 | 7 => {
                // Noisy trends with either sign of slope and a nonzero intercept.
                let slope = rng.uniform_f32(-4.0, 4.0);
                let intercept = rng.uniform_f32(-10.0, 10.0) * xscale;
                for i in 0..n {
                    y[i] = slope * x[i] + intercept + 0.01 * y[i];
                }
            }
            8 => {
                // Small variation about a large mean causes cancellation in C.
                for i in 0..n {
                    x[i] = xscale + x[i] * 1e-5;
                    y[i] = -yscale + y[i] * 1e-5;
                }
            }
            9 => {
                for i in 0..n {
                    x[i] = -x[i].abs();
                    y[i] = -y[i].abs();
                }
            }
            _ => (),
        }
        compare(&y, &x);
    }
}

#[test]
fn linear_regression_degenerate_and_extreme_matches_c() {
    assert_eq!(compare(&[0.0], &[0.0])[2], 999.0);
    let two = compare(&[1.0, -1.0], &[-1.0, 1.0]);
    assert_eq!(&two[..3], &[-1.0, 0.0, 1.0]);
    assert!(two[3].is_nan()); // C computes 0/0 for the standard error at n=2.
    assert_eq!(
        compare(&[7.0, 5.0, 3.0, 1.0], &[-1.0, 0.0, 1.0, 2.0]),
        [-2.0, 5.0, 1.0, 0.0]
    );
    assert_eq!(
        compare(&[4.0, 4.0, 4.0], &[-2.0, 0.0, 2.0]),
        [0.0, 4.0, 999.0, 0.0]
    );
    let values = [
        0.0,
        -0.0,
        f32::from_bits(1),
        -f32::from_bits(1),
        f32::MIN_POSITIVE,
        -f32::MIN_POSITIVE,
        1e-20,
        -1e-20,
        1.0,
        -1.0,
        1e20,
        -1e20,
        f32::MAX,
        -f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ];
    for n in [1, 2, 3, 17, 128, 4093, 65536] {
        for &v in &values {
            // Test both constant and alternating signs, including underflow,
            // overflow, and nonfinite inputs; every case is sent to C.
            let x = vec![v; n];
            let y: Vec<f32> = (0..n).map(|i| if i % 2 == 0 { v } else { -v }).collect();
            compare(&y, &x);
            compare(&y, &y);
        }
    }
}

#[test]
fn linear_regression_random_float_bits_matches_c() {
    let mut rng = XorShift::new(0xbeef_0123_4567_89ab);
    for case in 0..3000 {
        let n = [1, 2, 3, 7, 32, 127, 257][case % 7];
        let x: Vec<f32> = (0..n).map(|_| f32::from_bits(rng.next_u32())).collect();
        let y: Vec<f32> = (0..n).map(|_| f32::from_bits(rng.next_u32())).collect();
        compare(&y, &x);
    }
}
