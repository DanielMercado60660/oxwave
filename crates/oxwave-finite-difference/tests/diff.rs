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

//! Differential tests for cwp/lib/holbergd1.c, mkdiff.c, and mkhdiff.c.
//! Author of the original routines: Dave Hale, Colorado School of Mines.
use oxwave_cref::XorShift;
use oxwave_finite_difference as rust;
extern "C" {
    fn holbergd1(e: f32, n: i32, d: *mut f32);
    fn mkdiff(n: i32, a: f32, h: f32, l: i32, m: i32, d: *mut f32);
    fn mkhdiff(h: f32, l: i32, d: *mut f32);
}
fn compare(a: &[f32], b: &[f32], context: &str) {
    let scale = b
        .iter()
        .filter(|x| x.is_finite())
        .fold(0.0_f32, |m, x| m.max(x.abs()));
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        assert!(
            (x.is_nan() && y.is_nan())
                || x == y
                || (x.is_finite() && y.is_finite() && (x - y).abs() <= 1e-5 * scale),
            "{context}: index {i}, Rust {x:?}, C {y:?}, scale {scale}"
        );
    }
}
// The C loop reads beyond a four-row coefficient subarray for n >= 6.
// Evaluate only its retained coefficients for those filter lengths.
fn holberg_reference(e: f32, n: usize) -> Vec<f32> {
    let rows: &[[f32; 6]] = if n == 6 {
        &[
            [1.172, 0.2742, -0.3006, 0.2637, -0.2391, 0.0],
            [-0.06510, -0.1371, -0.0100, 0.1077, -0.0086, 0.0],
            [0.004688, 0.0274, 0.0661, 0.0826, 0.0530, 0.0],
        ]
    } else {
        &[
            [1.1963, 0.1987, -0.2195, 0.2096, -0.2038, 0.1779],
            [-0.07975, -0.1192, 0.03929, 0.04049, -0.05841, 0.04389],
            [0.009570, 0.03975, 0.04853, -0.00967, -0.05507, -0.01117],
            [-0.0006975, -0.00568, -0.02014, -0.04039, -0.04937, -0.03625],
        ]
    };
    let mut result = vec![0.0; n];
    for i in 0..n / 2 {
        let mut coefficient = 0.0_f32;
        for j in 0..6 {
            coefficient = (coefficient as f64
                + rows[i][j] as f64 * (e as f64).powf(2.0 * j as f64 / n as f64))
                as f32;
        }
        result[n / 2 + i] = -coefficient;
    }
    for i in 0..n / 2 {
        result[i] = -result[n - i - 1];
    }
    result
}

#[test]
fn holberg_differential() {
    let mut rng = XorShift::new(123);
    for n in [2, 4, 6, 8] {
        for k in 0..2500 {
            let e = match k % 10 {
                0 => 0.0,
                1 => -rng.uniform_f32(0.001, 2.0),
                2 => 1.0,
                _ => rng.uniform_f32(0.000001, 2.0),
            };
            let mut a = vec![0.0; n];
            let mut b = a.clone();
            rust::holbergd1(e, &mut a);
            if n <= 4 {
                unsafe { holbergd1(e, n as i32, b.as_mut_ptr()) };
            } else {
                b = holberg_reference(e, n);
            }
            compare(&a, &b, &format!("n={n} e={e}"));
        }
    }
}
#[test]
fn mkdiff_differential() {
    let mut rng = XorShift::new(456);
    for size in [1, 2, 3, 5, 8, 13, 16, 31, 64, 100] {
        for k in 0..400 {
            let n = (rng.next_u32() as usize % size) as i32;
            let l = (rng.next_u32() % 21) as i32 - 10;
            let m = l + size as i32 - 1;
            let a = if k % 7 == 0 {
                0.0
            } else {
                rng.uniform_f32(-1.0, 1.0)
            };
            let h = match k % 9 {
                0 => 0.0,
                1 => -1.0,
                _ => rng.uniform_f32(0.2, 3.0),
            };
            let mut x = vec![0.0; size];
            let mut y = x.clone();
            rust::mkdiff(n, a, h, l, m, &mut x);
            unsafe { mkdiff(n, a, h, l, m, y.as_mut_ptr()) };
            compare(&x, &y, &format!("size={size} n={n} a={a} h={h} l={l}"));
        }
    }
}
#[test]
fn mkhdiff_differential() {
    let mut rng = XorShift::new(789);
    for l in [0, 1, 2, 3, 5, 8, 13, 32, 64] {
        for k in 0..300 {
            let h = match k % 10 {
                0 => 0.0,
                1 => -1.0,
                _ => rng.uniform_f32(0.001, 10.0),
            };
            let mut x = vec![0.0; 1 + 2 * l as usize];
            let mut y = x.clone();
            rust::mkhdiff(h, l, &mut x);
            unsafe { mkhdiff(h, l, y.as_mut_ptr()) };
            compare(&x, &y, &format!("l={l} h={h}"));
        }
    }
}
