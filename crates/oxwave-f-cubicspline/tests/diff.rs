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

//! Differential tests for `cwp/lib/cubicspline.c`.
//! Author: Dave Hale, Colorado School of Mines c. 1989, 1990, 1991.
use oxwave_cref::XorShift;
use oxwave_f_cubicspline as rust;

extern "C" {
    fn cakima(n: i32, x: *mut f32, y: *mut f32, yd: *mut [f32; 4]);
    fn chermite(n: i32, x: *mut f32, y: *mut f32, yd: *mut [f32; 4]);
    fn cmonot(n: i32, x: *mut f32, y: *mut f32, yd: *mut [f32; 4]);
    fn csplin(n: i32, x: *mut f32, y: *mut f32, yd: *mut [f32; 4]);
}

type CFn = unsafe extern "C" fn(i32, *mut f32, *mut f32, *mut [f32; 4]);
type RustFn = fn(&[f32], &[f32], &mut [[f32; 4]]);

fn check_case(c: CFn, f: RustFn, hermite: bool, x: &[f32], y: &[f32], case: usize) {
    let n = x.len();
    assert_eq!(y.len(), n);
    assert!(x.windows(2).all(|w| w[0] < w[1]) || x.windows(2).all(|w| w[0] > w[1]));
    let mut cx = x.to_vec();
    let mut cy = y.to_vec();
    // Extra ordinates exercise the slice API's prefix semantics.
    let mut ry = y.to_vec();
    ry.extend_from_slice(&[f32::NAN, f32::INFINITY]);
    let sentinel = [17.0, -19.0, f32::from_bits(0x7fc1_2345), -0.0];
    let mut expected = vec![sentinel; n + 2];
    let mut actual = expected.clone();
    unsafe {
        c(
            n as i32,
            cx.as_mut_ptr(),
            cy.as_mut_ptr(),
            expected.as_mut_ptr(),
        )
    };
    f(x, &ry, &mut actual);
    for i in 0..n {
        assert_eq!(cx[i].to_bits(), x[i].to_bits(), "C changed x");
        assert_eq!(cy[i].to_bits(), y[i].to_bits(), "C changed y");
        assert_eq!(ry[i].to_bits(), y[i].to_bits(), "Rust changed y");
    }
    for i in 0..n + 2 {
        for j in 0..4 {
            let a = actual[i][j];
            let b = expected[i][j];
            if i >= n || (hermite && n >= 3 && i == n - 1 && j >= 2) {
                assert_eq!(b.to_bits(), sentinel[j].to_bits(), "C changed sentinel");
            } else {
                // Do not let NaN sentinels or overflow make a comparison vacuous.
                assert!(b.is_finite(), "case={case} n={n} [{i}][{j}] C={b:?}");
            }
            // Stronger than the required 1e-5 relative tolerance, including zeros.
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "case={case} n={n} [{i}][{j}] rust={a:?} c={b:?}"
            );
        }
    }
}

fn check(c: CFn, f: RustFn, hermite: bool) {
    // Explicit endpoint limiter, constant-slope, asymmetric Hermite, plateau,
    // sign-change, signed-zero, and tiny-ordinate cases in both x directions.
    let fixtures: &[(&[f32], &[f32])] = &[
        (&[0.0], &[-0.0]),
        (&[-2.0, 3.0], &[0.0, -0.0]),
        (&[0.0, 1.0, 2.0], &[0.0, 1.0, -9.0]),
        (&[0.0, 1.0, 2.0], &[0.0, -10.0, -9.0]),
        (&[0.0, 1.0, 2.0], &[0.0, 1.0, 11.0]),
        (&[0.0, 1.0, 2.0], &[0.0, 10.0, 11.0]),
        (&[-3.0, -2.0, 0.0, 4.0, 5.0], &[-9.0, -7.0, -3.0, 5.0, 7.0]),
        (&[-3.0, -2.0, 0.0, 4.0, 5.0], &[2.0, 2.0, -1.0, -1.0, 3.0]),
        (&[0.0, 0.001, 1.0, 100.0], &[0.0, -1.0, 2.0, -3.0]),
        (&[0.0, 1.0, 2.0, 3.0], &[0.0, -0.0, 0.0, -0.0]),
        (&[0.0, 1.0, 2.0], &[0.0, 1e-40, -1e-40]),
    ];
    for (case, &(x, y)) in fixtures.iter().enumerate() {
        check_case(c, f, hermite, x, y, case);
        let descending: Vec<f32> = x.iter().map(|v| -v).collect();
        check_case(c, f, hermite, &descending, y, case);
    }

    let mut rng = XorShift::new(12345);
    let sizes = [1, 2, 3, 4, 5, 7, 16, 31, 64, 127, 257, 1024, 4096];
    let xscales = [1e-4f32, 0.01, 1.0, 100.0, 1e4];
    let yscales = [1e-6f32, 0.001, 1.0, 1000.0, 1e6];
    for case in 0..6500 {
        let n = sizes[case % sizes.len()];
        let xscale = xscales[rng.next_u32() as usize % xscales.len()];
        let yscale = yscales[rng.next_u32() as usize % yscales.len()];
        let dir = if rng.next_u32() & 1 == 0 { 1.0 } else { -1.0 };
        let shape = (case / sizes.len()) % 10;
        let mut x = vec![0.0; n];
        let mut y = vec![0.0; n];
        x[0] = rng.uniform_f32(-10.0, 10.0) * xscale;
        for i in 0..n {
            if i > 0 {
                let spacing = if case % 3 == 0 {
                    // A ratio of up to 512 between neighboring intervals.
                    [0.0625, 1.0, 32.0][rng.next_u32() as usize % 3]
                } else {
                    rng.uniform_f32(0.25, 2.25)
                };
                x[i] = x[i - 1] + dir * xscale * spacing;
            }
            y[i] = match shape {
                0 => 0.0,
                1 => -0.0,
                2 => -7.0 * yscale,
                3 => (x[i] / xscale * 0.5 - 2.0) * yscale,
                4 => {
                    if i == 0 {
                        0.0
                    } else {
                        y[i - 1] + rng.uniform_f32(0.1, 2.0) * yscale
                    }
                }
                5 => {
                    if i == 0 {
                        0.0
                    } else {
                        y[i - 1] - rng.uniform_f32(0.1, 2.0) * yscale
                    }
                }
                6 => (i / 3) as f32 * yscale,
                7 => {
                    if i % 2 == 0 {
                        yscale
                    } else {
                        -yscale
                    }
                }
                _ => rng.uniform_f32(-10.0, 10.0) * yscale,
            };
        }
        check_case(c, f, hermite, &x, &y, case + fixtures.len());
    }
}

#[test]
fn cakima_matches_c() {
    check(cakima, rust::cakima, false);
}
#[test]
fn chermite_matches_c() {
    check(chermite, rust::chermite, true);
}
#[test]
fn cmonot_matches_c() {
    check(cmonot, rust::cmonot, false);
}
#[test]
fn csplin_matches_c() {
    check(csplin, rust::csplin, false);
}
