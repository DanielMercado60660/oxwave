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

//! Differential tests for cwp/lib/dblas.c, adapted from LINPACK FORTRAN.
//! Author: Dave Hale, Colorado School of Mines, 10/01/89.
use oxwave_cref::XorShift;
use oxwave_f_dblas as rust;

extern "C" {
    fn dasum(n: i32, sx: *const f64, incx: i32) -> f64;
    fn daxpy(n: i32, sa: f64, sx: *const f64, incx: i32, sy: *mut f64, incy: i32);
    fn dcopy(n: i32, sx: *const f64, incx: i32, sy: *mut f64, incy: i32);
    fn ddot(n: i32, sx: *const f64, incx: i32, sy: *const f64, incy: i32) -> f64;
    fn dnrm2(n: i32, sx: *const f64, incx: i32) -> f64;
    fn dscal(n: i32, sa: f64, sx: *mut f64, incx: i32);
    fn dswap(n: i32, sx: *mut f64, incx: i32, sy: *mut f64, incy: i32);
    fn idamax(n: i32, sx: *const f64, incx: i32) -> i32;
}

struct Case {
    n: i32,
    ix: i32,
    iy: i32,
    a: f64,
    x: Vec<f64>,
    y: Vec<f64>,
}
fn cases(mut f: impl FnMut(Case)) {
    let mut rng = XorShift::new(429);
    let sizes = [
        -2, 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 31, 64, 127, 256, 1021, 4096, 8191,
    ];
    let strides = [-5i32, -3, -1, 0, 1, 2, 5];
    let scales = [0.0, 1e-200, 1e-150, 1e-20, 1.0, 1e20, 1e150, 1e200];
    let combinations = sizes.len() * strides.len() * strides.len();
    // 7,448 random cases per function. Every size/stride pair is exercised at
    // every scale, including deliberate square/product overflow and underflow.
    for k in 0..combinations * scales.len() {
        let n = sizes[k % sizes.len()];
        let ix = strides[(k / sizes.len()) % strides.len()];
        let iy = strides[(k / (sizes.len() * strides.len())) % strides.len()];
        let scale = scales[k / combinations];
        // Two trailing guard elements also get compared after mutations.
        let x = (0..storage(n, ix))
            .map(|_| rng.uniform(-1.0, 1.0) * scale)
            .collect();
        let y = (0..storage(n, iy))
            .map(|_| rng.uniform(-1.0, 1.0) * scale)
            .collect();
        let a = if k % 9 == 0 {
            0.0
        } else {
            rng.uniform(-5.0, 5.0)
        };
        f(Case { n, ix, iy, a, x, y });
    }
    // Exercise every cleanup remainder, ties, signed zeros, subnormals,
    // cancellation, nonfinite values, and NaN payloads on every stride pair.
    let patterns: &[&[f64]] = &[
        &[0.0, -0.0, 0.0, -0.0],
        &[1.0, -1.0, 7.0, -7.0, 2.0, -7.0],
        &[1e16, 1.0, -1e16, 2.0, -2.0, -1.0],
        &[
            f64::from_bits(1),
            -f64::from_bits(1),
            f64::MIN_POSITIVE,
            -1e-200,
        ],
        &[f64::MAX, -f64::MAX, 1e200, -1e200, 1.0],
        &[
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8_0000_0000_0042),
            f64::from_bits(0xfff8_0000_0000_0123),
            1.0,
            -0.0,
        ],
    ];
    let scalars = [0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 0.5];
    for pattern in patterns {
        for ix in strides {
            for iy in strides {
                for n in 1..=18 {
                    let x = (0..storage(n, ix))
                        .map(|i| pattern[(i + n as usize) % pattern.len()])
                        .collect();
                    let y = (0..storage(n, iy))
                        .map(|i| pattern[(i + 1) % pattern.len()])
                        .collect();
                    let a = scalars[n as usize % scalars.len()];
                    f(Case { n, ix, iy, a, x, y });
                }
            }
        }
    }
}
fn storage(n: i32, inc: i32) -> usize {
    3 + (n.max(1) - 1) as usize * inc.unsigned_abs() as usize
}
fn unary_offset(n: i32, inc: i32) -> usize {
    if n > 0 && inc < 0 {
        ((n - 1) * -inc) as usize
    } else {
        0
    }
}
fn close(a: f64, b: f64) {
    if a.is_nan() || b.is_nan() {
        assert!(a.is_nan() && b.is_nan(), "Rust {a:?}, C {b:?}");
    } else if a.is_infinite() || b.is_infinite() {
        assert_eq!(a, b, "infinite results must have the same sign");
    } else {
        assert!(
            (a - b).abs() <= 1e-5 * a.abs().max(b.abs()),
            "Rust {a:?}, C {b:?}"
        );
    }
}
fn exact(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (i, (a, b)) in a.iter().zip(b).enumerate() {
        assert_eq!(a.to_bits(), b.to_bits(), "element {i}: Rust {a:?}, C {b:?}");
    }
}

#[test]
fn test_dasum() {
    cases(|c| {
        // Unary routines do not rebase negative increments inside C. Supply
        // the high-address pointer while Rust receives the whole storage.
        let p = unsafe { c.x.as_ptr().add(unary_offset(c.n, c.ix)) };
        close(rust::dasum(c.n, &c.x, c.ix), unsafe { dasum(c.n, p, c.ix) });
    });
}
#[test]
fn test_dnrm2() {
    cases(|c| {
        let p = unsafe { c.x.as_ptr().add(unary_offset(c.n, c.ix)) };
        close(rust::dnrm2(c.n, &c.x, c.ix), unsafe { dnrm2(c.n, p, c.ix) });
    });
}
#[test]
fn test_idamax() {
    cases(|c| {
        let p = unsafe { c.x.as_ptr().add(unary_offset(c.n, c.ix)) };
        assert_eq!(rust::idamax(c.n, &c.x, c.ix), unsafe {
            idamax(c.n, p, c.ix)
        });
    });
}
#[test]
fn test_ddot() {
    cases(|c| {
        close(rust::ddot(c.n, &c.x, c.ix, &c.y, c.iy), unsafe {
            ddot(c.n, c.x.as_ptr(), c.ix, c.y.as_ptr(), c.iy)
        })
    });
}
#[test]
fn test_daxpy() {
    cases(|c| {
        let mut r = c.y.clone();
        let mut v = c.y.clone();
        rust::daxpy(c.n, c.a, &c.x, c.ix, &mut r, c.iy);
        unsafe { daxpy(c.n, c.a, c.x.as_ptr(), c.ix, v.as_mut_ptr(), c.iy) };
        for (a, b) in r.iter().zip(&v) {
            close(*a, *b);
        }
    });
}
#[test]
fn test_dcopy() {
    cases(|c| {
        let mut r = c.y.clone();
        let mut v = c.y.clone();
        rust::dcopy(c.n, &c.x, c.ix, &mut r, c.iy);
        unsafe { dcopy(c.n, c.x.as_ptr(), c.ix, v.as_mut_ptr(), c.iy) };
        exact(&r, &v);
    });
}
#[test]
fn test_dscal() {
    cases(|c| {
        let mut r = c.x.clone();
        let mut v = c.x.clone();
        rust::dscal(c.n, c.a, &mut r, c.ix);
        let p = unsafe { v.as_mut_ptr().add(unary_offset(c.n, c.ix)) };
        unsafe { dscal(c.n, c.a, p, c.ix) };
        exact(&r, &v);
    });
}
#[test]
fn test_dswap() {
    cases(|c| {
        let (mut x, mut y) = (c.x.clone(), c.y.clone());
        let (mut a, mut b) = (c.x.clone(), c.y.clone());
        rust::dswap(c.n, &mut x, c.ix, &mut y, c.iy);
        unsafe { dswap(c.n, a.as_mut_ptr(), c.ix, b.as_mut_ptr(), c.iy) };
        exact(&x, &a);
        exact(&y, &b);
    });
}

#[test]
fn nonpositive_sizes_and_zero_axpy_do_not_access_arrays() {
    for n in [i32::MIN, -1, 0] {
        for inc in [i32::MIN, -1, 0, 1, i32::MAX] {
            let (mut x, mut y) = ([], []);
            // C's early returns allow null pointers for these seven routines.
            let p = std::ptr::null();
            let q = std::ptr::null_mut();
            unsafe {
                assert_eq!(
                    rust::dasum(n, &x, inc).to_bits(),
                    dasum(n, p, inc).to_bits()
                );
                assert_eq!(
                    rust::ddot(n, &x, inc, &y, inc).to_bits(),
                    ddot(n, p, inc, p, inc).to_bits()
                );
                assert_eq!(
                    rust::dnrm2(n, &x, inc).to_bits(),
                    dnrm2(n, p, inc).to_bits()
                );
                rust::daxpy(n, f64::NAN, &x, inc, &mut y, inc);
                daxpy(n, f64::NAN, p, inc, q, inc);
                rust::dcopy(n, &x, inc, &mut y, inc);
                dcopy(n, p, inc, q, inc);
                rust::dscal(n, f64::NAN, &mut x, inc);
                dscal(n, f64::NAN, q, inc);
                rust::dswap(n, &mut x, inc, &mut y, inc);
                dswap(n, q, inc, q, inc);
            }
            // Unlike the other routines, idamax reads sx[0] even for n <= 0.
            let one = [f64::NAN];
            assert_eq!(rust::idamax(n, &one, inc), unsafe {
                idamax(n, one.as_ptr(), inc)
            });
        }
    }
    for a in [0.0, -0.0] {
        rust::daxpy(i32::MAX, a, &[], i32::MIN, &mut [], i32::MAX);
        unsafe {
            daxpy(
                i32::MAX,
                a,
                std::ptr::null(),
                i32::MIN,
                std::ptr::null_mut(),
                i32::MAX,
            );
        }
    }
}

#[test]
fn idamax_ties_and_logical_indices() {
    let x = [1.0, -9.0, 2.0, 9.0, -3.0];
    for (inc, expected) in [(1, 1), (-1, 1), (0, 0)] {
        let p = unsafe { x.as_ptr().add(unary_offset(5, inc)) };
        assert_eq!(rust::idamax(5, &x, inc), expected);
        assert_eq!(unsafe { idamax(5, p, inc) }, expected);
    }
    // The maximum is at physical index 4 but logical index 2.
    let x = [1.0, 100.0, -3.0, 100.0, -9.0];
    assert_eq!(rust::idamax(3, &x, 2), 2);
    assert_eq!(unsafe { idamax(3, x.as_ptr(), 2) }, 2);
    let x = [f64::NAN, 1e100, f64::INFINITY];
    assert_eq!(rust::idamax(3, &x, 1), 0);
    assert_eq!(unsafe { idamax(3, x.as_ptr(), 1) }, 0);
}

#[test]
fn reductions_keep_the_c_unrolled_grouping() {
    // These lose different small terms if a block is combined with the
    // running sum before the C expression has finished evaluating.
    let mut x = [1.0; 13];
    x[0] = 1e16;
    let mut norm_x = [1.0; 13];
    norm_x[0] = 1e8;
    for n in [6, 7, 12, 13] {
        assert_eq!(
            rust::dasum(n, &x, 1).to_bits(),
            unsafe { dasum(n, x.as_ptr(), 1) }.to_bits()
        );
        assert_eq!(
            rust::dnrm2(n, &norm_x, 1).to_bits(),
            unsafe { dnrm2(n, norm_x.as_ptr(), 1) }.to_bits()
        );
    }
    let x = [1e16, 1.0, 1.0, -1e16, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let y = [1.0; 11];
    for n in [5, 6, 10, 11] {
        assert_eq!(
            rust::ddot(n, &x, 1, &y, 1).to_bits(),
            unsafe { ddot(n, x.as_ptr(), 1, y.as_ptr(), 1) }.to_bits()
        );
    }
}
