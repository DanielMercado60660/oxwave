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

//! Differential tests for cwp/lib/complexf.c.
//! Credits: Dave Hale, original version in C++; Gabriel Alvarez, translation to C.

use oxwave_cref::XorShift;
use oxwave_f_complexf::{self as rust, Complex};

mod c {
    use super::Complex;
    extern "C" {
        pub fn cipow(a: Complex, p: i32) -> Complex;
        pub fn crpow(a: Complex, p: f32) -> Complex;
        pub fn rcpow(a: f32, p: Complex) -> Complex;
        pub fn ccpow(a: Complex, p: Complex) -> Complex;
        pub fn cwp_ccos(a: Complex) -> Complex;
        pub fn cwp_csin(a: Complex) -> Complex;
        pub fn cwp_ccosh(a: Complex) -> Complex;
        pub fn cwp_csinh(a: Complex) -> Complex;
        pub fn cwp_cexp1(a: Complex) -> Complex;
        pub fn cwp_clog(a: Complex) -> Complex;
    }
}

fn compare(a: Complex, b: Complex, input: &str) {
    // Use only the C result to set the tolerance. Nonfinite components must
    // agree in classification (and infinity sign), independently of the scale.
    let scale = [b.r, b.i]
        .into_iter()
        .filter(|x| x.is_finite())
        .fold(0.0f64, |m, x| m.max((x as f64).abs()));
    for (x, y) in [(a.r, b.r), (a.i, b.i)] {
        assert!(
            (x.is_nan() && y.is_nan())
                || x == y
                || (x.is_finite()
                    && y.is_finite()
                    && ((x as f64) - (y as f64)).abs() <= 1e-5 * scale),
            "{input}: Rust {a:?}, C {b:?}"
        );
    }
}

fn compare_bits(a: Complex, b: Complex, input: &str) {
    for (x, y) in [(a.r, b.r), (a.i, b.i)] {
        // NaN payload/sign propagation is not specified by either language.
        assert!(
            (x.is_nan() && y.is_nan()) || x.to_bits() == y.to_bits(),
            "{input}: Rust {a:?}, C {b:?}; bits {:08x} vs {:08x}",
            x.to_bits(),
            y.to_bits()
        );
    }
}

fn edges() -> [f32; 25] {
    [
        0.0,
        -0.0,
        1.0,
        -1.0,
        0.5,
        -0.5,
        2.0,
        -2.0,
        f32::EPSILON,
        f32::MIN_POSITIVE,
        -f32::MIN_POSITIVE,
        f32::from_bits(1),
        -f32::from_bits(1),
        1e-40,
        -1e-40,
        1e-20,
        -1e-20,
        1e20,
        -1e20,
        f32::MAX,
        -f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
        -f32::NAN,
    ]
}

fn samples() -> Vec<Complex> {
    let mut rng = XorShift::new(42);
    let mut v = Vec::new();
    // These APIs are scalar. Exercise singleton, prime, power-of-two and large
    // batches of calls rather than inventing array arguments absent from C.
    for scale in [1e-40, 1e-20, 1e-5, 1.0, 10.0, 80.0, 1e20] {
        for n in [1, 2, 3, 17, 32, 127, 256, 1024, 4096] {
            for _ in 0..n {
                v.push(Complex {
                    r: rng.uniform_f32(-scale, scale),
                    i: rng.uniform_f32(-scale, scale),
                });
            }
        }
    }
    // Raw bit patterns cover mixed component magnitudes, subnormals and the
    // entire f32 exponent range, without discarding nonfinite inputs.
    for _ in 0..8192 {
        v.push(Complex {
            r: rng.any_f32(),
            i: rng.any_f32(),
        });
    }
    for r in edges() {
        for i in edges() {
            v.push(Complex { r, i });
        }
    }
    // Approach the negative-real-axis branch cut and square-magnitude and
    // exp/sinh/cosh underflow/overflow boundaries from both directions.
    for x in [
        1.0f32,
        3.743392e-23,
        1.8446743e19,
        88.72284,
        89.415985,
        103.97208,
        710.0,
    ] {
        for bits in [x.to_bits() - 1, x.to_bits(), x.to_bits() + 1] {
            for r in [f32::from_bits(bits), -f32::from_bits(bits)] {
                for i in [0.0, -0.0, f32::from_bits(1), -f32::from_bits(1), 1.0, -1.0] {
                    v.push(Complex { r, i });
                    v.push(Complex { r: i, i: r });
                }
            }
        }
    }
    v
}
#[test]
fn cwp_ccos_matches() {
    for a in samples() {
        compare(
            rust::cwp_ccos(a),
            unsafe { c::cwp_ccos(a) },
            &format!("a={a:?}"),
        );
    }
}
#[test]
fn cwp_csin_matches() {
    for a in samples() {
        compare(
            rust::cwp_csin(a),
            unsafe { c::cwp_csin(a) },
            &format!("a={a:?}"),
        );
    }
}
#[test]
fn cwp_ccosh_matches() {
    for a in samples() {
        compare(
            rust::cwp_ccosh(a),
            unsafe { c::cwp_ccosh(a) },
            &format!("a={a:?}"),
        );
    }
}
#[test]
fn cwp_csinh_matches() {
    for a in samples() {
        compare(
            rust::cwp_csinh(a),
            unsafe { c::cwp_csinh(a) },
            &format!("a={a:?}"),
        );
    }
}
#[test]
fn cwp_cexp1_matches() {
    for a in samples() {
        compare(
            rust::cwp_cexp1(a),
            unsafe { c::cwp_cexp1(a) },
            &format!("a={a:?}"),
        );
    }
}
#[test]
fn cwp_clog_matches() {
    for a in samples() {
        compare(
            rust::cwp_clog(a),
            unsafe { c::cwp_clog(a) },
            &format!("a={a:?}"),
        );
    }
}

#[test]
fn cipow_matches() {
    let mut rng = XorShift::new(21);
    for a in samples() {
        let p = rng.range_i32(-32, 32);
        for p in [
            p,
            -3,
            -1,
            0,
            1,
            2,
            17,
            -65537,
            65537,
            i32::MIN + 1,
            i32::MAX,
        ] {
            compare_bits(
                rust::cipow(a, p),
                unsafe { c::cipow(a, p) },
                &format!("a={a:?}, p={p}"),
            );
        }
    }
    // INT_MIN negation is undefined in C for nonzero bases. Zero bases return
    // before that operation and are defined for the entire i32 exponent range.
    for r in [0.0, -0.0] {
        for i in [0.0, -0.0] {
            let a = Complex { r, i };
            compare_bits(
                rust::cipow(a, i32::MIN),
                unsafe { c::cipow(a, i32::MIN) },
                &format!("a={a:?}, p=INT_MIN"),
            );
        }
    }
    // Exact unit bases allow huge exponents without overflow masking the
    // repeated-squaring result.
    for a in [
        Complex { r: 1.0, i: 0.0 },
        Complex { r: -1.0, i: 0.0 },
        Complex { r: 0.0, i: 1.0 },
    ] {
        for p in [-1073741824, -1024, -127, -2, 3, 127, 1024, 1073741824] {
            compare_bits(
                rust::cipow(a, p),
                unsafe { c::cipow(a, p) },
                &format!("a={a:?}, p={p}"),
            );
        }
    }
}
#[test]
fn crpow_matches() {
    let mut rng = XorShift::new(22);
    for a in samples() {
        let p = rng.uniform_f32(-8.0, 8.0);
        for p in [p, rng.any_f32(), 0.0, -0.0, -1.0, 0.5, 2.0] {
            compare(
                rust::crpow(a, p),
                unsafe { c::crpow(a, p) },
                &format!("a={a:?}, p={p}"),
            );
        }
    }
    for a in [
        Complex { r: 0.0, i: -0.0 },
        Complex { r: -2.0, i: 0.0 },
        Complex { r: -2.0, i: -0.0 },
        Complex { r: 0.5, i: -0.5 },
    ] {
        for p in edges() {
            compare(
                rust::crpow(a, p),
                unsafe { c::crpow(a, p) },
                &format!("a={a:?}, p={p}"),
            );
        }
    }
}
#[test]
fn rcpow_matches() {
    let mut rng = XorShift::new(23);
    for p in samples() {
        let a = rng.uniform_f32(-10.0, 10.0);
        for a in [a, rng.any_f32(), 0.0, -0.0, -1.0, 1e-40, f32::MAX] {
            compare(
                rust::rcpow(a, p),
                unsafe { c::rcpow(a, p) },
                &format!("a={a:?}, p={p:?}"),
            );
        }
    }
    for a in edges() {
        for r in edges() {
            for i in edges() {
                let p = Complex { r, i };
                compare(
                    rust::rcpow(a, p),
                    unsafe { c::rcpow(a, p) },
                    &format!("a={a:?}, p={p:?}"),
                );
            }
        }
    }
}
#[test]
fn ccpow_matches() {
    let mut rng = XorShift::new(24);
    for a in samples() {
        let p = Complex {
            r: rng.uniform_f32(-5.0, 5.0),
            i: rng.uniform_f32(-5.0, 5.0),
        };
        let raw_p = Complex {
            r: rng.any_f32(),
            i: rng.any_f32(),
        };
        for p in [
            p,
            raw_p,
            Complex { r: 0.0, i: 0.0 },
            Complex { r: -0.0, i: -0.0 },
            Complex { r: -1.0, i: 2.0 },
        ] {
            compare(
                rust::ccpow(a, p),
                unsafe { c::ccpow(a, p) },
                &format!("a={a:?}, p={p:?}"),
            );
        }
    }
    for a in [
        Complex { r: 0.0, i: -0.0 },
        Complex { r: -2.0, i: 0.0 },
        Complex { r: -2.0, i: -0.0 },
        Complex { r: 0.5, i: -0.5 },
    ] {
        for r in edges() {
            for i in edges() {
                let p = Complex { r, i };
                compare(
                    rust::ccpow(a, p),
                    unsafe { c::ccpow(a, p) },
                    &format!("a={a:?}, p={p:?}"),
                );
            }
        }
    }
}
