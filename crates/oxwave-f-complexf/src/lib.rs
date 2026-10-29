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

//! Port of cwp/lib/complexf.c.
//! Credits: Dave Hale, original version in C++; Gabriel Alvarez, translation to C.
//!
//! Float intermediates and double-precision math calls follow the C evaluation
//! order, including overflow and underflow in the unscaled squared magnitudes.

use oxwave_cwp::complex::cmul;
pub use oxwave_cwp::complex::Complex;

// complex.c uses an unscaled inverse, including its float overflow behavior.
fn cinv(z: Complex) -> Complex {
    let s = (1.0 / (z.r * z.r + z.i * z.i) as f64) as f32;
    Complex {
        r: z.r * s,
        i: -z.i * s,
    }
}

/// C prototype: `complex cipow(complex a, int p);`
///
/// Raises `a` to an integer power by repeated squaring. A zero exponent returns
/// one, including for a zero base; any other exponent of zero returns zero.
/// Like the original C, a nonzero base with `p == i32::MIN` is outside the
/// supported domain: signed negation overflows in C and its loop may not finish.
pub fn cipow(a: Complex, mut p: i32) -> Complex {
    if p == 0 {
        return Complex { r: 1.0, i: 0.0 };
    }
    if a.r == 0.0 && a.i == 0.0 {
        return Complex { r: 0.0, i: 0.0 };
    }
    let mut res = Complex { r: 1.0, i: 0.0 };
    let mut b = a;
    if p < 0 {
        p = p.wrapping_neg();
        b = cinv(b);
    }
    loop {
        if p & 1 != 0 {
            res = cmul(res, b);
        }
        p >>= 1;
        if p == 0 {
            return res;
        }
        b = cmul(b, b);
    }
}

/// C prototype: `complex crpow(complex a, float p);`
///
/// Raises `a` to a real power using its principal argument.
pub fn crpow(a: Complex, p: f32) -> Complex {
    if p == 0.0 {
        return Complex { r: 1.0, i: 0.0 };
    }
    if a.r == 0.0 && a.i == 0.0 {
        return Complex { r: 0.0, i: 0.0 };
    }
    let (ar, ai) = (a.r, a.i);
    let amp = (0.5 * p as f64 * ((ar * ar + ai * ai) as f64).ln()).exp() as f32;
    let phs = (p as f64 * (ai as f64).atan2(ar as f64)) as f32;
    Complex {
        r: (amp as f64 * (phs as f64).cos()) as f32,
        i: (amp as f64 * (phs as f64).sin()) as f32,
    }
}

/// C prototype: `complex rcpow(float a, complex p);`
///
/// Uses `0.5 * log(a*a)`, so negative real bases behave like their magnitudes,
/// as in the C routine.
pub fn rcpow(a: f32, p: Complex) -> Complex {
    if p.r == 0.0 && p.i == 0.0 {
        return Complex { r: 1.0, i: 0.0 };
    }
    if a == 0.0 {
        return Complex { r: 0.0, i: 0.0 };
    }
    let (pr, pi) = (p.r, p.i);
    let loga = (0.5 * ((a * a) as f64).ln()) as f32;
    let amp = ((pr * loga) as f64).exp() as f32;
    let phs = pi * loga;
    Complex {
        r: (amp as f64 * (phs as f64).cos()) as f32,
        i: (amp as f64 * (phs as f64).sin()) as f32,
    }
}

/// C prototype: `complex ccpow(complex a, complex p);`
///
/// Raises `a` to a complex power using its principal argument.
pub fn ccpow(a: Complex, p: Complex) -> Complex {
    if p.r == 0.0 && p.i == 0.0 {
        return Complex { r: 1.0, i: 0.0 };
    }
    if a.r == 0.0 && a.i == 0.0 {
        return Complex { r: 0.0, i: 0.0 };
    }
    let (pr, pi, ar, ai) = (p.r, p.i, a.r, a.i);
    let loga = (0.5 * ((ar * ar + ai * ai) as f64).ln()) as f32;
    let arga = (ai as f64).atan2(ar as f64) as f32;
    let amp = ((pr * loga - pi * arga) as f64).exp() as f32;
    let phs = pr * arga + pi * loga;
    Complex {
        r: (amp as f64 * (phs as f64).cos()) as f32,
        i: (amp as f64 * (phs as f64).sin()) as f32,
    }
}

/// C prototype: `complex cwp_ccos(complex a);`
/// Computes the complex cosine.
pub fn cwp_ccos(a: Complex) -> Complex {
    Complex {
        r: ((a.r as f64).cos() * (a.i as f64).cosh()) as f32,
        i: (-(a.r as f64).sin() * (a.i as f64).sinh()) as f32,
    }
}
/// C prototype: `complex cwp_csin(complex a);`
/// Computes the complex sine.
pub fn cwp_csin(a: Complex) -> Complex {
    Complex {
        r: ((a.r as f64).sin() * (a.i as f64).cosh()) as f32,
        i: ((a.r as f64).cos() * (a.i as f64).sinh()) as f32,
    }
}
/// C prototype: `complex cwp_ccosh(complex a);`
/// Computes the complex hyperbolic cosine.
pub fn cwp_ccosh(a: Complex) -> Complex {
    Complex {
        r: ((a.i as f64).cos() * (a.r as f64).cosh()) as f32,
        i: ((a.i as f64).sin() * (a.r as f64).sinh()) as f32,
    }
}
/// C prototype: `complex cwp_csinh(complex a);`
/// Computes the complex hyperbolic sine.
pub fn cwp_csinh(a: Complex) -> Complex {
    Complex {
        r: ((a.i as f64).cos() * (a.r as f64).sinh()) as f32,
        i: ((a.i as f64).sin() * (a.r as f64).cosh()) as f32,
    }
}
/// C prototype: `complex cwp_cexp1(complex a);`
/// Computes `exp(a)`, rounding the real exponential to `f32` before multiplying.
pub fn cwp_cexp1(a: Complex) -> Complex {
    let r = (a.r as f64).exp() as f32;
    Complex {
        r: (r as f64 * (a.i as f64).cos()) as f32,
        i: (r as f64 * (a.i as f64).sin()) as f32,
    }
}
/// C prototype: `complex cwp_clog(complex a);`
/// Computes the principal complex logarithm with the C routine's unscaled norm.
pub fn cwp_clog(a: Complex) -> Complex {
    let (ar, ai) = (a.r, a.i);
    let h = ((ar * ar + ai * ai) as f64).sqrt() as f32;
    Complex {
        r: (h as f64).ln() as f32,
        i: (ai as f64).atan2(ar as f64) as f32,
    }
}
