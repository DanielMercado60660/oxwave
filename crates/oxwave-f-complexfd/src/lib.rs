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

//! Port of `cwp/lib/complexfd.c`.
//! Credits: Dave Hale, original version in C++; Gabriel Alvarez, translation to C.
//!
//! The C power and logarithm routines use `float` intermediates
//! despite their double complex arguments. Those roundings, including their
//! underflow and overflow behavior, are preserved here.

pub use oxwave_cwp::complex::DComplex;
use oxwave_cwp::complex::{dcmplx, dcmul};

/// C prototype: `dcomplex dcipow(dcomplex a, int p);`
///
/// Raises `a` to an integer power by repeated squaring. Like C, zero to a
/// negative power returns zero, and zero to the zeroth power returns one.
/// For `i32::MIN`, whose negation is undefined in C, Rust uses a wider magnitude.
pub fn dcipow(a: DComplex, p: i32) -> DComplex {
    if p == 0 {
        return dcmplx(1.0, 0.0);
    }
    if a.r == 0.0 && a.i == 0.0 {
        return dcmplx(0.0, 0.0);
    }
    let mut res = dcmplx(1.0, 0.0);
    let mut b = a;
    // INT_MIN negation is undefined in C; a wider magnitude handles it.
    let mut magnitude = if p < 0 {
        // Original complexd.c dcinv formula; the frozen core's scaled inverse
        // has different rounding and overflow behavior.
        let s = 1.0 / (b.r * b.r + b.i * b.i);
        b = dcmplx(b.r * s, -b.i * s);
        -(p as i64)
    } else {
        p as i64
    };
    loop {
        if magnitude & 1 != 0 {
            res = dcmul(res, b);
        }
        magnitude >>= 1;
        if magnitude == 0 {
            return res;
        }
        b = dcmul(b, b);
    }
}

/// C prototype: `dcomplex dcrpow(dcomplex a, float p);`
///
/// Raises `a` to a real power using the C routine's single precision amplitude
/// and phase intermediates.
pub fn dcrpow(a: DComplex, p: f32) -> DComplex {
    if p == 0.0 {
        return dcmplx(1.0, 0.0);
    }
    if a.r == 0.0 && a.i == 0.0 {
        return dcmplx(0.0, 0.0);
    }
    let ar = a.r as f32;
    let ai = a.i as f32;
    let amp = (0.5 * p as f64 * ((ar * ar + ai * ai) as f64).ln()).exp() as f32;
    let phs = (p as f64 * (ai as f64).atan2(ar as f64)) as f32;
    dcmplx(
        amp as f64 * (phs as f64).cos(),
        amp as f64 * (phs as f64).sin(),
    )
}

/// C prototype: `dcomplex rdcpow(float a, dcomplex p);`
///
/// Raises a real base to a complex power. The original uses `log(a*a)`, so
/// negative bases have the same result as their positive magnitude.
pub fn rdcpow(a: f32, p: DComplex) -> DComplex {
    if p.r == 0.0 && p.i == 0.0 {
        return dcmplx(1.0, 0.0);
    }
    if a == 0.0 {
        return dcmplx(0.0, 0.0);
    }
    let pr = p.r as f32;
    let pi = p.i as f32;
    let loga = (0.5 * ((a * a) as f64).ln()) as f32;
    let amp = ((pr * loga) as f64).exp() as f32;
    let phs = pi * loga;
    dcmplx(
        amp as f64 * (phs as f64).cos(),
        amp as f64 * (phs as f64).sin(),
    )
}

/// C prototype: `dcomplex dcdcpow(dcomplex a, dcomplex p);`
///
/// Raises `a` to a complex power using single precision intermediates.
pub fn dcdcpow(a: DComplex, p: DComplex) -> DComplex {
    if p.r == 0.0 && p.i == 0.0 {
        return dcmplx(1.0, 0.0);
    }
    if a.r == 0.0 && a.i == 0.0 {
        return dcmplx(0.0, 0.0);
    }
    let pr = p.r as f32;
    let pi = p.i as f32;
    let ar = a.r as f32;
    let ai = a.i as f32;
    let loga = (0.5 * ((ar * ar + ai * ai) as f64).ln()) as f32;
    let arga = (ai as f64).atan2(ar as f64) as f32;
    let amp = ((pr * loga - pi * arga) as f64).exp() as f32;
    let phs = pr * arga + pi * loga;
    dcmplx(
        amp as f64 * (phs as f64).cos(),
        amp as f64 * (phs as f64).sin(),
    )
}

/// C prototype: `dcomplex dccos(dcomplex a);`
///
/// Computes the complex cosine.
pub fn dccos(a: DComplex) -> DComplex {
    dcmplx(a.r.cos() * a.i.cosh(), -a.r.sin() * a.i.sinh())
}
/// C prototype: `dcomplex dcsin(dcomplex a);`
///
/// Computes the complex sine.
pub fn dcsin(a: DComplex) -> DComplex {
    dcmplx(a.r.sin() * a.i.cosh(), a.r.cos() * a.i.sinh())
}
/// C prototype: `dcomplex dccosh(dcomplex a);`
///
/// Computes the complex hyperbolic cosine.
pub fn dccosh(a: DComplex) -> DComplex {
    dcmplx(a.i.cos() * a.r.cosh(), a.i.sin() * a.r.sinh())
}
/// C prototype: `dcomplex dcsinh(dcomplex a);`
///
/// Computes the complex hyperbolic sine.
pub fn dcsinh(a: DComplex) -> DComplex {
    dcmplx(a.i.cos() * a.r.sinh(), a.i.sin() * a.r.cosh())
}
/// C prototype: `dcomplex dcexp1(dcomplex a);`
///
/// Computes the complex exponential with the C routine's single precision
/// amplitude, which can overflow before a double precision exponential would.
pub fn dcexp1(a: DComplex) -> DComplex {
    let r = a.r.exp() as f32;
    dcmplx(r as f64 * a.i.cos(), r as f64 * a.i.sin())
}
/// C prototype: `dcomplex dclog(dcomplex a);`
///
/// Computes the complex logarithm after rounding the components and magnitude
/// to single precision, as in C.
pub fn dclog(a: DComplex) -> DComplex {
    let ar = a.r as f32;
    let ai = a.i as f32;
    let h = ((ar * ar + ai * ai) as f64).sqrt() as f32;
    dcmplx((h as f64).ln(), (ai as f64).atan2(ar as f64))
}
