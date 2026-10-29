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

//! `complex` and `dcomplex` (cwp.h) with the cwp/lib complex arithmetic
//! API. Written from scratch with textbook formulas (Smith's division,
//! scaled modulus); not a translation of cwp/lib/complex.c. Results can
//! differ from the C library in the last bit for division, modulus, sqrt
//! and transcendental functions.
//!
//! | C                         | Rust                                   |
//! |---------------------------|----------------------------------------|
//! | `complex z; z.r, z.i`     | `Complex { r, i }` (`#[repr(C)]`, 8 bytes) |
//! | `dcomplex`                | `DComplex { r, i }` (16 bytes)          |
//! | `cmplx(r,i)`, `cadd`, `csub`, `cmul`, `cdiv`, `rcabs`, `conjg`, `cneg`, `cinv`, `crmul` | same names |
//! | `cwp_csqrt`, `cwp_cexp`, `cwp_clog`, `cipow`, `crpow`, `rcpow`, `ccpow`, `cwp_ccos`, `cwp_csin`, `cwp_ccosh`, `cwp_csinh`, `cwp_cexp1` | same names |
//! | `dcmplx`, `dcadd`, ... `dclog` | same names                        |
//!
//! The operators `+ - *` and `/` are also implemented (with the same
//! formulas as `cadd`, `csub`, `cmul`, `cdiv`).

use std::ops::{Add, Div, Mul, Neg, Sub};

/// Single precision complex number (C `complex`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    pub r: f32,
    pub i: f32,
}

/// Double precision complex number (C `dcomplex`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DComplex {
    pub r: f64,
    pub i: f64,
}

const _: () = assert!(std::mem::size_of::<Complex>() == 8);
const _: () = assert!(std::mem::size_of::<DComplex>() == 16);

macro_rules! complex_impl {
    ($C:ident, $T:ty,
     $cmplx:ident, $cadd:ident, $csub:ident, $cmul:ident, $cdiv:ident, $rcabs:ident,
     $conjg:ident, $cneg:ident, $cinv:ident, $csqrt:ident, $cexp:ident, $crmul:ident,
     $cipow:ident, $clog:ident, $ccos:ident, $csin:ident, $ccosh:ident, $csinh:ident,
     $cexp1:ident, $ccpow:ident) => {
        /// Make a complex number.
        #[inline]
        pub fn $cmplx(re: $T, im: $T) -> $C {
            $C { r: re, i: im }
        }
        /// a + b
        #[inline]
        pub fn $cadd(a: $C, b: $C) -> $C {
            $C { r: a.r + b.r, i: a.i + b.i }
        }
        /// a - b
        #[inline]
        pub fn $csub(a: $C, b: $C) -> $C {
            $C { r: a.r - b.r, i: a.i - b.i }
        }
        /// a * b
        #[inline]
        pub fn $cmul(a: $C, b: $C) -> $C {
            $C { r: a.r * b.r - a.i * b.i, i: a.i * b.r + a.r * b.i }
        }
        /// a / b (Smith's algorithm); b == 0 gives IEEE inf/NaN, no error.
        pub fn $cdiv(a: $C, b: $C) -> $C {
            if b.r.abs() >= b.i.abs() {
                let ratio = b.i / b.r;
                let den = b.r + ratio * b.i;
                $C { r: (a.r + ratio * a.i) / den, i: (a.i - ratio * a.r) / den }
            } else {
                let ratio = b.r / b.i;
                let den = b.i + ratio * b.r;
                $C { r: (a.r * ratio + a.i) / den, i: (a.i * ratio - a.r) / den }
            }
        }
        /// |z| (scaled to avoid overflow).
        pub fn $rcabs(z: $C) -> $T {
            let x = z.r.abs();
            let y = z.i.abs();
            if x == 0.0 {
                y
            } else if y == 0.0 {
                x
            } else if x > y {
                let q = (y / x) as f64;
                (x as f64 * (1.0 + q * q).sqrt()) as $T
            } else {
                let q = (x / y) as f64;
                (y as f64 * (1.0 + q * q).sqrt()) as $T
            }
        }
        /// complex conjugate
        #[inline]
        pub fn $conjg(z: $C) -> $C {
            $C { r: z.r, i: -z.i }
        }
        /// -z
        #[inline]
        pub fn $cneg(z: $C) -> $C {
            $C { r: -z.r, i: -z.i }
        }
        /// 1/z
        pub fn $cinv(z: $C) -> $C {
            $cdiv($C { r: 1.0, i: 0.0 }, z)
        }
        /// principal square root
        pub fn $csqrt(z: $C) -> $C {
            if z.r == 0.0 && z.i == 0.0 {
                return $C { r: 0.0, i: 0.0 };
            }
            let m = $rcabs(z) as f64;
            let x = z.r as f64;
            if x >= 0.0 {
                let t = ((m + x) * 0.5).sqrt();
                $C { r: t as $T, i: (z.i as f64 / (2.0 * t)) as $T }
            } else {
                let t = ((m - x) * 0.5).sqrt();
                let t = if z.i < 0.0 { -t } else { t };
                $C { r: (z.i as f64 / (2.0 * t)) as $T, i: t as $T }
            }
        }
        /// e^z
        pub fn $cexp(z: $C) -> $C {
            // SU stores exp(z.r) in a variable of the element type before multiplying
            // (a float for `complex`); keep that rounding or results drift by 1 ulp.
            let a = (z.r as f64).exp() as $T as f64;
            $C { r: (a * (z.i as f64).cos()) as $T, i: (a * (z.i as f64).sin()) as $T }
        }
        /// a * x (real x)
        #[inline]
        pub fn $crmul(a: $C, x: $T) -> $C {
            $C { r: a.r * x, i: a.i * x }
        }
        /// a^p for integer p (repeated squaring).
        pub fn $cipow(a: $C, p: i32) -> $C {
            let mut result = $C { r: 1.0, i: 0.0 };
            let mut base = a;
            let mut n = p.unsigned_abs();
            while n > 0 {
                if n & 1 == 1 {
                    result = $cmul(result, base);
                }
                base = $cmul(base, base);
                n >>= 1;
            }
            if p < 0 {
                $cinv(result)
            } else {
                result
            }
        }
        /// principal natural logarithm
        pub fn $clog(a: $C) -> $C {
            let m = $rcabs(a) as f64;
            $C { r: m.ln() as $T, i: (a.i as f64).atan2(a.r as f64) as $T }
        }
        /// cos(a)
        pub fn $ccos(a: $C) -> $C {
            let (x, y) = (a.r as f64, a.i as f64);
            $C { r: (x.cos() * y.cosh()) as $T, i: (-x.sin() * y.sinh()) as $T }
        }
        /// sin(a)
        pub fn $csin(a: $C) -> $C {
            let (x, y) = (a.r as f64, a.i as f64);
            $C { r: (x.sin() * y.cosh()) as $T, i: (x.cos() * y.sinh()) as $T }
        }
        /// cosh(a)
        pub fn $ccosh(a: $C) -> $C {
            let (x, y) = (a.r as f64, a.i as f64);
            $C { r: (x.cosh() * y.cos()) as $T, i: (x.sinh() * y.sin()) as $T }
        }
        /// sinh(a)
        pub fn $csinh(a: $C) -> $C {
            let (x, y) = (a.r as f64, a.i as f64);
            $C { r: (x.sinh() * y.cos()) as $T, i: (x.cosh() * y.sin()) as $T }
        }
        /// e^(i a) for complex a: exp(-a.i) * (cos a.r, sin a.r)
        pub fn $cexp1(a: $C) -> $C {
            $cexp($C { r: -a.i, i: a.r })
        }
        /// a^p for complex p (principal branch); 0^p = 0.
        pub fn $ccpow(a: $C, p: $C) -> $C {
            if a.r == 0.0 && a.i == 0.0 {
                return $C { r: 0.0, i: 0.0 };
            }
            $cexp($cmul(p, $clog(a)))
        }

        impl Add for $C {
            type Output = $C;
            fn add(self, b: $C) -> $C {
                $cadd(self, b)
            }
        }
        impl Sub for $C {
            type Output = $C;
            fn sub(self, b: $C) -> $C {
                $csub(self, b)
            }
        }
        impl Mul for $C {
            type Output = $C;
            fn mul(self, b: $C) -> $C {
                $cmul(self, b)
            }
        }
        impl Div for $C {
            type Output = $C;
            fn div(self, b: $C) -> $C {
                $cdiv(self, b)
            }
        }
        impl Neg for $C {
            type Output = $C;
            fn neg(self) -> $C {
                $cneg(self)
            }
        }
    };
}

complex_impl!(Complex, f32, cmplx, cadd, csub, cmul, cdiv, rcabs, conjg, cneg, cinv,
    cwp_csqrt, cwp_cexp, crmul, cipow, cwp_clog, cwp_ccos, cwp_csin, cwp_ccosh, cwp_csinh,
    cwp_cexp1, ccpow);
complex_impl!(DComplex, f64, dcmplx, dcadd, dcsub, dcmul, dcdiv, drcabs, dconjg, dcneg, dcinv,
    dcsqrt, dcexp, dcrmul, dcipow, dclog, dccos, dcsin, dccosh, dcsinh, dcexp1, dcdcpow);

/// a^p for real p.
pub fn crpow(a: Complex, p: f32) -> Complex {
    ccpow(a, cmplx(p, 0.0))
}
/// a^p for real a, complex p.
pub fn rcpow(a: f32, p: Complex) -> Complex {
    ccpow(cmplx(a, 0.0), p)
}
/// a^p for real p (dcomplex; p is float in cwp.h).
pub fn dcrpow(a: DComplex, p: f32) -> DComplex {
    dcdcpow(a, dcmplx(p as f64, 0.0))
}
/// a^p for real a (float in cwp.h), complex p.
pub fn rdcpow(a: f32, p: DComplex) -> DComplex {
    dcdcpow(dcmplx(a as f64, 0.0), p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic() {
        let a = cmplx(1.0, 2.0);
        let b = cmplx(3.0, -1.0);
        assert_eq!(cmul(a, b), cmplx(5.0, 5.0));
        let q = cdiv(cmplx(5.0, 5.0), b);
        assert!((q.r - 1.0).abs() < 1e-6 && (q.i - 2.0).abs() < 1e-6);
        assert_eq!(rcabs(cmplx(3.0, 4.0)), 5.0);
        let s = cwp_csqrt(cmplx(-4.0, 0.0));
        assert_eq!(s, cmplx(0.0, 2.0));
        let e = cwp_cexp(cmplx(0.0, std::f32::consts::PI));
        assert!((e.r + 1.0).abs() < 1e-6);
        assert_eq!(cipow(a, 2), cmul(a, a));
        assert_eq!(a + b, cadd(a, b));
        let l = dclog(dcmplx(-1.0, 0.0));
        assert!((l.i - std::f64::consts::PI).abs() < 1e-15);
    }
}
