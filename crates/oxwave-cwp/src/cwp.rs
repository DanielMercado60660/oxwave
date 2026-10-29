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

//! Constants and macro equivalents from `cwp/include/cwp.h`.
//!
//! The C macros are type-generic and evaluate with C's usual arithmetic
//! conversions; these functions reproduce the exact results:
//!
//! | C macro          | Rust                     | note |
//! |------------------|--------------------------|------|
//! | `PI`, `D_PI`     | `PI`, `D_PI` (f64)        | `PI` is a *double* literal: `float*PI` is computed in double |
//! | `ABS(x)`         | `abs(x)`                 | `x < 0 ? -x : x` (so `ABS(-0.0) == -0.0`, NaN unchanged) |
//! | `MAX(x,y)`       | `max(x, y)`              | `x > y ? x : y` (NaN in x -> y) -- NOT `f32::max` |
//! | `MIN(x,y)`       | `min(x, y)`              | `x < y ? x : y` -- NOT `f32::min` |
//! | `SGN(x)`         | `sgn(x) -> f64`          | `x < 0 ? -1.0 : 1.0` (a double) |
//! | `NINT(x)`        | `nint(x) -> i32`         | `(int)(x > 0 ? x+0.5 : x-0.5)` in double |
//! | `CLOSETO(x,y)`   | `closeto(x, y)`          | `ABS(x-y) <= FLT_EPSILON*ABS(y)` |
//! | `ISODD(n)`       | `isodd(n)`               | `n & 1` |
//! | `SQR(x)`         | `x * x`                  | |
//! | `FSIZE` etc.     | `FSIZE`, `ISIZE`, `DSIZE`, `CSIZE` | |
//!
//! Note `nint`'s argument: a C `float` is promoted to double *before* the
//! 0.5 is added (the literal 0.5 is a double), so `nint` takes `f64`; pass
//! `x as f64` for floats -- this is exact and matches C.

use std::ops::Neg;

pub const PI: f64 = 3.141592653589793;
pub const D_PI: f64 = 3.1415926535897932385;
pub const GOLDEN_RATIO: f64 = 1.618034;
pub const TRUE: i32 = 1;
pub const FALSE: i32 = 0;
pub const YES: i32 = 1;
pub const NO: i32 = 0;
pub const FLT_EPSILON: f32 = f32::EPSILON;
pub const FLT_MAX: f32 = f32::MAX;
pub const FLT_MIN: f32 = f32::MIN_POSITIVE;
pub const DBL_EPSILON: f64 = f64::EPSILON;
pub const DBL_MAX: f64 = f64::MAX;
pub const DBL_MIN: f64 = f64::MIN_POSITIVE;
pub const INT_MAX: i32 = i32::MAX;
pub const INT_MIN: i32 = i32::MIN;
pub const ISIZE: usize = 4;
pub const FSIZE: usize = 4;
pub const DSIZE: usize = 8;
pub const CSIZE: usize = 8;

/// `ABS(x)`: `x < 0 ? -x : x`.
#[inline]
pub fn abs<T: PartialOrd + Neg<Output = T> + Default + Copy>(x: T) -> T {
    if x < T::default() {
        -x
    } else {
        x
    }
}

/// `MAX(x,y)`: `x > y ? x : y`.
#[inline]
pub fn max<T: PartialOrd + Copy>(x: T, y: T) -> T {
    if x > y {
        x
    } else {
        y
    }
}

/// `MIN(x,y)`: `x < y ? x : y`.
#[inline]
pub fn min<T: PartialOrd + Copy>(x: T, y: T) -> T {
    if x < y {
        x
    } else {
        y
    }
}

/// `SGN(x)`: `x < 0 ? -1.0 : 1.0` (a double).
#[inline]
pub fn sgn<T: PartialOrd + Default>(x: T) -> f64 {
    if x < T::default() {
        -1.0
    } else {
        1.0
    }
}

/// `NINT(x)`: `(int)((x)>0.0?(x)+0.5:(x)-0.5)`, evaluated in double.
/// Out-of-range values saturate (as the arm64 conversion instruction does).
#[inline]
pub fn nint(x: f64) -> i32 {
    (if x > 0.0 { x + 0.5 } else { x - 0.5 }) as i32
}

/// `CLOSETO(x,y)` for floats: `ABS(x-y) <= FLT_EPSILON*ABS(y)`.
#[inline]
pub fn closeto(x: f32, y: f32) -> bool {
    abs(x - y) <= FLT_EPSILON * abs(y)
}

/// `ISODD(n)`.
#[inline]
pub fn isodd(n: i32) -> bool {
    n & 1 != 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn macros() {
        assert_eq!(nint(2.5), 3);
        assert_eq!(nint(-2.5), -3);
        assert_eq!(nint(0.49999997f32 as f64), 0);
        assert_eq!(nint(-0.0), 0);
        assert_eq!(max(f32::NAN, 1.0), 1.0);
        assert!(max(1.0, f32::NAN).is_nan());
        assert_eq!(min(3, -4), -4);
        assert_eq!(abs(-0.0f32).to_bits(), (-0.0f32).to_bits());
        assert_eq!(sgn(-3i32), -1.0);
        assert_eq!(sgn(0.0f32), 1.0);
        assert!(closeto(1.0, 1.0 + f32::EPSILON / 2.0));
        assert!(isodd(-3));
    }
}
