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

//! Double precision BLAS, adapted from LINPACK FORTRAN; port of cwp/lib/dblas.c.
//! Author: Dave Hale, Colorado School of Mines, 10/01/89.
//!
//! `n` counts visited elements; increments are in elements. Slices contain all
//! visited elements. With negative increments, slices begin at the lowest
//! address. Single-vector routines start at `(n-1)*-incx` (pass a C pointer
//! to that element); two-vector routines apply C's internal starting offset.
//! Zero increments repeatedly access the same element. Insufficient slice
//! storage panics. Except for `idamax`, nonpositive `n` accesses no elements.

// Unary C routines take a pointer to this starting element. Two-vector C
// routines calculate the same offset internally, using int arithmetic.
fn start(n: i32, inc: i32) -> i32 {
    if inc < 0 && n > 0 {
        1i32.wrapping_sub(n).wrapping_mul(inc)
    } else {
        0
    }
}
// Match cwp.h ABS, including negative zero and NaNs.
fn abs(x: f64) -> f64 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

/// C: `double dasum (int n, double *sx, int incx);`
pub fn dasum(n: i32, sx: &[f64], incx: i32) -> f64 {
    let mut sum = 0.0;
    if n <= 0 {
        return sum;
    }
    if incx != 1 {
        let mut j = start(n, incx);
        for _ in 0..n {
            sum += abs(sx[j as usize]);
            j = j.wrapping_add(incx);
        }
    } else {
        let m = n % 6;
        for i in 0..m as usize {
            sum += abs(sx[i]);
        }
        for i in (m as usize..n as usize).step_by(6) {
            sum += abs(sx[i])
                + abs(sx[i + 1])
                + abs(sx[i + 2])
                + abs(sx[i + 3])
                + abs(sx[i + 4])
                + abs(sx[i + 5]);
        }
    }
    sum
}

/// C: `void daxpy (int n, double sa, double *sx, int incx, double *sy, int incy);`
pub fn daxpy(n: i32, sa: f64, sx: &[f64], incx: i32, sy: &mut [f64], incy: i32) {
    if n <= 0 || sa == 0.0 {
        return;
    }
    if incx != incy || incx < 0 || incy < 0 {
        let (mut ix, mut iy) = (start(n, incx), start(n, incy));
        for _ in 0..n {
            sy[iy as usize] += sa * sx[ix as usize];
            ix = ix.wrapping_add(incx);
            iy = iy.wrapping_add(incy);
        }
    } else if incx != 1 {
        let mut j = 0i32;
        for _ in 0..n {
            sy[j as usize] += sa * sx[j as usize];
            j = j.wrapping_add(incx);
        }
    } else {
        let m = n % 4;
        for i in 0..m as usize {
            sy[i] += sa * sx[i];
        }
        for i in (m as usize..n as usize).step_by(4) {
            sy[i] += sa * sx[i];
            sy[i + 1] += sa * sx[i + 1];
            sy[i + 2] += sa * sx[i + 2];
            sy[i + 3] += sa * sx[i + 3];
        }
    }
}

/// C: `void dcopy (int n, double *sx, int incx, double *sy, int incy);`
pub fn dcopy(n: i32, sx: &[f64], incx: i32, sy: &mut [f64], incy: i32) {
    if n <= 0 {
        return;
    }
    if incx != incy || incx < 0 || incy < 0 {
        let (mut ix, mut iy) = (start(n, incx), start(n, incy));
        for _ in 0..n {
            sy[iy as usize] = sx[ix as usize];
            ix = ix.wrapping_add(incx);
            iy = iy.wrapping_add(incy);
        }
    } else if incx != 1 {
        let mut j = 0i32;
        for _ in 0..n {
            sy[j as usize] = sx[j as usize];
            j = j.wrapping_add(incx);
        }
    } else {
        let m = n % 7;
        for i in 0..m as usize {
            sy[i] = sx[i];
        }
        for i in (m as usize..n as usize).step_by(7) {
            sy[i] = sx[i];
            sy[i + 1] = sx[i + 1];
            sy[i + 2] = sx[i + 2];
            sy[i + 3] = sx[i + 3];
            sy[i + 4] = sx[i + 4];
            sy[i + 5] = sx[i + 5];
            sy[i + 6] = sx[i + 6];
        }
    }
}

/// C: `double ddot (int n, double *sx, int incx, double *sy, int incy);`
pub fn ddot(n: i32, sx: &[f64], incx: i32, sy: &[f64], incy: i32) -> f64 {
    let mut sum = 0.0;
    if n <= 0 {
        return sum;
    }
    if incx != incy || incx < 0 || incy < 0 {
        let (mut ix, mut iy) = (start(n, incx), start(n, incy));
        for _ in 0..n {
            sum += sx[ix as usize] * sy[iy as usize];
            ix = ix.wrapping_add(incx);
            iy = iy.wrapping_add(incy);
        }
    } else if incx != 1 {
        let mut j = 0i32;
        for _ in 0..n {
            sum += sx[j as usize] * sy[j as usize];
            j = j.wrapping_add(incx);
        }
    } else {
        let m = n % 5;
        for i in 0..m as usize {
            sum += sx[i] * sy[i];
        }
        for i in (m as usize..n as usize).step_by(5) {
            sum += sx[i] * sy[i]
                + sx[i + 1] * sy[i + 1]
                + sx[i + 2] * sy[i + 2]
                + sx[i + 3] * sy[i + 3]
                + sx[i + 4] * sy[i + 4];
        }
    }
    sum
}

/// C: `double dnrm2 (int n, double *sx, int incx);`
/// As in C, the unscaled sum of squares can overflow or underflow.
pub fn dnrm2(n: i32, sx: &[f64], incx: i32) -> f64 {
    let mut sum = 0.0;
    if n <= 0 {
        return sum;
    }
    if incx != 1 {
        let mut j = start(n, incx);
        for _ in 0..n {
            sum += sx[j as usize] * sx[j as usize];
            j = j.wrapping_add(incx);
        }
    } else {
        let m = n % 6;
        for i in 0..m as usize {
            sum += sx[i] * sx[i];
        }
        for i in (m as usize..n as usize).step_by(6) {
            sum += sx[i] * sx[i]
                + sx[i + 1] * sx[i + 1]
                + sx[i + 2] * sx[i + 2]
                + sx[i + 3] * sx[i + 3]
                + sx[i + 4] * sx[i + 4]
                + sx[i + 5] * sx[i + 5];
        }
    }
    sum.sqrt()
}

/// C: `void dscal (int n, double sa, double *sx, int incx);`
pub fn dscal(n: i32, sa: f64, sx: &mut [f64], incx: i32) {
    if n <= 0 {
        return;
    }
    if incx != 1 {
        let mut j = start(n, incx);
        for _ in 0..n {
            sx[j as usize] *= sa;
            j = j.wrapping_add(incx);
        }
    } else {
        let m = n % 5;
        for i in 0..m as usize {
            sx[i] *= sa;
        }
        for i in (m as usize..n as usize).step_by(5) {
            sx[i] *= sa;
            sx[i + 1] *= sa;
            sx[i + 2] *= sa;
            sx[i + 3] *= sa;
            sx[i + 4] *= sa;
        }
    }
}

/// C: `void dswap (int n, double *sx, int incx, double *sy, int incy);`
pub fn dswap(n: i32, sx: &mut [f64], incx: i32, sy: &mut [f64], incy: i32) {
    if n <= 0 {
        return;
    }
    if incx != incy || incx < 0 || incy < 0 {
        let (mut ix, mut iy) = (start(n, incx), start(n, incy));
        for _ in 0..n {
            let stemp1 = sx[ix as usize];
            sx[ix as usize] = sy[iy as usize];
            sy[iy as usize] = stemp1;
            ix = ix.wrapping_add(incx);
            iy = iy.wrapping_add(incy);
        }
    } else if incx != 1 {
        let mut j = 0i32;
        for _ in 0..n {
            let stemp1 = sx[j as usize];
            sx[j as usize] = sy[j as usize];
            sy[j as usize] = stemp1;
            j = j.wrapping_add(incx);
        }
    } else {
        let m = n % 3;
        for i in 0..m as usize {
            let stemp1 = sx[i];
            sx[i] = sy[i];
            sy[i] = stemp1;
        }
        for i in (m as usize..n as usize).step_by(3) {
            let stemp1 = sx[i];
            let stemp2 = sx[i + 1];
            let stemp3 = sx[i + 2];
            sx[i] = sy[i];
            sx[i + 1] = sy[i + 1];
            sx[i + 2] = sy[i + 2];
            sy[i] = stemp1;
            sy[i + 1] = stemp2;
            sy[i + 2] = stemp3;
        }
    }
}

/// C: `int idamax (int n, double *sx, int incx);`
/// Returns the zero-based index of the first maximum. C reads the first element
/// even when `n <= 0`, so the slice must then still contain one element.
/// A NaN at the starting element prevents subsequent maximum updates.
pub fn idamax(n: i32, sx: &[f64], incx: i32) -> i32 {
    let mut imax = 0;
    if incx != 1 {
        let mut j = start(n, incx);
        let mut smax = abs(sx[j as usize]);
        for i in 0..n {
            let xmag = abs(sx[j as usize]);
            if xmag > smax {
                imax = i;
                smax = xmag;
            }
            j = j.wrapping_add(incx);
        }
    } else {
        let mut smax = abs(sx[0]);
        for i in 1..n {
            let xmag = abs(sx[i as usize]);
            if xmag > smax {
                imax = i;
                smax = xmag;
            }
        }
    }
    imax
}
