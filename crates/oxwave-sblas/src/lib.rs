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

//! Port of cwp/lib/sblas.c. Credits: Dave Hale, Colorado School of Mines, 10/01/89.
//!
//! Single-vector routines with negative increments start at index
//! `(1-n)*incx`, corresponding to a C pointer positioned there. Two-vector
//! routines position negative increments internally. Slices must cover every
//! access; `isamax` reads one element even when n is nonpositive.

fn start(n: i32, inc: i32) -> i32 {
    if inc < 0 {
        (1 - n) * inc
    } else {
        0
    }
}
fn abs(x: f32) -> f32 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

/// C prototype: `int isamax (int n, float *sx, int incx);`
pub fn isamax(n: i32, sx: &[f32], incx: i32) -> i32 {
    let base = if n > 0 { start(n, incx) } else { 0 };
    let mut imax = 0;
    let mut smax = abs(sx[base as usize]);
    let first = if incx == 1 { 1 } else { 0 };
    for i in first..n {
        let xmag = abs(sx[(base + i * incx) as usize]);
        if xmag > smax {
            imax = i;
            smax = xmag;
        }
    }
    imax
}

/// C prototype: `float sasum (int n, float *sx, int incx);`
pub fn sasum(n: i32, sx: &[f32], incx: i32) -> f32 {
    let mut sum = 0.0;
    if n <= 0 {
        return sum;
    }
    if incx != 1 {
        let base = start(n, incx);
        for i in 0..n {
            sum += abs(sx[(base + i * incx) as usize]);
        }
    } else {
        let m = n % 6;
        for i in 0..m {
            sum += abs(sx[i as usize]);
        }
        for i in (m..n).step_by(6) {
            let i = i as usize;
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

/// C prototype: `void saxpy (int n, float sa, float *sx, int incx, float *sy, int incy);`
pub fn saxpy(n: i32, sa: f32, sx: &[f32], incx: i32, sy: &mut [f32], incy: i32) {
    if n <= 0 || sa == 0.0 {
        return;
    }
    let (mut ix, mut iy) = (start(n, incx), start(n, incy));
    for _ in 0..n {
        sy[iy as usize] += sa * sx[ix as usize];
        ix += incx;
        iy += incy;
    }
}

/// C prototype: `void scopy (int n, float *sx, int incx, float *sy, int incy);`
pub fn scopy(n: i32, sx: &[f32], incx: i32, sy: &mut [f32], incy: i32) {
    if n <= 0 {
        return;
    }
    let (mut ix, mut iy) = (start(n, incx), start(n, incy));
    for _ in 0..n {
        sy[iy as usize] = sx[ix as usize];
        ix += incx;
        iy += incy;
    }
}

/// C prototype: `float sdot (int n, float *sx, int incx, float *sy, int incy);`
pub fn sdot(n: i32, sx: &[f32], incx: i32, sy: &[f32], incy: i32) -> f32 {
    let mut sum = 0.0;
    if n <= 0 {
        return sum;
    }
    if incx != 1 || incy != 1 {
        let (mut ix, mut iy) = (start(n, incx), start(n, incy));
        for _ in 0..n {
            sum += sx[ix as usize] * sy[iy as usize];
            ix += incx;
            iy += incy;
        }
    } else {
        let m = n % 5;
        for i in 0..m {
            let i = i as usize;
            sum += sx[i] * sy[i];
        }
        for i in (m..n).step_by(5) {
            let i = i as usize;
            sum += sx[i] * sy[i]
                + sx[i + 1] * sy[i + 1]
                + sx[i + 2] * sy[i + 2]
                + sx[i + 3] * sy[i + 3]
                + sx[i + 4] * sy[i + 4];
        }
    }
    sum
}

/// C prototype: `float snrm2 (int n, float *sx, int incx);`
pub fn snrm2(n: i32, sx: &[f32], incx: i32) -> f32 {
    let mut sum = 0.0f32;
    if n <= 0 {
        return sum;
    }
    if incx != 1 {
        let base = start(n, incx);
        for i in 0..n {
            let x = sx[(base + i * incx) as usize];
            sum += x * x;
        }
    } else {
        let m = n % 6;
        for i in 0..m {
            let x = sx[i as usize];
            sum += x * x;
        }
        for i in (m..n).step_by(6) {
            let i = i as usize;
            sum += sx[i] * sx[i]
                + sx[i + 1] * sx[i + 1]
                + sx[i + 2] * sx[i + 2]
                + sx[i + 3] * sx[i + 3]
                + sx[i + 4] * sx[i + 4]
                + sx[i + 5] * sx[i + 5];
        }
    }
    (sum as f64).sqrt() as f32
}

/// C prototype: `void sscal (int n, float sa, float *sx, int incx);`
pub fn sscal(n: i32, sa: f32, sx: &mut [f32], incx: i32) {
    if n <= 0 {
        return;
    }
    let base = start(n, incx);
    for i in 0..n {
        sx[(base + i * incx) as usize] *= sa;
    }
}

/// C prototype: `void sswap (int n, float *sx, int incx, float *sy, int incy);`
pub fn sswap(n: i32, sx: &mut [f32], incx: i32, sy: &mut [f32], incy: i32) {
    if n <= 0 {
        return;
    }
    let (mut ix, mut iy) = (start(n, incx), start(n, incy));
    for _ in 0..n {
        std::mem::swap(&mut sx[ix as usize], &mut sy[iy as usize]);
        ix += incx;
        iy += incy;
    }
}
