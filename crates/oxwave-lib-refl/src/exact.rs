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

//! Exact float intermediates from cwp/lib/complex.c, required by the
//! cancellation-sensitive Filon integrals in Refl/lib/psvreflect.c and shreflect.c.
//! Credits: Copyright (c) Colorado School of Mines, 2011.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.
//! Modified: Dave Hale, Colorado School of Mines, 04/26/90; added cinv().
//! Reference: Adapted from Press et al., 1988, Numerical Recipes in C, Appendix E.
use oxwave_cwp::complex::{cmplx, Complex};

pub(super) fn cwp_csqrt(z: Complex) -> Complex {
    if z.r == 0.0 && z.i == 0.0 {
        return cmplx(0.0, 0.0);
    }
    let x = (z.r as f64).abs() as f32;
    let y = (z.i as f64).abs() as f32;
    let w = if x >= y {
        let r = y / x;
        ((x as f64).sqrt() * (0.5 * (1.0 + (1.0 + (r * r) as f64).sqrt())).sqrt()) as f32
    } else {
        let r = x / y;
        ((y as f64).sqrt() * (0.5 * (r as f64 + (1.0 + (r * r) as f64).sqrt())).sqrt()) as f32
    };
    if z.r >= 0.0 {
        cmplx(w, (z.i as f64 / (2.0 * w as f64)) as f32)
    } else {
        let i = if z.i >= 0.0 { w } else { -w };
        cmplx((z.i as f64 / (2.0 * i as f64)) as f32, i)
    }
}

pub(super) fn rcabs(z: Complex) -> f32 {
    let x = (z.r as f64).abs() as f32;
    let y = (z.i as f64).abs() as f32;
    if x == 0.0 {
        y
    } else if y == 0.0 {
        x
    } else if x > y {
        let temp = y / x;
        (x as f64 * (1.0 + (temp * temp) as f64).sqrt()) as f32
    } else {
        let temp = x / y;
        (y as f64 * (1.0 + (temp * temp) as f64).sqrt()) as f32
    }
}

pub(super) fn cinv(z: Complex) -> Complex {
    let s = (1.0 / (z.r * z.r + z.i * z.i) as f64) as f32;
    cmplx(z.r * s, -z.i * s)
}
