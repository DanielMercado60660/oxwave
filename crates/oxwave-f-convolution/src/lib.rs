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

//! Port of cwp/lib/convolution.c.
//! Author: Dave Hale, Colorado School of Mines, 11/23/91.

/// Convolve sequences with first sample indices `ifx`, `ify`, and `ifz`.
///
/// C prototype: `void convolve_cwp(int lx, int ifx, float *x, int ly,
/// int ify, float *y, int lz, int ifz, float *z);`
/// Slice lengths replace the C array lengths.
pub fn convolve_cwp<'a>(
    mut ifx: i32,
    mut x: &'a [f32],
    mut ify: i32,
    mut y: &'a [f32],
    ifz: i32,
    z: &mut [f32],
) {
    if x.len() > y.len() {
        std::mem::swap(&mut x, &mut y);
        std::mem::swap(&mut ifx, &mut ify);
    }
    if z.is_empty() {
        return;
    }
    if x.is_empty() {
        z.fill(0.0);
        return;
    }
    let (ifx, ify, ifz) = (ifx as i64, ify as i64, ifz as i64);
    let (ilx, ily, ilz) = (
        ifx + x.len() as i64 - 1,
        ify + y.len() as i64 - 1,
        ifz + z.len() as i64 - 1,
    );
    let xv = |j: i64| x[(j - ifx) as usize];
    let yv = |j: i64| y[(j - ify) as usize];
    let scalar = |i: i64, lo: i64, hi: i64| {
        let mut s = 0.0f32;
        for j in lo..=hi {
            s += xv(j) * yv(i - j);
        }
        s
    };
    if x.len() <= 30 {
        let midlo = (ify + ilx).max(ifz);
        let midhi = (ily + ifx).min(ilz);
        for (k, out) in z.iter_mut().enumerate() {
            let i = ifz + k as i64;
            if i < ifx + ify || i > ilx + ily {
                *out = 0.0;
            } else if i >= midlo && i <= midhi && (i < midhi || (midhi - midlo) % 2 == 1) {
                let mut s = xv(ifx) * yv(i - ifx);
                for j in ifx + 1..=ilx {
                    s += xv(j) * yv(i - j);
                }
                *out = s;
            } else {
                *out = scalar(i, ifx.max(i - ily), ilx.min(i - ify));
            }
        }
        return;
    }
    let mut i = ifz;
    while i <= (ify + ifx - 1).min(ilz) {
        z[(i - ifz) as usize] = 0.0;
        i += 1;
    }
    if i < ify + ilx && i <= ilz && (i - ify - ifx) % 2 != 0 {
        z[(i - ifz) as usize] = scalar(i, ifx, i - ify);
        i += 1;
    }
    let high = (ilx + ify - 1).min(ilz);
    let mut jhigh = i - ify;
    while i < high {
        let (mut sa, mut sb) = (0.0f32, 0.0f32);
        let (mut xb, mut yb, mut j) = (xv(jhigh + 1), 0.0f32, jhigh);
        while j >= ifx {
            sa += xb * yb;
            let ya = yv(i - j);
            sb += xb * ya;
            let xa = xv(j);
            sa += xa * ya;
            yb = yv(i + 1 - j);
            sb += xa * yb;
            // The final C preload lies outside x and is unused.
            if j - 1 >= ifx {
                xb = xv(j - 1);
            }
            j -= 2;
        }
        z[(i - ifz) as usize] = sa;
        z[(i + 1 - ifz) as usize] = sb;
        i += 2;
        jhigh += 2;
    }
    if i == high {
        z[(i - ifz) as usize] = scalar(i, ifx, i - ify);
        i += 1;
    }
    let high = (ily + ifx).min(ilz);
    while i < high {
        let (mut sa, mut sb) = (0.0f32, 0.0f32);
        let (mut yb, mut xa, mut j) = (yv(i + 1 - ifx), xv(ifx), ifx);
        while j < ilx {
            sb += xa * yb;
            let ya = yv(i - j);
            sa += xa * ya;
            let xb = xv(j + 1);
            sb += xb * ya;
            yb = yv(i - 1 - j);
            sa += xb * yb;
            // C sometimes preloads one element beyond x, without using it.
            if j + 2 <= ilx {
                xa = xv(j + 2);
            }
            j += 2;
        }
        if (ilx - ifx) % 2 == 0 {
            sa += xv(ilx) * yv(i - ilx);
            sb += xv(ilx) * yv(i + 1 - ilx);
        }
        z[(i - ifz) as usize] = sa;
        z[(i + 1 - ifz) as usize] = sb;
        i += 2;
    }
    if i == high {
        z[(i - ifz) as usize] = scalar(i, ifx, ilx);
        i += 1;
    }
    if i <= ily + ilx && i <= ilz && (ilx - (i - ily)) % 2 == 0 {
        z[(i - ifz) as usize] = scalar(i, i - ily, ilx);
        i += 1;
    }
    let high = (ily + ilx).min(ilz);
    while i < high {
        let (mut sa, mut sb) = (0.0f32, 0.0f32);
        let mut j = i - ily;
        let mut xa = xv(j);
        let mut yb = 0.0f32;
        while j < ilx - 2 {
            sb += xa * yb;
            let ya = yv(i - j);
            sa += xa * ya;
            let xb = xv(j + 1);
            sb += xb * ya;
            yb = yv(i - 1 - j);
            sa += xb * yb;
            xa = xv(j + 2);
            j += 2;
        }
        sb += xa * yb;
        let ya = yv(i - j);
        sa += xa * ya;
        let xb = xv(j + 1);
        sb += xb * ya;
        yb = yv(i - 1 - j);
        sa += xb * yb;
        z[(i - ifz) as usize] = sa;
        z[(i + 1 - ifz) as usize] = sb;
        i += 2;
    }
    if i == high {
        z[(i - ifz) as usize] = scalar(i, i - ily, ilx);
        i += 1;
    }
    while i <= ilz {
        z[(i - ifz) as usize] = 0.0;
        i += 1;
    }
}
