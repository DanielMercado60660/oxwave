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

//! Port of cwp/lib/polar.c.
//! Author: Dave Hale, Colorado School of Mines, 06/15/90.

// C truncates the fractional index before checking the lower boundary.
fn index(v: f32, n: usize) -> (usize, f32) {
    let mut i = v as i32;
    let mut v = v;
    if i < 0 {
        i = 0;
        v = 0.0;
    }
    if i > n as i32 - 2 {
        i = n as i32 - 2;
        v = (n - 1) as f32;
    }
    (i as usize, v - i as f32)
}

fn interpolate(grid: &[&[f32]], i: usize, j: usize, s: f32, t: f32) -> f32 {
    ((1.0 - t as f64) * ((1.0 - s as f64) * grid[j][i] as f64 + (s * grid[j][i + 1]) as f64)
        + t as f64 * ((1.0 - s as f64) * grid[j + 1][i] as f64 + (s * grid[j + 1][i + 1]) as f64))
        as f32
}

/// Convert a rectangular grid `p[y][x]` to a polar grid `q[r][a]`.
///
/// C prototype: `void recttopolar(int nx, float dx, float fx, int ny, float dy, float fy, float **p, int na, float da, float fa, int nr, float dr, float fr, float **q);`
///
/// Grid dimensions come from the slices. Input rows must have equal lengths
/// and both input dimensions must be at least two. Output rows must be equal.
#[allow(clippy::too_many_arguments)]
pub fn recttopolar(
    dx: f32,
    fx: f32,
    dy: f32,
    fy: f32,
    p: &[&[f32]],
    da: f32,
    fa: f32,
    dr: f32,
    fr: f32,
    q: &mut [&mut [f32]],
) {
    assert!(p.len() >= 2 && p[0].len() >= 2 && p.iter().all(|row| row.len() == p[0].len()));
    assert!(q
        .first()
        .map_or(true, |row| q.iter().all(|r| r.len() == row.len())));
    let mut r = fr;
    for row in q.iter_mut() {
        let mut a = fa;
        for output in row.iter_mut() {
            let x = (r as f64 * (a as f64).cos()) as f32;
            let y = (r as f64 * (a as f64).sin()) as f32;
            let (ix, sx) = index((x - fx) / dx, p[0].len());
            let (iy, sy) = index((y - fy) / dy, p.len());
            *output = interpolate(p, ix, iy, sx, sy);
            a += da;
        }
        r += dr;
    }
}

/// Convert a polar grid `q[r][a]` to a rectangular grid `p[y][x]`.
///
/// C prototype: `void polartorect(int na, float da, float fa, int nr, float dr, float fr, float **q, int nx, float dx, float fx, int ny, float dy, float fy, float **p);`
///
/// Grid dimensions come from the slices. Input rows must have equal lengths
/// with at least two radii and one angle. Output rows must be equal.
/// A single angle is independent of angle, matching C's zero-weight
/// preceding-angle term without reading before the input row.
#[allow(clippy::too_many_arguments)]
pub fn polartorect(
    da: f32,
    fa: f32,
    dr: f32,
    fr: f32,
    q: &[&[f32]],
    dx: f32,
    fx: f32,
    dy: f32,
    fy: f32,
    p: &mut [&mut [f32]],
) {
    assert!(q.len() >= 2 && !q[0].is_empty() && q.iter().all(|row| row.len() == q[0].len()));
    assert!(p
        .first()
        .map_or(true, |row| p.iter().all(|r| r.len() == row.len())));
    let mut a = 0.0f32;
    let mut y = fy;
    for row in p.iter_mut() {
        let mut x = fx;
        for output in row.iter_mut() {
            if x != 0.0 {
                a = (y as f64).atan2(x as f64) as f32;
            } else if y > 0.0 {
                a = (std::f64::consts::PI / 2.0) as f32;
            } else if y < 0.0 {
                a = (-std::f64::consts::PI / 2.0) as f32;
            } else if y == 0.0 {
                a = 0.0;
            }
            let r = ((x * x + y * y) as f64).sqrt() as f32;
            let (ia, sa) = if q[0].len() == 1 {
                // C clamps ia to -1 and ai to 0, giving sa == 1.
                (0, 1.0)
            } else {
                index((a - fa) / da, q[0].len())
            };
            let (ir, sr) = index((r - fr) / dr, q.len());
            *output = if q[0].len() == 1 {
                // (1.0-sa)*q[ir][-1] contributes zero. Keep the radial
                // interpolation in double, as in C's bilinear expression.
                ((1.0 - sr as f64) * (0.0 + q[ir][0] as f64)
                    + sr as f64 * (0.0 + q[ir + 1][0] as f64)) as f32
            } else {
                interpolate(q, ia, ir, sa, sr)
            };
            x += dx;
        }
        y += dy;
    }
}
