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

//! Port of `cwp/lib/cubicspline.c`.
//! Author: Dave Hale, Colorado School of Mines c. 1989, 1990, 1991.
//! Modified: Dave Hale, Colorado School of Mines, 02/28/91: changed to work for n=1.
//! Modified: Dave Hale, Colorado School of Mines, 08/04/91: fixed bug in computation
//! of left end derivative (cmonot and csplin).

fn setup(x: &[f32], y: &[f32], yd: &mut [[f32; 4]]) -> bool {
    let n = x.len();
    assert!(n > 0 && y.len() >= n && yd.len() >= n);
    for i in 0..n {
        yd[i][0] = y[i];
    }
    if n == 1 {
        yd[0][1..].fill(0.0);
        return true;
    }
    if n == 2 {
        let s = (y[1] - y[0]) / (x[1] - x[0]);
        yd[0][1] = s;
        yd[1][1] = s;
        yd[0][2..].fill(0.0);
        yd[1][2..].fill(0.0);
        return true;
    }
    false
}
fn derivatives(x: &[f32], y: &[f32], yd: &mut [[f32; 4]]) {
    let n = x.len();
    for i in 0..n - 1 {
        let h = x[i + 1] - x[i];
        let d = (y[i + 1] - y[i]) / h;
        // C adds the float slopes before the double literal promotes the result.
        let q = ((yd[i][1] + yd[i + 1][1]) as f64 - 2.0 * d as f64) as f32;
        yd[i][2] = (2.0 * (d - yd[i][1] - q) as f64 / h as f64) as f32;
        yd[i][3] = ((q / h) as f64 * (6.0 / h as f64)) as f32;
    }
    yd[n - 1][2] = yd[n - 2][2] + (x[n - 1] - x[n - 2]) * yd[n - 2][3];
    yd[n - 1][3] = yd[n - 2][3];
}
/// Akima interpolation. `x` must be strictly monotonic and nonempty; `y` and
/// `yd` must have at least `x.len()` elements. Extra elements are untouched.
/// Rows contain `[y, y', y'', y''']`, without factorial scaling. For an interval
/// offset `h`, evaluate `yd[i][0] + h*(yd[i][1] + h*(yd[i][2]/2 + h*yd[i][3]/6))`.
/// One point gives constant interpolation; two points give linear interpolation.
///
/// C prototype: `void cakima(int n, float x[], float y[], float yd[][4]);`
pub fn cakima(x: &[f32], y: &[f32], yd: &mut [[f32; 4]]) {
    if setup(x, y, yd) {
        return;
    }
    let n = x.len();
    for i in 1..n {
        yd[i][2] = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    }
    for i in 1..n - 1 {
        yd[i][3] = (yd[i + 1][2] - yd[i][2]).abs();
    }
    yd[0][3] = yd[1][3];
    yd[n - 1][3] = yd[n - 2][3];
    let mut w = yd[1][3] + yd[0][3];
    let a = (2.0 * yd[1][2] as f64 - yd[2][2] as f64) as f32;
    yd[0][1] = if w != 0.0 {
        (yd[1][3] * a + yd[0][3] * yd[1][2]) / w
    } else {
        (0.5 * (a + yd[1][2]) as f64) as f32
    };
    for i in 1..n - 1 {
        w = yd[i + 1][3] + yd[i - 1][3];
        yd[i][1] = if w != 0.0 {
            (yd[i + 1][3] * yd[i][2] + yd[i - 1][3] * yd[i + 1][2]) / w
        } else {
            (0.5 * (yd[i][2] + yd[i + 1][2]) as f64) as f32
        };
    }
    w = yd[n - 2][3] + yd[n - 1][3];
    let a = (2.0 * yd[n - 1][2] as f64 - yd[n - 2][2] as f64) as f32;
    yd[n - 1][1] = if w != 0.0 {
        (yd[n - 2][3] * a + yd[n - 1][3] * yd[n - 1][2]) / w
    } else {
        (0.5 * (yd[n - 1][2] + a) as f64) as f32
    };
    for i in 1..n {
        let h = x[i] - x[i - 1];
        let q = ((yd[i - 1][1] + yd[i][1]) as f64 - 2.0 * yd[i][2] as f64) as f32;
        yd[i - 1][2] = (2.0 * (yd[i][2] - yd[i - 1][1] - q) as f64 / h as f64) as f32;
        yd[i - 1][3] = ((q / h) as f64 * (6.0 / h as f64)) as f32;
    }
    yd[n - 1][2] = yd[n - 2][2] + (x[n - 1] - x[n - 2]) * yd[n - 2][3];
    yd[n - 1][3] = yd[n - 2][3];
}
fn endpoint(h1: f32, h2: f32, d1: f32, d2: f32, left: bool) -> f32 {
    let h = h1 + h2;
    let (w1, w2) = if left {
        ((h1 + h) / h, -h1 / h)
    } else {
        (-h2 / h, (h2 + h) / h)
    };
    let mut s = w1 * d1 + w2 * d2;
    let d = if left { d1 } else { d2 };
    if s * d <= 0.0 {
        s = 0.0;
    } else if d1 * d2 < 0.0 {
        let m = (3.0 * d as f64) as f32;
        if s.abs() > m.abs() {
            s = m;
        }
    }
    s
}
/// Shape-preserving Fritsch–Carlson interpolation. See [`cakima`] for inputs.
///
/// C prototype: `void cmonot(int n, float x[], float y[], float yd[][4]);`
pub fn cmonot(x: &[f32], y: &[f32], yd: &mut [[f32; 4]]) {
    if setup(x, y, yd) {
        return;
    }
    let n = x.len();
    let h1 = x[1] - x[0];
    let h2 = x[2] - x[1];
    yd[0][1] = endpoint(h1, h2, (y[1] - y[0]) / h1, (y[2] - y[1]) / h2, true);
    for i in 1..n - 1 {
        let h1 = x[i] - x[i - 1];
        let h2 = x[i + 1] - x[i];
        let h = h1 + h2;
        let d1 = (y[i] - y[i - 1]) / h1;
        let d2 = (y[i + 1] - y[i]) / h2;
        yd[i][1] = if d1 * d2 <= 0.0 {
            0.0
        } else {
            let h3 = h + h + h;
            let w1 = (h + h1) / h3;
            let w2 = (h + h2) / h3;
            let min = if d1.abs() < d2.abs() {
                d1.abs()
            } else {
                d2.abs()
            };
            let max = if d1.abs() > d2.abs() {
                d1.abs()
            } else {
                d2.abs()
            };
            let r1 = d1 / max;
            let r2 = d2 / max;
            min / (w1 * r1 + w2 * r2)
        };
    }
    let h1 = x[n - 2] - x[n - 3];
    let h2 = x[n - 1] - x[n - 2];
    yd[n - 1][1] = endpoint(
        h1,
        h2,
        (y[n - 2] - y[n - 3]) / h1,
        (y[n - 1] - y[n - 2]) / h2,
        false,
    );
    derivatives(x, y, yd);
}
/// Cubic spline with continuous second derivatives. See [`cakima`] for inputs.
///
/// C prototype: `void csplin(int n, float x[], float y[], float yd[][4]);`
pub fn csplin(x: &[f32], y: &[f32], yd: &mut [[f32; 4]]) {
    if setup(x, y, yd) {
        return;
    }
    let n = x.len();
    let h1 = x[1] - x[0];
    let h2 = x[2] - x[1];
    let left = endpoint(h1, h2, (y[1] - y[0]) / h1, (y[2] - y[1]) / h2, true);
    let h1 = x[n - 2] - x[n - 3];
    let h2 = x[n - 1] - x[n - 2];
    let right = endpoint(
        h1,
        h2,
        (y[n - 2] - y[n - 3]) / h1,
        (y[n - 1] - y[n - 2]) / h2,
        false,
    );
    yd[0][0] = 1.0;
    yd[0][2] = (2.0 * left as f64) as f32;
    for i in 1..n - 1 {
        let h1 = x[i] - x[i - 1];
        let h2 = x[i + 1] - x[i];
        let d1 = (y[i] - y[i - 1]) / h1;
        let d2 = (y[i + 1] - y[i]) / h2;
        let a = h2 / (h1 + h2);
        yd[i][0] = a;
        yd[i][2] = (3.0 * ((a * d1) as f64 + (1.0 - a as f64) * d2 as f64)) as f32;
    }
    yd[n - 1][0] = 0.0;
    yd[n - 1][2] = (2.0 * right as f64) as f32;
    let mut t = 2.0f32;
    yd[0][1] = yd[0][2] / t;
    for i in 1..n {
        yd[i][3] = ((1.0 - yd[i - 1][0] as f64) / t as f64) as f32;
        t = (2.0 - (yd[i][0] * yd[i][3]) as f64) as f32;
        yd[i][1] = (yd[i][2] - yd[i][0] * yd[i - 1][1]) / t;
    }
    for i in (0..n - 1).rev() {
        yd[i][1] -= yd[i + 1][3] * yd[i + 1][1];
    }
    for i in 0..n {
        yd[i][0] = y[i];
    }
    derivatives(x, y, yd);
}
/// Hermite interpolation. See [`cakima`] for inputs. For three or more points,
/// the final row's second and third derivative slots are left unchanged.
///
/// C prototype: `void chermite(int n, float x[], float y[], float yd[][4]);`
pub fn chermite(x: &[f32], y: &[f32], yd: &mut [[f32; 4]]) {
    if setup(x, y, yd) {
        return;
    }
    let n = x.len();
    yd[0][1] = (y[2] - y[0]) / (x[2] - x[0]);
    yd[n - 1][1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
    for i in 1..n - 1 {
        yd[i][1] = (0.5 * (y[i] - y[i - 1]) as f64 / (x[i] - x[i - 1]) as f64
            + 0.5 * (y[i + 1] - y[i]) as f64 / (x[i + 1] - x[i]) as f64) as f32;
    }
    for i in 0..n - 1 {
        let dx = x[i + 1] - x[i];
        let f1 = (y[i + 1] - y[i]) / dx;
        let f3 = ((yd[i][1] + yd[i + 1][1]) as f64 - 2.0 * f1 as f64) as f32;
        yd[i][2] = (2.0 * (f1 - yd[i][1] - f3) as f64 / dx as f64) as f32;
        yd[i][3] = (6.0 * f3 as f64 / (dx * dx) as f64) as f32;
    }
}
