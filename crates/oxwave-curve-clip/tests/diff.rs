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

//! Differential tests for cwp/lib/yclip.c.
//! Author: Dave Hale, Colorado School of Mines, 07/03/89.

use oxwave_cref::XorShift;

extern "C" {
    fn yclip(
        nx: i32,
        dx: f32,
        fx: f32,
        y: *const f32,
        ymin: f32,
        ymax: f32,
        xc: *mut f32,
        yc: *mut f32,
    ) -> i32;
}

fn compare(y: &[f32], dx: f32, fx: f32, ymin: f32, ymax: f32) {
    let mut cx = vec![-12345.0; 2 * y.len()];
    let mut cy = cx.clone();
    let mut rx = cx.clone();
    let mut ry = cx.clone();
    let cn = unsafe {
        yclip(
            y.len() as i32,
            dx,
            fx,
            y.as_ptr(),
            ymin,
            ymax,
            cx.as_mut_ptr(),
            cy.as_mut_ptr(),
        )
    };
    let rn = oxwave_curve_clip::yclip(dx, fx, y, ymin, ymax, &mut rx, &mut ry);
    assert_eq!(rn, cn, "y={y:?}, bounds={ymin:?}..{ymax:?}");
    for i in 0..cx.len() {
        assert_eq!(
            rx[i].to_bits(),
            cx[i].to_bits(),
            "x[{i}], y={y:?}, dx={dx}, fx={fx}, bounds={ymin}..{ymax}"
        );
        assert_eq!(ry[i].to_bits(), cy[i].to_bits(), "y[{i}]");
    }
}

#[test]
fn yclip_matches_c() {
    let mut rng = XorShift::new(0x912b_437a);
    for &n in &[1, 2, 3, 7, 16, 31, 64, 127, 256, 1024] {
        for case in 0..400 {
            let scale = [0.00001, 1.0, 1000.0][case % 3];
            let lo = rng.uniform_f32(-3.0, 0.0) * scale;
            let hi = if case % 11 == 0 {
                lo
            } else {
                rng.uniform_f32(0.0, 3.0) * scale
            };
            let mut y = rng.vec_f32(n, -5.0 * scale, 5.0 * scale);
            for (i, value) in y.iter_mut().enumerate() {
                match (case + i) % 13 {
                    0 => *value = 0.0,
                    1 => *value = lo,
                    2 => *value = hi,
                    _ => (),
                }
            }
            let dx = if case % 7 == 0 {
                0.0
            } else {
                rng.uniform_f32(-10.0, 10.0)
            };
            let fx = rng.uniform_f32(-100.0, 100.0);
            compare(&y, dx, fx, lo, hi);
        }
    }
    for y in [
        &[-2.0][..],
        &[2.0],
        &[-2.0, 2.0],
        &[2.0, -2.0],
        &[-2.0, -2.0, 0.0, 2.0, 2.0],
        &[0.0, -0.0, 0.0],
    ] {
        compare(y, -1.0, 0.0, -1.0, 1.0);
    }
}
