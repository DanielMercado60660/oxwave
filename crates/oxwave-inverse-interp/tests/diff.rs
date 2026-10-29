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

//! Differential tests for cwp/lib/yxtoxy.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

use oxwave_cref::XorShift;

extern "C" {
    fn yxtoxy(
        nx: i32,
        dx: f32,
        fx: f32,
        y: *const f32,
        ny: i32,
        dy: f32,
        fy: f32,
        xylo: f32,
        xyhi: f32,
        x: *mut f32,
    );
}

fn compare(
    y: &[f32],
    dx: f32,
    fx: f32,
    ny: usize,
    dy: f32,
    fy: f32,
    xylo: f32,
    xyhi: f32,
) -> Vec<f32> {
    let mut c = vec![12345.0; ny];
    let mut rust = c.clone();
    unsafe {
        yxtoxy(
            y.len() as i32,
            dx,
            fx,
            y.as_ptr(),
            ny as i32,
            dy,
            fy,
            xylo,
            xyhi,
            c.as_mut_ptr(),
        )
    };
    oxwave_inverse_interp::yxtoxy(dx, fx, y, dy, fy, xylo, xyhi, &mut rust);
    for (i, (&a, &b)) in c.iter().zip(&rust).enumerate() {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "index {i}, y={y:?}, dx={dx}, fx={fx}, dy={dy}, fy={fy}, xylo={xylo}, xyhi={xyhi}"
        );
    }
    rust
}

#[test]
fn yxtoxy_endpoints_match_c() {
    // Exact knots and their adjacent floats exercise both inclusive comparisons.
    for y in [
        &[0.0][..],
        &[-0.0][..],
        &[-2.0, 0.0][..],
        &[-2.0, 0.0, 3.0][..],
    ] {
        for ny in [0, 1, 2, 3, 7, 16] {
            for dy in [0.125, 1.0, 4.0] {
                for fy in [-4.0, -2.0, -0.0, 0.0, 3.0, 4.0] {
                    compare(y, 0.5, -1.0, ny, dy, fy, -9.0, 9.0);
                }
                for &knot in y {
                    for fy in [knot.next_down(), knot, knot.next_up()] {
                        compare(y, 0.1, 0.0, ny, dy, fy, 0.0, -0.0);
                    }
                }
            }
        }
    }
}

#[test]
fn yxtoxy_single_sample_matches_c() {
    // Preserve C's special handling of equality at only the final output slot.
    for yi in [-8.0, -0.0, 0.0, 8.0] {
        for ny in [1, 2, 3, 7, 16, 127] {
            let fy = yi - (ny - 1) as f32 * 0.5;
            let last_equal = compare(&[yi], 0.25, -3.0, ny, 0.5, fy, -9.0, 9.0);
            assert_eq!(&last_equal[..ny - 1], vec![-9.0; ny - 1]);
            assert_eq!(last_equal[ny - 1], -3.0);

            let first_equal = compare(&[yi], 0.25, -3.0, ny, 0.5, yi, -9.0, 9.0);
            let expected = if ny == 1 { -3.0 } else { 9.0 };
            assert_eq!(first_equal, vec![expected; ny]);
        }
    }
}

#[test]
fn yxtoxy_accumulated_coordinates_match_c() {
    // Recomputing coordinates as origin + index * step would change these.
    for nx in [257, 4096, 8191] {
        let mut y = Vec::with_capacity(nx);
        let mut yi = -100.0;
        for i in 0..nx {
            y.push(yi);
            yi += 0.1 + (i % 7) as f32 * 0.01;
        }
        for ny in [1, 2, 257, 4096, 8191] {
            for (dx, fx) in [(0.1, -10.0), (0.1, 1e6), (1e-8, -1.0)] {
                let dy = (y[nx - 1] - y[0]) / ny as f32;
                compare(&y, dx, fx, ny, dy, y[0] - dy, -123.0, 456.0);
                compare(&y, dx, fx, ny, 0.1, y[0], -123.0, 456.0);
            }
        }
    }
}

#[test]
fn yxtoxy_random_matches_c() {
    let mut rng = XorShift::new(0x98123);
    let sizes = [1, 2, 3, 7, 16, 31, 64, 127, 256, 1021, 4096, 8191];
    let scales = [1e-12, 0.001, 1.0, 1000.0, 1e12];
    for case in 0..10_000 {
        let nx = sizes[case % sizes.len()];
        let ny = sizes[(rng.next_u32() as usize) % sizes.len()];
        let scale = scales[case % scales.len()];
        let mut y = Vec::with_capacity(nx);
        let mut value = rng.uniform_f32(-20.0, 0.0) * scale;
        for _ in 0..nx {
            y.push(value);
            value += rng.uniform_f32(0.1, 2.0) * scale;
        }
        let dx = rng.uniform_f32(0.01, 3.0) * scale;
        let fx = rng.uniform_f32(-100.0, 100.0) * scale;
        let dy = rng.uniform_f32(0.05, 3.0) * scale;
        let fy = match case % 7 {
            0 => y[0],
            1 => y[nx - 1],
            2 => y[0] - dy * 3.0,
            3 => y[nx - 1] + dy,
            4 => 0.0,
            5 => y[(rng.next_u32() as usize) % nx],
            _ => rng.uniform_f32(-30.0, nx as f32) * scale,
        };
        let xylo = rng.uniform_f32(-100.0, 100.0) * scale;
        let xyhi = rng.uniform_f32(-100.0, 100.0) * scale;
        compare(&y, dx, fx, ny, dy, fy, xylo, xyhi);
    }
}
