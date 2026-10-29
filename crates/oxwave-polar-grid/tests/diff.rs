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

//! Differential tests for cwp/lib/polar.c.
//! Author: Dave Hale, Colorado School of Mines, 06/15/90.

use oxwave_cref::XorShift;
use oxwave_polar_grid::{polartorect as rust_polar, recttopolar as rust_rect};

extern "C" {
    fn recttopolar(
        nx: i32,
        dx: f32,
        fx: f32,
        ny: i32,
        dy: f32,
        fy: f32,
        p: *mut *mut f32,
        na: i32,
        da: f32,
        fa: f32,
        nr: i32,
        dr: f32,
        fr: f32,
        q: *mut *mut f32,
    );
    fn polartorect(
        na: i32,
        da: f32,
        fa: f32,
        nr: i32,
        dr: f32,
        fr: f32,
        q: *mut *mut f32,
        nx: i32,
        dx: f32,
        fx: f32,
        ny: i32,
        dy: f32,
        fy: f32,
        p: *mut *mut f32,
    );
}

fn check(c: &[Vec<f32>], r: &[Vec<f32>], case: usize) {
    let scale = c.iter().flatten().fold(0.0f32, |m, v| m.max(v.abs()));
    for (i, (a, b)) in c.iter().flatten().zip(r.iter().flatten()).enumerate() {
        assert!(
            a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-5 * scale,
            "case {case}, element {i}: C {a}, Rust {b}, scale {scale}"
        );
    }
}

fn run(polar: bool) {
    let mut rng = XorShift::new(if polar { 130 } else { 131 });
    let sizes = [1usize, 2, 3, 5, 7, 8, 16, 31, 64];
    let edge_sizes = [1usize, 2, 3, 5, 7, 8, 16, 31, 64, 257, 1024];
    for case in 0..2750 {
        // Retain the original 2,500 cases and add larger primes and powers
        // of two for both source and destination dimensions.
        let sizes = if case < 2500 {
            &sizes[..]
        } else {
            &edge_sizes[..]
        };
        // C needs both source dimensions >= 2 to access four neighbors.
        let width = sizes[1 + (rng.next_u32() as usize % (sizes.len() - 1))];
        let height = sizes[1 + (rng.next_u32() as usize % (sizes.len() - 1))];
        let ow = sizes[rng.next_u32() as usize % sizes.len()];
        let oh = sizes[rng.next_u32() as usize % sizes.len()];
        let mut input: Vec<Vec<f32>> = (0..height)
            .map(|_| {
                (0..width)
                    .map(|_| {
                        if case % 17 == 0 {
                            0.0
                        } else {
                            rng.uniform_f32(-1000.0, 1000.0)
                        }
                    })
                    .collect()
            })
            .collect();
        let mut c = vec![vec![0.0; ow]; oh];
        let mut result = c.clone();
        let mut ip: Vec<*mut f32> = input.iter_mut().map(|row| row.as_mut_ptr()).collect();
        let mut op: Vec<*mut f32> = c.iter_mut().map(|row| row.as_mut_ptr()).collect();
        let src = [0.01, 0.25, 1.0, -0.1, -1.0, 3.0];
        let dst = [0.0, 0.01, 0.25, 1.0, -0.5, 3.0];
        let d1 = src[rng.next_u32() as usize % src.len()];
        let d2 = src[rng.next_u32() as usize % src.len()];
        let d3 = dst[rng.next_u32() as usize % dst.len()];
        let d4 = dst[rng.next_u32() as usize % dst.len()];
        let f1 = if case % 3 == 0 {
            0.0
        } else {
            rng.uniform_f32(-5.0, 5.0)
        };
        let f2 = if case % 4 == 0 {
            0.0
        } else {
            rng.uniform_f32(-5.0, 5.0)
        };
        let f3 = if case % 5 == 0 {
            0.0
        } else {
            rng.uniform_f32(-5.0, 5.0)
        };
        let f4 = if case % 6 == 0 {
            0.0
        } else {
            rng.uniform_f32(-5.0, 5.0)
        };
        unsafe {
            if polar {
                polartorect(
                    width as i32,
                    d1,
                    f1,
                    height as i32,
                    d2,
                    f2,
                    ip.as_mut_ptr(),
                    ow as i32,
                    d3,
                    f3,
                    oh as i32,
                    d4,
                    f4,
                    op.as_mut_ptr(),
                );
            } else {
                recttopolar(
                    width as i32,
                    d1,
                    f1,
                    height as i32,
                    d2,
                    f2,
                    ip.as_mut_ptr(),
                    ow as i32,
                    d3,
                    f3,
                    oh as i32,
                    d4,
                    f4,
                    op.as_mut_ptr(),
                );
            }
        }
        let ins: Vec<&[f32]> = input.iter().map(Vec::as_slice).collect();
        let mut outs: Vec<&mut [f32]> = result.iter_mut().map(Vec::as_mut_slice).collect();
        if polar {
            rust_polar(d1, f1, d2, f2, &ins, d3, f3, d4, f4, &mut outs);
        } else {
            rust_rect(d1, f1, d2, f2, &ins, d3, f3, d4, f4, &mut outs);
        }
        check(&c, &result, case);
    }
}

#[test]
fn polartorect_matches_c() {
    run(true);
}

#[test]
fn recttopolar_matches_c() {
    run(false);
}

#[test]
fn polartorect_single_angle_matches_c() {
    let mut rng = XorShift::new(132);
    let sizes = [1usize, 2, 3, 5, 7, 8, 31, 64, 257, 1024];
    for case in 0..2500 {
        let nr = sizes[1 + case % (sizes.len() - 1)];
        let nx = sizes[(case / (sizes.len() - 1)) % sizes.len()];
        let ny = sizes[rng.next_u32() as usize % 6];
        // With na == 1, C uses ia == -1 and sa == 1. Back each row
        // with a finite guard so its zero-weight q[ir][-1] read is valid.
        let mut input: Vec<[f32; 2]> = (0..nr)
            .map(|_| {
                let guard = rng.uniform_f32(-1000.0, 1000.0);
                let sample = match case % 19 {
                    0 => 0.0,
                    1 => -0.0,
                    _ => rng.uniform_f32(-1000.0, 1000.0),
                };
                [guard, sample]
            })
            .collect();
        let mut c = vec![vec![f32::NAN; nx]; ny];
        let mut result = c.clone();
        let mut ip: Vec<*mut f32> = input
            .iter_mut()
            .map(|row| unsafe { row.as_mut_ptr().add(1) })
            .collect();
        let mut op: Vec<*mut f32> = c.iter_mut().map(Vec::as_mut_ptr).collect();
        let steps = [0.01, 0.25, 1.0, -0.1, -1.0, 3.0];
        let da = steps[case % steps.len()];
        let dr = steps[(case / steps.len()) % steps.len()];
        let fa = rng.uniform_f32(-5.0, 5.0);
        let fr = rng.uniform_f32(-5.0, 5.0);
        let dx = if case % 7 == 0 { 0.0 } else { 0.25 };
        let dy = if case % 11 == 0 { 0.0 } else { -0.5 };
        let fx = if case % 3 == 0 { 0.0 } else { -2.0 };
        let fy = match case % 4 {
            0 => 0.0,
            1 => 2.0,
            _ => -2.0,
        };
        unsafe {
            polartorect(
                1,
                da,
                fa,
                nr as i32,
                dr,
                fr,
                ip.as_mut_ptr(),
                nx as i32,
                dx,
                fx,
                ny as i32,
                dy,
                fy,
                op.as_mut_ptr(),
            );
        }
        let ins: Vec<&[f32]> = input.iter().map(|row| &row[1..]).collect();
        let mut outs: Vec<&mut [f32]> = result.iter_mut().map(Vec::as_mut_slice).collect();
        rust_polar(da, fa, dr, fr, &ins, dx, fx, dy, fy, &mut outs);
        check(&c, &result, case);
    }
}
