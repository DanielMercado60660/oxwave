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

//! Differential tests for `par/lib/fractal.c`.
//! Author: Balazs Nemeth, Potash Corporation of Saskatoon Saskatchewan;
//! code donated to CWP c. 2008.

use oxwave_cref::XorShift;
use oxwave_fractal as rust;

extern "C" {
    fn distance(x0: f32, y0: f32, x1: f32, y1: f32) -> f32;
    fn hausdorff_dimension(ar: *mut f32, n: i32, minl: i32, maxl: i32, dl: i32) -> f32;
}

fn same(actual: f32, expected: f32, context: &str) {
    if expected.is_nan() {
        assert!(
            actual.is_nan(),
            "{context}: Rust {actual:?}, C {expected:?}"
        );
    } else if !expected.is_finite() || expected == 0.0 {
        assert_eq!(actual.to_bits(), expected.to_bits(), "{context}");
    } else {
        let error = (actual as f64 - expected as f64).abs();
        let scale = (actual as f64).abs().max((expected as f64).abs());
        assert!(
            error <= 1e-5 * scale,
            "{context}: Rust {actual:?}, C {expected:?}"
        );
    }
}

#[test]
fn distance_matches_c() {
    let mut rng = XorShift::new(0xfacade);
    let values = [0.0, -0.0, 1.0, -1.0, 1e-30, -1e-30, 1e20, -1e20];
    for i in 0..12000 {
        let mut p = [0.0_f32; 4];
        for v in &mut p {
            *v = if i < 512 {
                values[(rng.next_u64() as usize) % values.len()]
            } else {
                rng.uniform(-10000.0, 10000.0) as f32
            };
        }
        let expected = unsafe { distance(p[0], p[1], p[2], p[3]) };
        same(
            rust::distance(p[0], p[1], p[2], p[3]),
            expected,
            &format!("points {p:?}"),
        );
    }
}

#[test]
fn hausdorff_dimension_matches_c() {
    let mut rng = XorShift::new(0xdecafbad);
    let sizes = [1, 2, 3, 7, 16, 31, 64, 127, 256, 1024];
    for case in 0..4000 {
        let n = sizes[(case / 5) % sizes.len()];
        let mut ar = vec![0.0_f32; n];
        for (i, v) in ar.iter_mut().enumerate() {
            *v = match case % 5 {
                0 => 0.0,
                1 => -3.0,
                2 => i as f32 * 0.25 - 10.0,
                3 => rng.uniform(-0.001, 0.001) as f32,
                _ => rng.uniform(-1000.0, 1000.0) as f32,
            };
        }
        let minl = if case % 2 == 0 {
            1 + (rng.next_u64() % 3) as i32
        } else {
            1 + (rng.next_u64() % (n as u64 + 2)) as i32
        };
        let dl = 1 + (rng.next_u64() % 4) as i32;
        let nm = (rng.next_u64() % 12) as i32;
        let maxl = minl + dl * nm + (rng.next_u64() % dl as u64) as i32;
        let expected = unsafe { hausdorff_dimension(ar.as_mut_ptr(), n as i32, minl, maxl, dl) };
        let actual = rust::hausdorff_dimension(&ar, minl, maxl, dl);
        same(
            actual,
            expected,
            &format!("case {case}, n={n}, minl={minl}, maxl={maxl}, dl={dl}"),
        );
    }
}
