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

//! Differential tests for tetra/lib/tetra.c.
//! Credits: CWP: Zhaobo Meng, October 1997.

use oxwave_cref::XorShift;

#[link(name = "tetra", kind = "static")]
extern "C" {
    fn area3d(x0: *const f32, x1: *const f32, x2: *const f32) -> f32;
    fn tetra_volume(x0: *const f32, x1: *const f32, x2: *const f32, x3: *const f32) -> f32;
}

fn check(rust: f32, c: f32, points: &[[f32; 3]]) {
    if rust.is_nan() || c.is_nan() {
        assert!(rust.is_nan() && c.is_nan(), "rust={rust:?}, C={c:?}, points={points:?}");
    } else if rust.is_infinite() || c.is_infinite() {
        assert_eq!(rust, c, "points={points:?}");
    } else {
        let delta = (rust as f64 - c as f64).abs();
        let scale = (rust as f64).abs().max((c as f64).abs());
        assert!(delta <= 1e-5 * scale, "rust={rust:?}, C={c:?}, points={points:?}");
    }
}

fn cases(mut f: impl FnMut([[f32; 3]; 4])) {
    let mut rng = XorShift::new(0x1248_ace5);
    for scale in [1e-20_f32, 1e-10, 1.0, 2.0, 7.0, 1024.0, 1e10, 1e20] {
        for i in 0..4000 {
            let mut p = [[0.0; 3]; 4];
            for point in &mut p {
                for coord in point { *coord = rng.uniform_f32(-scale, scale); }
            }
            match i % 8 {
                0 => p[0] = [0.0; 3],
                1 => p[1] = p[0],
                2 => p[2] = p[1],
                3 => p[3] = p[0],
                4 => for point in &mut p { point[2] = 0.0; },
                _ => (),
            }
            f(p);
        }
    }
    for p in [
        [[0.0; 3]; 4],
        [[-0.0; 3]; 4],
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[-1.0, -2.0, -3.0], [2.0, 4.0, 6.0], [3.0, 6.0, 9.0], [4.0, 8.0, 12.0]],
    ] { f(p); }
}

#[test]
fn area3d_matches_c() {
    cases(|p| {
        let c = unsafe { area3d(p[0].as_ptr(), p[1].as_ptr(), p[2].as_ptr()) };
        check(oxwave_lib_tetra::area3d(&p[0], &p[1], &p[2]), c, &p);
    });
}

#[test]
fn tetra_volume_matches_c() {
    cases(|p| {
        let c = unsafe { tetra_volume(p[0].as_ptr(), p[1].as_ptr(), p[2].as_ptr(), p[3].as_ptr()) };
        check(oxwave_lib_tetra::tetra_volume(&p[0], &p[1], &p[2], &p[3]), c, &p);
    });
}
