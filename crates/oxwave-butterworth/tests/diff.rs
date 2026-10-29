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

//! Differential tests for cwp/lib/antialias.c; Dave Hale, Colorado School of Mines.

use oxwave_cref::XorShift;
extern "C" {
    fn antialias(frac: f32, phase: i32, n: i32, p: *const f32, q: *mut f32);
}

#[test]
fn antialias_matches_c() {
    let mut rng = XorShift::new(0x12345678);
    let sizes = [1, 2, 3, 7, 16, 31, 64, 127, 256, 1024];
    for case in 0..5000 {
        let n = sizes[case % sizes.len()];
        let frac = match case % 16 {
            0 => 0.0, 1 => -0.0, 2 => 1.0, 3 => -1.0,
            4 => 1.2, 5 => -1.2, 6 => 0.001, 7 => 0.999,
            _ => rng.uniform_f32(-0.999, 0.999),
        };
        let phase = [0, 1, -1, 2][(case / 16) % 4];
        let input = if case % 19 == 0 { vec![0.0; n] } else { rng.vec_f32(n, -100.0, 100.0) };
        let mut expected = vec![0.0; n];
        let mut actual = vec![0.0; n];
        unsafe { antialias(frac, phase, n as i32, input.as_ptr(), expected.as_mut_ptr()); }
        oxwave_butterworth::antialias(frac, phase, &input, &mut actual);
        let scale = expected.iter().fold(0.0f32, |m, v| m.max(v.abs()));
        for (i, (&a, &b)) in actual.iter().zip(&expected).enumerate() {
            assert!((a - b).abs() <= 1e-5 * scale || (a.is_nan() && b.is_nan()),
                "case={case} n={n} frac={frac} phase={phase} i={i}: rust={a} C={b} scale={scale}");
        }
    }
}
