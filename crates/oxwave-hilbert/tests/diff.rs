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

//! Differential tests for cwp/lib/hilbert.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

use oxwave_cref::XorShift;

extern "C" {
    fn hilbert(n: i32, x: *const f32, y: *mut f32);
}

#[test]
fn hilbert_matches_c() {
    let mut rng = XorShift::new(0x4961_82a5);
    let sizes = [1, 2, 3, 7, 29, 30, 31, 32, 59, 60, 61, 62, 63, 64, 97, 127, 128, 257, 1024, 4096];
    for case in 0..2400 {
        let n = sizes[case % sizes.len()];
        let scale = [0.0, 1e-15, 1e-3, 1.0, 1e6][(case / sizes.len()) % 5];
        let mut x: Vec<f32> = (0..n).map(|_| rng.uniform_f32(-scale, scale)).collect();
        if case % 11 == 0 { x.fill(-scale); }
        if case % 13 == 0 { x.fill(0.0); x[case % n] = scale; }
        // The optimized C convolution speculatively loads a sample immediately
        // outside the logical input, even though that value is not used. Keep
        // those loads in allocated memory, including for tiny signals.
        const PAD: usize = 64;
        let mut padded_x = vec![0.0f32; n + 2 * PAD];
        padded_x[PAD..PAD + n].copy_from_slice(&x);
        let mut padded_y = vec![f32::NAN; n + 2 * PAD];
        let mut actual = vec![f32::NAN; n];
        unsafe { hilbert(n as i32, padded_x.as_ptr().add(PAD), padded_y.as_mut_ptr().add(PAD)); }
        assert!(padded_y[..PAD].iter().chain(&padded_y[PAD + n..]).all(|v| v.is_nan()));
        let expected = &padded_y[PAD..PAD + n];
        oxwave_hilbert::hilbert(&x, &mut actual);
        let magnitude = expected.iter().fold(0.0f32, |m, &v| m.max(v.abs()));
        for i in 0..n {
            assert!((actual[i] - expected[i]).abs() <= 1e-5 * magnitude,
                "case {case}, n {n}, index {i}: rust {} C {}, max {magnitude}", actual[i], expected[i]);
        }
    }
}
