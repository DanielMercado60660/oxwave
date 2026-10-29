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

//! Differential tests for cwp/lib/orthpolynomials.c.
//! Author: Werner M. Heigl, Apache Corporation, E&P Technology, December 2006.

use oxwave_cref::XorShift;

extern "C" {
    fn hermite_n_polynomial(h: *mut f64, h0: *mut f64, h1: *mut f64,
                            t: *const f64, nt: i32, n: i32, sigma: f64);
}

#[test]
fn hermite_n_polynomial_matches_c() {
    let mut rng = XorShift::new(0xa5be_2739);
    let sizes = [0, 1, 2, 3, 7, 31, 32, 127, 128, 257, 1024, 4096];
    let sigmas = [0.0, -0.0, 0.001, -0.001, 0.5, -1.0, 1.0, 10.0, 1e6];
    for case in 0..3000 {
        let len = sizes[case % sizes.len()];
        let n = match case % 10 { 0 => -3, 1 => 0, 2 => 1, 3 => 2, _ => rng.range_i32(3, 25) };
        let sigma = sigmas[(case / sizes.len()) % sigmas.len()];
        let scale = [0.0, 1e-6, 1.0, 10.0, 1e6][(case / 7) % 5];
        let t: Vec<f64> = (0..len).map(|_| rng.uniform(-scale, scale)).collect();
        let mut h: Vec<f64> = (0..len).map(|_| rng.uniform(-scale, scale)).collect();
        let mut h0: Vec<f64> = (0..len).map(|_| rng.uniform(-scale, scale)).collect();
        let mut h1: Vec<f64> = (0..len).map(|_| rng.uniform(-scale, scale)).collect();
        let (mut ch, mut ch0, mut ch1) = (h.clone(), h0.clone(), h1.clone());
        unsafe { hermite_n_polynomial(ch.as_mut_ptr(), ch0.as_mut_ptr(), ch1.as_mut_ptr(), t.as_ptr(), len as i32, n, sigma); }
        oxwave_orthogonal_polynomials::hermite_n_polynomial(&mut h, &mut h0, &mut h1, &t, n, sigma);
        for (name, actual, expected) in [("h", &h, &ch), ("h0", &h0, &ch0), ("h1", &h1, &ch1)] {
            for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
                assert!(a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
                        "case {case}, {name}[{i}], n {n}, sigma {sigma}: Rust {a:?}, C {b:?}");
            }
        }
    }
}
