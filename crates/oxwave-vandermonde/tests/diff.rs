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

//! Differential tests for cwp/lib/vandermonde.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

use oxwave_cref::XorShift;

extern "C" {
    fn vanded(n: i32, v: *const f64, b: *const f64, x: *mut f64);
    fn vandef(n: i32, v: *const f32, b: *const f32, x: *mut f32);
}

fn inputs(r: &mut XorShift, n: usize, case: usize) -> (Vec<f64>, Vec<f64>) {
    let mut v = Vec::with_capacity(n);
    let mut b = Vec::with_capacity(n);
    for i in 0..n {
        let a = (r.next_u32() as f64 / u32::MAX as f64) * 2.0 - 1.0;
        let c = (r.next_u32() as f64 / u32::MAX as f64) * 2.0 - 1.0;
        v.push(match case % 5 {
            0 => (i as f64 - n as f64 / 2.0) / n.max(1) as f64,
            1 => a,
            2 => a * 0.001,
            3 => {
                if i % 3 == 0 {
                    0.0
                } else {
                    a
                }
            }
            _ => a * 4.0,
        });
        b.push(if (case + i) % 7 == 0 { 0.0 } else { c * 10.0 });
    }
    (v, b)
}

#[test]
fn vanded_matches_c() {
    let mut rng = XorShift::new(17891);
    for case in 0..3000 {
        let n = [0, 1, 2, 3, 7, 8, 17, 32, 61, 128][case % 10];
        let (v, b) = inputs(&mut rng, n, case / 10);
        let mut expected = vec![0.0; n + 1];
        let mut actual = vec![0.0; n + 1];
        expected[n] = 123.0;
        actual[n] = 123.0;
        unsafe { vanded(n as i32, v.as_ptr(), b.as_ptr(), expected.as_mut_ptr()) };
        oxwave_vandermonde::vanded(&v, &b, &mut actual);
        for i in 0..=n {
            assert!(
                actual[i].to_bits() == expected[i].to_bits()
                    || (actual[i].is_nan() && expected[i].is_nan()),
                "case {case} n {n} index {i}: {} != {}",
                actual[i],
                expected[i]
            );
        }
    }
}

#[test]
fn vandef_matches_c() {
    let mut rng = XorShift::new(21783);
    for case in 0..3000 {
        let n = [0, 1, 2, 3, 7, 8, 17, 32, 61, 128][case % 10];
        let (v, b) = inputs(&mut rng, n, case / 10);
        let v: Vec<f32> = v.into_iter().map(|x| x as f32).collect();
        let b: Vec<f32> = b.into_iter().map(|x| x as f32).collect();
        let mut expected = vec![0.0; n + 1];
        let mut actual = vec![0.0; n + 1];
        expected[n] = 123.0;
        actual[n] = 123.0;
        unsafe { vandef(n as i32, v.as_ptr(), b.as_ptr(), expected.as_mut_ptr()) };
        oxwave_vandermonde::vandef(&v, &b, &mut actual);
        for i in 0..=n {
            assert!(
                actual[i].to_bits() == expected[i].to_bits()
                    || (actual[i].is_nan() && expected[i].is_nan()),
                "case {case} n {n} index {i}: {} != {}",
                actual[i],
                expected[i]
            );
        }
    }
}
