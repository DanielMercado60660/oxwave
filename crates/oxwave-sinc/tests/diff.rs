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

//! Differential tests for cwp/lib/sinc.c.
//! Author of the C source: Dave Hale, Colorado School of Mines, 06/02/89.

use oxwave_cref::XorShift;

extern "C" {
    fn dsinc(x: f64) -> f64;
    fn fsinc(x: f32) -> f32;
}

fn compare(actual: f64, expected: f64, input: f64) {
    if actual.is_nan() || expected.is_nan() {
        assert!(actual.is_nan() && expected.is_nan(), "x={input:?}, Rust={actual:?}, C={expected:?}");
    } else if actual != expected {
        let scale = actual.abs().max(expected.abs());
        assert!(actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= 1e-5 * scale,
            "x={input:?}, Rust={actual:?}, C={expected:?}");
    }
}

#[test]
fn dsinc_matches_c() {
    let mut rng = XorShift::new(0x51ac);
    let edges = [0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 3.0, 7.0, 16.0,
        f64::MIN_POSITIVE, -f64::MIN_POSITIVE, f64::from_bits(1),
        -f64::from_bits(1), f64::MAX, -f64::MAX, f64::INFINITY,
        f64::NEG_INFINITY, f64::NAN];
    for x in edges {
        compare(oxwave_sinc::dsinc(x), unsafe { dsinc(x) }, x);
    }
    for _ in 0..50_000 {
        for x in [rng.uniform(-1000.0, 1000.0), rng.uniform(-1e-8, 1e-8),
            rng.uniform(-1e12, 1e12), f64::from_bits(rng.next_u64())] {
            compare(oxwave_sinc::dsinc(x), unsafe { dsinc(x) }, x);
        }
    }
}

#[test]
fn fsinc_matches_c() {
    let mut rng = XorShift::new(0xf51ac);
    let edges = [0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 3.0, 7.0, 16.0,
        f32::MIN_POSITIVE, -f32::MIN_POSITIVE, f32::from_bits(1),
        -f32::from_bits(1), f32::MAX, -f32::MAX, f32::INFINITY,
        f32::NEG_INFINITY, f32::NAN];
    for x in edges {
        compare(oxwave_sinc::fsinc(x) as f64, unsafe { fsinc(x) } as f64, x as f64);
    }
    for _ in 0..50_000 {
        for x in [rng.uniform_f32(-1000.0, 1000.0), rng.uniform_f32(-1e-8, 1e-8),
            rng.uniform_f32(-1e12, 1e12), rng.any_f32()] {
            compare(oxwave_sinc::fsinc(x) as f64, unsafe { fsinc(x) } as f64, x as f64);
        }
    }
}
