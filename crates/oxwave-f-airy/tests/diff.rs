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

//! Differential tests against cwp/lib/airy.c.
//! Author:  Dave Hale, Colorado School of Mines, 06/06/89
//! Author:  Dave Hale, Colorado School of Mines, 07/05/89
//! Authors:  Dave Hale, Craig Artley, Colorado School of Mines, 11/13/90

use oxwave_cref::XorShift;

extern "C" {
    // Original C prototypes: float airya(float x); float airyap(float x);
    // float airyb(float x); float airybp(float x);
    fn airya(x: f32) -> f32;
    fn airyap(x: f32) -> f32;
    fn airyb(x: f32) -> f32;
    fn airybp(x: f32) -> f32;
}

fn neighbors(inputs: &mut Vec<f32>, x: f32) {
    inputs.extend([x.next_down(), x, x.next_up()]);
}

fn inputs(seed: u64) -> Vec<f32> {
    let mut rng = XorShift::new(seed);
    let mut inputs = vec![
        0.0,
        -0.0,
        f32::from_bits(1),
        -f32::from_bits(1),
        f32::MIN_POSITIVE,
        -f32::MIN_POSITIVE,
        f32::MAX,
        -f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
        f32::from_bits(0xffc0_0001),
        f32::from_bits(0x7f80_0001),
    ];
    // These are scalar APIs. Sizes exercise batches of calls, not array lengths.
    // Cover every branch with singletons, primes, powers of two and large batches.
    for size in [1, 2, 3, 7, 31, 64, 127, 256, 1024, 4093, 16_384] {
        for (lo, hi) in [
            (-1.0, 0.0),
            (0.0, 1.0),
            (-10.0, -1.0),
            (-1000.0, -10.0),
            (1.0, 2.1),
            (2.1, 30.0),
            (30.0, 100.0),
        ] {
            for _ in 0..size {
                inputs.push(rng.uniform_f32(lo, hi));
            }
        }
        // Random encodings include subnormals, huge values, infinities and NaNs.
        // Every generated case is compared; there is no filtering.
        for _ in 0..size {
            inputs.push(rng.any_f32());
        }
    }

    // Explicitly cover both signs of every finite exponent, including the
    // subnormal/normal boundary and encodings immediately below each power of 2.
    for exponent in 0..=254 {
        for mantissa in [0, 1, 0x003f_ffff, 0x0040_0000, 0x007f_fffe, 0x007f_ffff] {
            let x = f32::from_bits((exponent << 23) | mantissa);
            inputs.extend([x, -x]);
        }
    }

    // Both sides of every linear interpolation knot, including the branch joins.
    for i in 0..=100 {
        neighbors(&mut inputs, i as f32 / 100.0);
        neighbors(&mut inputs, -(i as f32) / 100.0);
    }
    for i in 0..=90 {
        neighbors(&mut inputs, -1.0 - i as f32 / 10.0);
    }
    // Knots in y = 1.5 / (x * sqrt(x)), for both positive tables.
    for i in 0..=10 {
        for y in [1.5 - i as f64 / 10.0, 0.5 - i as f64 / 20.0] {
            if y > 0.0 {
                neighbors(&mut inputs, (1.5 / y).powf(2.0 / 3.0) as f32);
            }
        }
    }
    for i in 0..5 {
        let y = (5 - i) as f64 / 100.0;
        neighbors(&mut inputs, -((1.5 / y).powf(2.0 / 3.0) as f32));
    }
    for x in [
        -10.0,
        -1.0,
        0.0,
        1.0,
        3.0f64.powf(2.0 / 3.0) as f32,
        25.0,
        26.0,
        27.0,
        28.0,
        29.0,
        30.0,
        32.0,
        40.0,
        100.0,
        1000.0,
        100_000.0,
        1_000_000.0,
        1e20,
        1e30,
        1e38,
    ] {
        neighbors(&mut inputs, x);
        neighbors(&mut inputs, -x);
    }
    // Dense coverage of Ai underflow, Bi overflow, and the negative-table
    // endpoint where the original six-entry table is read with a zero weight.
    for _ in 0..10_000 {
        inputs.push(rng.uniform_f32(24.0, 32.0));
        inputs.push(rng.uniform_f32(-1_000_000.0, -100_000.0));
    }
    // The C tables include apparent transcription anomalies. Keep those values.
    for x in [-0.08, -0.12, -0.58] {
        neighbors(&mut inputs, x);
    }
    inputs
}

fn check(name: &str, seed: u64, rust: fn(f32) -> f32, c: unsafe extern "C" fn(f32) -> f32) {
    let inputs = inputs(seed);
    let mut different_bits = 0usize;
    let mut max_relative = 0.0f64;
    for &x in &inputs {
        let actual = rust(x);
        let expected = unsafe { c(x) };
        if expected.is_nan() {
            assert!(actual.is_nan(), "{name}({x:?}): Rust={actual:?}, C=NaN");
        } else if expected.is_infinite() {
            assert_eq!(actual, expected, "{name}({x:?}): infinity/sign mismatch");
        } else {
            assert!(
                actual.is_finite(),
                "{name}({x:?}): Rust={actual:?}, C={expected:?}"
            );
            let error = (actual as f64 - expected as f64).abs();
            let scale = (actual as f64).abs().max((expected as f64).abs());
            let relative = if scale == 0.0 { 0.0 } else { error / scale };
            max_relative = max_relative.max(relative);
            assert!(
                error <= 1e-5 * scale,
                "{name}({x:?}): Rust={actual:?}, C={expected:?}, relative={relative:e}"
            );
            if expected == 0.0 {
                // Relative error cannot distinguish +0 from -0 after underflow.
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "{name}({x:?}): zero/sign mismatch"
                );
            }
            if x > -10.0 && x < 1.0 {
                // Interpolation needs no transcendental library calls.
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "{name}({x:?}): table interpolation must be bit-exact"
                );
            }
            different_bits += usize::from(actual.to_bits() != expected.to_bits());
        }
    }
    eprintln!(
        "{name}: {} comparisons, {different_bits} finite results differ in bits, max relative error {max_relative:e}",
        inputs.len()
    );
}

#[test]
fn airya_matches_c() {
    check("airya", 12345, oxwave_f_airy::airya, airya);
}

#[test]
fn airyap_matches_c() {
    check("airyap", 12346, oxwave_f_airy::airyap, airyap);
}

#[test]
fn airyb_matches_c() {
    check("airyb", 12347, oxwave_f_airy::airyb, airyb);
}

#[test]
fn airybp_matches_c() {
    check("airybp", 12348, oxwave_f_airy::airybp, airybp);
}
