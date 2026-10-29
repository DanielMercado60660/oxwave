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

//! Differential tests for `cwp/lib/differentiate.c`.
//! Author: John Stockwell, CWP, 1994.

use oxwave_cref::XorShift;
use std::panic::{catch_unwind, AssertUnwindSafe};

extern "C" {
    fn differentiate(n: i32, h: f32, f: *mut f32, fprime: *mut f32);
    fn ddifferentiate(n: i32, h: f64, f: *mut f64, fprime: *mut f64);
}

const SIZES: [usize; 12] = [2, 3, 5, 7, 16, 31, 64, 127, 256, 1009, 4096, 8191];

// Run the same coverage independently for both precisions.
macro_rules! differential_suite {
    ($module:ident, $float:ty, $c:ident, $rust:path, $tiny:expr, $large:expr, $seed:expr) => {
        mod $module {
            use super::*;

            fn check(h: $float, f: &[$float]) {
                let n = f.len();
                let mut input = f.to_vec();
                let mut actual = vec![42.0 as $float; n + 2];
                let mut expected = actual.clone();
                // n >= 2: every C input access is valid, with n writable output
                // samples and guards on both sides. The routine is reentrant.
                unsafe {
                    $c(n as i32, h, input.as_mut_ptr(), expected[1..].as_mut_ptr());
                }
                $rust(h, f, &mut actual[1..]);
                // Require bit equality, stronger than 1e-5 relative error. C
                // does not specify NaN payloads, so compare their classification.
                for i in 0..n + 2 {
                    assert!(
                        actual[i].to_bits() == expected[i].to_bits()
                            || (actual[i].is_nan() && expected[i].is_nan()),
                        "n={n}, h={h:?}, index={i}: Rust {:?}, C {:?}",
                        actual[i],
                        expected[i]
                    );
                }
                assert_eq!(actual[0], 42.0);
                assert_eq!(actual[n + 1], 42.0);
                for i in 0..n {
                    assert_eq!(input[i].to_bits(), f[i].to_bits(), "C changed input {i}");
                }
            }

            #[test]
            fn random_matches_c() {
                let mut rng = XorShift::new($seed);
                let scales: [$float; 9] = [
                    0.0,
                    <$float>::from_bits(1),
                    <$float>::MIN_POSITIVE,
                    $tiny,
                    1e-4,
                    1.0,
                    1e10,
                    $large,
                    <$float>::MAX / 2.0,
                ];
                let intervals: [$float; 12] = [
                    0.0,
                    -0.0,
                    <$float>::from_bits(1),
                    <$float>::MIN_POSITIVE,
                    $tiny,
                    -1e-4,
                    0.004,
                    0.5,
                    -1.0,
                    1e10,
                    <$float>::MAX,
                    -<$float>::MAX,
                ];
                // Full Cartesian coverage avoids correlating signs, input
                // scales, lengths, and intervals. Generate three arrays for
                // each combination, totaling 3,888 C calls per function.
                let mut cases = 0;
                for n in SIZES {
                    for scale in scales {
                        for h in intervals {
                            for _ in 0..3 {
                                let f: Vec<$float> = (0..n)
                                    .map(|_| {
                                        ((rng.next_f64() * 2.0 - 1.0) * scale as f64) as $float
                                    })
                                    .collect();
                                check(h, &f);
                                cases += 1;
                            }
                        }
                    }
                }
                assert_eq!(cases, 3888);
            }

            #[test]
            fn random_bits_match_c() {
                let mut rng = XorShift::new($seed ^ 0x9e37_79b9_7f4a_7c15);
                // Vary the interval, length, sign, exponent, and mantissa
                // independently, including non-finite and subnormal values.
                for _ in 0..4096 {
                    let n = rng.range_i32(2, 257) as usize;
                    let h = <$float>::from_bits(rng.next_u64() as _);
                    let f: Vec<$float> = (0..n)
                        .map(|_| <$float>::from_bits(rng.next_u64() as _))
                        .collect();
                    check(h, &f);
                }
            }

            #[test]
            fn long_arrays_match_c() {
                let mut rng = XorShift::new($seed ^ 0xd1b5_4a32_d192_ed03);
                // Include large primes and a power of two with its neighbors.
                for n in [32749, 65535, 65536, 65537] {
                    for sign in [-1.0, 1.0] {
                        let h = (sign * rng.uniform(0.0001, 0.016)) as $float;
                        let f: Vec<$float> =
                            (0..n).map(|_| rng.uniform(-1e4, 1e4) as $float).collect();
                        check(h, &f);
                    }
                }
            }

            #[test]
            fn structured_edges_match_c() {
                for n in SIZES {
                    let mut impulse = vec![0.0 as $float; n];
                    impulse[n / 2] = -1.0;
                    let patterns: [Vec<$float>; 6] = [
                        vec![0.0; n],
                        vec![-0.0; n],
                        vec![-17.0; n],
                        (0..n).map(|i| i as $float * 0.25 - 17.0).collect(),
                        (0..n)
                            .map(|i| {
                                if i % 2 == 0 {
                                    <$float>::MAX
                                } else {
                                    -<$float>::MAX
                                }
                            })
                            .collect(),
                        impulse,
                    ];
                    for f in &patterns {
                        for h in [
                            0.0,
                            -0.0,
                            0.004,
                            -1.0,
                            <$float>::INFINITY,
                            <$float>::NEG_INFINITY,
                            <$float>::NAN,
                        ] {
                            check(h, f);
                        }
                    }
                }
                let special: [$float; 5] = [
                    0.0,
                    -0.0,
                    <$float>::INFINITY,
                    <$float>::NEG_INFINITY,
                    <$float>::NAN,
                ];
                for a in special {
                    for b in special {
                        check(1.0, &[a, b]);
                        check(-1.0, &[a, b, a]);
                    }
                }
                // Also exercise an output slice of exactly n elements.
                let mut expected: [$float; 2] = [0.0; 2];
                let mut input: [$float; 2] = [3.0, -1.0];
                unsafe {
                    $c(2, 0.5, input.as_mut_ptr(), expected.as_mut_ptr());
                }
                let mut actual: [$float; 2] = [0.0; 2];
                $rust(0.5, &input, &mut actual);
                assert_eq!(actual, expected);
            }

            #[test]
            fn rejects_invalid_lengths() {
                // For n=0 or n=1 the C endpoint formulas access f[1] and
                // f[n-2] outside the array. There is no defined C result;
                // exercise Rust's rejection without calling C out of bounds.
                let cases: [(Vec<$float>, usize); 3] =
                    [(vec![], 0), (vec![1.0], 1), (vec![1.0, 2.0], 1)];
                for (input, output_len) in cases {
                    let mut output = vec![42.0 as $float; output_len];
                    assert!(catch_unwind(AssertUnwindSafe(|| {
                        $rust(1.0, &input, &mut output);
                    }))
                    .is_err());
                    assert!(output.iter().all(|&x| x == 42.0));
                }
            }
        }
    };
}

differential_suite!(
    single,
    f32,
    differentiate,
    oxwave_f_differentiate::differentiate,
    1e-20,
    1e30,
    0x12345
);
differential_suite!(
    double,
    f64,
    ddifferentiate,
    oxwave_f_differentiate::ddifferentiate,
    1e-200,
    1e200,
    0x54321
);
