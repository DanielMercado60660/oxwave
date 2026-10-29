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

//! Differential tests for cwp/lib/maxmin.c.
//! Author: Balasz Nemeth: Potash Corporation, given to CWP in 2008.

use oxwave_array_extrema::{max_index as rust_max_index, min_index as rust_min_index};
use oxwave_cref::XorShift;

extern "C" {
    fn max_index(n: i32, a: *mut f32, inc: i32) -> i32;
    fn min_index(n: i32, a: *mut f32, inc: i32) -> i32;
}

fn check(
    c: unsafe extern "C" fn(i32, *mut f32, i32) -> i32,
    rust: fn(&[f32], i32) -> i32,
    seed: u64,
) {
    let mut rng = XorShift::new(seed);
    let sizes = [
        0, 1, 2, 3, 7, 16, 31, 64, 127, 256, 1021, 4096, 65521, 65536,
    ];
    let increments = [1, 2, 3, 7, 17, 4097, 65537];
    for case in 0..6000 {
        let n = sizes[case % sizes.len()];
        let inc = increments[(case / sizes.len()) % increments.len()];
        let mut a = Vec::with_capacity(n);
        for _ in 0..n {
            let bits = rng.next_u32();
            let value = match (case / (sizes.len() * increments.len())) % 6 {
                0 => (bits as i32 as f32) / 2147483648.0 * 10000.0,
                1 => (bits % 9) as f32 - 4.0,
                2 => 0.0,
                3 => -((bits % 10000) as f32),
                4 => f32::from_bits(bits),
                _ => (bits as i32 as f32) * 1e-30,
            };
            a.push(value);
        }
        // Positive increments here cannot overflow C's int loop counter.
        // With n == 0, C returns 0 without dereferencing the pointer.
        let expected = unsafe { c(n as i32, a.as_mut_ptr(), inc) };
        assert_eq!(rust(&a, inc), expected, "case {case}, n={n}, inc={inc}");
    }
    for mut a in [
        vec![f32::NAN, 1.0, -1.0],
        vec![f32::NAN; 7],
        vec![0.0, f32::NAN, f32::INFINITY, f32::NEG_INFINITY],
        vec![-0.0, 0.0, -0.0],
        vec![2.0, 2.0, -3.0, -3.0],
        vec![f32::INFINITY, f32::NEG_INFINITY, f32::INFINITY],
        vec![f32::NEG_INFINITY, f32::INFINITY, f32::NEG_INFINITY],
        vec![
            0.0,
            f32::from_bits(1),
            -f32::from_bits(1),
            f32::MIN_POSITIVE,
        ],
        vec![0.0, f32::MAX, -f32::MAX],
    ] {
        for inc in increments {
            let expected = unsafe { c(a.len() as i32, a.as_mut_ptr(), inc) };
            assert_eq!(rust(&a, inc), expected, "edge values {a:?}, inc={inc}");
        }
    }

    // Exhaust every possible extreme position and tie pair in small arrays.
    // Off-stride extremes must be ignored, including at multiples of inc:
    // the C scan starts at index 1, rather than at index inc.
    for n in 1..=19 {
        for inc in 1..=n + 1 {
            for first in 0..n {
                for second in first..n {
                    for extreme in [-10.0, 10.0] {
                        let mut a = vec![0.0; n];
                        a[first] = extreme;
                        a[second] = extreme;
                        let expected = unsafe { c(n as i32, a.as_mut_ptr(), inc as i32) };
                        assert_eq!(
                            rust(&a, inc as i32),
                            expected,
                            "n={n}, inc={inc}, positions={first},{second}, extreme={extreme}"
                        );
                    }
                }
            }
        }
    }

    // For n <= 1 the loop never executes, even for nonpositive increments.
    for mut a in [vec![], vec![f32::NAN], vec![-0.0]] {
        for inc in [i32::MIN, -1, 0, 1, i32::MAX] {
            let expected = unsafe { c(a.len() as i32, a.as_mut_ptr(), inc) };
            assert_eq!(rust(&a, inc), expected);
        }
    }
}

#[test]
fn max_index_matches_c() {
    check(max_index, rust_max_index, 1);
}

#[test]
fn min_index_matches_c() {
    check(min_index, rust_min_index, 2);
}
