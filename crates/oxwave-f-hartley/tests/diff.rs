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

//! Differential tests for cwp/lib/hartley.c.
//! Credits: CENPET: Werner M. Heigl, April 2006.

use oxwave_cref::XorShift;
use oxwave_f_hartley as rust;

extern "C" {
    fn srfht(n: *mut i32, m: *mut i32, f: *mut f32);
    fn dsrfht(n: *mut i32, m: *mut i32, f: *mut f64);
    fn r4fht(n: i32, m: i32, f: *mut f32);
    fn nextpow2(n: i32) -> i32;
    fn nextpow4(n: i32) -> i32;
}

const LIMIT: i32 = 16_777_216;

#[test]
fn nextpow2_matches_c() {
    let mut rng = XorShift::new(11);
    for n in [i32::MIN, -1_000_003, -1, 0, 1, 2, 3, 5, 7, 17, 997, 65_537] {
        assert_eq!(rust::nextpow2(n), unsafe { nextpow2(n) }, "n={n}");
    }
    for exponent in 0..=24 {
        let power = 1i32 << exponent;
        for n in [power - 1, power] {
            assert_eq!(rust::nextpow2(n), unsafe { nextpow2(n) }, "n={n}");
        }
        if power < LIMIT {
            let n = power + 1;
            assert_eq!(rust::nextpow2(n), unsafe { nextpow2(n) }, "n={n}");
        }
    }
    for trial in 0..20_000 {
        let n = if trial % 2 == 0 {
            rng.range_i32(0, LIMIT)
        } else {
            rng.range_i32(i32::MIN, -1)
        };
        assert_eq!(rust::nextpow2(n), unsafe { nextpow2(n) }, "n={n}");
    }
}

#[test]
fn nextpow4_matches_c() {
    let mut rng = XorShift::new(12);
    for n in [i32::MIN, -1_000_003, -1, 0, 1, 2, 3, 5, 7, 17, 997, 65_537] {
        assert_eq!(rust::nextpow4(n), unsafe { nextpow4(n) }, "n={n}");
    }
    for exponent in 0..=12 {
        let power = 1i32 << (2 * exponent);
        for n in [power - 1, power] {
            assert_eq!(rust::nextpow4(n), unsafe { nextpow4(n) }, "n={n}");
        }
        if power < LIMIT {
            let n = power + 1;
            assert_eq!(rust::nextpow4(n), unsafe { nextpow4(n) }, "n={n}");
        }
    }
    for trial in 0..20_000 {
        let n = if trial % 2 == 0 {
            rng.range_i32(0, LIMIT)
        } else {
            rng.range_i32(i32::MIN, -1)
        };
        assert_eq!(rust::nextpow4(n), unsafe { nextpow4(n) }, "n={n}");
    }
}

// Structured cases exercise cancellation, DC, negative values, impulses, signed
// zeros, and mixed magnitudes in addition to uniform random sequences.
fn sample(rng: &mut XorShift, trial: usize, i: usize, impulse: usize, dc: f64) -> f64 {
    match trial % 16 {
        0 => 0.0,
        1 => -0.0,
        2 => dc,
        3 => -dc.abs(),
        4 => {
            if i == impulse {
                -1.0
            } else {
                0.0
            }
        }
        5 => {
            if i % 2 == 0 {
                1.0
            } else {
                -1.0
            }
        }
        6 => -rng.next_f64(),
        7 => rng.uniform(-1.0, 1.0) * if i % 2 == 0 { 1.0 } else { 1e-8 },
        _ => rng.uniform(-1.0, 1.0),
    }
}

fn input_f32(rng: &mut XorShift, trial: usize, len: usize) -> Vec<f32> {
    let scale = [1e-30f32, 1e-12, 1.0, 1e12, 1e30][rng.next_u32() as usize % 5];
    let impulse = rng.next_u32() as usize % len;
    let dc = rng.uniform(-1.0, 1.0);
    (0..len)
        .map(|i| sample(rng, trial, i, impulse, dc) as f32 * scale)
        .collect()
}

fn input_f64(rng: &mut XorShift, trial: usize, len: usize) -> Vec<f64> {
    let scale = [1e-250, 1e-100, 1.0, 1e100, 1e250][rng.next_u32() as usize % 5];
    let impulse = rng.next_u32() as usize % len;
    let dc = rng.uniform(-1.0, 1.0);
    (0..len)
        .map(|i| sample(rng, trial, i, impulse, dc) * scale)
        .collect()
}

fn check_f32(actual: &[f32], expected: &[f32], exact: bool, n: i32, trial: usize) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected
        .iter()
        .fold(0.0f64, |s, &x| s.max((x as f64).abs()));
    for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
        if exact {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "n={n} trial={trial} i={i}: {a} != {b}"
            );
        } else {
            assert!(a.is_finite() && b.is_finite());
            assert!(
                (a as f64 - b as f64).abs() <= 1e-5 * scale,
                "n={n} trial={trial} i={i}: {a} != {b}, scale={scale}"
            );
        }
    }
}

fn check_f64(actual: &[f64], expected: &[f64], exact: bool, n: i32, trial: usize) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().fold(0.0f64, |s, &x| s.max(x.abs()));
    for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
        if exact {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "n={n} trial={trial} i={i}: {a} != {b}"
            );
        } else {
            assert!(a.is_finite() && b.is_finite());
            assert!(
                (a - b).abs() <= 1e-5 * scale,
                "n={n} trial={trial} i={i}: {a} != {b}, scale={scale}"
            );
        }
    }
}

fn single_case(rng: &mut XorShift, trial: usize, m: i32) {
    let n = 1i32 << m;
    // C's n=1 last stage accesses two samples. Both allocations explicitly
    // contain them, plus guards that must remain unchanged for every length.
    let used = n.max(2) as usize;
    let mut a = input_f32(rng, trial, used + 3);
    let tail = a[used..].to_vec();
    let mut b = a.clone();
    let (mut cn, mut cm) = (n, m);
    unsafe { srfht(&mut cn, &mut cm, b.as_mut_ptr()) };
    assert_eq!((cn, cm), (n, m));
    rust::srfht(n, m, &mut a);
    check_f32(&a[..used], &b[..used], m <= 3, n, trial);
    for ((&a, &b), &original) in a[used..].iter().zip(&b[used..]).zip(&tail) {
        assert_eq!(a.to_bits(), original.to_bits());
        assert_eq!(b.to_bits(), original.to_bits());
    }
}

#[test]
fn srfht_matches_c() {
    let mut rng = XorShift::new(13);
    for trial in 0..3_000 {
        let m = rng.range_i32(0, 12);
        single_case(&mut rng, trial, m);
    }
    for m in 13..=18 {
        for trial in [0, 4, 15] {
            single_case(&mut rng, trial, m);
        }
    }
}

fn double_case(rng: &mut XorShift, trial: usize, m: i32) {
    let n = 1i32 << m;
    let used = n.max(2) as usize;
    let mut a = input_f64(rng, trial, used + 3);
    let tail = a[used..].to_vec();
    let mut b = a.clone();
    let (mut cn, mut cm) = (n, m);
    unsafe { dsrfht(&mut cn, &mut cm, b.as_mut_ptr()) };
    assert_eq!((cn, cm), (n, m));
    rust::dsrfht(n, m, &mut a);
    check_f64(&a[..used], &b[..used], m <= 3, n, trial);
    for ((&a, &b), &original) in a[used..].iter().zip(&b[used..]).zip(&tail) {
        assert_eq!(a.to_bits(), original.to_bits());
        assert_eq!(b.to_bits(), original.to_bits());
    }
}

#[test]
fn dsrfht_matches_c() {
    let mut rng = XorShift::new(14);
    for trial in 0..3_000 {
        let m = rng.range_i32(0, 12);
        double_case(&mut rng, trial, m);
    }
    for m in 13..=18 {
        for trial in [0, 4, 15] {
            double_case(&mut rng, trial, m);
        }
    }
}

fn radix4_case(rng: &mut XorShift, trial: usize, m: i32) {
    let n = 1i32 << (2 * m);
    // C's first radix-4 stage unconditionally accesses four samples for n=1.
    let used = n.max(4) as usize;
    let mut a = input_f32(rng, trial, used + 3);
    let tail = a[used..].to_vec();
    let mut b = a.clone();
    unsafe { r4fht(n, m, b.as_mut_ptr()) };
    rust::r4fht(n, m, &mut a);
    check_f32(&a[..used], &b[..used], m <= 1, n, trial);
    for ((&a, &b), &original) in a[used..].iter().zip(&b[used..]).zip(&tail) {
        assert_eq!(a.to_bits(), original.to_bits());
        assert_eq!(b.to_bits(), original.to_bits());
    }
}

#[test]
fn r4fht_matches_c() {
    let mut rng = XorShift::new(15);
    for trial in 0..3_000 {
        let m = rng.range_i32(0, 6);
        radix4_case(&mut rng, trial, m);
    }
    for m in 7..=9 {
        for trial in [0, 4, 15] {
            radix4_case(&mut rng, trial, m);
        }
    }
}

#[test]
fn invalid_parameters_are_rejected() {
    use std::panic::catch_unwind;
    // Non-power lengths can loop forever or access outside the C buffers.
    // Exercise the Rust preconditions without invoking undefined C behavior.
    for n in [-1, 0, 3, 5, 7, 17, 997] {
        assert!(catch_unwind(|| rust::srfht(n, 1, &mut [0.0; 1024])).is_err());
        assert!(catch_unwind(|| rust::dsrfht(n, 1, &mut [0.0; 1024])).is_err());
        assert!(catch_unwind(|| rust::r4fht(n, 1, &mut [0.0; 1024])).is_err());
    }
    assert!(catch_unwind(|| rust::r4fht(2, 1, &mut [0.0; 4])).is_err());
    assert!(catch_unwind(|| rust::srfht(1, 0, &mut [0.0; 1])).is_err());
    assert!(catch_unwind(|| rust::dsrfht(1, 0, &mut [0.0; 1])).is_err());
    assert!(catch_unwind(|| rust::r4fht(1, 0, &mut [0.0; 1])).is_err());
    assert!(catch_unwind(|| rust::srfht(8, 2, &mut [0.0; 8])).is_err());
    assert!(catch_unwind(|| rust::dsrfht(8, 2, &mut [0.0; 8])).is_err());
    assert!(catch_unwind(|| rust::r4fht(16, 1, &mut [0.0; 16])).is_err());
    assert!(catch_unwind(|| rust::srfht(8, 3, &mut [0.0; 7])).is_err());
    assert!(catch_unwind(|| rust::dsrfht(8, 3, &mut [0.0; 7])).is_err());
    assert!(catch_unwind(|| rust::r4fht(16, 2, &mut [0.0; 15])).is_err());
    for n in [LIMIT + 1, i32::MAX] {
        assert!(catch_unwind(|| rust::nextpow2(n)).is_err());
        assert!(catch_unwind(|| rust::nextpow4(n)).is_err());
    }
}
