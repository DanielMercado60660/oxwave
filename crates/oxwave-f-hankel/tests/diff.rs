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

//! Differential tests for cwp/lib/hankel.c.
//! Credits: Dave Hale, Colorado School of Mines, 06/04/90.

use oxwave_cref::XorShift;
use oxwave_f_hankel as rust;
use std::ffi::c_void;

extern "C" {
    fn hankelalloc(nfft: i32) -> *mut c_void;
    fn hankelfree(ht: *mut c_void);
    fn hankel0(ht: *mut c_void, f: *mut f32, h: *mut f32);
    fn hankel1(ht: *mut c_void, f: *mut f32, h: *mut f32);
}

// nfft/2 must be a valid PFA length. Include the smallest lengths, primes
// among both half-lengths and sample counts, every valid power of two, and
// large mixed-radix transforms. nfft=1 (one sample) is outside the C contract:
// Abel reads uninitialized coefficients and pfarc writes past gk's allocation.
const SIZES: [i32; 26] = [
    2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 60, 130, 240, 480, 720, 1008, 2002,
    2310, 8190, 32760,
];

struct CTransformer(*mut c_void);

impl CTransformer {
    fn new(nfft: i32) -> Self {
        let ht = unsafe { hankelalloc(nfft) };
        assert!(!ht.is_null());
        Self(ht)
    }
}

impl Drop for CTransformer {
    fn drop(&mut self) {
        unsafe { hankelfree(self.0) };
    }
}

fn input(rng: &mut XorShift, n: usize, trial: usize) -> Vec<f32> {
    let amplitude = [1.0, 1e-25, 1e-15, 1e15, 100.0][trial % 5];
    let mut f = rng.vec_f32(n, -amplitude, amplitude);
    match trial % 16 {
        0 => f.fill(0.0),
        1 => f.fill(-amplitude),
        2 => f.fill(amplitude),
        3..=5 => {
            f.fill(0.0);
            let index = match trial % 16 {
                3 => 0,
                4 => n / 2,
                _ => n - 1,
            };
            f[index] = -amplitude;
        }
        6 => {
            for x in &mut f {
                *x = -x.abs();
            }
        }
        7 => {
            for (i, x) in f.iter_mut().enumerate() {
                *x = if i % 2 == 0 { amplitude } else { -amplitude };
            }
        }
        _ => {}
    }
    f
}

fn assert_output(actual: &[f32], expected: &[f32], nfft: i32, trial: usize) {
    let mut magnitude = 0.0f64;
    for &x in expected {
        assert!(
            x.is_finite(),
            "nonfinite C output: nfft={nfft}, trial={trial}"
        );
        magnitude = magnitude.max((x as f64).abs());
    }
    for (i, (&a, &e)) in actual.iter().zip(expected).enumerate() {
        assert!(
            a.is_finite() && (a as f64 - e as f64).abs() <= 1e-5 * magnitude,
            "nfft={nfft} trial={trial} i={i} actual={a} expected={e} max={magnitude}",
        );
    }
}

fn check_transform(
    c: &CTransformer,
    r: &mut rust::HankelTransformer,
    nfft: i32,
    f: &[f32],
    first_order: bool,
    trial: usize,
) {
    let n = nfft as usize / 2 + 1;
    // Extra slots verify that neither implementation overwrites the tail.
    let sentinel = f32::from_bits(0x4abc_def0);
    let mut c_input = f.to_vec();
    let mut expected = vec![sentinel; n + 3];
    let mut actual = expected.clone();
    if first_order {
        unsafe { hankel1(c.0, c_input.as_mut_ptr(), expected.as_mut_ptr()) };
        rust::hankel1(r, f, &mut actual);
    } else {
        unsafe { hankel0(c.0, c_input.as_mut_ptr(), expected.as_mut_ptr()) };
        rust::hankel0(r, f, &mut actual);
    }
    for (&before, &after) in f.iter().zip(&c_input) {
        assert_eq!(before.to_bits(), after.to_bits(), "C changed its input");
    }
    for i in n..n + 3 {
        assert_eq!(expected[i].to_bits(), sentinel.to_bits());
        assert_eq!(actual[i].to_bits(), sentinel.to_bits());
    }
    assert_output(&actual[..n], &expected[..n], nfft, trial);
}

fn initialize(c: &CTransformer, r: &mut rust::HankelTransformer, nfft: i32) {
    // C's first hankel1 reads an uninitialized gx[0]. Initialize it through
    // hankel0 so comparisons exercise defined C behavior, including reuse.
    let f = vec![0.0; nfft as usize / 2 + 1];
    check_transform(c, r, nfft, &f, false, 0);
}

fn compare(first_order: bool, seed: u64) {
    let mut rng = XorShift::new(seed);
    for &nfft in &SIZES {
        let n = nfft as usize / 2 + 1;
        let c = CTransformer::new(nfft);
        let mut r = rust::hankelalloc(nfft);
        initialize(&c, &mut r, nfft);
        // 6,656 cases per function, reusing each workspace across inputs.
        for trial in 0..256 {
            let f = input(&mut rng, n, trial);
            check_transform(&c, &mut r, nfft, &f, first_order, trial);
        }
        rust::hankelfree(r);
    }
}

#[test]
fn hankel0_matches_c() {
    compare(false, 123);
}

#[test]
fn hankel1_matches_c() {
    compare(true, 456);
}

#[test]
fn fresh_hankel1_matches_initialized_c() {
    let mut rng = XorShift::new(0x3c6ef372fe94f82b);
    for trial in 0..2000 {
        let nfft = SIZES[rng.next_u32() as usize % SIZES.len()];
        let n = nfft as usize / 2 + 1;
        let c = CTransformer::new(nfft);
        // Initialize only C: its uninitialized gx[0] otherwise prevents a
        // defined comparison. Rust must work on its very first hankel1 call.
        let mut zero = vec![0.0; n];
        let mut initialized = vec![f32::NAN; n];
        unsafe { hankel0(c.0, zero.as_mut_ptr(), initialized.as_mut_ptr()) };
        assert!(initialized.iter().all(|&x| x == 0.0));

        let mut r = rust::hankelalloc(nfft);
        let f = input(&mut rng, n, trial);
        check_transform(&c, &mut r, nfft, &f, true, trial);
        rust::hankelfree(r);
    }
}

#[test]
fn alternating_transforms_match_c() {
    let mut rng = XorShift::new(0x6a09e667f3bcc909);
    for &nfft in &SIZES {
        let n = nfft as usize / 2 + 1;
        let c = CTransformer::new(nfft);
        let mut r = rust::hankelalloc(nfft);
        for trial in 0..256 {
            let f = input(&mut rng, n, trial);
            check_transform(&c, &mut r, nfft, &f, trial % 2 != 0, trial);
        }
        rust::hankelfree(r);
    }
}

#[test]
fn hankelalloc_matches_c() {
    let mut rng = XorShift::new(789);
    for trial in 0..2000 {
        let nfft = SIZES[rng.next_u32() as usize % SIZES.len()];
        let c = CTransformer::new(nfft);
        let mut r = rust::hankelalloc(nfft);
        // Compare precomputed coefficients through an observable transform,
        // rather than merely checking that the allocations succeeded.
        let f = input(&mut rng, nfft as usize / 2 + 1, trial);
        check_transform(&c, &mut r, nfft, &f, false, trial);
        rust::hankelfree(r);
    }
}

#[test]
fn hankelfree_matches_c() {
    let mut rng = XorShift::new(987);
    for trial in 0..2000 {
        let nfft = SIZES[rng.next_u32() as usize % SIZES.len()];
        let survivor_c = CTransformer::new(nfft);
        let mut survivor_r = rust::hankelalloc(nfft);
        let other_nfft = SIZES[rng.next_u32() as usize % SIZES.len()];
        let other_c = CTransformer::new(other_nfft);
        let mut other_r = rust::hankelalloc(other_nfft);
        let f = input(&mut rng, nfft as usize / 2 + 1, trial);
        let other_f = input(&mut rng, other_nfft as usize / 2 + 1, trial);
        check_transform(&survivor_c, &mut survivor_r, nfft, &f, false, trial);
        check_transform(&other_c, &mut other_r, other_nfft, &other_f, false, trial);
        // Free a populated workspace while a distinct transformer remains live.
        drop(other_c);
        rust::hankelfree(other_r);
        check_transform(&survivor_c, &mut survivor_r, nfft, &f, true, trial);
        drop(survivor_c);
        rust::hankelfree(survivor_r);
    }
}

#[test]
fn c_in_place_transforms_match_rust() {
    let mut rng = XorShift::new(0xbb67ae8584caa73b);
    for &nfft in &SIZES {
        let n = nfft as usize / 2 + 1;
        let c = CTransformer::new(nfft);
        let mut r = rust::hankelalloc(nfft);
        initialize(&c, &mut r, nfft);
        for trial in 0..256 {
            let f = input(&mut rng, n, trial);
            let mut expected = f.clone();
            let mut actual = vec![0.0; n];
            let ptr = expected.as_mut_ptr();
            if trial % 2 == 0 {
                unsafe { hankel0(c.0, ptr, ptr) };
                rust::hankel0(&mut r, &f, &mut actual);
            } else {
                unsafe { hankel1(c.0, ptr, ptr) };
                rust::hankel1(&mut r, &f, &mut actual);
            }
            assert_output(&actual, &expected, nfft, trial);
        }
        rust::hankelfree(r);
    }
}
