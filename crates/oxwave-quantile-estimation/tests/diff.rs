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

//! Differential tests against cwp/lib/quest.c.
//! Author of the C routines: Dave Hale, Colorado School of Mines, 05/07/89.

use oxwave_cref::XorShift;
use oxwave_quantile_estimation as rust;
use std::ffi::c_void;

// These are the original, unmodified libcwp symbols linked by oxwave-cref.
extern "C" {
    fn quest(p: f32, n: i32, x: *mut f32) -> f32;
    fn questalloc(p: f32, n: i32, x: *mut f32) -> *mut c_void;
    fn questfree(q: *mut c_void);
    fn questupdate(q: *mut c_void, n: i32, x: *mut f32) -> f32;
}

// C owns each allocation. RAII also releases it if an assertion panics.
struct CQuest(*mut c_void);

impl CQuest {
    fn new(p: f32, x: &mut [f32]) -> Self {
        assert!(x.len() >= 5);
        let q = unsafe { questalloc(p, x.len() as i32, x.as_mut_ptr()) };
        assert!(!q.is_null());
        Self(q)
    }

    fn update(&mut self, x: &mut [f32]) -> f32 {
        unsafe { questupdate(self.0, x.len() as i32, x.as_mut_ptr()) }
    }
}

impl Drop for CQuest {
    fn drop(&mut self) {
        unsafe { questfree(self.0) };
    }
}

fn same(actual: f32, expected: f32, context: &str) {
    if expected.is_nan() {
        assert!(
            actual.is_nan(),
            "{context}: Rust={actual:?}, C={expected:?}"
        );
    } else {
        // Stronger than 1e-5 relative: includes sign of zero and infinities.
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{context}: Rust={actual:?}, C={expected:?}"
        );
    }
}

fn same_input(actual: &[f32], expected: &[f32]) {
    assert_eq!(
        actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        "input must remain unchanged"
    );
}

const SIZES: &[usize] = &[5, 6, 7, 11, 16, 31, 64, 127, 256, 1021, 4093];
const CHUNKS: &[usize] = &[0, 1, 2, 3, 7, 16, 31, 64, 127, 256, 1021, 4093];

fn probability(rng: &mut XorShift, case: usize) -> f32 {
    match case % 12 {
        0 => 0.0,
        1 => 1.0,
        2 => 0.5,
        3 => 0.01,
        4 => 0.99,
        5 => f32::from_bits(1),
        6 => f32::from_bits(1.0f32.to_bits() - 1),
        7 => 0.25,
        8 => 0.75,
        _ => rng.uniform_f32(0.0, 1.0),
    }
}

fn samples(rng: &mut XorShift, n: usize, case: usize) -> Vec<f32> {
    let scale = match (case / 12) % 5 {
        0 => 1.0,
        1 => 1e-30,
        2 => 1e30,
        3 => 1e-5,
        _ => 1e5,
    };
    let offset = rng.uniform_f32(-100.0, 100.0);
    (0..n)
        .map(|i| match case % 12 {
            0 => 0.0,
            1 => -0.0,
            2 => offset,
            3 => (i as f32 - n as f32 / 2.0) * scale,
            4 => (n as f32 / 2.0 - i as f32) * scale,
            5 => {
                if i % 2 == 0 {
                    -scale
                } else {
                    scale
                }
            }
            6 => rng.range_i32(-3, 3) as f32 * scale,
            7 => -rng.uniform_f32(0.0, 1000.0) * scale,
            8 => {
                if i % 17 == 0 {
                    scale * 1e6
                } else {
                    scale
                }
            }
            9 => rng.uniform_f32(-1000.0, 1000.0) * scale,
            10 => rng.uniform_f32(-1.0, 1.0) * scale + offset,
            _ => rng.uniform_f32(-1000.0, 1000.0),
        })
        .collect()
}

#[test]
fn quest_matches_c() {
    let mut rng = XorShift::new(0x51e5_7a11);
    for case in 0..4096 {
        let p = probability(&mut rng, case);
        let mut x = samples(&mut rng, SIZES[case % SIZES.len()], case);
        let original = x.clone();
        let expected = unsafe { quest(p, x.len() as i32, x.as_mut_ptr()) };
        same(
            rust::quest(p, &x),
            expected,
            &format!("quest case={case}, p={p}"),
        );
        same_input(&x, &original);
    }
}

#[test]
fn questalloc_matches_c() {
    let mut rng = XorShift::new(0xa110_ca7e);
    for case in 0..4096 {
        let p = probability(&mut rng, case);
        let mut x = samples(&mut rng, SIZES[case % SIZES.len()], case);
        let original = x.clone();
        let mut cq = CQuest::new(p, &mut x);
        let mut rq = rust::questalloc(p, &x);
        same(
            rust::questupdate(&mut rq, &[]),
            cq.update(&mut []),
            "questalloc initial estimate",
        );

        // Continuing after allocation exercises all hidden state, including
        // cases where questalloc had to process more than the first five.
        for step in 0..8 {
            let mut more = samples(&mut rng, CHUNKS[(case + step) % CHUNKS.len()], case + step);
            same(
                rust::questupdate(&mut rq, &more),
                cq.update(&mut more),
                &format!("questalloc case={case}, step={step}, p={p}"),
            );
        }
        same_input(&x, &original);
        rust::questfree(rq);
    }
}

#[test]
fn questupdate_matches_c() {
    let mut rng = XorShift::new(0x0ada_7e55);
    for case in 0..4096 {
        let p = probability(&mut rng, case);
        let mut initial = samples(&mut rng, 5, case);
        let mut cq = CQuest::new(p, &mut initial);
        let mut rq = rust::questalloc(p, &initial);
        // Every stream sees each chunk size, in a different order, and an
        // empty update after every chunk. Updates may change distribution.
        for step in 0..CHUNKS.len() {
            let n = CHUNKS[(case + step) % CHUNKS.len()];
            let mut x = samples(&mut rng, n, case + step);
            let original = x.clone();
            same(
                rust::questupdate(&mut rq, &x),
                cq.update(&mut x),
                &format!("questupdate case={case}, step={step}, n={n}, p={p}"),
            );
            same(
                rust::questupdate(&mut rq, &[]),
                cq.update(&mut []),
                "empty update",
            );
            same_input(&x, &original);
        }
        rust::questfree(rq);
    }
}

#[test]
fn questfree_matches_c() {
    let mut rng = XorShift::new(0xf4ee_5a7e);
    let mut initial = [-3.0, -1.0, 0.0, 1.0, 3.0];
    let mut csurvivor = CQuest::new(0.5, &mut initial);
    let mut rsurvivor = rust::questalloc(0.5, &initial);

    for case in 0..4096 {
        // Two live allocations with different contents, released in alternating
        // order. A third state must remain valid across all of those frees.
        let p = probability(&mut rng, case);
        let mut x = samples(&mut rng, SIZES[case % SIZES.len()], case);
        let mut cfirst = CQuest::new(p, &mut x);
        let rfirst = rust::questalloc(p, &x);
        let mut csecond = CQuest::new(1.0 - p, &mut x);
        let rsecond = rust::questalloc(1.0 - p, &x);

        if case % 2 == 0 {
            drop(cfirst);
            rust::questfree(rfirst);
            let mut rsecond = rsecond;
            same(
                rust::questupdate(&mut rsecond, &[]),
                csecond.update(&mut []),
                "second survives free",
            );
            drop(csecond);
            rust::questfree(rsecond);
        } else {
            drop(csecond);
            rust::questfree(rsecond);
            let mut rfirst = rfirst;
            same(
                rust::questupdate(&mut rfirst, &[]),
                cfirst.update(&mut []),
                "first survives free",
            );
            drop(cfirst);
            rust::questfree(rfirst);
        }

        let mut more = samples(&mut rng, case % 8, case);
        same(
            rust::questupdate(&mut rsurvivor, &more),
            csurvivor.update(&mut more),
            "independent estimator survives repeated frees",
        );
    }
    drop(csurvivor);
    rust::questfree(rsurvivor);
}

#[test]
fn large_arrays_and_partitioned_streams_match_c() {
    let mut rng = XorShift::new(0x1a49_e001);
    for n in [8192, 65521, 65536, 100003] {
        for case in 0..12 {
            let p = probability(&mut rng, case);
            let mut x = samples(&mut rng, n, case + 60);
            let expected = unsafe { quest(p, n as i32, x.as_mut_ptr()) };
            same(rust::quest(p, &x), expected, "large quest");

            let mut cq = CQuest::new(p, &mut x);
            let mut rq = rust::questalloc(p, &x);
            same(
                rust::questupdate(&mut rq, &[]),
                cq.update(&mut []),
                "large questalloc",
            );

            // Chunking must preserve exactly the streaming algorithm's state.
            let mut cpartitioned = CQuest::new(p, &mut x[..5]);
            let mut rpartitioned = rust::questalloc(p, &x[..5]);
            let mut pos = 5;
            while pos < n {
                let end = (pos + 1 + rng.next_u32() as usize % 4093).min(n);
                same(
                    rust::questupdate(&mut rpartitioned, &x[pos..end]),
                    cpartitioned.update(&mut x[pos..end]),
                    "large partitioned update",
                );
                pos = end;
            }
            same(
                rust::questupdate(&mut rpartitioned, &[]),
                cq.update(&mut []),
                "chunked versus whole C streaming result",
            );
            rust::questfree(rq);
            rust::questfree(rpartitioned);
        }
    }
}

#[test]
fn signed_zeros_extreme_floats_and_nonfinite_samples_match_c() {
    let mut rng = XorShift::new(0xed9e_ca5e);
    let values = [
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
        1.0,
        -1.0,
    ];
    for case in 0..2048 {
        let p = probability(&mut rng, case);
        let n = SIZES[case % SIZES.len()];
        let mut x: Vec<f32> = (0..n)
            .map(|_| values[rng.next_u32() as usize % values.len()])
            .collect();
        let original = x.clone();
        let expected = unsafe { quest(p, n as i32, x.as_mut_ptr()) };
        same(rust::quest(p, &x), expected, "extreme quest");
        let mut cq = CQuest::new(p, &mut x);
        let mut rq = rust::questalloc(p, &x);
        same(
            rust::questupdate(&mut rq, &[]),
            cq.update(&mut []),
            "extreme questalloc",
        );
        for more in x.chunks_mut(7) {
            same(
                rust::questupdate(&mut rq, more),
                cq.update(more),
                "extreme update",
            );
        }
        same_input(&x, &original);
        rust::questfree(rq);
    }
}

#[test]
fn initialization_requires_five_samples() {
    // C documents n>=5 and unconditionally reads x[0..5], so invoking C here
    // would be undefined behavior. Rust rejects these lengths explicitly.
    for n in 0..5 {
        let x = vec![0.0; n];
        assert!(std::panic::catch_unwind(|| rust::quest(0.5, &x)).is_err());
        assert!(std::panic::catch_unwind(|| rust::questalloc(0.5, &x)).is_err());
    }
}
