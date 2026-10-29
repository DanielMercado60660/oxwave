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

//! Differential tests for cwp/lib/franuni.c.
//! Author: Dave Hale, Colorado School of Mines, 12/30/89.
//! Adapted from subroutine uni in Kahaner, Moler, and Nash (1988),
//! referencing unpublished notes by Marsaglia.

use oxwave_cref::XorShift;
use std::sync::Mutex;

extern "C" {
    fn franuni() -> f32;
    fn sranuni(seed: i32);
}

static LOCK: Mutex<()> = Mutex::new(());

fn seed_both(seed: i32) {
    unsafe { sranuni(seed) };
    oxwave_f_franuni::sranuni(seed);
}

fn compare_draws(seed: i32, len: usize) {
    for n in 0..len {
        let c = unsafe { franuni() };
        let rust = oxwave_f_franuni::franuni();
        assert_eq!(c.to_bits(), rust.to_bits(), "seed {seed}, draw {n}");
        assert!((0.0..1.0).contains(&rust), "seed {seed}, draw {n}");
    }
}

#[test]
fn franuni_matches_c() {
    let _guard = LOCK.lock().unwrap();
    for seed in [0, 1, -1, 42, -42, i32::MAX, i32::MIN, i32::MIN + 1] {
        seed_both(seed);
        compare_draws(seed, 100_000);
    }

    // Exercise index wraparound, primes, powers of two, and random states.
    let mut rng = XorShift::new(0x4652_414e_554e_49);
    for case in 0..4096 {
        let seed = rng.next_u32() as i32;
        seed_both(seed);
        let len = [1, 2, 3, 16, 17, 18, 31, 64, 257, 1024, 4093][case % 11];
        compare_draws(seed, len);
    }
}

#[test]
fn sranuni_matches_c() {
    let _guard = LOCK.lock().unwrap();
    let mut rng = XorShift::new(437);
    let edges = [0, 1, -1, 42, -42, i32::MAX, i32::MIN, i32::MIN + 1];
    for (case, seed) in edges
        .into_iter()
        .chain((0..10_000).map(|_| rng.next_u32() as i32))
        .enumerate()
    {
        seed_both(seed);
        let len = [1, 2, 17, 31, 64, 257][case % 6];
        compare_draws(seed, len);
    }

    // Values around every modulus used to turn the seed into initial state.
    for modulus in [147, 157, 167, 177] {
        for multiple in [1, 2, 17, 1024, 1_000_000] {
            for delta in [-1, 0, 1] {
                let magnitude = modulus * multiple + delta;
                for seed in [magnitude, -magnitude] {
                    seed_both(seed);
                    compare_draws(seed, 257);
                }
            }
        }
    }
}

#[test]
fn sranuni_resets_after_arbitrary_draws() {
    let _guard = LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x5245_5345_4544);
    for _ in 0..2048 {
        // Avoid INT_MIN here so negating the seed is representable.
        let seed = (rng.next_u32() & 0x7fff_ffff) as i32;
        seed_both(seed);
        let expected: Vec<u32> = (0..31)
            .map(|n| {
                let c = unsafe { franuni() };
                let rust = oxwave_f_franuni::franuni();
                assert_eq!(c.to_bits(), rust.to_bits(), "seed {seed}, draw {n}");
                c.to_bits()
            })
            .collect();
        compare_draws(seed, (rng.next_u32() % 1024) as usize);
        // Each seed must reset every state component, including the carry and
        // indices, after both a partial cycle and an unrelated seed's draws.
        for reset_seed in [seed, -seed, seed] {
            seed_both(reset_seed);
            for (n, &bits) in expected.iter().enumerate() {
                assert_eq!(
                    unsafe { franuni() }.to_bits(),
                    bits,
                    "C reset {reset_seed}, draw {n}"
                );
                assert_eq!(
                    oxwave_f_franuni::franuni().to_bits(),
                    bits,
                    "Rust reset {reset_seed}, draw {n}"
                );
            }
            // Consecutive calls must replace the entire state even when no
            // values are drawn between them, including zero and INT_MIN seeds.
            for intermediate in [0, i32::MIN, rng.next_u32() as i32] {
                seed_both(intermediate);
            }
            seed_both(reset_seed);
            for (n, &bits) in expected.iter().enumerate() {
                assert_eq!(
                    unsafe { franuni() }.to_bits(),
                    bits,
                    "C consecutive reset {reset_seed}, draw {n}"
                );
                assert_eq!(
                    oxwave_f_franuni::franuni().to_bits(),
                    bits,
                    "Rust consecutive reset {reset_seed}, draw {n}"
                );
            }
            let unrelated = rng.next_u32() as i32;
            seed_both(unrelated);
            compare_draws(unrelated, (rng.next_u32() % 257) as usize);
        }
    }
}

#[test]
fn unseeded_initial_state_matches_c() {
    const CHILD: &str = "OXWAVE_FRANUNI_INITIAL_STATE_CHILD";
    if std::env::var_os(CHILD).is_none() {
        let status = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "unseeded_initial_state_matches_c"])
            .env(CHILD, "1")
            .status()
            .unwrap();
        assert!(status.success());
        return;
    }
    // The fresh process runs this test alone, before either generator is seeded.
    for n in 0..100_000 {
        let c = unsafe { franuni() };
        let rust = oxwave_f_franuni::franuni();
        assert_eq!(c.to_bits(), rust.to_bits(), "unseeded draw {n}");
    }
}
