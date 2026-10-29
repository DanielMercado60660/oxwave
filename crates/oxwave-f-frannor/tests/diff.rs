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

//! Differential tests for cwp/lib/frannor.c.
//! Author: Dave Hale, Colorado School of Mines, 01/21/89.
//! Adapted from Kahaner, Moler and Nash (1988), based on Marsaglia and Tsang (1984).

use oxwave_cref::XorShift;
use std::sync::Mutex;

static LOCK: Mutex<()> = Mutex::new(());

extern "C" {
    fn frannor() -> f32;
    fn srannor(seed: i32);
}

fn compare_sequence(count: usize, context: &str) {
    for index in 0..count {
        let expected = unsafe { frannor() };
        let actual = oxwave_f_frannor::frannor();
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{context}, sample {index}: Rust {actual:?}, C {expected:?}"
        );
    }
}

fn compare(seed: i32, count: usize) {
    unsafe { srannor(seed) };
    oxwave_f_frannor::srannor(seed);
    compare_sequence(count, &format!("seed {seed}"));
}

#[test]
fn frannor_default_state_matches_c() {
    // A separate process keeps the original static state intact regardless of
    // the order or concurrency of the seeded tests in the parent process.
    const CHILD: &str = "OXWAVE_FRANNOR_DEFAULT_STATE_CHILD";
    if std::env::var_os(CHILD).is_some() {
        compare_sequence(100_000, "default state");
    } else {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "frannor_default_state_matches_c", "--nocapture"])
            .env(CHILD, "1")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "child test failed:\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn frannor_matches_c() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    for (seed, count) in [
        (0, 1),
        (1, 2),
        (-1, 3),
        (32707, 127),
        (i32::MAX, 1024),
        (42, 100_000),
    ] {
        compare(seed, count);
    }
    let mut rng = XorShift::new(0x9e37_79b9_7f4a_7c15);
    // These routines take no array arguments: exercise sequence lengths instead.
    let lengths = [1, 2, 3, 7, 17, 31, 64, 127, 256, 1024];
    for trial in 0..4_000 {
        compare(rng.next_u32() as i32, lengths[trial % lengths.len()]);
    }
}

#[test]
fn srannor_matches_c() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0x12345678);
    for seed in [
        i32::MIN,
        i32::MIN + 1,
        i32::MAX,
        0,
        1,
        -1,
        32706,
        32707,
        32708,
        -32706,
        -32707,
        -32708,
    ] {
        compare(seed, 128);
    }
    for _ in 0..4_000 {
        let seed = rng.next_u32() as i32;
        compare(seed, 64);
        // Reinitialize after a different amount of state has been consumed.
        compare_sequence((rng.next_u32() % 131) as usize, "before reseeding");
        compare(seed, 64);
    }
    // ABS(seed) and the modulus collapse these pairs to the same initial state.
    for seed in [0, 1, 1111, 1947, 32706] {
        for equivalent in [seed, -seed, seed + 32707, -(seed + 32707)] {
            unsafe { srannor(seed) };
            oxwave_f_frannor::srannor(equivalent);
            compare_sequence(128, "equivalent seeds");
        }
    }
}

#[test]
fn srannor_all_seed_residues_match_c() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    // Except for INT_MIN's C signed overflow, ABS(seed) % 32707 reduces every
    // i32 seed to one of these states. Check every residue with both signs and
    // on both sides of the first modulus boundary. The 64 draws wrap the
    // 17-element generator state several times, including after each reseed.
    for seed in 0..32707 {
        for equivalent in [seed, -seed, seed + 32707, -(seed + 32707)] {
            compare(equivalent, 64);
        }
    }
}
