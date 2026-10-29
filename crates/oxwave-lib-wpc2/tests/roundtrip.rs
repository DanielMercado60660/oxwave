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

//! Round-trip accuracy of wpcCompress/wpcUncompress on data that compresses
//! poorly (ratio below 4), which forces the encoder's buffer-overflow retry.
//! SU re-encodes an already transformed block on that retry, so its error can
//! exceed 100%; OxWave restores the block first (see DEVIATIONS.md). These
//! tests have no C reference on purpose: matching SU here means matching the bug.

use oxwave_lib_wpc2::*;

/// Noisy synthetic gather: a few dipping wavelets plus uniform noise from a
/// fixed linear congruential generator, so the test is deterministic.
fn gather(n1: usize, n2: usize, noise: f32) -> Vec<f32> {
    let mut f = vec![0.0f32; n1 * n2];
    let mut seed: u32 = 12345;
    for j in 0..n2 {
        for (k, t0) in [60.0f32, 180.0, 330.0].iter().enumerate() {
            let tc = t0 + j as f32 * (0.4 + 0.3 * k as f32);
            for i in 0..n1 {
                let x = (i as f32 - tc) / 6.0;
                f[j * n1 + i] += (1.0 - 2.0 * x * x) * (-x * x).exp();
            }
        }
        for i in 0..n1 {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            f[j * n1 + i] += noise * (((seed >> 8) as f32 / 16777216.0) - 0.5);
        }
    }
    f
}

fn rel_rms_after_roundtrip(f: &[f32], n1: usize, n2: usize, error: f32) -> (f64, usize) {
    let w = wpcCompress(f, n1, n2, error).expect("compression failed");
    let nbytes = wpcSize(&w) as usize;
    let mut g = vec![0.0f32; n1 * n2];
    wpcUncompress(&w, &mut g);
    let (mut num, mut den) = (0.0f64, 0.0f64);
    for (a, b) in f.iter().zip(&g) {
        num += ((a - b) as f64).powi(2);
        den += (*a as f64).powi(2);
    }
    ((num / den).sqrt(), nbytes)
}

#[test]
fn poorly_compressible_data_meets_requested_error() {
    let (n1, n2) = (512, 128);
    let f = gather(n1, n2, 1.0);
    let (e, nbytes) = rel_rms_after_roundtrip(&f, n1, n2, 0.01);
    // The retry path must actually be exercised: encoded size above n1*n2 bytes.
    assert!(nbytes > n1 * n2, "test data compressed too well ({nbytes} bytes); retry path not reached");
    assert!(e < 0.03, "relative RMS error {e} for error=0.01");
}

#[test]
fn well_compressible_data_unchanged_by_fix() {
    let (n1, n2) = (512, 128);
    let f = gather(n1, n2, 0.0);
    let (e, _) = rel_rms_after_roundtrip(&f, n1, n2, 0.05);
    assert!(e < 0.1, "relative RMS error {e} for error=0.05");
}

/// Small and random inputs nearly always overflow the initial n1*n2-byte buffer.
/// Every round trip must stay within 3x the requested relative RMS error.
#[test]
fn random_small_inputs_meet_requested_error() {
    let mut seed: u32 = 878;
    let mut next = || {
        seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
        ((seed >> 8) as f32 / 16777216.0) * 200.0 - 100.0
    };
    let mut retried = 0;
    for i in 0..2000 {
        let n1 = [1, 2, 3, 7, 31, 64, 127, 256, 257][i % 9];
        let n2 = if i % 499 == 498 { 257 } else { [1, 2, 3, 7, 17][(i / 9) % 5] };
        let mut f: Vec<f32> = (0..n1 * n2).map(|_| next()).collect();
        if i % 31 == 0 {
            f.fill(0.0);
        }
        let error = [0.01, 0.1, 1.0][i % 3];
        let w = wpcCompress(&f, n1, n2, error).expect("compression failed");
        if wpcSize(&w) as usize > n1 * n2 {
            retried += 1;
        }
        let mut g = vec![0.0f32; n1 * n2];
        wpcUncompress(&w, &mut g);
        let num: f64 = f.iter().zip(&g).map(|(a, b)| ((a - b) as f64).powi(2)).sum();
        let den: f64 = f.iter().map(|a| (*a as f64).powi(2)).sum();
        if den > 0.0 {
            let e = (num / den).sqrt();
            assert!(e <= 3.0 * error as f64, "case {i} ({n1}x{n2}): error {e} for error={error}");
        } else {
            assert!(g.iter().all(|x| *x == 0.0), "case {i}: zero input must decode to zero");
        }
    }
    assert!(retried > 1000, "only {retried} cases reached the retry path");
}
