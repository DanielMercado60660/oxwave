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

//! Differential tests for cwp/lib/windows.c.
//! Author: Potash Corporation, Sascatchewan: Balasz Nemeth, given to CWP 2008.

use oxwave_cref::XorShift;

extern "C" {
    fn hanningnWindow(n: i32, w: *mut f32);
}

fn compare(initial: &[f32]) {
    // Include sentinels on both sides to check the C and Rust write bounds.
    let n = initial.len() - 2;
    let mut c = initial.to_vec();
    let mut rust = initial.to_vec();
    // The C loop overwrites every output element without reading its old value.
    unsafe { hanningnWindow(n as i32, c[1..].as_mut_ptr()) };
    oxwave_windows::hanningnWindow(&mut rust[1..n + 1]);
    assert_eq!(c[0].to_bits(), initial[0].to_bits());
    assert_eq!(c[n + 1].to_bits(), initial[n + 1].to_bits());
    assert_eq!(rust[0].to_bits(), initial[0].to_bits());
    assert_eq!(rust[n + 1].to_bits(), initial[n + 1].to_bits());
    let scale = c[1..n + 1].iter().fold(0.0_f32, |m, x| m.max(x.abs()));
    for i in 1..n + 1 {
        assert!(
            (c[i] - rust[i]).abs() <= 1e-5 * scale,
            "n={n}, index={}, C={}, Rust={}",
            i - 1,
            c[i],
            rust[i]
        );
        assert_eq!(
            c[i].to_bits(),
            rust[i].to_bits(),
            "bit mismatch: n={n}, index={}, C={}, Rust={}",
            i - 1,
            c[i],
            rust[i]
        );
    }
}

#[test]
fn hanningn_window_matches_c() {
    let mut rng = XorShift::new(0x1234_5678);
    for n in [
        0, 1, 2, 3, 5, 7, 16, 31, 64, 127, 256, 1021, 4093, 4096, 65521, 65536, 1_000_000,
    ] {
        let initial: Vec<f32> = (0..n + 2)
            .map(|_| rng.uniform_f32(-1000.0, 1000.0))
            .collect();
        compare(&initial);
    }
    for _ in 0..3000 {
        let n = (rng.next_u32() % 2049) as usize;
        let initial: Vec<f32> = (0..n + 2)
            .map(|_| rng.uniform_f32(-1000.0, 1000.0))
            .collect();
        compare(&initial);
    }
}

#[test]
fn hanningn_window_overwrites_zero_negative_and_nonfinite_values() {
    // There are no input samples: these are initial output-buffer contents.
    for n in [0, 1, 2, 3, 7, 31, 128, 1021, 4096] {
        for value in [
            0.0_f32,
            -0.0,
            -1.0,
            -f32::MAX,
            f32::MIN_POSITIVE,
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::NAN,
        ] {
            compare(&vec![value; n + 2]);
        }
    }
}
