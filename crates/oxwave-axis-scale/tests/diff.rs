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

//! Differential tests for cwp/lib/scaxis.c.
//! Author: Dave Hale, Colorado School of Mines, 01/13/89.

use oxwave_cref::XorShift;

extern "C" {
    fn scaxis(x1: f32, x2: f32, nxnum: *mut i32, dxnum: *mut f32, fxnum: *mut f32);
}

fn compare(x1: f32, x2: f32, count: i32) {
    let (mut cn, mut cd, mut cf) = (count, 0.0, 0.0);
    let (mut rn, mut rd, mut rf) = (count, 0.0, 0.0);
    unsafe { scaxis(x1, x2, &mut cn, &mut cd, &mut cf) };
    oxwave_axis_scale::scaxis(x1, x2, &mut rn, &mut rd, &mut rf);
    assert_eq!(rn, cn, "count: {x1}, {x2}, {count}");
    assert_eq!(rd.to_bits(), cd.to_bits(), "increment: {x1}, {x2}, {count}");
    assert_eq!(rf.to_bits(), cf.to_bits(), "first: {x1}, {x2}, {count}");
}

#[test]
fn scaxis_matches_c() {
    let counts = [-100, 0, 1, 2, 3, 7, 16, 31, 64, 101, 1024, 4096];
    for &(a, b) in &[(0.0, 0.0), (-0.0, 0.0), (1.0, 1.0), (-5.0, -5.0),
        (0.0, 1.0), (1.0, 0.0), (-1.0, 1.0), (-100.0, -1.0),
        (1e-20, 9e-20), (-1e20, 1e20), (0.001, 0.01)] {
        for &n in &counts { compare(a, b, n); }
    }
    let mut rng = XorShift::new(0x736361786973);
    for _ in 0..20_000 {
        let scale = 10.0f64.powf(rng.range_i32(-20, 20) as f64) as f32;
        let a = rng.uniform_f32(-10.0, 10.0) * scale;
        let b = rng.uniform_f32(-10.0, 10.0) * scale;
        let n = counts[(rng.next_u32() as usize) % counts.len()];
        compare(a, b, n);
    }
}

#[test]
fn scaxis_rejects_unrepresentable_scale() {
    if std::env::var_os("OXWAVE_SCAXIS_ERROR_CASE").is_some() {
        let mut n = 4096;
        let (mut d, mut f) = (0.0, 0.0);
        oxwave_axis_scale::scaxis(1.0, f32::from_bits(1.0f32.to_bits() + 1), &mut n, &mut d, &mut f);
        panic!("unrepresentable scale should exit");
    }
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "scaxis_rejects_unrepresentable_scale", "--nocapture"])
        .env("OXWAVE_SCAXIS_ERROR_CASE", "1")
        .output().unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&output.stderr).contains("scaxis: invalid or unrepresentable axis scale"));
}

#[test]
fn scaxis_rejects_non_progressing_tick() {
    if std::env::var_os("OXWAVE_SCAXIS_STALLED_CASE").is_some() {
        let mut n = 922;
        let (mut d, mut f) = (0.0, 0.0);
        oxwave_axis_scale::scaxis(
            f32::from_bits(0x2f22_d456),
            f32::from_bits(0x2f22_d46b),
            &mut n, &mut d, &mut f,
        );
        panic!("non-progressing tick should exit");
    }
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "scaxis_rejects_non_progressing_tick", "--nocapture"])
        .env("OXWAVE_SCAXIS_STALLED_CASE", "1")
        .output().unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&output.stderr).contains("scaxis: tick increment does not advance"));
}
