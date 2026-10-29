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

//! Differential tests for su/lib/elco_scalar.c.
//! Credits: Colorado School of Mines (copyright above); no individual credits in the C source.

use oxwave_cref::XorShift;
use oxwave_segy_scalar as rust;

extern "C" {
    fn elco_scalar(ncoords: i32, c: *mut f64) -> i16;
    fn from_segy_elco_multiplier(scalar: i16) -> f64;
    fn to_segy_elco_multiplier(scalar: i16) -> f64;
}

fn compare(coords: &[f64]) {
    let mut c = coords.to_vec();
    let expected = unsafe { elco_scalar(c.len() as i32, c.as_mut_ptr()) };
    assert_eq!(rust::elco_scalar(coords), expected, "coordinates: {coords:?}");
}

#[test]
fn elco_scalar_matches_c() {
    compare(&[]);
    compare(&[0.0, -0.0, f64::NAN]);
    compare(&[f64::INFINITY, f64::NEG_INFINITY]);
    let thresholds: [f64; 8] = [2147483647500.0, 214748364750.0, 21474836475.0,
        2147483647.5, 214748364.75, 21474836.475, 2147483.6475, 214748.36475];
    for x in thresholds {
        for bits in [x.to_bits() - 1, x.to_bits(), x.to_bits() + 1] {
            let value = f64::from_bits(bits);
            compare(&[value]);
            compare(&[-value, 0.0, f64::NAN]);
        }
    }
    let mut rng = XorShift::new(173);
    for size in [1, 2, 3, 7, 16, 31, 64, 127, 1024, 4096] {
        for _ in 0..400 {
            let scale = 10.0f64.powf((rng.next_u32() % 18) as f64 - 5.0);
            let coords: Vec<f64> = (0..size).map(|_| {
                match rng.next_u32() % 20 {
                    0 => 0.0,
                    1 => -0.0,
                    _ => rng.uniform(-1.0, 1.0) * scale,
                }
            }).collect();
            compare(&coords);
        }
    }
}

#[test]
fn from_segy_elco_multiplier_matches_c() {
    for scalar in i16::MIN..=i16::MAX {
        let c = unsafe { from_segy_elco_multiplier(scalar) };
        assert_eq!(rust::from_segy_elco_multiplier(scalar).to_bits(), c.to_bits(), "scalar {scalar}");
    }
}

#[test]
fn to_segy_elco_multiplier_matches_c() {
    for scalar in i16::MIN..=i16::MAX {
        let c = unsafe { to_segy_elco_multiplier(scalar) };
        assert_eq!(rust::to_segy_elco_multiplier(scalar).to_bits(), c.to_bits(), "scalar {scalar}");
    }
}

