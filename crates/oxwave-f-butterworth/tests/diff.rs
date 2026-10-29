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

//! Differential tests for cwp/lib/butterworth.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

use oxwave_cref::XorShift;
use oxwave_f_butterworth as rust;

extern "C" {
    fn bfdesign(fpass: f32, apass: f32, fstop: f32, astop: f32, npoles: *mut i32, f3db: *mut f32);
    fn bfhighpass(npoles: i32, f3db: f32, n: i32, p: *mut f32, q: *mut f32);
    fn bflowpass(npoles: i32, f3db: f32, n: i32, p: *mut f32, q: *mut f32);
}

fn compare_design(fpass: f32, apass: f32, fstop: f32, astop: f32) {
    let (mut cn, mut cf) = (-1, f32::NAN);
    let (mut rn, mut rf) = (-1, f32::NAN);
    unsafe { bfdesign(fpass, apass, fstop, astop, &mut cn, &mut cf) };
    rust::bfdesign(fpass, apass, fstop, astop, &mut rn, &mut rf);
    let context = format!("fpass={fpass} apass={apass} fstop={fstop} astop={astop}");
    assert_eq!(rn, cn, "{context}");
    assert!(cf.is_finite() && rf.is_finite(), "{context}: {rf} vs {cf}");
    assert!(
        (rf - cf).abs() <= 1e-5 * cf.abs(),
        "{context}: {rf} vs {cf}"
    );
    assert_eq!(rf.to_bits(), cf.to_bits(), "{context}: {rf} vs {cf}");
}

#[test]
fn design_matches_c() {
    let mut rng = XorShift::new(847);
    for _ in 0..20_000 {
        let fpass = (0.0001 + 0.4998 * rng.next_f64()) as f32;
        let mut fstop = (0.0001 + 0.4998 * rng.next_f64()) as f32;
        if fstop == fpass {
            fstop = if fpass < 0.25 { 0.49 } else { 0.01 };
        }
        let apass = (0.05 + 0.949 * rng.next_f64()) as f32;
        let astop = (apass as f64 * (0.001 + 0.998 * rng.next_f64())) as f32;
        compare_design(fpass, apass, fstop, astop);
    }
    // Valid constraints near DC/Nyquist, adjacent float frequencies, very
    // small stop-band amplitudes, and pass-band amplitudes approaching one.
    let frequencies = [
        1e-8,
        1e-5,
        0.01,
        0.1,
        0.25,
        f32::from_bits(0.25f32.to_bits() + 1),
        0.49,
        f32::from_bits(0.5f32.to_bits() - 1),
    ];
    for fpass in frequencies {
        for fstop in frequencies {
            if fpass != fstop {
                for apass in [0.01, 0.5, 0.95, f32::from_bits(1.0f32.to_bits() - 1)] {
                    for fraction in [1e-12, 0.01, 0.5, 0.9999] {
                        compare_design(fpass, apass, fstop, apass * fraction);
                    }
                }
            }
        }
    }
}

fn compare_filter(high: bool) {
    let mut rng = XorShift::new(if high { 21 } else { 22 });
    let sizes = [0, 1, 2, 3, 7, 16, 31, 64, 127, 256, 1021, 4096, 8191];
    let edges = [
        0.0,
        -0.0,
        0.5,
        -0.5,
        -0.25,
        1e-8,
        -1e-8,
        f32::from_bits(0.5f32.to_bits() - 1),
        0.25,
        f32::MIN_POSITIVE,
        f32::from_bits(1),
    ];
    for case in 0..6000 {
        let n = sizes[case % sizes.len()];
        // Cross all lengths and boundary frequencies with even/odd orders,
        // including the exact copy performed for zero poles.
        let poles = ((case / sizes.len()) % 33) as i32;
        let freq = if case % 2 == 0 {
            edges[(case / 2) % edges.len()]
        } else {
            (0.005 + 0.49 * rng.next_f64()) as f32
        };
        let amplitude = [1e-40, 1e-20, 1.0, 1e20][case % 4];
        let mut input: Vec<f32> = (0..n)
            .map(|j| match case % 17 {
                0 => 0.0,
                1 => -0.0,
                2 => {
                    if j == 0 {
                        amplitude
                    } else {
                        0.0
                    }
                }
                3 => -amplitude,
                4 => {
                    if j % 2 == 0 {
                        amplitude
                    } else {
                        -amplitude
                    }
                }
                5 => {
                    if j == n - 1 {
                        -amplitude
                    } else {
                        0.0
                    }
                }
                6 => {
                    if j % 2 == 0 {
                        0.0
                    } else {
                        -0.0
                    }
                }
                _ => (20.0 * rng.next_f64() - 10.0) as f32 * amplitude,
            })
            .collect();
        let saved_input = input.clone();
        // Leave a sentinel suffix to check that n (the input length) controls
        // both C and Rust writes even when the output buffer is larger.
        let sentinel = 12345.0f32;
        let mut c = vec![sentinel; n + 3];
        let mut r = vec![sentinel; n + 3];
        unsafe {
            if high {
                bfhighpass(poles, freq, n as i32, input.as_mut_ptr(), c.as_mut_ptr());
            } else {
                bflowpass(poles, freq, n as i32, input.as_mut_ptr(), c.as_mut_ptr());
            }
        }
        if high {
            rust::bfhighpass(poles, freq, &input, &mut r);
        } else {
            rust::bflowpass(poles, freq, &input, &mut r);
        }
        for j in 0..n {
            assert!(
                c[j].is_finite() && r[j].is_finite(),
                "nonfinite high={high} case={case} n={n} poles={poles} freq={freq} j={j}"
            );
            assert_eq!(
                r[j].to_bits(),
                c[j].to_bits(),
                "high={high} case={case} n={n} poles={poles} freq={freq} j={j}: {} vs {}",
                r[j],
                c[j]
            );
        }
        assert_eq!(&r[n..], &[sentinel; 3]);
        assert_eq!(&c[n..], &[sentinel; 3]);
        for j in 0..n {
            assert_eq!(input[j].to_bits(), saved_input[j].to_bits());
            if poles == 0 {
                assert_eq!(r[j].to_bits(), input[j].to_bits(), "zero poles must copy");
            }
        }
        // C explicitly permits p and q to alias; compare that behavior to the slice API.
        let mut inplace = input.clone();
        unsafe {
            let ptr = inplace.as_mut_ptr();
            if high {
                bfhighpass(poles, freq, n as i32, ptr, ptr);
            } else {
                bflowpass(poles, freq, n as i32, ptr, ptr);
            }
        }
        for j in 0..n {
            assert_eq!(
                r[j].to_bits(),
                inplace[j].to_bits(),
                "C aliasing case {case} j={j}"
            );
        }
    }
}

#[test]
fn highpass_matches_c() {
    compare_filter(true);
}
#[test]
fn lowpass_matches_c() {
    compare_filter(false);
}
