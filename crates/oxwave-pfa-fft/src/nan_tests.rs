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

use super::*;
use oxwave_cref::XorShift;

extern "C" {
    #[link_name = "pfacc"]
    fn c_pfacc(sign: i32, n: i32, z: *mut Complex);
    #[link_name = "pfacr"]
    fn c_pfacr(sign: i32, n: i32, cz: *mut Complex, rz: *mut f32);
    #[link_name = "pfarc"]
    fn c_pfarc(sign: i32, n: i32, rz: *mut f32, cz: *mut Complex);
    #[link_name = "pfamcc"]
    fn c_pfamcc(sign: i32, n: i32, nt: i32, k: i32, kt: i32, z: *mut Complex);
    #[link_name = "pfa2cr"]
    fn c_pfa2cr(sign: i32, idim: i32, n1: i32, n2: i32, cz: *mut Complex, rz: *mut f32);
    #[link_name = "pfa2rc"]
    fn c_pfa2rc(sign: i32, idim: i32, n1: i32, n2: i32, rz: *mut f32, cz: *mut Complex);
}

// Deliberately compare NaN bits, which the ordinary numerical differential
// tests accept regardless of sign/payload. Include signaling NaNs to check
// operand priority and quieting, and infinities/large finite values to create
// NaNs during the transform rather than only supplying them as input.
fn input(n: usize, case: usize) -> Vec<f32> {
    let mut rng = XorShift::new(1 + case as u64);
    (0..n)
        .map(|i| {
            let bits = match case {
                0 => 0x7fc00001 + i as u32,
                1 => 0xffc00100 + i as u32,
                2 => [0x7fc00123, 0xff800456, 0xffc00789, 0x7f800abc][i % 4],
                3 => [0x7f800000, 0xff800000, 0x7f7fffff, 0x80000000, 0][i % 5],
                4 => [0x7fc00123, 0xffc00456, 0x7fc00789, 0xffc00abc][i % 4],
                _ if i % 7 == 0 => rng.next_u32(),
                _ => return (i as f32 - 3.5) / 7.0,
            };
            f32::from_bits(bits)
        })
        .collect()
}

fn assert_bits(got: &[f32], expected: &[f32], context: &str) {
    assert_eq!(got.len(), expected.len());
    for (i, (a, b)) in got.iter().zip(expected).enumerate() {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "{context} slot={i}: {:08x} != {:08x}",
            a.to_bits(),
            b.to_bits()
        );
    }
}

#[test]
#[cfg(target_arch = "aarch64")]
fn scalar_nan_operand_order_and_negative_multipliers() {
    let a = Single(f32::from_bits(0xffc00123));
    let b = Single(f32::from_bits(0x7fc00456));
    let signaling = Single(f32::from_bits(0x7f800789));
    // ARM arithmetic chooses the first quiet NaN unless a signaling NaN is
    // present. Multiplying a NaN by -1 does not negate its sign bit.
    assert_eq!((a + b).0.to_bits(), 0xffc00123);
    assert_eq!((b + a).0.to_bits(), 0x7fc00456);
    assert_eq!((a - b).0.to_bits(), 0xffc00123);
    assert_eq!((a + signaling).0.to_bits(), 0x7fc00789);
    assert_eq!((Single(-1.0) * a).0.to_bits(), 0xffc00123);
    assert_eq!((Single(-1.0) * b).0.to_bits(), 0x7fc00456);
    assert_eq!((Double(-1.0) * a.wide()).narrow().0.to_bits(), 0xffc00123);
}

#[test]
fn complex_nan_sign_payload_and_quieting_match_c() {
    for n in [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 30, 240, 4095] {
        for sign in [-1, 1] {
            for case in 0..12 {
                let mut got = vec![Complex::default(); n];
                floats_mut(&mut got).copy_from_slice(&input(2 * n, case));
                let mut expected = got.clone();
                unsafe {
                    c_pfacc(sign, n as i32, expected.as_mut_ptr());
                }
                pfacc(sign, &mut got);
                assert_bits(
                    floats(&got),
                    floats(&expected),
                    &format!("pfacc n={n} sign={sign} case={case}"),
                );
            }
        }
    }
}

#[test]
fn real_nan_sign_payload_and_quieting_match_c() {
    for n in [
        2, 4, 6, 8, 10, 14, 16, 18, 22, 26, 32, 85, 120, 168, 8190, 8580,
    ] {
        for sign in [-1, 1] {
            for case in 0..12 {
                let mut rz = input(n, case);
                let mut got = vec![Complex::default(); (n + 3) / 2];
                let mut expected = got.clone();
                unsafe {
                    c_pfarc(sign, n as i32, rz.as_mut_ptr(), expected.as_mut_ptr());
                }
                pfarc(sign, &rz, &mut got);
                assert_bits(
                    floats(&got),
                    floats(&expected),
                    &format!("pfarc n={n} sign={sign} case={case}"),
                );

                // Independent spectra exercise the inverse twiddle, including
                // imaginary DC/Nyquist slots that C ignores for even lengths.
                floats_mut(&mut got).copy_from_slice(&input(2 * ((n + 3) / 2), case));
                let mut actual = vec![0.0; n];
                let mut expected = actual.clone();
                unsafe {
                    c_pfacr(sign, n as i32, got.as_mut_ptr(), expected.as_mut_ptr());
                }
                pfacr(sign, &got, &mut actual);
                assert_bits(
                    &actual,
                    &expected,
                    &format!("pfacr n={n} sign={sign} case={case}"),
                );
            }
        }
    }
}

#[test]
fn strided_and_2d_nan_bits_match_c() {
    for n in [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 30] {
        for sign in [-1, 1] {
            for case in 0..8 {
                let mut got = vec![Complex::default(); n * 3 + 5];
                floats_mut(&mut got).copy_from_slice(&input(2 * (n * 3 + 5), case));
                let mut expected = got.clone();
                unsafe {
                    c_pfamcc(sign, n as i32, 3, 3, 1, expected.as_mut_ptr());
                }
                pfamcc(sign, n as i32, 3, 3, 1, &mut got);
                assert_bits(floats(&got), floats(&expected), "pfamcc");

                for idim in [1, 2] {
                    let (n1, n2) = if idim == 1 { (2 * n, 3) } else { (3, 2 * n) };
                    let nc = if idim == 1 { (n + 1) * 3 } else { 3 * (n + 1) };
                    let mut rz = input(n1 * n2, case);
                    let mut got = vec![Complex::default(); nc];
                    let mut expected = got.clone();
                    unsafe {
                        c_pfa2rc(
                            sign,
                            idim,
                            n1 as i32,
                            n2 as i32,
                            rz.as_mut_ptr(),
                            expected.as_mut_ptr(),
                        );
                    }
                    pfa2rc(sign, idim, n1 as i32, n2 as i32, &rz, &mut got);
                    assert_bits(floats(&got), floats(&expected), "pfa2rc");
                    floats_mut(&mut got).copy_from_slice(&input(2 * nc, case));
                    let mut actual = vec![0.0; n1 * n2];
                    let mut expected = actual.clone();
                    unsafe {
                        c_pfa2cr(
                            sign,
                            idim,
                            n1 as i32,
                            n2 as i32,
                            got.as_mut_ptr(),
                            expected.as_mut_ptr(),
                        );
                    }
                    pfa2cr(sign, idim, n1 as i32, n2 as i32, &got, &mut actual);
                    assert_bits(&actual, &expected, "pfa2cr");
                }
            }
        }
    }
}
