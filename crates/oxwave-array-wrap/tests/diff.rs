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

//! Differential tests for `cwp/lib/wrap.c`.
//! Author: Potash Corporation: Balasz Nemeth, given to CWP 2008.

use oxwave_cref::XorShift;
use std::ffi::c_void;

extern "C" {
    fn wrapArray(base: *mut c_void, nmemb: usize, size: usize, f: i32);
}

fn compare<const WIDTH: usize>(rng: &mut XorShift, n: usize, f: i32) {
    let mut input = vec![[0u8; WIDTH]; n];
    for element in &mut input {
        for byte in element {
            *byte = rng.next_u32() as u8;
        }
    }
    let mut expected = input.clone();
    unsafe { wrapArray(expected.as_mut_ptr().cast(), n, WIDTH, f) };
    oxwave_array_wrap::wrapArray(&mut input, f);
    assert_eq!(input, expected, "n={n}, width={WIDTH}, f={f}");
}

#[test]
fn wrap_array_matches_c() {
    let mut rng = XorShift::new(0xa391_52bd);
    let sizes = [1, 2, 3, 7, 16, 31, 64, 127, 256, 1024, 4093];
    for &n in &sizes {
        let shifts = [0, 1, -1, n as i32, -(n as i32), n as i32 + 1,
            -(n as i32) - 1, i32::MIN, i32::MAX];
        for &f in &shifts {
            compare::<1>(&mut rng, n, f);
            compare::<3>(&mut rng, n, f);
            compare::<8>(&mut rng, n, f);
        }
        for _ in 0..300 {
            let f = rng.next_u32() as i32;
            compare::<1>(&mut rng, n, f);
            compare::<3>(&mut rng, n, f);
            compare::<8>(&mut rng, n, f);
        }
    }
    let mut zeros = vec![0i32; 17];
    let mut expected = zeros.clone();
    unsafe { wrapArray(expected.as_mut_ptr().cast(), expected.len(), std::mem::size_of::<i32>(), -9) };
    oxwave_array_wrap::wrapArray(&mut zeros, -9);
    assert_eq!(zeros, expected);
    let mut values = [-9i32, 0, 7, i32::MIN, i32::MAX];
    let mut expected = values;
    unsafe { wrapArray(expected.as_mut_ptr().cast(), expected.len(), std::mem::size_of::<i32>(), 2) };
    oxwave_array_wrap::wrapArray(&mut values, 2);
    assert_eq!(values, expected);
}
