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

//! Differential tests for cwp/lib/zasc.c.
//! Author: Stew Levin of Mobil, 1997; translated by f2c.

use oxwave_cref::XorShift;
use std::ffi::{c_char, c_int};

extern "C" {
    fn zasc(input: *mut c_char, output: *mut c_char, nchar: c_int) -> c_int;
}

#[test]
fn zasc_matches_c() {
    fn check(input: &[u8]) {
        let mut c_input = input.to_vec();
        let mut expected = vec![0xa5; input.len() + 5];
        let mut actual = expected.clone();
        let c_result = unsafe {
            zasc(c_input.as_mut_ptr().cast(), expected.as_mut_ptr().cast(), input.len() as c_int)
        };
        let result = oxwave_character_encoding::zasc(input, &mut actual);
        assert_eq!(result, c_result);
        assert_eq!(actual, expected, "input: {input:?}");
        assert_eq!(c_input, input);
    }
    check(&[]);
    check(&(0..=255).collect::<Vec<u8>>());
    for byte in 0..=255 {
        check(&[byte]);
    }
    let mut rng = XorShift::new(0x1234_5678);
    for n in [1, 2, 3, 7, 17, 31, 64, 127, 256, 1024, 4093, 8192] {
        check(&vec![0; n]);
        check(&vec![255; n]);
        for _ in 0..300 {
            let input: Vec<u8> = (0..n).map(|_| rng.next_u32() as u8).collect();
            check(&input);
        }
    }
    // A negative C count executes no iterations.
    let mut input = [0u8];
    let mut output = [0xa5u8];
    assert_eq!(unsafe { zasc(input.as_mut_ptr().cast(), output.as_mut_ptr().cast(), -1) }, 0);
    assert_eq!(output, [0xa5]);
}
