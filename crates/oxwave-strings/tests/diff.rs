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

//! Differential tests for cwp/lib/strstuff.c.
//! Credits: John Stockwell; John Stockwell and Jack K. Cohen; reversal based
//! on public-domain strrev by Bob Stout.

use oxwave_cref::XorShift;
use std::ffi::{c_char, c_void, CStr};

extern "C" {
    fn cwp_strdup(str: *mut c_char) -> *mut c_char;
    fn cwp_strrev(str: *mut c_char) -> *mut c_char;
    fn strchop(s: *mut c_char, t: *mut c_char);
    fn free(ptr: *mut c_void);
}

fn cases() -> Vec<Vec<u8>> {
    let mut rng = XorShift::new(0x12345678);
    let lengths = [0, 1, 2, 3, 7, 16, 31, 64, 127, 256, 1021, 4096];
    let mut cases = vec![vec![], vec![0], vec![b','], vec![b'a', 0, b',', b'z'],
        vec![b'a', b',', b'b', 0], (1..=255).collect()];
    for n in 0..6000 {
        let len = lengths[n % lengths.len()];
        let mut bytes = Vec::with_capacity(len + 1);
        for _ in 0..len { bytes.push(rng.next_u32() as u8); }
        if n % 3 == 0 {
            for b in &mut bytes { if *b == 0 { *b = 1; } }
        }
        if n % 5 == 0 && len > 0 { bytes[len / 2] = b','; }
        bytes.push(0);
        cases.push(bytes);
    }
    cases
}

fn c_input(input: &[u8]) -> Vec<u8> {
    let mut input = input.to_vec();
    input.push(0);
    input
}

#[test]
fn strdup_matches_c() {
    for input in cases() {
        let mut c_input = c_input(&input);
        let ptr = unsafe { cwp_strdup(c_input.as_mut_ptr().cast()) };
        assert!(!ptr.is_null());
        let actual = unsafe { CStr::from_ptr(ptr).to_bytes().to_vec() };
        unsafe { free(ptr.cast()) };
        assert_eq!(oxwave_strings::cwp_strdup(&input), actual);
    }
}

#[test]
fn strrev_matches_c() {
    for input in cases() {
        let mut c_input = c_input(&input);
        let ptr = unsafe { cwp_strrev(c_input.as_mut_ptr().cast()) };
        assert!(!ptr.is_null());
        let actual = unsafe { CStr::from_ptr(ptr).to_bytes().to_vec() };
        unsafe { free(ptr.cast()) };
        assert_eq!(oxwave_strings::cwp_strrev(&input), actual);
    }
}

#[test]
fn strchop_matches_c() {
    for input in cases() {
        let mut c_input = c_input(&input);
        let mut expected = vec![0xa5; input.len() + 1];
        let mut actual = expected.clone();
        unsafe { strchop(c_input.as_mut_ptr().cast(), expected.as_mut_ptr().cast()) };
        oxwave_strings::strchop(&input, &mut actual);
        assert_eq!(actual, expected);
    }
}
