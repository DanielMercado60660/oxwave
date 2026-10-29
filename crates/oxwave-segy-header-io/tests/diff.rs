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


//! Differential tests for su/lib/fgethdr.c.
//! Authors: zhiming li and j. dulac, unocal; modified for CWP/SU: R. Beardsley.

use oxwave_core::cio::{efwrite, erewind, etmpfile, ftello};
use oxwave_core::segy::{Bhed, EBCBYTES};
use oxwave_cref::XorShift;
use std::ffi::{c_char, c_int, c_long, c_void};

extern "C" {
    fn fgethdr(fp: *mut c_void, chdr: *mut c_char, bhdr: *mut Bhed);
    fn tmpfile() -> *mut c_void;
    fn fwrite(ptr: *const c_void, size: usize, count: usize, fp: *mut c_void) -> usize;
    fn rewind(fp: *mut c_void);
    fn ftell(fp: *mut c_void) -> c_long;
    fn fclose(fp: *mut c_void) -> c_int;
}

#[test]
fn fgethdr_matches_c() {
    let mut rng = XorShift::new(0x539a_6d23);
    let lengths = [3600, 3600, 0, 1, 2, 3, 17, 1024, 3199, 3200, 3201, 3202, 3571, 3599, 3600, 3601, 4096, 8192];
    for case in 0..2200 {
        let len = if case < lengths.len() { lengths[case] } else { (rng.next_u32() as usize) % 5000 };
        let input: Vec<u8> = match case {
            0 => vec![0; len],
            1 => vec![0xff; len],
            _ => (0..len).map(|_| rng.next_u32() as u8).collect(),
        };
        let mut ctext = [0u8; EBCBYTES + 7];
        for b in &mut ctext { *b = rng.next_u32() as u8; }
        let mut rtext = ctext;
        let mut cbinary = Bhed::zeroed();
        for b in cbinary.as_bytes_mut() { *b = rng.next_u32() as u8; }
        let mut rbinary = cbinary;
        let rust_file = etmpfile();
        efwrite(&input, &rust_file);
        erewind(&rust_file);
        unsafe {
            let c_file = tmpfile();
            assert!(!c_file.is_null());
            assert_eq!(fwrite(input.as_ptr().cast(), 1, input.len(), c_file), input.len());
            rewind(c_file);
            fgethdr(c_file, ctext.as_mut_ptr().cast(), &mut cbinary);
            oxwave_segy_header_io::fgethdr(&rust_file, &mut rtext, &mut rbinary);
            assert_eq!(rtext, ctext, "text case {case}, input length {len}");
            assert_eq!(rbinary.as_bytes(), cbinary.as_bytes(), "binary case {case}, input length {len}");
            assert_eq!(ftello(&rust_file), ftell(c_file) as i64, "position case {case}");
            assert_eq!(fclose(c_file), 0);
        }
    }
}
