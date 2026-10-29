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

//! Differential tests for cwp/lib/printerplot.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89 and 12/26/89.

use oxwave_core::cio::{etmpfile, erewind, efclose, fread as rust_fread, ftello};
use oxwave_cref::XorShift;
use std::ffi::{c_char, c_int, c_void, CString};

extern "C" {
    fn pp1d(fp: *mut c_void, title: *const c_char, lx: c_int, ifx: c_int, x: *const f32);
    fn pplot1(fp: *mut c_void, title: *const c_char, nx: c_int, ax: *const f32);
    fn tmpfile() -> *mut c_void;
    fn fclose(fp: *mut c_void) -> c_int;
    fn rewind(fp: *mut c_void);
    fn ftell(fp: *mut c_void) -> i64;
    fn fread(ptr: *mut c_void, size: usize, n: usize, fp: *mut c_void) -> usize;
}

fn compare(which: usize) {
    let mut rng = XorShift::new(0x37ac + which as u64);
    let cf = unsafe { tmpfile() };
    assert!(!cf.is_null());
    let rf = etmpfile();
    let sizes = [1, 2, 3, 7, 16, 31, 64, 127, 256, 1024];
    for case in 0..3000 {
        let n = sizes[case % sizes.len()];
        let mode = (case / sizes.len()) % 10;
        let x: Vec<f32> = (0..n).map(|i| match mode {
            0 => 0.0,
            1 => -0.0,
            2 => -rng.uniform_f32(0.01, 1000.0),
            3 => rng.uniform_f32(0.01, 1000.0),
            4 => rng.uniform_f32(-1e-20, 1e-20),
            5 => rng.uniform_f32(-1e20, 1e20),
            6 => 42.0,
            7 => if i % 2 == 0 { -1.0 } else { 1.0 },
            _ => rng.uniform_f32(-1000.0, 1000.0),
        }).collect();
        let ifx = rng.range_i32(-10000, 10000);
        let title = if case % 2 == 0 { "Plot: negative and positive" } else { "" };
        let ct = CString::new(title).unwrap();
        unsafe {
            rewind(cf);
            if which == 0 { pp1d(cf, ct.as_ptr(), n as i32, ifx, x.as_ptr()); }
            else { pplot1(cf, ct.as_ptr(), n as i32, x.as_ptr()); }
        }
        erewind(&rf);
        if which == 0 { oxwave_printer_plot::pp1d(&rf, title, ifx, &x); }
        else { oxwave_printer_plot::pplot1(&rf, title, &x); }
        let c_len = unsafe { ftell(cf) } as usize;
        let r_len = ftello(&rf) as usize;
        let mut c = vec![0u8; c_len];
        let mut r = vec![0u8; r_len];
        unsafe { rewind(cf); assert_eq!(fread(c.as_mut_ptr().cast(), 1, c_len, cf), c_len); }
        erewind(&rf);
        assert_eq!(rust_fread(&mut r, &rf), r_len);
        assert_eq!(r, c, "function {which}, case {case}, size {n}, mode {mode}");
    }
    unsafe { fclose(cf); }
    efclose(&rf);
}

#[test]
fn pp1d_matches_c() { compare(0); }
#[test]
fn pplot1_matches_c() { compare(1); }

#[test]
fn tiny_ranges_produce_bounded_plots() {
    let cases: [(&[f32], &[f32]); 5] = [
        (&[-1e-39, 0.0, 1e-39], &[-1.0, 0.0, 1.0]),
        (&[0.0, 1e-39], &[0.0, 1.0]),
        (&[-1e-39, 0.0], &[-1.0, 0.0]),
        (&[-f32::from_bits(1), 0.0, f32::from_bits(1)], &[-1.0, 0.0, 1.0]),
        (&[-f32::from_bits(1), 0.0, f32::from_bits(2)], &[-1.0, 0.0, 2.0]),
    ];
    for (tiny, normal) in cases {
        for which in 0..2 {
            let rows = |values: &[f32]| {
                let fp = etmpfile();
                if which == 0 { oxwave_printer_plot::pp1d(&fp, "", 0, values); }
                else { oxwave_printer_plot::pplot1(&fp, "", values); }
                let len = ftello(&fp) as usize;
                assert!(len < 1024);
                let mut bytes = vec![0; len];
                erewind(&fp);
                assert_eq!(rust_fread(&mut bytes, &fp), len);
                efclose(&fp);
                String::from_utf8(bytes).unwrap().lines().skip(4)
                    .map(|line| line[19..].to_owned()).collect::<Vec<_>>()
            };
            assert_eq!(rows(tiny), rows(normal));
        }
    }
}
