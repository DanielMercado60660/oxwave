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


//! Differential tests for su/lib/tabplot.c.
//! Authors: CWP: Brian Sumner, Jack K. Cohen.

use oxwave_cref::XorShift;
use std::fs::OpenOptions;
use std::io::{Read, Seek, SeekFrom};
use std::os::fd::AsRawFd;

#[repr(C)]
struct Segy {
    header: [u8; 240],
    data: [f32; 65535],
}

extern "C" {
    fn tabplot(tp: *mut Segy, itmin: i32, itmax: i32);
    fn dup(fd: i32) -> i32;
    fn dup2(old: i32, new: i32) -> i32;
    fn close(fd: i32) -> i32;
    fn fflush(stream: *mut std::ffi::c_void) -> i32;
}

#[test]
fn tabplot_matches_c() {
    let mut rng = XorShift::new(0x12345678);
    let mut cases = Vec::new();
    for case in 0..3000 {
        let n = [1, 2, 3, 7, 16, 31, 64, 127, 256, 1021][case % 10];
        let scale = [1.0, 1e-30, 1e30, 1e-10, 1e10][case % 5];
        let mut data: Vec<f32> = (0..n).map(|_| rng.uniform_f32(-scale, scale)).collect();
        if case % 11 == 0 { data.fill(0.0); }
        if case % 13 == 0 { data[0] = -0.0; }
        let lo = rng.range_i32(0, n as i32 - 1);
        let hi = rng.range_i32(lo, n as i32 - 1);
        cases.push((data, lo, hi));
    }
    cases.push((vec![-1.0, -0.5, -0.0, 0.0, 0.5, 1.0], 0, 5));
    cases.push((vec![-f32::from_bits(1), 0.0, f32::from_bits(1)], 0, 2));
    cases.push((vec![1.0, -2.0], 1, 0));
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(format!(".tabplot-diff-{}", std::process::id()));
    let mut file = OpenOptions::new().create_new(true).read(true).write(true).open(&path).unwrap();
    let saved = unsafe { dup(1) };
    assert!(saved >= 0);
    assert_eq!(unsafe { dup2(file.as_raw_fd(), 1) }, 1);
    let mut trace = Box::new(Segy { header: [0; 240], data: [0.0; 65535] });
    for (data, lo, hi) in &cases {
        trace.data[..data.len()].copy_from_slice(data);
        unsafe { tabplot(&mut *trace, *lo, *hi) };
    }
    assert_eq!(unsafe { fflush(std::ptr::null_mut()) }, 0);
    file.seek(SeekFrom::Start(0)).unwrap();
    let mut expected = Vec::new();
    file.read_to_end(&mut expected).unwrap();
    file.set_len(0).unwrap();
    file.seek(SeekFrom::Start(0)).unwrap();
    for (data, lo, hi) in &cases {
        oxwave_trace_display::tabplot(data, *lo, *hi);
    }
    assert_eq!(oxwave_core::cio::fflush(&oxwave_core::cio::stdout()), 0);
    assert_eq!(unsafe { dup2(saved, 1) }, 1);
    unsafe { close(saved) };
    file.seek(SeekFrom::Start(0)).unwrap();
    let mut actual = Vec::new();
    file.read_to_end(&mut actual).unwrap();
    std::fs::remove_file(path).unwrap();
    assert_eq!(actual, expected);
}
