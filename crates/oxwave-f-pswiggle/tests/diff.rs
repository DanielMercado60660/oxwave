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

//! Differential tests for psplot/lib/psWiggle.c; credits Dave Hale, Craig Artley, Paul Michaels.
use oxwave_cref::XorShift;
use std::{ffi::{c_char, c_void, CString}, fs::OpenOptions, io::{Read, Seek, SeekFrom}, os::fd::AsRawFd};
#[link(name = "psplot", kind = "static")]
extern "C" {
    fn psWiggle(n: i32, z: *const f32, zmin: f32, zmax: f32, zbase: f32,
        yzmin: f32, yzmax: f32, xfirst: f32, xlast: f32, ifill: i32, color: *const c_char);
    fn dup(fd: i32) -> i32;
    fn dup2(a: i32, b: i32) -> i32;
    fn close(fd: i32) -> i32;
    fn fflush(f: *mut c_void) -> i32;
}
#[test]
fn pswiggle_matches_c() {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(format!("pswiggle-diff-{}", std::process::id()));
    let mut file = OpenOptions::new().create(true).truncate(true).read(true).write(true).open(&path).unwrap();
    unsafe { fflush(std::ptr::null_mut()); }
    oxwave_core::cio::flush_all();
    let saved = unsafe { dup(1) };
    let mut outputs = Vec::new();
    for c in [false, true] {
        file.set_len(0).unwrap(); file.seek(SeekFrom::Start(0)).unwrap();
        unsafe { assert_eq!(dup2(file.as_raw_fd(), 1), 1); }
        let mut rng = XorShift::new(4123);
        let color = CString::new("red").unwrap();
        for case in 0..2400 {
            let n = [0, 1, 2, 3, 7, 16, 31, 199, 200, 201, 397, 512, 1021][case % 13];
            let mut z: Vec<f32> = (0..n).map(|_| rng.uniform_f32(-4.0, 4.0)).collect();
            if case % 7 == 0 { z.fill(0.0); }
            if case % 11 == 0 { z.fill(-0.25); }
            let zmin = if case % 19 == 0 { 0.0 } else { -1.0 };
            let zmax = if case % 19 == 0 { 0.0 } else { 1.0 };
            let base = if case % 3 == 0 { zmin } else if case % 3 == 1 { zmax } else { 0.0 };
            let a = rng.uniform_f32(-20.0, 20.0); let b = rng.uniform_f32(-20.0, 20.0);
            let x = rng.uniform_f32(-20.0, 20.0); let y = rng.uniform_f32(-20.0, 20.0);
            let ifill = ((case / 13) % 13) as i32 - 6;
            if c { unsafe { psWiggle(n as i32, z.as_ptr(), zmin, zmax, base, a, b, x, y, ifill, color.as_ptr()); } }
            else { oxwave_f_pswiggle::psWiggle(&z, zmin, zmax, base, a, b, x, y, ifill, "red"); }
        }
        oxwave_core::cio::flush_all();
        unsafe { fflush(std::ptr::null_mut()); dup2(saved, 1); }
        file.seek(SeekFrom::Start(0)).unwrap();
        let mut v = Vec::new(); file.read_to_end(&mut v).unwrap(); outputs.push(v);
    }
    unsafe { close(saved); }
    std::fs::remove_file(path).unwrap();
    if outputs[0] != outputs[1] {
        let p = outputs[0].iter().zip(&outputs[1]).position(|(a,b)| a != b).unwrap_or(outputs[0].len().min(outputs[1].len()));
        panic!("output differs at {p}: {:?} vs {:?}", String::from_utf8_lossy(&outputs[0][p.saturating_sub(40)..(p+100).min(outputs[0].len())]), String::from_utf8_lossy(&outputs[1][p.saturating_sub(40)..(p+100).min(outputs[1].len())]));
    }
}
