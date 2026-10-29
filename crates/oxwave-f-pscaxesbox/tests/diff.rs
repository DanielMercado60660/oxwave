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

//! Differential test for psplot/lib/psCAxesBox.c; Zhiming Li, Dave Hale, Craig Artley.
use oxwave_cref::XorShift;
use std::{
    ffi::{c_char, c_void, CString},
    fs::OpenOptions,
    io::{Read, Seek, SeekFrom},
    os::fd::AsRawFd,
};

#[link(name = "psplot", kind = "static")]
extern "C" {
    fn psCubeAxesBox(
        x: f32,
        y: f32,
        s1: f32,
        s2: f32,
        s3: f32,
        a: f32,
        b1: f32,
        e1: f32,
        p1: f32,
        q1: f32,
        d1: f32,
        f1: f32,
        n1: i32,
        g1: i32,
        l1: *const c_char,
        b2: f32,
        e2: f32,
        p2: f32,
        q2: f32,
        d2: f32,
        f2: f32,
        n2: i32,
        g2: i32,
        l2: *const c_char,
        b3: f32,
        e3: f32,
        p3: f32,
        q3: f32,
        d3: f32,
        f3: f32,
        n3: i32,
        g3: i32,
        l3: *const c_char,
        lf: *const c_char,
        ls: f32,
        t: *const c_char,
        tf: *const c_char,
        ts: f32,
        tc: *const c_char,
        ac: *const c_char,
        gc: *const c_char,
    );
    fn dup(fd: i32) -> i32;
    fn dup2(a: i32, b: i32) -> i32;
    fn close(fd: i32) -> i32;
    fn fflush(f: *mut c_void) -> i32;
}

#[test]
fn cube_axes_matches_c() {
    let path = std::env::var_os("CARGO_TARGET_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")))
        .join(format!("cube-diff-{}", std::process::id()));
    let mut file = OpenOptions::new()
        .create(true)
        .truncate(true)
        .read(true)
        .write(true)
        .open(&path)
        .unwrap();
    unsafe {
        fflush(std::ptr::null_mut());
    }
    oxwave_core::cio::flush_all();
    let saved = unsafe { dup(1) };
    let mut outputs = Vec::new();
    for c in [false, true] {
        file.set_len(0).unwrap();
        file.seek(SeekFrom::Start(0)).unwrap();
        unsafe { assert_eq!(dup2(file.as_raw_fd(), 1), 1) };
        let mut rng = XorShift::new(187);
        let font = CString::new("Helvetica").unwrap();
        let label = CString::new("Axis").unwrap();
        let title = CString::new("Title").unwrap();
        let color = CString::new("black").unwrap();
        for i in 0..2500 {
            let x = (rng.next_u32() % 2000) as f32 / 10.0 - 100.0;
            let y = (rng.next_u32() % 2000) as f32 / 10.0 - 100.0;
            let sizes = [1.0, 2.0, 7.0, 32.0, 127.0, 512.0];
            let s1 = sizes[i % 6];
            let s2 = sizes[(i / 6) % 6];
            let s3 = sizes[(i / 36) % 6];
            let angle = (rng.next_u32() % 1200) as f32 / 1000.0 - 0.6;
            let beg = -((rng.next_u32() % 100) as f32) / 10.0;
            let end = (rng.next_u32() % 100 + 1) as f32 / 10.0;
            let d = if i % 2 == 0 { 0.0 } else { 1.0 };
            let grid = (i % 4) as i32;
            let nt = [1, 2, 3, 7][i % 4];
            let pad = if i % 3 == 0 { 0.25 } else { 0.0 };
            if c {
                unsafe {
                    psCubeAxesBox(
                        x,
                        y,
                        s1,
                        s2,
                        s3,
                        angle,
                        beg,
                        end,
                        -pad,
                        pad,
                        d,
                        0.0,
                        nt,
                        grid,
                        label.as_ptr(),
                        end,
                        beg,
                        pad,
                        -pad,
                        d,
                        0.0,
                        nt,
                        (grid + 1) % 4,
                        label.as_ptr(),
                        beg,
                        end,
                        -pad,
                        pad,
                        d,
                        0.0,
                        nt,
                        (grid + 2) % 4,
                        label.as_ptr(),
                        font.as_ptr(),
                        10.0,
                        title.as_ptr(),
                        font.as_ptr(),
                        12.0,
                        color.as_ptr(),
                        color.as_ptr(),
                        color.as_ptr(),
                    );
                }
            } else {
                oxwave_f_pscaxesbox::psCubeAxesBox(
                    x,
                    y,
                    s1,
                    s2,
                    s3,
                    angle,
                    beg,
                    end,
                    -pad,
                    pad,
                    d,
                    0.0,
                    nt,
                    grid,
                    "Axis",
                    end,
                    beg,
                    pad,
                    -pad,
                    d,
                    0.0,
                    nt,
                    (grid + 1) % 4,
                    "Axis",
                    beg,
                    end,
                    -pad,
                    pad,
                    d,
                    0.0,
                    nt,
                    (grid + 2) % 4,
                    "Axis",
                    "Helvetica",
                    10.0,
                    "Title",
                    "Helvetica",
                    12.0,
                    "black",
                    "black",
                    "black",
                );
            }
        }
        oxwave_core::cio::flush_all();
        unsafe {
            fflush(std::ptr::null_mut());
            dup2(saved, 1);
        }
        file.seek(SeekFrom::Start(0)).unwrap();
        let mut bytes = Vec::new();
        file.read_to_end(&mut bytes).unwrap();
        outputs.push(bytes);
    }
    unsafe {
        close(saved);
    }
    std::fs::remove_file(path).unwrap();
    let a = String::from_utf8(outputs.remove(0)).unwrap();
    let b = String::from_utf8(outputs.remove(0)).unwrap();
    let aa: Vec<_> = a.split_whitespace().collect();
    let bb: Vec<_> = b.split_whitespace().collect();
    assert_eq!(aa.len(), bb.len());
    for (i, (x, y)) in aa.iter().zip(&bb).enumerate() {
        if x != y {
            match (x.parse::<f64>(), y.parse::<f64>()) {
                (Ok(x), Ok(y)) => assert!(
                    (x - y).abs() <= 1e-5 * x.abs().max(y.abs()),
                    "token {i}: {x} vs {y}"
                ),
                _ => panic!("token {i}: {x} vs {y}"),
            }
        }
    }
}
