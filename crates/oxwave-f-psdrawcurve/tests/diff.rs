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

//! Differential tests for psplot/lib/psDrawCurve.c; credits: Brian Macy, Dave Hale and others.
use oxwave_cref::XorShift;
use std::{
    ffi::{c_char, c_void, CString},
    fs::OpenOptions,
    io::{Read, Seek, SeekFrom},
    os::fd::AsRawFd,
};

#[link(name = "psplot", kind = "static")]
extern "C" {
    fn psDrawCurve(
        x: f32,
        y: f32,
        width: f32,
        height: f32,
        x1beg: f32,
        x1end: f32,
        p1beg: f32,
        p1end: f32,
        x2beg: f32,
        x2end: f32,
        p2beg: f32,
        p2end: f32,
        x1curve: *mut f32,
        x2curve: *mut f32,
        ncurve: i32,
        color: *mut c_char,
        curvewidth: f32,
        curvedash: i32,
        style: i32,
    );
    fn dup(fd: i32) -> i32;
    fn dup2(old: i32, new: i32) -> i32;
    fn close(fd: i32) -> i32;
    fn fflush(file: *mut c_void) -> i32;
}

fn run(c: bool) {
    let mut rng = XorShift::new(0x82454589);
    let sizes = [0, 1, 2, 3, 7, 16, 31, 64, 127, 256, 1024];
    let colors = ["black", "red", "blue", "white", "gray", "#abc", "unknown"];
    for case in 0..2400 {
        let n = sizes[case % sizes.len()];
        let mut xs = rng.vec_f32(n, -1000.0, 1000.0);
        let mut ys = rng.vec_f32(n, -1000.0, 1000.0);
        if n > 0 && case % 5 == 0 {
            xs[0] = 0.0;
            ys[0] = -0.0;
        }
        let mut a = [0.0; 12];
        for v in &mut a {
            *v = rng.uniform_f32(-100.0, 100.0);
        }
        if case % 17 == 0 {
            a[4..8].fill(0.0);
        }
        if case % 19 == 0 {
            a[8..12].fill(0.0);
        }
        if case % 23 == 0 {
            a[2] = 0.0;
            a[3] = 0.0;
        }
        let color = colors[case % colors.len()];
        let width = if case % 11 == 0 {
            0.0
        } else {
            rng.uniform_f32(-5.0, 5.0)
        };
        let dash = (case % 12) as i32;
        let style = [0, 1, -1, 2][case % 4];
        if c {
            let color = CString::new(color).unwrap();
            unsafe {
                psDrawCurve(
                    a[0],
                    a[1],
                    a[2],
                    a[3],
                    a[4],
                    a[5],
                    a[6],
                    a[7],
                    a[8],
                    a[9],
                    a[10],
                    a[11],
                    xs.as_mut_ptr(),
                    ys.as_mut_ptr(),
                    n as i32,
                    color.as_ptr() as *mut c_char,
                    width,
                    dash,
                    style,
                )
            };
        } else {
            oxwave_f_psdrawcurve::psDrawCurve(
                a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], &xs, &ys,
                color, width, dash, style,
            );
        }
    }
}

#[test]
fn differential_ps_draw_curve() {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(format!("psdrawcurve-diff-{}", std::process::id()));
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
    assert!(saved >= 0);
    let mut outputs = Vec::new();
    for c in [false, true] {
        file.set_len(0).unwrap();
        file.seek(SeekFrom::Start(0)).unwrap();
        assert_eq!(unsafe { dup2(file.as_raw_fd(), 1) }, 1);
        run(c);
        oxwave_core::cio::flush_all();
        unsafe {
            fflush(std::ptr::null_mut());
            assert_eq!(dup2(saved, 1), 1);
        }
        file.seek(SeekFrom::Start(0)).unwrap();
        let mut out = Vec::new();
        file.read_to_end(&mut out).unwrap();
        assert!(!out.is_empty(), "{} produced no PostScript", if c { "C" } else { "Rust" });
        outputs.push(out);
    }
    unsafe {
        close(saved);
    }
    std::fs::remove_file(path).unwrap();
    if outputs[0] != outputs[1] {
        let at = outputs[0]
            .iter()
            .zip(&outputs[1])
            .position(|(a, b)| a != b)
            .unwrap_or(outputs[0].len().min(outputs[1].len()));
        panic!(
            "PostScript differs at byte {at}: Rust {:?}, C {:?}",
            String::from_utf8_lossy(
                &outputs[0][at.saturating_sub(100)..(at + 150).min(outputs[0].len())]
            ),
            String::from_utf8_lossy(
                &outputs[1][at.saturating_sub(100)..(at + 150).min(outputs[1].len())]
            )
        );
    }
}
