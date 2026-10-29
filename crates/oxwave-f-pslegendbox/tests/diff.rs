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

//! Differential tests for psplot/lib/psLegendBox.c.
//! Credits: Dave Hale, Ken Larner, Torsten Schoenfelder.
use oxwave_cref::XorShift;
use oxwave_f_pslegendbox::psLegendBBox as rust_bbox;
use std::ffi::{c_char, CString};
use std::io::{Read, Seek, SeekFrom};
use std::os::fd::AsRawFd;

#[link(name = "psplot", kind = "static")]
extern "C" {
    #[link_name = "psLegendBBox"]
    fn psLegendBBox(
        x: f32,
        y: f32,
        width: f32,
        height: f32,
        font: *const c_char,
        size: f32,
        style: i32,
        bbox: *mut i32,
    );
    #[link_name = "psLegendBox"]
    fn psLegendBox(
        x: f32,
        y: f32,
        width: f32,
        height: f32,
        beg: f32,
        end: f32,
        pb: f32,
        pe: f32,
        dn: f32,
        first: f32,
        ticks: i32,
        grid: i32,
        label: *const c_char,
        font: *const c_char,
        size: f32,
        axes: *const c_char,
        grid_color: *const c_char,
        style: i32,
    );
    fn dup(fd: i32) -> i32;
    fn dup2(a: i32, b: i32) -> i32;
    fn close(fd: i32) -> i32;
    fn fflush(f: *mut std::ffi::c_void) -> i32;
}

#[test]
fn box_matches_c() {
    let path = std::path::PathBuf::from(std::env::var_os("CARGO_TARGET_DIR").unwrap_or_else(|| {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../.targets/lib-f_pslegendbox")
            .into_os_string()
    }))
        .join(format!("legend-diff-{}", std::process::id()));
    let mut file = std::fs::OpenOptions::new()
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
        unsafe {
            assert_eq!(dup2(file.as_raw_fd(), 1), 1);
        }
        let mut rng = XorShift(0x723458);
        let label = CString::new("Amplitude").unwrap();
        let font = CString::new("Helvetica").unwrap();
        let black = CString::new("black").unwrap();
        let red = CString::new("red").unwrap();
        for i in 0..2500 {
            let x = (rng.next_f64() * 300.0 - 150.0) as f32;
            let y = (rng.next_f64() * 300.0 - 150.0) as f32;
            let w = (rng.next_f64() * 250.0 + 10.0) as f32;
            let h = (rng.next_f64() * 250.0 + 10.0) as f32;
            let factor = [0.00001, 0.001, 1.0, 1000.0, 100000.0][i % 5];
            let lo = ((rng.next_u32() % 41) as i32 - 20) as f32 * factor;
            let hi = lo + (rng.next_u32() % 20 + 1) as f32 * factor;
            let (beg, end) = if i % 7 == 0 { (hi, lo) } else { (lo, hi) };
            let pad = if i % 4 == 0 { factor * 0.3 } else { 0.0 };
            let d = if i % 3 == 0 { 0.0 } else { factor };
            let style = (i % 3) as i32;
            let grid = (i % 4) as i32;
            let ticks = (i % 5 + 1) as i32;
            if c {
                unsafe {
                    psLegendBox(
                        x,
                        y,
                        w,
                        h,
                        beg,
                        end,
                        -pad,
                        pad,
                        d,
                        lo,
                        ticks,
                        grid,
                        label.as_ptr(),
                        font.as_ptr(),
                        12.0,
                        black.as_ptr(),
                        red.as_ptr(),
                        style,
                    );
                }
            } else {
                oxwave_f_pslegendbox::psLegendBox(
                    x,
                    y,
                    w,
                    h,
                    beg,
                    end,
                    -pad,
                    pad,
                    d,
                    lo,
                    ticks,
                    grid,
                    "Amplitude",
                    "Helvetica",
                    12.0,
                    "black",
                    "red",
                    style,
                );
            }
        }
        oxwave_core::cio::flush_all();
        unsafe {
            fflush(std::ptr::null_mut());
            dup2(saved, 1);
        }
        file.seek(SeekFrom::Start(0)).unwrap();
        let mut out = Vec::new();
        file.read_to_end(&mut out).unwrap();
        let clean = String::from_utf8(out)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with("test bbox_matches_c ..."))
            .collect::<Vec<_>>()
            .join("\n")
            .into_bytes();
        outputs.push(clean);
    }
    unsafe {
        close(saved);
    }
    std::fs::remove_file(path).unwrap();
    if outputs[0] != outputs[1] {
        let pos = outputs[0]
            .iter()
            .zip(&outputs[1])
            .position(|(a, b)| a != b)
            .unwrap_or(outputs[0].len().min(outputs[1].len()));
        panic!(
            "PostScript differs at byte {pos}: Rust {:?}, C {:?}",
            String::from_utf8_lossy(
                &outputs[0][pos.saturating_sub(100)..(pos + 150).min(outputs[0].len())]
            ),
            String::from_utf8_lossy(
                &outputs[1][pos.saturating_sub(100)..(pos + 150).min(outputs[1].len())]
            )
        );
    }
}

#[test]
fn bbox_matches_c() {
    let mut rng = XorShift(0x91827);
    let fonts = ["Helvetica", "Times-Roman", "Courier", "Symbol"];
    for i in 0..10000 {
        let mut rand = || (rng.next_f64() * 2000.0 - 1000.0) as f32;
        let (x, y, w, h, size) = (rand(), rand(), rand(), rand(), rand());
        let font = fonts[i % fonts.len()];
        let cfont = CString::new(font).unwrap();
        let style = (i % 5) as i32;
        let mut a = [0; 4];
        let mut b = [0; 4];
        rust_bbox(x, y, w, h, font, size, style, &mut a);
        unsafe { psLegendBBox(x, y, w, h, cfont.as_ptr(), size, style, b.as_mut_ptr()) };
        assert_eq!(a, b, "case {i}");
    }
}
