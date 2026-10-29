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

//! Differential tests for `psplot/lib/psContour.c`.
//! Author: Dave Hale, Colorado School of Mines, 06/28/89.
//! Contour labeling added by: Zhenyue Liu, August 1993.
use oxwave_cref::XorShift;
use std::ffi::{c_char, c_int, c_void, CString};
use std::fs::OpenOptions;
use std::io::{Read, Seek, SeekFrom};
use std::os::fd::AsRawFd;
use std::sync::atomic::{AtomicUsize, Ordering};
#[link(name = "psplot", kind = "static")]
extern "C" {
    fn psContour(
        c: f32,
        nx: c_int,
        x: *const f32,
        ny: c_int,
        y: *const f32,
        z: *mut f32,
        lcs: f32,
        lcf: *const c_char,
        lcc: *const c_char,
        w: *mut f32,
        nplaces: c_int,
    );
    fn dup(fd: c_int) -> c_int;
    fn dup2(fd: c_int, other: c_int) -> c_int;
    fn close(fd: c_int) -> c_int;
    fn fflush(f: *mut c_void) -> c_int;
}
fn capture(f: impl FnOnce()) -> Vec<u8> {
    // Cargo need not export CARGO_TARGET_DIR to test executables. The executable
    // itself always lives in the build directory, including with --target-dir.
    static NEXT_FILE: AtomicUsize = AtomicUsize::new(0);
    let path = std::env::current_exe().unwrap().with_file_name(format!(
        "pscontour-output-{}-{}.tmp",
        std::process::id(),
        NEXT_FILE.fetch_add(1, Ordering::Relaxed)
    ));
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    // On Unix the open descriptor remains usable; no artifact survives a panic.
    std::fs::remove_file(path).unwrap();
    unsafe {
        fflush(std::ptr::null_mut());
    }
    oxwave_core::cio::flush_all();
    let saved = unsafe { dup(1) };
    assert!(saved >= 0);
    struct RestoreStdout(c_int);
    impl Drop for RestoreStdout {
        fn drop(&mut self) {
            oxwave_core::cio::flush_all();
            unsafe {
                fflush(std::ptr::null_mut());
                assert!(dup2(self.0, 1) >= 0);
                close(self.0);
            }
        }
    }
    let restore = RestoreStdout(saved);
    assert!(unsafe { dup2(file.as_raw_fd(), 1) } >= 0);
    f();
    drop(restore);
    file.seek(SeekFrom::Start(0)).unwrap();
    let mut result = Vec::new();
    file.read_to_end(&mut result).unwrap();
    result
}

fn compare(
    description: &str,
    c: f32,
    x: &[f32],
    y: &[f32],
    z: &mut [f32],
    w: &mut [f32],
    lcs: f32,
    font: &str,
    color: &str,
    nplaces: i32,
) -> Vec<u8> {
    let mut cz = z.to_vec();
    let mut cw = w.to_vec();
    let cfont = CString::new(font).unwrap();
    let ccolor = CString::new(color).unwrap();
    let co = capture(|| unsafe {
        psContour(
            c,
            x.len() as i32,
            x.as_ptr(),
            y.len() as i32,
            y.as_ptr(),
            cz.as_mut_ptr(),
            lcs,
            cfont.as_ptr(),
            ccolor.as_ptr(),
            cw.as_mut_ptr(),
            nplaces,
        )
    });
    let ro = capture(|| oxwave_f_pscontour::psContour(c, x, y, z, lcs, font, color, w, nplaces));
    assert_eq!(
        cz.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        z.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        "z: {description}"
    );
    assert_eq!(
        cw.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        w.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        "w: {description}"
    );
    assert_eq!(co, ro, "PostScript bytes: {description}");
    ro
}

#[test]
fn pscontour_matches_c() {
    let mut rng = XorShift::new(92381);
    let sizes = [1, 2, 3, 5, 7, 8, 16, 31];
    for case in 0..3020 {
        let nx = if case >= 3000 {
            64
        } else {
            sizes[case % sizes.len()]
        };
        let ny = if case >= 3000 {
            64
        } else {
            sizes[(case / 8) % sizes.len()]
        };
        let x: Vec<f32> = (0..nx)
            .map(|i| (if case % 5 == 0 { -(i as f32) } else { i as f32 }) * 4.0 - 20.0)
            .collect();
        let y: Vec<f32> = (0..ny)
            .map(|i| (if case % 5 == 0 { -(i as f32) } else { i as f32 }) * 3.0 - 15.0)
            .collect();
        // Crossing contours in a one-dimensional grid cause out-of-bounds reads in C.
        let z: Vec<f32> = (0..nx * ny)
            .map(|_| {
                if case % 13 == 0 || nx == 1 || ny == 1 {
                    0.0
                } else {
                    rng.uniform_f32(-8.0, 8.0)
                }
            })
            .collect();
        let w: Vec<f32> = (0..nx * ny)
            .map(|_| {
                if case % 3 == 0 {
                    0.0
                } else {
                    rng.uniform_f32(-0.5, 1.0)
                }
            })
            .collect();
        let c = if case % 7 == 0 {
            0.0
        } else {
            rng.uniform_f32(-5.0, 5.0)
        };
        let lcs = match case % 4 {
            0 => 0.0,
            1 => 1.0,
            2 => 2.0,
            _ => 4.0,
        };
        let nplaces = 1 + case as i32 % 7;
        compare(
            &format!("uniform case {case}"),
            c,
            &x,
            &y,
            &mut z.clone(),
            &mut w.clone(),
            lcs,
            "Helvetica",
            "black",
            nplaces,
        );
    }

    // Both directions are independent, and spacing is deliberately nonuniform.
    // Reuse mutated z and w across levels, as callers plotting a map do.
    let sizes = [2, 3, 5, 7, 8, 16, 31, 32];
    for case in 0..2048 {
        let nx = sizes[case % sizes.len()];
        let ny = sizes[(case / sizes.len()) % sizes.len()];
        let mut x = vec![rng.uniform_f32(-200.0, 200.0); nx];
        let mut y = vec![rng.uniform_f32(-200.0, 200.0); ny];
        for i in 1..nx {
            x[i] = x[i - 1] + rng.uniform_f32(0.05, 30.0);
        }
        for i in 1..ny {
            y[i] = y[i - 1] + rng.uniform_f32(0.05, 30.0);
        }
        if case % 2 == 0 {
            x.reverse();
        }
        if case % 3 == 0 {
            y.reverse();
        }
        let scale = [1e-30_f32, 1e-6, 1.0, 1e6, 1e30][case % 5];
        let mut z: Vec<f32> = (0..nx * ny)
            .map(|k| {
                let ix = k % nx;
                let iy = k / nx;
                let a = ix as f32 / (nx - 1) as f32 - 0.5;
                let b = iy as f32 / (ny - 1) as f32 - 0.5;
                scale
                    * match case % 8 {
                        0 => a,
                        1 => b,
                        2 => a * a + b * b - 0.125, // Closed interior loops.
                        3 => a * b,                 // Saddles and exact contour vertices.
                        4 => {
                            if (ix + iy) % 2 == 0 {
                                -1.0
                            } else {
                                1.0
                            }
                        }
                        5 => {
                            if k % 2 == 0 {
                                -0.0
                            } else {
                                0.0
                            }
                        }
                        _ => rng.uniform_f32(-1.0, 1.0),
                    }
            })
            .collect();
        let mut w: Vec<f32> = (0..nx * ny)
            .map(|_| match case % 4 {
                0 => 0.0,
                1 => -1.0,
                2 => 1.0,
                _ => rng.uniform_f32(-2.0, 2.0),
            })
            .collect();
        let lcs = [-4.0, 0.0, 1.0, 1.0001, 2.5, 10.0][case % 6];
        let font = ["Helvetica", "Times-Roman", "Courier-Bold"][case % 3];
        let color = ["black", "red", "blue", "#369"][case % 4];
        let nplaces = [-1, 0, 1, 4, 6, 9][case % 6];
        for level in 0..3 {
            let c = match level {
                0 => 0.0,
                1 => scale * rng.uniform_f32(-0.5, 0.5),
                _ => scale * 2.0, // No intersections, after existing contours.
            };
            compare(
                &format!("nonuniform case {case}, level {level}"),
                c,
                &x,
                &y,
                &mut z,
                &mut w,
                lcs,
                font,
                color,
                nplaces,
            );
        }
    }

    // A single long contour must cross PATHMAX and resume after stroking.
    let x = [0.0, 100.0, 200.0];
    let y: Vec<f32> = (0..2048).map(|i| i as f32).collect();
    let mut z: Vec<f32> = (0..x.len() * y.len())
        .map(|k| (k % x.len()) as f32 - 0.5)
        .collect();
    let mut w = vec![0.0; z.len()];
    let output = compare(
        "long boundary contour",
        0.0,
        &x,
        &y,
        &mut z,
        &mut w,
        0.0,
        "Helvetica",
        "black",
        6,
    );
    assert!(String::from_utf8(output).unwrap().matches("\nS\n").count() >= 3);

    // Ensure the label drawing and mask update paths are actually exercised.
    let x: Vec<f32> = (0..65).map(|i| 10.0 * i as f32).collect();
    let y = x.clone();
    let mut z: Vec<f32> = (0..x.len() * y.len())
        .map(|k| {
            let a = (k % x.len()) as f32 - 32.0;
            let b = (k / x.len()) as f32 - 32.0;
            a * a + b * b
        })
        .collect();
    let mut w = vec![0.0; z.len()];
    let output = compare(
        "labeled closed contour",
        400.0,
        &x,
        &y,
        &mut z,
        &mut w,
        4.0,
        "Times-Roman",
        "red",
        6,
    );
    assert!(String::from_utf8(output).unwrap().contains("(400) SH\n"));
    assert!((10..55).any(|iy| (10..55).any(|ix| w[iy * x.len() + ix] > 0.0)));
}
