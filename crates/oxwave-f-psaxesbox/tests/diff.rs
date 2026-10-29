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

//! Differential tests for psplot/lib/psAxesBox.c.
//! Author: Dave Hale, Colorado School of Mines, 06/27/89.
//! Modified: Ken Larner, Colorado School of Mines, 08/30/90;
//! Dave Hale, Advance Geophysical, 10/18/92;
//! Morten Wendell Pedersen, Aarhus University, 23/3-97.

use oxwave_cref::XorShift;
use std::{
    ffi::{c_char, c_void, CString},
    fs::{File, OpenOptions},
    io::{Read, Seek, SeekFrom},
    os::fd::AsRawFd,
    path::PathBuf,
    process::Command,
};

#[link(name = "psplot", kind = "static")]
extern "C" {
    fn psAxesBBox(
        x: f32,
        y: f32,
        width: f32,
        height: f32,
        label_font: *const c_char,
        label_size: f32,
        title_font: *const c_char,
        title_size: f32,
        style: i32,
        bbox: *mut i32,
    );
    fn psAxesBox(
        x: f32,
        y: f32,
        width: f32,
        height: f32,
        x1_beg: f32,
        x1_end: f32,
        p1_beg: f32,
        p1_end: f32,
        d1_num: f32,
        f1_num: f32,
        n1_tic: i32,
        grid1: i32,
        label1: *const c_char,
        x2_beg: f32,
        x2_end: f32,
        p2_beg: f32,
        p2_end: f32,
        d2_num: f32,
        f2_num: f32,
        n2_tic: i32,
        grid2: i32,
        label2: *const c_char,
        label_font: *const c_char,
        label_size: f32,
        title: *const c_char,
        title_font: *const c_char,
        title_size: f32,
        title_color: *const c_char,
        axes_color: *const c_char,
        grid_color: *const c_char,
        ticwidth: f32,
        axeswidth: f32,
        gridwidth: f32,
        style: i32,
    );
}

#[test]
fn bbox_matches_c() {
    let fonts = FONTS;
    let mut rng = XorShift::new(12345);
    for i in 0..10000 {
        let x = if i % 11 == 0 {
            0.0
        } else {
            rng.uniform_f32(-10000.0, 10000.0)
        };
        let y = if i % 13 == 0 {
            -0.0
        } else {
            rng.uniform_f32(-10000.0, 10000.0)
        };
        let sizes = [
            0.0, 1.0, 2.0, 3.0, 7.0, 16.0, 31.0, 64.0, 127.0, 256.0, 1024.0, 65536.0,
        ];
        let w = if i % 3 == 0 {
            sizes[(i / 3) % sizes.len()]
        } else {
            rng.uniform_f32(-1000.0, 1000.0)
        };
        let h = if i % 5 == 0 {
            -sizes[(i / 5) % sizes.len()]
        } else {
            rng.uniform_f32(-1000.0, 1000.0)
        };
        let ls = if i % 7 == 0 {
            0.0
        } else {
            rng.uniform_f32(-100.0, 100.0)
        };
        let ts = rng.uniform_f32(-100.0, 100.0);
        let lf = fonts[i % fonts.len()];
        let tf = fonts[(i / 4) % fonts.len()];
        let style = [0, 1, -1, 2, i32::MIN, i32::MAX][i % 6];
        let mut expected = [0; 4];
        let mut actual = [0; 4];
        let lc = CString::new(lf).unwrap();
        let tc = CString::new(tf).unwrap();
        unsafe {
            psAxesBBox(
                x,
                y,
                w,
                h,
                lc.as_ptr(),
                ls,
                tc.as_ptr(),
                ts,
                style,
                expected.as_mut_ptr(),
            )
        };
        oxwave_f_psaxesbox::psAxesBBox(x, y, w, h, lf, ls, tf, ts, style, &mut actual);
        assert_eq!(actual, expected, "case {i}");
    }
}

const FONTS: &[&str] = &[
    "AvantGarde-Book",
    "Bookman-DemiItalic",
    "Courier",
    "Courier-BoldOblique",
    "Helvetica",
    "Helvetica-Bold",
    "Helvetica-Narrow-Oblique",
    "NewCenturySchlbk-Roman",
    "Palatino-Roman",
    "SanSerif-Roman",
    "Symbol",
    "Times-BoldItalic",
    "Times-Roman",
    "ZapfChancery-MediumItalic",
    "Unknown",
];

extern "C" {
    fn dup(fd: i32) -> i32;
    fn dup2(from: i32, to: i32) -> i32;
    fn close(fd: i32) -> i32;
    fn fflush(stream: *mut c_void) -> i32;
}

// Restore stdout even if a test panics.
struct Redirect(i32);
impl Redirect {
    fn to(file: &File) -> Self {
        oxwave_core::cio::flush_all();
        assert_eq!(unsafe { fflush(std::ptr::null_mut()) }, 0);
        let saved = unsafe { dup(1) };
        assert!(saved >= 0);
        assert_eq!(unsafe { dup2(file.as_raw_fd(), 1) }, 1);
        Self(saved)
    }
}
impl Drop for Redirect {
    fn drop(&mut self) {
        oxwave_core::cio::flush_all();
        unsafe {
            fflush(std::ptr::null_mut());
            dup2(self.0, 1);
            close(self.0);
        }
    }
}

struct Capture {
    file: File,
    path: PathBuf,
}
impl Capture {
    fn new() -> Self {
        let build = std::env::var_os("CARGO_TARGET_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| {
                PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../.targets/lib-f_psaxesbox")
            });
        std::fs::create_dir_all(&build).unwrap();
        let path = build.join(format!("psaxesbox-diff-{}.ps", std::process::id()));
        let file = OpenOptions::new()
            .create_new(true)
            .read(true)
            .write(true)
            .open(&path)
            .unwrap();
        Self { file, path }
    }

    fn output(&mut self, run: impl FnOnce()) -> Vec<u8> {
        self.file.set_len(0).unwrap();
        self.file.seek(SeekFrom::Start(0)).unwrap();
        {
            let _redirect = Redirect::to(&self.file);
            run();
        }
        self.file.seek(SeekFrom::Start(0)).unwrap();
        let mut bytes = Vec::new();
        self.file.read_to_end(&mut bytes).unwrap();
        bytes
    }
}
impl Drop for Capture {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

#[derive(Clone, Copy, Debug)]
struct Axis {
    beg: f32,
    end: f32,
    pbeg: f32,
    pend: f32,
    dnum: f32,
    fnum: f32,
    ntic: i32,
    grid: i32,
}
impl Axis {
    fn random(rng: &mut XorShift, i: usize) -> Self {
        let magnitudes = [
            1e-12, 1e-8, 1e-5, 1e-4, 0.1, 1.0, 10.0, 1e3, 1e4, 1e6, 1e10, 1e15,
        ];
        let magnitude = magnitudes[i % magnitudes.len()];
        let lo = rng.uniform_f32(-4.0, 4.0) * magnitude;
        let span = rng.uniform_f32(0.25, 5.0) * magnitude;
        let hi = lo + span;
        let (beg, end) = if i % 3 == 0 { (hi, lo) } else { (lo, hi) };
        let steps = [1.0, 2.0, 3.0, 7.0, 16.0, 31.0][(i / 3) % 6];
        let dnum = span / steps;
        let mut a = Self {
            beg,
            end,
            pbeg: rng.uniform_f32(-0.2, 0.2) * span,
            pend: rng.uniform_f32(-0.2, 0.2) * span,
            dnum: if i % 4 == 0 { 0.0 } else { dnum },
            fnum: lo + rng.uniform_f32(-0.5, 0.5) * dnum,
            ntic: [0, 1, 2, 3, 5, 7, 16][(i / 7) % 7],
            grid: [0, 1, 2, 3, -1, 4][(i / 11) % 6],
        };
        match i % 37 {
            0 => {
                // Equal limits, including zero, with and without padding.
                a.beg = if i % 2 == 0 { 0.0 } else { lo };
                a.end = a.beg;
                a.fnum = a.beg;
                // scaxis returns increment 1 for equal limits. At large
                // coordinates that cannot advance, so specify a finite step.
                a.dnum = if a.beg.abs() > 4096.0 {
                    a.beg.abs() * 0.01
                } else {
                    0.0
                };
                a.pbeg = 0.0;
                a.pend = if i % 2 == 0 { 0.0 } else { span };
            }
            1 => {
                // Crossing zero exercises near-zero label formatting.
                a.beg = -magnitude;
                a.end = magnitude;
                a.dnum = magnitude / 2.0;
                a.fnum = -magnitude;
                a.pbeg = 0.0;
                a.pend = 0.0;
            }
            2 => {
                // No numbered ticks; axis 2 inherits C's nexpmax.
                a.dnum = dnum;
                a.fnum = hi + 3.0 * span;
            }
            3 => {
                // Negative increments with empty, terminating loops.
                a.dnum = -dnum;
                a.fnum = hi + 3.0 * span;
                a.ntic = 0;
            }
            4 => {
                // Negative subdivision count, also with empty loops.
                a.dnum = dnum;
                a.fnum = hi + 3.0 * span;
                a.ntic = -7;
            }
            5 => {
                // Large subdivision count.
                a.dnum = span;
                a.fnum = lo;
                a.ntic = 257;
            }
            _ => {}
        }
        a
    }
}

#[derive(Debug)]
struct Plot {
    x: f32,
    y: f32,
    width: f32,
    height: f32,
    a: Axis,
    b: Axis,
    label_size: f32,
    title_size: f32,
    ticwidth: f32,
    axeswidth: f32,
    gridwidth: f32,
    style: i32,
    // label1, label2, labelFont, title, titleFont, titleColor, axesColor, gridColor
    text: [&'static str; 8],
}
impl Plot {
    fn random(rng: &mut XorShift, i: usize) -> Self {
        let sizes = [
            0.0, 1.0, 2.0, 3.0, 7.0, 16.0, 31.0, 64.0, 127.0, 256.0, 1024.0, 4096.0,
        ];
        let width = if i % 2 == 0 {
            sizes[(i / 2) % sizes.len()]
        } else {
            rng.uniform_f32(-720.0, 720.0)
        };
        let height = if i % 3 == 0 {
            -sizes[(i / 3) % sizes.len()]
        } else {
            rng.uniform_f32(-720.0, 720.0)
        };
        let mut a = Axis::random(rng, i);
        let mut b = Axis::random(rng, i * 17 + 13);
        let label_size = match i % 19 {
            0 => 0.0,
            1 => -12.0,
            _ => rng.uniform_f32(6.0, 30.0),
        };
        // Zero label width would ask scaxis for an unbounded tick count.
        // Still test zero-size fonts, using explicit finite increments.
        if label_size == 0.0 {
            for axis in [&mut a, &mut b] {
                if axis.dnum == 0.0 {
                    axis.dnum = if axis.beg == axis.end {
                        // Advance even for equal limits of large magnitude.
                        (axis.beg.abs() * 0.0001).max(1.0)
                    } else {
                        (axis.end - axis.beg).abs()
                    };
                }
            }
        }
        let labels = [
            "",
            "Time (s)",
            "Offset / km",
            "x\\y",
            "negative - values",
            "a(b)c",
            "depth",
        ];
        let colors = [
            "black",
            "red",
            "green",
            "blue",
            "gray",
            "cyan",
            "goldenrod",
            "Black",
            "#f00",
            "#123456",
            "#123456789",
            "#0123456789ab",
            "unrecognized",
        ];
        Self {
            x: if i % 13 == 0 {
                -0.0
            } else {
                rng.uniform_f32(-1000.0, 1000.0)
            },
            y: if i % 17 == 0 {
                0.0
            } else {
                rng.uniform_f32(-1000.0, 1000.0)
            },
            width,
            height,
            a,
            b,
            label_size,
            title_size: if i % 17 == 0 {
                0.0
            } else {
                rng.uniform_f32(-24.0, 36.0)
            },
            ticwidth: if i % 5 == 0 {
                0.0
            } else {
                rng.uniform_f32(-2.0, 4.0)
            },
            axeswidth: rng.uniform_f32(-2.0, 4.0),
            gridwidth: rng.uniform_f32(-2.0, 4.0),
            style: [0, 1, -1, 2, i32::MIN, i32::MAX][(i / 5) % 6],
            text: [
                labels[i % labels.len()],
                labels[(i / 3) % labels.len()],
                FONTS[i % FONTS.len()],
                labels[(i / 7) % labels.len()],
                FONTS[(i / FONTS.len()) % FONTS.len()],
                colors[i % colors.len()],
                colors[(i / 13) % colors.len()],
                colors[(i / 169) % colors.len()],
            ],
        }
    }

    fn rust(&self) {
        let (a, b, t) = (self.a, self.b, self.text);
        oxwave_f_psaxesbox::psAxesBox(
            self.x,
            self.y,
            self.width,
            self.height,
            a.beg,
            a.end,
            a.pbeg,
            a.pend,
            a.dnum,
            a.fnum,
            a.ntic,
            a.grid,
            t[0],
            b.beg,
            b.end,
            b.pbeg,
            b.pend,
            b.dnum,
            b.fnum,
            b.ntic,
            b.grid,
            t[1],
            t[2],
            self.label_size,
            t[3],
            t[4],
            self.title_size,
            t[5],
            t[6],
            t[7],
            self.ticwidth,
            self.axeswidth,
            self.gridwidth,
            self.style,
        );
    }

    fn c(&self) {
        let (a, b) = (self.a, self.b);
        // setcolor's C implementation lowercases in place, hence writable buffers.
        let mut text = self
            .text
            .map(|s| CString::new(s).unwrap().into_bytes_with_nul());
        let t = text.each_mut().map(|s| s.as_mut_ptr().cast::<c_char>());
        unsafe {
            psAxesBox(
                self.x,
                self.y,
                self.width,
                self.height,
                a.beg,
                a.end,
                a.pbeg,
                a.pend,
                a.dnum,
                a.fnum,
                a.ntic,
                a.grid,
                t[0],
                b.beg,
                b.end,
                b.pbeg,
                b.pend,
                b.dnum,
                b.fnum,
                b.ntic,
                b.grid,
                t[1],
                t[2],
                self.label_size,
                t[3],
                t[4],
                self.title_size,
                t[5],
                t[6],
                t[7],
                self.ticwidth,
                self.axeswidth,
                self.gridwidth,
                self.style,
            );
        }
    }
}

#[test]
fn axes_box_matches_c() {
    // Descriptor redirection is process-wide. Isolate it from other tests and
    // libtest's progress output, running this same test alone in the child.
    if std::env::var_os("OXWAVE_PSAXESBOX_DIFF_CHILD").is_none() {
        let result = Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "axes_box_matches_c", "--nocapture"])
            .env("OXWAVE_PSAXESBOX_DIFF_CHILD", "1")
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "child failed: {}
{}",
            String::from_utf8_lossy(&result.stdout),
            String::from_utf8_lossy(&result.stderr)
        );
        return;
    }
    let mut rng = XorShift::new(0x7861_7865_7362_6f78);
    let mut capture = Capture::new();
    for i in 0..5000 {
        let plot = Plot::random(&mut rng, i);
        let expected = capture.output(|| plot.c());
        let actual = capture.output(|| plot.rust());
        assert!(!expected.is_empty());
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(actual.len().min(expected.len()));
            let start = first.saturating_sub(100);
            panic!(
                "case {i}, first different byte {first}, Rust length {}, C length {}
{plot:#?}
Rust: {}
C: {}",
                actual.len(),
                expected.len(),
                String::from_utf8_lossy(&actual[start..(first + 200).min(actual.len())]),
                String::from_utf8_lossy(&expected[start..(first + 200).min(expected.len())])
            );
        }
    }
}
