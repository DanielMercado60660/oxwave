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

//! Differential tests of cwp/lib/rfwtva.c and cwp/lib/rfwtvaint.c.
//! Credits: Dave Hale, Colorado School of Mines, 07/01/89; Craig Artley,
//! Colorado School of Mines, 04/14/92; Tony Kocurko, Memorial University of
//! Newfoundland, Sept 1995 (sinc interpolation); Paul Michaels, Boise State
//! University, 29 December 2000 (shading); Kris Vanneste, Royal Observatory
//! of Belgium, July 2005 (amplitude-direction wiggle bits).

use oxwave_cref::XorShift;
use std::sync::Mutex;

// rfwtvaint and ints8r have mutable static storage in the original C.
static C_LOCK: Mutex<()> = Mutex::new(());

extern "C" {
    fn rfwtva(
        n: i32,
        z: *mut f32,
        zmin: f32,
        zmax: f32,
        zbase: f32,
        yzmin: i32,
        yzmax: i32,
        xfirst: i32,
        xlast: i32,
        wiggle: i32,
        nbpr: i32,
        bits: *mut u8,
        endian: i32,
    );
    fn rfwtvaint(
        n: i32,
        z: *mut f32,
        zmin: f32,
        zmax: f32,
        zbase: f32,
        yzmin: i32,
        yzmax: i32,
        xfirst: i32,
        xlast: i32,
        wiggle: i32,
        nbpr: i32,
        bits: *mut u8,
        endian: i32,
    );
}

#[derive(Debug)]
struct Case {
    z: Vec<f32>,
    zmin: f32,
    zmax: f32,
    zbase: f32,
    yzmin: i32,
    yzmax: i32,
    xfirst: i32,
    xlast: i32,
    wiggle: i32,
    nbpr: i32,
    endian: i32,
}

impl Case {
    fn rows(&self, interpolated: bool) -> usize {
        if interpolated {
            (self.xlast - self.xfirst + 1).abs() as usize
        } else {
            self.xfirst.max(self.xlast) as usize + 1
        }
    }

    fn c(&mut self, interpolated: bool, bits: &mut [u8]) {
        let f = if interpolated { rfwtvaint } else { rfwtva };
        unsafe {
            f(
                self.z.len() as i32,
                self.z.as_mut_ptr(),
                self.zmin,
                self.zmax,
                self.zbase,
                self.yzmin,
                self.yzmax,
                self.xfirst,
                self.xlast,
                self.wiggle,
                self.nbpr,
                bits.as_mut_ptr(),
                self.endian,
            );
        }
    }

    fn rust(&self, interpolated: bool, bits: &mut [u8]) {
        let f = if interpolated {
            oxwave_wiggle_raster::rfwtvaint
        } else {
            oxwave_wiggle_raster::rfwtva
        };
        f(
            &self.z,
            self.zmin,
            self.zmax,
            self.zbase,
            self.yzmin,
            self.yzmax,
            self.xfirst,
            self.xlast,
            self.wiggle,
            self.nbpr,
            bits,
            self.endian,
        );
    }
}

fn compare(case: &mut Case, interpolated: bool, rng: &mut XorShift, index: usize) {
    const GUARD: usize = 32;
    let rows = case.rows(interpolated);
    // Grey shading may advance into the next byte past the last row. Give
    // both implementations that byte explicitly; compare it too, rather
    // than hiding a C out-of-bounds access in allocator padding.
    let len = rows * case.nbpr as usize + 1;
    let mut c = vec![0xa5u8; GUARD + len + GUARD];
    for byte in &mut c[GUARD..GUARD + len] {
        *byte = match index % 3 {
            0 => 0,
            1 => rng.next_u32() as u8,
            _ => 0x55,
        };
    }
    // Start the extra byte at zero so any C write to it is observable.
    c[GUARD + len - 1] = 0;
    let original = c.clone();
    let mut rust = c.clone();
    let z_before: Vec<u32> = case.z.iter().map(|z| z.to_bits()).collect();
    case.c(interpolated, &mut c[GUARD..GUARD + len]);
    case.rust(interpolated, &mut rust[GUARD..GUARD + len]);
    assert_eq!(&c[..GUARD], &original[..GUARD], "C prefix guard");
    assert_eq!(
        &c[GUARD + len..],
        &original[GUARD + len..],
        "C suffix guard"
    );
    assert_eq!(
        z_before,
        case.z.iter().map(|z| z.to_bits()).collect::<Vec<_>>()
    );
    if let Some(byte) = c.iter().zip(&rust).position(|(a, b)| a != b) {
        panic!("case {index}, interpolated={interpolated}, byte={byte}: C={:02x}, Rust={:02x}; {case:?}",
            c[byte], rust[byte]);
    }
    for (&before, &after) in original.iter().zip(&rust) {
        assert_eq!(before & after, before, "rasterization only sets bits");
    }
    if c[GUARD + len - 1] == 0 {
        // Also exercise a slice of exactly rows * nbpr bytes whenever C did
        // not need the extra byte. Allocator padding must not hide a panic.
        let mut exact = original[GUARD..GUARD + len - 1].to_vec();
        case.rust(interpolated, &mut exact);
        assert_eq!(exact, c[GUARD..GUARD + len - 1]);
    }
}

fn random_cases(interpolated: bool) {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(if interpolated {
        0x1935_8ace
    } else {
        0x8931_cdef
    });
    let sizes = [
        1, 2, 3, 5, 7, 8, 16, 17, 31, 32, 63, 64, 127, 128, 257, 512, 1024, 4093,
    ];
    let spans = [0, 1, 2, 3, 7, 8, 16, 17, 31, 64, 127, 257, 1024, 2053];
    let shades = [-5, -1, 0, 1, 2, 3, 4, 5, 6, 10, i32::MAX];
    for index in 0..6500 {
        let n = sizes[index % sizes.len()];
        let nbpr = rng.range_i32(1, 40);
        let width = nbpr * 8;
        let mut yzmin = rng.range_i32(0, width - 1);
        let mut yzmax = rng.range_i32(0, width - 1);
        match index % 8 {
            0 => {
                yzmin = 0;
                yzmax = width - 1;
            }
            1 => {
                yzmin = width - 1;
                yzmax = 0;
            }
            2 => {
                yzmax = yzmin;
            }
            3 => {
                yzmin = -17;
                yzmax = width + 13;
            }
            4 => {
                yzmin = width + 9;
                yzmax = -9;
            }
            _ => {}
        }
        let span = spans[(index / sizes.len()) % spans.len()];
        let start = rng.range_i32(0, 9);
        let (xfirst, xlast) = if rng.next_u32() & 1 == 0 {
            (start, start + span)
        } else {
            (start + span, start)
        };
        // Scale inputs together, keeping projected coordinates in int range.
        let scale = [1e-18f32, 1e-6, 1.0, 1e6, 1e18][index % 5];
        let mut zmin = -scale * rng.uniform_f32(0.25, 3.0);
        let mut zmax = scale * rng.uniform_f32(0.25, 3.0);
        let mut z = rng.vec_f32(n, -4.0 * scale, 4.0 * scale);
        let mut zbase = rng.uniform_f32(-5.0 * scale, 5.0 * scale);
        match index % 13 {
            0 => z.fill(0.0),
            1 => z.fill(-scale),
            2 => z.fill(scale),
            3 => {
                for i in 0..n {
                    z[i] = if i % 2 == 0 {
                        -scale * 4.0
                    } else {
                        scale * 4.0
                    };
                }
            }
            4 => std::mem::swap(&mut zmin, &mut zmax),
            5 => {
                // Equal clip limits select unit yscale in C.
                zmin = -0.5;
                zmax = -0.5;
                zbase = 0.0;
                for sample in &mut z {
                    *sample /= scale;
                }
            }
            6 => zbase = zmin,
            7 => zbase = zmax,
            8 => {
                zmin = 0.0;
                zmax = 0.0;
                zbase = -0.0;
                z.fill(-0.0);
            }
            _ => {}
        }
        let mut case = Case {
            z,
            zmin,
            zmax,
            zbase,
            yzmin,
            yzmax,
            xfirst,
            xlast,
            wiggle: shades[(index / 3) % shades.len()],
            nbpr,
            endian: (index % 2) as i32,
        };
        compare(&mut case, interpolated, &mut rng, index);
    }
}

fn edge_cases(interpolated: bool) {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x2345);
    let mut index = 0;
    for z in [
        vec![0.0],
        vec![-0.0],
        vec![-3.0],
        vec![3.0],
        vec![-2.0, 2.0],
        vec![2.0, -2.0],
        vec![0.0; 17],
        vec![-2.0, 2.0, -2.0, 2.0, -2.0, 2.0, -2.0],
    ] {
        // Equal coordinates, reversed nxout=0/1, compressed and expanded
        // traces, nonzero origins, and both directions.
        for (xfirst, xlast) in [
            (0, 0),
            (5, 5),
            (1, 0),
            (2, 0),
            (0, 1),
            (0, 2),
            (0, 31),
            (31, 0),
            (4, 101),
            (101, 4),
        ] {
            for (yzmin, yzmax) in [
                (0, 0),
                (7, 7),
                (0, 7),
                (7, 0),
                (-3, 10),
                (10, -3),
                (6, 7),
                (7, 6),
            ] {
                for wiggle in [-1, 0, 1, 2, 3, 4, 5, 6] {
                    for endian in [0, 1] {
                        let mut case = Case {
                            z: z.clone(),
                            zmin: -1.0,
                            zmax: 1.0,
                            zbase: if index % 2 == 0 { 1.0 } else { -1.0 },
                            yzmin,
                            yzmax,
                            xfirst,
                            xlast,
                            wiggle,
                            nbpr: 1,
                            endian,
                        };
                        compare(&mut case, interpolated, &mut rng, index);
                        index += 1;
                    }
                }
            }
        }
    }
}

#[test]
fn rfwtva_random_matches_c() {
    random_cases(false);
}

#[test]
fn rfwtvaint_random_matches_c() {
    random_cases(true);
}

#[test]
fn rfwtva_edges_match_c() {
    edge_cases(false);
}

#[test]
fn rfwtvaint_edges_match_c() {
    edge_cases(true);
}

#[test]
fn rfwtva_empty_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    for nbpr in [1, 2, 17] {
        for (xfirst, xlast) in [(0, 0), (0, 31), (31, 0)] {
            for wiggle in [0, 1, 2, 3, 4, 5] {
                for endian in [0, 1] {
                    let mut c = vec![0xa5; 32 * nbpr as usize];
                    let mut rust = c.clone();
                    // rfwtva loads *z even with n=0: give C a valid pointer.
                    let mut dummy = -2.0;
                    unsafe {
                        rfwtva(
                            0,
                            &mut dummy,
                            -1.0,
                            1.0,
                            0.0,
                            0,
                            nbpr * 8 - 1,
                            xfirst,
                            xlast,
                            wiggle,
                            nbpr,
                            c.as_mut_ptr(),
                            endian,
                        );
                    }
                    oxwave_wiggle_raster::rfwtva(
                        &[],
                        -1.0,
                        1.0,
                        0.0,
                        0,
                        nbpr * 8 - 1,
                        xfirst,
                        xlast,
                        wiggle,
                        nbpr,
                        &mut rust,
                        endian,
                    );
                    assert_eq!(rust, c);
                    assert!(rust.iter().all(|&byte| byte == 0xa5));
                }
            }
        }
    }
    oxwave_wiggle_raster::rfwtva(&[], -1.0, 1.0, 0.0, 0, 7, 0, 0, 2, 1, &mut [], 1);
}

#[test]
fn grey_boundary_and_no_write_match_c() {
    let _guard = C_LOCK.lock().unwrap();
    for (interpolated, n) in [(false, 1), (false, 2), (true, 1), (true, 2)] {
        for endian in [0, 1] {
            let mut case = Case {
                z: vec![6.0; n],
                zmin: 0.0,
                zmax: 7.0,
                zbase: 7.0,
                yzmin: 0,
                yzmax: 7,
                xfirst: 0,
                xlast: if interpolated { 1 } else { 0 },
                wiggle: 2,
                nbpr: 1,
                endian,
            };
            let rows = case.rows(interpolated);
            let mut c = vec![0; rows + 1];
            let mut rust = c.clone();
            case.c(interpolated, &mut c);
            case.rust(interpolated, &mut rust);
            assert_eq!(rust, c);
            assert_eq!(rust[rows], if endian == 0 { 1 } else { 128 });

            // Grey fill stopping inside the last byte needs no padding.
            case.zbase = 6.0;
            c = vec![0; rows];
            rust = c.clone();
            case.c(interpolated, &mut c);
            case.rust(interpolated, &mut rust);
            assert_eq!(rust, c);

            // Valid VA-only input can perform no writes, even with no raster
            // capacity. Do not impose a conservative rows * nbpr requirement.
            case.wiggle = 0;
            case.z.fill(0.0);
            case.zbase = 7.0;
            let mut dummy = [0xa5];
            case.c(interpolated, &mut dummy);
            case.rust(interpolated, &mut []);
            assert_eq!(dummy, [0xa5]);
        }
    }
}

// These inputs cause C memory corruption or undefined signed arithmetic;
// exercise only Rust, in children so an err! exit can be checked directly.
#[test]
fn invalid_inputs_exit_cleanly() {
    const ENV: &str = "OXWAVE_WIGGLE_RASTER_INVALID_CHILD";
    if let Ok(which) = std::env::var(ENV) {
        let mut case = Case {
            z: vec![6.0, 6.0],
            zmin: 0.0,
            zmax: 7.0,
            zbase: 7.0,
            yzmin: 0,
            yzmax: 7,
            xfirst: 0,
            xlast: 0,
            wiggle: 2,
            nbpr: 1,
            endian: 1,
        };
        let interpolated = which.starts_with("sinc");
        let mut bits = vec![0; 1];
        match which.as_str() {
            "linear-grey" => {}
            "linear-grey-single" => case.z.truncate(1),
            "sinc-grey" => {
                case.xlast = 1;
                bits.resize(2, 0);
            }
            "linear-short" | "sinc-short" => {
                case.xlast = 1;
                bits.clear();
            }
            "linear-negative-row" => {
                case.xfirst = -1;
                case.xlast = -1;
            }
            "sinc-empty" => case.z.clear(),
            "linear-zero-width" | "sinc-zero-width" => case.nbpr = 0,
            "linear-negative-width" | "sinc-negative-width" => case.nbpr = -1,
            "sinc-count-add" => case.xlast = i32::MAX,
            "sinc-count-add-negative-origin" => {
                case.xfirst = i32::MIN;
                case.xlast = -1;
            }
            "sinc-count-sub-positive" => {
                case.xfirst = i32::MIN;
                case.xlast = i32::MAX;
            }
            "sinc-count-sub-negative" => {
                case.xfirst = i32::MAX;
                case.xlast = i32::MIN;
            }
            _ => panic!("unknown child case"),
        }
        case.rust(interpolated, &mut bits);
        panic!("invalid input returned instead of exiting through err!");
    }
    for (which, message) in [
        (
            "linear-grey",
            "wiggle raster: raster byte outside supplied buffer",
        ),
        (
            "linear-grey-single",
            "wiggle raster: raster byte outside supplied buffer",
        ),
        (
            "sinc-grey",
            "wiggle raster: raster byte outside supplied buffer",
        ),
        (
            "linear-short",
            "wiggle raster: raster byte outside supplied buffer",
        ),
        (
            "sinc-short",
            "wiggle raster: raster byte outside supplied buffer",
        ),
        (
            "linear-negative-row",
            "wiggle raster: raster byte outside supplied buffer",
        ),
        ("sinc-empty", "rfwtvaint: sample array must not be empty"),
        ("linear-zero-width", "rfwtva: nbpr must be positive"),
        ("sinc-zero-width", "rfwtvaint: nbpr must be positive"),
        ("linear-negative-width", "rfwtva: nbpr must be positive"),
        ("sinc-negative-width", "rfwtvaint: nbpr must be positive"),
        (
            "sinc-count-add",
            "rfwtvaint: scan-line count exceeds int range",
        ),
        (
            "sinc-count-add-negative-origin",
            "rfwtvaint: scan-line count exceeds int range",
        ),
        (
            "sinc-count-sub-positive",
            "rfwtvaint: scan-line count exceeds int range",
        ),
        (
            "sinc-count-sub-negative",
            "rfwtvaint: scan-line count exceeds int range",
        ),
    ] {
        let child = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "invalid_inputs_exit_cleanly", "--nocapture"])
            .env(ENV, which)
            .output()
            .unwrap();
        assert_eq!(child.status.code(), Some(1), "{which}: {child:?}");
        // No initargs was called, so err!'s argv0 prefix is empty.
        assert_eq!(
            child.stderr,
            format!("\n: {message}\n").as_bytes(),
            "{which}"
        );
    }
}

// Use child processes so C's fprintf and Rust's fprintf! can both be captured
// as actual stderr bytes, without a dependency or changing process-wide fds.
#[test]
fn invalid_endian_diagnostics_match_c() {
    const ENV: &str = "OXWAVE_WIGGLE_RASTER_DIAGNOSTIC_CHILD";
    if let Ok(which) = std::env::var(ENV) {
        let mut case = Case {
            z: vec![-1.0, 1.0],
            zmin: -1.0,
            zmax: 1.0,
            zbase: 0.0,
            yzmin: 0,
            yzmax: 7,
            xfirst: 0,
            xlast: 3,
            wiggle: 1,
            nbpr: 1,
            endian: -1,
        };
        for interpolated in [false, true] {
            for wiggle in [0, 1, 2, 3, 4, 5] {
                for endian in [-1, 2] {
                    case.wiggle = wiggle;
                    case.endian = endian;
                    let mut bits = vec![0x55; 16];
                    if which == "C" {
                        case.c(interpolated, &mut bits);
                    } else {
                        case.rust(interpolated, &mut bits);
                    }
                    assert_eq!(bits, vec![0x55; 16]);
                }
            }
        }
        return;
    }
    let run = |which| {
        std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "invalid_endian_diagnostics_match_c",
                "--nocapture",
            ])
            .env(ENV, which)
            .output()
            .unwrap()
    };
    let c = run("C");
    let rust = run("Rust");
    assert!(c.status.success(), "C child: {:?}", c);
    assert!(rust.status.success(), "Rust child: {:?}", rust);
    assert!(!c.stderr.is_empty());
    assert_eq!(rust.stderr, c.stderr);
    assert!(String::from_utf8(c.stderr)
        .unwrap()
        .lines()
        .all(|line| line == "endian must equal either 0 or 1"));
}
