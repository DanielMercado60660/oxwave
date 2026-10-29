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

//! Differential tests against cwp/lib/intl2b.c and cwp/lib/intlirr2b.c.
//! Author: Dave Hale, Colorado School of Mines, c. 1989-1991;
//! 07/02/89, modified 05/30/90 and 06/01/91.

use oxwave_cref::XorShift;
use std::sync::Mutex;

static C_LOCK: Mutex<()> = Mutex::new(());

mod c {
    extern "C" {
        pub fn intl2b(
            nxin: i32,
            dxin: f32,
            fxin: f32,
            nyin: i32,
            dyin: f32,
            fyin: f32,
            zin: *mut u8,
            nxout: i32,
            dxout: f32,
            fxout: f32,
            nyout: i32,
            dyout: f32,
            fyout: f32,
            zout: *mut u8,
        );
        pub fn intlirr2b(
            nxin: i32,
            xin: *mut f32,
            nyin: i32,
            dyin: f32,
            fyin: f32,
            zin: *mut u8,
            nxout: i32,
            dxout: f32,
            fxout: f32,
            nyout: i32,
            dyout: f32,
            fyout: f32,
            zout: *mut u8,
        );
    }
}

#[derive(Clone, Copy, Debug)]
struct Case {
    nx: i32,
    ny: i32,
    ox: i32,
    oy: i32,
    dx: f32,
    fx: f32,
    dy: f32,
    fy: f32,
    odx: f32,
    ofx: f32,
    ody: f32,
    ofy: f32,
}

impl Default for Case {
    fn default() -> Self {
        Self {
            nx: 2,
            ny: 2,
            ox: 7,
            oy: 7,
            dx: 1.0,
            fx: 0.0,
            dy: 1.0,
            fy: 0.0,
            odx: 0.25,
            ofx: -0.25,
            ody: 0.25,
            ofy: -0.25,
        }
    }
}

// All output bytes, including the write guard, must match exactly.
// Singleton C accesses need a byte before zin and one extra row after zin:
// the missing neighbors are read even though their table weights are zero.
// Keep those accesses inside a real initialized allocation. Rust receives
// only the logical grid, with no padding.
fn check(case: Case, xin: Option<&[f32]>, zin: &[u8]) -> Vec<u8> {
    assert_eq!(zin.len(), (case.nx * case.ny) as usize);
    let prefix = 8;
    let mut cinput = vec![0x6du8; prefix + zin.len() + case.nx as usize + 8];
    cinput[prefix..prefix + zin.len()].copy_from_slice(zin);
    let saved = cinput.clone();
    let nout = (case.ox * case.oy) as usize;
    let mut expected = vec![0xa5u8; nout + 16];
    let mut actual = expected.clone();
    if let Some(xin) = xin {
        assert_eq!(xin.len(), case.nx as usize);
        let mut cxin = xin.to_vec();
        // SAFETY: both grids and all C singleton neighbor reads are in bounds
        // of their allocations. Coordinates are ordered and dimensions valid.
        unsafe {
            c::intlirr2b(
                case.nx,
                cxin.as_mut_ptr(),
                case.ny,
                case.dy,
                case.fy,
                cinput.as_mut_ptr().add(prefix),
                case.ox,
                case.odx,
                case.ofx,
                case.oy,
                case.ody,
                case.ofy,
                expected.as_mut_ptr(),
            );
        }
        assert_eq!(cxin, xin, "C modified xin: {case:?}");
        oxwave_byte_interp::intlirr2b(
            xin,
            case.ny,
            case.dy,
            case.fy,
            zin,
            case.ox,
            case.odx,
            case.ofx,
            case.oy,
            case.ody,
            case.ofy,
            &mut actual,
        );
    } else {
        // SAFETY: see above. The C functions have the same grid accesses.
        unsafe {
            c::intl2b(
                case.nx,
                case.dx,
                case.fx,
                case.ny,
                case.dy,
                case.fy,
                cinput.as_mut_ptr().add(prefix),
                case.ox,
                case.odx,
                case.ofx,
                case.oy,
                case.ody,
                case.ofy,
                expected.as_mut_ptr(),
            );
        }
        oxwave_byte_interp::intl2b(
            case.nx,
            case.dx,
            case.fx,
            case.ny,
            case.dy,
            case.fy,
            zin,
            case.ox,
            case.odx,
            case.ofx,
            case.oy,
            case.ody,
            case.ofy,
            &mut actual,
        );
    }
    assert_eq!(saved, cinput, "C modified zin: {case:?}");
    if actual != expected {
        let i = actual
            .iter()
            .zip(&expected)
            .position(|(a, b)| a != b)
            .unwrap();
        panic!(
            "byte {i}: Rust {}, C {}; irregular={}, {case:?}",
            actual[i],
            expected[i],
            xin.is_some()
        );
    }
    assert!(actual[nout..].iter().all(|&v| v == 0xa5));
    actual.truncate(nout);
    actual
}

fn bytes(r: &mut XorShift, n: usize, case: usize) -> Vec<u8> {
    let constant = r.next_u32() as u8;
    (0..n)
        .map(|i| match case % 9 {
            0 => 0,
            1 => 255,
            2 => constant,
            3 => i as u8,
            4 => {
                if i % 2 == 0 {
                    0
                } else {
                    255
                }
            }
            _ => r.next_u32() as u8,
        })
        .collect()
}

fn output_axis(
    r: &mut XorShift,
    nin: i32,
    din: f32,
    fin: f32,
    nout: i32,
    mode: usize,
) -> (f32, f32) {
    let span = (nin + 3) as f32 / (nout.max(2) - 1) as f32;
    match mode % 10 {
        0 => (din * span, fin - din * 2.0),
        1 => (-din * span, fin + din * (nin + 1) as f32),
        2 => (din, fin),
        3 => (din * 0.0625, fin - din * 0.5),
        4 => (din * 2.25, fin - din),
        5 => (0.0, fin + din * r.uniform_f32(-3.0, nin as f32 + 3.0)),
        6 => (0.0, fin - din * 4.0),
        7 => (0.0, fin + din * (nin + 3) as f32),
        _ => (
            din * r.uniform_f32(-3.0, 3.0),
            fin + din * r.uniform_f32(-3.0, nin as f32 + 3.0),
        ),
    }
}

fn random_case(r: &mut XorShift, i: usize) -> Case {
    let sizes = [1, 2, 3, 5, 7, 16, 31, 32, 97, 128];
    let mut c = Case::default();
    c.nx = sizes[r.range_i32(0, 9) as usize];
    c.ny = sizes[r.range_i32(0, 9) as usize];
    c.ox = sizes[r.range_i32(0, 9) as usize];
    c.oy = sizes[r.range_i32(0, 9) as usize];
    // Long rows/columns supplement prime and power-of-two small grids.
    if i % 257 == 0 {
        c.nx = [257, 1024, 4096][i / 257 % 3];
        c.ny = 3;
        c.ox = 4096;
        c.oy = 5;
    } else if i % 263 == 0 {
        c.nx = 2;
        c.ny = 1024;
        c.ox = 3;
        c.oy = 2048;
    }
    if i % 71 == 0 {
        c.ox = 0;
    }
    if i % 73 == 0 {
        c.oy = 0;
    }
    let scale = [0.000244140625, 0.03125, 1.0, 128.0, 16384.0][i % 5];
    c.dx = r.uniform_f32(0.125, 2.0) * scale;
    c.dy = r.uniform_f32(0.125, 2.0) * scale;
    if r.next_u32() & 1 != 0 {
        c.dx = -c.dx;
    }
    if r.next_u32() & 1 != 0 {
        c.dy = -c.dy;
    }
    c.fx = r.uniform_f32(-100.0, 100.0) * scale;
    c.fy = r.uniform_f32(-100.0, 100.0) * scale;
    (c.odx, c.ofx) = output_axis(r, c.nx, c.dx, c.fx, c.ox, i);
    (c.ody, c.ofy) = output_axis(r, c.ny, c.dy, c.fy, c.oy, i / 10);
    c
}

#[test]
fn intl2b_randomized_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(0x1ab2_0042);
    for i in 0..6000 {
        let c = random_case(&mut r, i);
        let zin = bytes(&mut r, (c.nx * c.ny) as usize, i);
        let actual = check(c, None, &zin);
        if i % 9 == 0 {
            assert!(actual.iter().all(|&b| b == 0));
        }
        if i % 9 == 1 {
            assert!(actual.iter().all(|&b| b == 255));
        }
    }
}

#[test]
fn intlirr2b_randomized_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(0x1abb_0089);
    for i in 0..6000 {
        let mut c = random_case(&mut r, i);
        let step = c.dx.abs();
        let mut xin = vec![c.fx; c.nx as usize];
        for j in 1..xin.len() {
            xin[j] = xin[j - 1] + step * r.uniform_f32(0.125, 4.0);
        }
        let mean_step = if c.nx > 1 {
            (xin[xin.len() - 1] - xin[0]) / (c.nx - 1) as f32
        } else {
            step
        };
        (c.odx, c.ofx) = output_axis(&mut r, c.nx, mean_step, xin[0], c.ox, i);
        let zin = bytes(&mut r, (c.nx * c.ny) as usize, i);
        let actual = check(c, Some(&xin), &zin);
        if i % 9 == 0 {
            assert!(actual.iter().all(|&b| b == 0));
        }
        if i % 9 == 1 {
            assert!(actual.iter().all(|&b| b == 255));
        }
    }
}

fn coefficient_cases(irregular: bool) {
    // Check every byte value at every coefficient and immediately around all
    // rounding boundaries; exercise x and y tables independently.
    let xin = [0.0, 1.0];
    let horizontal: Vec<u8> = (0..=255u8).flat_map(|b| [b, 255 - b]).collect();
    let vertical: Vec<u8> = (0..=255u8).chain((0..=255u8).rev()).collect();
    let mut fractions: Vec<f32> = (0..=99).map(|i| i as f32 / 99.0).collect();
    for i in 0..99 {
        let f = ((i as f64 + 0.5) / 99.0) as f32;
        fractions.extend([
            f32::from_bits(f.to_bits() - 1),
            f,
            f32::from_bits(f.to_bits() + 1),
        ]);
    }
    for frac in fractions {
        let h = Case {
            nx: 2,
            ny: 256,
            ox: 1,
            oy: 256,
            odx: 0.0,
            ofx: frac,
            ody: 1.0,
            ofy: 0.0,
            ..Case::default()
        };
        check(h, irregular.then_some(xin.as_slice()), &horizontal);
        let v = Case {
            nx: 256,
            ny: 2,
            ox: 256,
            oy: 1,
            odx: 1.0,
            ofx: 0.0,
            ody: 0.0,
            ofy: frac,
            ..Case::default()
        };
        let vx: Vec<f32> = (0..256).map(|i| i as f32).collect();
        check(v, irregular.then_some(vx.as_slice()), &vertical);
    }
}

#[test]
fn intl2b_all_coefficients_and_rounding_boundaries_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    coefficient_cases(false);
}

#[test]
fn intlirr2b_all_coefficients_and_rounding_boundaries_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    coefficient_cases(true);
}

#[test]
fn singleton_zero_spacing_and_extreme_coordinates_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(0x1eed);
    for nx in [1, 2, 7] {
        for ny in [1, 2, 7] {
            for dx in [0.0, -0.0, 1.0, -2.0, f32::MIN_POSITIVE, f32::from_bits(1)] {
                for dy in [0.0, -0.0, 1.0, -2.0, f32::MIN_POSITIVE, f32::from_bits(1)] {
                    let zin = bytes(&mut r, (nx * ny) as usize, 8);
                    for (ofx, ofy) in [(-8.0, -8.0), (0.0, 0.0), (8.0, 8.0)] {
                        let c = Case {
                            nx,
                            ny,
                            dx,
                            dy,
                            ox: 13,
                            oy: 17,
                            ofx,
                            ofy,
                            odx: 0.5,
                            ody: -0.5,
                            ..Case::default()
                        };
                        check(c, None, &zin);
                        let xin: Vec<f32> = (0..nx).map(|i| i as f32 - 3.0).collect();
                        check(c, Some(&xin), &zin);
                    }
                }
            }
        }
    }
    for origin in [1.0e20, -1.0e20, 16777216.0, -16777216.0] {
        let c = Case {
            fx: origin,
            fy: origin,
            ofx: origin,
            ofy: origin,
            odx: 0.25,
            ody: -0.25,
            ..Case::default()
        };
        // Small increments round back to the original f32 coordinate.
        check(c, None, &[0, 255, 127, 42]);
        check(c, Some(&[-1.0, 1.0]), &[0, 255, 127, 42]);
    }
}

#[test]
fn large_grids_and_y_row_reuse_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(0x1234_5678);
    for (nx, ny, ox, oy) in [(257, 257, 512, 512), (32, 31, 97, 128), (2, 2, 17, 19)] {
        let zin = bytes(&mut r, (nx * ny) as usize, 8);
        let xin: Vec<f32> = (0..nx).map(|i| (i * i + i) as f32 * 0.25 - 10.0).collect();
        for direction in [-1.0, 1.0] {
            for step in [0.0, 0.0625, 1.0, 2.5] {
                let c = Case {
                    nx,
                    ny,
                    ox,
                    oy,
                    odx: 0.5,
                    ofx: -1.0,
                    ody: direction * step,
                    ofy: if direction > 0.0 {
                        -2.0
                    } else {
                        ny as f32 + 1.0
                    },
                    ..Case::default()
                };
                check(c, None, &zin);
                check(c, Some(&xin), &zin);
            }
        }
    }
}

#[test]
fn irregular_exact_knots_match_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let xin = [
        -100.0, -12.0, -11.999, -1.0, 0.0, 0.0001, 9.0, 10000.0, 1.0e10,
    ];
    let zin: Vec<u8> = (0..xin.len() * 7).map(|i| (i * 37) as u8).collect();
    for x in xin {
        let c = Case {
            nx: xin.len() as i32,
            ny: 7,
            ox: 1,
            oy: 7,
            odx: 0.0,
            ofx: x,
            ody: 1.0,
            ofy: 0.0,
            ..Case::default()
        };
        let actual = check(c, Some(&xin), &zin);
        let ix = xin.iter().position(|&v| v == x).unwrap();
        for iy in 0..7 {
            assert_eq!(actual[iy], zin[iy * xin.len() + ix]);
        }
    }
}
