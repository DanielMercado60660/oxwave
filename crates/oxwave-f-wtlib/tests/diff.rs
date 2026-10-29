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

//! Differential tests for par/lib/wtlib.c (Revision 1.4).
//! CWP: Zhaobo Meng, Colorado School of Mines, Sept 1996.
//! Modified: Carlos E. Theodoro, Colorado School of Mines, Jun 97.

use oxwave_cref::XorShift;
use oxwave_f_wtlib::*;
use std::process::Command;
use std::ptr;
use std::sync::Mutex;

// wto1dset writes the original C library's static reverse-filter arrays.
static LOCK: Mutex<()> = Mutex::new(());

#[repr(C)]
#[derive(Default)]
struct CF {
    order: i32,
    ioff: i32,
    joff: i32,
    cc: *mut f32,
    cr: *mut f32,
    phi: *mut f32,
    psi: *mut f32,
    weight: *mut f32,
}

#[repr(C)]
struct CS {
    maxlevel: i32,
    ndim: i32,
    order: i32,
    npoints: *mut i32,
    levels: *mut i32,
}

mod c {
    use super::{CF, CS};
    extern "C" {
        pub fn fhierfromcd(
            f: *mut f32,
            cd: *mut f32,
            wfilt: *mut CF,
            wtsizes: *mut CS,
            just_dc: i32,
        );
        pub fn wt1(
            cd: *mut f32,
            tocord: i32,
            npoints: i32,
            wfilt: *mut CF,
            wtsizes: *mut CS,
            idim: i32,
        );
        pub fn wt_cascade(wfilt: *mut CF, wtsizes: *mut CS);
        pub fn wtn(cd: *mut f32, tocord: i32, wfilt: *mut CF, wtsizes: *mut CS, dconly: i32);
        pub fn wto1d(cd: *mut f32, npoints: i32, tocord: i32, wfilt: *mut CF);
        pub fn wto1dset(wfilt: *mut CF, wtsizes: *mut CS);
    }
}

extern "C" {
    fn free(p: *mut std::ffi::c_void);
}

struct Pair {
    f: WtFilter,
    s: WtSizes,
    cf: CF,
    cs: CS,
}

impl Pair {
    fn new(order: i32, maxlevel: i32, points: &[i32], levels: &[i32]) -> Self {
        assert_eq!(points.len(), levels.len());
        let mut s = WtSizes {
            MaxLevel: maxlevel,
            NDim: points.len() as i32,
            order,
            NPointsn: [vec![0], points.to_vec()].concat(),
            Mraleveln: [vec![0], levels.to_vec()].concat(),
        };
        let mut cs = CS {
            maxlevel,
            ndim: s.NDim,
            order,
            npoints: s.NPointsn.as_mut_ptr(),
            levels: s.Mraleveln.as_mut_ptr(),
        };
        let mut cf = CF::default();
        unsafe { c::wto1dset(&mut cf, &mut cs) };
        let mut f = WtFilter::default();
        wto1dset(&mut f, &s);
        Self { f, s, cf, cs }
    }

    fn free_wavelets(&mut self) {
        unsafe {
            free(self.cf.phi.cast());
            free(self.cf.psi.cast());
        }
        self.cf.phi = ptr::null_mut();
        self.cf.psi = ptr::null_mut();
    }

    fn compare_wavelets(&self) {
        let len = (self.s.order << self.s.MaxLevel) as usize;
        assert_eq!(self.f.phi.len(), len);
        assert_eq!(self.f.psi.len(), len);
        unsafe {
            compare(&self.f.phi, std::slice::from_raw_parts(self.cf.phi, len));
            compare(&self.f.psi, std::slice::from_raw_parts(self.cf.psi, len));
        }
    }
}

impl Drop for Pair {
    fn drop(&mut self) {
        self.free_wavelets();
    }
}

#[track_caller]
fn compare(rust: &[f32], c: &[f32]) {
    assert_eq!(rust.len(), c.len());
    let scale = c.iter().fold(0.0f64, |m, &v| m.max((v as f64).abs()));
    for (i, (&x, &y)) in rust.iter().zip(c).enumerate() {
        assert!(
            x.is_finite()
                && y.is_finite()
                && ((x as f64) - (y as f64)).abs() <= 1e-5 * scale + 1e-30,
            "sample {i}: Rust {x} != C {y}, max C magnitude {scale}"
        );
    }
}

fn exact(rust: &[f32], c: &[f32]) {
    assert_eq!(rust.len(), c.len());
    for (i, (&x, &y)) in rust.iter().zip(c).enumerate() {
        assert_eq!(x.to_bits(), y.to_bits(), "sample {i}: {x} != {y}");
    }
}

fn choose(r: &mut XorShift, n: usize) -> usize {
    r.next_u32() as usize % n
}

fn order(r: &mut XorShift) -> i32 {
    4 + 2 * choose(r, 9) as i32
}

fn direction(r: &mut XorShift) -> ToCorD {
    if choose(r, 2) == 0 {
        ToCorD::ToC
    } else {
        ToCorD::ToD
    }
}

fn samples(r: &mut XorShift, n: usize, case: usize) -> Vec<f32> {
    let scale = [1e-12, 1e-6, 1.0, 1e4, 1e12][case % 5];
    (0..n)
        .map(|i| match case % 23 {
            0 => 0.0,
            1 => -0.0,
            2 => -r.uniform_f32(0.0, 1.0) * scale,
            3 => scale,
            4 => {
                if i % 2 == 0 {
                    scale
                } else {
                    -scale
                }
            }
            5 => {
                if i == n / 2 {
                    -scale
                } else {
                    0.0
                }
            }
            6 => (i as f32 / n as f32 - 0.5) * scale,
            _ => r.uniform_f32(-1.0, 1.0) * scale,
        })
        .collect()
}

// Execute fatal C/Rust branches in a child process, so neither can terminate the
// differential-test parent. This is called only by the two corresponding tests.
fn error_child() {
    let Ok(case) = std::env::var("OXWAVE_WTLIB_ERROR") else {
        return;
    };
    let fields: Vec<&str> = case.split(':').collect();
    let value: i32 = fields[2].parse().unwrap();
    let mut p = Pair::new(4, 0, &[16], &[1]);
    if fields[1] == "setup" {
        p.s.order = value;
        p.cs.order = value;
        if fields[0] == "c" {
            unsafe { c::wto1dset(&mut p.cf, &mut p.cs) };
        } else {
            wto1dset(&mut p.f, &p.s);
        }
    } else {
        // A level of two gives coarsest length 0, 1, 2 or 3 in these cases.
        p.s.Mraleveln[1] = 2;
        let mut cd = vec![0.0; value as usize + 1];
        if fields[0] == "c" {
            unsafe { c::wt1(cd.as_mut_ptr(), 1, value, &mut p.cf, &mut p.cs, 1) };
        } else {
            wt1(&mut cd, ToCorD::ToD, value, &p.f, &p.s, 1);
        }
    }
    panic!("fatal branch unexpectedly returned");
}

fn compare_error(test: &str, function: &str, value: i32, expected: &[u8]) {
    for implementation in ["c", "rust"] {
        let output = Command::new(std::env::current_exe().unwrap())
            .args(["--exact", test, "--nocapture", "--test-threads=1"])
            .env(
                "OXWAVE_WTLIB_ERROR",
                format!("{implementation}:{function}:{value}"),
            )
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{implementation}: {output:?}"
        );
        assert_eq!(output.stderr, expected, "{implementation}");
    }
}

#[test]
fn wto1dset_matches() {
    error_child();
    let _guard = LOCK.lock().unwrap();
    let mut r = XorShift::new(1);
    for case in 0..3000 {
        // All orders at every level 0..=6, then thousands of random settings.
        let (order, maxlevel) = if case < 63 {
            (4 + 2 * (case % 9) as i32, (case / 9) as i32)
        } else {
            (order(&mut r), choose(&mut r, 5) as i32)
        };
        let p = Pair::new(order, maxlevel, &[1, 2], &[0, 0]);
        assert_eq!(p.f.order, p.cf.order);
        assert_eq!(p.f.ioff, p.cf.ioff);
        assert_eq!(p.f.joff, p.cf.joff);
        unsafe {
            exact(
                &p.f.cc,
                std::slice::from_raw_parts(p.cf.cc, order as usize + 1),
            );
            exact(
                &p.f.cr,
                std::slice::from_raw_parts(p.cf.cr, order as usize + 1),
            );
        }
        p.compare_wavelets();
    }
    for order in [-4, 0, 1, 2, 3, 5, 21, 22] {
        compare_error(
            "wto1dset_matches",
            "setup",
            order,
            b"unimplemented order in wto1dset",
        );
    }
}

#[test]
fn wt_cascade_matches() {
    let _guard = LOCK.lock().unwrap();
    let mut r = XorShift::new(2);
    for case in 0..3000 {
        let mut p = Pair::new(order(&mut r), 0, &[8, 8], &[1, 1]);
        // The public cascade also accepts caller-provided coefficient arrays:
        // cover one/two taps, primes, zero filters, and negative coefficients.
        let taps = [1, 2, 3, 4, 5, 7, 8, 12, 16, 20, 31][choose(&mut r, 11)];
        p.s.order = taps;
        p.cs.order = taps;
        p.s.MaxLevel = choose(&mut r, 5) as i32;
        p.cs.maxlevel = p.s.MaxLevel;
        let mut cc = (0..=taps)
            .map(|_| {
                if case % 17 == 0 {
                    0.0
                } else {
                    r.uniform_f32(-1.0, 1.0)
                }
            })
            .collect::<Vec<_>>();
        p.f.cc = cc.clone();
        p.cf.cc = cc.as_mut_ptr();
        p.free_wavelets();
        unsafe { c::wt_cascade(&mut p.cf, &mut p.cs) };
        wt_cascade(&mut p.f, &p.s);
        p.compare_wavelets();
        exact(&p.f.cc, &cc);
    }
}

#[test]
fn wto1d_matches() {
    let _guard = LOCK.lock().unwrap();
    let mut r = XorShift::new(3);
    let sizes = [0, 1, 2, 3, 5, 7, 8, 16, 31, 64, 127, 256, 1024, 4096, 16384];
    for case in 0..4000 {
        let mut p = Pair::new(order(&mut r), 0, &[8], &[1]);
        let n = if case < sizes.len() {
            sizes[case]
        } else {
            sizes[choose(&mut r, sizes.len())]
        };
        // Independent offsets and arbitrary filters exercise the operator beyond
        // the built-in Daubechies initialization. Keep both backing arrays alive.
        let mut cc = p.f.cc.clone();
        let mut cr = p.f.cr.clone();
        if case % 3 == 0 {
            for x in cc.iter_mut().chain(cr.iter_mut()).skip(1) {
                *x = r.uniform_f32(-1.0, 1.0);
            }
            p.f.cc.clone_from(&cc);
            p.f.cr.clone_from(&cr);
            // Exercise the C build's wrapping index arithmetic, including
            // overflow before applying its npoints - 1 bit mask.
            let offsets = [
                i32::MIN,
                i32::MIN + 1,
                -1024,
                -1,
                0,
                1,
                1024,
                i32::MAX - 1,
                i32::MAX,
            ];
            p.f.ioff = offsets[choose(&mut r, offsets.len())];
            p.f.joff = offsets[choose(&mut r, offsets.len())];
            p.cf.cc = cc.as_mut_ptr();
            p.cf.cr = cr.as_mut_ptr();
            p.cf.ioff = p.f.ioff;
            p.cf.joff = p.f.joff;
        }
        let d = direction(&mut r);
        let mut a = samples(&mut r, n + 3, case);
        let mut b = a.clone();
        let before = a.clone();
        wto1d(&mut a, n as i32, d, &p.f);
        unsafe { c::wto1d(b.as_mut_ptr(), n as i32, d as i32, &mut p.cf) };
        compare(&a[1..=n], &b[1..=n]);
        exact(&a[..1], &before[..1]);
        exact(&b[..1], &before[..1]);
        exact(&a[n + 1..], &before[n + 1..]);
        exact(&b[n + 1..], &before[n + 1..]);
    }
}

#[test]
fn wt1_matches() {
    error_child();
    let _guard = LOCK.lock().unwrap();
    let mut r = XorShift::new(4);
    for case in 0..4000 {
        let n = 1usize << (2 + choose(&mut r, 13));
        let level = choose(&mut r, n.ilog2() as usize - 1) as i32;
        let idim = 1 + choose(&mut r, 3);
        let mut points = vec![8; 3];
        let mut levels = vec![0; 3];
        points[idim - 1] = n as i32;
        levels[idim - 1] = level;
        let mut p = Pair::new(order(&mut r), 0, &points, &levels);
        let mut a = samples(&mut r, n + 3, case);
        let mut b = a.clone();
        let before = a.clone();
        let d = direction(&mut r);
        wt1(&mut a, d, n as i32, &p.f, &p.s, idim as i32);
        unsafe {
            c::wt1(
                b.as_mut_ptr(),
                d as i32,
                n as i32,
                &mut p.cf,
                &mut p.cs,
                idim as i32,
            )
        };
        compare(&a[1..=n], &b[1..=n]);
        exact(&a[..1], &before[..1]);
        exact(&b[..1], &before[..1]);
        exact(&a[n + 1..], &before[n + 1..]);
        exact(&b[n + 1..], &before[n + 1..]);
        // Also compare an inverse immediately following every forward result.
        if d == ToCorD::ToD {
            wt1(&mut a, ToCorD::ToC, n as i32, &p.f, &p.s, idim as i32);
            unsafe {
                c::wt1(
                    b.as_mut_ptr(),
                    0,
                    n as i32,
                    &mut p.cf,
                    &mut p.cs,
                    idim as i32,
                )
            };
            compare(&a[1..=n], &b[1..=n]);
            // Compare reconstruction to C as well: its tabulated filters are
            // preserved verbatim, including their reconstruction error.
        }
    }
    for n in [1, 2, 3, 4, 7, 8, 12, 15] {
        compare_error(
            "wt1_matches",
            "wt1",
            n,
            format!(
                "lenlastc={}\nStLevel too large or npoints too small.\n",
                n >> 2
            )
            .as_bytes(),
        );
    }
}

#[test]
fn wtn_matches() {
    let _guard = LOCK.lock().unwrap();
    let mut r = XorShift::new(5);
    for case in 0..5000 {
        let ndim = 1 + choose(&mut r, 4);
        let mut points = Vec::new();
        let mut levels = Vec::new();
        let mut total = 1usize;
        for _ in 0..ndim {
            let mut n = 1 << choose(&mut r, 5);
            if total * n > 4096 {
                n = 1;
            }
            points.push(n as i32);
            levels.push(choose(&mut r, n.ilog2() as usize + 1) as i32);
            total *= n;
        }
        // Odd and prime axes exercise both unchanged level-zero axes and the
        // original bit-mask operator at positive levels, even though those
        // lengths need not give a mathematically invertible wavelet transform.
        if case % 11 == 0 {
            points[0] = [1, 2, 3, 7, 31, 127][choose(&mut r, 6)];
            levels[0] = if case % 22 == 0 {
                choose(&mut r, (points[0] as u32).ilog2() as usize + 1) as i32
            } else {
                0
            };
            total = points.iter().map(|&n| n as usize).product();
        }
        let mut p = Pair::new(order(&mut r), 0, &points, &levels);
        let dc = if ndim >= 2 {
            [-1, 0, 1, 2][choose(&mut r, 4)]
        } else {
            0
        };
        let d = direction(&mut r);
        let mut a = samples(&mut r, total + 3, case);
        let mut b = a.clone();
        let before = a.clone();
        wtn(&mut a, d, &p.f, &p.s, dc);
        unsafe { c::wtn(b.as_mut_ptr(), d as i32, &mut p.cf, &mut p.cs, dc) };
        compare(&a[1..=total], &b[1..=total]);
        exact(&a[..1], &before[..1]);
        exact(&b[..1], &before[..1]);
        exact(&a[total + 1..], &before[total + 1..]);
        exact(&b[total + 1..], &before[total + 1..]);
        if d == ToCorD::ToD {
            wtn(&mut a, ToCorD::ToC, &p.f, &p.s, dc);
            unsafe { c::wtn(b.as_mut_ptr(), 0, &mut p.cf, &mut p.cs, dc) };
            compare(&a[1..=total], &b[1..=total]);
        }
    }
}

#[test]
fn fhierfromcd_matches() {
    let _guard = LOCK.lock().unwrap();
    let mut r = XorShift::new(6);
    for case in 0..3000 {
        let (n, m) = if case % 500 == 0 {
            (128usize, 64usize)
        } else {
            (1 << choose(&mut r, 6), 1 << (1 + choose(&mut r, 5)))
        };
        let l1 = choose(&mut r, n.ilog2() as usize + 1) as i32;
        let mut l2 = choose(&mut r, m.ilog2() as usize + 1) as i32;
        // Equal zero levels make C shift by -1 (undefined behavior). Generate
        // a valid nonzero second level directly, rather than skipping a case.
        if l1 == 0 && l2 == 0 {
            l2 = 1;
        }
        let mut p = Pair::new(order(&mut r), 0, &[n as i32, m as i32], &[l1, l2]);
        let mut cd = samples(&mut r, n * m + 3, case);
        let cd_before = cd.clone();
        let mut a = samples(&mut r, cd.len(), case + 1);
        let mut b = a.clone();
        let before = a.clone();
        let dc = [-1, 0, 1, 2][choose(&mut r, 4)];
        fhierfromcd(&mut a, &cd, &p.f, &p.s, dc);
        unsafe { c::fhierfromcd(b.as_mut_ptr(), cd.as_mut_ptr(), &mut p.cf, &mut p.cs, dc) };
        compare(&a[..=n * m], &b[..=n * m]);
        exact(&cd, &cd_before);
        exact(&a[n * m + 1..], &before[n * m + 1..]);
        exact(&b[n * m + 1..], &before[n * m + 1..]);
    }
}
