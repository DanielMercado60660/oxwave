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

//! Differential tests for cwp/lib/sqr.c, adapted from LINPACK FORTRAN.
//! Credits: Dave Hale, Colorado School of Mines, 12/29/89 (sqrdc, sqrsl)
//! and 12/31/89 (sqrst).

use oxwave_cref::XorShift;
use oxwave_qr::{sqrdc as rust_sqrdc, sqrsl as rust_sqrsl, sqrst as rust_sqrst};

extern "C" {
    fn sqrdc(
        x: *mut *mut f32,
        n: i32,
        p: i32,
        qraux: *mut f32,
        jpvt: *mut i32,
        work: *mut f32,
        job: i32,
    );
    fn sqrsl(
        x: *mut *mut f32,
        n: i32,
        k: i32,
        qraux: *mut f32,
        y: *mut f32,
        qy: *mut f32,
        qty: *mut f32,
        b: *mut f32,
        rsd: *mut f32,
        xb: *mut f32,
        job: i32,
        info: *mut i32,
    );
    fn sqrst(
        x: *mut *mut f32,
        n: i32,
        p: i32,
        y: *mut f32,
        tol: f32,
        b: *mut f32,
        rsd: *mut f32,
        k: *mut i32,
        jpvt: *mut i32,
        qraux: *mut f32,
        work: *mut f32,
    );
}

// Tall, wide, square; 1, 2, primes, powers of two, and larger dimensions.
const SHAPES: &[(usize, usize)] = &[
    (1, 1),
    (1, 7),
    (2, 1),
    (2, 2),
    (2, 7),
    (3, 5),
    (5, 3),
    (7, 7),
    (8, 13),
    (13, 8),
    (16, 16),
    (17, 31),
    (31, 17),
    (32, 32),
    (64, 7),
    (7, 64),
    (97, 5),
    (5, 97),
    (127, 8),
    (8, 127),
    (256, 3),
    (3, 256),
];
const LARGE: &[(usize, usize)] = &[(128, 128), (257, 17), (17, 257), (1024, 5), (5, 1024)];

fn pointers(x: &mut [Vec<f32>]) -> Vec<*mut f32> {
    x.iter_mut().map(|c| c.as_mut_ptr()).collect()
}
fn slices(x: &mut [Vec<f32>]) -> Vec<&mut [f32]> {
    x.iter_mut().map(|c| c.as_mut_slice()).collect()
}

// Require bit equality for every output, stronger than the requested 1e-5
// relative-to-maximum bound. This also checks untouched entries and signed zero.
fn same(a: &[f32], b: &[f32], label: &str, case: usize) {
    assert_eq!(a.len(), b.len());
    for i in 0..a.len() {
        assert_eq!(
            a[i].to_bits(),
            b[i].to_bits(),
            "{label} case {case}, index {i}: Rust {:?}, C {:?}",
            a[i],
            b[i]
        );
    }
}
fn same_matrix(a: &[Vec<f32>], b: &[Vec<f32>], case: usize) {
    for (a, b) in a.iter().zip(b) {
        same(a, b, "x", case);
    }
}

fn matrix(rng: &mut XorShift, n: usize, p: usize, mode: usize) -> Vec<Vec<f32>> {
    let scale = [1.0f32, 1e-12, 1e12, 1e-4, 1e4][mode % 5];
    // Two guard entries on every column verify that swaps touch only n rows.
    let mut x: Vec<Vec<f32>> = (0..p).map(|_| rng.vec_f32(n + 2, -scale, scale)).collect();
    match mode % 12 {
        0 => {
            for col in &mut x {
                col[..n].fill(0.0);
            }
        }
        1 => {
            for col in &mut x {
                col[..n].fill(-0.0);
            }
        }
        2 => {
            // Diagonal, with exact zero columns when p > n.
            for j in 0..p {
                x[j][..n].fill(0.0);
                if j < n {
                    x[j][j] = if j % 2 == 0 { scale } else { -scale };
                }
            }
        }
        3 => {
            // Repeated columns: deficient rank and norm recomputation.
            for j in 1..p {
                for i in 0..n {
                    x[j][i] = x[0][i];
                }
            }
        }
        4 => {
            // Zero columns mixed with random columns.
            for j in (0..p).step_by(2) {
                x[j][..n].fill(0.0);
            }
        }
        5 => {
            // Nearly dependent columns.
            for j in 1..p {
                for i in 0..n {
                    x[j][i] = x[0][i] + x[j][i] * 1e-6;
                }
            }
        }
        6 => {
            for col in &mut x {
                for v in &mut col[..n] {
                    *v = -v.abs();
                }
            }
        }
        7 => {
            // Mixed column scales.
            for j in 0..p {
                let s = [1e-6, 1.0, 1e6][j % 3];
                for v in &mut x[j][..n] {
                    *v *= s;
                }
            }
        }
        8 => {
            // Equal-norm columns and pivot ties.
            for j in 0..p {
                x[j][..n].fill(0.0);
                x[j][j % n] = if j % 2 == 0 { scale } else { -scale };
            }
        }
        _ => {}
    }
    x
}

#[test]
fn sqrdc_matches_c() {
    let mut rng = XorShift::new(0x7172_0101);
    for case in 0..4000 + LARGE.len() {
        let (n, p) = if case < 4000 {
            SHAPES[case % SHAPES.len()]
        } else {
            LARGE[case - 4000]
        };
        let mode = if case < 4000 {
            case / SHAPES.len()
        } else {
            11 + 12 * (case - 4000)
        };
        let mut c = matrix(&mut rng, n, p, mode);
        let mut r = c.clone();
        let mut ca = rng.vec_f32(p + 2, -2.0, 2.0);
        let mut ra = ca.clone();
        let mut cw = rng.vec_f32(p + 2, -2.0, 2.0);
        let mut rw = cw.clone();
        let mut cp: Vec<i32> = (0..p + 2)
            .map(|_| match case % 7 {
                0 => 0,
                1 => -1,
                2 => 1,
                _ => rng.range_i32(-5, 5),
            })
            .collect();
        let mut rp = cp.clone();
        let job = [0, 1, -1, 7][case % 4];
        unsafe {
            sqrdc(
                pointers(&mut c).as_mut_ptr(),
                n as i32,
                p as i32,
                ca.as_mut_ptr(),
                cp.as_mut_ptr(),
                cw.as_mut_ptr(),
                job,
            );
        }
        rust_sqrdc(
            &mut slices(&mut r),
            n as i32,
            p as i32,
            &mut ra,
            &mut rp,
            &mut rw,
            job,
        );
        same_matrix(&r, &c, case);
        same(&ra, &ca, "qraux", case);
        same(&rw, &cw, "work", case);
        assert_eq!(rp, cp, "jpvt case {case}");
    }
    // sqrdc itself supports empty dimensions without dereferencing x rows.
    for (n, p) in [(0, 0), (0, 3), (3, 0)] {
        for job in [0, 1] {
            let mut c = vec![vec![]; p as usize];
            let mut r = c.clone();
            let mut ca = vec![7.0; p as usize];
            let mut ra = ca.clone();
            let mut cp = vec![0; p as usize];
            let mut rp = cp.clone();
            let mut cw = ca.clone();
            let mut rw = cw.clone();
            unsafe {
                sqrdc(
                    pointers(&mut c).as_mut_ptr(),
                    n,
                    p,
                    ca.as_mut_ptr(),
                    cp.as_mut_ptr(),
                    cw.as_mut_ptr(),
                    job,
                );
            }
            rust_sqrdc(&mut slices(&mut r), n, p, &mut ra, &mut rp, &mut rw, job);
            same_matrix(&r, &c, 4005);
            same(&ra, &ca, "empty qraux", 4005);
            same(&rw, &cw, "empty work", 4005);
            assert_eq!(rp, cp);
        }
    }
}

#[test]
fn sqrsl_matches_c() {
    let mut rng = XorShift::new(0x7172_0202);
    for case in 0..4096 + LARGE.len() {
        let (n, p) = if case < 4096 {
            SHAPES[case % SHAPES.len()]
        } else {
            LARGE[case - 4096]
        };
        let mut c = matrix(&mut rng, n, p, case / SHAPES.len());
        let mut qraux = vec![0.0; p];
        let mut jpvt = vec![0; p];
        let mut work = vec![0.0; p];
        // Factor with C to test sqrsl independently of our sqrdc.
        unsafe {
            sqrdc(
                pointers(&mut c).as_mut_ptr(),
                n as i32,
                p as i32,
                qraux.as_mut_ptr(),
                jpvt.as_mut_ptr(),
                work.as_mut_ptr(),
                (case % 2) as i32,
            );
        }
        let k = match (case / 32) % 4 {
            0 => n.min(p),
            1 => 1,
            2 => rng.range_i32(1, n.min(p) as i32) as usize,
            _ => 0,
        };
        // Exercise exactly singular R at every diagonal, including index zero.
        if case % 11 == 0 && k > 0 {
            c[case % k][case % k] = 0.0;
        }
        let mut r = c.clone();
        let original = c.clone();
        let mut y = rng.vec_f32(n + 2, -100.0, 100.0);
        if case % 13 == 0 {
            y[..n].fill(0.0);
        }
        let y0 = y.clone();
        let qraux0 = qraux.clone();
        let mask = case % 32;
        let mut job = 0;
        for (bit, digit) in [10000, 1000, 100, 10, 1].iter().enumerate() {
            if mask & (1 << bit) != 0 {
                job += digit;
            }
        }
        if case % 37 == 0 {
            job *= -3;
        } // Nonbinary and negative nonzero digits.
        let requested = [
            job / 10000 != 0,
            job % 10000 != 0,
            (job % 1000) / 100 != 0,
            (job % 100) / 10 != 0,
            job % 10 != 0,
        ];
        let mut co: Vec<Vec<f32>> = (0..5)
            .map(|i| rng.vec_f32(if i == 2 { k.max(1) + 2 } else { n + 2 }, -10.0, 10.0))
            .collect();
        let mut ro = co.clone();
        // Alternate sentinel-backed unrequested arrays with NULL/empty arrays.
        let empty_unrequested = case % 2 == 0;
        let outptr: Vec<*mut f32> = co
            .iter_mut()
            .enumerate()
            .map(|(i, v)| {
                if empty_unrequested && !requested[i] {
                    std::ptr::null_mut()
                } else {
                    v.as_mut_ptr()
                }
            })
            .collect();
        let mut ci = -999;
        let mut ri = -999;
        unsafe {
            sqrsl(
                pointers(&mut c).as_mut_ptr(),
                n as i32,
                k as i32,
                qraux.as_mut_ptr(),
                y.as_mut_ptr(),
                outptr[0],
                outptr[1],
                outptr[2],
                outptr[3],
                outptr[4],
                job,
                &mut ci,
            );
        }
        let mut outs: Vec<&mut [f32]> = ro
            .iter_mut()
            .enumerate()
            .map(|(i, v)| {
                if empty_unrequested && !requested[i] {
                    &mut v[..0]
                } else {
                    v.as_mut_slice()
                }
            })
            .collect();
        let (qy, tail) = outs.split_at_mut(1);
        let (qty, tail) = tail.split_at_mut(1);
        let (b, tail) = tail.split_at_mut(1);
        let (rsd, xb) = tail.split_at_mut(1);
        rust_sqrsl(
            &mut slices(&mut r),
            n as i32,
            k as i32,
            &qraux0,
            &y0,
            qy[0],
            qty[0],
            b[0],
            rsd[0],
            xb[0],
            job,
            &mut ri,
        );
        assert_eq!(ri, ci, "info case {case} n={n} p={p} k={k} job={job}");
        for i in 0..5 {
            same(&ro[i], &co[i], ["qy", "qty", "b", "rsd", "xb"][i], case);
        }
        same_matrix(&r, &c, case);
        same_matrix(&r, &original, case); // Temporary diagonal replacement is restored.
        same(&y, &y0, "y unchanged", case);
        same(&qraux, &qraux0, "qraux unchanged", case);
    }
}

#[test]
fn sqrst_matches_c() {
    let mut rng = XorShift::new(0x7172_0303);
    let tolerances = [0.0, 1e-7, 1e-5, 0.01, 0.2, 0.9, 1.0, 2.0, -0.1];
    for case in 0..4000 + LARGE.len() {
        let (n, p) = if case < 4000 {
            SHAPES[case % SHAPES.len()]
        } else {
            LARGE[case - 4000]
        };
        let mode = if case < 4000 {
            case / SHAPES.len()
        } else {
            11 + 12 * (case - 4000)
        };
        let mut c = matrix(&mut rng, n, p, mode);
        let mut r = c.clone();
        let mut cy = rng.vec_f32(n + 2, -100.0, 100.0);
        if case % 17 == 0 {
            cy[..n].fill(0.0);
        }
        let ry = cy.clone();
        let mut cb = rng.vec_f32(p + 2, -7.0, 7.0);
        let mut rb = cb.clone();
        let mut cr = rng.vec_f32(n + 2, -7.0, 7.0);
        let mut rr = cr.clone();
        let mut cp = vec![-99; p + 2];
        let mut rp = cp.clone();
        let mut ca = rng.vec_f32(p + 2, -7.0, 7.0);
        let mut ra = ca.clone();
        let mut cw = rng.vec_f32(p + 2, -7.0, 7.0);
        let mut rw = cw.clone();
        let mut ck = -99;
        let mut rk = -99;
        let tol = tolerances[case % tolerances.len()];
        unsafe {
            sqrst(
                pointers(&mut c).as_mut_ptr(),
                n as i32,
                p as i32,
                cy.as_mut_ptr(),
                tol,
                cb.as_mut_ptr(),
                cr.as_mut_ptr(),
                &mut ck,
                cp.as_mut_ptr(),
                ca.as_mut_ptr(),
                cw.as_mut_ptr(),
            );
        }
        rust_sqrst(
            &mut slices(&mut r),
            n as i32,
            p as i32,
            &ry,
            tol,
            &mut rb,
            &mut rr,
            &mut rk,
            &mut rp,
            &mut ra,
            &mut rw,
        );
        assert_eq!(rk, ck, "rank case {case} tol={tol}");
        assert_eq!(rp, cp, "jpvt case {case}");
        same_matrix(&r, &c, case);
        same(&rb, &cb, "b", case);
        same(&rr, &cr, "rsd", case);
        same(&ra, &ca, "qraux", case);
        same(&rw, &cw, "work", case);
        same(&ry, &cy, "y", case);
    }
}
