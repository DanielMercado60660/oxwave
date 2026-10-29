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

//! Differential tests for cwp/lib/dge.c.
//! Author: Dave Hale, Colorado School of Mines, 10/01/89.
use oxwave_cref::XorShift;
use oxwave_f_dge as rust;

extern "C" {
    fn dgefa(a: *mut *mut f64, n: i32, ipvt: *mut i32, info: *mut i32);
    fn dgeco(a: *mut *mut f64, n: i32, ipvt: *mut i32, rcond: *mut f64, z: *mut f64);
    fn dgesl(a: *mut *mut f64, n: i32, ipvt: *mut i32, b: *mut f64, job: i32);
}

fn cases(mut f: impl FnMut(Vec<Vec<f64>>, usize)) {
    let mut rng = XorShift::new(7319);
    let sizes = [1, 2, 3, 4, 5, 7, 8, 11, 16, 31, 32, 64, 97];
    let scales = [0.0, 1e-100, 1e-10, 1.0, 1e10, 1e100];
    for case in 0..2600 {
        let n = sizes[case % sizes.len()];
        let scale = scales[(case / sizes.len()) % scales.len()];
        let mut a: Vec<Vec<f64>> = (0..n)
            .map(|_| (0..n).map(|_| rng.uniform(-2.0, 2.0) * scale).collect())
            .collect();
        match (case / (sizes.len() * scales.len())) % 6 {
            0 => {
                for i in 0..n {
                    a[i][i] += n as f64 * scale;
                }
            }
            1 => {
                for col in &mut a {
                    col[0] = 0.0;
                }
            }
            2 => {
                if n > 1 {
                    a[1] = a[0].clone();
                }
            }
            3 => {
                for i in 0..n {
                    for j in 0..n {
                        if i != j {
                            a[i][j] = 0.0;
                        }
                    }
                }
            }
            _ => {}
        }
        f(a, case);
    }
}
fn pointers(a: &mut [Vec<f64>]) -> Vec<*mut f64> {
    a.iter_mut().map(|c| c.as_mut_ptr()).collect()
}
fn cols(a: &mut [Vec<f64>]) -> Vec<&mut [f64]> {
    a.iter_mut().map(Vec::as_mut_slice).collect()
}
fn compare(x: &[f64], y: &[f64], label: &str) {
    let magnitude = x
        .iter()
        .chain(y)
        .filter(|v| v.is_finite())
        .fold(0.0f64, |a, b| a.max(b.abs()));
    for (i, (&a, &b)) in x.iter().zip(y).enumerate() {
        assert!(
            (a.is_nan() && b.is_nan())
                || a == b
                || (a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-5 * magnitude),
            "{label}[{i}] rust={a:?} C={b:?}, scale={magnitude}"
        );
    }
}
fn matrices(x: &[Vec<f64>], y: &[Vec<f64>]) {
    for (a, b) in x.iter().zip(y) {
        compare(a, b, "matrix");
    }
}

#[test]
fn dgefa_matches_c() {
    cases(|mut a, _| {
        let n = a.len();
        let mut c = a.clone();
        let mut p = vec![-7; n];
        let mut q = p.clone();
        let mut info = 77;
        let mut ci = info;
        unsafe {
            dgefa(
                pointers(&mut c).as_mut_ptr(),
                n as i32,
                q.as_mut_ptr(),
                &mut ci,
            );
        }
        rust::dgefa(&mut cols(&mut a), n as i32, &mut p, &mut info);
        assert_eq!(p, q);
        assert_eq!(info, ci);
        matrices(&a, &c);
    });
}

#[test]
fn dgeco_matches_c() {
    cases(|mut a, _| {
        let n = a.len();
        let mut c = a.clone();
        let mut p = vec![-7; n];
        let mut q = p.clone();
        let mut z = vec![9.0; n];
        let mut cz = z.clone();
        let mut r = 9.0;
        let mut cr = r;
        unsafe {
            dgeco(
                pointers(&mut c).as_mut_ptr(),
                n as i32,
                q.as_mut_ptr(),
                &mut cr,
                cz.as_mut_ptr(),
            );
        }
        rust::dgeco(&mut cols(&mut a), n as i32, &mut p, &mut r, &mut z);
        assert_eq!(p, q);
        matrices(&a, &c);
        compare(&z, &cz, "z");
        compare(&[r], &[cr], "rcond");
    });
}

#[test]
fn dgesl_matches_c() {
    let mut rng = XorShift::new(894);
    cases(|mut a, case| {
        let n = a.len();
        let mut p = vec![0; n];
        let mut info = 0;
        rust::dgefa(&mut cols(&mut a), n as i32, &mut p, &mut info);
        let mut c = a.clone();
        for job in [0, 1, -1] {
            let mut b: Vec<f64> = (0..n)
                .map(|_| {
                    if case % 7 == 0 {
                        0.0
                    } else {
                        rng.uniform(-3.0, 3.0)
                    }
                })
                .collect();
            let mut cb = b.clone();
            let mut cp = p.clone();
            unsafe {
                dgesl(
                    pointers(&mut c).as_mut_ptr(),
                    n as i32,
                    cp.as_mut_ptr(),
                    cb.as_mut_ptr(),
                    job,
                );
            }
            let ac: Vec<&[f64]> = a.iter().map(Vec::as_slice).collect();
            rust::dgesl(&ac, n as i32, &p, &mut b, job);
            compare(&b, &cb, "solution");
            matrices(&a, &c);
            assert_eq!(p, cp);
        }
    });
}
