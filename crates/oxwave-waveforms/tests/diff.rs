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

//! Differential tests for cwp/lib/waveforms.c, credited to Tong Fei,
//! Ken Larner, Nils Maercklin, and Werner M. Heigl, with Gaussian derivative
//! code credited to Dominique Rousset.
use oxwave_cref::XorShift;
use oxwave_waveforms as rust;
extern "C" {
    fn akb_wavelet(n: i32, dt: f32, f: f32, w: *mut f32);
    fn ricker1_wavelet(n: i32, dt: f32, f: f32, w: *mut f32);
    fn gaussian_wavelet(n: i32, dt: f32, f: f32, w: *mut f32);
    fn gaussderiv_wavelet(n: i32, dt: f32, f: f32, w: *mut f32);
    fn ricker2_wavelet(n: i32, dt: f32, p: f32, a: f32, d: f32, w: *mut f32);
    fn berlage_wavelet(n: i32, dt: f32, f: f32, a: f32, tn: f32, decay: f32, ipa: f32, w: *mut f32);
    fn spike_wavelet(n: i32, index: i32, w: *mut f32);
    fn unit_wavelet(n: i32, w: *mut f32);
    fn zero_wavelet(n: i32, w: *mut f32);
    fn deriv_n_gauss(
        dt: f64,
        nt: i32,
        t0: f64,
        f: f32,
        n: i32,
        w: *mut f64,
        sign: i32,
        verbose: i32,
    );
}
const SIZES: [usize; 9] = [1, 2, 3, 7, 16, 31, 64, 127, 1024];
fn check(a: &[f32], b: &[f32]) {
    let scale = b
        .iter()
        .filter(|x| x.is_finite())
        .fold(0f64, |m, x| m.max(x.abs() as f64));
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        assert!(
            (x.is_nan() && y.is_nan())
                || x == y
                || (x.is_finite()
                    && y.is_finite()
                    && ((x as f64 - y as f64).abs() <= 1e-5 * scale + 1e-30)),
            "index {i}: rust={x:?} C={y:?}, scale={scale}"
        );
    }
}
fn check_simple(
    rust_fn: fn(f32, f32, &mut [f32]),
    c_fn: unsafe extern "C" fn(i32, f32, f32, *mut f32),
) {
    let mut r = XorShift::new(14);
    for k in 0..3000 {
        let n = SIZES[k % SIZES.len()];
        let dt = if k % 19 == 0 {
            0.0
        } else {
            r.uniform_f32(-0.01, 0.01)
        };
        let f = if k % 23 == 0 {
            0.0
        } else {
            r.uniform_f32(-100.0, 100.0)
        };
        let mut a = vec![0.0; n];
        let mut b = a.clone();
        rust_fn(dt, f, &mut a);
        unsafe {
            c_fn(n as i32, dt, f, b.as_mut_ptr());
        }
        check(&a, &b);
    }
}
#[test]
fn akb_matches() {
    check_simple(rust::akb_wavelet, akb_wavelet);
}
#[test]
fn ricker1_matches() {
    check_simple(rust::ricker1_wavelet, ricker1_wavelet);
}
#[test]
fn gaussian_matches() {
    check_simple(rust::gaussian_wavelet, gaussian_wavelet);
}
#[test]
fn gaussderiv_matches() {
    check_simple(rust::gaussderiv_wavelet, gaussderiv_wavelet);
}
#[test]
fn ricker2_matches() {
    let mut r = XorShift::new(22);
    for k in 0..3000 {
        let n = SIZES[k % SIZES.len()];
        let dt = r.uniform_f32(-0.02, 0.02);
        let p = if k % 23 == 0 {
            0.0
        } else {
            r.uniform_f32(-1.0, 1.0)
        };
        let a = r.uniform_f32(-5.0, 5.0);
        let d = r.uniform_f32(-3.0, 3.0);
        let mut x = vec![0.0; 2 * n - 1];
        let mut y = x.clone();
        rust::ricker2_wavelet(n, dt, p, a, d, &mut x);
        unsafe { ricker2_wavelet(n as i32, dt, p, a, d, y.as_mut_ptr()) };
        check(&x, &y);
    }
}
#[test]
fn berlage_matches() {
    let mut r = XorShift::new(32);
    for k in 0..3000 {
        let n = SIZES[k % SIZES.len()];
        let dt = if k % 17 == 0 {
            0.0
        } else {
            r.uniform_f32(-0.01, 0.01)
        };
        let f = r.uniform_f32(-60.0, 60.0);
        let a = r.uniform_f32(-3.0, 3.0);
        let tn = if k % 2 == 0 {
            r.range_i32(0, 4) as f32
        } else {
            r.uniform_f32(-2.0, 3.0)
        };
        let decay = r.uniform_f32(-3.0, 3.0);
        let ipa = r.uniform_f32(-3.0, 3.0);
        let mut x = vec![0.0; n];
        let mut y = x.clone();
        rust::berlage_wavelet(dt, f, a, tn, decay, ipa, &mut x);
        unsafe { berlage_wavelet(n as i32, dt, f, a, tn, decay, ipa, y.as_mut_ptr()) };
        check(&x, &y);
    }
}
#[test]
fn spike_matches() {
    let mut r = XorShift::new(33);
    for k in 0..3000 {
        let n = SIZES[k % SIZES.len()];
        let i = r.range_i32(0, n as i32 - 1) as usize;
        let mut a = vec![-5.0; n];
        let mut b = a.clone();
        rust::spike_wavelet(i, &mut a);
        unsafe { spike_wavelet(n as i32, i as i32, b.as_mut_ptr()) };
        assert_eq!(a, b);
    }
}
fn check_constant(rust_fn: fn(&mut [f32]), c_fn: unsafe extern "C" fn(i32, *mut f32)) {
    let mut r = XorShift::new(99);
    for k in 0..3000 {
        let n = SIZES[k % SIZES.len()];
        let mut a: Vec<f32> = (0..n).map(|_| r.any_f32()).collect();
        let mut b = a.clone();
        rust_fn(&mut a);
        unsafe { c_fn(n as i32, b.as_mut_ptr()) };
        assert_eq!(a, b);
    }
}
#[test]
fn unit_matches() {
    check_constant(rust::unit_wavelet, unit_wavelet);
}
#[test]
fn zero_matches() {
    check_constant(rust::zero_wavelet, zero_wavelet);
}
#[test]
fn deriv_matches() {
    let mut r = XorShift::new(47);
    for k in 0..3000 {
        let len = SIZES[k % SIZES.len()];
        let dt = if k % 23 == 0 {
            0.0
        } else {
            r.uniform(-0.01, 0.01)
        };
        let t0 = r.uniform(-0.1, 0.1);
        let f = if k % 29 == 0 {
            0.0
        } else {
            r.uniform_f32(-60.0, 60.0)
        };
        let n = r.range_i32(-1, 8);
        let sign = r.range_i32(-1, 1);
        let mut a = vec![0.0; len];
        let mut b = a.clone();
        rust::deriv_n_gauss(dt, t0, f, n, &mut a, sign, 0);
        unsafe { deriv_n_gauss(dt, len as i32, t0, f, n, b.as_mut_ptr(), sign, 0) };
        let scale = b
            .iter()
            .filter(|x| x.is_finite())
            .fold(0f64, |m, x| m.max(x.abs()));
        for (i, (&x, &y)) in a.iter().zip(&b).enumerate() {
            assert!(
                (x.is_nan() && y.is_nan())
                    || x == y
                    || (x.is_finite() && y.is_finite() && (x - y).abs() <= 1e-5 * scale + 1e-30),
                "index {i} rust={x:?} C={y:?} n={n} dt={dt} t0={t0} f={f}"
            );
        }
    }
}
