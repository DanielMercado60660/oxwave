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

//! Differential tests against original `cwp/lib/tridiagonal.c`.
//! Credits: tridif/tridid: Dave Hale, Colorado School of Mines, 10/03/89;
//! tripd/tripp: Zhenyue Liu, Colorado School of Mines, 1992-1993.

use oxwave_cref::XorShift;
use oxwave_tridiag as rust;

extern "C" {
    fn tridid(n: i32, a: *mut f64, b: *mut f64, c: *mut f64, r: *mut f64, u: *mut f64);
    fn tridif(n: i32, a: *mut f32, b: *mut f32, c: *mut f32, r: *mut f32, u: *mut f32);
    fn tripd(d: *mut f32, e: *mut f32, b: *mut f32, n: i32);
    fn tripp(n: i32, d: *mut f32, e: *mut f32, c: *mut f32, b: *mut f32);
}

const SIZES: [usize; 16] = [
    1, 2, 3, 4, 5, 7, 16, 31, 32, 64, 97, 128, 257, 1021, 1024, 4096,
];
const CASES_PER_SIZE: usize = 256;

// Bit equality is stronger than the requested 1e-5 relative tolerance.
// NaN payloads/signs are not part of the comparison, but their positions are.
fn equal32(rust: &[f32], c: &[f32]) {
    assert_eq!(rust.len(), c.len());
    for (i, (&r, &c)) in rust.iter().zip(c).enumerate() {
        assert!(
            r.to_bits() == c.to_bits() || (r.is_nan() && c.is_nan()),
            "index {i}: Rust {r:?} ({:08x}), C {c:?} ({:08x})",
            r.to_bits(),
            c.to_bits()
        );
    }
}

fn equal64(rust: &[f64], c: &[f64]) {
    assert_eq!(rust.len(), c.len());
    for (i, (&r, &c)) in rust.iter().zip(c).enumerate() {
        assert!(
            r.to_bits() == c.to_bits() || (r.is_nan() && c.is_nan()),
            "index {i}: Rust {r:?} ({:016x}), C {c:?} ({:016x})",
            r.to_bits(),
            c.to_bits()
        );
    }
}

type System64 = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

// a[1..] is the lower diagonal. Most cases are diagonally dominant;
// modes 3, 5, and 6 also exercise general, nearly singular, and pivoting cases.
fn system(rng: &mut XorShift, n: usize, case: usize, scale: f64, rhs: f64) -> System64 {
    let mut a = vec![0.0; n];
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    let mut b = vec![0.0; n];
    for k in 0..n {
        a[k] = rng.uniform(-0.8, 0.8) * scale;
        e[k] = rng.uniform(-0.8, 0.8) * scale;
        d[k] = rng.uniform(2.0, 4.0) * scale;
        b[k] = rng.uniform(-4.0, 4.0) * rhs;
        match case % 8 {
            0 => b[k] = 0.0,
            1 => {
                a[k] = 0.0;
                e[k] = -0.0;
            }
            2 => d[k] *= if rng.next_u32() & 1 == 0 { -1.0 } else { 1.0 },
            3 => d[k] = rng.uniform(-1.0, 1.0) * scale,
            4 => {
                if k % 3 == 0 {
                    a[k] = 0.0;
                }
                if k % 5 == 0 {
                    e[k] = 0.0;
                }
            }
            5 => {
                a[k] = -scale;
                d[k] = 2.0 * scale;
                e[k] = -scale;
            }
            _ => {}
        }
    }
    if case % 8 == 6 && n > 1 {
        d[0] = 0.125 * scale;
        a[1] = 2.0 * scale;
        e[0] = 0.25 * scale;
    }
    (a, d, e, b)
}

fn f32vec(v: Vec<f64>) -> Vec<f32> {
    v.into_iter().map(|x| x as f32).collect()
}

fn check_tridid(mut a: Vec<f64>, mut d: Vec<f64>, mut e: Vec<f64>, mut b: Vec<f64>) {
    let n = d.len();
    let saved = (a.clone(), d.clone(), e.clone(), b.clone());
    // The output suffix must remain untouched.
    let mut c_out = vec![12345.0; n + 2];
    let mut r_out = c_out.clone();
    rust::tridid(&a, &d, &e, &b, &mut r_out);
    unsafe {
        tridid(
            n as i32,
            a.as_mut_ptr(),
            d.as_mut_ptr(),
            e.as_mut_ptr(),
            b.as_mut_ptr(),
            c_out.as_mut_ptr(),
        );
    }
    equal64(&r_out, &c_out);
    equal64(&a, &saved.0);
    equal64(&d, &saved.1);
    equal64(&e, &saved.2);
    equal64(&b, &saved.3);
}

fn check_tridif(mut a: Vec<f32>, mut d: Vec<f32>, mut e: Vec<f32>, mut b: Vec<f32>) {
    let n = d.len();
    let saved = (a.clone(), d.clone(), e.clone(), b.clone());
    let mut c_out = vec![12345.0; n + 2];
    let mut r_out = c_out.clone();
    rust::tridif(&a, &d, &e, &b, &mut r_out);
    unsafe {
        tridif(
            n as i32,
            a.as_mut_ptr(),
            d.as_mut_ptr(),
            e.as_mut_ptr(),
            b.as_mut_ptr(),
            c_out.as_mut_ptr(),
        );
    }
    equal32(&r_out, &c_out);
    equal32(&a, &saved.0);
    equal32(&d, &saved.1);
    equal32(&e, &saved.2);
    equal32(&b, &saved.3);
}

#[test]
fn tridid_matches_c() {
    let mut rng = XorShift::new(0x747269646964);
    for n in SIZES {
        for case in 0..CASES_PER_SIZE {
            let scale = 2.0_f64.powi([-400, -100, 0, 100, 400][case % 5]);
            let rhs = 2.0_f64.powi([-400, -100, 0, 100, 400][(case / 5) % 5]);
            let (mut a, d, mut e, b) = system(&mut rng, n, case, scale, rhs);
            a[0] = f64::NAN; // C ignores this entry, even for n=1.
            e[n - 1] = f64::NAN;
            check_tridid(a, d, e, b);
        }
    }
    for x in [
        0.0,
        -0.0,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ] {
        check_tridid(vec![0.0], vec![x], vec![0.0], vec![-1.0]);
        check_tridid(
            vec![0.0, -1.0, x],
            vec![1.0, 1.0, 2.0],
            vec![1.0, x, 0.0],
            vec![x, -0.0, 1.0],
        );
    }
}

#[test]
fn tridif_matches_c() {
    let mut rng = XorShift::new(0x747269646966);
    for n in SIZES {
        for case in 0..CASES_PER_SIZE {
            let scale = 2.0_f64.powi([-60, -20, 0, 20, 60][case % 5]);
            let rhs = 2.0_f64.powi([-60, -20, 0, 20, 60][(case / 5) % 5]);
            let (a, d, e, b) = system(&mut rng, n, case, scale, rhs);
            let mut a = f32vec(a);
            let mut e = f32vec(e);
            a[0] = f32::NAN;
            e[n - 1] = f32::NAN;
            check_tridif(a, f32vec(d), e, f32vec(b));
        }
    }
    for x in [
        0.0,
        -0.0,
        f32::MIN_POSITIVE,
        f32::from_bits(1),
        f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ] {
        check_tridif(vec![0.0], vec![x], vec![0.0], vec![-1.0]);
        check_tridif(
            vec![0.0, -1.0, x],
            vec![1.0, 1.0, 2.0],
            vec![1.0, x, 0.0],
            vec![x, -0.0, 1.0],
        );
    }
}

fn check_tripd(mut d: Vec<f32>, mut e: Vec<f32>, mut b: Vec<f32>) {
    let n = d.len();
    // Exercise both the minimal n-1 superdiagonal and extra unused storage.
    b.extend_from_slice(&[8765.0, -8765.0]);
    let (mut rd, mut re, mut rb) = (d.clone(), e.clone(), b.clone());
    rust::tripd(&mut rd, &mut re, &mut rb);
    unsafe {
        tripd(d.as_mut_ptr(), e.as_mut_ptr(), b.as_mut_ptr(), n as i32);
    }
    equal32(&rd, &d);
    equal32(&re, &e);
    equal32(&rb, &b);
}

#[test]
fn tripd_matches_c() {
    let mut rng = XorShift::new(0x7472697064);
    for n in SIZES {
        for case in 0..CASES_PER_SIZE {
            let scale = 2.0_f32.powi([-60, -20, 0, 20, 60][case % 5]);
            let rhs = 2.0_f32.powi([-60, -20, 0, 20, 60][(case / 5) % 5]);
            let mut d = vec![0.0; n];
            let mut e = vec![0.0; n];
            let mut b = vec![0.0; n];
            for k in 0..n {
                // Strict diagonal dominance with d>0 ensures positive definiteness.
                d[k] = rng.uniform_f32(2.0, 4.0) * scale;
                e[k] = rng.uniform_f32(-0.8, 0.8) * scale;
                b[k] = rng.uniform_f32(-4.0, 4.0) * rhs;
                match case % 4 {
                    0 => b[k] = 0.0,
                    1 => e[k] = -0.0,
                    2 => {
                        d[k] = 2.0 * scale;
                        e[k] = -scale;
                    }
                    _ => {}
                }
            }
            if case % 2 == 0 {
                e.truncate(n - 1);
            } else {
                e[n - 1] = f32::NAN; // This last superdiagonal is untouched.
            }
            check_tripd(d, e, b);
        }
    }
    // Also match the C arithmetic outside the positive definite contract.
    for x in [
        0.0,
        -0.0,
        -2.0,
        f32::MIN_POSITIVE,
        f32::from_bits(1),
        f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ] {
        check_tripd(vec![x], vec![], vec![-1.0]);
        check_tripd(vec![1.0, 1.0, x], vec![1.0, x], vec![x, -0.0, 1.0]);
    }
}

fn check_tripp(mut d: Vec<f32>, mut e: Vec<f32>, mut c: Vec<f32>, mut b: Vec<f32>) {
    let n = d.len();
    // C unconditionally accesses b[-1], d[-1], e[-1] when n==1.
    // Provide initialized prefix cells inside the same allocations, so every
    // C access is valid. Compare all logical entries and the trailing guards;
    // Rust intentionally never exposes or writes C's extra prefix cell.
    d.insert(0, 1.0);
    e.insert(0, 0.0);
    c.insert(0, 0.0);
    b.insert(0, 123.0);
    e.push(5432.0);
    c.push(-5432.0);
    b.push(6789.0);
    let (mut rd, mut re, mut rc, mut rb) = (
        d[1..].to_vec(),
        e[1..].to_vec(),
        c[1..].to_vec(),
        b[1..].to_vec(),
    );
    rust::tripp(&mut rd, &mut re, &mut rc, &mut rb);
    unsafe {
        tripp(
            n as i32,
            d.as_mut_ptr().add(1),
            e.as_mut_ptr().add(1),
            c.as_mut_ptr().add(1),
            b.as_mut_ptr().add(1),
        );
    }
    equal32(&rd, &d[1..]);
    equal32(&re, &e[1..]);
    equal32(&rc, &c[1..]);
    equal32(&rb, &b[1..]);
    assert_eq!(d[0], 1.0);
    assert_eq!(e[0], 0.0);
    assert_eq!(c[0], 0.0);
}

#[test]
fn tripp_matches_c() {
    let mut rng = XorShift::new(0x7472697070);
    for n in SIZES {
        for case in 0..CASES_PER_SIZE {
            let scale = 2.0_f64.powi([-60, -20, 0, 20, 60][case % 5]);
            let rhs = 2.0_f64.powi([-60, -20, 0, 20, 60][(case / 5) % 5]);
            let (c, d, e, b) = system(&mut rng, n, case, scale, rhs);
            check_tripp(f32vec(d), f32vec(e), f32vec(c), f32vec(b));
        }
    }
    for x in [
        0.0,
        -0.0,
        f32::MIN_POSITIVE,
        f32::from_bits(1),
        f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ] {
        check_tripp(vec![x], vec![0.0], vec![0.0], vec![-1.0]);
        // A single zero pivot (x==0) or a nonfinite arithmetic path.
        check_tripp(vec![x, 2.0], vec![1.0, 0.0], vec![7.0, -1.0], vec![1.0, x]);
    }
}

#[test]
fn tripp_retains_noop_pivot_and_last_superdiagonal_write() {
    let (mut d, mut e, mut c, mut b) = (
        vec![1.0, 3.0],
        vec![2.0, -0.0],
        vec![99.0, 4.0],
        vec![5.0, 6.0],
    );
    // |d[0]| < |c[1]| takes the C exch branch, but no values are swapped.
    check_tripp(d.clone(), e.clone(), c.clone(), b.clone());
    rust::tripp(&mut d, &mut e, &mut c, &mut b);
    equal32(&d, &[1.0, -5.0]);
    equal32(&e, &[2.0, -0.0]);
    equal32(&c, &[0.0, 4.0]);
    equal32(&b, &[-0.5999999, 2.8]);

    // Even the unused last superdiagonal must be updated: subtracting a
    // negative zero from negative zero changes its sign to positive zero.
    let (mut d, mut e, mut c, mut b) = (
        vec![1.0, 3.0],
        vec![2.0, -0.0],
        vec![99.0, -4.0],
        vec![5.0, 6.0],
    );
    check_tripp(d.clone(), e.clone(), c.clone(), b.clone());
    rust::tripp(&mut d, &mut e, &mut c, &mut b);
    assert_eq!(e[1].to_bits(), 0.0_f32.to_bits());
}
