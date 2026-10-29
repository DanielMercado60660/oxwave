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

//! Differential tests for cwp/lib/dpfafft.c.
//! Author: Dave Hale, Colorado School of Mines, 04/27/89.
//! Revised by Baoniu Han to handle double precision, 12/14/98.
//! Cost table rebuilt by Jack Cohen, 05/03/94, on advice of David and
//! Gregory Chudnovsky.

use oxwave_cref::XorShift;
use oxwave_cwp::complex::DComplex;
use oxwave_f_dpfafft as port;

extern "C" {
    fn npfa_d(nmin: i32) -> i32;
    fn npfao_d(nmin: i32, nmax: i32) -> i32;
    fn npfar_d(nmin: i32) -> i32;
    fn npfaro_d(nmin: i32, nmax: i32) -> i32;
    fn pfacc_d(isign: i32, n: i32, z: *mut DComplex);
    fn pfamcc_d(isign: i32, n: i32, nt: i32, k: i32, kt: i32, z: *mut DComplex);
    fn pfacr_d(isign: i32, n: i32, cz: *mut DComplex, rz: *mut f64);
    fn pfarc_d(isign: i32, n: i32, rz: *mut f64, cz: *mut DComplex);
    fn pfa2cc_d(isign: i32, idim: i32, n1: i32, n2: i32, z: *mut DComplex);
    fn pfa2cr_d(isign: i32, idim: i32, n1: i32, n2: i32, cz: *mut DComplex, rz: *mut f64);
    fn pfa2rc_d(isign: i32, idim: i32, n1: i32, n2: i32, rz: *mut f64, cz: *mut DComplex);
}

// Enumerate all legal lengths independently of the C and Rust lookup tables.
// Higher powers of two or three are not valid prime-factor transform lengths.
fn lengths() -> Vec<usize> {
    let mut sizes = Vec::new();
    for a in [1, 2, 4, 8, 16] {
        for b in [1, 3, 9] {
            for c in [1, 5] {
                for d in [1, 7] {
                    for e in [1, 11] {
                        for f in [1, 13] {
                            sizes.push(a * b * c * d * e * f);
                        }
                    }
                }
            }
        }
    }
    sizes.sort_unstable();
    assert_eq!(sizes.len(), 240);
    sizes
}

// Each supported length and sign in both 2-D orientations, then thousands
// of random lengths/signs. Largest real length is 2*720720 = 1441440.
fn cases(r: &mut XorShift) -> Vec<(usize, i32)> {
    let sizes = lengths();
    let mut result = Vec::new();
    for &n in &sizes {
        for sign in [-1, -1, 1, 1] {
            result.push((n, sign));
        }
    }
    let small: Vec<_> = sizes.into_iter().filter(|&n| n <= 8192).collect();
    for _ in 0..2500 {
        let n = small[r.next_u32() as usize % small.len()];
        let sign = [-7, -1, 0, 1, 9][r.next_u32() as usize % 5];
        result.push((n, sign));
    }
    result
}

fn value(r: &mut XorShift, iteration: usize) -> f64 {
    if iteration % 19 == 0 {
        0.0
    } else {
        (r.next_u32() as f64 / u32::MAX as f64 - 0.5)
            * 200.0
            * [1e-120, 1e-6, 1.0, 1e6, 1e120][iteration % 5]
    }
}
fn complex(r: &mut XorShift, len: usize, i: usize) -> Vec<DComplex> {
    (0..len)
        .map(|_| DComplex {
            r: value(r, i),
            i: value(r, i),
        })
        .collect()
}
fn real(r: &mut XorShift, len: usize, i: usize) -> Vec<f64> {
    (0..len).map(|_| value(r, i)).collect()
}
fn check(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        assert!(
            (x - y).abs() <= 1e-5 * scale,
            "index {i}: rust={x}, c={y}, scale={scale}"
        );
    }
}
fn check_complex(a: &[DComplex], b: &[DComplex]) {
    assert_eq!(a.len(), b.len());
    let scale = b
        .iter()
        .fold(0.0_f64, |m, z| m.max(z.r.abs()).max(z.i.abs()));
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        for (x, y) in [(x.r, y.r), (x.i, y.i)] {
            assert!(x.is_finite() && y.is_finite());
            assert!(
                (x - y).abs() <= 1e-5 * scale,
                "index {i}: rust={x}, c={y}, scale={scale}"
            );
        }
    }
}

// No libm calls occur in the complex kernels. Require bit equality,
// including signed zeros, untouched guard elements and strided gaps.
fn check_complex_exact(a: &[DComplex], b: &[DComplex]) {
    assert_eq!(a.len(), b.len());
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        assert_eq!(
            x.r.to_bits(),
            y.r.to_bits(),
            "real index {i}: {x:?} vs {y:?}"
        );
        assert_eq!(
            x.i.to_bits(),
            y.i.to_bits(),
            "imag index {i}: {x:?} vs {y:?}"
        );
    }
}

#[test]
fn pfacc_d_matches_c() {
    let mut r = XorShift::new(177);
    for (i, (n, sign)) in cases(&mut r).into_iter().enumerate() {
        let mut a = complex(&mut r, n + 2, i);
        let mut b = a.clone();
        port::pfacc_d(sign, &mut a[1..n + 1]);
        unsafe { pfacc_d(sign, n as i32, b.as_mut_ptr().add(1)) };
        check_complex_exact(&a, &b);
    }
}
#[test]
fn pfamcc_d_matches_c() {
    let mut r = XorShift::new(178);
    for (i, (n, sign)) in cases(&mut r).into_iter().enumerate() {
        let nt = if n > 8192 { 1 + i % 2 } else { i % 6 };
        let (k, kt) = match i % 6 {
            0 => (nt + 1, 1),
            1 => (1, n + 2),
            2 => (2, 2 * n + 3),
            3 => (3 * nt + 1, 3),
            4 => (0, 2),
            _ => (1, 0),
        };
        let len = if nt == 0 {
            0
        } else {
            (n - 1) * k + (nt - 1) * kt + 1
        };
        let mut a = complex(&mut r, len + 2, i);
        let mut b = a.clone();
        port::pfamcc_d(
            sign,
            n as i32,
            nt as i32,
            k as i32,
            kt as i32,
            &mut a[1..len + 1],
        );
        unsafe {
            pfamcc_d(
                sign,
                n as i32,
                nt as i32,
                k as i32,
                kt as i32,
                b.as_mut_ptr().add(1),
            )
        };
        check_complex_exact(&a, &b);
    }
}
#[test]
fn pfarc_d_matches_c() {
    let mut r = XorShift::new(179);
    for (i, (half, sign)) in cases(&mut r).into_iter().enumerate() {
        let n = 2 * half;
        let mut input = real(&mut r, n, i);
        let mut a = complex(&mut r, n / 2 + 1, i);
        let mut b = a.clone();
        port::pfarc_d(sign, &input, &mut a);
        unsafe { pfarc_d(sign, n as i32, input.as_mut_ptr(), b.as_mut_ptr()) };
        check_complex(&a, &b);
    }
}
#[test]
fn pfacr_d_matches_c() {
    let mut r = XorShift::new(180);
    for (i, (half, sign)) in cases(&mut r).into_iter().enumerate() {
        let n = 2 * half;
        let mut input = complex(&mut r, n / 2 + 1, i);
        let mut a = real(&mut r, n, i);
        let mut b = a.clone();
        port::pfacr_d(sign, &input, &mut a);
        unsafe { pfacr_d(sign, n as i32, input.as_mut_ptr(), b.as_mut_ptr()) };
        check(&a, &b);
    }
}

#[test]
fn pfa2cc_d_matches_c() {
    let mut r = XorShift::new(181);
    for (i, (n, sign)) in cases(&mut r).into_iter().enumerate() {
        let idim = 1 + (i % 2) as i32;
        // The other dimension need not be a valid FFT length.
        let other = if n > 8192 {
            2
        } else {
            1 + r.next_u32() as usize % 17
        };
        let (n1, n2) = if idim == 1 { (n, other) } else { (other, n) };
        let mut a = complex(&mut r, n1 * n2 + 2, i);
        let mut b = a.clone();
        port::pfa2cc_d(sign, idim, n1 as i32, n2 as i32, &mut a[1..n1 * n2 + 1]);
        unsafe { pfa2cc_d(sign, idim, n1 as i32, n2 as i32, b.as_mut_ptr().add(1)) };
        check_complex_exact(&a, &b);
    }
}
#[test]
fn pfa2rc_d_matches_c() {
    let mut r = XorShift::new(182);
    for (i, (half, sign)) in cases(&mut r).into_iter().enumerate() {
        let idim = 1 + (i % 2) as i32;
        let n = 2 * half;
        let other = if half > 8192 {
            2
        } else {
            1 + r.next_u32() as usize % 17
        };
        let (n1, n2) = if idim == 1 { (n, other) } else { (other, n) };
        let mut input = real(&mut r, n1 * n2, i);
        let len = if idim == 1 {
            (n1 + 2) * n2 / 2
        } else {
            n1 * (n2 + 2) / 2
        };
        let mut a = complex(&mut r, len, i);
        let mut b = a.clone();
        port::pfa2rc_d(sign, idim, n1 as i32, n2 as i32, &input, &mut a);
        unsafe {
            pfa2rc_d(
                sign,
                idim,
                n1 as i32,
                n2 as i32,
                input.as_mut_ptr(),
                b.as_mut_ptr(),
            )
        };
        check_complex(&a, &b);
    }
}
#[test]
fn pfa2cr_d_matches_c() {
    let mut r = XorShift::new(183);
    for (i, (half, sign)) in cases(&mut r).into_iter().enumerate() {
        let idim = 1 + (i % 2) as i32;
        let n = 2 * half;
        let other = if half > 8192 {
            2
        } else {
            1 + r.next_u32() as usize % 17
        };
        let (n1, n2) = if idim == 1 { (n, other) } else { (other, n) };
        let len = if idim == 1 {
            (n1 + 2) * n2 / 2
        } else {
            n1 * (n2 + 2) / 2
        };
        let mut input = complex(&mut r, len, i);
        let mut a = real(&mut r, n1 * n2, i);
        let mut b = a.clone();
        port::pfa2cr_d(sign, idim, n1 as i32, n2 as i32, &input, &mut a);
        unsafe {
            pfa2cr_d(
                sign,
                idim,
                n1 as i32,
                n2 as i32,
                input.as_mut_ptr(),
                b.as_mut_ptr(),
            )
        };
        check(&a, &b);
    }
}

// Requests on either side of every supported length, full-width signed
// integers, and bounds above the transform limit. Integer results are exact.
fn requests(r: &mut XorShift, real_lengths: bool) -> Vec<i32> {
    let mut values = vec![i32::MIN, i32::MIN + 1, -1, 0, 1, 2, i32::MAX - 1, i32::MAX];
    for n in lengths() {
        let n = (n * if real_lengths { 2 } else { 1 }) as i32;
        values.extend_from_slice(&[n - 1, n, n + 1]);
    }
    for _ in 0..30_000 {
        values.push(r.next_u32() as i32);
        values.push((r.next_u32() % 2_000_000) as i32 - 100_000);
    }
    values
}

#[test]
fn npfa_d_matches_c() {
    let mut r = XorShift::new(0x12345678);
    for n in requests(&mut r, false) {
        assert_eq!(port::npfa_d(n), unsafe { npfa_d(n) }, "n={n}");
    }
}

#[test]
fn npfao_d_matches_c() {
    let mut r = XorShift::new(0x23456789);
    for n in requests(&mut r, false) {
        for m in [
            n,
            n.saturating_sub(1),
            n.saturating_add(1),
            r.next_u32() as i32,
        ] {
            assert_eq!(
                port::npfao_d(n, m),
                unsafe { npfao_d(n, m) },
                "n={n}, m={m}"
            );
        }
    }
}

#[test]
fn npfar_d_matches_c() {
    let mut r = XorShift::new(0x34567890);
    for n in requests(&mut r, true) {
        assert_eq!(port::npfar_d(n), unsafe { npfar_d(n) }, "n={n}");
    }
}

#[test]
fn npfaro_d_matches_c() {
    let mut r = XorShift::new(0x45678901);
    for n in requests(&mut r, true) {
        for m in [
            n,
            n.saturating_sub(1),
            n.saturating_add(1),
            r.next_u32() as i32,
        ] {
            assert_eq!(
                port::npfaro_d(n, m),
                unsafe { npfaro_d(n, m) },
                "n={n}, m={m}"
            );
        }
    }
}

// Odd real lengths are outside C's documented even-length recommendation.
// Padded, initialized allocations make the C indexing defined and let us
// verify its actual behavior (including dependence on workspace contents).
// Each selected n has a supported floor(n/2), including zero for n=1.
const ODD: &[usize] = &[1, 3, 5, 7, 9, 11, 15, 17, 19, 23, 27, 31, 33];

#[test]
fn pfarc_d_odd_and_unit_matches_padded_c() {
    let mut r = XorShift::new(190);
    for i in 0..2600 {
        let n = ODD[i % ODD.len()];
        let sign = [-3, -1, 0, 1, 5][i % 5];
        let mut input = real(&mut r, n, i);
        // C accesses double slots n and n+1, beyond its documented n/2+1
        // complex elements when n is odd.
        let nc = (n + 3) / 2;
        let mut a = complex(&mut r, nc + 2, i);
        let mut b = a.clone();
        port::pfarc_d(sign, &input, &mut a[1..nc + 1]);
        unsafe { pfarc_d(sign, n as i32, input.as_mut_ptr(), b.as_mut_ptr().add(1)) };
        check_complex(&a[1..nc + 1], &b[1..nc + 1]);
        check_complex_exact(&[a[0], a[nc + 1]], &[b[0], b[nc + 1]]);
    }
}

#[test]
fn pfacr_d_odd_and_unit_matches_padded_c() {
    let mut r = XorShift::new(191);
    for i in 0..2600 {
        let n = ODD[i % ODD.len()];
        let sign = [-3, -1, 0, 1, 5][i % 5];
        let mut input = complex(&mut r, n / 2 + 1, i);
        // C writes rz[1] even when n=1. Give C a spare slot; the Rust API
        // preserves the in-bounds result without this out-of-slice write.
        let nr = n.max(2);
        let mut a = real(&mut r, nr + 2, i);
        let mut b = a.clone();
        port::pfacr_d(sign, &input, &mut a[1..n + 1]);
        unsafe { pfacr_d(sign, n as i32, input.as_mut_ptr(), b.as_mut_ptr().add(1)) };
        check(&a[1..n + 1], &b[1..n + 1]);
        assert_eq!(a[0].to_bits(), b[0].to_bits());
        assert_eq!(a[nr + 1].to_bits(), b[nr + 1].to_bits());
    }
}

#[test]
fn pfa2rc_d_odd_and_unit_matches_padded_c() {
    let mut r = XorShift::new(192);
    for i in 0..2600 {
        let n = ODD[(i / 2) % ODD.len()];
        let idim = 1 + (i % 2) as i32;
        let other = 1 + r.next_u32() as usize % 17;
        let (n1, n2) = if idim == 1 { (n, other) } else { (other, n) };
        let sign = [-3, -1, 0, 1, 5][i % 5];
        let mut input = real(&mut r, n1 * n2, i);
        let slots = if idim == 1 {
            (n1 + 2) * n2
        } else {
            n1 * (n2 + 2)
        };
        let nc = (slots + 1) / 2;
        let mut a = complex(&mut r, nc + 2, i);
        let mut b = a.clone();
        port::pfa2rc_d(sign, idim, n1 as i32, n2 as i32, &input, &mut a[1..nc + 1]);
        unsafe {
            pfa2rc_d(
                sign,
                idim,
                n1 as i32,
                n2 as i32,
                input.as_mut_ptr(),
                b.as_mut_ptr().add(1),
            )
        };
        check_complex(&a[1..nc + 1], &b[1..nc + 1]);
        check_complex_exact(&[a[0], a[nc + 1]], &[b[0], b[nc + 1]]);
    }
}

#[test]
fn pfa2cr_d_odd_and_unit_matches_padded_c() {
    let mut r = XorShift::new(193);
    for i in 0..2600 {
        let n = ODD[(i / 2) % ODD.len()];
        let idim = 1 + (i % 2) as i32;
        let other = 1 + r.next_u32() as usize % 17;
        let (n1, n2) = if idim == 1 { (n, other) } else { (other, n) };
        let sign = [-3, -1, 0, 1, 5][i % 5];
        let slots = if idim == 1 {
            (n1 + 2) * n2
        } else {
            n1 * (n2 + 2)
        };
        let mut input = complex(&mut r, slots / 2, i);
        // Unit fast dimension writes one extra real slot; an odd slow
        // dimension unmerges one extra row, including caller workspace.
        let nr = if idim == 1 {
            n1 * n2 + usize::from(n1 == 1)
        } else {
            n1 * (n2 + 1)
        };
        let mut a = real(&mut r, nr + 2, i);
        let mut b = a.clone();
        port::pfa2cr_d(sign, idim, n1 as i32, n2 as i32, &input, &mut a[1..nr + 1]);
        unsafe {
            pfa2cr_d(
                sign,
                idim,
                n1 as i32,
                n2 as i32,
                input.as_mut_ptr(),
                b.as_mut_ptr().add(1),
            )
        };
        check(&a[1..nr + 1], &b[1..nr + 1]);
        assert_eq!(a[0].to_bits(), b[0].to_bits());
        assert_eq!(a[nr + 1].to_bits(), b[nr + 1].to_bits());
    }
}
