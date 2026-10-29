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

//! Differential tests against cwp/lib/pfafft.c.
//! Credits: Dave Hale, Colorado School of Mines, 1989–1991;
//! FFT cost table rebuilt by Jack Cohen, 05/03/94, on advice of
//! David and Gregory Chudnovsky. FFT algorithms reference C. Temperton.

use oxwave_cref::XorShift;
use oxwave_pfa_fft::{self as rust, Complex};

mod c {
    use super::Complex;
    extern "C" {
        pub fn npfa(nmin: i32) -> i32;
        pub fn npfao(nmin: i32, nmax: i32) -> i32;
        pub fn npfar(nmin: i32) -> i32;
        pub fn npfaro(nmin: i32, nmax: i32) -> i32;
        pub fn pfacc(isign: i32, n: i32, z: *mut Complex);
        pub fn pfacr(isign: i32, n: i32, cz: *mut Complex, rz: *mut f32);
        pub fn pfarc(isign: i32, n: i32, rz: *mut f32, cz: *mut Complex);
        pub fn pfamcc(isign: i32, n: i32, nt: i32, k: i32, kt: i32, z: *mut Complex);
        pub fn pfa2cc(isign: i32, idim: i32, n1: i32, n2: i32, z: *mut Complex);
        pub fn pfa2cr(isign: i32, idim: i32, n1: i32, n2: i32, cz: *mut Complex, rz: *mut f32);
        pub fn pfa2rc(isign: i32, idim: i32, n1: i32, n2: i32, rz: *mut f32, cz: *mut Complex);
    }
}

// Derive the admissible lengths independently of either implementation.
fn lengths() -> Vec<usize> {
    let mut out = Vec::new();
    for p2 in [1, 2, 4, 8, 16] {
        for p3 in [1, 3, 9] {
            for p5 in [1, 5] {
                for p7 in [1, 7] {
                    for p11 in [1, 11] {
                        for p13 in [1, 13] {
                            out.push(p2 * p3 * p5 * p7 * p11 * p13);
                        }
                    }
                }
            }
        }
    }
    out.sort_unstable();
    assert_eq!(out.len(), 240);
    out
}

fn bounds(rng: &mut XorShift) -> Vec<i32> {
    let mut v = vec![
        i32::MIN,
        i32::MIN + 1,
        -100,
        -2,
        -1,
        0,
        1,
        2,
        3,
        7,
        11,
        13,
        17,
        31,
        32,
        65536,
        720719,
        720720,
        720721,
        1441439,
        1441440,
        1441441,
        i32::MAX - 1,
        i32::MAX,
    ];
    for n in lengths() {
        for k in [-1, 0, 1] {
            v.push(n as i32 + k);
            v.push(2 * n as i32 + k);
        }
    }
    for _ in 0..20_000 {
        v.push(rng.range_i32(-2_000_000, 3_000_000));
        v.push(rng.next_u32() as i32);
    }
    v
}

#[test]
fn npfa_matches_c() {
    let mut rng = XorShift::new(1);
    for nmin in bounds(&mut rng) {
        assert_eq!(rust::npfa(nmin), unsafe { c::npfa(nmin) }, "nmin={nmin}");
    }
}

#[test]
fn npfar_matches_c() {
    let mut rng = XorShift::new(2);
    for nmin in bounds(&mut rng) {
        // INT_MAX+1 is undefined in C. Exercise that boundary against the
        // original complex selector using the mathematically rounded bound.
        let expected = unsafe {
            if nmin == i32::MAX {
                2 * c::npfa(((i64::from(nmin) + 1) / 2) as i32)
            } else {
                c::npfar(nmin)
            }
        };
        assert_eq!(rust::npfar(nmin), expected, "nmin={nmin}");
    }
}

#[test]
fn npfao_matches_c() {
    let mut rng = XorShift::new(3);
    for nmin in bounds(&mut rng) {
        // Includes reversed/empty intervals, equality, saturation, cost ties.
        for nmax in [
            nmin,
            nmin.wrapping_sub(1),
            nmin.wrapping_mul(2),
            rng.next_u32() as i32,
        ] {
            assert_eq!(
                rust::npfao(nmin, nmax),
                unsafe { c::npfao(nmin, nmax) },
                "nmin={nmin} nmax={nmax}"
            );
        }
    }
}

#[test]
fn npfaro_matches_c() {
    let mut rng = XorShift::new(4);
    for nmin in bounds(&mut rng) {
        for nmax in [
            nmin,
            nmin.wrapping_sub(1),
            nmin.wrapping_mul(2),
            rng.next_u32() as i32,
        ] {
            let expected = unsafe {
                if nmin == i32::MAX || nmax == i32::MAX {
                    // Avoid undefined signed addition in either C bound.
                    2 * c::npfao(
                        ((i64::from(nmin) + 1) / 2) as i32,
                        ((i64::from(nmax) + 1) / 2) as i32,
                    )
                } else {
                    c::npfaro(nmin, nmax)
                }
            };
            assert_eq!(
                rust::npfaro(nmin, nmax),
                expected,
                "nmin={nmin} nmax={nmax}"
            );
        }
    }
}

#[test]
fn real_selectors_saturate_at_integer_limit() {
    let mut rng = XorShift::new(18);
    for case in 0..4096 {
        let high = i32::MAX - case;
        assert_eq!(rust::npfar(high), 1_441_440);
        for other in [i32::MIN, -1, 0, 1, i32::MAX, rng.next_u32() as i32] {
            assert_eq!(rust::npfaro(high, other), 1_441_440);
        }
        // A huge upper bound must still allow cost optimization for ordinary
        // lower bounds. INT_MAX-1 is a defined, table-equivalent C bound.
        let low = rng.range_i32(-100, 1_441_441);
        assert_eq!(rust::npfaro(low, high), unsafe {
            c::npfaro(low, i32::MAX - 1)
        });
    }
}

#[derive(Default)]
struct Check {
    samples: usize,
    different_bits: usize,
    worst_relative: f64,
}

impl Check {
    fn real(&mut self, got: &[f32], expected: &[f32], context: &str) {
        assert_eq!(got.len(), expected.len());
        let mut scale = 0.0f64;
        for &v in expected {
            assert!(v.is_finite(), "nonfinite C result: {context}");
            if (v as f64).abs() > scale {
                scale = (v as f64).abs();
            }
        }
        for (i, (&a, &b)) in got.iter().zip(expected).enumerate() {
            assert!(a.is_finite(), "nonfinite Rust result: {context} slot={i}");
            let delta = (a as f64 - b as f64).abs();
            // A zero C vector must match exactly; no absolute floor masks errors.
            assert!(
                delta <= 1e-5 * scale,
                "{context} slot={i}: Rust={a:e} C={b:e} scale={scale:e}"
            );
            if scale > 0.0 {
                self.worst_relative = self.worst_relative.max(delta / scale);
            }
            self.samples += 1;
            self.different_bits += usize::from(a.to_bits() != b.to_bits());
        }
    }

    fn complex(&mut self, got: &[Complex], expected: &[Complex], context: &str) {
        let a: Vec<_> = got.iter().flat_map(|z| [z.r, z.i]).collect();
        let b: Vec<_> = expected.iter().flat_map(|z| [z.r, z.i]).collect();
        self.real(&a, &b, context);
    }

    fn report(&self, name: &str) {
        eprintln!(
            "{name}: {} float components, {} differing bits, worst error/max={:.3e}",
            self.samples, self.different_bits, self.worst_relative
        );
    }
}

fn sign(case: usize) -> i32 {
    [1, -1, 0, 3, -7][case % 5]
}

fn input(rng: &mut XorShift, n: usize, case: usize) -> Vec<f32> {
    let scale = [1.0, 1e-30, 1e20, 1e-10, 1000.0, 1e-35, 1e25][case % 7];
    (0..n)
        .map(|i| match case % 19 {
            0 => 0.0,
            1 => {
                if i % 2 == 0 {
                    -0.0
                } else {
                    0.0
                }
            }
            2 => {
                if i == n / 2 {
                    -scale
                } else {
                    0.0
                }
            }
            3 => -scale,
            4 => {
                if i % 2 == 0 {
                    scale
                } else {
                    -scale
                }
            }
            5 => -rng.uniform_f32(0.0, 1.0) * scale,
            _ => rng.uniform_f32(-1.0, 1.0) * scale,
        })
        .collect()
}

fn complex_input(rng: &mut XorShift, n: usize, case: usize) -> Vec<Complex> {
    input(rng, 2 * n, case)
        .chunks_exact(2)
        .map(|z| Complex { r: z[0], i: z[1] })
        .collect()
}

fn select_length(rng: &mut XorShift, all: &[usize], case: usize) -> usize {
    // Every admissible length, including 720720, then thousands of varied small/medium
    // transforms. These exercise all rotation factors and factor combinations.
    if case < all.len() {
        all[case]
    } else {
        all[rng.next_u32() as usize % 185]
    }
}

#[test]
fn pfacc_matches_c() {
    let mut rng = XorShift::new(5);
    let all = lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let n = select_length(&mut rng, &all, case);
        let s = sign(case);
        let mut got = complex_input(&mut rng, n, case);
        let mut expected = got.clone();
        unsafe {
            c::pfacc(s, n as i32, expected.as_mut_ptr());
        }
        rust::pfacc(s, &mut got);
        check.complex(
            &got,
            &expected,
            &format!("pfacc n={n} sign={s} case={case}"),
        );
    }
    check.report("pfacc");
}

#[test]
fn pfarc_matches_c() {
    let mut rng = XorShift::new(6);
    let all = lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let n = 2 * select_length(&mut rng, &all, case);
        let s = sign(case);
        let mut rz = input(&mut rng, n, case);
        let before = rz.clone();
        let sentinel = Complex {
            r: -987.25,
            i: 321.5,
        };
        let mut got = vec![sentinel; n / 2 + 4];
        let mut expected = got.clone();
        unsafe {
            c::pfarc(s, n as i32, rz.as_mut_ptr(), expected.as_mut_ptr());
        }
        rust::pfarc(s, &before, &mut got);
        check.complex(
            &got[..n / 2 + 1],
            &expected[..n / 2 + 1],
            &format!("pfarc n={n} sign={s} case={case}"),
        );
        assert_eq!(got[n / 2 + 1..], expected[n / 2 + 1..]);
        for (a, b) in rz.iter().zip(&before) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
    check.report("pfarc");
}

#[test]
fn pfacr_matches_c() {
    let mut rng = XorShift::new(7);
    let all = lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let n = 2 * select_length(&mut rng, &all, case);
        let s = sign(case);
        // Arbitrary spectra also cover ignored imaginary DC/Nyquist values.
        let mut cz = complex_input(&mut rng, n / 2 + 1, case);
        let before = cz.clone();
        let mut got = vec![-987.25; n + 3];
        let mut expected = got.clone();
        unsafe {
            c::pfacr(s, n as i32, cz.as_mut_ptr(), expected.as_mut_ptr());
        }
        rust::pfacr(s, &before, &mut got[..n]);
        check.real(
            &got[..n],
            &expected[..n],
            &format!("pfacr n={n} sign={s} case={case}"),
        );
        assert_eq!(got[n..], expected[n..]);
        for (a, b) in cz.iter().zip(&before) {
            assert_eq!(a.r.to_bits(), b.r.to_bits());
            assert_eq!(a.i.to_bits(), b.i.to_bits());
        }
    }
    check.report("pfacr");
}

#[test]
fn pfamcc_matches_c() {
    let mut rng = XorShift::new(8);
    let all = lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let n = select_length(&mut rng, &all, case);
        let s = sign(case);
        let nt = if case < 240 {
            1 + case % 2
        } else {
            rng.range_i32(0, 5) as usize
        };
        // Contiguous, interleaved, and genuinely gapped layouts.
        let gap = rng.range_i32(1, 3) as usize;
        let (k, kt) = match case % 3 {
            0 => (1, n + gap),
            1 => (nt.max(1) + gap, 1),
            _ => (gap, n * gap + 3),
        };
        let len = if nt == 0 {
            0
        } else {
            (n - 1) * k + (nt - 1) * kt + 1
        };
        let mut got = complex_input(&mut rng, len + 3, case);
        let original = got.clone();
        let mut expected = got.clone();
        unsafe {
            c::pfamcc(
                s,
                n as i32,
                nt as i32,
                k as i32,
                kt as i32,
                expected.as_mut_ptr(),
            );
        }
        rust::pfamcc(s, n as i32, nt as i32, k as i32, kt as i32, &mut got);
        check.complex(
            &got,
            &expected,
            &format!("pfamcc n={n} nt={nt} k={k} kt={kt} sign={s}"),
        );
        let mut touched = vec![false; got.len()];
        for t in 0..nt {
            for j in 0..n {
                touched[t * kt + j * k] = true;
            }
        }
        for (j, touched) in touched.into_iter().enumerate() {
            if !touched {
                assert_eq!(got[j].r.to_bits(), original[j].r.to_bits());
                assert_eq!(got[j].i.to_bits(), original[j].i.to_bits());
                assert_eq!(expected[j].r.to_bits(), original[j].r.to_bits());
                assert_eq!(expected[j].i.to_bits(), original[j].i.to_bits());
            }
        }
    }
    check.report("pfamcc");
}

#[test]
fn pfamcc_zero_count_ignores_strides() {
    let mut rng = XorShift::new(19);
    for case in 0..4096 {
        let n = [1, 2, 16, 720720, i32::MAX, rng.next_u32() as i32][case % 6];
        let k = [0, 1, i32::MAX / 2, i32::MAX, rng.next_u32() as i32][case % 5];
        let kt = [i32::MAX, 0, 1, i32::MAX / 2, rng.next_u32() as i32][case % 5];
        let mut got = complex_input(&mut rng, if case % 2 == 0 { 0 } else { 17 }, case);
        let before = got.clone();
        // There is no C oracle for overflowing stride arithmetic. No
        // transforms must leave even a supplied nonempty buffer untouched.
        rust::pfamcc(sign(case), n, 0, k, kt, &mut got);
        for (a, b) in got.iter().zip(&before) {
            assert_eq!(a.r.to_bits(), b.r.to_bits());
            assert_eq!(a.i.to_bits(), b.i.to_bits());
        }
    }
    rust::pfamcc(1, 2, 0, i32::MAX, i32::MAX, &mut []);
}

#[test]
fn pfamcc_wide_unused_strides_match_c() {
    let mut rng = XorShift::new(20);
    let all = lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let s = sign(case);
        let wide = [i32::MAX / 2, i32::MAX / 2 + 1, i32::MAX - 1, i32::MAX][case % 4];
        let n = all[case % 80];
        let k = (case % 3) as i32;
        let len = (n - 1) * k as usize + 1;
        let mut got = complex_input(&mut rng, len + 3, case);
        let mut expected = got.clone();
        // kt is unused by a single transform, but C doubles it and adds it
        // after the final iteration. Use kt=0 for a defined C oracle.
        unsafe { c::pfamcc(s, n as i32, 1, k, 0, expected.as_mut_ptr()) };
        rust::pfamcc(s, n as i32, 1, k, wide, &mut got[..len]);
        check.complex(&got, &expected, "pfamcc wide kt");

        // A unit transform also never uses k, even when several transforms
        // are requested. Both large strides are harmless for n=nt=1.
        let before = got.clone();
        rust::pfamcc(s, 1, 1, wide, wide, &mut got[..1]);
        rust::pfamcc(s, 1, len as i32, wide, 1, &mut got[..len]);
        for (a, b) in got.iter().zip(&before) {
            assert_eq!(a.r.to_bits(), b.r.to_bits());
            assert_eq!(a.i.to_bits(), b.i.to_bits());
        }

        // Even without a radix inverse, k=0 makes every access safe. Keep
        // C's aliasing behavior instead of rejecting these lengths outright.
        let n = [32, 49, 121, 169][case % 4];
        let mut got = [Complex {
            r: rng.range_i32(-16, 16) as f32,
            i: rng.range_i32(-16, 16) as f32,
        }; 2];
        let mut expected = got;
        unsafe { c::pfamcc(s, n, 2, 0, 1, expected.as_mut_ptr()) };
        rust::pfamcc(s, n, 2, 0, 1, &mut got);
        check.complex(&got, &expected, "pfamcc zero sample stride");
    }
    check.report("pfamcc wide strides");
}

#[test]
fn pfamcc_invalid_inputs_exit_cleanly() {
    const CHILD: &str = "OXWAVE_PFAMCC_INVALID_CASE";
    if let Ok(case) = std::env::var(CHILD) {
        let mut z = vec![Complex { r: 1.0, i: -2.0 }; if case == "radix" { 49 } else { 2 }];
        let (n, nt, k, kt) = match case.as_str() {
            "n" => (0, 1, 1, 1),
            "nt" => (2, -1, 1, 1),
            "k" => (2, 1, -1, 1),
            "kt" => (2, 1, 1, -1),
            "storage" => (3, 1, 1, 1),
            "wide_storage" => (i32::MAX, i32::MAX, i32::MAX, i32::MAX),
            "radix" => (49, 1, 1, 1),
            _ => panic!("unknown child case"),
        };
        rust::pfamcc(1, n, nt, k, kt, &mut z);
        panic!("invalid pfamcc input returned");
    }
    for case in ["n", "nt", "k", "kt", "storage", "wide_storage", "radix"] {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "pfamcc_invalid_inputs_exit_cleanly",
                "--nocapture",
            ])
            .env(CHILD, case)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "case={case}");
        let message = if case.ends_with("storage") {
            "pfafft.c: pfamcc: insufficient complex storage"
        } else if case == "radix" {
            "pfafft.c: pfamcc: invalid transform length"
        } else {
            "pfafft.c: pfamcc: invalid dimensions or strides"
        };
        assert_eq!(
            String::from_utf8(output.stderr).unwrap(),
            format!("\n: {message}\n")
        );
    }
}

#[test]
fn fft_invalid_inputs_exit_cleanly() {
    const CHILD: &str = "OXWAVE_PFA_INVALID_CASE";
    if let Ok(case) = std::env::var(CHILD) {
        let mut z = [Complex::default(); 50];
        let mut r = [0.0; 98];
        match case.as_str() {
            "pfacc_empty" => rust::pfacc(1, &mut []),
            "pfacc_radix" => rust::pfacc(1, &mut z[..49]),
            "pfacr_empty" => rust::pfacr(1, &z, &mut []),
            "pfacr_storage" => rust::pfacr(1, &[], &mut r[..2]),
            "pfacr_radix" => rust::pfacr(1, &z, &mut r),
            "pfarc_empty" => rust::pfarc(1, &[], &mut z),
            "pfarc_storage" => rust::pfarc(1, &r[..2], &mut []),
            "pfarc_radix" => rust::pfarc(1, &r, &mut z),
            "pfa2cc_dims" => rust::pfa2cc(1, 1, -1, 2, &mut z),
            "pfa2cc_storage" => rust::pfa2cc(1, 2, i32::MAX, i32::MAX, &mut z),
            "pfa2cc_radix" => rust::pfa2cc(1, 1, 49, 1, &mut z),
            "pfa2cr_dims" => rust::pfa2cr(1, 1, 0, 2, &z, &mut r),
            "pfa2cr_storage" => rust::pfa2cr(1, 2, 1, 1, &z, &mut r[..1]),
            "pfa2cr_radix" => rust::pfa2cr(1, 1, 98, 1, &z, &mut r),
            "pfa2rc_dims" => rust::pfa2rc(1, 2, 2, -1, &r, &mut z),
            "pfa2rc_storage" => rust::pfa2rc(1, 1, 1, 1, &r, &mut z[..1]),
            "pfa2rc_radix" => rust::pfa2rc(1, 1, 98, 1, &r, &mut z),
            _ => panic!("unknown child case"),
        }
        panic!("invalid FFT input returned");
    }
    for (case, message) in [
        ("pfacc_empty", "pfacc: invalid transform length"),
        ("pfacc_radix", "pfacc: invalid transform length"),
        ("pfacr_empty", "pfacr: invalid transform length"),
        ("pfacr_storage", "pfacr: insufficient complex storage"),
        ("pfacr_radix", "pfacc: invalid transform length"),
        ("pfarc_empty", "pfarc: invalid transform length"),
        ("pfarc_storage", "pfarc: insufficient complex storage"),
        ("pfarc_radix", "pfacc: invalid transform length"),
        ("pfa2cc_dims", "pfa2cc: invalid dimensions"),
        ("pfa2cc_storage", "pfa2cc: insufficient complex storage"),
        ("pfa2cc_radix", "pfamcc: invalid transform length"),
        ("pfa2cr_dims", "pfa2cr: invalid dimensions"),
        ("pfa2cr_storage", "pfa2cr: insufficient storage"),
        ("pfa2cr_radix", "pfamcc: invalid transform length"),
        ("pfa2rc_dims", "pfa2rc: invalid dimensions"),
        ("pfa2rc_storage", "pfa2rc: insufficient storage"),
        ("pfa2rc_radix", "pfamcc: invalid transform length"),
    ] {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "fft_invalid_inputs_exit_cleanly", "--nocapture"])
            .env(CHILD, case)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "case={case}");
        assert_eq!(
            String::from_utf8(output.stderr).unwrap(),
            format!("\n: pfafft.c: {message}\n"),
            "case={case}"
        );
    }
}

// Every length along both axes, plus thousands of random rectangles, including
// odd and unit passive dimensions. idim alone selects the transformed axis.
fn shape(rng: &mut XorShift, all: &[usize], case: usize, real: bool) -> (i32, usize, usize) {
    let idim = 1 + (case % 2) as i32;
    let (n, other) = if case < 480 {
        (all[case / 2], 1 + (case / 2) % 2)
    } else {
        (
            all[rng.next_u32() as usize % 140],
            rng.range_i32(1, 9) as usize,
        )
    };
    let n = if real { 2 * n } else { n };
    if idim == 1 {
        (idim, n, other)
    } else {
        (idim, other, n)
    }
}

#[test]
fn pfa2cc_matches_c() {
    let mut rng = XorShift::new(9);
    let all = lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let (idim, n1, n2) = shape(&mut rng, &all, case, false);
        let s = sign(case);
        let mut got = complex_input(&mut rng, n1 * n2 + 3, case);
        let mut expected = got.clone();
        unsafe {
            c::pfa2cc(s, idim, n1 as i32, n2 as i32, expected.as_mut_ptr());
        }
        rust::pfa2cc(s, idim, n1 as i32, n2 as i32, &mut got);
        check.complex(
            &got[..n1 * n2],
            &expected[..n1 * n2],
            &format!("pfa2cc {n1}x{n2} idim={idim} sign={s}"),
        );
        assert_eq!(got[n1 * n2..], expected[n1 * n2..]);
    }
    check.report("pfa2cc");
}

#[test]
fn pfa2rc_matches_c() {
    let mut rng = XorShift::new(10);
    let all = lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let (idim, n1, n2) = shape(&mut rng, &all, case, true);
        let s = sign(case);
        let nc = if idim == 1 {
            (n1 / 2 + 1) * n2
        } else {
            n1 * (n2 / 2 + 1)
        };
        let mut rz = input(&mut rng, n1 * n2, case);
        let original = rz.clone();
        let sentinel = Complex {
            r: -987.25,
            i: 321.5,
        };
        let mut got = vec![sentinel; nc + 3];
        let mut expected = got.clone();
        unsafe {
            c::pfa2rc(
                s,
                idim,
                n1 as i32,
                n2 as i32,
                rz.as_mut_ptr(),
                expected.as_mut_ptr(),
            );
        }
        rust::pfa2rc(s, idim, n1 as i32, n2 as i32, &original, &mut got);
        check.complex(
            &got[..nc],
            &expected[..nc],
            &format!("pfa2rc {n1}x{n2} idim={idim} sign={s}"),
        );
        assert_eq!(got[nc..], expected[nc..]);
        for (a, b) in rz.iter().zip(&original) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
    check.report("pfa2rc");
}

#[test]
fn pfa2cr_matches_c() {
    let mut rng = XorShift::new(11);
    let all = lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let (idim, n1, n2) = shape(&mut rng, &all, case, true);
        let s = sign(case);
        let nc = if idim == 1 {
            (n1 / 2 + 1) * n2
        } else {
            n1 * (n2 / 2 + 1)
        };
        let mut cz = complex_input(&mut rng, nc, case);
        let original = cz.clone();
        let mut got = vec![-987.25; n1 * n2 + 3];
        let mut expected = got.clone();
        unsafe {
            c::pfa2cr(
                s,
                idim,
                n1 as i32,
                n2 as i32,
                cz.as_mut_ptr(),
                expected.as_mut_ptr(),
            );
        }
        rust::pfa2cr(s, idim, n1 as i32, n2 as i32, &original, &mut got);
        check.real(
            &got[..n1 * n2],
            &expected[..n1 * n2],
            &format!("pfa2cr {n1}x{n2} idim={idim} sign={s}"),
        );
        assert_eq!(got[n1 * n2..], expected[n1 * n2..]);
        for (a, b) in cz.iter().zip(&original) {
            assert_eq!(a.r.to_bits(), b.r.to_bits());
            assert_eq!(a.i.to_bits(), b.i.to_bits());
        }
    }
    check.report("pfa2cr");
}

// C can also process some undocumented kernel lengths (e.g. 17 is a no-op).
// A modular inverse must exist for every selected radix to keep C's indices
// inside the n-element buffer. Repeated unsupported factors such as 7*7
// have no inverse and cause an out-of-bounds C access, so cannot be oracles.
fn c_kernel_has_bounded_indices(n: usize) -> bool {
    if n == 0 {
        return true;
    }
    let mut nleft = n;
    for factor in [16, 13, 11, 9, 8, 7, 5, 4, 3, 2] {
        if nleft % factor == 0 {
            nleft /= factor;
            if !(1..=factor).any(|j| (j * (n / factor)) % factor == 1) {
                return false;
            }
        }
    }
    true
}

// Odd lengths with admissible floor(n/2), npfa/npfar sizes, and the
// nf=nfft/2+1 lengths used by cepstrum programs, wherever C stays in bounds.
fn compatibility_lengths() -> Vec<usize> {
    let all = lengths();
    let mut out = vec![1];
    for &n in &all {
        out.push(2 * n + 1);
        out.push(unsafe { c::npfar((2 * n) as i32) } as usize);
        for candidate in [unsafe { c::npfa(n as i32) } as usize, n + 1] {
            if c_kernel_has_bounded_indices(candidate / 2) {
                out.push(candidate);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    assert!(out.contains(&85)); // suicepstrum: nfft=168, newnf=85.
    assert!(out.contains(&35)); // npfa(35), with undocumented no-op kernel 17.
    out
}

fn same_bits(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn pfacr_odd_and_program_lengths_match_c() {
    let mut rng = XorShift::new(12);
    let all = compatibility_lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let n = select_length(&mut rng, &all, case);
        let s = sign(case);
        let mut cz = complex_input(&mut rng, n / 2 + 1, case);
        let original = cz.clone();
        let mut got = input(&mut rng, n.max(2) + 3, case + 7);
        let before = got.clone();
        let mut expected = got.clone();
        unsafe { c::pfacr(s, n as i32, cz.as_mut_ptr(), expected.as_mut_ptr()) };
        rust::pfacr(s, &original, &mut got[..n]);
        check.real(
            &got[..n],
            &expected[..n],
            &format!("pfacr compatibility n={n} sign={s}"),
        );
        // At n=1, C writes an unused rz[1] beyond its documented output.
        // Rust exposes only the in-bounds result and must preserve its guard.
        same_bits(&got[n..], &before[n..]);
        same_bits(&expected[n.max(2)..], &before[n.max(2)..]);
        assert_eq!(cz, original);
    }
    check.report("pfacr odd/program");
}

#[test]
fn pfarc_odd_and_program_lengths_match_c() {
    let mut rng = XorShift::new(13);
    let all = compatibility_lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let n = select_length(&mut rng, &all, case);
        let s = sign(case);
        let mut rz = input(&mut rng, n, case);
        let original = rz.clone();
        // C writes float slot n+1 even for odd n. Pass precisely the
        // required output extent to Rust, with guard storage behind it.
        let nc = (n + 3) / 2;
        let mut got = complex_input(&mut rng, nc + 3, case + 7);
        let before = got.clone();
        let mut expected = got.clone();
        unsafe { c::pfarc(s, n as i32, rz.as_mut_ptr(), expected.as_mut_ptr()) };
        rust::pfarc(s, &original, &mut got[..nc]);
        check.complex(
            &got[..nc],
            &expected[..nc],
            &format!("pfarc compatibility n={n} sign={s}"),
        );
        assert_eq!(got[nc..], before[nc..]);
        assert_eq!(expected[nc..], before[nc..]);
        same_bits(&rz, &original);
    }
    check.report("pfarc odd/program");
}

fn compatibility_shape(rng: &mut XorShift, all: &[usize], case: usize) -> (i32, usize, usize) {
    let idim = 1 + (case % 2) as i32;
    let (n, other) = if case < 2 * all.len() {
        (all[case / 2], 1 + (case / 2) % 2)
    } else {
        (
            all[rng.next_u32() as usize % 185],
            rng.range_i32(1, 9) as usize,
        )
    };
    if idim == 1 {
        (idim, n, other)
    } else {
        (idim, other, n)
    }
}

#[test]
fn pfa2cr_odd_and_program_lengths_match_c() {
    let mut rng = XorShift::new(14);
    let all = compatibility_lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let (idim, n1, n2) = compatibility_shape(&mut rng, &all, case);
        let s = sign(case);
        // Actual C float-slot extents, including the odd-length row stride
        // and final unmerge row. No uninitialized or out-of-bounds C reads.
        let (nc, nr) = if idim == 1 {
            (((n1 + 2) * n2) / 2, n1 * n2 + usize::from(n1 == 1))
        } else {
            ((n1 * (n2 + 2)) / 2, n1 * (n2 + n2 % 2))
        };
        let mut cz = complex_input(&mut rng, nc, case);
        let original = cz.clone();
        let mut got = input(&mut rng, nr + 3, case + 7);
        let before = got.clone();
        let mut expected = got.clone();
        unsafe {
            c::pfa2cr(
                s,
                idim,
                n1 as i32,
                n2 as i32,
                cz.as_mut_ptr(),
                expected.as_mut_ptr(),
            );
        }
        rust::pfa2cr(s, idim, n1 as i32, n2 as i32, &original, &mut got[..nr]);
        check.real(
            &got[..nr],
            &expected[..nr],
            &format!("pfa2cr compatibility {n1}x{n2} idim={idim} sign={s}"),
        );
        same_bits(&got[nr..], &before[nr..]);
        same_bits(&expected[nr..], &before[nr..]);
        assert_eq!(cz, original);
    }
    check.report("pfa2cr odd/program");
}

#[test]
fn pfa2rc_odd_and_program_lengths_match_c() {
    let mut rng = XorShift::new(15);
    let all = compatibility_lengths();
    let mut check = Check::default();
    for case in 0..4096 {
        let (idim, n1, n2) = compatibility_shape(&mut rng, &all, case);
        let s = sign(case);
        let slots = if idim == 1 {
            (n1 + 2) * n2
        } else {
            n1 * (n2 + 2)
        };
        let nc = (slots + 1) / 2;
        let mut rz = input(&mut rng, n1 * n2, case);
        let original = rz.clone();
        let mut got = complex_input(&mut rng, nc + 3, case + 7);
        let before = got.clone();
        let mut expected = got.clone();
        unsafe {
            c::pfa2rc(
                s,
                idim,
                n1 as i32,
                n2 as i32,
                rz.as_mut_ptr(),
                expected.as_mut_ptr(),
            );
        }
        rust::pfa2rc(s, idim, n1 as i32, n2 as i32, &original, &mut got[..nc]);
        check.complex(
            &got[..nc],
            &expected[..nc],
            &format!("pfa2rc compatibility {n1}x{n2} idim={idim} sign={s}"),
        );
        assert_eq!(got[nc..], before[nc..]);
        assert_eq!(expected[nc..], before[nc..]);
        same_bits(&rz, &original);
    }
    check.report("pfa2rc odd/program");
}

#[test]
fn pfacr_suicepstrum_168_regression() {
    let nfft = rust::npfar(168) as usize;
    assert_eq!(nfft, 168);
    let nf = nfft / 2 + 1;
    assert_eq!(nf, 85);
    let mut rng = XorShift::new(16);
    let mut check = Check::default();
    for case in 0..2048 {
        let mut cz = complex_input(&mut rng, nf, case);
        // The program allocates and clears nfft output slots, but transforms
        // only newnf slots. The remaining output must stay zero.
        let mut got = vec![0.0; nfft];
        let mut expected = got.clone();
        let s = sign(case);
        unsafe { c::pfacr(s, nf as i32, cz.as_mut_ptr(), expected.as_mut_ptr()) };
        rust::pfacr(s, &cz, &mut got[..nf]);
        check.real(
            &got,
            &expected,
            &format!("suicepstrum nfft=168 newnf=85 sign={s}"),
        );
        same_bits(&got[nf..], &vec![0.0; nfft - nf]);
    }
    check.report("suicepstrum regression");
}
