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

//! Differential tests against par/lib/refaniso.c.
//! Author: CWP: Andreas Rueger, 1994-1996, Colorado School of Mines
#![allow(non_snake_case)]
use oxwave_cref::{self, rel_diff, XorShift};
use oxwave_f_refaniso::{self as rust, Stiff2D};
extern "C" {
    fn graebner2D(
        _: *mut Stiff2D,
        _: f64,
        _: *mut Stiff2D,
        _: f64,
        _: f64,
        _: i32,
        _: i32,
        _: i32,
        _: *mut f64,
    ) -> i32;
    fn gvelpolSH(
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
    );
    fn gvelpolTI(
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
    ) -> i32;
    fn p_hor2DTI(_: *mut Stiff2D, _: f64, _: i32, _: *mut f64) -> i32;
    fn p_vert2DVTI(_: *mut Stiff2D, _: f64, _: i32, _: *mut f64) -> i32;
    fn rottens2D(_: *mut Stiff2D, _: f64);
    fn stiff2thomVTI(
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
    ) -> i32;
    fn stiff2tv(
        _: *mut Stiff2D,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
    ) -> i32;
    fn thom2stiffTI(_: f64, _: f64, _: f64, _: f64, _: f64, _: f64, _: *mut Stiff2D, _: i32)
        -> i32;
    fn thom2tv(
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
        _: *mut f64,
    ) -> i32;
    fn v_phase2DVTI(_: *mut Stiff2D, _: f64, _: i32, _: *mut f64) -> i32;
}

const CASES: usize = 8192;
const SENTINEL: f64 = 123.456789;
const PI: f64 = std::f64::consts::PI;

#[track_caller]
fn compare<const N: usize>(rust: [f64; N], c: [f64; N], case: usize) {
    for (j, (a, b)) in rust.into_iter().zip(c).enumerate() {
        assert!(
            a == b
                || (a.is_nan() && b.is_nan())
                || (a.is_finite() && b.is_finite() && rel_diff(a, b) <= 1e-5),
            "case {case}, output {j}: Rust {a:?} != C {b:?}",
        );
    }
}

// Distinct initial values check outputs left untouched on early returns.
fn five<T>(
    call: impl FnOnce(&mut f64, &mut f64, &mut f64, &mut f64, &mut f64) -> T,
) -> (T, [f64; 5]) {
    let mut out = [SENTINEL, -29.0, 17.0, -41.0, 83.0];
    let [a, b, c, d, e] = &mut out;
    (call(a, b, c, d, e), out)
}

fn tensor(s: Stiff2D) -> [f64; 9] {
    [
        s.a1111, s.a3333, s.a1133, s.a1313, s.a1113, s.a3313, s.a1212, s.a2323, s.a1223,
    ]
}

// Generate shared inputs with the reference rather than the port under test.
fn physical(rng: &mut XorShift, tilted: bool) -> Stiff2D {
    let vp = rng.uniform(1.5, 7.0);
    let vs = vp * rng.uniform(0.15, 0.6);
    let delta = rng.uniform(-0.05, 0.15);
    let eps = delta + rng.uniform(0.01, 0.3);
    let gamma = rng.uniform(-0.2, 0.5);
    let phi = if tilted { rng.uniform(-PI, PI) } else { 0.0 };
    let mut s = Stiff2D::default();
    assert_eq!(
        unsafe { thom2stiffTI(vp, vs, eps, delta, gamma, phi, &mut s, 1) },
        1
    );
    s
}

fn signed_tensor(rng: &mut XorShift) -> Stiff2D {
    let scale = 10.0_f64.powf(rng.uniform(-12.0, 12.0));
    Stiff2D {
        a1111: rng.uniform(-20.0, 20.0) * scale,
        a3333: rng.uniform(-20.0, 20.0) * scale,
        a1133: rng.uniform(-20.0, 20.0) * scale,
        a1313: rng.uniform(-20.0, 20.0) * scale,
        a1113: rng.uniform(-20.0, 20.0) * scale,
        a3313: rng.uniform(-20.0, 20.0) * scale,
        a1212: rng.uniform(-20.0, 20.0) * scale,
        a2323: rng.uniform(-20.0, 20.0) * scale,
        a1223: rng.uniform(-20.0, 20.0) * scale,
    }
}

fn mixed_tensor(rng: &mut XorShift, i: usize) -> Stiff2D {
    let mut s = if i % 2 == 0 {
        physical(rng, false)
    } else {
        signed_tensor(rng)
    };
    match i % 37 {
        0 => s = Stiff2D::default(),
        1 => s.a3333 = 0.0,
        2 => s.a1313 = 0.0,
        3 => s.a3333 = s.a1313,
        4 => s.a2323 = 0.0,
        5 => s.a1111 = -0.0,
        6 => s.a1133 = -s.a1313,
        _ => {}
    }
    s
}

fn angle(rng: &mut XorShift, i: usize) -> f64 {
    // Straddle the FLT_EPSILON clamp and the physical angle range.
    let edge = [
        0.0,
        -0.0,
        1.0,
        -1.0,
        1.0 + f32::EPSILON as f64 / 4.0,
        1.0 + f32::EPSILON as f64,
        1.0 - f32::EPSILON as f64 / 4.0,
        -1.0 - f32::EPSILON as f64 / 4.0,
        1.1,
        -1.1,
    ];
    if i % 31 < edge.len() {
        edge[i % 31]
    } else {
        rng.uniform(-1.0, 1.0)
    }
}

#[test]
fn graebner2D_matches_c() {
    let mut rng = XorShift::new(451);
    let mut successes = [0; 8];
    let mut errors = 0;
    for i in 0..CASES {
        let mut a = physical(&mut rng, false);
        let mut b = physical(&mut rng, false);
        if i % 17 == 0 {
            b = a;
        }
        let mut rho1 = rng.uniform(1.0, 4.0);
        let mut rho2 = rng.uniform(1.0, 4.0);
        let mut pl = rng.uniform(-0.85, 0.85) / a.a1111.max(b.a1111).sqrt();
        let mut modei = (i & 1) as i32;
        let mut modet = ((i >> 1) & 1) as i32;
        let mut rort = ((i >> 2) & 1) as i32;
        // Singular media, evanescence, bad flags and signed densities.
        match i {
            0..=7 => pl = 0.0,
            8..=15 => pl = -0.0,
            16..=23 => pl = 10.0,
            24..=31 => a = Stiff2D::default(),
            32..=39 => b = Stiff2D::default(),
            40..=47 => a.a1111 = -a.a1111,
            48..=55 => rho1 = 0.0,
            56..=63 => rho2 = -rho2,
            64..=71 => {
                rho1 = 0.0;
                rho2 = 0.0;
            }
            72..=79 => a.a3333 = a.a1313,
            80..=87 => modei = -1,
            88..=95 => modet = 2,
            96..=103 => rort = 2,
            _ => {}
        }
        let mut rr = SENTINEL;
        let mut cc = SENTINEL;
        let r = rust::graebner2D(&a, rho1, &b, rho2, pl, modei, modet, rort, &mut rr);
        let c = unsafe { graebner2D(&mut a, rho1, &mut b, rho2, pl, modei, modet, rort, &mut cc) };
        assert_eq!(r, c, "case {i}");
        compare([rr], [cc], i);
        if r == 1 {
            successes[i & 7] += 1;
        } else {
            errors += 1;
        }
    }
    assert!(successes.iter().all(|&n| n > 500), "{successes:?}");
    assert!(errors > 0);
}

#[test]
fn gvelpolSH_matches_c() {
    let mut rng = XorShift::new(452);
    for i in 0..CASES {
        let s = mixed_tensor(&mut rng, i);
        let px = if i % 17 == 0 {
            -0.0
        } else {
            rng.uniform(-2.0, 2.0)
        };
        let pz = if i % 19 == 0 {
            0.0
        } else {
            rng.uniform(-2.0, 2.0)
        };
        let (_, rr) =
            five(|a, b, c, d, e| rust::gvelpolSH(s.a1212, s.a2323, s.a1223, px, pz, a, b, c, d, e));
        let (_, cc) = five(|a, b, c, d, e| unsafe {
            gvelpolSH(s.a1212, s.a2323, s.a1223, px, pz, a, b, c, d, e)
        });
        for j in 0..5 {
            assert_eq!(rr[j].to_bits(), cc[j].to_bits(), "case {i}, output {j}");
        }
    }
}

#[test]
fn gvelpolTI_matches_c() {
    let mut rng = XorShift::new(453);
    let mut successes = 0;
    let mut errors = 0;
    for i in 0..CASES {
        let mut s = if i % 2 == 0 {
            physical(&mut rng, true)
        } else {
            mixed_tensor(&mut rng, i)
        };
        let (mut px, mut pz) = (rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0));
        if i % 2 == 0 {
            // On-shell slowness for both P and SV in tilted media.
            let theta = rng.uniform(-1.5, 1.5);
            let sin = theta.sin();
            assert_eq!(
                unsafe { p_hor2DTI(&mut s, sin, ((i >> 1) & 1) as i32, &mut px) },
                1
            );
            pz = theta.cos() * (px / sin);
        }
        if i % 31 == 0 {
            px = 0.0;
            pz = 0.0;
        }
        if i == 1 {
            // Exact zero polarization denominator.
            s = Stiff2D {
                a1111: 1.0,
                a1313: 1.0,
                ..Stiff2D::default()
            };
            px = 1.0;
            pz = 0.0;
        }
        let (r, rr) = five(|a, b, c, d, e| {
            rust::gvelpolTI(
                s.a1111, s.a3333, s.a1133, s.a1313, s.a1113, s.a3313, px, pz, a, b, c, d, e,
            )
        });
        let (c, cc) = five(|a, b, c, d, e| unsafe {
            gvelpolTI(
                s.a1111, s.a3333, s.a1133, s.a1313, s.a1113, s.a3313, px, pz, a, b, c, d, e,
            )
        });
        assert_eq!(r, c, "case {i}");
        compare(rr, cc, i);
        if r == 1 {
            successes += 1;
        } else {
            errors += 1;
        }
    }
    assert!(successes > 3000 && errors > 100);
}

#[test]
fn p_hor2DTI_matches_c() {
    // Reference err() exits. Child processes let us compare its fatal path.
    if let Ok(child) = std::env::var("OXWAVE_REFANISO_FATAL") {
        let mut s = Stiff2D {
            a1111: -4.0,
            a3333: -4.0,
            a1313: -1.0,
            ..Stiff2D::default()
        };
        let mode = if child.ends_with('0') { 0 } else { 1 };
        let mut out = SENTINEL;
        if child.starts_with("rust") {
            oxwave_core::par::initargs(&["refaniso-diff".to_owned()]);
            rust::p_hor2DTI(&s, 0.0, mode, &mut out);
        } else {
            let mut argv = oxwave_cref::CArgv::new(&["refaniso-diff"]);
            unsafe {
                oxwave_cref::initargs(1, argv.argv());
                p_hor2DTI(&mut s, 0.0, mode, &mut out);
            }
        }
        panic!("expected ray initialization to exit");
    }
    let mut rng = XorShift::new(454);
    let mut successes = 0;
    let mut errors = 0;
    for i in 0..CASES {
        let mut s = physical(&mut rng, true);
        if i % 41 == 0 {
            s = Stiff2D::default();
        }
        let sin = angle(&mut rng, i);
        let mode = if i % 43 == 0 { -1 } else { (i & 1) as i32 };
        let mut rr = SENTINEL;
        let mut cc = SENTINEL;
        let r = rust::p_hor2DTI(&s, sin, mode, &mut rr);
        let c = unsafe { p_hor2DTI(&mut s, sin, mode, &mut cc) };
        assert_eq!(r, c, "case {i}");
        compare([rr], [cc], i);
        if r == 1 {
            successes += 1;
        } else {
            errors += 1;
        }
    }
    assert!(successes > 5000 && errors > 100);
    for mode in 0..=1 {
        let run = |implementation| {
            std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "p_hor2DTI_matches_c", "--nocapture"])
                .env("OXWAVE_REFANISO_FATAL", format!("{implementation}{mode}"))
                .output()
                .unwrap()
        };
        let r = run("rust");
        let c = run("c");
        assert_eq!(r.status.code(), Some(1));
        assert_eq!(c.status.code(), Some(1));
        assert_eq!(r.stderr, c.stderr);
        assert_eq!(
            r.stderr,
            b"\nrefaniso-diff: \n ERROR: ray initialization \n\n"
        );
    }
}

#[test]
fn p_vert2DVTI_matches_c() {
    let mut rng = XorShift::new(455);
    let mut successes = 0;
    let mut errors = 0;
    for i in 0..CASES {
        let mut s = mixed_tensor(&mut rng, i);
        let pl = if i % 17 == 0 {
            0.0
        } else {
            rng.uniform(-1.5, 1.5)
        };
        let mode = [-1, 0, 1, 2][i % 4];
        let mut rr = SENTINEL;
        let mut cc = SENTINEL;
        let r = rust::p_vert2DVTI(&s, pl, mode, &mut rr);
        let c = unsafe { p_vert2DVTI(&mut s, pl, mode, &mut cc) };
        assert_eq!(r, c, "case {i}");
        compare([rr], [cc], i);
        if r == 1 {
            successes += 1;
        } else {
            errors += 1;
        }
    }
    assert!(successes > 1000 && errors > 1000);
}

#[test]
fn rottens2D_matches_c() {
    let mut rng = XorShift::new(456);
    for i in 0..CASES {
        let mut rr = if i % 2 == 0 {
            physical(&mut rng, true)
        } else {
            mixed_tensor(&mut rng, i)
        };
        let mut cc = rr;
        let phi = match i % 17 {
            0 => 0.0,
            1 => -0.0,
            2 => PI / 2.0,
            3 => -PI / 2.0,
            4 => PI,
            5 => 2.0 * PI,
            _ => rng.uniform(-1e6, 1e6),
        };
        rust::rottens2D(&mut rr, phi);
        unsafe {
            rottens2D(&mut cc, phi);
        }
        compare(tensor(rr), tensor(cc), i);
    }
}

#[test]
fn stiff2thomVTI_matches_c() {
    let mut rng = XorShift::new(457);
    let mut successes = 0;
    let mut errors = 0;
    for i in 0..CASES {
        let s = mixed_tensor(&mut rng, i);
        let (r, rr) = five(|a, b, c, d, e| {
            rust::stiff2thomVTI(s.a1111, s.a3333, s.a1133, s.a1313, s.a1212, a, b, c, d, e)
        });
        let (c, cc) = five(|a, b, c, d, e| unsafe {
            stiff2thomVTI(s.a1111, s.a3333, s.a1133, s.a1313, s.a1212, a, b, c, d, e)
        });
        assert_eq!(r, c, "case {i}");
        compare(rr, cc, i);
        if r == 1 {
            successes += 1;
        } else {
            errors += 1;
        }
    }
    assert!(successes > 1000 && errors > 1000);
}

#[test]
fn stiff2tv_matches_c() {
    let mut rng = XorShift::new(458);
    for i in 0..CASES {
        let mut s = mixed_tensor(&mut rng, i);
        let (r, rr) = five(|a, b, c, d, e| rust::stiff2tv(&s, a, b, c, d, e));
        let (c, cc) = five(|a, b, c, d, e| unsafe { stiff2tv(&mut s, a, b, c, d, e) });
        assert_eq!(r, c, "case {i}");
        compare(rr, cc, i);
    }
}

#[test]
fn thom2stiffTI_matches_c() {
    let mut rng = XorShift::new(459);
    let mut successes = 0;
    let mut errors = 0;
    for i in 0..CASES {
        let mut vp = rng.uniform(0.1, 8.0);
        let mut vs = vp * rng.uniform(0.0, 1.2);
        let eps = rng.uniform(-0.6, 0.5);
        let mut delta = rng.uniform(-0.5, 0.5);
        let gamma = rng.uniform(-0.6, 0.5);
        let phi = if i % 17 == 0 {
            0.0
        } else {
            rng.uniform(-PI, PI)
        };
        let sign = [1, -1, 0, 2][i % 4];
        match i % 31 {
            0 => vp = 0.0,
            1 => vp = -vp,
            2 => vs = -vs,
            3 => vs = 0.0,
            4 => vs = vp,
            5 => delta = -2.0,
            6 => {
                vp *= 1000.0;
                vs *= 1000.0;
            }
            _ => {}
        }
        let mut rr = signed_tensor(&mut rng);
        let mut cc = rr;
        let r = rust::thom2stiffTI(vp, vs, eps, delta, gamma, phi, &mut rr, sign);
        let c = unsafe { thom2stiffTI(vp, vs, eps, delta, gamma, phi, &mut cc, sign) };
        assert_eq!(r, c, "case {i}");
        compare(tensor(rr), tensor(cc), i);
        if r == 1 {
            successes += 1;
        } else {
            errors += 1;
        }
    }
    assert!(successes > 1000 && errors > 1000);
}

#[test]
fn thom2tv_matches_c() {
    let mut rng = XorShift::new(460);
    for i in 0..CASES {
        let mut vp = rng.uniform(-8.0, 8.0);
        let mut vs = rng.uniform(-4.0, 4.0);
        let mut eps = rng.uniform(-1.0, 1.0);
        let delta = rng.uniform(-1.0, 1.0);
        let mut gamma = rng.uniform(-1.0, 1.0);
        match i % 23 {
            0 => vp = 0.0,
            1 => vs = 0.0,
            2 => vs = vp,
            3 => vs = -vp,
            4 => eps = -0.5,
            5 => gamma = -0.5,
            6 => {
                vp = 0.0;
                vs = 0.0;
                eps = 0.0;
                gamma = 0.0;
            }
            7 => eps = -0.5 + f64::EPSILON,
            8 => gamma = -0.5 - f64::EPSILON,
            _ => {}
        }
        let (r, rr) = five(|a, b, c, d, e| rust::thom2tv(vp, vs, eps, delta, gamma, a, b, c, d, e));
        let (c, cc) =
            five(|a, b, c, d, e| unsafe { thom2tv(vp, vs, eps, delta, gamma, a, b, c, d, e) });
        assert_eq!(r, c, "case {i}");
        compare(rr, cc, i);
    }
}

#[test]
fn v_phase2DVTI_matches_c() {
    let mut rng = XorShift::new(461);
    let mut successes = 0;
    let mut errors = 0;
    for i in 0..CASES {
        let mut s = mixed_tensor(&mut rng, i);
        let sin = angle(&mut rng, i);
        let mode = [-1, 0, 1, 2][i % 4];
        let mut rr = SENTINEL;
        let mut cc = SENTINEL;
        let r = rust::v_phase2DVTI(&s, sin, mode, &mut rr);
        let c = unsafe { v_phase2DVTI(&mut s, sin, mode, &mut cc) };
        assert_eq!(r, c, "case {i}");
        compare([rr], [cc], i);
        if r == 1 {
            successes += 1;
        } else {
            errors += 1;
        }
    }
    assert!(successes > 5000 && errors > 0);
}

#[test]
fn nonfatal_diagnostics_match_c() {
    const DIAGNOSTIC_CASES: usize = 4096;
    const CHILD_ENV: &str = "OXWAVE_REFANISO_DIAGNOSTICS";
    // Separate processes capture the actual stderr bytes without redirecting
    // file descriptors shared by the other concurrently running tests.
    if let Ok(implementation) = std::env::var(CHILD_ENV) {
        let mut rng = XorShift::new(462);
        for i in 0..DIAGNOSTIC_CASES {
            let mut a = mixed_tensor(&mut rng, i);
            let mut b = mixed_tensor(&mut rng, i + 1);
            let rho1 = rng.uniform(-4.0, 4.0);
            let rho2 = rng.uniform(-4.0, 4.0);
            let pl = if i % 17 == 0 {
                0.0
            } else {
                rng.uniform(-2.0, 2.0)
            };
            let sin = angle(&mut rng, i);
            let modei = (i & 1) as i32;
            let modet = ((i >> 1) & 1) as i32;
            let rort = ((i >> 2) & 1) as i32;
            let (mut coeff, mut vertical, mut phase) = (SENTINEL, SENTINEL, SENTINEL);
            let codes = if implementation == "rust" {
                [
                    rust::graebner2D(&a, rho1, &b, rho2, pl, modei, modet, rort, &mut coeff),
                    rust::p_vert2DVTI(&a, pl, modei, &mut vertical),
                    rust::v_phase2DVTI(&a, sin, modei, &mut phase),
                ]
            } else {
                assert_eq!(implementation, "c");
                unsafe {
                    [
                        graebner2D(
                            &mut a, rho1, &mut b, rho2, pl, modei, modet, rort, &mut coeff,
                        ),
                        p_vert2DVTI(&mut a, pl, modei, &mut vertical),
                        v_phase2DVTI(&mut a, sin, modei, &mut phase),
                    ]
                }
            };
            // Round-trip float formatting lets the parent compare outputs and
            // return codes as well as diagnostics for every generated case.
            println!(
                "\ncase {i} {} {} {} {coeff:?} {vertical:?} {phase:?}",
                codes[0], codes[1], codes[2]
            );
        }
        return;
    }
    let run = |implementation| {
        std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "nonfatal_diagnostics_match_c", "--nocapture"])
            .env(CHILD_ENV, implementation)
            .output()
            .unwrap()
    };
    let r = run("rust");
    let c = run("c");
    assert!(r.status.success(), "Rust child failed: {:?}", r.status);
    assert!(c.status.success(), "C child failed: {:?}", c.status);
    let r_stderr = String::from_utf8(r.stderr).unwrap();
    let c_stderr = String::from_utf8(c.stderr).unwrap();
    assert_eq!(r_stderr, c_stderr);
    for diagnostic in [
        "ERROR in slowness medium1 ",
        "ERROR slowness medium 2 ",
        " evanescent energy generated  ",
        "evanescent energy generated  ",
        " error in v_phase ",
    ] {
        assert!(
            c_stderr.lines().any(|line| line == diagnostic),
            "missing {diagnostic:?}"
        );
    }
    let r_stdout = String::from_utf8(r.stdout).unwrap();
    let c_stdout = String::from_utf8(c.stdout).unwrap();
    let cases = |s: &str| {
        s.lines()
            .filter(|line| line.starts_with("case "))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    let r_cases = cases(&r_stdout);
    let c_cases = cases(&c_stdout);
    assert_eq!(r_cases.len(), DIAGNOSTIC_CASES);
    assert_eq!(c_cases.len(), DIAGNOSTIC_CASES);
    for (i, (r, c)) in r_cases.iter().zip(&c_cases).enumerate() {
        let r: Vec<_> = r.split_whitespace().collect();
        let c: Vec<_> = c.split_whitespace().collect();
        assert_eq!(r.len(), 8);
        assert_eq!(c.len(), 8);
        assert_eq!(r[..5], c[..5], "case {i}: return codes");
        compare::<3>(
            std::array::from_fn(|j| r[5 + j].parse().unwrap()),
            std::array::from_fn(|j| c[5 + j].parse().unwrap()),
            i,
        );
    }
}
