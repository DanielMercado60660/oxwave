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

//! Differential tests for `par/lib/upweik.c`.
//! Authors: CWP: Zhenuye Liu, based on code by Dave Hale, 1992.
use oxwave_cref::XorShift;
use oxwave_f_upweik as rust;
extern "C" {
    fn eikpex(
        n: i32,
        da: f32,
        r: f32,
        dr: f32,
        sc: *const f32,
        uc: *const f32,
        wc: *const f32,
        tc: *const f32,
        sn: *const f32,
        un: *mut f32,
        wn: *mut f32,
        tn: *mut f32,
    );
    fn ray_theoretic_sigma(
        n: i32,
        da: f32,
        r: f32,
        dr: f32,
        uc: *const f32,
        wc: *const f32,
        sc: *const f32,
        un: *const f32,
        wn: *const f32,
        sn: *mut f32,
    );
    fn ray_theoretic_beta(
        n: i32,
        da: f32,
        r: f32,
        dr: f32,
        uc: *const f32,
        wc: *const f32,
        bc: *const f32,
        un: *const f32,
        wn: *const f32,
        bn: *mut f32,
    );
    fn eiktam(
        xs: f32,
        zs: f32,
        nz: i32,
        dz: f32,
        fz: f32,
        nx: i32,
        dx: f32,
        fx: f32,
        vel: *mut *mut f32,
        time: *mut *mut f32,
        angle: *mut *mut f32,
        sig: *mut *mut f32,
        bet: *mut *mut f32,
    );
}
fn rand(r: &mut XorShift) -> f32 {
    (r.next_u32() as f64 / u32::MAX as f64) as f32
}
fn check(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    let scale = b
        .iter()
        .filter(|y| y.is_finite())
        .fold(0f32, |x, y| x.max(y.abs()))
        .max(1e-30);
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        if y.is_nan() {
            assert!(x.is_nan(), "index {i}: {x} != NaN");
            continue;
        }
        if y.is_infinite() {
            assert_eq!(x, y, "index {i}: infinity sign/classification differs");
            continue;
        }
        assert!(
            x.is_finite() && (x - y).abs() <= 1e-5 * scale,
            "index {i}: {x} != {y}, scale {scale}"
        );
    }
}
#[test]
fn eikpex_matches() {
    let mut g = XorShift::new(3);
    for k in 0..4000 {
        let n = if k % 37 == 0 {
            1
        } else {
            [2, 3, 7, 16, 31, 64, 127, 256, 1021, 1024][k % 10]
        };
        let da = if k % 3 == 0 {
            (2.0 * std::f64::consts::PI / n as f64) as f32
        } else {
            0.05 + rand(&mut g)
        };
        let r = [0.01, 1., 100., 1000.][(k / 10) % 4] * (1. + rand(&mut g));
        let dr = if k % 37 == 0 {
            -1.
        } else if k % 19 == 0 {
            0.
        } else {
            r * (0.01 + rand(&mut g) * 0.5)
        };
        let speed = [0.0002, 0.01, 1., 10.][(k / 40) % 4] * (1. + rand(&mut g));
        let sc: Vec<_> = (0..n)
            .map(|i| speed * (1. + 0.001 * (2. * std::f32::consts::PI * i as f32 / n as f32).cos()))
            .collect();
        let sn: Vec<_> = sc
            .iter()
            .map(|&s| s * (1. + 0.001 * rand(&mut g)))
            .collect();
        let wc: Vec<_> = (0..n)
            .map(|i| {
                r * speed * 0.2 * (2. * std::f32::consts::PI * i as f32 / n as f32 + 0.3).sin()
            })
            .collect();
        let uc: Vec<_> = (0..n)
            .map(|i| {
                if k % 13 == 0 {
                    0.
                } else {
                    ((sc[i] * sc[i] - (wc[i] / r) * (wc[i] / r)) as f64).sqrt() as f32
                }
            })
            .collect();
        let tc: Vec<_> = (0..n).map(|_| r * speed * (rand(&mut g) - 0.5)).collect();
        let (mut u, mut w, mut t) = (vec![0.; n], vec![0.; n], vec![0.; n]);
        let (mut cu, mut cw, mut ct) = (u.clone(), w.clone(), t.clone());
        unsafe {
            eikpex(
                n as i32,
                da,
                r,
                dr,
                sc.as_ptr(),
                uc.as_ptr(),
                wc.as_ptr(),
                tc.as_ptr(),
                sn.as_ptr(),
                cu.as_mut_ptr(),
                cw.as_mut_ptr(),
                ct.as_mut_ptr(),
            );
        }
        rust::eikpex(da, r, dr, &sc, &uc, &wc, &tc, &sn, &mut u, &mut w, &mut t);
        check(&u, &cu);
        check(&w, &cw);
        check(&t, &ct);
        if dr <= 0. {
            // This path only copies: require bit identity, including signed zeros.
            for (actual, expected) in [(&u, &uc), (&w, &wc), (&t, &tc)] {
                assert_eq!(
                    actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                );
            }
        }
    }
}
fn transport(beta: bool) {
    let mut g = XorShift::new(57);
    // C allocates na-2 entries, then tripp accesses index -1 when na == 3;
    // na < 3 also accesses outside the workspaces. Four is the C minimum
    // with defined memory accesses; sizes one and two are covered by eikpex.
    for k in 0..4000 {
        let n = [4, 5, 7, 16, 31, 64, 127, 256, 1021, 1024][k % 10];
        let da = (0.05 + rand(&mut g)) * if k % 7 == 0 { -1. } else { 1. };
        let r = [0.01, 1., 100.][(k / 10) % 3] * (1. + rand(&mut g));
        let dr = r * (0.1 + rand(&mut g)) * if k % 9 == 0 { -0.2 } else { 1. };
        let magnitude = [0.0002, 1., 10.][(k / 30) % 3];
        let sign = if k % 5 == 0 { -1. } else { 1. };
        let uc: Vec<_> = (0..n)
            .map(|_| sign * magnitude * (1. + rand(&mut g)))
            .collect();
        let un: Vec<_> = (0..n)
            .map(|_| sign * magnitude * (1. + rand(&mut g)))
            .collect();
        let wc: Vec<_> = (0..n)
            .map(|_| {
                if k % 11 == 0 {
                    0.
                } else {
                    (rand(&mut g) - 0.5) * magnitude * r * r * da / dr.abs()
                }
            })
            .collect();
        let wn: Vec<_> = wc.iter().map(|w| w * (0.5 + rand(&mut g))).collect();
        let val: Vec<_> = (0..n)
            .map(|_| {
                if k % 17 == 0 {
                    0.
                } else {
                    (rand(&mut g) - 0.5) * [0.001, 1., 1000.][(k / 90) % 3]
                }
            })
            .collect();
        let mut out = vec![0.; n];
        let mut c = out.clone();
        unsafe {
            if beta {
                ray_theoretic_beta(
                    n as i32,
                    da,
                    r,
                    dr,
                    uc.as_ptr(),
                    wc.as_ptr(),
                    val.as_ptr(),
                    un.as_ptr(),
                    wn.as_ptr(),
                    c.as_mut_ptr(),
                )
            } else {
                ray_theoretic_sigma(
                    n as i32,
                    da,
                    r,
                    dr,
                    uc.as_ptr(),
                    wc.as_ptr(),
                    val.as_ptr(),
                    un.as_ptr(),
                    wn.as_ptr(),
                    c.as_mut_ptr(),
                )
            }
        }
        if beta {
            rust::ray_theoretic_beta(da, r, dr, &uc, &wc, &val, &un, &wn, &mut out)
        } else {
            rust::ray_theoretic_sigma(da, r, dr, &uc, &wc, &val, &un, &wn, &mut out)
        }
        check(&out, &c);
    }
}
#[test]
fn sigma_matches() {
    transport(false)
}
#[test]
fn beta_matches() {
    transport(true)
}
fn ptrs(v: &mut [Vec<f32>]) -> Vec<*mut f32> {
    v.iter_mut().map(|x| x.as_mut_ptr()).collect()
}
#[test]
fn eiktam_matches() {
    let mut g = XorShift::new(987);
    for k in 0..3012 {
        let nx = if k >= 3000 {
            127
        } else {
            [2, 3, 5, 8, 17, 32][k % 6]
        };
        let nz = if k >= 3000 {
            128
        } else {
            [2, 3, 7, 8, 16, 31][(k / 6) % 6]
        };
        let dx = [0.1, 1., 10., 25., 100.][(k / 36) % 5];
        let dz = dx * (0.5 + 1.5 * rand(&mut g));
        let fx = -100. + 200. * rand(&mut g);
        let fz = -100. + 200. * rand(&mut g);
        let xf = match k % 4 {
            0 | 2 => 0.,
            1 => 1.,
            _ => rand(&mut g),
        };
        let zf = match k % 4 {
            0 | 1 => 0.,
            2 => 1.,
            _ => rand(&mut g),
        };
        let xs = fx + (nx - 1) as f32 * dx * xf;
        let zs = fz + (nz - 1) as f32 * dz * zf;
        let mut vel = vec![vec![0.; nz]; nx];
        let speed = [1., 1500., 6000.][(k / 180) % 3] * (1. + rand(&mut g));
        for (ix, row) in vel.iter_mut().enumerate() {
            for (iz, x) in row.iter_mut().enumerate() {
                *x = if k % 1000 == 0 {
                    // Division by zero is defined floating point behavior in C;
                    // require the same infinities/NaNs, without a tolerance escape.
                    0.
                } else if k % 11 == 0 {
                    speed
                } else {
                    speed
                        * (1.
                            + 0.03 * ix as f32 / nx as f32
                            + 0.02 * iz as f32 / nz as f32
                            + 0.0001 * (rand(&mut g) - 0.5))
                };
            }
        }
        let mut c = vec![vec![vec![0.; nz]; nx]; 4];
        let mut vp = ptrs(&mut vel);
        let mut cp: Vec<_> = c.iter_mut().map(|x| ptrs(x)).collect();
        unsafe {
            eiktam(
                xs,
                zs,
                nz as i32,
                dz,
                fz,
                nx as i32,
                dx,
                fx,
                vp.as_mut_ptr(),
                cp[0].as_mut_ptr(),
                cp[1].as_mut_ptr(),
                cp[2].as_mut_ptr(),
                cp[3].as_mut_ptr(),
            );
        }
        let mut o = vec![vec![vec![0.; nz]; nx]; 4];
        let [t, a, s, b]: &mut [_; 4] = o.as_mut_slice().try_into().unwrap();
        rust::eiktam(
            xs,
            zs,
            dz,
            fz,
            dx,
            fx,
            &vel.iter().map(|x| x.as_slice()).collect::<Vec<_>>(),
            &mut t.iter_mut().map(|x| x.as_mut_slice()).collect::<Vec<_>>(),
            &mut a.iter_mut().map(|x| x.as_mut_slice()).collect::<Vec<_>>(),
            &mut s.iter_mut().map(|x| x.as_mut_slice()).collect::<Vec<_>>(),
            &mut b.iter_mut().map(|x| x.as_mut_slice()).collect::<Vec<_>>(),
        );
        for j in 0..4 {
            for i in 0..nx {
                check(&o[j][i], &c[j][i]);
            }
        }
    }
}

#[test]
fn eikpex_nonpositive_step_copies_bits() {
    let mut g = XorShift::new(123);
    for k in 0..2000 {
        let n = [1, 2, 7, 32, 257, 1024][k % 6];
        let input: Vec<_> = (0..n)
            .map(|i| match i % 5 {
                0 => 0.,
                1 => -0.,
                2 => -rand(&mut g),
                3 => f32::INFINITY,
                _ => f32::from_bits(0x7fc01234),
            })
            .collect();
        let dr = if k % 2 == 0 { 0. } else { -rand(&mut g) };
        let mut c = [vec![1.; n], vec![1.; n], vec![1.; n]];
        let [cu, cw, ct] = &mut c;
        unsafe {
            eikpex(
                n as i32,
                0.,
                0.,
                dr,
                input.as_ptr(),
                input.as_ptr(),
                input.as_ptr(),
                input.as_ptr(),
                input.as_ptr(),
                cu.as_mut_ptr(),
                cw.as_mut_ptr(),
                ct.as_mut_ptr(),
            );
        }
        let mut o = [vec![1.; n], vec![1.; n], vec![1.; n]];
        let [u, w, t] = &mut o;
        rust::eikpex(0., 0., dr, &input, &input, &input, &input, &input, u, w, t);
        for j in 0..3 {
            for i in 0..n {
                assert_eq!(o[j][i].to_bits(), c[j][i].to_bits());
                assert_eq!(o[j][i].to_bits(), input[i].to_bits());
            }
        }
    }
}

#[test]
fn eikpex_signed_constant_fields_match_bits() {
    let mut g = XorShift::new(0x51_6e_ed);
    for k in 0..2048 {
        let n = [2, 3, 7, 16, 31, 64, 127, 256, 1021, 1024][k % 10];
        let da = if k % 2 == 0 {
            (2.0 * std::f64::consts::PI / n as f64) as f32
        } else {
            0.05 + rand(&mut g)
        };
        let r = [0.0001, 0.01, 1., 100., 10000.][(k / 10) % 5] * (1. + rand(&mut g));
        let dr = r * (0.01 + 0.5 * rand(&mut g));
        let magnitude = [0.0001, 0.01, 1., 100., 10000.][(k / 50) % 5] * (1. + rand(&mut g));
        let sc = vec![if k % 3 == 0 { -magnitude } else { magnitude }; n];
        let next = magnitude * (0.5 + rand(&mut g));
        let sn = vec![if k % 5 == 0 { -next } else { next }; n];
        let uc = vec![if k % 7 == 0 { magnitude } else { -magnitude }; n];
        let wc = vec![if k % 2 == 0 { -0.0 } else { 0.0 }; n];
        let tc: Vec<_> = (0..n)
            .map(|_| r * magnitude * (rand(&mut g) - 0.5))
            .collect();
        // Constant angular fields give a single CFL step without curvature
        // failure, even for negative slowness or a negative initial dt/dr.
        // Only sqrt and the original scalar arithmetic remain: require bits.
        let mut c = [vec![f32::NAN; n], vec![f32::NAN; n], vec![f32::NAN; n]];
        let [cu, cw, ct] = &mut c;
        unsafe {
            eikpex(
                n as i32,
                da,
                r,
                dr,
                sc.as_ptr(),
                uc.as_ptr(),
                wc.as_ptr(),
                tc.as_ptr(),
                sn.as_ptr(),
                cu.as_mut_ptr(),
                cw.as_mut_ptr(),
                ct.as_mut_ptr(),
            );
        }
        let mut o = [vec![f32::NAN; n], vec![f32::NAN; n], vec![f32::NAN; n]];
        let [u, w, t] = &mut o;
        rust::eikpex(da, r, dr, &sc, &uc, &wc, &tc, &sn, u, w, t);
        for j in 0..3 {
            for i in 0..n {
                assert!(c[j][i].is_finite(), "case {k}, output {j}, index {i}");
                assert_eq!(
                    o[j][i].to_bits(),
                    c[j][i].to_bits(),
                    "case {k}, output {j}, index {i}: {} != {}",
                    o[j][i],
                    c[j][i]
                );
            }
        }
    }
}

#[test]
fn eikpex_curvature_error_matches() {
    const CHILD: &str = "OXWAVE_UPWEIK_CURVATURE_CHILD";
    if let Ok(implementation) = std::env::var(CHILD) {
        // Both paths exit, so run each in a child of this same test binary.
        let mut argv = oxwave_cref::CArgv::new(&["upweik-diff"]);
        unsafe {
            oxwave_cref::initargs(argv.argc(), argv.argv());
        }
        oxwave_core::par::initargs(&["upweik-diff".into()]);
        let (mut u, mut w, mut t) = ([0.; 4], [0.; 4], [0.; 4]);
        // At the left boundary a negative angular derivative is retained;
        // with zero slowness it forces uu < 0 on the first intermediate step.
        if implementation == "c" {
            unsafe {
                eikpex(
                    4,
                    0.1,
                    1.,
                    0.1,
                    [0.; 4].as_ptr(),
                    [1.; 4].as_ptr(),
                    [-1.; 4].as_ptr(),
                    [0.; 4].as_ptr(),
                    [0.; 4].as_ptr(),
                    u.as_mut_ptr(),
                    w.as_mut_ptr(),
                    t.as_mut_ptr(),
                );
            }
        } else {
            rust::eikpex(
                0.1, 1., 0.1, &[0.; 4], &[1.; 4], &[-1.; 4], &[0.; 4], &[0.; 4], &mut u, &mut w,
                &mut t,
            );
        }
        panic!("curvature error should exit");
    }
    let run = |implementation| {
        std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "eikpex_curvature_error_matches", "--nocapture"])
            .env(CHILD, implementation)
            .output()
            .unwrap()
    };
    let c = run("c");
    let rust = run("rust");
    assert_eq!(c.status.code(), Some(1));
    assert_eq!(rust.status.code(), c.status.code());
    assert_eq!(rust.stdout, c.stdout);
    assert_eq!(rust.stderr, c.stderr);
    assert_eq!(c.stderr, b"\nupweik-diff: \tRaypath has a too large curvature!\n\t A smoother velocity is required. \n\n");
}
