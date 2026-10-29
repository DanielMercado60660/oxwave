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

//! Differential tests for par/lib/taup.c, by John Anderson and Gabriel Alvarez.
use oxwave_cref::{self, XorShift};
use oxwave_cwp::complex::Complex;
use oxwave_f_taup as rust;

extern "C" {
    fn forward_p_transform(
        nx: i32,
        nt: i32,
        dt: f32,
        pmax: f32,
        pmin: f32,
        dp: f32,
        depthref: f32,
        f1: f32,
        f2: f32,
        freq1: f32,
        freq2: f32,
        lagc: i32,
        lent: i32,
        lenx: i32,
        xopt: i32,
        ninterp: i32,
        nwin: i32,
        prewhite: f32,
        interoff: f32,
        offref: f32,
        igopt: i32,
        dx: f32,
        fx: f32,
        pmula: f32,
        pmulb: f32,
        input: *const *const f32,
        output: *mut *mut f32,
    );
    fn ga_xinterpolate(
        input: *const *const f32,
        output: *mut *mut f32,
        ninterp: i32,
        nt: i32,
        nx: i32,
        freq1: f32,
        freq2: f32,
        lagc: i32,
        lent: i32,
        lenx: i32,
        xopt: i32,
        dt: f32,
        iopt: i32,
    );
    fn fwd_FK_sstack(
        dt: f32,
        nt: i32,
        nx: i32,
        xmin: f32,
        dx: f32,
        np: i32,
        pmin: f32,
        dp: f32,
        fmin: f32,
        input: *const *const f32,
        output: *mut *mut f32,
    );
    fn inv_FK_sstack(
        dt: f32,
        nt: i32,
        nx: i32,
        xmin: f32,
        dx: f32,
        np: i32,
        pmin: f32,
        dp: f32,
        fmin: f32,
        input: *const *const f32,
        output: *mut *mut f32,
    );
    fn inverse_p_transform(
        nx: i32,
        nt: i32,
        dt: f32,
        pmax: f32,
        pmin: f32,
        dp: f32,
        depthref: f32,
        f1: f32,
        f2: f32,
        interoff: f32,
        offref: f32,
        igopt: i32,
        dx: f32,
        fx: f32,
        input: *const *const f32,
        output: *mut *mut f32,
    );
    fn gofx(opt: i32, x: f32, off: f32, depth: f32) -> f32;
    fn freqweight(j: i32, df: f32, f1: f32, f2: f32) -> f32;
    fn compute_r(w: f32, nx: i32, g: *const f32, np: i32, dp: f32, r: *mut Complex);
    fn compute_rhs(
        w: f32,
        nx: i32,
        g: *const f32,
        d: *const Complex,
        np: i32,
        p: f32,
        dp: f32,
        r: *mut Complex,
    );
    fn rcdot(n: i32, a: *const Complex, b: *const Complex) -> f32;
    fn htmul(n: i32, a: *const Complex, x: *const Complex, y: *mut Complex);
    fn runav(n: i32, len: i32, a: *const f32, b: *mut f32);
    fn taupmute(ip: i32, ipa: i32, ipb: i32, nt: i32, ltap: i32, rt: *mut f32);
    fn rho_filter(npoints: i32, nt: i32, dt: f32, rho: *mut f32);
    fn fwd_tx_sstack(
        dt: f32,
        nt: i32,
        nx: i32,
        xmin: f32,
        dx: f32,
        np: i32,
        pmin: f32,
        dp: f32,
        tr: *const *const f32,
        out: *mut *mut f32,
    );
    fn inv_tx_sstack(
        dt: f32,
        nt: i32,
        nx: i32,
        npoints: i32,
        xmin: f32,
        dx: f32,
        np: i32,
        pmin: f32,
        dp: f32,
        tr: *const *const f32,
        out: *mut *mut f32,
    );
    fn ctoephcg(
        it: i32,
        n: i32,
        a: *const Complex,
        x: *mut Complex,
        y: *const Complex,
        s: *mut Complex,
        ss: *mut Complex,
        g: *mut Complex,
        rr: *mut Complex,
    ) -> i32;
    fn ctoep(
        n: i32,
        r: *const Complex,
        a: *mut Complex,
        b: *const Complex,
        f: *mut Complex,
        g: *mut Complex,
    ) -> i32;
}

#[test]
fn test_forward_p_transform() {
    let mut rng = XorShift::new(534);
    for k in 0..2000 {
        let nt = [7, 16, 31, 64][k % 4];
        let nx = k % 4 + 2;
        let ninterp = k % 2;
        let np = 3;
        let input: Vec<Vec<f32>> = (0..nx)
            .map(|_| (0..nt).map(|_| val(&mut rng)).collect())
            .collect();
        let mut a = vec![vec![0.0; nt]; np];
        let mut b = a.clone();
        let inp: Vec<_> = input.iter().map(|x| x.as_ptr()).collect();
        let mut out: Vec<_> = b.iter_mut().map(|x| x.as_mut_ptr()).collect();
        rust::forward_p_transform(
            nt, 0.004, 2.0, 0.0, 1.0, 50.0, 30.0, 100.0, 3.0, 20.0, 4, 3, 1, 1, ninterp, 1, 0.1,
            0.0, 100.0, 3, 10.0, -20.0, 3.0, 4.0, &input, &mut a,
        );
        unsafe {
            forward_p_transform(
                nx as i32,
                nt as i32,
                0.004,
                2.0,
                0.0,
                1.0,
                50.0,
                30.0,
                100.0,
                3.0,
                20.0,
                4,
                3,
                1,
                1,
                ninterp as i32,
                1,
                0.1,
                0.0,
                100.0,
                3,
                10.0,
                -20.0,
                3.0,
                4.0,
                inp.as_ptr(),
                out.as_mut_ptr(),
            );
        }
        let m = b.iter().flatten().fold(1e-20f32, |v, x| v.max(x.abs()));
        for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
            assert!((x - y).abs() <= 1e-5 * m, "case {k}: {x} != {y}, max {m}");
        }
    }
}

#[test]
fn test_ga_xinterpolate() {
    let mut rng = XorShift::new(533);
    for k in 0..2000 {
        let nt = [7, 16, 31, 64][k % 4];
        let nx = [1, 2, 3, 7, 16][k % 5];
        let ninterp = (k / 3) % 3;
        let iopt = (k % 3) as i32;
        let nout = if iopt == 0 {
            (nx - 1) * (ninterp + 1) + 1
        } else {
            nx
        };
        let input: Vec<Vec<f32>> = (0..nx)
            .map(|_| (0..nt).map(|_| val(&mut rng)).collect())
            .collect();
        let mut a = vec![vec![0.0; nt]; nout];
        let mut b = a.clone();
        let inp: Vec<_> = input.iter().map(|x| x.as_ptr()).collect();
        let mut out: Vec<_> = b.iter_mut().map(|x| x.as_mut_ptr()).collect();
        let xopt = (k % 2) as i32;
        let lenx = if k % 2 == 0 { 1 } else { 3 };
        rust::ga_xinterpolate(
            &input, &mut a, ninterp, nt, 3.0, 20.0, 4, 3, lenx, xopt, 0.004, iopt,
        );
        unsafe {
            ga_xinterpolate(
                inp.as_ptr(),
                out.as_mut_ptr(),
                ninterp as i32,
                nt as i32,
                nx as i32,
                3.0,
                20.0,
                4,
                3,
                lenx,
                xopt,
                0.004,
                iopt,
            );
        }
        let m = b.iter().flatten().fold(1e-20f32, |v, x| v.max(x.abs()));
        for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
            assert!((x - y).abs() <= 1e-5 * m, "case {k}: {x} != {y}, max {m}");
        }
    }
}

#[test]
fn test_fwd_fk_sstack() {
    let mut rng = XorShift::new(532);
    for k in 0..2000 {
        let nt = [2, 3, 7, 16, 31, 64][k % 6];
        let np = k % 5 + 1;
        let nx = k % 7 + 1;
        let input: Vec<Vec<f32>> = (0..nx)
            .map(|_| (0..nt).map(|_| val(&mut rng)).collect())
            .collect();
        let mut a = vec![vec![0.0; nt]; np];
        let mut b = a.clone();
        let inp: Vec<_> = input.iter().map(|x| x.as_ptr()).collect();
        let mut out: Vec<_> = b.iter_mut().map(|x| x.as_mut_ptr()).collect();
        rust::fwd_FK_sstack(
            0.004, nt, -10.0, 10.0, -0.00001, 0.00001, 0.0, &input, &mut a,
        );
        unsafe {
            fwd_FK_sstack(
                0.004,
                nt as i32,
                nx as i32,
                -10.0,
                10.0,
                np as i32,
                -0.00001,
                0.00001,
                0.0,
                inp.as_ptr(),
                out.as_mut_ptr(),
            );
        }
        let m = b.iter().flatten().fold(1e-20f32, |v, x| v.max(x.abs()));
        for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
            assert!((x - y).abs() <= 1e-5 * m, "{x} != {y}");
        }
    }
}

#[test]
fn test_inv_fk_sstack() {
    let mut rng = XorShift::new(531);
    for k in 0..2000 {
        let nt = [2, 3, 7, 16, 31, 64][k % 6];
        let np = k % 5 + 1;
        let nx = k % (4 * np) + 1;
        let input: Vec<Vec<f32>> = (0..np)
            .map(|_| (0..nt).map(|_| val(&mut rng)).collect())
            .collect();
        let mut a = vec![vec![0.0; nt]; nx];
        let mut b = a.clone();
        let inp: Vec<_> = input.iter().map(|x| x.as_ptr()).collect();
        let mut out: Vec<_> = b.iter_mut().map(|x| x.as_mut_ptr()).collect();
        rust::inv_FK_sstack(0.004, nt, -10.0, 10.0, 0.0, 0.00001, 0.0, &input, &mut a);
        unsafe {
            inv_FK_sstack(
                0.004,
                nt as i32,
                nx as i32,
                -10.0,
                10.0,
                np as i32,
                0.0,
                0.00001,
                0.0,
                inp.as_ptr(),
                out.as_mut_ptr(),
            );
        }
        let m = b.iter().flatten().fold(1e-20f32, |v, x| v.max(x.abs()));
        for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
            assert!((x - y).abs() <= 1e-5 * m, "{x} != {y}");
        }
    }
}

#[test]
fn test_inverse_p_transform() {
    let mut rng = XorShift::new(12345);
    for k in 0..2000 {
        let nt = [2, 3, 7, 16, 31, 64][k % 6];
        let nx = k % 5 + 1;
        let np = 3;
        let nfft = oxwave_pfa_fft::npfar(nt as i32) as usize;
        let input: Vec<Vec<f32>> = (0..np)
            .map(|_| {
                let mut a = vec![0.0; nfft + 2];
                for x in &mut a[..nt] {
                    *x = val(&mut rng);
                }
                a
            })
            .collect();
        let mut a = vec![vec![0.0; nt]; nx];
        let mut b = a.clone();
        let inp: Vec<_> = input.iter().map(|x| x.as_ptr()).collect();
        let mut out: Vec<_> = b.iter_mut().map(|x| x.as_mut_ptr()).collect();
        rust::inverse_p_transform(
            nt, 0.004, 2.0, 0.0, 1.0, 50.0, 30.0, 100.0, 0.0, 100.0, 3, 10.0, -20.0, &input, &mut a,
        );
        unsafe {
            inverse_p_transform(
                nx as i32,
                nt as i32,
                0.004,
                2.0,
                0.0,
                1.0,
                50.0,
                30.0,
                100.0,
                0.0,
                100.0,
                3,
                10.0,
                -20.0,
                inp.as_ptr(),
                out.as_mut_ptr(),
            );
        }
        let m = b.iter().flatten().fold(1e-20f32, |v, x| v.max(x.abs()));
        for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
            assert!((x - y).abs() <= 1e-5 * m, "{x} != {y}");
        }
    }
}
fn val(r: &mut XorShift) -> f32 {
    (r.next_u32() as i32 as f32) / (i32::MAX as f32)
}
fn close(a: f32, b: f32) {
    assert!(
        (a - b).abs() <= 1e-5 * a.abs().max(b.abs()).max(1e-20),
        "{a} != {b}"
    );
}
fn cv(a: &[Complex], b: &[Complex]) {
    let m = b
        .iter()
        .map(|x| x.r.abs().max(x.i.abs()))
        .fold(1e-20f32, f32::max);
    for (x, y) in a.iter().zip(b) {
        assert!(
            (x.r - y.r).abs() <= 1e-5 * m && (x.i - y.i).abs() <= 1e-5 * m,
            "{x:?} != {y:?}"
        );
    }
}
#[test]
fn test_gofx() {
    let mut r = XorShift::new(1);
    for i in 0..5000 {
        let (x, o, d) = (val(&mut r) * 100., val(&mut r) * 30., val(&mut r) * 100.);
        close(rust::gofx(i % 7, x, o, d), unsafe { gofx(i % 7, x, o, d) });
    }
}
#[test]
fn test_freqweight() {
    let mut r = XorShift::new(2);
    for i in 0..5000 {
        let (df, a, b) = (val(&mut r), val(&mut r), val(&mut r));
        close(rust::freqweight(i % 30 - 15, df, a, b), unsafe {
            freqweight(i % 30 - 15, df, a, b)
        });
    }
}
#[test]
fn test_compute_r() {
    let mut r = XorShift::new(3);
    for k in 0..2000 {
        let n = [1, 2, 3, 7, 16, 31][k % 6];
        let g: Vec<_> = (0..n).map(|_| val(&mut r)).collect();
        let mut a = vec![Complex::default(); n];
        let mut b = a.clone();
        let w = val(&mut r) * 10.;
        let dp = val(&mut r);
        rust::compute_r(w, &g, dp, &mut a);
        unsafe { compute_r(w, n as i32, g.as_ptr(), n as i32, dp, b.as_mut_ptr()) };
        cv(&a, &b);
    }
}
#[test]
fn test_compute_rhs() {
    let mut r = XorShift::new(4);
    for k in 0..2000 {
        let n = [1, 2, 3, 7, 16, 31][k % 6];
        let g: Vec<_> = (0..n).map(|_| val(&mut r)).collect();
        let d: Vec<_> = (0..n)
            .map(|_| Complex {
                r: val(&mut r),
                i: val(&mut r),
            })
            .collect();
        let mut a = vec![Complex::default(); n];
        let mut b = a.clone();
        let (w, p, dp) = (val(&mut r) * 10., val(&mut r), val(&mut r));
        rust::compute_rhs(w, &g, &d, p, dp, &mut a);
        unsafe {
            compute_rhs(
                w,
                n as i32,
                g.as_ptr(),
                d.as_ptr(),
                n as i32,
                p,
                dp,
                b.as_mut_ptr(),
            )
        };
        cv(&a, &b);
    }
}
#[test]
fn test_rcdot() {
    let mut r = XorShift::new(5);
    for k in 0..2000 {
        let n = k % 31 + 1;
        let a: Vec<_> = (0..n)
            .map(|_| Complex {
                r: val(&mut r),
                i: val(&mut r),
            })
            .collect();
        let b: Vec<_> = (0..n)
            .map(|_| Complex {
                r: val(&mut r),
                i: val(&mut r),
            })
            .collect();
        close(rust::rcdot(&a, &b), unsafe {
            rcdot(n as i32, a.as_ptr(), b.as_ptr())
        });
    }
}
#[test]
fn test_htmul() {
    let mut r = XorShift::new(6);
    for k in 0..2000 {
        let n = k % 31 + 1;
        let a: Vec<_> = (0..n)
            .map(|_| Complex {
                r: val(&mut r),
                i: val(&mut r),
            })
            .collect();
        let x = a.clone();
        let mut y = vec![Complex::default(); n];
        let mut z = y.clone();
        rust::htmul(&a, &x, &mut y);
        unsafe { htmul(n as i32, a.as_ptr(), x.as_ptr(), z.as_mut_ptr()) };
        cv(&y, &z);
    }
}
#[test]
fn test_runav() {
    let mut r = XorShift::new(7);
    for k in 0..2000 {
        let n = k % 31 + 1;
        let len = (k % 65) as i32;
        let a: Vec<_> = (0..n).map(|_| val(&mut r)).collect();
        let mut b = vec![0.; n];
        let mut c = b.clone();
        rust::runav(len, &a, &mut b);
        unsafe { runav(n as i32, len, a.as_ptr(), c.as_mut_ptr()) };
        for (x, y) in b.iter().zip(c) {
            close(*x, y);
        }
    }
}
#[test]
fn test_taupmute() {
    let mut r = XorShift::new(8);
    for k in 0..2000 {
        let n = k % 31 + 1;
        let mut a: Vec<_> = (0..n).map(|_| val(&mut r)).collect();
        let mut b = a.clone();
        let ip = k as i32 % 25 - 5;
        rust::taupmute(ip, 0, 10, 7, &mut a);
        unsafe { taupmute(ip, 0, 10, n as i32, 7, b.as_mut_ptr()) };
        assert_eq!(a, b);
    }
}
#[test]
fn test_rho_filter() {
    let mut rng = XorShift::new(87);
    for k in 0..2000 {
        let nt = [3, 5, 16, 31, 64, 127][k % 6];
        let n = [1, 3, 5][k % 3];
        let dt = 0.001 + val(&mut rng).abs() * 0.02;
        let mut a = vec![0.; n];
        let mut b = a.clone();
        rust::rho_filter(nt, dt, &mut a);
        unsafe { rho_filter(n as i32, nt as i32, dt, b.as_mut_ptr()) };
        for (x, y) in a.iter().zip(b) {
            close(*x, y);
        }
    }
}
#[test]
fn test_fwd_tx_sstack() {
    let mut rng = XorShift::new(9);
    for k in 0..2000 {
        let nt = k % 31 + 2;
        let nx = k % 7 + 1;
        let np = k % 5 + 1;
        let tr: Vec<Vec<f32>> = (0..nx)
            .map(|_| (0..nt).map(|_| val(&mut rng)).collect())
            .collect();
        let mut a = vec![vec![0.; nt]; np];
        let mut b = a.clone();
        let t: Vec<_> = tr.iter().map(|x| x.as_ptr()).collect();
        let mut o: Vec<_> = b.iter_mut().map(|x| x.as_mut_ptr()).collect();
        rust::fwd_tx_sstack(1., nt, -1., 0.2, -0.1, 0.05, &tr, &mut a);
        unsafe {
            fwd_tx_sstack(
                1.,
                nt as i32,
                nx as i32,
                -1.,
                0.2,
                np as i32,
                -0.1,
                0.05,
                t.as_ptr(),
                o.as_mut_ptr(),
            )
        };
        for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
            close(*x, *y);
        }
    }
}
#[test]
fn test_inv_tx_sstack() {
    let mut rng = XorShift::new(12);
    for k in 0..2000 {
        let nt = k % 31 + 3;
        let nx = k % 7 + 1;
        let np = k % 5 + 1;
        let tr: Vec<Vec<f32>> = (0..np)
            .map(|_| (0..nt).map(|_| val(&mut rng)).collect())
            .collect();
        let mut a = vec![vec![0.; nt]; nx];
        let mut b = a.clone();
        let t: Vec<_> = tr.iter().map(|x| x.as_ptr()).collect();
        let mut o: Vec<_> = b.iter_mut().map(|x| x.as_mut_ptr()).collect();
        rust::inv_tx_sstack(1., nt, 3, -1., 0.2, -0.1, 0.05, &tr, &mut a);
        unsafe {
            inv_tx_sstack(
                1.,
                nt as i32,
                nx as i32,
                3,
                -1.,
                0.2,
                np as i32,
                -0.1,
                0.05,
                t.as_ptr(),
                o.as_mut_ptr(),
            )
        };
        let m = b.iter().flatten().fold(1e-20f32, |a, x| a.max(x.abs()));
        for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
            assert!((x - y).abs() < 1e-5 * m);
        }
    }
}
#[test]
fn test_ctoephcg() {
    let mut rng = XorShift::new(10);
    for k in 0..2000 {
        let n = k % 12 + 1;
        let mut a: Vec<_> = (0..n)
            .map(|_| Complex {
                r: val(&mut rng) * 0.01,
                i: val(&mut rng) * 0.01,
            })
            .collect();
        a[0] = Complex { r: 2., i: 0. };
        let y: Vec<_> = (0..n)
            .map(|_| Complex {
                r: val(&mut rng),
                i: val(&mut rng),
            })
            .collect();
        let mut w = vec![vec![Complex::default(); n]; 5];
        let mut c = w.clone();
        let [x, s, ss, g, rr]: &mut [Vec<Complex>; 5] = w.as_mut_slice().try_into().unwrap();
        let it = rust::ctoephcg(5, &a, x, &y, s, ss, g, rr);
        let [x, s, ss, g, rr]: &mut [Vec<Complex>; 5] = c.as_mut_slice().try_into().unwrap();
        let ci = unsafe {
            ctoephcg(
                5,
                n as i32,
                a.as_ptr(),
                x.as_mut_ptr(),
                y.as_ptr(),
                s.as_mut_ptr(),
                ss.as_mut_ptr(),
                g.as_mut_ptr(),
                rr.as_mut_ptr(),
            )
        };
        assert_eq!(it, ci);
        for (a, b) in w.iter().zip(c.iter()) {
            cv(a, b);
        }
    }
}
#[test]
fn test_ctoep() {
    let mut rng = XorShift::new(11);
    for k in 0..2000 {
        let n = k % 16 + 1;
        let mut r: Vec<_> = (0..n)
            .map(|_| Complex {
                r: val(&mut rng) * 0.01,
                i: val(&mut rng) * 0.01,
            })
            .collect();
        r[0] = Complex { r: 2., i: 0. };
        let b: Vec<_> = (0..n)
            .map(|_| Complex {
                r: val(&mut rng),
                i: val(&mut rng),
            })
            .collect();
        let mut a = vec![Complex::default(); n];
        let mut f = a.clone();
        let mut g = a.clone();
        let mut aa = a.clone();
        let mut ff = a.clone();
        let mut gg = a.clone();
        let x = rust::ctoep(&r, &mut a, &b, &mut f, &mut g);
        let y = unsafe {
            ctoep(
                n as i32,
                r.as_ptr(),
                aa.as_mut_ptr(),
                b.as_ptr(),
                ff.as_mut_ptr(),
                gg.as_mut_ptr(),
            )
        };
        assert_eq!(x, y);
        cv(&a, &aa);
        cv(&f, &ff);
        cv(&g, &gg);
    }
}
