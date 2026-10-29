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

//! Differential tests for su/lib/bilinear.c. Credits: none listed in the C source.
use oxwave_bilinear_interp as port;
use oxwave_cref::XorShift;

extern "C" {
    fn binterpfind(
        ki: i32,
        ai: *const i32,
        ni: i32,
        ei: i32,
        kc: i32,
        ac: *const i32,
        nc: i32,
        ec: i32,
        oi: *mut i32,
        oc: *mut i32,
        wi: *mut f64,
        wc: *mut f64,
    );
    fn binterpapply(
        a: *const f64,
        b: *const f64,
        ni: i32,
        wi: f64,
        c: *const f64,
        d: *const f64,
        nc: i32,
        wc: f64,
        n: i32,
        out: *mut f64,
    );
    fn binterpvalue(
        x: f64,
        mode: i32,
        ax: *const f64,
        ay: *const f64,
        an: i32,
        bx: *const f64,
        by: *const f64,
        bn: i32,
        ni: i32,
        wi: f64,
        cx: *const f64,
        cy: *const f64,
        cn: i32,
        dx: *const f64,
        dy: *const f64,
        dn: i32,
        nc: i32,
        wc: f64,
        out: *mut f64,
    );
}
const SIZES: [usize; 9] = [1, 2, 3, 7, 16, 31, 64, 127, 1024];
fn close(a: f64, b: f64) {
    assert!(
        a.to_bits() == b.to_bits()
            || (a.is_nan() && b.is_nan())
            || (a - b).abs() <= 1e-5 * a.abs().max(b.abs()),
        "{a:?} != {b:?}"
    );
}
fn weight(r: &mut XorShift, k: usize) -> f64 {
    match k % 5 {
        0 => 0.0,
        1 => 1.0,
        2 => -0.0,
        _ => r.uniform(-2.0, 3.0),
    }
}
#[test]
fn find_matches_c() {
    let mut r = XorShift::new(101);
    for &n in &SIZES {
        for &m in &SIZES {
            let a: Vec<i32> = (0..n).map(|i| i as i32 * 3 - 150).collect();
            let b: Vec<i32> = (0..m).map(|i| i as i32 * 5 - 100).collect();
            for k in 0..160 {
                let ki = if k % 3 == 0 {
                    a[k % n]
                } else {
                    r.range_i32(-300, 4000)
                };
                let kc = if k % 4 == 0 {
                    b[k % m]
                } else {
                    r.range_i32(-300, 6000)
                };
                let ei = (k % 4) as i32;
                let ec = ((k / 4) % 4) as i32;
                let (mut oi, mut oc, mut wi, mut wc) = (0, 0, 0.0, 0.0);
                let (mut ci, mut cc, mut cw, mut cd) = (0, 0, 0.0, 0.0);
                port::binterpfind(ki, &a, ei, kc, &b, ec, &mut oi, &mut oc, &mut wi, &mut wc);
                unsafe {
                    binterpfind(
                        ki,
                        a.as_ptr(),
                        n as i32,
                        ei,
                        kc,
                        b.as_ptr(),
                        m as i32,
                        ec,
                        &mut ci,
                        &mut cc,
                        &mut cw,
                        &mut cd,
                    );
                }
                assert_eq!((oi, oc), (ci, cc));
                close(wi, cw);
                close(wc, cd);
            }
        }
    }
}
#[test]
fn apply_matches_c() {
    let mut r = XorShift::new(102);
    for &n in &SIZES {
        for k in 0..500 {
            let arrays: Vec<Vec<f64>> = (0..4)
                .map(|_| {
                    (0..n)
                        .map(|_| {
                            if k % 13 == 0 {
                                0.0
                            } else {
                                r.uniform(-1000.0, 1000.0)
                            }
                        })
                        .collect()
                })
                .collect();
            let ni = if k % 2 == 0 { 1 } else { 7 };
            let nc = if k % 4 < 2 { 1 } else { 16 };
            let wi = weight(&mut r, k);
            let wc = weight(&mut r, k / 5);
            let mut out = vec![99.0; n];
            let mut c = out.clone();
            port::binterpapply(
                &arrays[0], &arrays[1], ni, wi, &arrays[2], &arrays[3], nc, wc, &mut out,
            );
            unsafe {
                binterpapply(
                    arrays[0].as_ptr(),
                    arrays[1].as_ptr(),
                    ni,
                    wi,
                    arrays[2].as_ptr(),
                    arrays[3].as_ptr(),
                    nc,
                    wc,
                    n as i32,
                    c.as_mut_ptr(),
                )
            };
            for i in 0..n {
                close(out[i], c[i]);
            }
        }
    }
    let mut out = [];
    port::binterpapply(&[], &[], 2, 0.5, &[], &[], 2, 0.5, &mut out);
    unsafe {
        binterpapply(
            out.as_ptr(),
            out.as_ptr(),
            2,
            0.5,
            out.as_ptr(),
            out.as_ptr(),
            2,
            0.5,
            0,
            out.as_mut_ptr(),
        )
    };
}
#[test]
fn value_matches_c() {
    let mut r = XorShift::new(103);
    for &n in &SIZES {
        for k in 0..500 {
            let sizes = [
                n,
                SIZES[(k + 1) % 9],
                SIZES[(k + 3) % 9],
                SIZES[(k + 5) % 9],
            ];
            let xs: Vec<Vec<f64>> = sizes
                .iter()
                .enumerate()
                .map(|(j, &s)| (0..s).map(|i| i as f64 * 2.5 - 20.0 + j as f64).collect())
                .collect();
            let ys: Vec<Vec<f64>> = sizes
                .iter()
                .map(|&s| {
                    (0..s)
                        .map(|_| {
                            if k % 11 == 0 {
                                0.0
                            } else {
                                r.uniform(-1000.0, 1000.0)
                            }
                        })
                        .collect()
                })
                .collect();
            let x = if k % 4 == 0 {
                xs[0][k % n]
            } else {
                r.uniform(-100.0, 2700.0)
            };
            let mode = (k % 4) as i32;
            let ni = if k % 2 == 0 { 1 } else { 7 };
            let nc = if k % 4 < 2 { 1 } else { 16 };
            let wi = weight(&mut r, k / 4);
            let wc = weight(&mut r, k / 20);
            let (mut out, mut c) = (99.0, 99.0);
            port::binterpvalue(
                x, mode, &xs[0], &ys[0], &xs[1], &ys[1], ni, wi, &xs[2], &ys[2], &xs[3], &ys[3],
                nc, wc, &mut out,
            );
            unsafe {
                binterpvalue(
                    x,
                    mode,
                    xs[0].as_ptr(),
                    ys[0].as_ptr(),
                    sizes[0] as i32,
                    xs[1].as_ptr(),
                    ys[1].as_ptr(),
                    sizes[1] as i32,
                    ni,
                    wi,
                    xs[2].as_ptr(),
                    ys[2].as_ptr(),
                    sizes[2] as i32,
                    xs[3].as_ptr(),
                    ys[3].as_ptr(),
                    sizes[3] as i32,
                    nc,
                    wc,
                    &mut c,
                )
            };
            close(out, c);
        }
    }
}
