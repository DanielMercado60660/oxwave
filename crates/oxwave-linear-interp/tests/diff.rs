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

//! Differential tests for cwp/lib/intlin.c, intlinc.c and su/lib/linterpd.c.
//! Credits: Dave Hale, 1989; complex adaptation Hans Ecke, 2002; linterpd lists none.
use oxwave_cref::XorShift;
use oxwave_cwp::complex::Complex;
use oxwave_linear_interp as port;
use std::sync::Mutex;
static LOCK: Mutex<()> = Mutex::new(());
extern "C" {
    fn bhighd(a: *const f64, n: i32, x: f64) -> i32;
    fn bhighi(a: *const i32, n: i32, x: i32) -> i32;
    fn intlin(
        n: i32,
        x: *const f32,
        y: *const f32,
        l: f32,
        r: f32,
        m: i32,
        q: *const f32,
        o: *mut f32,
    );
    fn intlinc(
        n: i32,
        x: *const f32,
        y: *const Complex,
        l: Complex,
        r: Complex,
        m: i32,
        q: *const f32,
        o: *mut Complex,
    );
    fn linterpd(x: f64, a: *const f64, y: *const f64, n: i32, mode: i32, o: *mut f64);
}
const SIZES: [usize; 9] = [1, 2, 3, 7, 16, 31, 64, 127, 1024];
const NAN_BITS: [u64; 4] = [
    0x7ff8_0000_0000_0000,
    0xfff8_0000_0000_0000,
    0x7ff8_0000_0000_0042,
    0xfff8_0000_0000_0042,
];
#[test]
fn bhighd_matches() {
    let mut rng = XorShift::new(1);
    for n in SIZES {
        let a: Vec<f64> = (0..n).map(|i| (i as i32 / 3 - 100) as f64).collect();
        for k in 0..2000 {
            let x = if k % 3 == 0 {
                a[k % n]
            } else {
                rng.uniform(-150.0, 400.0)
            };
            assert_eq!(port::bhighd(&a, x), unsafe {
                bhighd(a.as_ptr(), n as i32, x)
            });
        }
    }
    assert_eq!(port::bhighd(&[], 0.0), unsafe {
        bhighd(std::ptr::null(), 0, 0.0)
    });
}
#[test]
fn bhighi_matches() {
    let mut rng = XorShift::new(2);
    for n in SIZES {
        let a: Vec<i32> = (0..n).map(|i| i as i32 / 3 - 100).collect();
        for k in 0..2000 {
            let x = if k % 3 == 0 {
                a[k % n]
            } else {
                rng.range_i32(-150, 400)
            };
            assert_eq!(port::bhighi(&a, x), unsafe {
                bhighi(a.as_ptr(), n as i32, x)
            });
        }
    }
    assert_eq!(port::bhighi(&[], 0), unsafe {
        bhighi(std::ptr::null(), 0, 0)
    });
}
fn inputs(r: &mut XorShift, n: usize, reverse: bool) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let mut x: Vec<f32> = (0..n).map(|i| i as f32 * 0.75 - 100.0).collect();
    if reverse {
        x.reverse();
    }
    let y = (0..n).map(|_| r.uniform_f32(-1000.0, 1000.0)).collect();
    let mut q: Vec<f32> = (0..400)
        .map(|_| r.uniform_f32(-110.0, n as f32 * 0.75 - 90.0))
        .collect();
    q.extend_from_slice(&x);
    q.extend_from_slice(&[0.0, -100.0, -1000.0, 10000.0]);
    (x, y, q)
}
#[test]
fn intlin_matches() {
    let _g = LOCK.lock().unwrap();
    let mut r = XorShift::new(3);
    for n in SIZES {
        for rev in [false, true] {
            for _ in 0..5 {
                let (x, y, q) = inputs(&mut r, n, rev);
                let mut c = vec![0.0; q.len()];
                let mut a = c.clone();
                unsafe {
                    intlin(
                        n as i32,
                        x.as_ptr(),
                        y.as_ptr(),
                        -3.0,
                        7.0,
                        q.len() as i32,
                        q.as_ptr(),
                        c.as_mut_ptr(),
                    )
                };
                port::intlin(&x, &y, -3.0, 7.0, &q, &mut a);
                for (v, w) in a.iter().zip(&c) {
                    assert_eq!(v.to_bits(), w.to_bits());
                }
            }
        }
    }
}
#[test]
fn intlinc_matches() {
    let _g = LOCK.lock().unwrap();
    let mut r = XorShift::new(4);
    for n in SIZES {
        for rev in [false, true] {
            for _ in 0..5 {
                let (x, y, q) = inputs(&mut r, n, rev);
                let y: Vec<_> = y
                    .iter()
                    .map(|&v| Complex {
                        r: v,
                        i: r.uniform_f32(-500.0, 500.0),
                    })
                    .collect();
                let l = Complex { r: -3.0, i: 0.0 };
                let rr = Complex { r: 7.0, i: -4.0 };
                let mut c = vec![Complex::default(); q.len()];
                let mut a = c.clone();
                unsafe {
                    intlinc(
                        n as i32,
                        x.as_ptr(),
                        y.as_ptr(),
                        l,
                        rr,
                        q.len() as i32,
                        q.as_ptr(),
                        c.as_mut_ptr(),
                    )
                };
                port::intlinc(&x, &y, l, rr, &q, &mut a);
                for (v, w) in a.iter().zip(&c) {
                    assert_eq!(v.r.to_bits(), w.r.to_bits());
                    assert_eq!(v.i.to_bits(), w.i.to_bits());
                }
            }
        }
    }
}
#[test]
fn linterpd_matches() {
    let mut r = XorShift::new(5);
    for n in SIZES {
        let x: Vec<f64> = (0..n).map(|i| i as f64 * 0.75 - 100.0).collect();
        let y: Vec<f64> = (0..n).map(|_| r.uniform(-1000.0, 1000.0)).collect();
        for mode in 0..=4 {
            for k in 0..2000 {
                let q = if k % 4 == 0 {
                    x[k % n]
                } else {
                    r.uniform(-120.0, n as f64 * 0.75 - 80.0)
                };
                let mut c = 0.0;
                let mut a = 0.0;
                unsafe { linterpd(q, x.as_ptr(), y.as_ptr(), n as i32, mode, &mut c) };
                port::linterpd(q, &x, &y, mode, &mut a);
                assert_eq!(a.to_bits(), c.to_bits(), "n={n} mode={mode} q={q}");
            }
        }
    }
}

#[test]
fn linterpd_nan_single_entry_matches() {
    // C returns the sole time before inspecting the offset, including NaNs.
    for bits in NAN_BITS {
        let offset = f64::from_bits(bits);
        for value in [-0.0, 0.0, -7.0, 42.0] {
            for mode in [-1, 0, 1, 2, 3, 4] {
                let offs = [0.0];
                let tims = [value];
                let mut c = 0.0;
                let mut a = 0.0;
                unsafe { linterpd(offset, offs.as_ptr(), tims.as_ptr(), 1, mode, &mut c) };
                port::linterpd(offset, &offs, &tims, mode, &mut a);
                assert_eq!(a.to_bits(), c.to_bits());
            }
        }
    }
}

#[test]
fn linterpd_nan_offset_exits_cleanly() {
    const CHILD: &str = "OXWAVE_LINEAR_INTERP_NAN_CHILD";
    if let Ok(case) = std::env::var(CHILD) {
        let mut fields = case.split(',');
        let n = fields.next().unwrap();
        let mode = fields.next().unwrap();
        let bits = fields.next().unwrap().parse::<u64>().unwrap();
        let n: usize = n.parse().unwrap();
        let mode: i32 = mode.parse().unwrap();
        oxwave_core::par::initargs(&["linear-interp-test".to_owned()]);
        let offs: Vec<_> = (0..n).map(|i| i as f64 - 1.0).collect();
        let tims = vec![0.0; n];
        let mut time = 0.0;
        // Do not invoke C here: its out-of-bounds accesses are undefined.
        port::linterpd(f64::from_bits(bits), &offs, &tims, mode, &mut time);
        panic!("NaN offset must exit via err!");
    }

    for n in [2, 7, 16] {
        for mode in [-1, 0, 1, 2, 3, 4] {
            for bits in NAN_BITS {
                let output = std::process::Command::new(std::env::current_exe().unwrap())
                    .args([
                        "--exact",
                        "linterpd_nan_offset_exits_cleanly",
                        "--nocapture",
                    ])
                    .env(CHILD, format!("{n},{mode},{bits}"))
                    .output()
                    .unwrap();
                assert_eq!(
                    output.status.code(),
                    Some(1),
                    "n={n} mode={mode} bits={bits:#x}"
                );
                assert_eq!(
                    output.stderr, b"\nlinear-interp-test: linterpd: offset is NaN\n",
                    "n={n} mode={mode} bits={bits:#x}"
                );
            }
        }
    }
}
