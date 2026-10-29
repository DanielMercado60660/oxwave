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

//! Differential tests against cwp/lib/intsinc8.c, inttable8.c, mksinc.c,
//! ressinc8.c, and shfs8r.c.
//! Credits: Dave Hale, Colorado School of Mines, 06/02/89;
//! ressinc8.c: Dave Hale, Colorado School of Mines, 06/06/90.

use oxwave_cref::XorShift;
use oxwave_cwp::complex::Complex;
use oxwave_sinc_interp as rust;
use std::sync::Mutex;

// The original C functions lazily initialize mutable static tables.
static C_LOCK: Mutex<()> = Mutex::new(());
const CASES: usize = 4096;
const SIZES: &[usize] = &[
    1, 2, 3, 5, 7, 8, 9, 16, 17, 31, 64, 127, 256, 509, 1024, 4093,
];
const GUARD: f32 = 1234567.0;

extern "C" {
    fn mksinc(d: f32, lsinc: i32, sinc: *mut f32);
    fn ints8r(
        nxin: i32,
        dxin: f32,
        fxin: f32,
        yin: *mut f32,
        yinl: f32,
        yinr: f32,
        nxout: i32,
        xout: *mut f32,
        yout: *mut f32,
    );
    fn ints8c(
        nxin: i32,
        dxin: f32,
        fxin: f32,
        yin: *mut Complex,
        yinl: Complex,
        yinr: Complex,
        nxout: i32,
        xout: *mut f32,
        yout: *mut Complex,
    );
    fn intt8r(
        ntable: i32,
        table: *mut [f32; 8],
        nxin: i32,
        dxin: f32,
        fxin: f32,
        yin: *mut f32,
        yinl: f32,
        yinr: f32,
        nxout: i32,
        xout: *mut f32,
        yout: *mut f32,
    );
    fn intt8c(
        ntable: i32,
        table: *mut [f32; 8],
        nxin: i32,
        dxin: f32,
        fxin: f32,
        yin: *mut Complex,
        yinl: Complex,
        yinr: Complex,
        nxout: i32,
        xout: *mut f32,
        yout: *mut Complex,
    );
    fn ress8r(
        nxin: i32,
        dxin: f32,
        fxin: f32,
        yin: *mut f32,
        yinl: f32,
        yinr: f32,
        nxout: i32,
        dxout: f32,
        fxout: f32,
        yout: *mut f32,
    );
    fn ress8c(
        nxin: i32,
        dxin: f32,
        fxin: f32,
        yin: *mut Complex,
        yinl: Complex,
        yinr: Complex,
        nxout: i32,
        dxout: f32,
        fxout: f32,
        yout: *mut Complex,
    );
    fn shfs8r(
        dx: f32,
        nxin: i32,
        fxin: f32,
        yin: *mut f32,
        yinl: f32,
        yinr: f32,
        nxout: i32,
        fxout: f32,
        yout: *mut f32,
    );
}

// This checks 1e-5 relative to the largest magnitude in the C output,
// with no absolute tolerance floor. All generated samples are finite.
// Also count bit differences to report the actual precision.
fn compare(a: &[f32], b: &[f32], case: usize) -> usize {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().fold(0.0f64, |s, &x| s.max((x as f64).abs()));
    let mut different = 0;
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        assert!(x.is_finite() && y.is_finite());
        assert!(
            (x as f64 - y as f64).abs() <= 1e-5 * scale,
            "case={case} sample={i} Rust={x:?} C={y:?} scale={scale:e}"
        );
        different += usize::from(x.to_bits() != y.to_bits());
    }
    different
}

fn compare_complex(a: &[Complex], b: &[Complex], case: usize) -> usize {
    let a: Vec<f32> = a.iter().flat_map(|z| [z.r, z.i]).collect();
    let b: Vec<f32> = b.iter().flat_map(|z| [z.r, z.i]).collect();
    compare(&a, &b, case)
}

struct Data {
    yin: Vec<f32>,
    cyin: Vec<Complex>,
    left: f32,
    right: f32,
    cleft: Complex,
    cright: Complex,
    dx: f32,
    fx: f32,
    xout: Vec<f32>,
}

impl Data {
    fn new(rng: &mut XorShift, case: usize) -> Self {
        let n = SIZES[case % SIZES.len()];
        // Include no-output calls, then all small/prime/power-of-two/large sizes.
        let m = if case % 37 == 0 {
            0
        } else {
            SIZES[(case / SIZES.len()) % SIZES.len()]
        };
        // Vary scales independently of input and output sizes.
        let amplitude = match rng.next_u32() % 8 {
            0 => 0.0,
            1 => 1e-30,
            2 => 1e30,
            _ => 100.0,
        };
        let mut yin = Vec::with_capacity(n);
        let mut cyin = Vec::with_capacity(n);
        for i in 0..n {
            let (r, im) = match case % 11 {
                0 => (0.0, -0.0),
                1 => (if i == n / 2 { amplitude } else { 0.0 }, 0.0),
                2 => (-amplitude, amplitude),
                3 => (
                    (i as f64 * 0.37).sin() as f32 * amplitude,
                    (i as f64 * 0.53).cos() as f32 * amplitude,
                ),
                _ => (
                    rng.uniform_f32(-1.0, 1.0) * amplitude,
                    rng.uniform_f32(-1.0, 1.0) * amplitude,
                ),
            };
            yin.push(r);
            cyin.push(Complex { r, i: im });
        }
        let left = rng.uniform_f32(-1.0, 1.0) * amplitude;
        let right = rng.uniform_f32(-1.0, 1.0) * amplitude;
        let cleft = Complex {
            r: left,
            i: rng.uniform_f32(-1.0, 1.0) * amplitude,
        };
        let cright = Complex {
            r: right,
            i: rng.uniform_f32(-1.0, 1.0) * amplitude,
        };
        let dx = match (case + case / SIZES.len()) % 12 {
            0 => 1.0,
            1 => -1.0,
            2 => 0.004,
            3 => -0.004,
            4 => 1e-5,
            5 => 1e4,
            _ => rng.uniform_f32(0.001, 10.0) * if case % 2 == 0 { 1.0 } else { -1.0 },
        };
        // Keep tiny intervals resolvable at nonzero origins.
        let fx = if case % 5 == 0 {
            0.0
        } else {
            rng.uniform_f32(-20.0, 20.0) * dx
        };
        let mut xout = Vec::with_capacity(m);
        for i in 0..m {
            let sample = match i % 20 {
                0 => -1_000_000.25,
                1 => 1_000_000.25,
                2 => -8.5,
                3 => -8.0,
                4 => -7.5,
                5 => -1.0,
                6 => 0.0,
                7 => 0.5,
                8 => (n - 1) as f32,
                9 => n as f32,
                10 => 3.0,
                11 => n as f32 - 4.0,
                // Table midpoint ties and their immediate f32 neighbours.
                12 => 4.0 + 0.5 / 512.0,
                13 => f32::from_bits((4.0f32 + 0.5 / 512.0).to_bits() - 1),
                14 => f32::from_bits((4.0f32 + 0.5 / 512.0).to_bits() + 1),
                _ => rng.uniform_f32(-16.0, n as f32 + 16.0),
            };
            xout.push(fx + sample * dx);
        }
        Self {
            yin,
            cyin,
            left,
            right,
            cleft,
            cright,
            dx,
            fx,
            xout,
        }
    }

    fn resample_grid(&self, rng: &mut XorShift, case: usize) -> (f32, f32) {
        // Each mode gets every pair of input/output sizes.
        let (ratio, offset) = match (case / (SIZES.len() * SIZES.len())) % 16 {
            0 => (1.0, 0.0),
            1 => (2.0, -4.0),
            2 => (-1.0, self.yin.len() as f32),
            3 => (0.0, 0.0099),
            4 => (0.0, 0.0101),
            5 => (0.0, -0.0099),
            6 => (0.0, -0.0101),
            7 => (0.0001, 0.0), // initially on grid, then requires full interpolation
            8 => (0.0, 0.01),
            9 => (0.0, -0.01),
            10 => (1.0, -1_000_000.25),
            11 => (1.0, 1_000_000.25),
            12 => (0.5, -0.5),
            _ => (
                rng.uniform_f32(-3.0, 3.0),
                rng.uniform_f32(-20.0, self.yin.len() as f32 + 20.0),
            ),
        };
        (ratio * self.dx, self.fx + offset * self.dx)
    }
}

// Exercise generic (not necessarily sinc) coefficients too. Generate sinc
// rows through the ORIGINAL C implementation so intt8 tests are independent
// of the Rust mksinc. Endpoint rows follow the documented table contract.
fn make_table(rng: &mut XorShift, case: usize) -> Vec<[f32; 8]> {
    let nt = [2, 3, 7, 16, 31, 64, 127, 513, 1024][case % 9];
    let mut table = vec![[0.0; 8]; nt];
    for j in 1..nt - 1 {
        if case % 3 == 0 {
            for i in 0..8 {
                table[j][i] = rng.uniform_f32(-2.0, 2.0);
            }
        } else {
            unsafe { mksinc(j as f32 / (nt - 1) as f32, 8, table[j].as_mut_ptr()) };
        }
    }
    table[0][3] = 1.0;
    table[nt - 1][4] = 1.0;
    table
}

#[test]
fn mksinc_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x51ac8);
    let mut different = 0;
    // Every supported length, including 1, odd/prime sizes and powers of two;
    // beyond 20 the original C writes past its fixed-size work arrays.
    for case in 0..12000 {
        let n = case % 20 + 1;
        let d = match (case / 20) % 12 {
            0 => 0.0,
            1 => -0.0,
            2 => 1.0,
            3 => 0.5,
            4 => f32::from_bits(1),
            5 => f32::from_bits(1.0f32.to_bits() - 1),
            6 => rng.uniform_f32(-2.0, 0.0),
            7 => rng.uniform_f32(1.0, 3.0),
            _ => rng.uniform_f32(0.0, 1.0),
        };
        let mut expected = vec![GUARD; n + 2];
        let mut actual = expected.clone();
        unsafe { mksinc(d, n as i32, expected.as_mut_ptr()) };
        rust::mksinc(d, &mut actual[..n]);
        different += compare(&actual[..n], &expected[..n], case);
        assert_eq!(&actual[n..], &[GUARD; 2]);
        assert_eq!(&expected[n..], &[GUARD; 2]);
    }
    eprintln!("mksinc: {different} scalar outputs differ in bits");
    assert_eq!(
        different, 0,
        "mksinc must be bit-exact on these finite inputs"
    );
}

#[test]
fn ints8r_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(8442);
    let mut different = 0;
    for case in 0..CASES {
        let d = Data::new(&mut rng, case);
        let m = d.xout.len();
        let mut expected = vec![GUARD; m + 2];
        let mut actual = expected.clone();
        unsafe {
            ints8r(
                d.yin.len() as i32,
                d.dx,
                d.fx,
                d.yin.as_ptr() as *mut f32,
                d.left,
                d.right,
                m as i32,
                d.xout.as_ptr() as *mut f32,
                expected.as_mut_ptr(),
            )
        };
        rust::ints8r(d.dx, d.fx, &d.yin, d.left, d.right, &d.xout, &mut actual);
        different += compare(&actual[..m], &expected[..m], case);
        assert_eq!(&actual[m..], &[GUARD; 2]);
        assert_eq!(&expected[m..], &[GUARD; 2]);
    }
    eprintln!("ints8r: {different} scalar outputs differ in bits");
    assert_eq!(
        different, 0,
        "ints8r must be bit-exact on these finite inputs"
    );
}

#[test]
fn ints8c_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(9257);
    let mut different = 0;
    for case in 0..CASES {
        let d = Data::new(&mut rng, case);
        let m = d.xout.len();
        let mut expected = vec![
            Complex {
                r: GUARD,
                i: -GUARD
            };
            m + 2
        ];
        let mut actual = expected.clone();
        unsafe {
            ints8c(
                d.cyin.len() as i32,
                d.dx,
                d.fx,
                d.cyin.as_ptr() as *mut Complex,
                d.cleft,
                d.cright,
                m as i32,
                d.xout.as_ptr() as *mut f32,
                expected.as_mut_ptr(),
            )
        };
        rust::ints8c(d.dx, d.fx, &d.cyin, d.cleft, d.cright, &d.xout, &mut actual);
        different += compare_complex(&actual[..m], &expected[..m], case);
        assert_eq!(
            &actual[m..],
            &[Complex {
                r: GUARD,
                i: -GUARD
            }; 2]
        );
        assert_eq!(
            &expected[m..],
            &[Complex {
                r: GUARD,
                i: -GUARD
            }; 2]
        );
    }
    eprintln!("ints8c: {different} scalar outputs differ in bits");
    assert_eq!(
        different, 0,
        "ints8c must be bit-exact on these finite inputs"
    );
}

#[test]
fn intt8r_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(10169);
    let mut different = 0;
    for case in 0..CASES {
        let d = Data::new(&mut rng, case);
        let m = d.xout.len();
        let mut expected = vec![GUARD; m + 2];
        let mut actual = expected.clone();
        let table = make_table(&mut rng, case);
        unsafe {
            intt8r(
                table.len() as i32,
                table.as_ptr() as *mut [f32; 8],
                d.yin.len() as i32,
                d.dx,
                d.fx,
                d.yin.as_ptr() as *mut f32,
                d.left,
                d.right,
                m as i32,
                d.xout.as_ptr() as *mut f32,
                expected.as_mut_ptr(),
            )
        };
        rust::intt8r(
            &table,
            d.dx,
            d.fx,
            &d.yin,
            d.left,
            d.right,
            &d.xout,
            &mut actual,
        );
        different += compare(&actual[..m], &expected[..m], case);
        assert_eq!(&actual[m..], &[GUARD; 2]);
        assert_eq!(&expected[m..], &[GUARD; 2]);
    }
    eprintln!("intt8r: {different} scalar outputs differ in bits");
    assert_eq!(
        different, 0,
        "intt8r must be bit-exact on these finite inputs"
    );
}

#[test]
fn intt8c_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(11094);
    let mut different = 0;
    for case in 0..CASES {
        let d = Data::new(&mut rng, case);
        let m = d.xout.len();
        let mut expected = vec![
            Complex {
                r: GUARD,
                i: -GUARD
            };
            m + 2
        ];
        let mut actual = expected.clone();
        let table = make_table(&mut rng, case);
        unsafe {
            intt8c(
                table.len() as i32,
                table.as_ptr() as *mut [f32; 8],
                d.cyin.len() as i32,
                d.dx,
                d.fx,
                d.cyin.as_ptr() as *mut Complex,
                d.cleft,
                d.cright,
                m as i32,
                d.xout.as_ptr() as *mut f32,
                expected.as_mut_ptr(),
            )
        };
        rust::intt8c(
            &table,
            d.dx,
            d.fx,
            &d.cyin,
            d.cleft,
            d.cright,
            &d.xout,
            &mut actual,
        );
        different += compare_complex(&actual[..m], &expected[..m], case);
        assert_eq!(
            &actual[m..],
            &[Complex {
                r: GUARD,
                i: -GUARD
            }; 2]
        );
        assert_eq!(
            &expected[m..],
            &[Complex {
                r: GUARD,
                i: -GUARD
            }; 2]
        );
    }
    eprintln!("intt8c: {different} scalar outputs differ in bits");
    assert_eq!(
        different, 0,
        "intt8c must be bit-exact on these finite inputs"
    );
}

#[test]
fn ress8r_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(12116);
    let mut different = 0;
    for case in 0..CASES {
        let d = Data::new(&mut rng, case);
        let m = d.xout.len();
        let mut expected = vec![GUARD; m + 2];
        let mut actual = expected.clone();
        let (dxout, fxout) = d.resample_grid(&mut rng, case);
        unsafe {
            ress8r(
                d.yin.len() as i32,
                d.dx,
                d.fx,
                d.yin.as_ptr() as *mut f32,
                d.left,
                d.right,
                m as i32,
                dxout,
                fxout,
                expected.as_mut_ptr(),
            )
        };
        rust::ress8r(
            d.dx,
            d.fx,
            &d.yin,
            d.left,
            d.right,
            dxout,
            fxout,
            &mut actual[..m],
        );
        different += compare(&actual[..m], &expected[..m], case);
        assert_eq!(&actual[m..], &[GUARD; 2]);
        assert_eq!(&expected[m..], &[GUARD; 2]);
        // Guaranteed exact grid copying, including signed zero, must
        // preserve every sample bit (also tests the all-on-grid fast path).
        if (case / (SIZES.len() * SIZES.len())) % 16 == 0 && d.dx.abs() == 1.0 {
            for i in 0..m {
                let wanted = if i < d.yin.len() { d.yin[i] } else { d.right };
                assert_eq!(actual[i].to_bits(), wanted.to_bits());
            }
        }
    }
    eprintln!("ress8r: {different} scalar outputs differ in bits");
    assert_eq!(
        different, 0,
        "ress8r must be bit-exact on these finite inputs"
    );
}

#[test]
fn ress8c_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(13383);
    let mut different = 0;
    for case in 0..CASES {
        let d = Data::new(&mut rng, case);
        let m = d.xout.len();
        let mut expected = vec![
            Complex {
                r: GUARD,
                i: -GUARD
            };
            m + 2
        ];
        let mut actual = expected.clone();
        let (dxout, fxout) = d.resample_grid(&mut rng, case);
        unsafe {
            ress8c(
                d.cyin.len() as i32,
                d.dx,
                d.fx,
                d.cyin.as_ptr() as *mut Complex,
                d.cleft,
                d.cright,
                m as i32,
                dxout,
                fxout,
                expected.as_mut_ptr(),
            )
        };
        rust::ress8c(
            d.dx,
            d.fx,
            &d.cyin,
            d.cleft,
            d.cright,
            dxout,
            fxout,
            &mut actual[..m],
        );
        different += compare_complex(&actual[..m], &expected[..m], case);
        assert_eq!(
            &actual[m..],
            &[Complex {
                r: GUARD,
                i: -GUARD
            }; 2]
        );
        assert_eq!(
            &expected[m..],
            &[Complex {
                r: GUARD,
                i: -GUARD
            }; 2]
        );
        // Guaranteed exact grid copying, including signed zero, must
        // preserve every sample bit (also tests the all-on-grid fast path).
        if (case / (SIZES.len() * SIZES.len())) % 16 == 0 && d.dx.abs() == 1.0 {
            for i in 0..m {
                let wanted = if i < d.cyin.len() {
                    d.cyin[i]
                } else {
                    d.cright
                };
                assert_eq!(actual[i].r.to_bits(), wanted.r.to_bits());
                assert_eq!(actual[i].i.to_bits(), wanted.i.to_bits());
            }
        }
    }
    eprintln!("ress8c: {different} scalar outputs differ in bits");
    assert_eq!(
        different, 0,
        "ress8c must be bit-exact on these finite inputs"
    );
}

#[test]
fn shfs8r_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(14825);
    let mut different = 0;
    for case in 0..CASES {
        let d = Data::new(&mut rng, case);
        let m = d.xout.len();
        let mut expected = vec![GUARD; m + 2];
        let mut actual = expected.clone();
        let shift = match case % 14 {
            0 => 0.0,
            1 => -1.0,
            2 => 1.0,
            3 => -4.0,
            4 => -0.5,
            5 => 0.5,
            6 => -2_000_000.3,
            7 => 2_000_000.3,
            8 => 0.5 / 512.0,
            9 => f32::from_bits((0.5f32 / 512.0).to_bits() - 1),
            10 => f32::from_bits((0.5f32 / 512.0).to_bits() + 1),
            _ => rng.uniform_f32(-20.0, d.yin.len() as f32 + 20.0),
        };
        let fxout = d.fx + shift * d.dx;
        unsafe {
            shfs8r(
                d.dx,
                d.yin.len() as i32,
                d.fx,
                d.yin.as_ptr() as *mut f32,
                d.left,
                d.right,
                m as i32,
                fxout,
                expected.as_mut_ptr(),
            )
        };
        rust::shfs8r(d.dx, d.fx, &d.yin, d.left, d.right, fxout, &mut actual[..m]);
        different += compare(&actual[..m], &expected[..m], case);
        assert_eq!(&actual[m..], &[GUARD; 2]);
        assert_eq!(&expected[m..], &[GUARD; 2]);
    }
    eprintln!("shfs8r: {different} scalar outputs differ in bits");
    assert_eq!(
        different, 0,
        "shfs8r must be bit-exact on these finite inputs"
    );
}
