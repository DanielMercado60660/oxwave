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

//! Differential tests for cwp/lib/sblas.c. Credits: Dave Hale, Colorado School of Mines, 10/01/89.
use oxwave_cref::XorShift;
use oxwave_sblas as rust;

extern "C" {
    fn isamax(n: i32, x: *const f32, inc: i32) -> i32;
    fn sasum(n: i32, x: *const f32, inc: i32) -> f32;
    fn saxpy(n: i32, a: f32, x: *const f32, ix: i32, y: *mut f32, iy: i32);
    fn scopy(n: i32, x: *const f32, ix: i32, y: *mut f32, iy: i32);
    fn sdot(n: i32, x: *const f32, ix: i32, y: *const f32, iy: i32) -> f32;
    fn snrm2(n: i32, x: *const f32, inc: i32) -> f32;
    fn sscal(n: i32, a: f32, x: *mut f32, inc: i32);
    fn sswap(n: i32, x: *mut f32, ix: i32, y: *mut f32, iy: i32);
}

fn cases(mut f: impl FnMut(i32, i32, i32, f32, Vec<f32>, Vec<f32>)) {
    let mut rng = XorShift::new(0x1234567);
    let sizes = [-2, 0, 1, 2, 3, 5, 7, 11, 16, 31, 64, 127, 256, 1024];
    let strides = [-3, -1, 0, 1, 2, 3];
    for t in 0..3500 {
        let n = sizes[t % sizes.len()];
        let ix = strides[(t / sizes.len()) % strides.len()];
        let iy = strides[(t / (sizes.len() * strides.len())) % strides.len()];
        let len = 1 + (n.max(1) as usize - 1) * 3;
        let scale = match t % 5 {
            0 => 0.0,
            1 => 1e-15,
            2 => 1e15,
            _ => 100.0,
        };
        let mut v = || (rng.next_f64() as f32 * 2.0 - 1.0) * scale;
        let x = (0..len).map(|_| v()).collect();
        let y = (0..len).map(|_| v()).collect();
        let a = if t % 7 == 0 {
            0.0
        } else {
            (rng.next_f64() as f32 * 2.0 - 1.0) * 3.0
        };
        f(n, ix, iy, a, x, y);
    }
}
fn offset(n: i32, inc: i32) -> usize {
    if n > 0 && inc < 0 {
        ((1 - n) * inc) as usize
    } else {
        0
    }
}
fn same(a: f32, b: f32) {
    assert!(
        a.to_bits() == b.to_bits()
            || (a.is_nan() && b.is_nan())
            || (a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-5 * a.abs().max(b.abs())),
        "{a:?} != {b:?}"
    );
}
fn array(a: &[f32], b: &[f32]) {
    for (&a, &b) in a.iter().zip(b) {
        same(a, b);
    }
}

#[test]
fn isamax_matches() {
    cases(|n, ix, _, _, x, _| unsafe {
        assert_eq!(
            rust::isamax(n, &x, ix),
            isamax(n, x.as_ptr().add(offset(n, ix)), ix)
        );
    });
}
#[test]
fn sasum_matches() {
    cases(|n, ix, _, _, x, _| unsafe {
        same(
            rust::sasum(n, &x, ix),
            sasum(n, x.as_ptr().add(offset(n, ix)), ix),
        );
    });
}
#[test]
fn snrm2_matches() {
    cases(|n, ix, _, _, x, _| unsafe {
        same(
            rust::snrm2(n, &x, ix),
            snrm2(n, x.as_ptr().add(offset(n, ix)), ix),
        );
    });
}
#[test]
fn sscal_matches() {
    cases(|n, ix, _, a, mut x, _| unsafe {
        let mut c = x.clone();
        sscal(n, a, c.as_mut_ptr().add(offset(n, ix)), ix);
        rust::sscal(n, a, &mut x, ix);
        array(&x, &c);
    });
}
#[test]
fn saxpy_matches() {
    cases(|n, ix, iy, a, x, mut y| unsafe {
        let mut c = y.clone();
        saxpy(n, a, x.as_ptr(), ix, c.as_mut_ptr(), iy);
        rust::saxpy(n, a, &x, ix, &mut y, iy);
        array(&y, &c);
    });
}
#[test]
fn scopy_matches() {
    cases(|n, ix, iy, _, x, mut y| unsafe {
        let mut c = y.clone();
        scopy(n, x.as_ptr(), ix, c.as_mut_ptr(), iy);
        rust::scopy(n, &x, ix, &mut y, iy);
        for (a, b) in y.iter().zip(c) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    });
}
#[test]
fn sdot_matches() {
    cases(|n, ix, iy, _, x, y| unsafe {
        same(
            rust::sdot(n, &x, ix, &y, iy),
            sdot(n, x.as_ptr(), ix, y.as_ptr(), iy),
        );
    });
}
#[test]
fn sswap_matches() {
    cases(|n, ix, iy, _, mut x, mut y| unsafe {
        let mut cx = x.clone();
        let mut cy = y.clone();
        sswap(n, cx.as_mut_ptr(), ix, cy.as_mut_ptr(), iy);
        rust::sswap(n, &mut x, ix, &mut y, iy);
        for (a, b) in x.iter().chain(&y).zip(cx.iter().chain(&cy)) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    });
}
