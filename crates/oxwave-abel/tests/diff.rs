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

//! Differential tests for cwp/lib/abel.c.
//! Credits: Dave Hale and Lydia Deng, Colorado School of Mines, 06/01/90.
use oxwave_cref::XorShift;
use std::ffi::c_void;
extern "C" {
    fn abelalloc(n: i32) -> *mut c_void;
    fn abelfree(at: *mut c_void);
    fn abel(at: *mut c_void, f: *mut f32, g: *mut f32);
}
fn check(n: usize, input: &[f32]) {
    let at = oxwave_abel::abelalloc(n as i32);
    let cat = unsafe { abelalloc(n as i32) };
    assert!(!cat.is_null());
    let mut f = input.to_vec();
    let mut cg = vec![0.0; n];
    unsafe { abel(cat, f.as_mut_ptr(), cg.as_mut_ptr()) };
    let mut rg = vec![0.0; n];
    oxwave_abel::abel(&at, input, &mut rg);
    let scale = cg
        .iter()
        .chain(rg.iter())
        .fold(0.0f32, |m, x| m.max(x.abs()));
    for i in 0..n {
        assert!(
            (cg[i] - rg[i]).abs() <= 1e-5 * scale,
            "n={n} i={i} C={} Rust={} scale={scale}",
            cg[i],
            rg[i]
        );
    }
    let mut in_place = input.to_vec();
    unsafe { abel(cat, in_place.as_mut_ptr(), in_place.as_mut_ptr()) };
    assert_eq!(cg, in_place);
    unsafe { abelfree(cat) };
    oxwave_abel::abelfree(at);
}
#[test]
fn abel_matches_c() {
    let mut r = XorShift::new(519);
    let sizes = [2, 3, 5, 7, 16, 31, 64, 127, 256, 1024];
    for case in 0..3000 {
        let n = sizes[case % sizes.len()];
        let input: Vec<f32> = (0..n)
            .map(|i| match case % 8 {
                0 => 0.0,
                1 => -1.0,
                2 => {
                    if i == 0 {
                        1.0
                    } else {
                        0.0
                    }
                }
                3 => r.uniform_f32(-1e-10, 1e-10),
                4 => r.uniform_f32(-1e10, 1e10),
                _ => r.uniform_f32(-100.0, 100.0),
            })
            .collect();
        check(n, &input);
    }
}
#[test]
fn abelalloc_matches_c() {
    let mut r = XorShift::new(72);
    for _ in 0..3000 {
        let n = 1 + (r.next_u32() % 257) as i32;
        let cat = unsafe { abelalloc(n) };
        assert!(!cat.is_null());
        let at = oxwave_abel::abelalloc(n);
        if n > 1 {
            let input = vec![1.0; n as usize];
            let mut cg = vec![0.0; n as usize];
            let mut rg = cg.clone();
            let mut f = input.clone();
            unsafe { abel(cat, f.as_mut_ptr(), cg.as_mut_ptr()) };
            oxwave_abel::abel(&at, &input, &mut rg);
            let scale = cg.iter().fold(0.0f32, |m, x| m.max(x.abs()));
            for (c, rust) in cg.iter().zip(rg.iter()) {
                assert!((c - rust).abs() <= 1e-5 * scale);
            }
        }
        unsafe { abelfree(cat) };
        oxwave_abel::abelfree(at);
    }
}
#[test]
fn abelfree_matches_c() {
    let mut r = XorShift::new(35);
    for _ in 0..3000 {
        let n = 1 + (r.next_u32() % 257) as i32;
        let cat = unsafe { abelalloc(n) };
        let at = oxwave_abel::abelalloc(n);
        unsafe { abelfree(cat) };
        oxwave_abel::abelfree(at);
    }
}
