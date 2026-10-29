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

//! Differential tests for cwp/lib/sort.c, credited to Dave Hale.
use oxwave_cref::XorShift;
use std::sync::Mutex;

static LOCK: Mutex<()> = Mutex::new(());
extern "C" {
    fn hpsort(n: i32, a: *mut f32);
    fn qksort(n: i32, a: *mut f32);
    fn qkfind(m: i32, n: i32, a: *mut f32);
    fn qkisort(n: i32, a: *mut f32, i: *mut i32);
    fn qkifind(m: i32, n: i32, a: *mut f32, i: *mut i32);
}
fn input(r: &mut XorShift, case: usize) -> Vec<f32> {
    let sizes = [1, 2, 3, 7, 8, 9, 11, 17, 31, 32, 63, 64, 127, 128, 257, 512, 1024];
    let n = sizes[case % sizes.len()];
    (0..n).map(|k| {
        let x = r.next_u32();
        match case % 8 {
            0 => 0.0,
            1 => if x & 1 == 0 { 0.0 } else { -0.0 },
            2 => (x % 7) as f32 - 3.0,
            3 => (x as i32 as f32) * 1e-20,
            4 => (x as i32 as f32) * 1e10,
            5 => k as f32 - (n / 2) as f32,
            6 => (n - k) as f32,
            _ => (x as i32 as f32) / 1000.0,
        }
    }).collect()
}
fn bits(a: &[f32]) -> Vec<u32> { a.iter().map(|x| x.to_bits()).collect() }
fn indices(n: usize, r: &mut XorShift) -> Vec<i32> {
    let mut v: Vec<i32> = (0..n as i32).collect();
    for k in (1..n).rev() { v.swap(k, r.next_u32() as usize % (k + 1)); }
    v
}
#[test]
fn hpsort_matches_c() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(1);
    for case in 0..3000 {
        let mut rust = input(&mut r, case);
        let mut c = rust.clone();
        unsafe { hpsort(c.len() as i32, c.as_mut_ptr()); }
        oxwave_sort::hpsort(&mut rust);
        assert_eq!(bits(&rust), bits(&c), "case {case}");
    }
}
#[test]
fn qksort_matches_c() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(2);
    for case in 0..3000 {
        let mut rust = input(&mut r, case);
        let mut c = rust.clone();
        unsafe { qksort(c.len() as i32, c.as_mut_ptr()); }
        oxwave_sort::qksort(&mut rust);
        assert_eq!(bits(&rust), bits(&c), "case {case}");
    }
}
#[test]
fn qkfind_matches_c() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(3);
    for case in 0..3000 {
        let mut rust = input(&mut r, case);
        let mut c = rust.clone();
        let m = match case % 3 { 0 => 0, 1 => c.len() - 1, _ => r.next_u32() as usize % c.len() };
        unsafe { qkfind(m as i32, c.len() as i32, c.as_mut_ptr()); }
        oxwave_sort::qkfind(m, &mut rust);
        assert_eq!(bits(&rust), bits(&c), "case {case}, m {m}");
    }
}
#[test]
fn qkisort_matches_c() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(4);
    for case in 0..3000 {
        let mut a = input(&mut r, case);
        let mut rust = indices(a.len(), &mut r);
        let mut c = rust.clone();
        unsafe { qkisort(c.len() as i32, a.as_mut_ptr(), c.as_mut_ptr()); }
        oxwave_sort::qkisort(&a, &mut rust);
        assert_eq!(rust, c, "case {case}");
    }
}
#[test]
fn qkifind_matches_c() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(5);
    for case in 0..3000 {
        let mut a = input(&mut r, case);
        let mut rust = indices(a.len(), &mut r);
        let mut c = rust.clone();
        let m = match case % 3 { 0 => 0, 1 => c.len() - 1, _ => r.next_u32() as usize % c.len() };
        unsafe { qkifind(m as i32, c.len() as i32, a.as_mut_ptr(), c.as_mut_ptr()); }
        oxwave_sort::qkifind(m, &a, &mut rust);
        assert_eq!(rust, c, "case {case}, m {m}");
    }
}

