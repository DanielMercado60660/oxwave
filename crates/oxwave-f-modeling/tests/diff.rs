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

//! Differential tests for par/lib/modeling.c; credits: CWP Dave Hale, 09/17/91.
use oxwave_cref::XorShift;
use oxwave_f_modeling as rust;
use std::ffi::{c_char, c_void, CString};
static PAR_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
#[repr(C)]
struct Wavelet {
    lw: i32,
    iw: i32,
    wv: *mut f32,
}
#[repr(C)]
struct Segment {
    x: f32,
    z: f32,
    s: f32,
    c: f32,
}
#[repr(C)]
struct Reflector {
    ns: i32,
    ds: f32,
    a: f32,
    rs: *mut Segment,
}
extern "C" {
    fn addsinc(t: f32, a: f32, n: i32, dt: f32, ft: f32, tr: *mut f32);
    fn makericker(f: f32, dt: f32, w: *mut *mut Wavelet);
    fn raylv2(
        v: f32,
        dx: f32,
        dz: f32,
        x0: f32,
        z0: f32,
        x: f32,
        z: f32,
        c: *mut f32,
        s: *mut f32,
        t: *mut f32,
        q: *mut f32,
    );
    fn decodeReflector(
        s: *mut c_char,
        a: *mut f32,
        n: *mut i32,
        x: *mut *mut f32,
        z: *mut *mut f32,
    ) -> i32;
    fn decodeReflectors(
        n: *mut i32,
        a: *mut *mut f32,
        nu: *mut *mut i32,
        x: *mut *mut *mut f32,
        z: *mut *mut *mut f32,
    );
    fn breakReflectors(
        n: *mut i32,
        a: *mut *mut f32,
        nu: *mut *mut i32,
        x: *mut *mut *mut f32,
        z: *mut *mut *mut f32,
    );
    fn makeref(
        ds: f32,
        n: i32,
        a: *mut f32,
        nu: *mut i32,
        x: *mut *mut f32,
        z: *mut *mut f32,
        r: *mut *mut Reflector,
    );
    fn free(p: *mut c_void);
}
fn close(a: f32, b: f32, scale: f32) {
    assert!(
        (a - b).abs() <= 1e-5 * scale.max(1e-30),
        "{a} != {b}, scale {scale}"
    )
}
#[test]
fn sinc() {
    let mut rng = XorShift::new(1);
    for i in 0..4000 {
        let n = [1, 2, 7, 8, 31, 64, 257][i % 7];
        let dt = 0.004;
        let time = (rng.next_f64() as f32) * (n + 5) as f32 * dt;
        let amp = rng.next_f64() as f32 * 4.0 - 2.0;
        let mut a = vec![0.0; n];
        let mut b = a.clone();
        unsafe { addsinc(time, amp, n as i32, dt, 0.0, a.as_mut_ptr()) };
        rust::addsinc(time, amp, dt, 0.0, &mut b);
        for (x, y) in a.iter().zip(b) {
            close(*x, y, amp.abs())
        }
    }
}
#[test]
fn ricker() {
    let mut rng = XorShift::new(2);
    for _ in 0..3000 {
        let f = 5.0 + rng.next_f64() as f32 * 95.0;
        let dt = 0.001 + rng.next_f64() as f32 * 0.009;
        let mut p = std::ptr::null_mut();
        unsafe {
            makericker(f, dt, &mut p);
            let w = rust::makericker(f, dt);
            assert_eq!((w.lw, w.iw), ((*p).lw, (*p).iw));
            for (x, y) in
                w.wv.iter()
                    .zip(std::slice::from_raw_parts((*p).wv, (*p).lw as usize))
            {
                close(*x, *y, 1.0)
            }
            free((*p).wv.cast());
            free(p.cast())
        }
    }
}
#[test]
fn ray() {
    let mut rng = XorShift::new(3);
    for i in 0..5000 {
        let dx = if i % 3 == 0 {
            0.0
        } else {
            rng.next_f64() as f32 * 0.2
        };
        let dz = if i % 5 == 0 {
            0.0
        } else {
            rng.next_f64() as f32 * 0.2
        };
        let x = rng.next_f64() as f32 * 1000.0;
        let z = rng.next_f64() as f32 * 1000.0;
        let (mut c, mut s, mut t, mut q) = (0.0, 0.0, 0.0, 0.0);
        unsafe {
            raylv2(
                2000.0, dx, dz, 0.0, 0.0, x, z, &mut c, &mut s, &mut t, &mut q,
            )
        }
        let r = rust::raylv2(2000.0, dx, dz, 0.0, 0.0, x, z);
        for (a, b) in [(c, r.0), (s, r.1), (t, r.2), (q, r.3)] {
            close(a, b, a.abs().max(b.abs()))
        }
    }
}
#[test]
fn decode() {
    let _guard = PAR_LOCK.lock().unwrap();
    let mut rng = XorShift::new(4);
    for _ in 0..3000 {
        let spec = format!(
            "{}:{},{};{},{}",
            rng.next_f64(),
            rng.next_f64(),
            rng.next_f64(),
            rng.next_f64(),
            rng.next_f64()
        );
        let mut bytes = CString::new(spec.clone()).unwrap().into_bytes_with_nul();
        let (mut a, mut n, mut x, mut z) = (0.0, 0, std::ptr::null_mut(), std::ptr::null_mut());
        unsafe {
            assert_eq!(
                decodeReflector(bytes.as_mut_ptr().cast(), &mut a, &mut n, &mut x, &mut z),
                1
            );
            let r = rust::decodeReflector(&spec).unwrap();
            assert_eq!(r.a, a);
            assert_eq!(r.x, std::slice::from_raw_parts(x, n as usize));
            assert_eq!(r.z, std::slice::from_raw_parts(z, n as usize));
            free(x.cast());
            free(z.cast())
        }
    }
}
extern "C" {
    fn malloc(n: usize) -> *mut c_void;
}
unsafe fn copy_alloc<T: Copy>(x: &[T]) -> *mut T {
    let p = malloc(std::mem::size_of_val(x)).cast::<T>();
    std::ptr::copy_nonoverlapping(x.as_ptr(), p, x.len());
    p
}
unsafe fn input(r: &rust::ReflectorInput) -> (*mut f32, *mut i32, *mut *mut f32, *mut *mut f32) {
    let x = copy_alloc(&r.x);
    let z = copy_alloc(&r.z);
    (
        copy_alloc(&[r.a]),
        copy_alloc(&[r.x.len() as i32]),
        copy_alloc(&[x]),
        copy_alloc(&[z]),
    )
}
#[test]
fn break_reflectors() {
    let mut rng = XorShift::new(5);
    for i in 0..2000 {
        let n = [2, 3, 7, 16, 31][i % 5];
        let r = rust::ReflectorInput {
            a: rng.next_f64() as f32,
            x: (0..n).map(|_| rng.next_f64() as f32).collect(),
            z: (0..n).map(|_| rng.next_f64() as f32).collect(),
        };
        unsafe {
            let (mut a, mut nu, mut x, mut z) = input(&r);
            let mut nr = 1;
            breakReflectors(&mut nr, &mut a, &mut nu, &mut x, &mut z);
            let out = rust::breakReflectors(&[r]);
            assert_eq!(nr as usize, out.len());
            for j in 0..nr as usize {
                assert_eq!(*a.add(j), out[j].a);
                assert_eq!(*nu.add(j), 2);
                assert_eq!(std::slice::from_raw_parts(*x.add(j), 2), out[j].x);
                assert_eq!(std::slice::from_raw_parts(*z.add(j), 2), out[j].z);
                free((*x.add(j)).cast());
                free((*z.add(j)).cast())
            }
            free(a.cast());
            free(nu.cast());
            free(x.cast());
            free(z.cast())
        }
    }
}
#[test]
fn make_ref() {
    let mut rng = XorShift::new(6);
    for i in 0..2000 {
        let n = [2, 3, 7, 16][i % 4];
        let r = rust::ReflectorInput {
            a: rng.next_f64() as f32,
            x: (0..n).map(|j| j as f32 * 2.0).collect(),
            z: (0..n).map(|_| rng.next_f64() as f32 * 2.0).collect(),
        };
        unsafe {
            let (a, nu, x, z) = input(&r);
            let mut p = std::ptr::null_mut();
            makeref(1.0, 1, a, nu, x, z, &mut p);
            let out = rust::makeref(1.0, &[r]);
            assert_eq!((*p).ns, out[0].ns);
            close((*p).ds, out[0].ds, (*p).ds.abs());
            for j in 0..(*p).ns as usize {
                let c = &*(*p).rs.add(j);
                let v = &out[0].rs[j];
                for (a, b) in [(c.x, v.x), (c.z, v.z), (c.s, v.s), (c.c, v.c)] {
                    close(a, b, (n as f32 * 2.0).max(2.0))
                }
            }
            free((*p).rs.cast());
            free(p.cast());
            free(a.cast())
        }
    }
}
#[test]
fn decode_reflectors() {
    let _guard = PAR_LOCK.lock().unwrap();
    unsafe {
        let mut rng = XorShift::new(7);
        let mut args = vec!["test".to_string()];
        for _ in 0..2000 {
            args.push(format!(
                "ref={}:0,{};1,{}",
                rng.next_f64(),
                rng.next_f64(),
                rng.next_f64()
            ));
        }
        let cargs: Vec<_> = args
            .iter()
            .map(|s| CString::new(s.as_str()).unwrap())
            .collect();
        let mut argv: Vec<_> = cargs.iter().map(|s| s.as_ptr() as *mut c_char).collect();
        oxwave_cref::initargs(argv.len() as i32, argv.as_mut_ptr());
        oxwave_core::par::initargs(&args);
        let (mut n, mut a, mut nu, mut x, mut z) = (
            0,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        );
        decodeReflectors(&mut n, &mut a, &mut nu, &mut x, &mut z);
        let r = rust::decodeReflectors();
        assert_eq!(n as usize, r.len());
        for i in 0..n as usize {
            assert_eq!(*a.add(i), r[i].a);
            assert_eq!(
                std::slice::from_raw_parts(*x.add(i), *nu.add(i) as usize),
                r[i].x
            );
            assert_eq!(
                std::slice::from_raw_parts(*z.add(i), *nu.add(i) as usize),
                r[i].z
            );
            free((*x.add(i)).cast());
            free((*z.add(i)).cast());
        }
        free(a.cast());
        free(nu.cast());
        free(x.cast());
        free(z.cast())
    }
}
