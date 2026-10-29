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

//! Differential tests: oxwave-core / oxwave-cwp vs the original SU C code.

use oxwave_cref::*;
use std::ffi::{CStr, CString};
use std::sync::Mutex;

static C_LOCK: Mutex<()> = Mutex::new(());

extern "C" {
    /// macOS: address of the thread's errno.
    fn __error() -> *mut std::ffi::c_int;
}

#[test]
fn ibm_to_float_matches_c() {
    let _g = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(1);
    let n = 200_000;
    let mut words: Vec<i32> = (0..n).map(|_| r.next_u32() as i32).collect();
    words.extend_from_slice(&[0, 1, -1, i32::MIN, i32::MAX, 0x4110_0000, 0x0010_0000, 0x7f00_0000u32 as i32]);
    for endian in [0, 1] {
        let mut from = words.clone();
        let mut to = vec![0i32; words.len()];
        unsafe { oxcref_ibm_to_float(from.as_mut_ptr(), to.as_mut_ptr(), words.len() as i32, endian, 0) };
        let mut mine: Vec<f32> = words.iter().map(|&w| f32::from_bits(w as u32)).collect();
        oxwave_core::segy::ibm_to_float(&mut mine, words.len(), endian);
        for k in 0..words.len() {
            assert_eq!(mine[k].to_bits(), to[k] as u32, "word {:08x} endian {}", words[k], endian);
        }
    }
}

#[test]
fn float_to_ibm_matches_c() {
    let _g = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(2);
    let n = 200_000;
    let mut words: Vec<i32> = (0..n).map(|_| r.next_u32() as i32).collect();
    for x in [0.0f32, -0.0, 1.0, -1.0, f32::MAX, f32::MIN_POSITIVE, 1e-45, f32::INFINITY, f32::NAN] {
        words.push(x.to_bits() as i32);
    }
    for endian in [0, 1] {
        let mut from = words.clone();
        let mut to = vec![0i32; words.len()];
        unsafe { oxcref_float_to_ibm(from.as_mut_ptr(), to.as_mut_ptr(), words.len() as i32, endian) };
        let mut mine: Vec<f32> = words.iter().map(|&w| f32::from_bits(w as u32)).collect();
        oxwave_core::segy::float_to_ibm(&mut mine, words.len(), endian);
        for k in 0..words.len() {
            assert_eq!(mine[k].to_bits(), to[k] as u32, "word {:08x} endian {}", words[k], endian);
        }
    }
}

#[test]
fn zebc_matches_c() {
    let _g = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let input: Vec<u8> = (0..=255u8).collect();
    let mut out = vec![0u8; 256];
    unsafe { zebc(input.as_ptr() as *const _, out.as_mut_ptr() as *mut _, 256) };
    let mut mine = vec![0u8; 256];
    oxwave_core::segy::zebc(&input, &mut mine, 256);
    assert_eq!(mine, out);
}

#[test]
fn header_table_matches_c() {
    let _g = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    for (i, k) in oxwave_core::header::HDR.iter().enumerate() {
        let ck = CString::new(k.key).unwrap();
        assert_eq!(unsafe { getindex(ck.as_ptr()) }, i as i32);
        let t = unsafe { CStr::from_ptr(hdtype(ck.as_ptr())) }.to_str().unwrap();
        assert_eq!(t, k.type_);
    }
    let ck = CString::new("nokey").unwrap();
    assert_eq!(unsafe { getindex(ck.as_ptr()) }, -1);
}

fn c_getparint(name: &str) -> (i32, Vec<i32>) {
    let n = CString::new(name).unwrap();
    let mut buf = vec![-999i32; 32];
    let c = unsafe { getparint(n.as_ptr(), buf.as_mut_ptr()) };
    buf.truncate(c.max(0) as usize);
    (c, buf)
}
fn c_getparfloat(name: &str) -> (i32, Vec<f32>) {
    let n = CString::new(name).unwrap();
    let mut buf = vec![-999f32; 32];
    let c = unsafe { getparfloat(n.as_ptr(), buf.as_mut_ptr()) };
    buf.truncate(c.max(0) as usize);
    (c, buf)
}
fn c_getpardouble(name: &str) -> (i32, Vec<f64>) {
    let n = CString::new(name).unwrap();
    let mut buf = vec![-999f64; 32];
    let c = unsafe { getpardouble(n.as_ptr(), buf.as_mut_ptr()) };
    buf.truncate(c.max(0) as usize);
    (c, buf)
}
fn c_getparstring(name: &str) -> (i32, String) {
    let n = CString::new(name).unwrap();
    let mut p: *mut std::ffi::c_char = std::ptr::null_mut();
    let c = unsafe { getparstring(n.as_ptr(), &mut p) };
    let s = if c != 0 { unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned() } else { String::new() };
    (c, s)
}

#[test]
fn getpar_matches_c() {
    let _g = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let dir = std::env::temp_dir().join(format!("oxw_cref_par_{}", std::process::id()));
    std::fs::write(&dir, "a=1 b=2,3\n  c=\"x y\"  junk d=\t9\na=7 e='q r' f=\n").unwrap();
    let pf = format!("par={}", dir.display());
    let argsets: Vec<Vec<&str>> = vec![
        vec!["prog", "n=5", "v=1,,2", "w=3,", "x=", "f=1e-5,2.5x,abc,-0", "d=0x1p3,nan,-.5e2,4e-39", "s=a,b", "n=-7"],
        vec!["prog", "n= 12", "v=+4,  5", "p=4294967295", "n=1,2,3", "g=1e38,3.4e38"],
        vec!["prog", "a=3", &pf, "b=4"],
        vec!["prog", "par=", "a=1", "par=", "x=,"],
        vec!["prog"],
        vec!["prog", "noequals", "=5", "a==2", "a=b=c"],
    ];
    let names = ["n", "v", "w", "x", "f", "d", "s", "p", "g", "a", "b", "c", "e", "par", "", "q"];
    for set in &argsets {
        let mut cargv = CArgv::new(set);
        unsafe { oxwave_cref::initargs(cargv.argc(), cargv.argv()) };
        let rargv: Vec<String> = set.iter().map(|s| s.to_string()).collect();
        oxwave_core::par::initargs(&rargv);
        for name in names {
            let cn = CString::new(name).unwrap();
            assert_eq!(unsafe { countparval(cn.as_ptr()) }, oxwave_core::par::countparval(name), "countparval {} {:?}", name, set);
            assert_eq!(unsafe { countparname(cn.as_ptr()) }, oxwave_core::par::countparname(name), "countparname {} {:?}", name, set);
            if name == "p" || name == "g" {
                continue; // values overflow int / float on purpose
            }
            let (c, cv) = c_getparint(name);
            let mut rv: Vec<i32> = Vec::new();
            let r = oxwave_core::par::getparint(name, &mut rv);
            assert_eq!((c, &cv), (r, &rv), "getparint {} {:?}", name, set);
            let (c, cv) = c_getparfloat(name);
            let mut rv: Vec<f32> = Vec::new();
            let r = oxwave_core::par::getparfloat(name, &mut rv);
            assert_eq!(c, r, "getparfloat {} {:?}", name, set);
            for (x, y) in cv.iter().zip(rv.iter()) {
                assert!(x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()), "getparfloat {} {:?}", name, set);
            }
            let (c, cv) = c_getpardouble(name);
            let mut rv: Vec<f64> = Vec::new();
            let r = oxwave_core::par::getpardouble(name, &mut rv);
            assert_eq!(c, r);
            for (x, y) in cv.iter().zip(rv.iter()) {
                assert!(x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()), "getpardouble {} {:?}", name, set);
            }
            let (c, cs) = c_getparstring(name);
            let mut rs = String::new();
            let r = oxwave_core::par::getparstring(name, &mut rs);
            assert_eq!((c, cs), (r, rs), "getparstring {} {:?}", name, set);
            for k in 1..4 {
                let mut ci = -999;
                let c = unsafe { getnparint(k, cn.as_ptr(), &mut ci) };
                let mut ri = -999;
                let r = oxwave_core::par::getnparint(k, name, &mut ri);
                assert_eq!((c, ci), (r, ri), "getnparint {} {} {:?}", k, name, set);
            }
        }
    }
    std::fs::remove_file(&dir).unwrap();
}

#[test]
fn eatof_eatod_match_c_on_random_strings() {
    let _g = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(7);
    for _ in 0..20000 {
        let x = f64::from_bits(r.next_u64());
        if !x.is_finite() || x.abs() > 1e300 || (x != 0.0 && x.abs() < 1e-300) {
            continue;
        }
        let s = match r.next_u32() % 4 {
            0 => format!("{:e}", x),
            1 => format!("{}", x as f32),
            2 => format!("{:.3}", x / 1e300 * 1e3),
            _ => format!("{}x", r.uniform(-1e6, 1e6)),
        };
        let cs = CString::new(s.clone()).unwrap();
        oxwave_core::par::set_errno(0);
        let d = oxwave_core::par::strtod(&s).0;
        if oxwave_core::par::errno() != 0 || !d.is_finite() {
            continue; // C would err() and exit
        }
        unsafe { *__error() = 0 };
        let cd = unsafe { eatod(cs.as_ptr()) };
        assert_eq!(d.to_bits(), cd.to_bits(), "{}", s);
        if (d as f32).is_finite() {
            let f = oxwave_core::par::eatof(&s);
            unsafe { *__error() = 0 };
            let cf = unsafe { eatof(cs.as_ptr()) };
            assert_eq!(f.to_bits(), cf.to_bits(), "{}", s);
        }
    }
}

#[test]
fn complex_vs_c() {
    use oxwave_cwp::complex as cx;
    let _g = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut r = XorShift::new(3);
    let (mut wdiv, mut wabs, mut wsqrt, mut wexp) = (0f64, 0f64, 0f64, 0f64);
    let (mut ediv, mut eabs, mut esqrt, mut eexp) = (0u32, 0u32, 0u32, 0u32);
    let n = 100_000;
    for _ in 0..n {
        let a = CComplex { r: r.uniform_f32(-100.0, 100.0), i: r.uniform_f32(-100.0, 100.0) };
        let b = CComplex { r: r.uniform_f32(-100.0, 100.0), i: r.uniform_f32(-100.0, 100.0) };
        let ra = cx::cmplx(a.r, a.i);
        let rb = cx::cmplx(b.r, b.i);
        let c = unsafe { cadd(a, b) };
        assert_eq!((c.r, c.i), (cx::cadd(ra, rb).r, cx::cadd(ra, rb).i));
        let c = unsafe { csub(a, b) };
        assert_eq!((c.r, c.i), (cx::csub(ra, rb).r, cx::csub(ra, rb).i));
        let c = unsafe { cmul(a, b) };
        assert_eq!((c.r, c.i), (cx::cmul(ra, rb).r, cx::cmul(ra, rb).i));
        let c = unsafe { cdiv(a, b) };
        let m = cx::cdiv(ra, rb);
        let d = rel_diff(c.r as f64, m.r as f64).max(rel_diff(c.i as f64, m.i as f64));
        wdiv = wdiv.max(d);
        ediv += (c.r.to_bits() != m.r.to_bits() || c.i.to_bits() != m.i.to_bits()) as u32;
        let c = unsafe { rcabs(a) };
        let m = cx::rcabs(ra);
        wabs = wabs.max(rel_diff(c as f64, m as f64));
        eabs += (c.to_bits() != m.to_bits()) as u32;
        let c = unsafe { cwp_csqrt(a) };
        let m = cx::cwp_csqrt(ra);
        wsqrt = wsqrt.max(rel_diff(c.r as f64, m.r as f64).max(rel_diff(c.i as f64, m.i as f64)));
        esqrt += (c.r.to_bits() != m.r.to_bits() || c.i.to_bits() != m.i.to_bits()) as u32;
        let s = CComplex { r: a.r / 20.0, i: a.i };
        let c = unsafe { cwp_cexp(s) };
        let m = cx::cwp_cexp(cx::cmplx(s.r, s.i));
        wexp = wexp.max(rel_diff(c.r as f64, m.r as f64).max(rel_diff(c.i as f64, m.i as f64)));
        eexp += (c.r.to_bits() != m.r.to_bits() || c.i.to_bits() != m.i.to_bits()) as u32;
    }
    eprintln!(
        "complex vs C over {} cases: cdiv max rel {:.2e} ({} not bit-equal), rcabs {:.2e} ({}), csqrt {:.2e} ({}), cexp {:.2e} ({})",
        n, wdiv, ediv, wabs, eabs, wsqrt, esqrt, wexp, eexp
    );
    assert!(wdiv < 1e-5 && wabs < 1e-6 && wsqrt < 1e-5 && wexp < 1e-5);
}

#[test]
fn pfa_and_random_link() {
    let _g = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    unsafe {
        assert_eq!(npfa(1000), 1001); // 7*11*13
        let n = npfar(100);
        assert!(n >= 100);
        let mut z = vec![CComplex { r: 1.0, i: 0.0 }; 16];
        pfacc(1, 16, z.as_mut_ptr());
        assert!((z[0].r - 16.0).abs() < 1e-5 && z[1].r.abs() < 1e-5);
        sranuni(42);
        let a = franuni();
        sranuni(42);
        assert_eq!(a, franuni());
        assert!((0.0..1.0).contains(&a));
    }
}
