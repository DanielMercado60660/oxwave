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

//! Differential tests for cwp/lib/complexfd.c.
//! Credits: Dave Hale, original version in C++; Gabriel Alvarez, translation to C.

use oxwave_cref::XorShift;
use oxwave_f_complexfd as rust;
use oxwave_f_complexfd::DComplex;

mod c {
    use super::DComplex;
    use std::ffi::{c_char, c_void, CString};
    use std::process::Command;
    use std::sync::OnceLock;

    // libcwp.a in the supplied SU build omits complexfd.o. Declare its
    // original prototypes here and load the same source compiled against that
    // archive, so these tests remain runnable without changing shared crates.
    #[allow(dead_code)]
    extern "C" {
        fn dccos(a: DComplex) -> DComplex;
        fn dccosh(a: DComplex) -> DComplex;
        fn dcdcpow(a: DComplex, p: DComplex) -> DComplex;
        fn dcexp1(a: DComplex) -> DComplex;
        fn dcipow(a: DComplex, p: i32) -> DComplex;
        fn dclog(a: DComplex) -> DComplex;
        fn dcrpow(a: DComplex, p: f32) -> DComplex;
        fn dcsin(a: DComplex) -> DComplex;
        fn dcsinh(a: DComplex) -> DComplex;
        fn rdcpow(a: f32, p: DComplex) -> DComplex;
        fn dlopen(path: *const c_char, mode: i32) -> *mut c_void;
        fn dlsym(handle: *mut c_void, name: *const c_char) -> *mut c_void;
    }

    fn symbol(name: &str) -> *mut c_void {
        static HANDLE: OnceLock<usize> = OnceLock::new();
        let handle = *HANDLE.get_or_init(|| {
            let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
            let executable = std::env::current_exe().expect("test executable path");
            let lib = executable
                .parent()
                .unwrap()
                .join(format!("complexfd-ref-{}.dylib", std::process::id()));
            let status = Command::new("/usr/bin/cc")
                .args([
                    "-dynamiclib",
                    "-std=c90",
                    "-m64",
                    "-ansi",
                    "-ffp-contract=off",
                    "-DCWP_LITTLE_ENDIAN",
                    "-DSUXDR",
                ])
                .arg(format!(
                    "-I{}",
                    root.join("vendor/su-build/include").display()
                ))
                .arg(root.join("vendor/su-original/src/cwp/lib/complexfd.c"))
                .arg(root.join("vendor/su-build/lib/libcwp.a"))
                .arg("-o")
                .arg(&lib)
                .status()
                .expect("compile original complexfd.c");
            assert!(status.success(), "compiling original complexfd.c failed");
            let path = CString::new(lib.to_str().unwrap()).unwrap();
            let handle = unsafe { dlopen(path.as_ptr(), 2) };
            assert!(!handle.is_null(), "loading original complexfd.c failed");
            handle as usize
        });
        let name = CString::new(name).unwrap();
        let ptr = unsafe { dlsym(handle as *mut c_void, name.as_ptr()) };
        assert!(!ptr.is_null());
        ptr
    }

    macro_rules! wrapper {
        ($rust:ident, $name:literal, ($($arg:ident : $ty:ty),*)) => {
            pub unsafe fn $rust($($arg: $ty),*) -> DComplex {
                let f: unsafe extern "C" fn($($ty),*) -> DComplex = std::mem::transmute(symbol($name));
                f($($arg),*)
            }
        };
    }
    wrapper!(call_dccos, "dccos", (a: DComplex));
    wrapper!(call_dccosh, "dccosh", (a: DComplex));
    wrapper!(call_dcdcpow, "dcdcpow", (a: DComplex, p: DComplex));
    wrapper!(call_dcexp1, "dcexp1", (a: DComplex));
    wrapper!(call_dcipow, "dcipow", (a: DComplex, p: i32));
    wrapper!(call_dclog, "dclog", (a: DComplex));
    wrapper!(call_dcrpow, "dcrpow", (a: DComplex, p: f32));
    wrapper!(call_dcsin, "dcsin", (a: DComplex));
    wrapper!(call_dcsinh, "dcsinh", (a: DComplex));
    wrapper!(call_rdcpow, "rdcpow", (a: f32, p: DComplex));
}

fn z(r: f64, i: f64) -> DComplex {
    DComplex { r, i }
}

fn compare(a: DComplex, b: DComplex, context: &str) {
    // Compare nonfinite components by classification, without allowing an
    // infinite component to give the other component an infinite tolerance.
    let scale = [a.r, a.i, b.r, b.i]
        .into_iter()
        .filter(|x| x.is_finite())
        .fold(0.0f64, |scale, x| scale.max(x.abs()));
    for (x, y) in [(a.r, b.r), (a.i, b.i)] {
        if x == 0.0 && y == 0.0 {
            assert_eq!(x.to_bits(), y.to_bits(), "{context}: signed zero");
        }
        assert!(
            x == y
                || (x.is_nan() && y.is_nan())
                || (x.is_finite() && y.is_finite() && (x - y).abs() <= 1e-5 * scale),
            "{context}: Rust {a:?}, C {b:?}"
        );
    }
}

fn inputs() -> Vec<DComplex> {
    let mut rng = XorShift::new(0x142928);
    let mut xs = Vec::new();
    // Scalar functions: these sizes exercise runs of varied lengths.
    for n in [1, 2, 7, 16, 101, 1024, 8191] {
        for _ in 0..n {
            let scale = match rng.next_u32() % 5 {
                0 => 1e-10,
                1 => 0.01,
                2 => 1.0,
                3 => 10.0,
                _ => 80.0,
            };
            xs.push(z(rng.uniform(-scale, scale), rng.uniform(-scale, scale)));
        }
    }
    for v in [
        0.0,
        -0.0,
        1.0,
        -1.0,
        // Rounding ties when the original narrows double components to float.
        1.0 + 2.0f64.powi(-24),
        1.0 + 3.0 * 2.0f64.powi(-24),
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        f64::MAX,
        -f64::MAX,
        f32::from_bits(1) as f64,
        f32::MIN_POSITIVE as f64,
        f32::MAX as f64,
        // Float squares underflow near 1e-23 and overflow near 1e19.
        1e-24,
        1e-22,
        1e-40,
        -1e-40,
        1e19,
        2e19,
        1e20,
        -1e20,
        // Float exp amplitude underflows/overflows near -104 and 89.
        -104.0,
        -103.0,
        88.0,
        89.0,
        1000.0,
        -1000.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ] {
        for w in [0.0, -0.0, 1.0, -1.0, v] {
            xs.push(z(v, w));
        }
    }
    // Exercise independent exponents/signs across the full IEEE range,
    // including subnormals and nonfinite values, without filtering any cases.
    for _ in 0..4096 {
        xs.push(z(
            f64::from_bits(rng.next_u64()),
            f64::from_bits(rng.next_u64()),
        ));
        xs.push(z(rng.any_f32() as f64, rng.any_f32() as f64));
    }
    xs
}

#[test]
fn dccos_matches_c() {
    for a in inputs() {
        compare(
            rust::dccos(a),
            unsafe { c::call_dccos(a) },
            &format!("dccos({a:?})"),
        );
    }
}
#[test]
fn dccosh_matches_c() {
    for a in inputs() {
        compare(
            rust::dccosh(a),
            unsafe { c::call_dccosh(a) },
            &format!("dccosh({a:?})"),
        );
    }
}
#[test]
fn dcsin_matches_c() {
    for a in inputs() {
        compare(
            rust::dcsin(a),
            unsafe { c::call_dcsin(a) },
            &format!("dcsin({a:?})"),
        );
    }
}
#[test]
fn dcsinh_matches_c() {
    for a in inputs() {
        compare(
            rust::dcsinh(a),
            unsafe { c::call_dcsinh(a) },
            &format!("dcsinh({a:?})"),
        );
    }
}
#[test]
fn dcexp1_matches_c() {
    for a in inputs() {
        compare(
            rust::dcexp1(a),
            unsafe { c::call_dcexp1(a) },
            &format!("dcexp1({a:?})"),
        );
    }
}
#[test]
fn dclog_matches_c() {
    for a in inputs() {
        compare(
            rust::dclog(a),
            unsafe { c::call_dclog(a) },
            &format!("dclog({a:?})"),
        );
    }
}

#[test]
fn dcipow_matches_c() {
    let mut rng = XorShift::new(7);
    // Zero bases return before C negates p, so INT_MIN is defined here.
    // Check all four signed-zero combinations without entering its undefined
    // nonzero-base branch.
    for r in [0.0, -0.0] {
        for i in [0.0, -0.0] {
            let a = z(r, i);
            compare(
                rust::dcipow(a, i32::MIN),
                unsafe { c::call_dcipow(a, i32::MIN) },
                &format!("dcipow({a:?}, INT_MIN)"),
            );
        }
    }
    for a in inputs() {
        // INT_MIN has undefined negation in C and is outside its defined
        // nonzero-base domain. The adjacent value exercises all 31 bits.
        for p in [
            0,
            1,
            -1,
            2,
            -2,
            i32::MIN + 1,
            i32::MAX,
            rng.range_i32(-12, 12),
            rng.range_i32(-1024, 1024),
            (rng.next_u32() & 0x7fff_ffff) as i32,
            -((rng.next_u32() & 0x7fff_ffff) as i32),
        ] {
            let got = rust::dcipow(a, p);
            let expected = unsafe { c::call_dcipow(a, p) };
            for (x, y) in [(got.r, expected.r), (got.i, expected.i)] {
                assert!(
                    x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()),
                    "dcipow({a:?}, {p}): Rust {got:?}, C {expected:?}"
                );
            }
        }
    }
}

#[test]
fn dcrpow_matches_c() {
    let mut rng = XorShift::new(8);
    for a in inputs() {
        for p in [
            0.0,
            -0.0,
            -1.0,
            0.5,
            rng.uniform_f32(-4.0, 4.0),
            rng.any_f32(),
        ] {
            compare(
                rust::dcrpow(a, p),
                unsafe { c::call_dcrpow(a, p) },
                &format!("dcrpow({a:?}, {p})"),
            );
        }
    }
}

#[test]
fn dcdcpow_matches_c() {
    let mut rng = XorShift::new(9);
    for a in inputs() {
        for p in [
            z(0.0, 0.0),
            z(-1.0, 0.0),
            z(0.0, -1.0),
            z(rng.uniform(-4.0, 4.0), rng.uniform(-4.0, 4.0)),
            z(
                f64::from_bits(rng.next_u64()),
                f64::from_bits(rng.next_u64()),
            ),
        ] {
            compare(
                rust::dcdcpow(a, p),
                unsafe { c::call_dcdcpow(a, p) },
                &format!("dcdcpow({a:?}, {p:?})"),
            );
        }
    }
}

#[test]
fn rdcpow_matches_c() {
    let mut rng = XorShift::new(10);
    for p in inputs() {
        for a in [
            0.0,
            -0.0,
            1.0,
            -1.0,
            rng.uniform_f32(-10.0, 10.0),
            rng.any_f32(),
        ] {
            compare(
                rust::rdcpow(a, p),
                unsafe { c::call_rdcpow(a, p) },
                &format!("rdcpow({a}, {p:?})"),
            );
        }
    }
}
