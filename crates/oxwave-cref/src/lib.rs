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

//! Test-only bindings to the ORIGINAL SU C libraries (libsu, libpar,
//! libcwp from vendor/su-build/lib), for differential tests of the Rust
//! ports: call the C function and the Rust port on the same random inputs
//! and compare. Never a dependency of a program; use it from
//! `[dev-dependencies]` only.
//!
//! Adding a binding: declare it in the `extern "C"` block below (or in your
//! own test file) with the exact C prototype. `complex` is `Complex` from
//! oxwave-cwp (or any `#[repr(C)] struct { r: f32, i: f32 }`); `float*`
//! is `*mut f32`; `int` is `i32`. Functions that are `static` inside SU
//! programs are not in any library: copy them verbatim into
//! `csrc/shim.c` with an `oxcref_` prefix.
//!
//! Note: libpar's `err()` calls `exit(1)`, which kills the test process;
//! only call C functions on inputs where they succeed.

use std::ffi::{c_char, c_int, c_uint};

/// `#[repr(C)]` single-precision complex, layout-identical to C `complex`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CComplex {
    pub r: f32,
    pub i: f32,
}

extern "C" {
    // ---- cwp/lib/pfafft.c
    pub fn npfa(nmin: c_int) -> c_int;
    pub fn npfao(nmin: c_int, nmax: c_int) -> c_int;
    pub fn npfar(nmin: c_int) -> c_int;
    pub fn npfaro(nmin: c_int, nmax: c_int) -> c_int;
    pub fn pfacc(isign: c_int, n: c_int, z: *mut CComplex);
    pub fn pfarc(isign: c_int, n: c_int, rz: *mut f32, cz: *mut CComplex);
    pub fn pfacr(isign: c_int, n: c_int, cz: *mut CComplex, rz: *mut f32);
    // ---- cwp/lib/franuni.c, frannor.c
    pub fn franuni() -> f32;
    pub fn sranuni(seed: c_int);
    pub fn frannor() -> f32;
    pub fn srannor(seed: c_int);
    // ---- cwp/lib/complex.c
    pub fn cadd(a: CComplex, b: CComplex) -> CComplex;
    pub fn csub(a: CComplex, b: CComplex) -> CComplex;
    pub fn cmul(a: CComplex, b: CComplex) -> CComplex;
    pub fn cdiv(a: CComplex, b: CComplex) -> CComplex;
    pub fn rcabs(z: CComplex) -> f32;
    pub fn cwp_csqrt(z: CComplex) -> CComplex;
    pub fn cwp_cexp(z: CComplex) -> CComplex;
    // ---- par/lib/getpars.c, atopkge.c
    pub fn initargs(argc: c_int, argv: *mut *mut c_char);
    pub fn getparint(name: *const c_char, p: *mut c_int) -> c_int;
    pub fn getparfloat(name: *const c_char, p: *mut f32) -> c_int;
    pub fn getpardouble(name: *const c_char, p: *mut f64) -> c_int;
    pub fn getparstring(name: *const c_char, p: *mut *mut c_char) -> c_int;
    pub fn getnparint(n: c_int, name: *const c_char, p: *mut c_int) -> c_int;
    pub fn countparval(name: *const c_char) -> c_int;
    pub fn countparname(name: *const c_char) -> c_int;
    pub fn eatoi(s: *const c_char) -> c_int;
    pub fn eatof(s: *const c_char) -> f32;
    pub fn eatod(s: *const c_char) -> f64;
    pub fn eatop(s: *const c_char) -> c_uint;
    // ---- su/lib/hdrpkge.c
    pub fn getindex(key: *const c_char) -> c_int;
    pub fn hdtype(key: *const c_char) -> *const c_char;
    // ---- cwp/lib/zebc.c
    pub fn zebc(ainput: *const c_char, aoutput: *mut c_char, nchar: c_int) -> c_int;
    // ---- csrc/shim.c (static functions copied from SU programs)
    pub fn oxcref_ibm_to_float(from: *mut c_int, to: *mut c_int, n: c_int, endian: c_int, verbose: c_int);
    pub fn oxcref_float_to_ibm(from: *mut c_int, to: *mut c_int, n: c_int, endian: c_int);
}

/// Deterministic xorshift64 generator for differential tests.
#[derive(Clone, Debug)]
pub struct XorShift(pub u64);

impl XorShift {
    /// Seeded generator (seed 0 is remapped to a fixed nonzero value).
    pub fn new(seed: u64) -> XorShift {
        XorShift(if seed == 0 { 0x9E37_79B9_7F4A_7C15 } else { seed })
    }
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    pub fn next_u32(&mut self) -> u32 {
        (self.next_u64() >> 32) as u32
    }
    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Uniform in [lo, hi).
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
    /// Uniform f32 in [lo, hi).
    pub fn uniform_f32(&mut self, lo: f32, hi: f32) -> f32 {
        self.uniform(lo as f64, hi as f64) as f32
    }
    /// Integer in [lo, hi] inclusive.
    pub fn range_i32(&mut self, lo: i32, hi: i32) -> i32 {
        lo + (self.next_u64() % ((hi as i64 - lo as i64 + 1) as u64)) as i32
    }
    /// Random f32 bit pattern (includes NaN, inf, denormals).
    pub fn any_f32(&mut self) -> f32 {
        f32::from_bits(self.next_u32())
    }
    /// Vector of uniform f32 values.
    pub fn vec_f32(&mut self, n: usize, lo: f32, hi: f32) -> Vec<f32> {
        (0..n).map(|_| self.uniform_f32(lo, hi)).collect()
    }
}

/// Relative difference helper: |a-b| / max(|a|,|b|,tiny).
pub fn rel_diff(a: f64, b: f64) -> f64 {
    if a == b || (a.is_nan() && b.is_nan()) {
        return 0.0;
    }
    (a - b).abs() / a.abs().max(b.abs()).max(1e-300)
}

/// Owned C argv for calling `initargs` (keeps the strings alive).
pub struct CArgv {
    _strings: Vec<std::ffi::CString>,
    ptrs: Vec<*mut c_char>,
}

impl CArgv {
    pub fn new(args: &[&str]) -> CArgv {
        let strings: Vec<std::ffi::CString> =
            args.iter().map(|s| std::ffi::CString::new(*s).unwrap()).collect();
        let mut ptrs: Vec<*mut c_char> = strings.iter().map(|s| s.as_ptr() as *mut c_char).collect();
        ptrs.push(std::ptr::null_mut());
        CArgv { _strings: strings, ptrs }
    }
    pub fn argc(&self) -> c_int {
        (self.ptrs.len() - 1) as c_int
    }
    pub fn argv(&mut self) -> *mut *mut c_char {
        self.ptrs.as_mut_ptr()
    }
}
