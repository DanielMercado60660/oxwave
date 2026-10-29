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

//! Differential tests for Trielas/lib/{findPiso,findSiso,findqSH,rottensh,rottensor,
//! rotvector,rt_SHa_real,rt_ani_real,rt_iso_cmplx,rt_iso_real,solveForSlowqP,
//! solveForSlowqSV}.c. Credits: Andreas Rueger and Tariq Alkalifah, CSM.
use oxwave_cref as _;
use oxwave_lib_trielas as rust;

#[allow(dead_code)]
mod c {
    #[link(name = "trielas", kind = "static")]
    unsafe extern "C" {
        pub fn solveForSlowqP(
            s: f32,
            c: f32,
            p: f32,
            a: f32,
            b: f32,
            f: f32,
            l: f32,
            m: f32,
            n: f32,
            x: *mut f32,
            z: *mut f32,
            rt: f32,
        ) -> i32;
        pub fn solveForSlowqSV(
            s: f32,
            c: f32,
            p: f32,
            a: f32,
            b: f32,
            f: f32,
            l: f32,
            m: f32,
            n: f32,
            x: *mut f32,
            z: *mut f32,
            rt: f32,
        ) -> i32;
        pub fn rt_iso_real(
            vp1: f32,
            vp2: f32,
            vs1: f32,
            vs2: f32,
            rho1: f32,
            rho2: f32,
            pl: f32,
            modei: i32,
            modet: i32,
            rort: i32,
            coeff: *mut f32,
        ) -> i32;
        pub fn rt_iso_cmplx(
            vp1: f32,
            vp2: f32,
            vs1: f32,
            vs2: f32,
            rho1: f32,
            rho2: f32,
            pl: f32,
            modei: i32,
            modet: i32,
            rort: i32,
            coeff: *mut f32,
            phase: *mut f32,
        ) -> i32;
        pub fn rt_ani_real(
            gz: f32,
            gx: f32,
            a1111i: f32,
            a3333i: f32,
            a1133i: f32,
            a1313i: f32,
            a1113i: f32,
            a3313i: f32,
            rhoi: f32,
            a1111t: f32,
            a3333t: f32,
            a1133t: f32,
            a1313t: f32,
            a1113t: f32,
            a3313t: f32,
            rhot: f32,
            pxi: f32,
            pzi: f32,
            g11i: f32,
            g13i: f32,
            g33i: f32,
            pxg: f32,
            pzg: f32,
            g11g: f32,
            g13g: f32,
            g33g: f32,
            modei: i32,
            modeg: i32,
            rort: i32,
            coeff: *mut f32,
            ifp: *mut std::ffi::c_void,
        ) -> i32;
        pub fn rotvector(x: *mut f32, z: *mut f32, s: f32, c: f32);
        pub fn rottensh(a: *mut f32, b: *mut f32, c: *mut f32, s: f32, k: f32);
        pub fn rottensor(
            a: *mut f32,
            b: *mut f32,
            c: *mut f32,
            d: *mut f32,
            e: *mut f32,
            f: *mut f32,
            s: f32,
            k: f32,
        );
        pub fn findPiso(
            a: f32,
            p: f32,
            x: f32,
            z: f32,
            px: *mut f32,
            pz: *mut f32,
            vx: *mut f32,
            vz: *mut f32,
            g11: *mut f32,
            g13: *mut f32,
            g33: *mut f32,
            r: i32,
        ) -> i32;
        pub fn findSiso(
            a: f32,
            p: f32,
            x: f32,
            z: f32,
            px: *mut f32,
            pz: *mut f32,
            vx: *mut f32,
            vz: *mut f32,
            g11: *mut f32,
            g13: *mut f32,
            g33: *mut f32,
            r: i32,
        ) -> i32;
        pub fn findqSH(
            s: f32,
            c: f32,
            p: f32,
            a: f32,
            b: f32,
            d: f32,
            px: *mut f32,
            pz: *mut f32,
            vx: *mut f32,
            vz: *mut f32,
            r: i32,
        ) -> i32;
        pub fn rt_SHa_real(
            a: f32,
            b: f32,
            c: f32,
            d: f32,
            e: f32,
            f: f32,
            g: f32,
            h: f32,
            i: f32,
            j: f32,
            k: f32,
            l: f32,
            r: i32,
            out: *mut f32,
            s: f32,
            co: f32,
        ) -> i32;
    }
}
use oxwave_cref::XorShift as Rng;
trait Samples {
    fn f(&mut self) -> f32;
}
impl Samples for Rng {
    fn f(&mut self) -> f32 {
        self.uniform_f32(-2.0, 2.0)
    }
}
// Compare each output independently, with no tolerance for finite values:
// a large phase must never mask even a one-ulp coefficient error. NaN payloads
// are unspecified, but zeros (including their signs) and infinities are exact.
#[track_caller]
fn exact(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        assert!(
            x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()),
            "element {i}: {x} ({:#010x}) != {y} ({:#010x})",
            x.to_bits(),
            y.to_bits()
        );
    }
}

fn sample(rng: &mut Rng, n: usize) -> f32 {
    match n % 19 {
        0 => 0.0,
        1 => -0.0,
        2 => 1.0,
        3 => -1.0,
        4 => rng.f() * 1e-12,
        5 => rng.f() * 1e12,
        _ => rng.f(),
    }
}

fn rotation(rng: &mut Rng, n: usize) -> (f32, f32) {
    match n % 13 {
        0 => (0.0, 1.0),
        1 => (0.0, -1.0),
        2 => (1.0, 0.0),
        3 => (-1.0, 0.0),
        4 => (0.0, 0.0),
        5 => (rng.f(), rng.f()), // The C API also accepts non-unit rotations.
        _ => {
            let angle = rng.uniform(-std::f64::consts::PI, std::f64::consts::PI);
            (angle.sin() as f32, angle.cos() as f32)
        }
    }
}
#[test]
fn rotvector() {
    let mut r = Rng::new(7);
    for n in 0..5000 {
        let (s, c) = rotation(&mut r, n);
        let (mut x, mut z) = (sample(&mut r, n), sample(&mut r, n + 7));
        let (mut a, mut b) = (x, z);
        unsafe { c::rotvector(&mut a, &mut b, s, c) };
        rust::rotvector(&mut x, &mut z, s, c);
        exact(&[x, z], &[a, b]);
    }
}
#[test]
fn rottensh() {
    let mut r = Rng::new(9);
    for n in 0..5000 {
        let (s, c) = rotation(&mut r, n);
        let (mut x, mut y, mut z) = (
            sample(&mut r, n),
            sample(&mut r, n + 7),
            sample(&mut r, n + 11),
        );
        let (mut a, mut b, mut d) = (x, y, z);
        unsafe { c::rottensh(&mut a, &mut b, &mut d, s, c) };
        rust::rottensh(&mut x, &mut y, &mut z, s, c);
        exact(&[x, y, z], &[a, b, d]);
    }
}
#[test]
fn rottensor() {
    let mut r = Rng::new(11);
    for n in 0..5000 {
        let (s, k) = rotation(&mut r, n);
        let (mut a, mut b, mut c, mut d, mut e, mut f) = (
            sample(&mut r, n),
            sample(&mut r, n + 1),
            sample(&mut r, n + 3),
            sample(&mut r, n + 7),
            sample(&mut r, n + 9),
            sample(&mut r, n + 11),
        );
        let (mut aa, mut bb, mut cc, mut dd, mut ee, mut ff) = (a, b, c, d, e, f);
        unsafe { c::rottensor(&mut aa, &mut bb, &mut cc, &mut dd, &mut ee, &mut ff, s, k) };
        rust::rottensor(&mut a, &mut b, &mut c, &mut d, &mut e, &mut f, s, k);
        exact(&[a, b, c, d, e, f], &[aa, bb, cc, dd, ee, ff]);
    }
}
fn iso(shear: bool) {
    let mut rng = Rng::new(13);
    for n in 0..5000 {
        let a = sample(&mut rng, n);
        let p = match n % 11 {
            0 => 0.0,
            1 => 0.5,
            2 => -0.5,
            3 => 1.0,
            4 => 1e6,
            _ => rng.f(),
        };
        let (x, z) = rotation(&mut rng, n);
        let mut q = [
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
        ];
        let mut c = q;
        let [px, pz, vx, vz, g11, g13, g33] = &mut q;
        let [cp, cz, cvx, cvz, c11, c13, c33] = &mut c;
        let r = (n % 2) as i32;
        let rc = unsafe {
            if shear {
                c::findSiso(a, p, x, z, cp, cz, cvx, cvz, c11, c13, c33, r)
            } else {
                c::findPiso(a, p, x, z, cp, cz, cvx, cvz, c11, c13, c33, r)
            }
        };
        let rr = if shear {
            rust::findSiso(a, p, x, z, px, pz, vx, vz, g11, g13, g33, r)
        } else {
            rust::findPiso(a, p, x, z, px, pz, vx, vz, g11, g13, g33, r)
        };
        assert_eq!(rc, rr);
        exact(&q, &c);
    }
}
#[test]
fn find_piso() {
    iso(false)
}
#[test]
fn find_siso() {
    iso(true)
}
#[test]
fn find_sh() {
    let mut rng = Rng::new(17);
    for n in 0..5000 {
        let (s, c) = rotation(&mut rng, n);
        let (p, a, b, d) = if n % 17 == 0 {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            (
                sample(&mut rng, n),
                sample(&mut rng, n + 3),
                sample(&mut rng, n + 7),
                sample(&mut rng, n + 11),
            )
        };
        let mut q = [123.0, -456.0, 789.0, -12.0];
        let mut cq = q;
        let [x, z, v, w] = &mut q;
        let [cx, cz, cv, cw] = &mut cq;
        let rc = unsafe { c::findqSH(s, c, p, a, b, d, cx, cz, cv, cw, (n % 2) as i32) };
        let rr = rust::findqSH(s, c, p, a, b, d, x, z, v, w, (n % 2) as i32);
        assert_eq!(rc, rr);
        exact(&q, &cq);
    }
}
#[test]
fn rt_sh() {
    let mut rng = Rng::new(19);
    for n in 0..5000 {
        let (a, b, c, d, e, f, g, h, i, j, k, l, s, co) = (
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
            rng.f(),
        );
        let (mut x, mut y) = (0.0, 0.0);
        let rc =
            unsafe { c::rt_SHa_real(a, b, c, d, e, f, g, h, i, j, k, l, n % 3, &mut x, s, co) };
        let rr = rust::rt_SHa_real(a, b, c, d, e, f, g, h, i, j, k, l, n % 3, &mut y, s, co);
        assert_eq!(rc, rr);
        exact(&[x], &[y]);
    }
    // Zero traction, near-singular denominators, equal media and normal incidence.
    for stiffness in [0.0, f32::EPSILON * 0.25, 1.0, -1.0, 1e12] {
        for pz in [0.0, -0.0, 1.0, -1.0] {
            for rort in 0..3 {
                let (mut x, mut y) = (123.0, 123.0);
                let rc = unsafe {
                    c::rt_SHa_real(
                        stiffness, stiffness, 0.0, stiffness, stiffness, 0.0, 0.0, pz, 0.0, pz,
                        0.0, -pz, rort, &mut x, 0.0, 1.0,
                    )
                };
                let rr = rust::rt_SHa_real(
                    stiffness, stiffness, 0.0, stiffness, stiffness, 0.0, 0.0, pz, 0.0, pz, 0.0,
                    -pz, rort, &mut y, 0.0, 1.0,
                );
                assert_eq!(rc, rr);
                exact(&[x], &[y]);
            }
        }
    }
}

fn slow_test(sv: bool) {
    let mut rng = Rng::new(33);
    for i in 0..5000 {
        let s = match i % 17 {
            0 => 0.0,
            1 => -0.0,
            2 => 1.0,
            3 => -1.0,
            4 => 0.7,
            5 => -0.7,
            6 => f32::from_bits(0.7f32.to_bits() + 1),
            7 => f32::from_bits(0.7f32.to_bits() - 1),
            8 => (0.5f32).sqrt(),
            9 => f32::from_bits((0.5f32).sqrt().to_bits() - 1),
            _ => rng.f() / 2.0,
        };
        let c = (1.0 - s * s).sqrt();
        let p = match i % 13 {
            0 => 0.0,
            1 => 1.0,
            2 => -1.0,
            3 => 1e4,
            _ => rng.f() * 0.2,
        };
        let [a, b, f, l, m, n] = match i % 31 {
            0 => [0.0; 6],
            1 => [-4.0, -3.0, -1.0, -1.0, 0.0, 0.0],
            2 => [4.0, 4.0, 2.0, 1.0, 0.0, 0.0],
            3 => [4.0, 3.0, 1.0, 0.0, 0.0, 0.0],
            _ => [
                4.0 + rng.f() * 0.2,
                3.0 + rng.f() * 0.2,
                1.0 + rng.f() * 0.2,
                1.0 + rng.f() * 0.1,
                rng.f() * 0.03,
                rng.f() * 0.03,
            ],
        };
        let rt = if i % 2 == 0 { 1.0 } else { -1.0 };
        let (mut x, mut z, mut cx, mut cz) = (0.0, 0.0, 0.0, 0.0);
        let rc = unsafe {
            if sv {
                c::solveForSlowqSV(s, c, p, a, b, f, l, m, n, &mut cx, &mut cz, rt)
            } else {
                c::solveForSlowqP(s, c, p, a, b, f, l, m, n, &mut cx, &mut cz, rt)
            }
        };
        let rr = if sv {
            rust::solveForSlowqSV(s, c, p, a, b, f, l, m, n, &mut x, &mut z, rt)
        } else {
            rust::solveForSlowqP(s, c, p, a, b, f, l, m, n, &mut x, &mut z, rt)
        };
        assert_eq!(rc, rr);
        exact(&[x, z], &[cx, cz]);
    }
}
#[test]
fn solve_for_slow_qp() {
    slow_test(false)
}
#[test]
fn solve_for_slow_qsv() {
    slow_test(true)
}

#[test]
fn rt_iso_real() {
    let mut r = Rng::new(107);
    for i in 0..10000 {
        let (v1, v2, s1, s2, d1, d2, p) = (
            3.0 + r.f() * 0.5,
            3.0 + r.f() * 0.5,
            1.5 + r.f() * 0.3,
            1.5 + r.f() * 0.3,
            2.0 + r.f() * 0.3,
            2.0 + r.f() * 0.3,
            r.f() * 0.3,
        );
        let mi = if i % 5 == 0 { 2 } else { i % 2 };
        let mt = if mi == 2 { 2 } else { (i / 2) % 2 };
        let rt = (i / 4) % 3;
        let (mut x, mut y) = (0.0, 0.0);
        let a = unsafe { c::rt_iso_real(v1, v2, s1, s2, d1, d2, p, mi, mt, rt, &mut x) };
        let b = rust::rt_iso_real(v1, v2, s1, s2, d1, d2, p, mi, mt, rt, &mut y);
        assert_eq!(a, b);
        exact(&[x], &[y]);
    }
    iso_edges(false);
}

#[test]
fn rt_iso_cmplx() {
    let mut r = Rng::new(109);
    for i in 0..400000 {
        let (v1, v2, s1, s2, d1, d2, p) = (
            3.0 + r.f() * 0.5,
            3.0 + r.f() * 0.5,
            1.5 + r.f() * 0.3,
            1.5 + r.f() * 0.3,
            2.0 + r.f() * 0.3,
            2.0 + r.f() * 0.3,
            r.f() * 0.5,
        );
        let mi = if i % 5 == 0 { 2 } else { i % 2 };
        let mt = if mi == 2 { 2 } else { (i / 2) % 2 };
        let rt = (i / 4) % 2;
        let (mut x, mut y, mut ph, mut qh) = (0.0, 0.0, 0.0, 0.0);
        let a = unsafe { c::rt_iso_cmplx(v1, v2, s1, s2, d1, d2, p, mi, mt, rt, &mut x, &mut ph) };
        let b = rust::rt_iso_cmplx(v1, v2, s1, s2, d1, d2, p, mi, mt, rt, &mut y, &mut qh);
        assert_eq!(a, b);
        exact(&[x, ph], &[y, qh]);
    }
    iso_edges(true);
}

fn iso_edges(complex: bool) {
    for [vp1, vp2, vs1, vs2, rho1, rho2] in [
        [3.0, 3.0, 1.5, 1.5, 2.0, 2.0],
        [2.0, 4.0, 1.0, 2.0, 2.0, 3.0],
        [0.0; 6],
        [3.0, 4.0, 0.0, 0.0, 2.0, 3.0],
        [3.0, 4.0, 1.0, 2.0, 0.0, 0.0],
        [-3.0, 4.0, -1.0, 2.0, -2.0, 3.0],
        [3.0, 4.0, 1.0, 2.0, 1e-12, 1e-12],
        [3.0, 4.0, 1.0, 2.0, 1e12, 1e12],
    ] {
        for p in [
            0.0,
            -0.0,
            0.25,
            f32::from_bits(0.25f32.to_bits() - 1),
            f32::from_bits(0.25f32.to_bits() + 1),
            1.0 / 3.0,
            -0.25,
            0.5,
            1.0,
            -1.0,
            1e-12,
            1e3,
        ] {
            for (mi, mt) in [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2), (0, 2), (2, 0)] {
                for rt in 0..if complex { 2 } else { 3 } {
                    let (mut x, mut y, mut ph, mut qh) = (123.0, 123.0, -456.0, -456.0);
                    let (rc, rr) = if complex {
                        (
                            unsafe {
                                c::rt_iso_cmplx(
                                    vp1, vp2, vs1, vs2, rho1, rho2, p, mi, mt, rt, &mut x, &mut ph,
                                )
                            },
                            rust::rt_iso_cmplx(
                                vp1, vp2, vs1, vs2, rho1, rho2, p, mi, mt, rt, &mut y, &mut qh,
                            ),
                        )
                    } else {
                        (
                            unsafe {
                                c::rt_iso_real(
                                    vp1, vp2, vs1, vs2, rho1, rho2, p, mi, mt, rt, &mut x,
                                )
                            },
                            rust::rt_iso_real(
                                vp1, vp2, vs1, vs2, rho1, rho2, p, mi, mt, rt, &mut y,
                            ),
                        )
                    };
                    assert_eq!(rc, rr);
                    exact(&[x, ph], &[y, qh]);
                }
            }
        }
    }
}

// rt_ani_real's polar, findqPqSV and gvelreal helpers live in elaray.c,
// rather than libtrielas.a. Compile those ORIGINAL definitions alongside the
// original coefficient routine, without altering vendor or the frozen cref shim.
// All compiler artifacts stay alongside this test binary in CARGO_TARGET_DIR.
fn ani_reference() -> usize {
    use std::{ffi::CString, path::Path, process::Command, sync::OnceLock};
    static HANDLE: OnceLock<usize> = OnceLock::new();
    *HANDLE.get_or_init(|| {
        let ws = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let out = std::env::current_exe().unwrap().parent().unwrap()
            .join(format!("trielas-reference-{}", std::process::id()));
        std::fs::create_dir_all(&out).unwrap();
        let elaray = include_str!("../../../vendor/su-original/src/Trielas/main/elaray.c");
        let mut source = String::from("#include \"par.h\"\n#include \"elastic.h\"\nchar *sdoc[] = { NULL };\n");
        for marker in ["static void gvelreal (", "void polar(", "int findqPqSV("] {
            let tail = &elaray[elaray.rfind(marker).expect("original helper missing")..];
            let open = tail.find('{').unwrap();
            let mut depth = 1;
            let mut end = open + 1;
            for (offset, ch) in tail[open+1..].char_indices() {
                if ch == '{' { depth += 1; }
                if ch == '}' { depth -= 1; }
                if depth == 0 { end = open+1+offset+1; break; }
            }
            assert_eq!(depth, 0);
            source.push_str(&tail[..end]);
            source.push('\n');
        }
        // elastic.h has no include guard; it is already included above.
        source.push_str(&include_str!("../../../vendor/su-original/src/Trielas/lib/rt_ani_real.c")
            .replace("#include \"elastic.h\"", ""));
        source.push_str("\nint oxcref_ani_case(float *v, int *m, float *out, FILE *fp) { static int initialized=0; static char *argv[] = { \"trielas-reference\", NULL }; if (!initialized) { initargs(1,argv); initialized=1; } return rt_ani_real(");
        for i in 0..26 { source.push_str(&format!("v[{i}],")); }
        source.push_str("m[0],m[1],m[2],out,fp); }\n");
        source.push_str("int oxcref_root(float s,float c,float p,float *a,int mode,float rt,float *o,FILE *fp) { return findqPqSV(s,c,p,a[0],a[1],a[2],a[3],a[4],a[5],mode,&o[0],&o[1],&o[2],&o[3],&o[4],&o[5],&o[6],rt,fp); }\n");
        let src = out.join("reference.c");
        std::fs::write(&src, source).unwrap();
        let lib = out.join("reference.dylib");
        let mut cmd = Command::new("/usr/bin/cc");
        cmd.args(["-std=c90", "-ffp-contract=off", "-fPIC", "-Wno-implicit-function-declaration"]);
        cmd.arg(if cfg!(target_os="macos") { "-dynamiclib" } else { "-shared" });
        cmd.arg("-I").arg(ws.join("vendor/su-build/include"));
        cmd.arg("-I").arg(ws.join("vendor/su-original/src/Trielas/include"));
        cmd.arg(&src).arg("-L").arg(ws.join("vendor/su-build/lib"));
        cmd.args(["-ltrielas", "-lpar", "-lcwp", "-lm", "-o"]).arg(&lib);
        let compiled = cmd.output().expect("compile original C anisotropic reference");
        assert!(compiled.status.success(), "{}", String::from_utf8_lossy(&compiled.stderr));
        let path = CString::new(lib.to_str().unwrap()).unwrap();
        let handle = unsafe { dlopen(path.as_ptr(), 2) };
        assert!(!handle.is_null(), "dlopen: {}", unsafe { std::ffi::CStr::from_ptr(dlerror()) }.to_string_lossy());
        handle as usize
    })
}

unsafe extern "C" {
    fn dlopen(path: *const std::ffi::c_char, mode: i32) -> *mut std::ffi::c_void;
    fn dlsym(handle: *mut std::ffi::c_void, name: *const std::ffi::c_char)
        -> *mut std::ffi::c_void;
    fn dlerror() -> *const std::ffi::c_char;
    fn tmpfile() -> *mut std::ffi::c_void;
    fn fclose(fp: *mut std::ffi::c_void) -> i32;
    fn rewind(fp: *mut std::ffi::c_void);
    fn ftell(fp: *mut std::ffi::c_void) -> std::ffi::c_long;
    fn fread(dst: *mut u8, size: usize, n: usize, fp: *mut std::ffi::c_void) -> usize;
}

fn ani_rust(
    v: &[f32; 26],
    m: [i32; 3],
    out: &mut f32,
    fp: Option<&oxwave_core::cio::CFile>,
) -> i32 {
    rust::rt_ani_real(
        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13],
        v[14], v[15], v[16], v[17], v[18], v[19], v[20], v[21], v[22], v[23], v[24], v[25], m[0],
        m[1], m[2], out, fp,
    )
}

#[test]
fn rt_ani_real() {
    type Coeff =
        unsafe extern "C" fn(*const f32, *const i32, *mut f32, *mut std::ffi::c_void) -> i32;
    type Root = unsafe extern "C" fn(
        f32,
        f32,
        f32,
        *const f32,
        i32,
        f32,
        *mut f32,
        *mut std::ffi::c_void,
    ) -> i32;
    let handle = ani_reference() as *mut std::ffi::c_void;
    let coeff_ptr = unsafe { dlsym(handle, c"oxcref_ani_case".as_ptr()) };
    let root_ptr = unsafe { dlsym(handle, c"oxcref_root".as_ptr()) };
    assert!(!coeff_ptr.is_null() && !root_ptr.is_null());
    let cref: Coeff = unsafe { std::mem::transmute(coeff_ptr) };
    let root: Root = unsafe { std::mem::transmute(root_ptr) };
    let fp = unsafe { tmpfile() };
    assert!(!fp.is_null());
    let mut rng = oxwave_cref::XorShift::new(0x1917);
    let mut counts = [[0; 2]; 3];
    for i in 0..12000 {
        let s = if i % 19 == 0 {
            0.0
        } else {
            rng.uniform_f32(-0.65, 0.65)
        };
        let c = ((1.0 - s * s) as f64).sqrt() as f32;
        let pl = if i % 17 == 0 {
            0.0
        } else {
            rng.uniform_f32(-0.15, 0.15)
        };
        let ai = [
            rng.uniform_f32(8.5, 9.5),
            rng.uniform_f32(8.0, 9.0),
            rng.uniform_f32(3.0, 4.0),
            rng.uniform_f32(1.8, 2.2),
            rng.uniform_f32(-0.08, 0.08),
            rng.uniform_f32(-0.08, 0.08),
        ];
        let at = if i % 23 == 0 {
            ai
        } else {
            [
                rng.uniform_f32(9.0, 10.0),
                rng.uniform_f32(8.0, 9.0),
                rng.uniform_f32(3.0, 4.0),
                rng.uniform_f32(1.8, 2.2),
                rng.uniform_f32(-0.08, 0.08),
                rng.uniform_f32(-0.08, 0.08),
            ]
        };
        let mi = [0, 1, 3, 4][i % 4];
        let mg = [0, 1, 3, 4][(i / 4) % 4];
        let rort = (i / 16 % 3) as i32;
        let mut inc = [0.0; 7];
        let mut generated = [0.0; 7];
        assert_eq!(
            unsafe {
                root(
                    s,
                    c,
                    pl,
                    ai.as_ptr(),
                    if mi == 0 || mi == 3 { 3 } else { 4 },
                    1.0,
                    inc.as_mut_ptr(),
                    fp,
                )
            },
            1,
            "incident root, case {i}"
        );
        assert_eq!(
            unsafe {
                root(
                    s,
                    c,
                    pl,
                    if rort == 0 { at.as_ptr() } else { ai.as_ptr() },
                    if mg == 0 || mg == 3 { 3 } else { 4 },
                    if rort == 0 { 1.0 } else { -1.0 },
                    generated.as_mut_ptr(),
                    fp,
                )
            },
            1,
            "generated root, case {i}"
        );
        let mut v = [0.0; 26];
        v[0] = s;
        v[1] = c;
        v[2..8].copy_from_slice(&ai);
        v[8] = rng.uniform_f32(1.8, 3.2);
        v[9..15].copy_from_slice(&at);
        v[15] = rng.uniform_f32(1.8, 3.2);
        if i % 31 == 0 && rort != 2 {
            v[8] = 0.0;
        }
        if i % 37 == 0 && rort != 2 {
            v[15] = 0.0;
        }
        if i % 41 == 0 {
            v[15] = -v[15];
        }
        v[16..18].copy_from_slice(&inc[..2]);
        v[18..21].copy_from_slice(&inc[4..]);
        v[21..23].copy_from_slice(&generated[..2]);
        v[23..26].copy_from_slice(&generated[4..]);
        let modes = [mi, mg, rort];
        let (mut x, mut y) = (123.0, 123.0);
        let diagnostics = i % 997 == 0;
        if diagnostics {
            unsafe {
                rewind(fp);
            }
        }
        let rc = unsafe { cref(v.as_ptr(), modes.as_ptr(), &mut x, fp) };
        let rr = ani_rust(&v, modes, &mut y, None);
        assert_eq!(rc, rr, "case {i}, {modes:?}, {v:?}");
        assert_eq!(rc, 1, "physical roots must reach the coefficient solve");
        exact(&[x], &[y]);
        if diagnostics {
            let len = unsafe { ftell(fp) };
            assert!(len >= 0);
            unsafe {
                rewind(fp);
            }
            let mut c_log = vec![0u8; len as usize];
            assert_eq!(
                unsafe { fread(c_log.as_mut_ptr(), 1, c_log.len(), fp) },
                c_log.len()
            );
            let path = std::env::current_exe()
                .unwrap()
                .parent()
                .unwrap()
                .join(format!(
                    "trielas-rust-diagnostics-{}.txt",
                    std::process::id()
                ));
            let rust_fp = oxwave_core::cio::CFile::from_file(std::fs::File::create(&path).unwrap());
            let mut logged = 123.0;
            assert_eq!(ani_rust(&v, modes, &mut logged, Some(&rust_fp)), rc);
            oxwave_core::cio::flush_all();
            assert_eq!(std::fs::read(&path).unwrap(), c_log, "diagnostics case {i}");
            exact(&[x], &[logged]);
            drop(rust_fp);
            std::fs::remove_file(path).unwrap();
        }
        counts[rort as usize][(mg == 1 || mg == 4) as usize] += 1;
    }
    assert!(counts.iter().flatten().all(|&n| n >= 1900));

    // Preserve the C free-surface pole check at normal incidence, including
    // its unchanged coefficient on error (the P traction divisor is zero).
    let mut edge = [0.0; 26];
    edge[1] = 1.0;
    edge[2..8].copy_from_slice(&[9.0, 9.0, 5.0, 2.0, 0.0, 0.0]);
    edge[8] = 2.0;
    edge[9..15].copy_from_slice(&[9.0, 9.0, 5.0, 2.0, 0.0, 0.0]);
    edge[15] = 2.0;
    edge[17] = 1.0 / 3.0;
    edge[20] = 1.0;
    edge[22] = -1.0 / 3.0;
    edge[25] = 1.0;
    let (mut x, mut y) = (123.0, 123.0);
    assert_eq!(
        unsafe { cref(edge.as_ptr(), [0, 0, 2].as_ptr(), &mut x, fp) },
        -1
    );
    assert_eq!(ani_rust(&edge, [0, 0, 2], &mut y, None), -1);
    exact(&[x, y], &[123.0, 123.0]);
    // No propagating P root at these tangential slownesses. Use a real FILE
    // for C because its failed-retry diagnostic dereferences ifp unconditionally.
    for tangential in [2.0, 10.0, 1e3] {
        edge[16] = tangential;
        let (mut x, mut y) = (123.0, 123.0);
        assert_eq!(
            unsafe { cref(edge.as_ptr(), [0, 0, 0].as_ptr(), &mut x, fp) },
            0
        );
        assert_eq!(ani_rust(&edge, [0, 0, 0], &mut y, None), 0);
        exact(&[x, y], &[123.0, 123.0]);
    }
    assert_eq!(unsafe { fclose(fp) }, 0);
}

#[test]
fn rt_iso_cmplx_rcabs_rounding() {
    // Audit M3 regression (audits/claude-wave3-tri-libs.md): squaring rcabs's
    // ratio in f64 instead of f32 produces 0.11859904, one ulp below the
    // original C coefficient.
    let (vp1, vp2, vs1, vs2, rho1, rho2, pl) = (
        3.5854886, 2.0740638, 1.0834757, 1.3912363, 2.14403, 2.0636044, 0.47742823,
    );
    let (mut c_coeff, mut c_phase) = (123.0, -456.0);
    let (mut coeff, mut phase) = (123.0, -456.0);
    assert_eq!(
        unsafe {
            c::rt_iso_cmplx(
                vp1, vp2, vs1, vs2, rho1, rho2, pl, 1, 0, 1, &mut c_coeff, &mut c_phase,
            )
        },
        1
    );
    assert_eq!(
        rust::rt_iso_cmplx(
            vp1, vp2, vs1, vs2, rho1, rho2, pl, 1, 0, 1, &mut coeff, &mut phase,
        ),
        1
    );
    exact(&[c_coeff], &[0.11859905]);
    exact(&[coeff, phase], &[c_coeff, c_phase]);
}
