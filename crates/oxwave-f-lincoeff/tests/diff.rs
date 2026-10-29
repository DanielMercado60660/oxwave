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

//! Differential tests against the original `par/lib/lincoeff.c` in libpar.
//! Credits from the C source: Author: CWP: Petr Jilek, December 1999.
//! These scalar routines and fixed-size ErrorFlag arrays have no length argument.

#![allow(non_snake_case, clippy::too_many_arguments)]

use oxwave_cref::XorShift;
use oxwave_f_lincoeff::{self as rust, ErrorFlag, Model};
use std::sync::Mutex;

// libpar expects these symbols from its host application. Every test that calls
// C takes this lock, so global medium parameters cannot race between tests.
static C_LOCK: Mutex<()> = Mutex::new(());
mod c {
    use super::{ErrorFlag, Model};
    #[no_mangle]
    pub static mut vp1: f32 = 0.0;
    #[no_mangle]
    pub static mut vp2: f32 = 0.0;
    #[no_mangle]
    pub static mut vs1: f32 = 0.0;
    #[no_mangle]
    pub static mut vs2: f32 = 0.0;
    #[no_mangle]
    pub static mut rho1: f32 = 0.0;
    #[no_mangle]
    pub static mut rho2: f32 = 0.0;
    #[no_mangle]
    pub static mut eps1_1: f32 = 0.0;
    #[no_mangle]
    pub static mut eps1_2: f32 = 0.0;
    #[no_mangle]
    pub static mut delta1_1: f32 = 0.0;
    #[no_mangle]
    pub static mut delta1_2: f32 = 0.0;
    #[no_mangle]
    pub static mut delta1_3: f32 = 0.0;
    #[no_mangle]
    pub static mut gamma1_1: f32 = 0.0;
    #[no_mangle]
    pub static mut gamma1_44: f32 = 0.0;
    #[no_mangle]
    pub static mut gamma1_55: f32 = 0.0;
    #[no_mangle]
    pub static mut eps2_1: f32 = 0.0;
    #[no_mangle]
    pub static mut eps2_2: f32 = 0.0;
    #[no_mangle]
    pub static mut delta2_1: f32 = 0.0;
    #[no_mangle]
    pub static mut delta2_2: f32 = 0.0;
    #[no_mangle]
    pub static mut delta2_3: f32 = 0.0;
    #[no_mangle]
    pub static mut gamma2_44: f32 = 0.0;
    #[no_mangle]
    pub static mut gamma2_55: f32 = 0.0;
    pub unsafe fn set_model(m: &Model) {
        vp1 = m.vp1;
        vp2 = m.vp2;
        vs1 = m.vs1;
        vs2 = m.vs2;
        rho1 = m.rho1;
        rho2 = m.rho2;
        eps1_1 = m.eps1_1;
        eps1_2 = m.eps1_2;
        delta1_1 = m.delta1_1;
        delta1_2 = m.delta1_2;
        delta1_3 = m.delta1_3;
        gamma1_1 = m.gamma1_1;
        gamma1_44 = m.gamma1_44;
        gamma1_55 = m.gamma1_55;
        eps2_1 = m.eps2_1;
        eps2_2 = m.eps2_2;
        delta2_1 = m.delta2_1;
        delta2_2 = m.delta2_2;
        delta2_3 = m.delta2_3;
        gamma2_44 = m.gamma2_44;
        gamma2_55 = m.gamma2_55;
    }
    extern "C" {
        pub fn P_err_2nd_order(
            rp_1st: *mut ErrorFlag,
            rp_2nd: *mut ErrorFlag,
            true_kappa: f32,
            index: i32,
        ) -> i32;
        pub fn P_term(
            Da_a: f32,
            Dz_z: f32,
            DG_G: f32,
            Delta1_1: f32,
            Delta1_2: f32,
            Delta1_3: f32,
            Delta2_1: f32,
            Delta2_2: f32,
            Delta2_3: f32,
            Eps1_1: f32,
            Eps1_2: f32,
            Eps2_1: f32,
            Eps2_2: f32,
            Gamma1_44: f32,
            Gamma2_44: f32,
            alpha: f32,
            beta: f32,
            ang: f32,
            azim: f32,
            a_k: f32,
        ) -> f32;
        pub fn Phi_rot(
            rs1: *mut f32,
            rs2: *mut f32,
            iso_plane: i32,
            pb_x: f32,
            pb_y: f32,
            pb_z: f32,
            gs1_x: f32,
            gs1_y: f32,
            gs1_z: f32,
            gs2_x: f32,
            gs2_y: f32,
            gs2_z: f32,
            CPhi1: *mut f32,
            SPhi1: *mut f32,
            CPhi2: *mut f32,
            SPhi2: *mut f32,
        ) -> i32;
        pub fn SH_term(
            Delta1_1: f32,
            Delta1_2: f32,
            Delta1_3: f32,
            Delta2_1: f32,
            Delta2_2: f32,
            Delta2_3: f32,
            Eps1_1: f32,
            Eps1_2: f32,
            Eps2_1: f32,
            Eps2_2: f32,
            Gamma1_44: f32,
            Gamma2_44: f32,
            alpha: f32,
            beta: f32,
            ang: f32,
            azim: f32,
            kappa: f32,
        ) -> f32;
        pub fn SV_term(
            DG_G: f32,
            Dr_r: f32,
            Delta1_1: f32,
            Delta1_2: f32,
            Delta1_3: f32,
            Delta2_1: f32,
            Delta2_2: f32,
            Delta2_3: f32,
            Eps1_1: f32,
            Eps1_2: f32,
            Eps2_1: f32,
            Eps2_2: f32,
            Gamma1_44: f32,
            Gamma2_44: f32,
            alpha: f32,
            beta: f32,
            ang: f32,
            azim: f32,
            kappa: f32,
        ) -> f32;
        pub fn S_err_2nd_order(
            rsv_1st: *mut ErrorFlag,
            rsv_2nd: *mut ErrorFlag,
            rsh_1st: *mut ErrorFlag,
            rsh_2nd: *mut ErrorFlag,
            true_kappa: f32,
            index: i32,
        ) -> i32;
        pub fn lincoef_Rp(
            ang: f32,
            azim: f32,
            kappa: f32,
            rpp: *mut f32,
            rp_1st: *mut ErrorFlag,
            rp_2nd: *mut ErrorFlag,
            count: i32,
        ) -> f32;
        pub fn lincoef_Rs(
            ang: f32,
            azim: f32,
            kappa: f32,
            rps1: *mut f32,
            rps2: *mut f32,
            sv: *mut f32,
            sh: *mut f32,
            cphi: *mut f32,
            sphi: *mut f32,
            i_hsp: i32,
            rsv_1st: *mut ErrorFlag,
            rsv_2nd: *mut ErrorFlag,
            rsh_1st: *mut ErrorFlag,
            rsh_2nd: *mut ErrorFlag,
            count: i32,
        ) -> f32;
        pub fn polar(ang: f32, azim: f32, cphi: *mut f32, sphi: *mut f32, i_hsp: i32) -> i32;
    }
}

fn same(c: &[f32], r: &[f32], context: &str, case: usize) {
    assert_eq!(c.len(), r.len());
    let mut scale = 0.0_f64;
    for (&a, &b) in c.iter().zip(r) {
        if a.is_finite() {
            scale = scale.max((a as f64).abs());
        }
        if b.is_finite() {
            scale = scale.max((b as f64).abs());
        }
    }
    for (j, (&a, &b)) in c.iter().zip(r).enumerate() {
        if a.is_nan() {
            assert!(b.is_nan(), "{context} case {case} [{j}] C={a} Rust={b}");
        } else if a.is_infinite() || b.is_infinite() {
            assert_eq!(a, b, "{context} case {case} [{j}]");
        } else {
            assert!(
                b.is_finite() && (a as f64 - b as f64).abs() <= 1e-5 * scale + 1e-30,
                "{context} case {case} [{j}] C={a:?} Rust={b:?} scale={scale}\nC={c:?}\nRust={r:?}"
            );
        }
    }
}

fn same_flag(c: &ErrorFlag, r: &ErrorFlag, case: usize) {
    // Compare each physical quantity separately: angle units must not inflate
    // the tolerance for small reflection coefficients.
    same(&c.iso, &r.iso, "iso", case);
    same(&c.upper, &r.upper, "upper", case);
    same(&c.lower, &r.lower, "lower", case);
    same(&c.global, &r.global, "global", case);
    same(&c.angle, &r.angle, "angle", case);
}

fn flag_bits(e: &ErrorFlag) -> Vec<u32> {
    e.iso
        .iter()
        .chain(&e.upper)
        .chain(&e.lower)
        .chain(&e.global)
        .chain(&e.angle)
        .map(|v| v.to_bits())
        .collect()
}

fn unchanged(c: &[ErrorFlag], before: &[ErrorFlag]) {
    for (a, b) in c.iter().zip(before) {
        assert_eq!(flag_bits(a), flag_bits(b));
    }
}

fn initial(rng: &mut XorShift) -> ErrorFlag {
    let mut e = ErrorFlag::default();
    for x in e
        .iso
        .iter_mut()
        .chain(&mut e.upper)
        .chain(&mut e.lower)
        .chain(&mut e.global)
        .chain(&mut e.angle)
    {
        *x = rng.uniform_f32(-0.25, 0.25);
    }
    e
}

fn model(rng: &mut XorShift, i: usize) -> Model {
    let mut m = Model {
        vp1: rng.uniform_f32(1.5, 8.0),
        vp2: 0.0,
        vs1: 0.0,
        vs2: 0.0,
        rho1: rng.uniform_f32(1.0, 3.5),
        rho2: 0.0,
        eps1_1: rng.uniform_f32(-0.25, 0.25),
        eps1_2: rng.uniform_f32(-0.25, 0.25),
        delta1_1: rng.uniform_f32(-0.25, 0.25),
        delta1_2: rng.uniform_f32(-0.25, 0.25),
        delta1_3: rng.uniform_f32(-0.25, 0.25),
        gamma1_1: rng.uniform_f32(-0.25, 0.25),
        gamma1_44: rng.uniform_f32(-0.25, 0.25),
        gamma1_55: rng.uniform_f32(-0.25, 0.25),
        eps2_1: rng.uniform_f32(-0.25, 0.25),
        eps2_2: rng.uniform_f32(-0.25, 0.25),
        delta2_1: rng.uniform_f32(-0.25, 0.25),
        delta2_2: rng.uniform_f32(-0.25, 0.25),
        delta2_3: rng.uniform_f32(-0.25, 0.25),
        gamma2_44: rng.uniform_f32(-0.25, 0.25),
        gamma2_55: rng.uniform_f32(-0.25, 0.25),
    };
    m.vp2 = m.vp1 * rng.uniform_f32(0.7, 1.4);
    m.vs1 = m.vp1 * rng.uniform_f32(0.35, 0.65);
    m.vs2 = m.vs1 * rng.uniform_f32(0.7, 1.4);
    m.rho2 = m.rho1 * rng.uniform_f32(0.7, 1.4);
    // Exact zeros, equal media, mixed signs and tiny anisotropy exercise all
    // sign-dependent empirical corrections and interpolation fallbacks.
    if i % 11 == 0 {
        m.eps1_1 = 0.0;
        m.eps1_2 = 0.0;
        m.delta1_1 = 0.0;
        m.delta1_2 = 0.0;
        m.delta1_3 = 0.0;
        m.gamma1_1 = 0.0;
        m.gamma1_44 = 0.0;
        m.gamma1_55 = 0.0;
        m.eps2_1 = 0.0;
        m.eps2_2 = 0.0;
        m.delta2_1 = 0.0;
        m.delta2_2 = 0.0;
        m.delta2_3 = 0.0;
        m.gamma2_44 = 0.0;
        m.gamma2_55 = 0.0;
    }
    if i % 13 == 0 {
        m.vp2 = m.vp1;
        m.vs2 = m.vs1;
        m.rho2 = m.rho1;
    }
    if i % 17 == 0 {
        m.delta1_1 = 0.0;
        m.delta2_2 = 0.0;
        m.gamma1_44 = 0.0;
    }
    if i % 19 == 0 {
        m.gamma1_55 = 1e-12;
        m.gamma1_44 = -1e-12;
    }
    if i % 127 == 0 {
        m.vp1 = -m.vp1;
        m.vp2 = -m.vp2;
        m.vs1 = -m.vs1;
        m.vs2 = -m.vs2;
    }
    if i % 131 == 0 {
        m.rho1 = -m.rho1;
        m.rho2 = -m.rho2;
    }
    if i % 137 == 0 {
        m.vp2 = 2.0 * m.vp1;
    } // exact isotropic helper's postcritical sentinel
    if i % 139 == 0 {
        m.vs1 = 0.0;
    }
    if i % 149 == 0 {
        m.rho2 = 0.0;
    }
    if i % 151 == 0 {
        m.vp2 = -m.vp1;
    } // zero background alpha
    if i % 157 == 0 {
        m.vs1 = m.vp1;
        m.vs2 = m.vp2;
    } // alpha == beta
    if i % 163 == 0 {
        m.vp1 = 0.0;
    } // undefined physical model, defined IEEE operations
    if i % 167 == 0 {
        m.vp2 = m.vp1 * (1.0 + 1e-6);
        m.vs2 = m.vs1;
        m.rho2 = m.rho1;
    }
    if i % 173 == 0 {
        // Also exercise m/s and kg/m^3 instead of km/s and g/cm^3.
        m.vp1 *= 1000.0;
        m.vp2 *= 1000.0;
        m.vs1 *= 1000.0;
        m.vs2 *= 1000.0;
        m.rho1 *= 1000.0;
        m.rho2 *= 1000.0;
    }
    m
}

fn angles(rng: &mut XorShift, i: usize) -> (f32, f32, f32) {
    let ang = match i % 23 {
        0 => 0.0,
        1 => -0.0,
        2 => std::f32::consts::FRAC_PI_2,
        _ => rng.uniform_f32(-1.5, 1.5),
    };
    let azim = if i % 7 == 0 {
        0.0
    } else {
        rng.uniform_f32(-12.0, 12.0)
    };
    let kappa = match i % 7 {
        0 => 0.0,
        1 => std::f32::consts::FRAC_PI_2,
        2 => -std::f32::consts::FRAC_PI_2,
        _ => rng.uniform_f32(-6.0, 6.0),
    };
    (ang, azim, kappa)
}

fn term_inputs<const N: usize>(rng: &mut XorShift, i: usize, alpha: usize) -> [f32; N] {
    let mut x = [0.0; N];
    for v in &mut x {
        *v = rng.uniform_f32(-0.5, 0.5);
    }
    if i % 7 == 0 {
        x[..alpha].fill(0.0);
    }
    if i % 11 == 0 {
        x[i % alpha] = 0.0;
    }
    x[alpha] = rng.uniform_f32(1.5, 8.0);
    x[alpha + 1] = x[alpha] * rng.uniform_f32(0.3, 0.7);
    if i % 13 == 0 {
        x[alpha] = -x[alpha];
    }
    if i % 17 == 0 {
        x[alpha + 1] = -x[alpha + 1];
    }
    if i % 19 == 0 {
        x[alpha + 1] = 0.0;
    }
    if i % 211 == 0 {
        x[alpha] = 0.0;
    }
    if i % 223 == 0 {
        x[alpha + 1] = x[alpha];
    }
    let (a, z, k) = angles(rng, i);
    x[alpha + 2] = a;
    x[alpha + 3] = z;
    x[alpha + 4] = k;
    x
}

#[test]
fn p_term_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(17);
    for i in 0..12000 {
        let x = term_inputs::<20>(&mut rng, i, 15);
        let c = unsafe {
            c::P_term(
                x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12],
                x[13], x[14], x[15], x[16], x[17], x[18], x[19],
            )
        };
        let r = rust::P_term(
            x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13],
            x[14], x[15], x[16], x[17], x[18], x[19],
        );
        same(&[c], &[r], "P_term", i);
    }
}

#[test]
fn sv_term_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(127);
    for i in 0..12000 {
        let x = term_inputs::<19>(&mut rng, i, 14);
        let c = unsafe {
            c::SV_term(
                x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12],
                x[13], x[14], x[15], x[16], x[17], x[18],
            )
        };
        let r = rust::SV_term(
            x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13],
            x[14], x[15], x[16], x[17], x[18],
        );
        same(&[c], &[r], "SV_term", i);
    }
}

#[test]
fn sh_term_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(491);
    for i in 0..12000 {
        let x = term_inputs::<17>(&mut rng, i, 12);
        let c = unsafe {
            c::SH_term(
                x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12],
                x[13], x[14], x[15], x[16],
            )
        };
        let r = rust::SH_term(
            x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13],
            x[14], x[15], x[16],
        );
        same(&[c], &[r], "SH_term", i);
    }
}

#[test]
fn polar_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(759);
    for i in 0..12000 {
        let m = model(&mut rng, i);
        let (ang, azim, _) = angles(&mut rng, i);
        let kind = [-1, 0, 1, 2, 3, 4][i % 6];
        unsafe {
            c::set_model(&m);
        }
        let (mut cc, mut cs) = (0.123, -0.456);
        let (mut rc, mut rs) = (cc, cs);
        let cr = unsafe { c::polar(ang, azim, &mut cc, &mut cs, kind) };
        let rr = rust::polar(&m, ang, azim, &mut rc, &mut rs, kind);
        assert_eq!(
            cr, rr,
            "case {i} model={m:?} kind={kind} ang={ang} azim={azim}"
        );
        same(&[cc, cs], &[rc, rs], "polar", i);
        if cr == -2 || cr == -100 {
            assert_eq!(
                (cc.to_bits(), cs.to_bits()),
                (0.123_f32.to_bits(), (-0.456_f32).to_bits())
            );
        }
    }
    // At azimuth zero B12 is exactly zero. The two signs of B11-B22 select
    // opposite axes; all-zero anisotropy terminates before writing outputs.
    for (delta, expected) in [
        (0.2, (0.0, 1.0)),
        (-0.2, (1.0, 0.0)),
        (0.0, (0.123, -0.456)),
    ] {
        let m = Model {
            vp1: 3.0,
            vs1: 1.5,
            delta1_2: delta,
            ..Model::default()
        };
        unsafe {
            c::set_model(&m);
        }
        let (mut cc, mut cs, mut rc, mut rs) = (0.123, -0.456, 0.123, -0.456);
        let cr = unsafe { c::polar(0.5, 0.0, &mut cc, &mut cs, 3) };
        let rr = rust::polar(&m, 0.5, 0.0, &mut rc, &mut rs, 3);
        assert_eq!(cr, rr);
        assert_eq!((cc, cs), expected);
        assert_eq!((rc, rs), expected);
    }
    // At normal incidence with vs1=1, D=2*gamma1_44 exactly. Check the
    // FLT_EPSILON singularity threshold and its two adjacent float values.
    let midpoint = f32::EPSILON / 2.0;
    for gamma in [
        f32::from_bits(midpoint.to_bits() - 1),
        midpoint,
        f32::from_bits(midpoint.to_bits() + 1),
    ] {
        let m = Model {
            vp1: 2.0,
            vs1: 1.0,
            gamma1_44: gamma,
            ..Model::default()
        };
        unsafe {
            c::set_model(&m);
        }
        for kind in [1, 2, 3] {
            let (mut cc, mut cs, mut rc, mut rs) = (0.123, -0.456, 0.123, -0.456);
            let cr = unsafe { c::polar(0.0, 0.0, &mut cc, &mut cs, kind) };
            let rr = rust::polar(&m, 0.0, 0.0, &mut rc, &mut rs, kind);
            assert_eq!(cr, rr);
            let expected = if gamma <= midpoint {
                if kind == 3 {
                    -2
                } else {
                    -1
                }
            } else {
                0
            };
            assert_eq!(cr, expected);
            same(
                &[cc, cs],
                &[rc, rs],
                "polar epsilon threshold",
                kind as usize,
            );
        }
    }
}

#[test]
fn p_err_2nd_order_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(771);
    for i in 0..4000 {
        let m = model(&mut rng, i);
        let (_, _, kappa) = angles(&mut rng, i);
        unsafe {
            c::set_model(&m);
        }
        let mut ce = [initial(&mut rng), initial(&mut rng)];
        let mut re = ce;
        for index in [1, 2, 0, -1, 3, i32::MAX] {
            let before = ce;
            let rust_before = re;
            let [c1, c2] = &mut ce;
            let [r1, r2] = &mut re;
            let cr = unsafe { c::P_err_2nd_order(c1, c2, kappa, index) };
            let rr = rust::P_err_2nd_order(&m, r1, r2, kappa, index);
            assert_eq!(cr, rr, "case {i} index={index}");
            for j in 0..2 {
                same_flag(&ce[j], &re[j], i);
            }
            if index != 1 && index != 2 {
                assert_eq!(cr, -1);
                unchanged(&ce, &before);
                unchanged(&re, &rust_before);
            }
        }
    }
}

#[test]
fn s_err_2nd_order_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(819);
    for i in 0..4000 {
        let m = model(&mut rng, i);
        let (_, _, kappa) = angles(&mut rng, i);
        unsafe {
            c::set_model(&m);
        }
        let mut ce = [
            initial(&mut rng),
            initial(&mut rng),
            initial(&mut rng),
            initial(&mut rng),
        ];
        let mut re = ce;
        for index in [1, 2, 0, -1, 3, i32::MAX] {
            let before = ce;
            let rust_before = re;
            let [c1, c2, c3, c4] = &mut ce;
            let [r1, r2, r3, r4] = &mut re;
            let cr = unsafe { c::S_err_2nd_order(c1, c2, c3, c4, kappa, index) };
            let rr = rust::S_err_2nd_order(&m, r1, r2, r3, r4, kappa, index);
            assert_eq!(cr, rr, "case {i} index={index}");
            for j in 0..4 {
                same_flag(&ce[j], &re[j], i);
            }
            if index != 1 && index != 2 {
                assert_eq!(cr, -1);
                unchanged(&ce, &before);
                unchanged(&re, &rust_before);
            }
        }
    }
}

#[test]
fn lincoef_rp_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(871);
    for i in 0..10000 {
        let m = model(&mut rng, i);
        let (ang, azim, kappa) = angles(&mut rng, i);
        let count = [0, 1, 2, -1, i32::MAX][i % 5];
        unsafe {
            c::set_model(&m);
        }
        let mut ce = [initial(&mut rng), initial(&mut rng)];
        let before = ce;
        let mut re = ce;
        let [c1, c2] = &mut ce;
        let [r1, r2] = &mut re;
        let (mut cv, mut rv) = (0.321, 0.321);
        let cr = unsafe { c::lincoef_Rp(ang, azim, kappa, &mut cv, c1, c2, count) };
        let rr = rust::lincoef_Rp(&m, ang, azim, kappa, &mut rv, r1, r2, count);
        assert_eq!(cr.to_bits(), rr.to_bits(), "Rp status case {i}");
        same(&[cv], &[rv], "Rp coefficient", i);
        for j in 0..2 {
            same_flag(&ce[j], &re[j], i);
        }
        if count != 1 {
            unchanged(&ce, &before);
            unchanged(&re, &before);
        }
    }
    let m = boundary_model();
    unsafe {
        c::set_model(&m);
    }
    let mut statuses = Vec::new();
    for (i, ang) in critical_angles(&m).into_iter().enumerate() {
        let (mut c1, mut c2, mut r1, mut r2) = (
            ErrorFlag::default(),
            ErrorFlag::default(),
            ErrorFlag::default(),
            ErrorFlag::default(),
        );
        let (mut cv, mut rv) = (0.0, 0.0);
        let cr = unsafe { c::lincoef_Rp(ang, 0.3, -0.2, &mut cv, &mut c1, &mut c2, 0) };
        let rr = rust::lincoef_Rp(&m, ang, 0.3, -0.2, &mut rv, &mut r1, &mut r2, 0);
        assert_eq!(cr.to_bits(), rr.to_bits());
        same(&[cv], &[rv], "Rp warning boundary", i);
        statuses.push(cr);
    }
    assert_eq!(statuses[0], 0.0);
    assert!(statuses[2] > 0.0);
}

fn boundary_model() -> Model {
    Model {
        vp1: 3.0,
        vp2: 4.0,
        vs1: 1.5,
        vs2: 2.0,
        rho1: 2.0,
        rho2: 2.5,
        ..Model::default()
    }
}

fn critical_angles(m: &Model) -> [f32; 3] {
    // The C assignment narrows asin to float before multiplying by 0.7.
    let critical = ((m.vp1 / m.vp2) as f64).asin() as f32;
    let threshold = (0.7 * critical as f64) as f32;
    [
        f32::from_bits(threshold.to_bits() - 1),
        threshold,
        f32::from_bits(threshold.to_bits() + 1),
    ]
}

#[test]
fn lincoef_rs_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(911);
    for i in 0..12000 {
        let m = model(&mut rng, i);
        let (ang, azim, kappa) = angles(&mut rng, i);
        let kind = [-1, 0, 1, 2, 3, 4][i % 6];
        let count = [0, 1, 2, -1, i32::MAX][(i / 6) % 5];
        unsafe {
            c::set_model(&m);
        }
        let mut ce = [
            initial(&mut rng),
            initial(&mut rng),
            initial(&mut rng),
            initial(&mut rng),
        ];
        let before = ce;
        let mut re = ce;
        let [c1, c2, c3, c4] = &mut ce;
        let [r1, r2, r3, r4] = &mut re;
        let mut cv = [0.123, -0.234, 0.345, -0.456, 0.567, -0.678];
        let mut rv = cv;
        let p = cv.as_mut_ptr();
        let cr = unsafe {
            c::lincoef_Rs(
                ang,
                azim,
                kappa,
                p,
                p.add(1),
                p.add(2),
                p.add(3),
                p.add(4),
                p.add(5),
                kind,
                c1,
                c2,
                c3,
                c4,
                count,
            )
        };
        let [v1, v2, v3, v4, v5, v6] = &mut rv;
        let rr = rust::lincoef_Rs(
            &m, ang, azim, kappa, v1, v2, v3, v4, v5, v6, kind, r1, r2, r3, r4, count,
        );
        assert_eq!(cr.to_bits(), rr.to_bits(), "Rs status case {i}");
        same(&cv[..4], &rv[..4], "Rs coefficients", i);
        same(&cv[4..], &rv[4..], "Rs polarization", i);
        for j in 0..4 {
            same_flag(&ce[j], &re[j], i);
        }
        if count != 1 || cr == -2.0 {
            unchanged(&ce, &before);
            unchanged(&re, &before);
        }
    }
    let m = boundary_model();
    unsafe {
        c::set_model(&m);
    }
    let mut statuses = Vec::new();
    for (i, ang) in critical_angles(&m).into_iter().enumerate() {
        let mut ce = [ErrorFlag::default(); 4];
        let mut re = ce;
        let [c1, c2, c3, c4] = &mut ce;
        let [r1, r2, r3, r4] = &mut re;
        let mut cv = [0.0; 6];
        let mut rv = cv;
        let p = cv.as_mut_ptr();
        let cr = unsafe {
            c::lincoef_Rs(
                ang,
                0.3,
                -0.2,
                p,
                p.add(1),
                p.add(2),
                p.add(3),
                p.add(4),
                p.add(5),
                0,
                c1,
                c2,
                c3,
                c4,
                0,
            )
        };
        let [v1, v2, v3, v4, v5, v6] = &mut rv;
        let rr = rust::lincoef_Rs(
            &m, ang, 0.3, -0.2, v1, v2, v3, v4, v5, v6, 0, r1, r2, r3, r4, 0,
        );
        assert_eq!(cr.to_bits(), rr.to_bits());
        same(&cv, &rv, "Rs warning boundary", i);
        statuses.push(cr);
    }
    assert_eq!(statuses[0], 0.0);
    assert!(statuses[2] > 0.0);
}

fn phi_call(c_backend: bool, x: [f32; 11], plane: i32) -> (i32, [f32; 6]) {
    let mut v = [x[0], x[1], 0.123, -0.234, 0.345, -0.456];
    let status = if c_backend {
        let p = v.as_mut_ptr();
        unsafe {
            c::Phi_rot(
                p,
                p.add(1),
                plane,
                x[2],
                x[3],
                x[4],
                x[5],
                x[6],
                x[7],
                x[8],
                x[9],
                x[10],
                p.add(2),
                p.add(3),
                p.add(4),
                p.add(5),
            )
        }
    } else {
        let [a, b, c, d, e, f] = &mut v;
        rust::Phi_rot(
            a, b, plane, x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], c, d, e, f,
        )
    };
    (status, v)
}

#[test]
fn phi_rot_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    // A subprocess captures the real diagnostic bytes for degenerate vectors.
    // It is the same test entry point, selected explicitly for each backend.
    if let Ok(backend) = std::env::var("OXWAVE_LINCOEFF_PHI_DIAGNOSTICS") {
        for x in [
            [0.5, -0.25, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.5, -0.25, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        ] {
            let (status, v) = phi_call(backend == "c", x, 1);
            assert_eq!(status, -1);
            assert_eq!(v, [0.5, -0.25, 0.123, -0.234, 0.345, -0.456]);
        }
        return;
    }
    let mut rng = XorShift::new(827);
    for i in 0..12000 {
        let mut x = [0.0; 11];
        for v in &mut x {
            *v = rng.uniform_f32(-1.0, 1.0);
        }
        let plane = [-1, 0, 1, 2][i % 4];
        if i % 31 == 0 {
            x[2] = 0.0;
            x[3] = 0.0;
        }
        if i % 37 == 0 {
            x[5] = 0.0;
            x[6] = 0.0;
            x[7] = 0.0;
        }
        if i % 41 == 0 {
            x[2] = 0.0;
        }
        if i % 43 == 0 {
            x[2] = 0.0;
            x[3] = 0.0;
            x[4] = 1.0;
        }
        let (cr, cv) = phi_call(true, x, plane);
        let (rr, rv) = phi_call(false, x, plane);
        assert_eq!(cr, rr, "case {i}, inputs={x:?}");
        assert_eq!(
            [cv[0].to_bits(), cv[1].to_bits()],
            [rv[0].to_bits(), rv[1].to_bits()],
            "Phi sign/swap case {i}"
        );
        same(&cv[2..], &rv[2..], "Phi polarization", i);
    }
    // Exact parallel normals and axis-aligned polarizations take the fallback
    // rotation axis and sign/swap paths, without relying on random coincidences.
    for plane in [-1, 0, 1, 2] {
        for sign in [-1.0, 1.0] {
            for x in [
                [0.5, -0.25, 0.0, 0.0, sign, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                [0.0, -0.0, 0.0, sign, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0],
            ] {
                let (cr, cv) = phi_call(true, x, plane);
                let (rr, rv) = phi_call(false, x, plane);
                assert_eq!(cr, rr);
                same(&cv, &rv, "Phi axes", 0);
            }
        }
        for epsilon in [
            1e-3,
            1e-6,
            1e-12,
            1e-30,
            f32::MIN_POSITIVE,
            f32::from_bits(1),
        ] {
            let x = [0.5, -0.25, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, epsilon, 0.0];
            let (cr, cv) = phi_call(true, x, plane);
            let (rr, rv) = phi_call(false, x, plane);
            assert_eq!(cr, rr);
            same(&cv, &rv, "Phi nearly parallel", 0);
        }
    }
    let run = |backend| {
        std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "phi_rot_matches_c", "--nocapture"])
            .env("OXWAVE_LINCOEFF_PHI_DIAGNOSTICS", backend)
            .output()
            .unwrap()
    };
    let c = run("c");
    let r = run("rust");
    assert!(c.status.success(), "{}", String::from_utf8_lossy(&c.stderr));
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    assert_eq!(c.stderr, r.stderr);
    assert_eq!(c.stderr, b"! Paralel S-wave polarizations !  Either wrong data or extreme anisotropy or a vicinity of S-wave singularity! \n ! S-wave polarizations in a plane paralel to the slowness vector!   Either wrong data or extreme anisotropy \n");
}
