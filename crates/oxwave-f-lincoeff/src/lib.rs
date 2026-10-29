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

//! Linearized reflection coefficients, ported from `par/lib/lincoeff.c`.
//! Credits from the C source: Author: CWP: Petr Jilek, December 1999.
//!
//! The C application's external medium parameters are explicit in [`Model`].
//! Angles are radians unless the C routine reports a critical angle in degrees.
//! Error-analysis index 2 uses the global results from an earlier index 1 call.

// Explicit casts and mutable C temporaries preserve intermediate precision and
// assignment order, including the original zero terms in empirical corrections.
#![allow(
    non_snake_case,
    unused_parens,
    unused_mut,
    unused_assignments,
    unused_variables,
    clippy::too_many_arguments
)]

/// Preprocessed medium parameters corresponding to the C application's extern floats.
/// Velocities and densities must use consistent units. Anisotropy parameters follow
/// the conventions in the original LinRORT program; no medium validation is added.
#[derive(Clone, Copy, Debug, Default)]
pub struct Model {
    pub vp1: f32,
    pub vp2: f32,
    pub vs1: f32,
    pub vs2: f32,
    pub rho1: f32,
    pub rho2: f32,
    pub eps1_1: f32,
    pub eps1_2: f32,
    pub delta1_1: f32,
    pub delta1_2: f32,
    pub delta1_3: f32,
    pub gamma1_1: f32,
    pub gamma1_44: f32,
    pub gamma1_55: f32,
    pub eps2_1: f32,
    pub eps2_2: f32,
    pub delta2_1: f32,
    pub delta2_2: f32,
    pub delta2_3: f32,
    pub gamma2_44: f32,
    pub gamma2_55: f32,
}

/// Error-analysis outputs, with the same layout as `ErrorFlag` in `par.h`.
/// Unwritten entries retain their incoming values, including on early returns.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ErrorFlag {
    pub iso: [f32; 5],
    pub upper: [f32; 2],
    pub lower: [f32; 2],
    pub global: [f32; 4],
    pub angle: [f32; 4],
}

/// Weak-contrast, weak-anisotropy P/P reflection coefficient.
/// `count == 1` also evaluates both error structures. Returns a critical-angle
/// warning in degrees, or zero.
///
/// C prototype: `float lincoef_Rp(float ang, float azim, float kappa, float * rpp, ErrorFlag * rp_1st, ErrorFlag * rp_2nd, int count);`
/// The Rust `model` argument supplies the C external medium parameters.
pub fn lincoef_Rp(
    model: &Model,
    ang: f32,
    azim: f32,
    kappa: f32,
    rpp: &mut f32,
    rp_1st: &mut ErrorFlag,
    rp_2nd: &mut ErrorFlag,
    count: i32,
) -> f32 {
    let vp1 = model.vp1;
    let vp2 = model.vp2;
    let vs1 = model.vs1;
    let vs2 = model.vs2;
    let rho1 = model.rho1;
    let rho2 = model.rho2;
    let eps1_1 = model.eps1_1;
    let eps1_2 = model.eps1_2;
    let delta1_1 = model.delta1_1;
    let delta1_2 = model.delta1_2;
    let delta1_3 = model.delta1_3;
    let gamma1_44 = model.gamma1_44;
    let gamma1_55 = model.gamma1_55;
    let eps2_1 = model.eps2_1;
    let eps2_2 = model.eps2_2;
    let delta2_1 = model.delta2_1;
    let delta2_2 = model.delta2_2;
    let delta2_3 = model.delta2_3;
    let gamma2_44 = model.gamma2_44;
    let gamma2_55 = model.gamma2_55;
    let mut err: f32 = 0.0_f32;
    let mut ANG_CR: f32;
    // to control thinks
    let mut Da_a: f32;
    let mut DG_G: f32;
    let mut alpha: f32;
    let mut beta: f32;
    let mut a_k: f32;
    // to simplify thinks
    let mut abs_lincoeff: f32;
    let mut S2: f32;
    let mut S2T2: f32;
    // to compute thinks
    // first, an approximate check for the crytical angle
    if vp1 <= vp2 {
        if (ang as f64)
            >= 0.7_f64
                * ({
                    ANG_CR = (((vp1 / vp2) as f64).asin() as f32);
                    ANG_CR
                } as f64)
        {
            err = ANG_CR;
        }
    }
    // some quantities in advance
    Da_a = ((((vp2 - vp1) as f64) / (0.5_f64 * ((vp2 + vp1) as f64))) as f32);
    DG_G = ((((vs2 as f64).powf(2.0_f64) * (rho2 as f64)
        - (vs1 as f64).powf(2.0_f64) * (rho1 as f64))
        / (0.5_f64
            * ((vs2 as f64).powf(2.0_f64) * (rho2 as f64)
                + (vs1 as f64).powf(2.0_f64) * (rho1 as f64)))) as f32);
    alpha = ((0.5_f64 * ((vp2 + vp1) as f64)) as f32);
    beta = ((0.5_f64 * ((vs2 + vs1) as f64)) as f32);
    if count == 1_i32 {
        P_err_2nd_order(model, rp_1st, rp_2nd, kappa, 1_i32);
        if (kappa as f64) != 0.0_f64 {
            P_err_2nd_order(model, rp_1st, rp_2nd, kappa, 2_i32);
        }
    }
    // now, the individual parts of the coefficient
    a_k = azim - kappa;
    abs_lincoeff = (vp2 * rho2 - vp1 * rho1) / (vp2 * rho2 + vp1 * rho1);
    S2 = ((0.5_f64
        * ((Da_a as f64) - ((2.0_f32 * beta / alpha) as f64).powf(2.0_f64) * (DG_G as f64)
            + (delta2_1 as f64) * (a_k as f64).sin().powf(2.0_f64)
            + (delta2_2 as f64) * (a_k as f64).cos().powf(2.0_f64)
            - (delta1_1 as f64) * (azim as f64).sin().powf(2.0_f64)
            - (delta1_2 as f64) * (azim as f64).cos().powf(2.0_f64)
            - 8.0_f64
                * ((beta / alpha) as f64).powf(2.0_f64)
                * ((gamma2_55 as f64) * (a_k as f64).cos().powf(2.0_f64)
                    + (gamma2_44 as f64) * (a_k as f64).sin().powf(2.0_f64)
                    - (gamma1_55 as f64) * (azim as f64).cos().powf(2.0_f64)
                    - (gamma1_44 as f64) * (azim as f64).sin().powf(2.0_f64)))) as f32);
    S2T2 = ((0.5_f64
        * ((Da_a as f64)
            + (eps2_1 as f64) * (a_k as f64).sin().powf(4.0_f64)
            + (eps2_2 as f64)
                * ((a_k as f64).cos().powf(4.0_f64)
                    + 2.0_f64
                        * (a_k as f64).cos().powf(2.0_f64)
                        * (a_k as f64).sin().powf(2.0_f64))
            + (delta2_3 as f64)
                * (a_k as f64).cos().powf(2.0_f64)
                * (a_k as f64).sin().powf(2.0_f64)
            - (eps1_1 as f64) * (azim as f64).sin().powf(4.0_f64)
            - (eps1_2 as f64)
                * ((azim as f64).cos().powf(4.0_f64)
                    + 2.0_f64
                        * (azim as f64).sin().powf(2.0_f64)
                        * (azim as f64).cos().powf(2.0_f64))
            - (delta1_3 as f64)
                * (azim as f64).cos().powf(2.0_f64)
                * (azim as f64).sin().powf(2.0_f64))) as f32);
    // finally, make the coefficient
    *rpp = (((abs_lincoeff as f64)
        + (S2 as f64) * (ang as f64).sin().powf(2.0_f64)
        + (S2T2 as f64) * (ang as f64).sin().powf(2.0_f64) * (ang as f64).tan().powf(2.0_f64))
        as f32);
    return ((((err * 180.0_f32) as f64) / 3.141592653589793_f64) as f32);
}

/// Weak-contrast, weak-anisotropy P/S reflection coefficients.
/// `count == 1` also evaluates all four error structures. Returns the C warning
/// or polarization status, preserving outputs untouched by early termination.
///
/// C prototype: `float lincoef_Rs(float ang, float azim, float kappa, float * rps1, float * rps2, float * sv, float * sh, float * cphi, float * sphi, int i_hsp, ErrorFlag * rsv_1st, ErrorFlag * rsv_2nd, ErrorFlag * rsh_1st, ErrorFlag * rsh_2nd, int count);`
/// The Rust `model` argument supplies the C external medium parameters.
pub fn lincoef_Rs(
    model: &Model,
    ang: f32,
    azim: f32,
    kappa: f32,
    rps1: &mut f32,
    rps2: &mut f32,
    sv: &mut f32,
    sh: &mut f32,
    cphi: &mut f32,
    sphi: &mut f32,
    i_hsp: i32,
    rsv_1st: &mut ErrorFlag,
    rsv_2nd: &mut ErrorFlag,
    rsh_1st: &mut ErrorFlag,
    rsh_2nd: &mut ErrorFlag,
    count: i32,
) -> f32 {
    let vp1 = model.vp1;
    let vp2 = model.vp2;
    let vs1 = model.vs1;
    let vs2 = model.vs2;
    let rho1 = model.rho1;
    let rho2 = model.rho2;
    let eps1_1 = model.eps1_1;
    let eps1_2 = model.eps1_2;
    let delta1_1 = model.delta1_1;
    let delta1_2 = model.delta1_2;
    let delta1_3 = model.delta1_3;
    let gamma1_44 = model.gamma1_44;
    let gamma1_55 = model.gamma1_55;
    let eps2_1 = model.eps2_1;
    let eps2_2 = model.eps2_2;
    let delta2_1 = model.delta2_1;
    let delta2_2 = model.delta2_2;
    let delta2_3 = model.delta2_3;
    let gamma2_44 = model.gamma2_44;
    let gamma2_55 = model.gamma2_55;
    let mut err: f32 = 0.0_f32;
    let mut ANG_CR: f32;
    let mut perr: i32;
    // float alpha, beta, Da_a, Dr_r, Db_b, CS, S_C, CS3, S3_C, S5_C, S, CS_C, CS3_C, S3;
    let mut alpha: f32;
    let mut beta: f32;
    let mut Dr_r: f32;
    let mut CS: f32;
    let mut S_C: f32;
    let mut CS3: f32;
    let mut S3_C: f32;
    let mut S5_C: f32;
    let mut S: f32;
    let mut CS_C: f32;
    let mut CS3_C: f32;
    let mut S3: f32;
    let mut cos2a: f32;
    let mut sin2a: f32;
    let mut cos2a_k: f32;
    let mut sin2a_k: f32;
    let mut DG_G: f32;
    // first, an approximate check for the crytical angle
    if vp1 <= vp2 {
        if (ang as f64)
            >= 0.7_f64
                * ({
                    ANG_CR = (((vp1 / vp2) as f64).asin() as f32);
                    ANG_CR
                } as f64)
        {
            err = ((((ANG_CR * 180.0_f32) as f64) / 3.141592653589793_f64) as f32);
        }
    }
    // compute cosine and sine of the rotation angle PHI
    perr = polar(model, ang, azim, cphi, sphi, i_hsp);
    if (perr != 0) {
        err = (perr as f32);
    }
    if err == (-2_i32 as f32) {
        return err;
    }
    // now, determine the SV and SH terms in Rps coefficients
    alpha = ((0.5_f64 * ((vp2 + vp1) as f64)) as f32);
    beta = ((0.5_f64 * ((vs2 + vs1) as f64)) as f32);
    Dr_r = ((((rho2 - rho1) as f64) / (0.5_f64 * ((rho2 + rho1) as f64))) as f32);
    cos2a_k = (((azim - kappa) as f64).cos().powf(2.0_f64) as f32);
    sin2a_k = (((azim - kappa) as f64).sin().powf(2.0_f64) as f32);
    cos2a = ((azim as f64).cos().powf(2.0_f64) as f32);
    sin2a = ((azim as f64).sin().powf(2.0_f64) as f32);
    // for DG_G insted of Dr_R + 2Db_b
    DG_G = ((((vs2 as f64).powf(2.0_f64) * (rho2 as f64)
        - (vs1 as f64).powf(2.0_f64) * (rho1 as f64))
        / (0.5_f64
            * ((vs2 as f64).powf(2.0_f64) * (rho2 as f64)
                + (vs1 as f64).powf(2.0_f64) * (rho1 as f64)))) as f32);
    if count == 1_i32 {
        S_err_2nd_order(model, rsv_1st, rsv_2nd, rsh_1st, rsh_2nd, kappa, 1_i32);
        if (kappa as f64) != 0.0_f64 {
            S_err_2nd_order(model, rsv_1st, rsv_2nd, rsh_1st, rsh_2nd, kappa, 2_i32);
        }
    }
    // for DG_G insted of Dr_R + 2Db_b:
    CS = ((((-beta / alpha * DG_G) as f64)
        - ((alpha * beta) as f64)
            / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
            * ((delta2_2 * cos2a_k + delta2_1 * sin2a_k - delta1_1 * sin2a - delta1_2 * cos2a)
                as f64)
        - ((2.0_f32 * beta / alpha
            * (gamma2_55 * cos2a_k + gamma2_44 * sin2a_k - gamma1_55 * cos2a - gamma1_44 * sin2a))
            as f64)) as f32);
    S_C = ((-0.5_f64 * (Dr_r as f64)
        + (alpha as f64).powf(2.0_f64)
            / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
            * ((delta2_2 * cos2a_k + delta2_1 * sin2a_k - delta1_2 * cos2a - delta1_1 * sin2a)
                as f64)) as f32);
    CS3 = ((((-alpha * beta) as f64) / ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64))
        * ((delta1_1 * sin2a + delta1_2 * cos2a
            - delta1_3 * cos2a * sin2a
            - eps1_2 * (cos2a * cos2a + 2.0_f32 * cos2a * sin2a)
            - eps1_1 * sin2a * sin2a
            + eps2_2 * (cos2a_k * cos2a_k + 2.0_f32 * sin2a_k * cos2a_k)
            + eps2_1 * sin2a_k * sin2a_k
            - delta2_1 * sin2a_k
            - delta2_2 * cos2a_k
            + delta2_3 * cos2a_k * sin2a_k) as f64)) as f32);
    S3_C = (((beta as f64).powf(2.0_f64) / (alpha as f64).powf(2.0_f64) * (DG_G as f64)
        + (beta as f64).powf(2.0_f64)
            / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
            * ((delta1_1 * sin2a + delta1_2 * cos2a - delta2_1 * sin2a_k - delta2_2 * cos2a_k)
                as f64)
        + (alpha as f64).powf(2.0_f64)
            / ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64))
            * ((delta1_1 * sin2a + delta1_2 * cos2a
                - delta1_3 * cos2a * sin2a
                - eps1_1 * sin2a * sin2a
                - eps1_2 * (cos2a * cos2a + 2.0_f32 * cos2a * sin2a)
                + eps2_1 * sin2a_k * sin2a_k
                + eps2_2 * (cos2a_k * cos2a_k + 2.0_f32 * sin2a_k * cos2a_k)
                - delta2_1 * sin2a_k
                - delta2_2 * cos2a_k
                + delta2_3 * cos2a_k * sin2a_k) as f64)
        + 2.0_f64
            * ((beta / alpha) as f64).powf(2.0_f64)
            * ((gamma2_55 * cos2a_k + gamma2_44 * sin2a_k - gamma1_55 * cos2a - gamma1_44 * sin2a)
                as f64)) as f32);
    S5_C = ((-(beta as f64).powf(2.0_f64)
        / ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64))
        * ((delta1_1 * sin2a + delta1_2 * cos2a
            - delta1_3 * cos2a * sin2a
            - eps1_1 * sin2a * sin2a
            - eps1_2 * (cos2a * cos2a + 2.0_f32 * cos2a * sin2a)
            + eps2_1 * sin2a_k * sin2a_k
            + eps2_2 * (cos2a_k * cos2a_k + 2.0_f32 * sin2a_k * cos2a_k)
            - delta2_1 * sin2a_k
            - delta2_2 * cos2a_k
            + delta2_3 * cos2a_k * sin2a_k) as f64)) as f32);
    S = (((alpha as f64).powf(2.0_f64)
        / (4.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
        * (((delta2_2 - delta2_1) as f64) * ((2.0_f32 * (azim - kappa)) as f64).sin()
            + ((delta1_1 - delta1_2) as f64) * ((2.0_f32 * azim) as f64).sin())) as f32);
    CS_C = ((((alpha * beta) as f64)
        / (4.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
        * (((delta2_1 - delta2_2) as f64) * ((2.0_f32 * (azim - kappa)) as f64).sin()
            + ((delta1_2 - delta1_1) as f64) * ((2.0_f32 * azim) as f64).sin())
        + ((beta / alpha) as f64)
            * (((gamma2_44 - gamma2_55) as f64) * ((2.0_f32 * (azim - kappa)) as f64).sin()
                - ((gamma1_44 - gamma1_55) as f64) * ((2.0_f32 * azim) as f64).sin()))
        as f32);
    CS3_C = ((((alpha * beta) as f64)
        / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
        * (((delta1_1 - delta1_2 - delta1_3 * (cos2a - sin2a) + 2.0_f32 * (eps1_2 - eps1_1) * sin2a)
            as f64)
            * (azim as f64).sin()
            * (azim as f64).cos()
            - ((delta2_1 - delta2_2 - delta2_3 * (cos2a_k - sin2a_k)
                + 2.0_f32 * (eps2_2 - eps2_1) * sin2a_k) as f64)
                * ((azim - kappa) as f64).sin()
                * ((azim - kappa) as f64).cos())) as f32);
    S3 = (((alpha as f64).powf(2.0_f64)
        / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
        * (((-delta1_1 + delta1_2 + delta1_3 * (cos2a - sin2a)
            - 2.0_f32 * (eps1_2 - eps1_1) * sin2a) as f64)
            * (azim as f64).sin()
            * (azim as f64).cos()
            + ((delta2_1 - delta2_2 - delta2_3 * (cos2a_k - sin2a_k)
                + 2.0_f32 * (eps2_2 - eps2_1) * sin2a_k) as f64)
                * ((azim - kappa) as f64).sin()
                * ((azim - kappa) as f64).cos())) as f32);
    *sv = (((CS as f64) * (ang as f64).cos() * (ang as f64).sin()
        + (S_C as f64) * (ang as f64).sin()
        + (CS3 as f64) * (ang as f64).cos() * (ang as f64).sin().powf(3.0_f64)
        + ((S3_C as f64) + (S_C as f64) * 0.5_f64 * ((beta / alpha) as f64).powf(2.0_f64))
            * (ang as f64).sin().powf(3.0_f64)
        + ((S5_C as f64) + (S3_C as f64) * 0.5_f64 * ((beta / alpha) as f64).powf(2.0_f64))
            * (ang as f64).sin().powf(5.0_f64)
        + (S5_C as f64)
            * 0.5_f64
            * ((beta / alpha) as f64).powf(2.0_f64)
            * (ang as f64).sin().powf(7.0_f64)) as f32);
    *sh = (((S as f64) * (ang as f64).sin()
        + (CS_C as f64) * (ang as f64).cos() * (ang as f64).sin()
        + ((CS3_C as f64) + (CS_C as f64) * 0.5_f64 * ((beta / alpha) as f64).powf(2.0_f64))
            * (ang as f64).cos()
            * (ang as f64).sin().powf(3.0_f64)
        + (S3 as f64) * (ang as f64).sin().powf(3.0_f64)
        + (CS3_C as f64)
            * 0.5_f64
            * ((beta / alpha) as f64).powf(2.0_f64)
            * (ang as f64).cos()
            * (ang as f64).sin().powf(5.0_f64)) as f32);
    // finnaly, complete the Rps1, Rps2
    *rps1 = *sv * *cphi + *sh * *sphi;
    *rps2 = -*sv * *sphi + *sh * *cphi;
    return err;
}

/// Evaluate first-order and empirical semi-second-order P/P error estimates.
/// Call with index 1, then index 2 for rotated media, retaining the same outputs.
/// Invalid indices return -1 without changing the outputs.
///
/// C prototype: `int P_err_2nd_order(ErrorFlag * rp_1st, ErrorFlag * rp_2nd, float true_kappa, int index);`
/// The Rust `model` argument supplies the C external medium parameters.
pub fn P_err_2nd_order(
    model: &Model,
    rp_1st: &mut ErrorFlag,
    rp_2nd: &mut ErrorFlag,
    true_kappa: f32,
    index: i32,
) -> i32 {
    let vp1 = model.vp1;
    let vp2 = model.vp2;
    let vs1 = model.vs1;
    let vs2 = model.vs2;
    let rho1 = model.rho1;
    let rho2 = model.rho2;
    let eps1_1 = model.eps1_1;
    let eps1_2 = model.eps1_2;
    let delta1_1 = model.delta1_1;
    let delta1_2 = model.delta1_2;
    let delta1_3 = model.delta1_3;
    let gamma1_44 = model.gamma1_44;
    let eps2_1 = model.eps2_1;
    let eps2_2 = model.eps2_2;
    let delta2_1 = model.delta2_1;
    let delta2_2 = model.delta2_2;
    let delta2_3 = model.delta2_3;
    let gamma2_44 = model.gamma2_44;
    let mut Da_a: f32;
    let mut Dz_z: f32;
    let mut DG_G: f32;
    let mut Db_b: f32;
    let mut Dr_r: f32;
    let mut alpha: f32;
    let mut beta: f32;
    let mut azim: f32;
    let mut Delta1_1: f32;
    let mut Delta1_2: f32;
    let mut Delta2_1: f32;
    let mut Delta2_2: f32;
    let mut Gamma1_44: f32;
    let mut Gamma2_44: f32;
    let mut D1_1: f32;
    let mut D1_2: f32;
    let mut D2_1: f32;
    let mut D2_2: f32;
    let mut GS_1: f32;
    let mut GS_2: f32;
    let mut ISO_2nd: f32;
    let mut d11_plus: f32 = 0.0_f32;
    let mut d11_min: f32 = 0.0_f32;
    let mut d12_plus: f32 = 0.0_f32;
    let mut d12_min: f32 = 0.0_f32;
    let mut gs1_plus: f32 = 0.0_f32;
    let mut gs1_min: f32 = 0.0_f32;
    let mut d21_plus: f32 = 0.0_f32;
    let mut d21_min: f32 = 0.0_f32;
    let mut d22_plus: f32 = 0.0_f32;
    let mut d22_min: f32 = 0.0_f32;
    let mut gs2_plus: f32 = 0.0_f32;
    let mut gs2_min: f32 = 0.0_f32;
    let mut a_plus: f32 = 0.0_f32;
    let mut a_min: f32 = 0.0_f32;
    let mut b_plus: f32 = 0.0_f32;
    let mut b_min: f32 = 0.0_f32;
    let mut r_plus: f32 = 0.0_f32;
    let mut r_min: f32 = 0.0_f32;
    let mut Fst0_0: f32;
    let mut Fst0_90: f32;
    let mut Fst90_0: f32;
    let mut Fst90_90: f32;
    let mut FstK_0: f32;
    let mut FstK_90: f32;
    let mut FstKK_0: f32;
    let mut FstKK_90: f32;
    let mut Fst0_kappa: f32;
    let mut Fst90_kappa: f32;
    let mut FstK_kappa: f32;
    let mut FstKK_kappa: f32;
    let mut Snd0_0: f32;
    let mut Snd0_90: f32;
    let mut Snd90_0: f32;
    let mut Snd90_90: f32;
    let mut SndK_0: f32;
    let mut SndK_90: f32;
    let mut SndKK_0: f32;
    let mut SndKK_90: f32;
    let mut Snd0_kappa: f32;
    let mut Snd90_kappa: f32;
    let mut SndK_kappa: f32;
    let mut SndKK_kappa: f32;
    let mut err_kappa: f32 = 0.0_f32;
    let lkappa: f32 = ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    let kkappa: f32 = ((60.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    let mut D: f32;
    let mut N: i32 = 0_i32;
    // some useful quantities in advance
    Da_a = ((((vp2 - vp1) as f64) / (0.5_f64 * ((vp2 + vp1) as f64))) as f32);
    Db_b = ((((vs2 - vs1) as f64) / (0.5_f64 * ((vs2 + vs1) as f64))) as f32);
    Dz_z = (vp2 * rho2 - vp1 * rho1) / (vp2 * rho2 + vp1 * rho1);
    DG_G = ((((vs2 as f64).powf(2.0_f64) * (rho2 as f64)
        - (vs1 as f64).powf(2.0_f64) * (rho1 as f64))
        / (0.5_f64
            * ((vs2 as f64).powf(2.0_f64) * (rho2 as f64)
                + (vs1 as f64).powf(2.0_f64) * (rho1 as f64)))) as f32);
    Dr_r = ((((rho2 - rho1) as f64) / (0.5_f64 * ((rho2 + rho1) as f64))) as f32);
    alpha = ((0.5_f64 * ((vp2 + vp1) as f64)) as f32);
    beta = ((0.5_f64 * ((vs2 + vs1) as f64)) as f32);
    if index == 1_i32 {
        err_kappa = 0.0_f32;
    }
    if index == 2_i32 {
        err_kappa = ((3.141592653589793_f64 / 2.0_f64) as f32);
    }
    if (index != 1_i32 && index != 2_i32) {
        return -1_i32;
    }
    if (delta1_1 != 0.0) {
        d11_plus = ((((delta1_1 as f64).abs() + (delta1_1 as f64))
            / (2.0_f64 * (delta1_1 as f64).abs())) as f32);
    }
    if (delta1_1 != 0.0) {
        d11_min = ((((delta1_1 as f64).abs() - (delta1_1 as f64))
            / (2.0_f64 * (delta1_1 as f64).abs())) as f32);
    }
    if (delta1_2 != 0.0) {
        d12_plus = ((((delta1_2 as f64).abs() + (delta1_2 as f64))
            / (2.0_f64 * (delta1_2 as f64).abs())) as f32);
    }
    if (delta1_2 != 0.0) {
        d12_min = ((((delta1_2 as f64).abs() - (delta1_2 as f64))
            / (2.0_f64 * (delta1_2 as f64).abs())) as f32);
    }
    if (gamma1_44 != 0.0) {
        gs1_plus = ((((gamma1_44 as f64).abs() + (gamma1_44 as f64))
            / (2.0_f64 * (gamma1_44 as f64).abs())) as f32);
    }
    if (gamma1_44 != 0.0) {
        gs1_min = ((((gamma1_44 as f64).abs() - (gamma1_44 as f64))
            / (2.0_f64 * (gamma1_44 as f64).abs())) as f32);
    }
    if (delta2_1 != 0.0) {
        d21_plus = ((((delta2_1 as f64).abs() + (delta2_1 as f64))
            / (2.0_f64 * (delta2_1 as f64).abs())) as f32);
    }
    if (delta2_1 != 0.0) {
        d21_min = ((((delta2_1 as f64).abs() - (delta2_1 as f64))
            / (2.0_f64 * (delta2_1 as f64).abs())) as f32);
    }
    if (delta2_2 != 0.0) {
        d22_plus = ((((delta2_2 as f64).abs() + (delta2_2 as f64))
            / (2.0_f64 * (delta2_2 as f64).abs())) as f32);
    }
    if (delta2_2 != 0.0) {
        d22_min = ((((delta2_2 as f64).abs() - (delta2_2 as f64))
            / (2.0_f64 * (delta2_2 as f64).abs())) as f32);
    }
    if (gamma2_44 != 0.0) {
        gs2_plus = ((((gamma2_44 as f64).abs() + (gamma2_44 as f64))
            / (2.0_f64 * (gamma2_44 as f64).abs())) as f32);
    }
    if (gamma2_44 != 0.0) {
        gs2_min = ((((gamma2_44 as f64).abs() - (gamma2_44 as f64))
            / (2.0_f64 * (gamma2_44 as f64).abs())) as f32);
    }
    if (Da_a != 0.0) {
        a_plus = ((((Da_a as f64).abs() + (Da_a as f64)) / (2.0_f64 * (Da_a as f64).abs())) as f32);
    }
    if (Da_a != 0.0) {
        a_min = ((((Da_a as f64).abs() - (Da_a as f64)) / (2.0_f64 * (Da_a as f64).abs())) as f32);
    }
    if (Db_b != 0.0) {
        b_plus = ((((Db_b as f64).abs() + (Db_b as f64)) / (2.0_f64 * (Db_b as f64).abs())) as f32);
    }
    if (Db_b != 0.0) {
        b_min = ((((Db_b as f64).abs() - (Db_b as f64)) / (2.0_f64 * (Db_b as f64).abs())) as f32);
    }
    if (Dr_r != 0.0) {
        r_plus = ((((Dr_r as f64).abs() + (Dr_r as f64)) / (2.0_f64 * (Dr_r as f64).abs())) as f32);
    }
    if (Dr_r != 0.0) {
        r_min = ((((Dr_r as f64).abs() - (Dr_r as f64)) / (2.0_f64 * (Dr_r as f64).abs())) as f32);
    }
    // the angles used for the anisotropic accuracy evaluation
    (*rp_2nd).angle[0_i32 as usize] = ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    (*rp_1st).angle[0_i32 as usize] = (*rp_2nd).angle[0_i32 as usize];
    // incidence angle
    (*rp_2nd).angle[1_i32 as usize] = ((0.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    (*rp_1st).angle[1_i32 as usize] = (*rp_2nd).angle[1_i32 as usize];
    // first azimuth
    (*rp_2nd).angle[2_i32 as usize] = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    (*rp_1st).angle[2_i32 as usize] = (*rp_2nd).angle[2_i32 as usize];
    // second azimuth
    (*rp_2nd).angle[3_i32 as usize] = true_kappa;
    (*rp_1st).angle[3_i32 as usize] = (*rp_2nd).angle[3_i32 as usize];
    // kappa
    // ISO part first
    N = 0_i32;
    while N < 5_i32 {
        (*rp_1st).iso[N as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            alpha,
            beta,
            ((((15_i32 + 5_i32 * N) as f64) * 3.141592653589793_f64 / 180.0_f64) as f32),
            0.0_f32,
            0.0_f32,
        );
        (*rp_2nd).iso[N as usize] = Iso_exact(
            0_i32,
            vp1,
            vs1,
            rho1,
            vp2,
            vs2,
            rho2,
            ((((15_i32 + 5_i32 * N) as f64) * 3.141592653589793_f64 / 180.0_f64) as f32),
        );
        N += 1;
    }
    // get the ISO semi-2nd order for the inc. angle err_ang
    ISO_2nd = Iso_exact(
        0_i32,
        vp1,
        vs1,
        rho1,
        vp2,
        vs2,
        rho2,
        ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
    ) - P_term(
        Da_a,
        Dz_z,
        DG_G,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        alpha,
        beta,
        ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
        0.0_f32,
        0.0_f32,
    );
    // UPPER halfspace contribution
    Delta1_1 = (((-0.45_f64 * (d11_plus as f64) - 0.5_f64 * (d11_min as f64))
        * (delta1_1 as f64).powf(2.0_f64)
        + (-0.22_f64 * (d11_plus as f64) * (gs1_plus as f64)
            + ((0.0_f32 * d11_plus * gs1_min) as f64)
            - 0.93_f64 * (d11_min as f64) * (gs1_plus as f64)
            - 0.33_f64 * (d11_min as f64) * (gs1_min as f64))
            * (delta1_1 as f64)
            * (gamma1_44 as f64)
        + (0.74_f64 * (d11_plus as f64) * (a_plus as f64)
            + 0.59_f64 * (d11_plus as f64) * (a_min as f64)
            + 2.56_f64 * (d11_min as f64) * (a_plus as f64)
            + 0.31_f64 * (d11_min as f64) * (a_min as f64))
            * (delta1_1 as f64)
            * (Da_a as f64)
        + (0.19_f64 * (d11_plus as f64) * (b_plus as f64)
            + ((0.0_f32 * d11_plus * b_min) as f64)
            + 1.27_f64 * (d11_min as f64) * (b_plus as f64)
            + 0.31_f64 * (d11_min as f64) * (b_min as f64))
            * (delta1_1 as f64)
            * (Db_b as f64)
        + (((0.0_f32 * d11_plus * r_plus) as f64)
            + 0.18_f64 * (d11_plus as f64) * (r_min as f64)
            + 1.27_f64 * (d11_min as f64) * (r_plus as f64)
            + 0.95_f64 * (d11_min as f64) * (r_min as f64))
            * (delta1_1 as f64)
            * (Dr_r as f64)) as f32);
    D1_1 = Delta1_1;
    Delta1_2 = (((-0.45_f64 * (d12_plus as f64) - 0.5_f64 * (d12_min as f64))
        * (delta1_2 as f64).powf(2.0_f64)
        + (0.74_f64 * (d12_plus as f64) * (a_plus as f64)
            + 0.59_f64 * (d12_plus as f64) * (a_min as f64)
            + 2.56_f64 * (d12_min as f64) * (a_plus as f64)
            + 0.31_f64 * (d12_min as f64) * (a_min as f64))
            * (delta1_2 as f64)
            * (Da_a as f64)
        + (0.19_f64 * (d12_plus as f64) * (b_plus as f64)
            + ((0.0_f32 * d12_plus * b_min) as f64)
            + 1.27_f64 * (d12_min as f64) * (b_plus as f64)
            + 0.31_f64 * (d12_min as f64) * (b_min as f64))
            * (delta1_2 as f64)
            * (Db_b as f64)
        + (((0.0_f32 * d12_plus * r_plus) as f64)
            + 0.18_f64 * (d12_plus as f64) * (r_min as f64)
            + 1.27_f64 * (d12_min as f64) * (r_plus as f64)
            + 0.95_f64 * (d12_min as f64) * (r_min as f64))
            * (delta1_2 as f64)
            * (Dr_r as f64)) as f32);
    D1_2 = Delta1_2;
    Delta2_1 = 0.0_f32;
    Delta2_2 = 0.0_f32;
    Gamma1_44 = (((1.0_f64 / 3.0_f64 * (gs1_plus as f64) + 0.55_f64 * (gs1_min as f64))
        * (gamma1_44 as f64).powf(2.0_f64)
        + (0.92_f64 * (gs1_plus as f64) * (a_plus as f64)
            + 0.9_f64 * (gs1_plus as f64) * (a_min as f64)
            + 1.27_f64 * (gs1_min as f64) * (a_plus as f64)
            + 0.98_f64 * (gs1_min as f64) * (a_min as f64))
            * (gamma1_44 as f64)
            * (Da_a as f64)
        + (-1.74_f64 * (gs1_plus as f64) * (b_plus as f64)
            - 2.1_f64 * (gs1_plus as f64) * (b_min as f64)
            - 1.84_f64 * (gs1_min as f64) * (b_plus as f64)
            - 1.91_f64 * (gs1_min as f64) * (b_min as f64))
            * (gamma1_44 as f64)
            * (Db_b as f64)
        + (-1.02_f64 * (gs1_plus as f64) * (r_plus as f64)
            - 1.04_f64 * (gs1_plus as f64) * (r_min as f64)
            - 0.9_f64 * (gs1_min as f64) * (r_plus as f64)
            - 0.9_f64 * (gs1_min as f64) * (r_min as f64))
            * (gamma1_44 as f64)
            * (Dr_r as f64)) as f32);
    GS_1 = Gamma1_44;
    Gamma2_44 = 0.0_f32;
    if index == 1_i32 {
        azim = 0.0_f32;
        (*rp_1st).upper[0_i32 as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            eps1_1,
            eps1_2,
            0.0_f32,
            0.0_f32,
            gamma1_44,
            0.0_f32,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rp_2nd).upper[0_i32 as usize] = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        ) + ISO_2nd;
        azim = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rp_1st).upper[1_i32 as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            eps1_1,
            eps1_2,
            0.0_f32,
            0.0_f32,
            gamma1_44,
            0.0_f32,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rp_2nd).upper[1_i32 as usize] = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        ) + ISO_2nd;
    }
    // LOWER halfspace contribution
    Delta1_1 = 0.0_f32;
    Delta1_2 = 0.0_f32;
    Delta2_1 = (((-0.44_f64 * (d21_plus as f64) - 0.89_f64 * (d21_min as f64))
        * (delta2_1 as f64).powf(2.0_f64)
        + (-(1.0_f64 / 3.0_f64) * (d21_plus as f64) * (gs2_plus as f64)
            + ((0.0_f32 * d21_plus * gs2_min) as f64)
            - 0.65_f64 * (d21_min as f64) * (gs2_plus as f64)
            + ((0.0_f32 * d21_min * gs2_min) as f64))
            * (delta2_1 as f64)
            * (gamma2_44 as f64)
        + (1.83_f64 * (d21_plus as f64) * (a_plus as f64)
            + 1.43_f64 * (d21_plus as f64) * (a_min as f64)
            + 2.56_f64 * (d21_min as f64) * (a_plus as f64)
            + 1.35_f64 * (d21_min as f64) * (a_min as f64))
            * (delta2_1 as f64)
            * (Da_a as f64)
        + (-0.37_f64 * (d21_plus as f64) * (b_plus as f64)
            - 0.59_f64 * (d21_plus as f64) * (b_min as f64)
            - 0.55_f64 * (d21_min as f64) * (b_plus as f64)
            - 0.07_f64 * (d21_min as f64) * (b_min as f64))
            * (delta2_1 as f64)
            * (Db_b as f64)
        + (((0.0_f32 * d21_plus * r_plus + 0.0_f32 * d21_plus * r_min) as f64)
            - 0.13_f64 * (d21_min as f64) * (r_plus as f64)
            + 0.1_f64 * (d21_min as f64) * (r_min as f64))
            * (delta2_1 as f64)
            * (Dr_r as f64)) as f32);
    D2_1 = Delta2_1;
    Delta2_2 = (((-0.44_f64 * (d22_plus as f64) - 0.89_f64 * (d22_min as f64))
        * (delta2_2 as f64).powf(2.0_f64)
        + (1.83_f64 * (d22_plus as f64) * (a_plus as f64)
            + 1.43_f64 * (d22_plus as f64) * (a_min as f64)
            + 2.56_f64 * (d22_min as f64) * (a_plus as f64)
            + 1.35_f64 * (d22_min as f64) * (a_min as f64))
            * (delta2_2 as f64)
            * (Da_a as f64)
        + (-0.37_f64 * (d22_plus as f64) * (b_plus as f64)
            - 0.59_f64 * (d22_plus as f64) * (b_min as f64)
            - 0.55_f64 * (d22_min as f64) * (b_plus as f64)
            - 0.07_f64 * (d22_min as f64) * (b_min as f64))
            * (delta2_2 as f64)
            * (Db_b as f64)
        + (((0.0_f32 * d22_plus * r_plus + 0.0_f32 * d22_plus * r_min) as f64)
            - 0.13_f64 * (d22_min as f64) * (r_plus as f64)
            + 0.1_f64 * (d22_min as f64) * (r_min as f64))
            * (delta2_2 as f64)
            * (Dr_r as f64)) as f32);
    D2_2 = Delta2_2;
    Gamma1_44 = 0.0_f32;
    Gamma2_44 = (((-0.37_f64 * (gs2_plus as f64) - 0.5_f64 * (gs2_min as f64))
        * (gamma2_44 as f64).powf(2.0_f64)
        + (((0.0_f32 * gs2_plus * a_plus) as f64)
            + 1.05_f64 * (gs2_plus as f64) * (a_min as f64)
            + 0.83_f64 * (gs2_min as f64) * (a_plus as f64)
            + 0.83_f64 * (gs2_min as f64) * (a_min as f64))
            * (gamma2_44 as f64)
            * (Da_a as f64)
        + (-0.46_f64 * (gs2_plus as f64) * (b_plus as f64)
            + 0.23_f64 * (gs2_plus as f64) * (b_min as f64)
            + 0.24_f64 * (gs2_min as f64) * (b_plus as f64)
            + 0.38_f64 * (gs2_min as f64) * (b_min as f64))
            * (gamma2_44 as f64)
            * (Db_b as f64)
        + (-0.13_f64 * (gs2_plus as f64) * (r_plus as f64)
            + 0.07_f64 * (gs2_plus as f64) * (r_min as f64)
            + ((0.0_f32 * gs2_min * r_plus) as f64)
            + 0.24_f64 * (gs2_min as f64) * (r_min as f64))
            * (gamma2_44 as f64)
            * (Dr_r as f64)) as f32);
    GS_2 = Gamma2_44;
    if index == 1_i32 {
        azim = 0.0_f32;
        (*rp_1st).lower[0_i32 as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            delta2_1,
            delta2_2,
            delta2_3,
            0.0_f32,
            0.0_f32,
            eps2_1,
            eps2_2,
            0.0_f32,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            0.0_f32,
            azim,
        );
        (*rp_2nd).lower[0_i32 as usize] = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            0.0_f32,
            azim,
        ) + ISO_2nd;
        azim = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rp_1st).lower[1_i32 as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            delta2_1,
            delta2_2,
            delta2_3,
            0.0_f32,
            0.0_f32,
            eps2_1,
            eps2_2,
            0.0_f32,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            0.0_f32,
            azim,
        );
        (*rp_2nd).lower[1_i32 as usize] = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            0.0_f32,
            azim,
        ) + ISO_2nd;
    }
    // GLOBAL contribution
    // The following terms are derived on the assumption that Rpp coefficients for
    // two identicaly anisotropic halfspaces are zero (with 0 velocity contrats)
    Delta1_1 = (((D1_1 as f64)
        + (0.45_f64 * (d11_plus as f64) * (d21_plus as f64)
            - 0.45_f64 * (d11_plus as f64) * (d21_min as f64)
            - 0.5_f64 * (d11_min as f64) * (d21_plus as f64)
            + 0.5_f64 * (d11_min as f64) * (d21_min as f64))
            * (delta1_1 as f64)
            * (delta2_1 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (0.45_f64 * (d11_plus as f64) * (d22_plus as f64)
            - 0.45_f64 * (d11_plus as f64) * (d22_min as f64)
            - 0.5_f64 * (d11_min as f64) * (d22_plus as f64)
            + 0.5_f64 * (d11_min as f64) * (d22_min as f64))
            * (delta1_1 as f64)
            * (delta2_2 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (0.22_f64 * (d11_plus as f64) * (gs2_plus as f64)
            + 0.0_f64 * (d11_plus as f64) * (gs2_min as f64)
            + 0.93_f64 * (d11_min as f64) * (gs2_plus as f64)
            + 0.33_f64 * (d11_min as f64) * (gs2_min as f64))
            * (delta1_1 as f64)
            * (gamma2_44 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)) as f32);
    Delta1_2 = (((D1_2 as f64)
        + (0.45_f64 * (d12_plus as f64) * (d21_plus as f64)
            - 0.45_f64 * (d12_plus as f64) * (d21_min as f64)
            - 0.5_f64 * (d12_min as f64) * (d21_plus as f64)
            + 0.5_f64 * (d12_min as f64) * (d21_min as f64))
            * (delta1_2 as f64)
            * (delta2_1 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (0.45_f64 * (d12_plus as f64) * (d22_plus as f64)
            - 0.45_f64 * (d12_plus as f64) * (d22_min as f64)
            - 0.5_f64 * (d12_min as f64) * (d22_plus as f64)
            + 0.5_f64 * (d12_min as f64) * (d22_min as f64))
            * (delta1_2 as f64)
            * (delta2_2 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (0.22_f64 * (d12_plus as f64) * (gs2_plus as f64)
            + 0.0_f64 * (d12_plus as f64) * (gs2_min as f64)
            + 0.93_f64 * (d12_min as f64) * (gs2_plus as f64)
            + 0.33_f64 * (d12_min as f64) * (gs2_min as f64))
            * (delta1_2 as f64)
            * (gamma2_44 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)) as f32);
    Delta2_1 = (((D2_1 as f64)
        + (0.44_f64 * (d21_plus as f64) * (d11_plus as f64)
            - 0.44_f64 * (d21_plus as f64) * (d11_min as f64)
            - 0.89_f64 * (d21_min as f64) * (d11_plus as f64)
            + 0.89_f64 * (d21_min as f64) * (d11_min as f64))
            * (delta2_1 as f64)
            * (delta1_1 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (0.44_f64 * (d21_plus as f64) * (d12_plus as f64)
            - 0.44_f64 * (d21_plus as f64) * (d12_min as f64)
            - 0.89_f64 * (d21_min as f64) * (d12_plus as f64)
            + 0.89_f64 * (d21_min as f64) * (d12_min as f64))
            * (delta2_1 as f64)
            * (delta1_2 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (1.0_f64 / 3.0_f64 * (d21_plus as f64) * (gs1_plus as f64)
            + 0.0_f64 * (d21_plus as f64) * (gs1_min as f64)
            + 0.65_f64 * (d21_min as f64) * (gs1_plus as f64)
            + 0.0_f64 * (d21_min as f64) * (gs1_min as f64))
            * (delta2_1 as f64)
            * (gamma1_44 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)) as f32);
    Delta2_2 = (((D2_2 as f64)
        + (0.44_f64 * (d22_plus as f64) * (d11_plus as f64)
            - 0.44_f64 * (d22_plus as f64) * (d11_min as f64)
            - 0.89_f64 * (d22_min as f64) * (d11_plus as f64)
            + 0.89_f64 * (d22_min as f64) * (d11_min as f64))
            * (delta2_2 as f64)
            * (delta1_1 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (0.44_f64 * (d22_plus as f64) * (d12_plus as f64)
            - 0.44_f64 * (d22_plus as f64) * (d12_min as f64)
            - 0.89_f64 * (d22_min as f64) * (d12_plus as f64)
            + 0.89_f64 * (d22_min as f64) * (d12_min as f64))
            * (delta2_2 as f64)
            * (delta1_2 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (1.0_f64 / 3.0_f64 * (d22_plus as f64) * (gs1_plus as f64)
            + 0.0_f64 * (d22_plus as f64) * (gs1_min as f64)
            + 0.65_f64 * (d22_min as f64) * (gs1_plus as f64)
            + 0.0_f64 * (d22_min as f64) * (gs1_min as f64))
            * (delta2_2 as f64)
            * (gamma1_44 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)) as f32);
    Gamma1_44 = (((GS_1 as f64)
        + (0.0_f64 * (gs1_plus as f64) * (d21_plus as f64)
            + 0.0_f64 * (gs1_plus as f64) * (d21_min as f64)
            + 0.0_f64 * (gs1_min as f64) * (d21_plus as f64)
            + 0.0_f64 * (gs1_min as f64) * (d21_min as f64))
            * (gamma1_44 as f64)
            * (delta2_1 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (0.0_f64 * (gs1_plus as f64) * (d22_plus as f64)
            + 0.0_f64 * (gs1_plus as f64) * (d22_min as f64)
            + 0.0_f64 * (gs1_min as f64) * (d22_plus as f64)
            + 0.0_f64 * (gs1_min as f64) * (d22_min as f64))
            * (gamma1_44 as f64)
            * (delta2_2 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (-(1.0_f64 / 3.0_f64) * (gs1_plus as f64) * (gs2_plus as f64)
            - 0.75_f64 * (gs1_plus as f64) * (gs2_min as f64)
            - 0.95_f64 * (gs1_min as f64) * (gs2_plus as f64)
            - 0.55_f64 * (gs1_min as f64) * (gs2_min as f64))
            * (gamma1_44 as f64)
            * (gamma2_44 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)) as f32);
    Gamma2_44 = (((GS_2 as f64)
        + (0.37_f64 * (gs2_plus as f64) * (gs1_plus as f64)
            + 0.0_f64 * (gs2_plus as f64) * (gs1_min as f64)
            + 0.0_f64 * (gs2_min as f64) * (gs1_plus as f64)
            + 0.5_f64 * (gs2_min as f64) * (gs1_min as f64))
            * (gamma2_44 as f64)
            * (gamma1_44 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)) as f32);
    if index == 1_i32 {
        azim = 0.0_f32;
        (*rp_1st).global[0_i32 as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        );
        (*rp_2nd).global[0_i32 as usize] = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        ) + ISO_2nd;
        azim = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rp_1st).global[1_i32 as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        );
        (*rp_2nd).global[1_i32 as usize] = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        ) + ISO_2nd;
        azim = lkappa;
        (*rp_1st).global[2_i32 as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        );
        (*rp_2nd).global[2_i32 as usize] = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        ) + ISO_2nd;
        azim = kkappa;
        (*rp_1st).global[3_i32 as usize] = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        );
        (*rp_2nd).global[3_i32 as usize] = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        ) + ISO_2nd;
    }
    Fst0_0 = (*rp_1st).global[0_i32 as usize];
    Snd0_0 = (*rp_2nd).global[0_i32 as usize];
    Fst90_0 = (*rp_1st).global[1_i32 as usize];
    Snd90_0 = (*rp_2nd).global[1_i32 as usize];
    FstK_0 = (*rp_1st).global[2_i32 as usize];
    SndK_0 = (*rp_2nd).global[2_i32 as usize];
    FstKK_0 = (*rp_1st).global[3_i32 as usize];
    SndKK_0 = (*rp_2nd).global[3_i32 as usize];
    if index == 2_i32 {
        azim = 0.0_f32;
        Fst0_90 = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        );
        Snd0_90 = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        ) + ISO_2nd;
        Fst0_kappa = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - true_kappa,
        );
        azim = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        Fst90_90 = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        );
        Snd90_90 = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        ) + ISO_2nd;
        Fst90_kappa = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - true_kappa,
        );
        azim = lkappa;
        FstK_90 = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        );
        SndK_90 = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        ) + ISO_2nd;
        FstK_kappa = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - true_kappa,
        );
        azim = kkappa;
        FstKK_90 = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        );
        SndKK_90 = P_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - err_kappa,
        ) + ISO_2nd;
        FstKK_kappa = P_term(
            Da_a,
            Dz_z,
            DG_G,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            azim - true_kappa,
        );
        // interpolate the error for the halfspaces rotated by the angle true_kappa
        if ({
            D = Fst0_90 - Fst0_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            Snd0_kappa = ((0.5_f64 * ((Snd0_0 + Snd0_90) as f64)) as f32);
        } else {
            Snd0_kappa = Snd0_0 + (Snd0_90 - Snd0_0) / D * (Fst0_kappa - Fst0_0);
        }
        if ({
            D = Fst90_90 - Fst90_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            Snd90_kappa = ((0.5_f64 * ((Snd90_0 + Snd90_90) as f64)) as f32);
        } else {
            Snd90_kappa = Snd90_0 + (Snd90_90 - Snd90_0) / D * (Fst90_kappa - Fst90_0);
        }
        if ({
            D = FstK_90 - FstK_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            SndK_kappa = ((0.5_f64 * ((SndK_90 + SndK_0) as f64)) as f32);
        } else {
            SndK_kappa = SndK_0 + (SndK_90 - SndK_0) / D * (FstK_kappa - FstK_0);
        }
        if ({
            D = FstKK_90 - FstKK_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            SndKK_kappa = ((0.5_f64 * ((SndKK_90 + SndKK_0) as f64)) as f32);
        } else {
            SndKK_kappa = SndKK_0 + (SndKK_90 - SndKK_0) / D * (FstKK_kappa - FstKK_0);
        }
        (*rp_1st).global[0_i32 as usize] = Fst0_kappa;
        (*rp_2nd).global[0_i32 as usize] = Snd0_kappa;
        (*rp_1st).global[1_i32 as usize] = Fst90_kappa;
        (*rp_2nd).global[1_i32 as usize] = Snd90_kappa;
        (*rp_1st).global[2_i32 as usize] = FstK_kappa;
        (*rp_2nd).global[2_i32 as usize] = SndK_kappa;
        (*rp_1st).global[3_i32 as usize] = FstKK_kappa;
        (*rp_2nd).global[3_i32 as usize] = SndKK_kappa;
    }
    return 0_i32;
}

/// First-order P-wave angular term.
///
/// C prototype: `float P_term(float Da_a, float Dz_z, float DG_G, float Delta1_1, float Delta1_2, float Delta1_3, float Delta2_1, float Delta2_2, float Delta2_3, float Eps1_1, float Eps1_2, float Eps2_1, float Eps2_2, float Gamma1_44, float Gamma2_44, float alpha, float beta, float ang, float azim, float a_k);`
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
) -> f32 {
    let mut absolute: f32;
    let mut sn_2: f32;
    let mut sn_2_tg_2: f32;
    absolute = Dz_z;
    sn_2 = ((0.5_f64
        * ((Da_a as f64) - ((2.0_f32 * beta / alpha) as f64).powf(2.0_f64) * (DG_G as f64)
            + (Delta2_1 as f64) * (a_k as f64).sin().powf(2.0_f64)
            + (Delta2_2 as f64) * (a_k as f64).cos().powf(2.0_f64)
            - (Delta1_1 as f64) * (azim as f64).sin().powf(2.0_f64)
            - (Delta1_2 as f64) * (azim as f64).cos().powf(2.0_f64)
            - 8.0_f64
                * ((beta / alpha) as f64).powf(2.0_f64)
                * ((Gamma2_44 as f64) * (a_k as f64).sin().powf(2.0_f64)
                    - (Gamma1_44 as f64) * (azim as f64).sin().powf(2.0_f64)))) as f32);
    sn_2_tg_2 = ((0.5_f64
        * ((Da_a as f64)
            + (Eps2_1 as f64) * (a_k as f64).sin().powf(4.0_f64)
            + (Eps2_2 as f64)
                * ((a_k as f64).cos().powf(4.0_f64)
                    + 2.0_f64
                        * (a_k as f64).cos().powf(2.0_f64)
                        * (a_k as f64).sin().powf(2.0_f64))
            + (Delta2_3 as f64)
                * (a_k as f64).cos().powf(2.0_f64)
                * (a_k as f64).sin().powf(2.0_f64)
            - (Eps1_1 as f64) * (azim as f64).sin().powf(4.0_f64)
            - (Eps1_2 as f64)
                * ((azim as f64).cos().powf(4.0_f64)
                    + 2.0_f64
                        * (azim as f64).sin().powf(2.0_f64)
                        * (azim as f64).cos().powf(2.0_f64))
            - (Delta1_3 as f64)
                * (azim as f64).cos().powf(2.0_f64)
                * (azim as f64).sin().powf(2.0_f64))) as f32);
    return (((absolute as f64)
        + (sn_2 as f64) * (ang as f64).sin().powf(2.0_f64)
        + (sn_2_tg_2 as f64) * (ang as f64).sin().powf(2.0_f64) * (ang as f64).tan().powf(2.0_f64))
        as f32);
}

/// Private exact isotropic reflection coefficient used by error analysis.
///
/// C prototype: `float Iso_exact(int type, float vp1, float vs1, float rho1, float vp2, float vs2, float rho2, float ang);`
fn Iso_exact(
    wave_type: i32,
    vp1: f32,
    vs1: f32,
    rho1: f32,
    vp2: f32,
    vs2: f32,
    rho2: f32,
    ang: f32,
) -> f32 {
    let mut p: f32;
    let mut P1: f32;
    let mut P2: f32;
    let mut P3: f32;
    let mut P4: f32;
    let mut q: f32;
    let mut X: f32;
    let mut Y: f32;
    let mut Z: f32;
    let mut D: f32;
    let mut COEF: f32;
    // some auxiliary quantities
    p = (((ang as f64).sin() / (vp1 as f64)) as f32);
    if ({
        P1 = ((1.0_f64 - ((vp1 * p) as f64).powf(2.0_f64)) as f32);
        P1
    } as f64)
        < 0.0_f64
    {
        return (-999.9_f64 as f32);
    }
    if ({
        P2 = ((1.0_f64 - ((vs1 * p) as f64).powf(2.0_f64)) as f32);
        P2
    } as f64)
        < 0.0_f64
    {
        return (-999.9_f64 as f32);
    }
    if ({
        P3 = ((1.0_f64 - ((vp2 * p) as f64).powf(2.0_f64)) as f32);
        P3
    } as f64)
        < 0.0_f64
    {
        return (-999.9_f64 as f32);
    }
    if ({
        P4 = ((1.0_f64 - ((vs2 * p) as f64).powf(2.0_f64)) as f32);
        P4
    } as f64)
        < 0.0_f64
    {
        return (-999.9_f64 as f32);
    }
    P1 = ((P1 as f64).sqrt() as f32);
    P2 = ((P2 as f64).sqrt() as f32);
    P3 = ((P3 as f64).sqrt() as f32);
    P4 = ((P4 as f64).sqrt() as f32);
    q = ((2.0_f64
        * ((rho2 as f64) * (vs2 as f64).powf(2.0_f64) - (rho1 as f64) * (vs1 as f64).powf(2.0_f64)))
        as f32);
    X = (((rho2 as f64) - (q as f64) * (p as f64).powf(2.0_f64)) as f32);
    Y = (((rho1 as f64) + (q as f64) * (p as f64).powf(2.0_f64)) as f32);
    Z = ((((rho2 - rho1) as f64) - (q as f64) * (p as f64).powf(2.0_f64)) as f32);
    D = ((((q * p) as f64).powf(2.0_f64) * (P1 as f64) * (P2 as f64) * (P3 as f64) * (P4 as f64)
        + ((rho1 * rho2 * (vp1 * vs2 * P2 * P3 + vs1 * vp2 * P1 * P4)) as f64)
        + ((vp1 * vs1 * P3 * P4) as f64) * (Y as f64).powf(2.0_f64)
        + ((vp2 * vs2 * P1 * P2) as f64) * (X as f64).powf(2.0_f64)
        + ((vp1 * vp2 * vs1 * vs2) as f64) * (p as f64).powf(2.0_f64) * (Z as f64).powf(2.0_f64))
        as f32);
    // the coefficient
    if wave_type == 0_i32 {
        COEF = (((((q * p) as f64).powf(2.0_f64)
            * (P1 as f64)
            * (P2 as f64)
            * (P3 as f64)
            * (P4 as f64)
            + ((rho1 * rho2 * (vs1 * vp2 * P1 * P4 - vp1 * vs2 * P2 * P3)) as f64)
            - ((vp1 * vs1 * P3 * P4) as f64) * (Y as f64).powf(2.0_f64)
            + ((vp2 * vs2 * P1 * P2) as f64) * (X as f64).powf(2.0_f64)
            - ((vp1 * vp2 * vs1 * vs2) as f64)
                * (p as f64).powf(2.0_f64)
                * (Z as f64).powf(2.0_f64))
            / (D as f64)) as f32);
    } else {
        if wave_type == 1_i32 {
            COEF = (-2_i32 as f32) * vp1 * p * P1 * (q * P3 * P4 * Y + vp2 * vs2 * X * Z) / D;
        } else {
            COEF = 1000.0_f32;
        }
    }
    return COEF;
}

/// Cosine and sine of the S-wave polarization rotation.
/// Halfspace kinds are 0 (ISO), 1 (VTI), 2 (HTI), and 3 (ORT).
/// Returns 0, -1 (completed at a singularity), -2 (ORT singularity), or
/// -100 (unknown kind), matching C.
///
/// C prototype: `int polar(float ang, float azim, float * cphi, float * sphi, int i_hsp);`
/// The Rust `model` argument supplies the C external medium parameters.
pub fn polar(
    model: &Model,
    ang: f32,
    azim: f32,
    cphi: &mut f32,
    sphi: &mut f32,
    i_hsp: i32,
) -> i32 {
    let vp1 = model.vp1;
    let vs1 = model.vs1;
    let eps1_1 = model.eps1_1;
    let eps1_2 = model.eps1_2;
    let delta1_1 = model.delta1_1;
    let delta1_2 = model.delta1_2;
    let delta1_3 = model.delta1_3;
    let gamma1_1 = model.gamma1_1;
    let gamma1_44 = model.gamma1_44;
    let gamma1_55 = model.gamma1_55;
    let mut c: i32;
    let mut perr: i32 = 0_i32;
    let mut B12: f32;
    let mut B11_B22: f32;
    let mut D: f32;
    let mut A55Gs: f32;
    let mut A55G1: f32;
    let mut A33: f32;
    let mut VS: f32;
    let mut sin_j: f32;
    let mut cos_j: f32;
    let mut sin_k: f32;
    let mut cos_k: f32;
    // do not do anything, if ISO inc. halfspace
    if i_hsp == 0_i32 {
        *cphi = 1.0_f32;
        *sphi = 0.0_f32;
        return {
            perr = 0_i32;
            perr
        };
    }
    // auxiliary quantities
    A55Gs = (((vs1 as f64).powf(2.0_f64) * ((gamma1_44 - gamma1_55) as f64)) as f32);
    A55G1 = (((vs1 as f64).powf(2.0_f64)
        * ((1.0_f32 + 2.0_f32 * gamma1_55) as f64)
        * (gamma1_1 as f64)) as f32);
    A33 = ((vp1 as f64).powf(2.0_f64) as f32);
    // get the reflection angle of the iso S-wave (in rad), and azimuth angles
    VS = (((2.0_f64 * (vs1 as f64).powf(2.0_f64) * ((1.0_f32 + 2.0_f32 * gamma1_44) as f64)
        + (vs1 as f64).powf(2.0_f64)
            * ((1.0_f32 + 2.0_f32 * gamma1_55) as f64)
            * ((1.0_f32 + gamma1_1) as f64))
        / 3.0_f64)
        .sqrt() as f32);
    sin_j = ((((VS / vp1) as f64) * (ang as f64).sin()) as f32);
    cos_j = ((1.0_f64 - 0.5_f64 * (sin_j as f64).powf(2.0_f64)) as f32);
    cos_k = ((azim as f64).cos() as f32);
    sin_k = ((azim as f64).sin() as f32);
    // get the polarization coefficients
    B11_B22 = ((((2.0_f32 * A33) as f64)
        * ((eps1_2 as f64) * (cos_k as f64).powf(2.0_f64) - (delta1_2 as f64))
        * (cos_k as f64).powf(2.0_f64)
        * (cos_j as f64).powf(2.0_f64)
        * (sin_j as f64).powf(2.0_f64)
        + ((2.0_f32 * A33) as f64)
            * ((eps1_1 as f64) * (sin_k as f64).powf(2.0_f64) - (delta1_1 as f64))
            * (sin_k as f64).powf(2.0_f64)
            * (cos_j as f64).powf(2.0_f64)
            * (sin_j as f64).powf(2.0_f64)
        + ((2.0_f32 * A33 * (1.0_f32 + 2.0_f32 * eps1_2) * delta1_3) as f64)
            * (cos_k as f64).powf(2.0_f64)
            * (sin_k as f64).powf(2.0_f64)
            * (sin_j as f64).powf(2.0_f64)
            * ((cos_j as f64).powf(2.0_f64) + 1.0_f64)
        + ((4.0_f32 * A33 * eps1_2) as f64)
            * (cos_k as f64).powf(2.0_f64)
            * (sin_k as f64).powf(2.0_f64)
            * (cos_j as f64).powf(2.0_f64)
            * (sin_j as f64).powf(2.0_f64)
        + ((2.0_f32 * A33 * (eps1_2 - eps1_1)) as f64)
            * (cos_k as f64).powf(2.0_f64)
            * (sin_k as f64).powf(2.0_f64)
            * (sin_j as f64).powf(2.0_f64)
        - ((2.0_f32 * A55G1) as f64) * (sin_j as f64).powf(2.0_f64)
        - ((2.0_f32 * A55Gs) as f64)
            * ((cos_k as f64).powf(2.0_f64) - (sin_k as f64).powf(2.0_f64))
        + ((2.0_f32 * A55Gs) as f64) * (cos_k as f64).powf(2.0_f64) * (sin_j as f64).powf(2.0_f64))
        as f32);
    B12 = ((((A33 * (1.0_f32 + 2.0_f32 * eps1_2) * delta1_3 * cos_k * sin_k) as f64)
        * ((cos_k as f64).powf(2.0_f64) - (sin_k as f64).powf(2.0_f64))
        * (cos_j as f64)
        * (sin_j as f64).powf(2.0_f64)
        + ((A33 * (delta1_2 - delta1_1) * cos_k * sin_k * cos_j) as f64)
            * (sin_j as f64).powf(2.0_f64)
        + ((2.0_f32 * A33 * (eps1_1 - eps1_2) * cos_k) as f64)
            * (sin_k as f64).powf(3.0_f64)
            * (cos_j as f64)
            * (sin_j as f64).powf(2.0_f64)
        + ((2.0_f32 * A55Gs * cos_k * sin_k * cos_j) as f64)) as f32);
    D = (((B11_B22 as f64).powf(2.0_f64) + 4.0_f64 * (B12 as f64).powf(2.0_f64)).sqrt() as f32);
    if D <= 1.1920929e-07_f32 {
        perr = -1_i32;
    }
    // do the logic for different halfspaces
    if i_hsp == 1_i32 {
        *cphi = 1.0_f32;
        *sphi = 0.0_f32;
        return perr;
    }
    if i_hsp == 2_i32 {
        *cphi = ((((cos_j * cos_k) as f64)
            / (1.0_f64 - ((sin_j * cos_k) as f64).powf(2.0_f64)).sqrt()) as f32);
        *sphi =
            (((-sin_k as f64) / (1.0_f64 - ((sin_j * cos_k) as f64).powf(2.0_f64)).sqrt()) as f32);
        return perr;
    }
    if i_hsp == 3_i32 {
        if perr == -1_i32 {
            return {
                perr = -2_i32;
                perr
            };
        }
        if (B12 as f64).abs() <= 1.1920929e-07_f64 {
            if (B11_B22 as f64).abs() <= 1.1920929e-07_f64 {
                return {
                    perr = -2_i32;
                    perr
                };
            }
            match {
                c = (((B11_B22 as f64) / (B11_B22 as f64).abs()).floor() as i32);
                c
            } {
                1 | 2 => {
                    *cphi = 1.0;
                    *sphi = 0.0;
                }
                0 | -1 | -2 => {
                    *cphi = 0.0;
                    *sphi = 1.0;
                }
                _ => {
                    perr = -2;
                    return perr;
                }
            }
        } else {
            *cphi = ((0.5_f64 * ((1.0_f32 + B11_B22 / D) as f64)).sqrt() as f32);
            *sphi = (((B12 as f64) / (B12 as f64).abs()
                * (0.5_f64 * ((1.0_f32 - B11_B22 / D) as f64)).sqrt()) as f32);
        }
        return 0_i32;
    }
    return -100_i32;
}

/// Evaluate first-order and empirical semi-second-order SV/SH error estimates.
/// Call with index 1, then index 2 for rotated media, retaining the same outputs.
/// Invalid indices return -1 without changing the outputs.
///
/// C prototype: `int S_err_2nd_order(ErrorFlag * rsv_1st, ErrorFlag * rsv_2nd, ErrorFlag * rsh_1st, ErrorFlag * rsh_2nd, float true_kappa, int index);`
/// The Rust `model` argument supplies the C external medium parameters.
pub fn S_err_2nd_order(
    model: &Model,
    rsv_1st: &mut ErrorFlag,
    rsv_2nd: &mut ErrorFlag,
    rsh_1st: &mut ErrorFlag,
    rsh_2nd: &mut ErrorFlag,
    true_kappa: f32,
    index: i32,
) -> i32 {
    let vp1 = model.vp1;
    let vp2 = model.vp2;
    let vs1 = model.vs1;
    let vs2 = model.vs2;
    let rho1 = model.rho1;
    let rho2 = model.rho2;
    let eps1_1 = model.eps1_1;
    let eps1_2 = model.eps1_2;
    let delta1_1 = model.delta1_1;
    let delta1_2 = model.delta1_2;
    let delta1_3 = model.delta1_3;
    let gamma1_44 = model.gamma1_44;
    let eps2_1 = model.eps2_1;
    let eps2_2 = model.eps2_2;
    let delta2_1 = model.delta2_1;
    let delta2_2 = model.delta2_2;
    let delta2_3 = model.delta2_3;
    let gamma2_44 = model.gamma2_44;
    let mut Da_a: f32;
    let mut Db_b: f32;
    let mut Dr_r: f32;
    let mut DG_G: f32;
    let mut alpha: f32;
    let mut beta: f32;
    let mut Delta1_1: f32;
    let mut Delta1_2: f32;
    let mut Delta2_1: f32;
    let mut Delta2_2: f32;
    let mut Gamma1_44: f32;
    let mut Gamma2_44: f32;
    let mut D1_1: f32;
    let mut D1_2: f32;
    let mut D2_1: f32;
    let mut D2_2: f32;
    let mut GS_1: f32;
    let mut GS_2: f32;
    let mut SHD1_1: f32;
    let mut SHD1_2: f32;
    let mut SHD2_1: f32;
    let mut SHD2_2: f32;
    let mut SHGS_1: f32;
    let mut SHGS_2: f32;
    let mut ISO_2nd: f32;
    let mut azim: f32;
    let mut d11_plus: f32 = 0.0_f32;
    let mut d11_min: f32 = 0.0_f32;
    let mut d12_plus: f32 = 0.0_f32;
    let mut d12_min: f32 = 0.0_f32;
    let mut gs1_plus: f32 = 0.0_f32;
    let mut gs1_min: f32 = 0.0_f32;
    let mut d21_plus: f32 = 0.0_f32;
    let mut d21_min: f32 = 0.0_f32;
    let mut d22_plus: f32 = 0.0_f32;
    let mut d22_min: f32 = 0.0_f32;
    let mut gs2_plus: f32 = 0.0_f32;
    let mut gs2_min: f32 = 0.0_f32;
    let mut a_plus: f32 = 0.0_f32;
    let mut a_min: f32 = 0.0_f32;
    let mut b_plus: f32 = 0.0_f32;
    let mut b_min: f32 = 0.0_f32;
    let mut r_plus: f32 = 0.0_f32;
    let mut r_min: f32 = 0.0_f32;
    let mut Fst0_0: f32;
    let mut Fst0_90: f32;
    let mut Fst90_0: f32;
    let mut Fst90_90: f32;
    let mut FstK_0: f32;
    let mut FstK_90: f32;
    let mut FstKK_0: f32;
    let mut FstKK_90: f32;
    let mut Fst0_kappa: f32;
    let mut Fst90_kappa: f32;
    let mut FstK_kappa: f32;
    let mut FstKK_kappa: f32;
    let mut Snd0_0: f32;
    let mut Snd0_90: f32;
    let mut Snd90_0: f32;
    let mut Snd90_90: f32;
    let mut SndK_0: f32;
    let mut SndK_90: f32;
    let mut SndKK_0: f32;
    let mut SndKK_90: f32;
    let mut Snd0_kappa: f32;
    let mut Snd90_kappa: f32;
    let mut SndK_kappa: f32;
    let mut SndKK_kappa: f32;
    let mut SHFst20_0: f32;
    let mut SHFst20_90: f32;
    let mut SHFst45_0: f32;
    let mut SHFst45_90: f32;
    let mut SHFst70_0: f32;
    let mut SHFst70_90: f32;
    let mut SHFst20_kappa: f32;
    let mut SHFst45_kappa: f32;
    let mut SHFst70_kappa: f32;
    let mut SHSnd20_0: f32;
    let mut SHSnd20_90: f32;
    let mut SHSnd45_0: f32;
    let mut SHSnd45_90: f32;
    let mut SHSnd70_0: f32;
    let mut SHSnd70_90: f32;
    let mut SHSnd20_kappa: f32;
    let mut SHSnd45_kappa: f32;
    let mut SHSnd70_kappa: f32;
    let mut err_kappa: f32 = 0.0_f32;
    let lkappa: f32 = ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    let kkappa: f32 = ((60.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    let mut D: f32;
    let mut N: i32 = 0_i32;
    // some quantities in advance
    Da_a = ((((vp2 - vp1) as f64) / (0.5_f64 * ((vp2 + vp1) as f64))) as f32);
    Db_b = ((((vs2 - vs1) as f64) / (0.5_f64 * ((vs2 + vs1) as f64))) as f32);
    Dr_r = ((((rho2 - rho1) as f64) / (0.5_f64 * ((rho2 + rho1) as f64))) as f32);
    DG_G = ((((vs2 as f64).powf(2.0_f64) * (rho2 as f64)
        - (vs1 as f64).powf(2.0_f64) * (rho1 as f64))
        / (0.5_f64
            * ((vs2 as f64).powf(2.0_f64) * (rho2 as f64)
                + (vs1 as f64).powf(2.0_f64) * (rho1 as f64)))) as f32);
    alpha = ((0.5_f64 * ((vp2 + vp1) as f64)) as f32);
    beta = ((0.5_f64 * ((vs2 + vs1) as f64)) as f32);
    if index == 1_i32 {
        err_kappa = 0.0_f32;
    }
    if index == 2_i32 {
        err_kappa = ((3.141592653589793_f64 / 2.0_f64) as f32);
    }
    if (index != 1_i32 && index != 2_i32) {
        return -1_i32;
    }
    if (delta1_1 != 0.0) {
        d11_plus = ((((delta1_1 as f64).abs() + (delta1_1 as f64))
            / (2.0_f64 * (delta1_1 as f64).abs())) as f32);
    }
    if (delta1_1 != 0.0) {
        d11_min = ((((delta1_1 as f64).abs() - (delta1_1 as f64))
            / (2.0_f64 * (delta1_1 as f64).abs())) as f32);
    }
    if (delta1_2 != 0.0) {
        d12_plus = ((((delta1_2 as f64).abs() + (delta1_2 as f64))
            / (2.0_f64 * (delta1_2 as f64).abs())) as f32);
    }
    if (delta1_2 != 0.0) {
        d12_min = ((((delta1_2 as f64).abs() - (delta1_2 as f64))
            / (2.0_f64 * (delta1_2 as f64).abs())) as f32);
    }
    if (gamma1_44 != 0.0) {
        gs1_plus = ((((gamma1_44 as f64).abs() + (gamma1_44 as f64))
            / (2.0_f64 * (gamma1_44 as f64).abs())) as f32);
    }
    if (gamma1_44 != 0.0) {
        gs1_min = ((((gamma1_44 as f64).abs() - (gamma1_44 as f64))
            / (2.0_f64 * (gamma1_44 as f64).abs())) as f32);
    }
    if (delta2_1 != 0.0) {
        d21_plus = ((((delta2_1 as f64).abs() + (delta2_1 as f64))
            / (2.0_f64 * (delta2_1 as f64).abs())) as f32);
    }
    if (delta2_1 != 0.0) {
        d21_min = ((((delta2_1 as f64).abs() - (delta2_1 as f64))
            / (2.0_f64 * (delta2_1 as f64).abs())) as f32);
    }
    if (delta2_2 != 0.0) {
        d22_plus = ((((delta2_2 as f64).abs() + (delta2_2 as f64))
            / (2.0_f64 * (delta2_2 as f64).abs())) as f32);
    }
    if (delta2_2 != 0.0) {
        d22_min = ((((delta2_2 as f64).abs() - (delta2_2 as f64))
            / (2.0_f64 * (delta2_2 as f64).abs())) as f32);
    }
    if (gamma2_44 != 0.0) {
        gs2_plus = ((((gamma2_44 as f64).abs() + (gamma2_44 as f64))
            / (2.0_f64 * (gamma2_44 as f64).abs())) as f32);
    }
    if (gamma2_44 != 0.0) {
        gs2_min = ((((gamma2_44 as f64).abs() - (gamma2_44 as f64))
            / (2.0_f64 * (gamma2_44 as f64).abs())) as f32);
    }
    if (Da_a != 0.0) {
        a_plus = ((((Da_a as f64).abs() + (Da_a as f64)) / (2.0_f64 * (Da_a as f64).abs())) as f32);
    }
    if (Da_a != 0.0) {
        a_min = ((((Da_a as f64).abs() - (Da_a as f64)) / (2.0_f64 * (Da_a as f64).abs())) as f32);
    }
    if (Db_b != 0.0) {
        b_plus = ((((Db_b as f64).abs() + (Db_b as f64)) / (2.0_f64 * (Db_b as f64).abs())) as f32);
    }
    if (Db_b != 0.0) {
        b_min = ((((Db_b as f64).abs() - (Db_b as f64)) / (2.0_f64 * (Db_b as f64).abs())) as f32);
    }
    if (Dr_r != 0.0) {
        r_plus = ((((Dr_r as f64).abs() + (Dr_r as f64)) / (2.0_f64 * (Dr_r as f64).abs())) as f32);
    }
    if (Dr_r != 0.0) {
        r_min = ((((Dr_r as f64).abs() - (Dr_r as f64)) / (2.0_f64 * (Dr_r as f64).abs())) as f32);
    }
    // the angles used for the anisotropic accuracy evaluation
    (*rsv_2nd).angle[0_i32 as usize] = ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    (*rsv_1st).angle[0_i32 as usize] = (*rsv_2nd).angle[0_i32 as usize];
    // incidence angle
    (*rsv_2nd).angle[1_i32 as usize] = 0.0_f32;
    (*rsv_1st).angle[1_i32 as usize] = (*rsv_2nd).angle[1_i32 as usize];
    // first azimuth
    (*rsv_2nd).angle[2_i32 as usize] = 90.0_f32;
    (*rsv_1st).angle[2_i32 as usize] = (*rsv_2nd).angle[2_i32 as usize];
    // second azimuth
    (*rsv_2nd).angle[3_i32 as usize] = true_kappa;
    (*rsv_1st).angle[3_i32 as usize] = (*rsv_2nd).angle[3_i32 as usize];
    // kappa
    (*rsh_2nd).angle[0_i32 as usize] = ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
    (*rsh_1st).angle[0_i32 as usize] = (*rsh_2nd).angle[0_i32 as usize];
    // incidence angle
    (*rsh_2nd).angle[1_i32 as usize] = 0.0_f32;
    (*rsh_1st).angle[1_i32 as usize] = (*rsh_2nd).angle[1_i32 as usize];
    // first azimuth
    (*rsh_2nd).angle[2_i32 as usize] = 90.0_f32;
    (*rsh_1st).angle[2_i32 as usize] = (*rsh_2nd).angle[2_i32 as usize];
    // second azimuth
    (*rsh_2nd).angle[3_i32 as usize] = true_kappa;
    (*rsh_1st).angle[3_i32 as usize] = (*rsh_2nd).angle[3_i32 as usize];
    // kappa
    // ISO part first
    // SV
    N = 0_i32;
    while N < 5_i32 {
        (*rsv_1st).iso[N as usize] = SV_term(
            DG_G,
            Dr_r,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            alpha,
            beta,
            ((((15_i32 + 5_i32 * N) as f64) * 3.141592653589793_f64 / 180.0_f64) as f32),
            0.0_f32,
            0.0_f32,
        );
        (*rsv_2nd).iso[N as usize] = Iso_exact(
            1_i32,
            vp1,
            vs1,
            rho1,
            vp2,
            vs2,
            rho2,
            ((((15_i32 + 5_i32 * N) as f64) * 3.141592653589793_f64 / 180.0_f64) as f32),
        );
        N += 1;
    }
    // evaluate the ISO 2nd order for the inc. angle err_ang
    ISO_2nd = Iso_exact(
        1_i32,
        vp1,
        vs1,
        rho1,
        vp2,
        vs2,
        rho2,
        ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
    ) - SV_term(
        DG_G,
        Dr_r,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        0.0_f32,
        alpha,
        beta,
        ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
        0.0_f32,
        0.0_f32,
    );
    // SH
    (*rsh_2nd).iso[0_i32 as usize] = 0.0_f32;
    (*rsh_1st).iso[0_i32 as usize] = (*rsh_2nd).iso[0_i32 as usize];
    (*rsh_2nd).iso[1_i32 as usize] = 0.0_f32;
    (*rsh_1st).iso[1_i32 as usize] = (*rsh_2nd).iso[1_i32 as usize];
    (*rsh_2nd).iso[2_i32 as usize] = 0.0_f32;
    (*rsh_1st).iso[2_i32 as usize] = (*rsh_2nd).iso[2_i32 as usize];
    (*rsh_2nd).iso[3_i32 as usize] = 0.0_f32;
    (*rsh_1st).iso[3_i32 as usize] = (*rsh_2nd).iso[3_i32 as usize];
    (*rsh_2nd).iso[4_i32 as usize] = 0.0_f32;
    (*rsh_1st).iso[4_i32 as usize] = (*rsh_2nd).iso[4_i32 as usize];
    // UPPER halfspace contribution
    // SV
    Delta1_1 = ((((((-1_i32 as f32) * d11_plus) as f64) - 2.5_f64 * (d11_min as f64))
        * (delta1_1 as f64).powf(2.0_f64)
        + (-0.22_f64 * (d11_plus as f64) * (gs1_plus as f64)
            + ((0.0_f32 * d11_plus * gs1_min) as f64)
            - 0.5_f64 * (d11_min as f64) * (gs1_plus as f64)
            - 0.39_f64 * (d11_min as f64) * (gs1_min as f64))
            * (delta1_1 as f64)
            * (gamma1_44 as f64)
        + (1.28_f64 * (d11_plus as f64) * (a_plus as f64)
            + 0.75_f64 * (d11_plus as f64) * (a_min as f64)
            + 1.93_f64 * (d11_min as f64) * (a_plus as f64)
            + 2.93_f64 * (d11_min as f64) * (a_min as f64))
            * (delta1_1 as f64)
            * (Da_a as f64)
        + (0.64_f64 * (d11_plus as f64) * (b_plus as f64)
            + 0.23_f64 * (d11_plus as f64) * (b_min as f64)
            + 0.18_f64 * (d11_min as f64) * (b_plus as f64)
            + 1.87_f64 * (d11_min as f64) * (b_min as f64))
            * (delta1_1 as f64)
            * (Db_b as f64)
        + (1.15_f64 * (d11_plus as f64) * (r_plus as f64)
            + 0.66_f64 * (d11_plus as f64) * (r_min as f64)
            + 3.51_f64 * (d11_min as f64) * (r_plus as f64)
            + 2.69_f64 * (d11_min as f64) * (r_min as f64))
            * (delta1_1 as f64)
            * (Dr_r as f64)) as f32);
    D1_1 = Delta1_1;
    Delta1_2 = ((((((-1_i32 as f32) * d12_plus) as f64) - 2.5_f64 * (d12_min as f64))
        * (delta1_2 as f64).powf(2.0_f64)
        + (1.28_f64 * (d12_plus as f64) * (a_plus as f64)
            + 0.75_f64 * (d12_plus as f64) * (a_min as f64)
            + 1.93_f64 * (d12_min as f64) * (a_plus as f64)
            + 2.93_f64 * (d12_min as f64) * (a_min as f64))
            * (delta1_2 as f64)
            * (Da_a as f64)
        + (0.64_f64 * (d12_plus as f64) * (b_plus as f64)
            + 0.23_f64 * (d12_plus as f64) * (b_min as f64)
            + 0.18_f64 * (d12_min as f64) * (b_plus as f64)
            + 1.87_f64 * (d12_min as f64) * (b_min as f64))
            * (delta1_2 as f64)
            * (Db_b as f64)
        + (1.15_f64 * (d12_plus as f64) * (r_plus as f64)
            + 0.66_f64 * (d12_plus as f64) * (r_min as f64)
            + 3.51_f64 * (d12_min as f64) * (r_plus as f64)
            + 2.69_f64 * (d12_min as f64) * (r_min as f64))
            * (delta1_2 as f64)
            * (Dr_r as f64)) as f32);
    D1_2 = Delta1_2;
    Delta2_1 = 0.0_f32;
    Delta2_2 = 0.0_f32;
    Gamma1_44 = (((-0.6_f64 * (gs1_plus as f64) - 1.1_f64 * (gs1_min as f64))
        * (gamma1_44 as f64).powf(2.0_f64)
        + (-0.385_f64 * (gs1_plus as f64) * (a_plus as f64)
            - 0.36_f64 * (gs1_plus as f64) * (a_min as f64)
            - 0.55_f64 * (gs1_min as f64) * (a_plus as f64)
            - 0.225_f64 * (gs1_min as f64) * (a_min as f64))
            * (gamma1_44 as f64)
            * (Da_a as f64)
        + (((0.0_f32 * gs1_plus * b_plus + 0.0_f32 * gs1_plus * b_min + 0.0_f32 * gs1_min * b_plus)
            as f64)
            - 0.3_f64 * (gs1_min as f64) * (b_min as f64))
            * (gamma1_44 as f64)
            * (Db_b as f64)
        + (-0.33_f64 * (gs1_plus as f64) * (r_plus as f64)
            - 0.32_f64 * (gs1_plus as f64) * (r_min as f64)
            - 0.64_f64 * (gs1_min as f64) * (r_plus as f64)
            - 0.57_f64 * (gs1_min as f64) * (r_min as f64))
            * (gamma1_44 as f64)
            * (Dr_r as f64)) as f32);
    GS_1 = Gamma1_44;
    Gamma2_44 = 0.0_f32;
    if index == 1_i32 {
        azim = 0.0_f32;
        (*rsv_1st).upper[0_i32 as usize] = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            eps1_1,
            eps1_2,
            0.0_f32,
            0.0_f32,
            gamma1_44,
            0.0_f32,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rsv_2nd).upper[0_i32 as usize] = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        ) + ISO_2nd;
        azim = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rsv_1st).upper[1_i32 as usize] = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            eps1_1,
            eps1_2,
            0.0_f32,
            0.0_f32,
            gamma1_44,
            0.0_f32,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rsv_2nd).upper[1_i32 as usize] = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        ) + ISO_2nd;
    }
    // SH
    Delta1_1 = (((-0.8_f64 * (d11_plus as f64) - 2.0_f64 * (d11_min as f64))
        * (delta1_1 as f64).powf(2.0_f64)
        + (-0.22_f64 * (d11_plus as f64) * (gs1_plus as f64)
            + ((0.0_f32 * d11_plus * gs1_min) as f64)
            - 0.4_f64 * (d11_min as f64) * (gs1_plus as f64)
            - 0.22_f64 * (d11_min as f64) * (gs1_min as f64))
            * (delta1_1 as f64)
            * (gamma1_44 as f64)
        + (((0.0_f32 * d11_plus * d12_plus
            + 0.0_f32 * d11_plus * d12_min
            + 0.0_f32 * d11_min * d12_plus
            + 0.0_f32 * d11_min * d12_min)
            * delta1_1
            * delta1_2) as f64)
        + (1.32_f64 * (d11_plus as f64) * (a_plus as f64)
            + 1.09_f64 * (d11_plus as f64) * (a_min as f64)
            + 2.93_f64 * (d11_min as f64) * (a_plus as f64)
            + 2.4_f64 * (d11_min as f64) * (a_min as f64))
            * (delta1_1 as f64)
            * (Da_a as f64)
        + (0.95_f64 * (d11_plus as f64) * (b_plus as f64)
            + 0.42_f64 * (d11_plus as f64) * (b_min as f64)
            + 0.73_f64 * (d11_min as f64) * (b_plus as f64)
            + 0.75_f64 * (d11_min as f64) * (b_min as f64))
            * (delta1_1 as f64)
            * (Db_b as f64)
        + (1.67_f64 * (d11_plus as f64) * (r_plus as f64)
            + 1.68_f64 * (d11_plus as f64) * (r_min as f64)
            + 3.71_f64 * (d11_min as f64) * (r_plus as f64)
            + 3.12_f64 * (d11_min as f64) * (r_min as f64))
            * (delta1_1 as f64)
            * (Dr_r as f64)) as f32);
    SHD1_1 = Delta1_1;
    Delta1_2 = (((-0.8_f64 * (d12_plus as f64) - 2.0_f64 * (d12_min as f64))
        * (delta1_2 as f64).powf(2.0_f64)
        + (0.22_f64 * (d12_plus as f64) * (gs1_plus as f64)
            - 0.33_f64 * (d12_plus as f64) * (gs1_min as f64)
            + 0.22_f64 * (d12_min as f64) * (gs1_plus as f64)
            + 0.32_f64 * (d12_min as f64) * (gs1_min as f64))
            * (delta1_2 as f64)
            * (gamma1_44 as f64)
        + (1.32_f64 * (d12_plus as f64) * (a_plus as f64)
            + 1.09_f64 * (d12_plus as f64) * (a_min as f64)
            + 2.93_f64 * (d12_min as f64) * (a_plus as f64)
            + 2.4_f64 * (d12_min as f64) * (a_min as f64))
            * (delta1_2 as f64)
            * (Da_a as f64)
        + (0.95_f64 * (d12_plus as f64) * (b_plus as f64)
            + 0.42_f64 * (d12_plus as f64) * (b_min as f64)
            + 0.73_f64 * (d12_min as f64) * (b_plus as f64)
            + 0.75_f64 * (d12_min as f64) * (b_min as f64))
            * (delta1_2 as f64)
            * (Db_b as f64)
        + (1.67_f64 * (d12_plus as f64) * (r_plus as f64)
            + 1.68_f64 * (d12_plus as f64) * (r_min as f64)
            + 3.71_f64 * (d12_min as f64) * (r_plus as f64)
            + 3.12_f64 * (d12_min as f64) * (r_min as f64))
            * (delta1_2 as f64)
            * (Dr_r as f64)) as f32);
    SHD1_2 = Delta1_2;
    Delta2_1 = 0.0_f32;
    Delta2_2 = 0.0_f32;
    Gamma1_44 = (((-0.45_f64 * (gs1_plus as f64) - 0.78_f64 * (gs1_min as f64))
        * (gamma1_44 as f64).powf(2.0_f64)
        + (-0.17_f64 * (gs1_plus as f64) * (a_plus as f64)
            - 0.16_f64 * (gs1_plus as f64) * (a_min as f64)
            - 0.28_f64 * (gs1_min as f64) * (a_plus as f64)
            - 0.09_f64 * (gs1_min as f64) * (a_min as f64))
            * (gamma1_44 as f64)
            * (Da_a as f64)
        + (-0.27_f64 * (gs1_plus as f64) * (b_plus as f64)
            - 0.45_f64 * (gs1_plus as f64) * (b_min as f64)
            - 0.55_f64 * (gs1_min as f64) * (b_plus as f64)
            - 0.52_f64 * (gs1_min as f64) * (b_min as f64))
            * (gamma1_44 as f64)
            * (Db_b as f64)
        + (-0.38_f64 * (gs1_plus as f64) * (r_plus as f64)
            - 0.48_f64 * (gs1_plus as f64) * (r_min as f64)
            - 0.57_f64 * (gs1_min as f64) * (r_plus as f64)
            - 0.61_f64 * (gs1_min as f64) * (r_min as f64))
            * (gamma1_44 as f64)
            * (Dr_r as f64)) as f32);
    SHGS_1 = Gamma1_44;
    Gamma2_44 = 0.0_f32;
    if index == 1_i32 {
        azim = ((45.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rsh_1st).upper[0_i32 as usize] = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            eps1_1,
            eps1_2,
            0.0_f32,
            0.0_f32,
            gamma1_44,
            0.0_f32,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rsh_2nd).upper[0_i32 as usize] = SH_term(
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rsh_2nd).upper[1_i32 as usize] = 0.0_f32;
        (*rsh_1st).upper[1_i32 as usize] = 0.0_f32;
    }
    // LOWER halfspace contribution
    // SV
    Delta1_1 = 0.0_f32;
    Delta1_2 = 0.0_f32;
    Delta2_1 = (((-0.72_f64 * (d21_plus as f64) - 1.89_f64 * (d21_min as f64))
        * (delta2_1 as f64).powf(2.0_f64)
        + (0.61_f64 * (d21_plus as f64) * (gs2_plus as f64)
            + 1.26_f64 * (d21_plus as f64) * (gs2_min as f64)
            + 0.83_f64 * (d21_min as f64) * (gs2_plus as f64)
            + 2.28_f64 * (d21_min as f64) * (gs2_min as f64))
            * (delta2_1 as f64)
            * (gamma2_44 as f64)
        + (-0.28_f64 * (d21_plus as f64) * (a_plus as f64) + ((0.0_f32 * d21_plus * a_min) as f64)
            - 1.1_f64 * (d21_min as f64) * (a_plus as f64)
            - 1.05_f64 * (d21_min as f64) * (a_min as f64))
            * (delta2_1 as f64)
            * (Da_a as f64)
        + (0.82_f64 * (d21_plus as f64) * (b_plus as f64)
            + 0.83_f64 * (d21_plus as f64) * (b_min as f64)
            + 2.93_f64 * (d21_min as f64) * (b_plus as f64)
            + 1.73_f64 * (d21_min as f64) * (b_min as f64))
            * (delta2_1 as f64)
            * (Db_b as f64)
        + (-0.07_f64 * (d21_plus as f64) * (r_plus as f64)
            + 0.1_f64 * (d21_plus as f64) * (r_min as f64)
            + ((0.0_f32 * d21_min * r_plus) as f64)
            + ((0.0_f32 * d21_min * r_min) as f64))
            * (delta2_1 as f64)
            * (Dr_r as f64)) as f32);
    D2_1 = Delta2_1;
    Delta2_2 = (((-0.72_f64 * (d22_plus as f64) - 1.89_f64 * (d22_min as f64))
        * (delta2_2 as f64).powf(2.0_f64)
        + (-0.28_f64 * (d22_plus as f64) * (a_plus as f64) + ((0.0_f32 * d22_plus * a_min) as f64)
            - 1.1_f64 * (d22_min as f64) * (a_plus as f64)
            - 1.05_f64 * (d22_min as f64) * (a_min as f64))
            * (delta2_2 as f64)
            * (Da_a as f64)
        + (0.82_f64 * (d22_plus as f64) * (b_plus as f64)
            + 0.83_f64 * (d22_plus as f64) * (b_min as f64)
            + 2.93_f64 * (d22_min as f64) * (b_plus as f64)
            + 1.73_f64 * (d22_min as f64) * (b_min as f64))
            * (delta2_2 as f64)
            * (Db_b as f64)
        + (-0.07_f64 * (d22_plus as f64) * (r_plus as f64)
            + 0.1_f64 * (d22_plus as f64) * (r_min as f64)
            + ((0.0_f32 * d22_min * r_plus) as f64)
            + ((0.0_f32 * d22_min * r_min) as f64))
            * (delta2_2 as f64)
            * (Dr_r as f64)) as f32);
    D2_2 = Delta2_2;
    Gamma1_44 = 0.0_f32;
    Gamma2_44 = (((-0.28_f64 * (gs2_plus as f64) - 2.0_f64 / 3.0_f64 * (gs2_min as f64))
        * (gamma2_44 as f64).powf(2.0_f64)
        + (-0.37_f64 * (gs2_plus as f64) * (a_plus as f64)
            - 0.3_f64 * (gs2_plus as f64) * (a_min as f64)
            - 0.55_f64 * (gs2_min as f64) * (a_plus as f64)
            - 0.45_f64 * (gs2_min as f64) * (a_min as f64))
            * (gamma2_44 as f64)
            * (Da_a as f64)
        + (1.1_f64 * (gs2_plus as f64) * (b_plus as f64)
            + 0.69_f64 * (gs2_plus as f64) * (b_min as f64)
            + 0.92_f64 * (gs2_min as f64) * (b_plus as f64)
            + 0.6_f64 * (gs2_min as f64) * (b_min as f64))
            * (gamma2_44 as f64)
            * (Db_b as f64)
        + (((0.0_f32 * gs2_plus * r_plus + 0.0_f32 * gs2_plus * r_min) as f64)
            + 0.06_f64 * (gs2_min as f64) * (r_plus as f64)
            + 0.2_f64 * (gs2_min as f64) * (r_min as f64))
            * (gamma2_44 as f64)
            * (Dr_r as f64)) as f32);
    GS_2 = Gamma2_44;
    if index == 1_i32 {
        azim = 0.0_f32;
        (*rsv_1st).lower[0_i32 as usize] = SV_term(
            DG_G,
            Dr_r,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            delta2_1,
            delta2_2,
            delta2_3,
            0.0_f32,
            0.0_f32,
            eps2_1,
            eps2_2,
            0.0_f32,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rsv_2nd).lower[0_i32 as usize] = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        ) + ISO_2nd;
        azim = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rsv_1st).lower[1_i32 as usize] = SV_term(
            DG_G,
            Dr_r,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            delta2_1,
            delta2_2,
            delta2_3,
            0.0_f32,
            0.0_f32,
            eps2_1,
            eps2_2,
            0.0_f32,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rsv_2nd).lower[1_i32 as usize] = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        ) + ISO_2nd;
    }
    // SH
    Delta1_1 = 0.0_f32;
    Delta1_2 = 0.0_f32;
    Delta2_1 = (((-0.61_f64 * (d21_plus as f64) - 1.56_f64 * (d21_min as f64))
        * (delta2_1 as f64).powf(2.0_f64)
        + (0.28_f64 * (d21_plus as f64) * (gs2_plus as f64)
            + 1.1_f64 * (d21_plus as f64) * (gs2_min as f64)
            + 0.44_f64 * (d21_min as f64) * (gs2_plus as f64)
            + 1.56_f64 * (d21_min as f64) * (gs2_min as f64))
            * (delta2_1 as f64)
            * (gamma2_44 as f64)
        + (((0.0_f32 * d21_plus * d22_plus
            + 0.0_f32 * d21_plus * d22_min
            + 0.0_f32 * d21_min * d22_plus
            + 0.0_f32 * d21_min * d22_min)
            * delta2_1
            * delta2_2) as f64)
        + (0.09_f64 * (d21_plus as f64) * (a_plus as f64)
            + 0.38_f64 * (d21_plus as f64) * (a_min as f64)
            - 0.47_f64 * (d21_min as f64) * (a_plus as f64)
            - 0.74_f64 * (d21_min as f64) * (a_min as f64))
            * (delta2_1 as f64)
            * (Da_a as f64)
        + (0.55_f64 * (d21_plus as f64) * (b_plus as f64)
            + 0.45_f64 * (d21_plus as f64) * (b_min as f64)
            + 2.01_f64 * (d21_min as f64) * (b_plus as f64)
            + 1.28_f64 * (d21_min as f64) * (b_min as f64))
            * (delta2_1 as f64)
            * (Db_b as f64)
        + (-0.07_f64 * (d21_plus as f64) * (r_plus as f64)
            + 0.05_f64 * (d21_plus as f64) * (r_min as f64)
            + ((0.0_f32 * d21_min * r_plus) as f64)
            + 0.14_f64 * (d21_min as f64) * (r_min as f64))
            * (delta2_1 as f64)
            * (Dr_r as f64)) as f32);
    SHD2_1 = Delta2_1;
    Delta2_2 = (((-0.61_f64 * (d22_plus as f64) - 1.56_f64 * (d22_min as f64))
        * (delta2_2 as f64).powf(2.0_f64)
        + (0.11_f64 * (d22_plus as f64) * (gs2_plus as f64)
            + 0.11_f64 * (d22_plus as f64) * (gs2_min as f64)
            + 0.0_f64 * (d22_min as f64) * (gs2_plus as f64)
            + 0.64_f64 * (d22_min as f64) * (gs2_min as f64))
            * (delta2_2 as f64)
            * (gamma2_44 as f64)
        + (0.09_f64 * (d22_plus as f64) * (a_plus as f64)
            + 0.38_f64 * (d22_plus as f64) * (a_min as f64)
            - 0.47_f64 * (d22_min as f64) * (a_plus as f64)
            - 0.74_f64 * (d22_min as f64) * (a_min as f64))
            * (delta2_2 as f64)
            * (Da_a as f64)
        + (0.55_f64 * (d22_plus as f64) * (b_plus as f64)
            + 0.45_f64 * (d22_plus as f64) * (b_min as f64)
            + 2.01_f64 * (d22_min as f64) * (b_plus as f64)
            + 1.28_f64 * (d22_min as f64) * (b_min as f64))
            * (delta2_2 as f64)
            * (Db_b as f64)
        + (-0.07_f64 * (d22_plus as f64) * (r_plus as f64)
            + 0.05_f64 * (d22_plus as f64) * (r_min as f64)
            + ((0.0_f32 * d22_min * r_plus) as f64)
            + 0.14_f64 * (d22_min as f64) * (r_min as f64))
            * (delta2_2 as f64)
            * (Dr_r as f64)) as f32);
    SHD2_2 = Delta2_2;
    Gamma1_44 = 0.0_f32;
    Gamma2_44 = (((-(1.0_f64 / 3.0_f64) * (gs2_plus as f64) - 2.0_f64 / 3.0_f64 * (gs2_min as f64))
        * (gamma2_44 as f64).powf(2.0_f64)
        + (-0.18_f64 * (gs2_plus as f64) * (a_plus as f64)
            - 0.08_f64 * (gs2_plus as f64) * (a_min as f64)
            - 0.18_f64 * (gs2_min as f64) * (a_plus as f64)
            - 0.15_f64 * (gs2_min as f64) * (a_min as f64))
            * (gamma2_44 as f64)
            * (Da_a as f64)
        + (0.64_f64 * (gs2_plus as f64) * (b_plus as f64)
            + 0.52_f64 * (gs2_plus as f64) * (b_min as f64)
            + 0.55_f64 * (gs2_min as f64) * (b_plus as f64)
            + 0.52_f64 * (gs2_min as f64) * (b_min as f64))
            * (gamma2_44 as f64)
            * (Db_b as f64)
        + (-0.06_f64 * (gs2_plus as f64) * (r_plus as f64)
            + ((0.0_f32 * gs2_plus * r_min) as f64)
            + 0.06_f64 * (gs2_min as f64) * (r_plus as f64)
            + 0.23_f64 * (gs2_min as f64) * (r_min as f64))
            * (gamma2_44 as f64)
            * (Dr_r as f64)) as f32);
    SHGS_2 = Gamma2_44;
    if index == 1_i32 {
        azim = ((45.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rsh_1st).lower[0_i32 as usize] = SH_term(
            0.0_f32,
            0.0_f32,
            0.0_f32,
            delta2_1,
            delta2_2,
            delta2_3,
            0.0_f32,
            0.0_f32,
            eps2_1,
            eps2_2,
            0.0_f32,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rsh_2nd).lower[0_i32 as usize] = SH_term(
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            0.0_f32,
        );
        (*rsh_1st).lower[1_i32 as usize] = 0.0_f32;
        (*rsh_2nd).lower[1_i32 as usize] = 0.0_f32;
    }
    // GLOBAL model contribution
    // SV
    // The following terms are derived on the assumption that Rps coefficients for
    // two identicaly anisotropic halfspaces are zero (with 0 velocity contrats)
    Delta1_1 = (((D1_1 as f64)
        + (1.0_f64 * (d11_plus as f64) * (d21_plus as f64)
            - 1.0_f64 * (d11_plus as f64) * (d21_min as f64)
            - 2.5_f64 * (d11_min as f64) * (d21_plus as f64)
            + 2.5_f64 * (d11_min as f64) * (d21_min as f64))
            * (delta1_1 as f64)
            * (delta2_1 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (1.0_f64 * (d11_plus as f64) * (d22_plus as f64)
            - 1.0_f64 * (d11_plus as f64) * (d22_min as f64)
            - 2.5_f64 * (d11_min as f64) * (d22_plus as f64)
            + 2.5_f64 * (d11_min as f64) * (d22_min as f64))
            * (delta1_1 as f64)
            * (delta2_2 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (0.22_f64 * (d11_plus as f64) * (gs2_plus as f64)
            + 0.0_f64 * (d11_plus as f64) * (gs2_min as f64)
            + 0.5_f64 * (d11_min as f64) * (gs2_plus as f64)
            + 0.39_f64 * (d11_min as f64) * (gs2_min as f64))
            * (delta1_1 as f64)
            * (gamma2_44 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)) as f32);
    Delta1_2 = (((D1_2 as f64)
        + (1.0_f64 * (d12_plus as f64) * (d21_plus as f64)
            - 1.0_f64 * (d12_plus as f64) * (d21_min as f64)
            - 2.5_f64 * (d12_min as f64) * (d21_plus as f64)
            + 2.5_f64 * (d12_min as f64) * (d21_min as f64))
            * (delta1_2 as f64)
            * (delta2_1 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (1.0_f64 * (d12_plus as f64) * (d22_plus as f64)
            - 1.0_f64 * (d12_plus as f64) * (d22_min as f64)
            - 2.5_f64 * (d12_min as f64) * (d22_plus as f64)
            + 2.5_f64 * (d12_min as f64) * (d22_min as f64))
            * (delta1_2 as f64)
            * (delta2_2 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (0.22_f64 * (d12_plus as f64) * (gs2_plus as f64)
            + 0.0_f64 * (d12_plus as f64) * (gs2_min as f64)
            + 0.5_f64 * (d12_min as f64) * (gs2_plus as f64)
            + 0.39_f64 * (d12_min as f64) * (gs2_min as f64))
            * (delta1_2 as f64)
            * (gamma2_44 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)) as f32);
    Delta2_1 = (((D2_1 as f64)
        + (0.72_f64 * (d21_plus as f64) * (d11_plus as f64)
            - 0.72_f64 * (d21_plus as f64) * (d11_min as f64)
            - 1.89_f64 * (d21_min as f64) * (d11_plus as f64)
            + 1.89_f64 * (d21_min as f64) * (d11_min as f64))
            * (delta2_1 as f64)
            * (delta1_1 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (0.72_f64 * (d21_plus as f64) * (d12_plus as f64)
            - 0.72_f64 * (d21_plus as f64) * (d12_min as f64)
            - 1.89_f64 * (d21_min as f64) * (d12_plus as f64)
            + 1.89_f64 * (d21_min as f64) * (d12_min as f64))
            * (delta2_1 as f64)
            * (delta1_2 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (-0.61_f64 * (d21_plus as f64) * (gs1_plus as f64)
            - 1.26_f64 * (d21_plus as f64) * (gs1_min as f64)
            - 0.83_f64 * (d21_min as f64) * (gs1_plus as f64)
            - 2.28_f64 * (d21_min as f64) * (gs1_min as f64))
            * (delta2_1 as f64)
            * (gamma1_44 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)) as f32);
    Delta2_2 = (((D2_2 as f64)
        + (0.72_f64 * (d22_plus as f64) * (d11_plus as f64)
            - 0.72_f64 * (d22_plus as f64) * (d11_min as f64)
            - 1.89_f64 * (d22_min as f64) * (d11_plus as f64)
            + 1.89_f64 * (d22_min as f64) * (d11_min as f64))
            * (delta2_2 as f64)
            * (delta1_1 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (0.72_f64 * (d22_plus as f64) * (d12_plus as f64)
            - 0.72_f64 * (d22_plus as f64) * (d12_min as f64)
            - 1.89_f64 * (d22_min as f64) * (d12_plus as f64)
            + 1.89_f64 * (d22_min as f64) * (d12_min as f64))
            * (delta2_2 as f64)
            * (delta1_2 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (-0.61_f64 * (d22_plus as f64) * (gs1_plus as f64)
            - 1.26_f64 * (d22_plus as f64) * (gs1_min as f64)
            - 0.83_f64 * (d22_min as f64) * (gs1_plus as f64)
            - 2.28_f64 * (d22_min as f64) * (gs1_min as f64))
            * (delta2_2 as f64)
            * (gamma1_44 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)) as f32);
    Gamma1_44 = (((GS_1 as f64)
        + (0.0_f64 * (gs1_plus as f64) * (d21_plus as f64)
            + 0.0_f64 * (gs1_plus as f64) * (d21_min as f64)
            + 0.0_f64 * (gs1_min as f64) * (d21_plus as f64)
            + 0.0_f64 * (gs1_min as f64) * (d21_min as f64))
            * (gamma1_44 as f64)
            * (delta2_1 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)
        + (0.0_f64 * (gs1_plus as f64) * (d22_plus as f64)
            + 0.0_f64 * (gs1_plus as f64) * (d22_min as f64)
            + 0.0_f64 * (gs1_min as f64) * (d22_plus as f64)
            + 0.0_f64 * (gs1_min as f64) * (d22_min as f64))
            * (gamma1_44 as f64)
            * (delta2_2 as f64)
            * (err_kappa as f64).sin().powf(2.0_f64)
        + (0.6_f64 * (gs1_plus as f64) * (gs2_plus as f64)
            + 0.3_f64 * (gs1_plus as f64) * (gs2_min as f64)
            + 0.3_f64 * (gs1_min as f64) * (gs2_plus as f64)
            + 1.1_f64 * (gs1_min as f64) * (gs2_min as f64))
            * (gamma1_44 as f64)
            * (gamma2_44 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)) as f32);
    Gamma2_44 = (((GS_2 as f64)
        + (0.28_f64 * (gs2_plus as f64) * (gs1_plus as f64)
            + 0.0_f64 * (gs2_plus as f64) * (gs1_min as f64)
            + 0.0_f64 * (gs2_min as f64) * (gs1_plus as f64)
            + 2.0_f64 / 3.0_f64 * (gs2_min as f64) * (gs1_min as f64))
            * (gamma2_44 as f64)
            * (gamma1_44 as f64)
            * (err_kappa as f64).cos().powf(2.0_f64)) as f32);
    if index == 1_i32 {
        azim = 0.0_f32;
        (*rsv_1st).global[0_i32 as usize] = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        (*rsv_2nd).global[0_i32 as usize] = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        ) + ISO_2nd;
        azim = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rsv_1st).global[1_i32 as usize] = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        (*rsv_2nd).global[1_i32 as usize] = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        ) + ISO_2nd;
        azim = lkappa;
        (*rsv_1st).global[2_i32 as usize] = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        (*rsv_2nd).global[2_i32 as usize] = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        ) + ISO_2nd;
        azim = kkappa;
        (*rsv_1st).global[3_i32 as usize] = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        (*rsv_2nd).global[3_i32 as usize] = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        ) + ISO_2nd;
    }
    Fst0_0 = (*rsv_1st).global[0_i32 as usize];
    Snd0_0 = (*rsv_2nd).global[0_i32 as usize];
    Fst90_0 = (*rsv_1st).global[1_i32 as usize];
    Snd90_0 = (*rsv_2nd).global[1_i32 as usize];
    FstK_0 = (*rsv_1st).global[2_i32 as usize];
    SndK_0 = (*rsv_2nd).global[2_i32 as usize];
    FstKK_0 = (*rsv_1st).global[3_i32 as usize];
    SndKK_0 = (*rsv_2nd).global[3_i32 as usize];
    if index == 2_i32 {
        azim = 0.0_f32;
        Fst0_90 = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        Snd0_90 = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        ) + ISO_2nd;
        Fst0_kappa = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            true_kappa,
        );
        azim = ((90.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        Fst90_90 = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        Snd90_90 = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        ) + ISO_2nd;
        Fst90_kappa = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            true_kappa,
        );
        azim = lkappa;
        FstK_90 = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        SndK_90 = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        ) + ISO_2nd;
        FstK_kappa = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            true_kappa,
        );
        azim = kkappa;
        FstKK_90 = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        SndKK_90 = SV_term(
            0.0_f32,
            0.0_f32,
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        ) + ISO_2nd;
        FstKK_kappa = SV_term(
            DG_G,
            Dr_r,
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            true_kappa,
        );
        // interpolate the error for the halfspaces rotated by the angle true_kappa
        if ({
            D = Fst0_90 - Fst0_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            Snd0_kappa = ((0.5_f64 * ((Snd0_0 + Snd0_90) as f64)) as f32);
        } else {
            Snd0_kappa = Snd0_0 + (Snd0_90 - Snd0_0) / D * (Fst0_kappa - Fst0_0);
        }
        if ({
            D = Fst90_90 - Fst90_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            Snd90_kappa = ((0.5_f64 * ((Snd90_0 + Snd90_90) as f64)) as f32);
        } else {
            Snd90_kappa = Snd90_0 + (Snd90_90 - Snd90_0) / D * (Fst90_kappa - Fst90_0);
        }
        if ({
            D = FstK_90 - FstK_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            SndK_kappa = ((0.5_f64 * ((SndK_90 + SndK_0) as f64)) as f32);
        } else {
            SndK_kappa = SndK_0 + (SndK_90 - SndK_0) / D * (FstK_kappa - FstK_0);
        }
        if ({
            D = FstKK_90 - FstKK_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            SndKK_kappa = ((0.5_f64 * ((SndKK_90 + SndKK_0) as f64)) as f32);
        } else {
            SndKK_kappa = SndKK_0 + (SndKK_90 - SndKK_0) / D * (FstKK_kappa - FstKK_0);
        }
        (*rsv_1st).global[0_i32 as usize] = Fst0_kappa;
        (*rsv_2nd).global[0_i32 as usize] = Snd0_kappa;
        (*rsv_1st).global[1_i32 as usize] = Fst90_kappa;
        (*rsv_2nd).global[1_i32 as usize] = Snd90_kappa;
        (*rsv_1st).global[2_i32 as usize] = FstK_kappa;
        (*rsv_2nd).global[2_i32 as usize] = SndK_kappa;
        (*rsv_1st).global[3_i32 as usize] = FstKK_kappa;
        (*rsv_2nd).global[3_i32 as usize] = SndKK_kappa;
    }
    // SH
    // The following terms are derived on the assumption that Rps coefficients for
    // two identicaly anisotropic halfspaces are zero (with 0 velocity contrats)
    Delta1_1 = (((SHD1_1 as f64)
        + (0.8_f64 * (d11_plus as f64) * (d21_plus as f64)
            + 0.33_f64 * (d11_plus as f64) * (d21_min as f64)
            + 0.33_f64 * (d11_min as f64) * (d21_plus as f64)
            + 2.0_f64 * (d11_min as f64) * (d21_min as f64))
            * (delta1_1 as f64)
            * (delta2_1 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)
        + (0.8_f64 * (d11_plus as f64) * (d22_plus as f64)
            + 0.33_f64 * (d11_plus as f64) * (d22_min as f64)
            + 0.33_f64 * (d11_min as f64) * (d22_plus as f64)
            + 2.0_f64 * (d11_min as f64) * (d22_min as f64))
            * (delta1_1 as f64)
            * (delta2_2 as f64)
            * ((1.0_f32 * err_kappa) as f64).sin().powf(2.0_f64)
        + (0.22_f64 * (d11_plus as f64) * (gs2_plus as f64)
            + 0.0_f64 * (d11_plus as f64) * (gs2_min as f64)
            + 0.4_f64 * (d11_min as f64) * (gs2_plus as f64)
            + 0.22_f64 * (d11_min as f64) * (gs2_min as f64))
            * (delta1_1 as f64)
            * (gamma2_44 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)) as f32);
    Delta1_2 = (((SHD1_2 as f64)
        + (0.8_f64 * (d12_plus as f64) * (d21_plus as f64)
            + 0.33_f64 * (d12_plus as f64) * (d21_min as f64)
            + 0.33_f64 * (d12_min as f64) * (d21_plus as f64)
            + 2.0_f64 * (d12_min as f64) * (d21_min as f64))
            * (delta1_2 as f64)
            * (delta2_1 as f64)
            * ((1.0_f32 * err_kappa) as f64).sin().powf(2.0_f64)
        + (0.8_f64 * (d12_plus as f64) * (d22_plus as f64)
            + 0.33_f64 * (d12_plus as f64) * (d22_min as f64)
            + 0.33_f64 * (d12_min as f64) * (d22_plus as f64)
            + 2.0_f64 * (d12_min as f64) * (d22_min as f64))
            * (delta1_2 as f64)
            * (delta2_2 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)
        + (-0.22_f64 * (d12_plus as f64) * (gs2_plus as f64)
            + 0.33_f64 * (d12_plus as f64) * (gs2_min as f64)
            - 0.22_f64 * (d12_min as f64) * (gs2_plus as f64)
            - 0.32_f64 * (d12_min as f64) * (gs2_min as f64))
            * (delta1_2 as f64)
            * (gamma2_44 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)) as f32);
    Delta2_1 = (((SHD2_1 as f64)
        + (0.61_f64 * (d21_plus as f64) * (d11_plus as f64)
            + 0.0_f64 * (d21_plus as f64) * (d11_min as f64)
            + 0.0_f64 * (d21_min as f64) * (d11_plus as f64)
            + 1.56_f64 * (d21_min as f64) * (d11_min as f64))
            * (delta2_1 as f64)
            * (delta1_1 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)
        + (0.61_f64 * (d21_plus as f64) * (d12_plus as f64)
            + 0.0_f64 * (d21_plus as f64) * (d12_min as f64)
            + 0.0_f64 * (d21_min as f64) * (d12_plus as f64)
            + 1.56_f64 * (d21_min as f64) * (d12_min as f64))
            * (delta2_1 as f64)
            * (delta1_2 as f64)
            * ((1.0_f32 * err_kappa) as f64).sin().powf(2.0_f64)
        + (-0.28_f64 * (d21_plus as f64) * (gs1_plus as f64)
            - 1.1_f64 * (d21_plus as f64) * (gs1_min as f64)
            - 0.44_f64 * (d21_min as f64) * (gs1_plus as f64)
            - 1.56_f64 * (d21_min as f64) * (gs1_min as f64))
            * (delta2_1 as f64)
            * (gamma1_44 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)) as f32);
    Delta2_2 = (((SHD2_2 as f64)
        + (0.61_f64 * (d22_plus as f64) * (d11_plus as f64)
            + 0.0_f64 * (d22_plus as f64) * (d11_min as f64)
            + 0.0_f64 * (d22_min as f64) * (d11_plus as f64)
            + 1.56_f64 * (d22_min as f64) * (d11_min as f64))
            * (delta2_2 as f64)
            * (delta1_1 as f64)
            * ((1.0_f32 * err_kappa) as f64).sin().powf(2.0_f64)
        + (0.61_f64 * (d22_plus as f64) * (d12_plus as f64)
            + ((0.0_f32 * d22_plus * d12_min) as f64)
            + ((0.0_f32 * d22_min * d12_plus) as f64)
            + 1.56_f64 * (d22_min as f64) * (d12_min as f64))
            * (delta2_2 as f64)
            * (delta1_2 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)
        + (-0.11_f64 * (d22_plus as f64) * (gs1_plus as f64)
            - 0.11_f64 * (d22_plus as f64) * (gs1_min as f64)
            - 0.0_f64 * (d22_min as f64) * (gs1_plus as f64)
            - 0.64_f64 * (d22_min as f64) * (gs1_min as f64))
            * (delta2_2 as f64)
            * (gamma1_44 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)) as f32);
    Gamma1_44 = (((SHGS_1 as f64)
        + (-0.0_f64 * (gs1_plus as f64) * (d21_plus as f64)
            - 0.0_f64 * (gs1_plus as f64) * (d21_min as f64)
            - 0.0_f64 * (gs1_min as f64) * (d21_plus as f64)
            - 0.0_f64 * (gs1_min as f64) * (d21_min as f64))
            * (gamma1_44 as f64)
            * (delta2_1 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)
        + (0.0_f64 * (gs1_plus as f64) * (d22_plus as f64)
            + 0.0_f64 * (gs1_plus as f64) * (d22_min as f64)
            + 0.0_f64 * (gs1_min as f64) * (d22_plus as f64)
            + 0.0_f64 * (gs1_min as f64) * (d22_min as f64))
            * (gamma1_44 as f64)
            * (delta2_2 as f64)
            * ((2.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)
        + (0.45_f64 * (gs1_plus as f64) * (gs2_plus as f64)
            + 0.0_f64 * (gs1_plus as f64) * (gs2_min as f64)
            + 0.1_f64 * (gs1_min as f64) * (gs2_plus as f64)
            + 0.78_f64 * (gs1_min as f64) * (gs2_min as f64))
            * (gamma1_44 as f64)
            * (gamma2_44 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)) as f32);
    Gamma2_44 = (((SHGS_2 as f64)
        + (1.0_f64 / 3.0_f64 * (gs2_plus as f64) * (gs1_plus as f64)
            + ((0.0_f32 * gs2_plus * gs1_min) as f64)
            + ((0.0_f32 * gs2_min * gs1_plus) as f64)
            + 2.0_f64 / 3.0_f64 * (gs2_min as f64) * (gs1_min as f64))
            * (gamma2_44 as f64)
            * (gamma1_44 as f64)
            * ((1.0_f32 * err_kappa) as f64).cos().powf(2.0_f64)) as f32);
    if index == 1_i32 {
        azim = ((20.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rsh_1st).global[0_i32 as usize] = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        (*rsh_2nd).global[0_i32 as usize] = SH_term(
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        azim = ((45.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rsh_1st).global[1_i32 as usize] = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        (*rsh_2nd).global[1_i32 as usize] = SH_term(
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        azim = ((70.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        (*rsh_1st).global[2_i32 as usize] = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        (*rsh_2nd).global[2_i32 as usize] = SH_term(
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        (*rsh_1st).global[3_i32 as usize] = 0.0_f32;
        (*rsh_2nd).global[3_i32 as usize] = 0.0_f32;
    }
    SHFst20_0 = (*rsh_1st).global[0_i32 as usize];
    SHSnd20_0 = (*rsh_2nd).global[0_i32 as usize];
    SHFst45_0 = (*rsh_1st).global[1_i32 as usize];
    SHSnd45_0 = (*rsh_2nd).global[1_i32 as usize];
    SHFst70_0 = (*rsh_1st).global[2_i32 as usize];
    SHSnd70_0 = (*rsh_2nd).global[2_i32 as usize];
    if index == 2_i32 {
        azim = ((20.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        SHFst20_90 = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        SHSnd20_90 = SH_term(
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        SHFst20_kappa = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            true_kappa,
        );
        azim = ((45.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        SHFst45_90 = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        SHSnd45_90 = SH_term(
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        SHFst45_kappa = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            true_kappa,
        );
        azim = ((70.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32);
        SHFst70_90 = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        SHSnd70_90 = SH_term(
            Delta1_1,
            Delta1_2,
            0.0_f32,
            Delta2_1,
            Delta2_2,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            0.0_f32,
            Gamma1_44,
            Gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            err_kappa,
        );
        SHFst70_kappa = SH_term(
            delta1_1,
            delta1_2,
            delta1_3,
            delta2_1,
            delta2_2,
            delta2_3,
            eps1_1,
            eps1_2,
            eps2_1,
            eps2_2,
            gamma1_44,
            gamma2_44,
            alpha,
            beta,
            ((30.0_f64 * 3.141592653589793_f64 / 180.0_f64) as f32),
            azim,
            true_kappa,
        );
        // interpolate the error for the halfspaces rotated by the angle true_kappa
        if ({
            D = SHFst20_90 - SHFst20_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            SHSnd20_kappa = ((0.5_f64 * ((SHSnd20_0 + SHSnd20_90) as f64)) as f32);
        } else {
            SHSnd20_kappa = SHSnd20_0 + (SHSnd20_90 - SHSnd20_0) / D * (SHFst20_kappa - SHFst20_0);
        }
        if ({
            D = SHFst45_90 - SHFst45_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            SHSnd45_kappa = ((0.5_f64 * ((SHSnd45_0 + SHSnd45_90) as f64)) as f32);
        } else {
            SHSnd45_kappa = SHSnd45_0 + (SHSnd45_90 - SHSnd45_0) / D * (SHFst45_kappa - SHFst45_0);
        }
        if ({
            D = SHFst70_90 - SHFst70_0;
            D
        } as f64)
            .abs()
            < 1.1920929e-07_f64
        {
            SHSnd70_kappa = ((0.5_f64 * ((SHSnd70_0 + SHSnd70_90) as f64)) as f32);
        } else {
            SHSnd70_kappa = SHSnd70_0 + (SHSnd70_90 - SHSnd70_0) / D * (SHFst70_kappa - SHFst70_0);
        }
        (*rsh_1st).global[0_i32 as usize] = SHFst20_kappa;
        (*rsh_2nd).global[0_i32 as usize] = SHSnd20_kappa;
        (*rsh_1st).global[1_i32 as usize] = SHFst45_kappa;
        (*rsh_2nd).global[1_i32 as usize] = SHSnd45_kappa;
        (*rsh_1st).global[2_i32 as usize] = SHFst70_kappa;
        (*rsh_2nd).global[2_i32 as usize] = SHSnd70_kappa;
    }
    return 0_i32;
}

/// First-order SV-wave angular term.
///
/// C prototype: `float SV_term(float DG_G, float Dr_r, float Delta1_1, float Delta1_2, float Delta1_3, float Delta2_1, float Delta2_2, float Delta2_3, float Eps1_1, float Eps1_2, float Eps2_1, float Eps2_2, float Gamma1_44, float Gamma2_44, float alpha, float beta, float ang, float azim, float kappa);`
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
) -> f32 {
    let mut cos2a: f32;
    let mut sin2a: f32;
    let mut cos2a_k: f32;
    let mut sin2a_k: f32;
    let mut CS: f32;
    let mut S_C: f32;
    let mut CS3: f32;
    let mut S3_C: f32;
    let mut S5_C: f32;
    let mut result: f32 = 0.0_f32;
    cos2a_k = (((azim - kappa) as f64).cos().powf(2.0_f64) as f32);
    sin2a_k = (((azim - kappa) as f64).sin().powf(2.0_f64) as f32);
    cos2a = ((azim as f64).cos().powf(2.0_f64) as f32);
    sin2a = ((azim as f64).sin().powf(2.0_f64) as f32);
    CS = ((((-beta / alpha * DG_G) as f64)
        - ((alpha * beta) as f64)
            / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
            * ((Delta2_2 * cos2a_k + Delta2_1 * sin2a_k - Delta1_1 * sin2a - Delta1_2 * cos2a)
                as f64)
        - ((2.0_f32 * beta / alpha * (Gamma2_44 * sin2a_k - Gamma1_44 * sin2a)) as f64))
        as f32);
    S_C = ((-0.5_f64 * (Dr_r as f64)
        + (alpha as f64).powf(2.0_f64)
            / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
            * ((Delta2_2 * cos2a_k + Delta2_1 * sin2a_k - Delta1_2 * cos2a - Delta1_1 * sin2a)
                as f64)) as f32);
    CS3 = ((((-alpha * beta) as f64) / ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64))
        * ((Delta1_1 * sin2a + Delta1_2 * cos2a
            - Delta2_1 * sin2a_k
            - Delta2_2 * cos2a_k
            - Delta1_3 * cos2a * sin2a
            - Eps1_2 * (cos2a * cos2a + 2.0_f32 * cos2a * sin2a)
            - Eps1_1 * sin2a * sin2a
            + Eps2_2 * (cos2a_k * cos2a_k + 2.0_f32 * sin2a_k * cos2a_k)
            + Eps2_1 * sin2a_k * sin2a_k
            + Delta2_3 * cos2a_k * sin2a_k) as f64)) as f32);
    S3_C = (((beta as f64).powf(2.0_f64) / (alpha as f64).powf(2.0_f64) * (DG_G as f64)
        + (beta as f64).powf(2.0_f64)
            / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
            * ((Delta1_1 * sin2a + Delta1_2 * cos2a - Delta2_1 * sin2a_k - Delta2_2 * cos2a_k)
                as f64)
        + (alpha as f64).powf(2.0_f64)
            / ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64))
            * ((Delta1_1 * sin2a + Delta1_2 * cos2a
                - Delta2_1 * sin2a_k
                - Delta2_2 * cos2a_k
                - Delta1_3 * cos2a * sin2a
                - Eps1_1 * sin2a * sin2a
                - Eps1_2 * (cos2a * cos2a + 2.0_f32 * cos2a * sin2a)
                + Eps2_1 * sin2a_k * sin2a_k
                + Eps2_2 * (cos2a_k * cos2a_k + 2.0_f32 * sin2a_k * cos2a_k)
                + Delta2_3 * cos2a_k * sin2a_k) as f64)
        + 2.0_f64
            * ((beta / alpha) as f64).powf(2.0_f64)
            * ((Gamma2_44 * sin2a_k - Gamma1_44 * sin2a) as f64)) as f32);
    S5_C = ((-(beta as f64).powf(2.0_f64)
        / ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64))
        * ((Delta1_1 * sin2a + Delta1_2 * cos2a
            - Delta2_1 * sin2a_k
            - Delta2_2 * cos2a_k
            - Delta1_3 * cos2a * sin2a
            - Eps1_1 * sin2a * sin2a
            - Eps1_2 * (cos2a * cos2a + 2.0_f32 * cos2a * sin2a)
            + Eps2_1 * sin2a_k * sin2a_k
            + Eps2_2 * (cos2a_k * cos2a_k + 2.0_f32 * sin2a_k * cos2a_k)
            + Delta2_3 * cos2a_k * sin2a_k) as f64)) as f32);
    result = (((CS as f64) * (ang as f64).cos() * (ang as f64).sin()
        + (S_C as f64) * (ang as f64).sin()
        + (CS3 as f64) * (ang as f64).cos() * (ang as f64).sin().powf(3.0_f64)
        + ((S3_C as f64) + (S_C as f64) * 0.5_f64 * ((beta / alpha) as f64).powf(2.0_f64))
            * (ang as f64).sin().powf(3.0_f64)
        + ((S5_C as f64) + (S3_C as f64) * 0.5_f64 * ((beta / alpha) as f64).powf(2.0_f64))
            * (ang as f64).sin().powf(5.0_f64)
        + (S5_C as f64)
            * 0.5_f64
            * ((beta / alpha) as f64).powf(2.0_f64)
            * (ang as f64).sin().powf(7.0_f64)) as f32);
    return result;
}

/// First-order SH-wave angular term.
///
/// C prototype: `float SH_term(float Delta1_1, float Delta1_2, float Delta1_3, float Delta2_1, float Delta2_2, float Delta2_3, float Eps1_1, float Eps1_2, float Eps2_1, float Eps2_2, float Gamma1_44, float Gamma2_44, float alpha, float beta, float ang, float azim, float kappa);`
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
) -> f32 {
    let mut cos2a: f32;
    let mut sin2a: f32;
    let mut cos2a_k: f32;
    let mut sin2a_k: f32;
    let mut S: f32;
    let mut CS_C: f32;
    let mut CS3_C: f32;
    let mut S3: f32;
    let mut result: f32;
    cos2a_k = (((azim - kappa) as f64).cos().powf(2.0_f64) as f32);
    sin2a_k = (((azim - kappa) as f64).sin().powf(2.0_f64) as f32);
    cos2a = ((azim as f64).cos().powf(2.0_f64) as f32);
    sin2a = ((azim as f64).sin().powf(2.0_f64) as f32);
    S = (((alpha as f64).powf(2.0_f64)
        / (4.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
        * (((Delta2_2 - Delta2_1) as f64) * ((2.0_f32 * (azim - kappa)) as f64).sin()
            + ((Delta1_1 - Delta1_2) as f64) * ((2.0_f32 * azim) as f64).sin())) as f32);
    CS_C = ((((alpha * beta) as f64)
        / (4.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
        * (((Delta2_1 - Delta2_2) as f64) * ((2.0_f32 * (azim - kappa)) as f64).sin()
            + ((Delta1_2 - Delta1_1) as f64) * ((2.0_f32 * azim) as f64).sin())
        + ((beta / alpha) as f64)
            * ((Gamma2_44 as f64) * ((2.0_f32 * (azim - kappa)) as f64).sin()
                - (Gamma1_44 as f64) * ((2.0_f32 * azim) as f64).sin())) as f32);
    CS3_C = ((((alpha * beta) as f64)
        / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
        * (((Delta1_1 - Delta1_2 - Delta1_3 * (cos2a - sin2a) + 2.0_f32 * (Eps1_2 - Eps1_1) * sin2a)
            as f64)
            * (azim as f64).sin()
            * (azim as f64).cos()
            - ((Delta2_1 - Delta2_2 - Delta2_3 * (cos2a_k - sin2a_k)
                + 2.0_f32 * (Eps2_2 - Eps2_1) * sin2a_k) as f64)
                * ((azim - kappa) as f64).sin()
                * ((azim - kappa) as f64).cos())) as f32);
    S3 = (((alpha as f64).powf(2.0_f64)
        / (2.0_f64 * ((alpha as f64).powf(2.0_f64) - (beta as f64).powf(2.0_f64)))
        * (((-Delta1_1 + Delta1_2 + Delta1_3 * (cos2a - sin2a)
            - 2.0_f32 * (Eps1_2 - Eps1_1) * sin2a) as f64)
            * (azim as f64).sin()
            * (azim as f64).cos()
            + ((Delta2_1 - Delta2_2 - Delta2_3 * (cos2a_k - sin2a_k)
                + 2.0_f32 * (Eps2_2 - Eps2_1) * sin2a_k) as f64)
                * ((azim - kappa) as f64).sin()
                * ((azim - kappa) as f64).cos())) as f32);
    result = (((S as f64) * (ang as f64).sin()
        + (CS_C as f64) * (ang as f64).cos() * (ang as f64).sin()
        + ((CS3_C as f64) + (CS_C as f64) * 0.5_f64 * ((beta / alpha) as f64).powf(2.0_f64))
            * (ang as f64).cos()
            * (ang as f64).sin().powf(3.0_f64)
        + (S3 as f64) * (ang as f64).sin().powf(3.0_f64)
        + (CS3_C as f64)
            * 0.5_f64
            * ((beta / alpha) as f64).powf(2.0_f64)
            * (ang as f64).cos()
            * (ang as f64).sin().powf(5.0_f64)) as f32);
    return result;
}

/// Rotate S-wave polarizations into the selected projection plane.
/// `iso_plane == 1` selects the plane normal to the polarization cross product;
/// other values select the background-isotropic SV/SH plane. Returns -1 for
/// degenerate polarizations, preserving C diagnostics and partial mutations.
///
/// C prototype: `int Phi_rot(float * rs1, float * rs2, int iso_plane, float pb_x, float pb_y, float pb_z, float gs1_x, float gs1_y, float gs1_z, float gs2_x, float gs2_y, float gs2_z, float * CPhi1, float * SPhi1, float * CPhi2, float * SPhi2);`
pub fn Phi_rot(
    rs1: &mut f32,
    rs2: &mut f32,
    iso_plane: i32,
    pb_x: f32,
    pb_y: f32,
    pb_z: f32,
    mut gs1_x: f32,
    mut gs1_y: f32,
    mut gs1_z: f32,
    mut gs2_x: f32,
    mut gs2_y: f32,
    mut gs2_z: f32,
    CPhi1: &mut f32,
    SPhi1: &mut f32,
    CPhi2: &mut f32,
    SPhi2: &mut f32,
) -> i32 {
    let mut gv_x: f32 = 0.0_f32;
    let mut gv_y: f32 = 0.0_f32;
    let mut gv_z: f32 = 0.0_f32;
    let mut gh_x: f32 = 0.0_f32;
    let mut gh_y: f32 = 0.0_f32;
    let mut gh_z: f32 = 0.0_f32;
    let mut p_h: f32;
    let mut norm: f32;
    let mut denom: f32;
    let mut alpha: f32;
    let mut temp: f32;
    let mut ax_x: f32;
    let mut ax_y: f32;
    let mut ax_z: f32;
    let mut n_x: f32;
    let mut n_y: f32;
    let mut n_z: f32;
    let mut prl_x: f32;
    let mut prl_y: f32;
    let mut prl_z: f32;
    let mut s1_prl: f32;
    let mut s1_ax: f32;
    let mut s2_prl: f32;
    let mut s2_ax: f32;
    let mut s1_v: f32;
    let mut s1_h: f32;
    let mut s2_v: f32;
    let mut s2_h: f32;
    let mut a11: f32;
    let mut a12: f32;
    let mut g10_v: f32;
    let mut g10_h: f32;
    let mut g20_v: f32;
    let mut g20_h: f32;
    let mut az_n: f32;
    let mut an_n: f32;
    let mut sv_n_x: f32;
    let mut sv_n_y: f32;
    let mut sv_n_z: f32;
    let mut sh_n_x: f32;
    let mut sh_n_y: f32;
    let mut sh_n_z: f32;
    // let's start: just to make sure
    norm = (((gs1_x as f64).powf(2.0_f64)
        + (gs1_y as f64).powf(2.0_f64)
        + (gs1_z as f64).powf(2.0_f64))
    .sqrt() as f32);
    gs1_x = gs1_x / norm;
    gs1_y = gs1_y / norm;
    gs1_z = gs1_z / norm;
    norm = (((gs2_x as f64).powf(2.0_f64)
        + (gs2_y as f64).powf(2.0_f64)
        + (gs2_z as f64).powf(2.0_f64))
    .sqrt() as f32);
    gs2_x = gs2_x / norm;
    gs2_y = gs2_y / norm;
    gs2_z = gs2_z / norm;
    if iso_plane != 1_i32 {
        // now, determine gsv and gsh polarization vectors in the background isotropy plane
        p_h = (((pb_x as f64).powf(2.0_f64) + (pb_y as f64).powf(2.0_f64)).sqrt() as f32);
        if (p_h as f64) != 0.0_f64 {
            gh_x = -pb_y / p_h;
            if (pb_x as f64) != 0.0_f64 {
                gh_y = (((pb_x as f64) / (pb_x as f64).abs()
                    * (1.0_f64 - (gh_x as f64).powf(2.0_f64)).sqrt())
                    as f32);
            }
            if (pb_x as f64) == 0.0_f64 {
                gh_y = 0.0_f32;
            }
        } else {
            gh_x = 0.0_f32;
            gh_y = 1.0_f32;
        }
        gh_z = 0.0_f32;
        gv_x = gh_y * pb_z - gh_z * pb_y;
        gv_y = gh_z * pb_x - gh_x * pb_z;
        gv_z = gh_x * pb_y - gh_y * pb_x;
        norm = (((gv_x as f64).powf(2.0_f64)
            + (gv_y as f64).powf(2.0_f64)
            + (gv_z as f64).powf(2.0_f64))
        .sqrt() as f32);
        gv_x = gv_x / norm;
        gv_y = gv_y / norm;
        gv_z = gv_z / norm;
    }
    // normal vector to the gs1-gs2 plane
    n_x = gs1_y * gs2_z - gs2_y * gs1_z;
    n_y = gs1_z * gs2_x - gs2_z * gs1_x;
    n_z = gs1_x * gs2_y - gs2_x * gs1_y;
    norm = (((n_x as f64).powf(2.0_f64) + (n_y as f64).powf(2.0_f64) + (n_z as f64).powf(2.0_f64))
        .sqrt() as f32);
    if norm == 0.0_f32 {
        eprint!("! Paralel S-wave polarizations !");
        eprint!(
            "  Either wrong data or extreme anisotropy or a vicinity of S-wave singularity! \n"
        );
        return -1_i32;
    }
    n_x = n_x / norm;
    n_y = n_y / norm;
    n_z = n_z / norm;
    if ((pb_x * n_x + pb_y * n_y + pb_z * n_z) as f64) < 0.0_f64 {
        n_x = -n_x;
        n_y = -n_y;
        n_z = -n_z;
        temp = gs1_x;
        gs1_x = gs2_x;
        gs2_x = temp;
        temp = gs1_y;
        gs1_y = gs2_y;
        gs2_y = temp;
        temp = gs1_z;
        gs1_z = gs2_z;
        gs2_z = temp;
        temp = *rs1;
        *rs1 = *rs2;
        *rs2 = temp;
    } else {
        if ((pb_x * n_x + pb_y * n_y + pb_z * n_z) as f64) == 0.0_f64 {
            eprint!(" ! S-wave polarizations in a plane paralel to the slowness vector!");
            eprint!("   Either wrong data or extreme anisotropy \n");
            return -1_i32;
        }
    }
    'projection: {
        // the last choice if the plane of projection: iso or aniso ?
        if iso_plane == 1_i32 {
            az_n = ((n_y as f64).atan2((n_x as f64)) as f32);
            an_n = ((n_z as f64)
                .atan2(((n_x as f64).powf(2.0_f64) + (n_y as f64).powf(2.0_f64)).sqrt())
                as f32);
            sv_n_x = (((an_n as f64).sin() * (az_n as f64).cos()) as f32);
            sv_n_y = (((an_n as f64).sin() * (az_n as f64).sin()) as f32);
            sv_n_z = (-(an_n as f64).cos() as f32);
            sh_n_x = (-(az_n as f64).sin() as f32);
            sh_n_y = ((az_n as f64).cos() as f32);
            sh_n_z = 0.0_f32;
            s1_v = gs1_x * sv_n_x + gs1_y * sv_n_y + gs1_z * sv_n_z;
            s1_h = gs1_x * sh_n_x + gs1_y * sh_n_y + gs1_z * sh_n_z;
            s2_v = gs2_x * sv_n_x + gs2_y * sv_n_y + gs2_z * sv_n_z;
            s2_h = gs2_x * sh_n_x + gs2_y * sh_n_y + gs2_z * sh_n_z;
            break 'projection;
        }
        // axis of rotation of the planes gs1-gs2 -> gs10-gs20
        ax_x = n_y * pb_z - n_z * pb_y;
        ax_y = n_z * pb_x - n_x * pb_z;
        ax_z = n_x * pb_y - n_y * pb_x;
        norm = (((ax_x as f64).powf(2.0_f64)
            + (ax_y as f64).powf(2.0_f64)
            + (ax_z as f64).powf(2.0_f64))
        .sqrt() as f32);
        if norm == 0.0_f32 {
            ax_x = gv_x;
            ax_y = gv_y;
            ax_z = gv_z;
        } else {
            ax_x = ax_x / norm;
            ax_y = ax_y / norm;
            ax_z = ax_z / norm;
        }
        // paralel vector with the plane gs1-gs2, perpendicular to ax and n:
        // (n,ax,prl) complete right-hand orthogonal basis
        prl_x = n_y * ax_z - n_z * ax_y;
        prl_y = n_z * ax_x - n_x * ax_z;
        prl_z = n_x * ax_y - n_y * ax_x;
        norm = (((prl_x as f64).powf(2.0_f64)
            + (prl_y as f64).powf(2.0_f64)
            + (prl_z as f64).powf(2.0_f64))
        .sqrt() as f32);
        prl_x = prl_x / norm;
        prl_y = prl_y / norm;
        prl_z = prl_z / norm;
        // coordinates of gs1, gs2 in (ax,prl) plane
        s1_prl = gs1_x * prl_x + gs1_y * prl_y + gs1_z * prl_z;
        s1_ax = gs1_x * ax_x + gs1_y * ax_y + gs1_z * ax_z;
        s2_prl = gs2_x * prl_x + gs2_y * prl_y + gs2_z * prl_z;
        s2_ax = gs2_x * ax_x + gs2_y * ax_y + gs2_z * ax_z;
        // rotation of gs1, gs2 onto sv-sh plane
        a11 = ax_x * gv_x + ax_y * gv_y + ax_z * gv_z;
        a12 = ax_x * gh_x + ax_y * gh_y + ax_z * gh_z;
        s1_v = a11 * s1_ax - a12 * s1_prl;
        s1_h = a12 * s1_ax + a11 * s1_prl;
        s2_v = a11 * s2_ax - a12 * s2_prl;
        s2_h = a12 * s2_ax + a11 * s2_prl;
    }
    // the closest stretched perpendicular polarization vectors g10, g20
    if (!({
        denom = s1_v * s2_h - s1_h * s2_v;
        denom
    } != 0.0))
    {
        eprint!("this should not happen \n");
        return -1_i32;
    }
    alpha = ((0.5_f64
        * (3.141592653589793_f64 / 2.0_f64 - ((s1_v * s2_v + s1_h * s2_h) as f64).acos()))
        as f32);
    // this is for stretching of vectors up to a common angle Phi
    if 0_i32 == 1_i32 {
        g10_v = ((((alpha as f64).cos() * (s2_h as f64) - (alpha as f64).sin() * (s1_h as f64))
            / (denom as f64)) as f32);
        g10_h = ((-((alpha as f64).cos() * (s2_v as f64) - (alpha as f64).sin() * (s1_v as f64))
            / (denom as f64)) as f32);
        g20_v = ((-((alpha as f64).cos() * (s1_h as f64) - (alpha as f64).sin() * (s2_h as f64))
            / (denom as f64)) as f32);
        g20_h = ((((alpha as f64).cos() * (s1_v as f64) - (alpha as f64).sin() * (s2_v as f64))
            / (denom as f64)) as f32);
        s1_v = g10_v;
        s1_h = g10_h;
        s2_v = g20_v;
        s2_h = g20_h;
    }
    // put s1 wave into 1 or 4 quadrant
    if (s1_v as f64) < 0.0_f64 {
        s1_v = -s1_v;
        s1_h = -s1_h;
        *rs1 = -*rs1;
    }
    // put s2 wave into 1 or 2 quadrant
    if (s2_h as f64) < 0.0_f64 {
        s2_v = -s2_v;
        s2_h = -s2_h;
        *rs2 = -*rs2;
    }
    // make sure s1->s2 is right-handed
    if ((s1_v * s2_h - s1_h * s2_v) as f64) < 0.0_f64 {
        temp = s1_v;
        s1_v = s2_v;
        s2_v = temp;
        temp = s1_h;
        s1_h = s2_h;
        s2_h = temp;
        temp = *rs1;
        *rs1 = *rs2;
        *rs2 = temp;
    }
    // cos and sin of Phi
    *CPhi1 = s1_v;
    *SPhi1 = s1_h;
    *CPhi2 = s2_v;
    *SPhi2 = s2_h;
    // O.K. this is it
    return 0_i32;
}
