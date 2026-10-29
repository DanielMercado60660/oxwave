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

//! Reflection coefficients for anisotropic media. Port of par/lib/refaniso.c.
//! Author: CWP: Andreas Rueger, 1994-1996, Colorado School of Mines
//!
//! C `double` inputs and outputs remain `f64`. Tensor pointers become borrowed
//! [`Stiff2D`] values and scalar output pointers become mutable references.
//! Return codes, partial output writes, and the C `FLT_EPSILON` thresholds are
//! preserved, including for singular or unphysical inputs.
#![allow(non_snake_case, clippy::too_many_arguments)]

use oxwave_core::{cio::stderr, fprintf};

/// Nine-component stiffness tensor with the C layout from par/include/anisotropy.h.
/// Axes 1 and 3 are in the propagation plane; axis 2 is out of plane.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct Stiff2D {
    pub a1111: f64,
    pub a3333: f64,
    pub a1133: f64,
    pub a1313: f64,
    pub a1113: f64,
    pub a3313: f64,
    pub a1212: f64,
    pub a2323: f64,
    pub a1223: f64,
}

/// Compute SH group velocity; the three in-plane polarization outputs are zero.
///
/// C: `void gvelpolSH(double a1212, double a2323, double a1223, double px, double pz, double *vgx, double *vgz, double *g11, double *g13, double *g33)`.
pub fn gvelpolSH(
    a: f64,
    b: f64,
    c: f64,
    px: f64,
    pz: f64,
    vgx: &mut f64,
    vgz: &mut f64,
    g11: &mut f64,
    g13: &mut f64,
    g33: &mut f64,
) {
    *vgx = a * px + c * pz;
    *vgz = c * px + b * pz;
    *g11 = 0.;
    *g13 = 0.;
    *g33 = 0.;
}
/// Compute P/SV group velocity and polarization from the supplied slowness.
/// Returns -1 for invalid polarization, after writing the group velocities but
/// before writing the polarization outputs.
///
/// C: `int gvelpolTI(double a1111, double a3333, double a1133, double a1313, double a1113, double a3313, double px, double pz, double *vgx, double *vgz, double *g11n, double *g13n, double *g33n)`.
pub fn gvelpolTI(
    a: f64,
    c: f64,
    f: f64,
    l: f64,
    m: f64,
    n: f64,
    px: f64,
    pz: f64,
    vgx: &mut f64,
    vgz: &mut f64,
    g11n: &mut f64,
    g13n: &mut f64,
    g33n: &mut f64,
) -> i32 {
    let px2 = px * px;
    let pz2 = pz * pz;
    let pxz = px * pz;
    let t11 = a * px2 + l * pz2 + 2. * m * pxz;
    let t33 = c * pz2 + l * px2 + 2. * n * pxz;
    let t13 = (f + l) * pxz + m * px2 + n * pz2;
    let den = 1. / (t11 + t33 - 2.);
    let mut g11 = (t33 - 1.) * den;
    let mut g33 = (t11 - 1.) * den;
    let g13 = -t13 * den;
    *vgx = a * px * g11
        + (f + l) * pz * g13
        + n * pz * g33
        + m * (pz * g11 + 2. * px * g13)
        + l * g33 * px;
    *vgz = c * pz * g33
        + (f + l) * px * g13
        + m * px * g11
        + n * (px * g33 + 2. * pz * g13)
        + l * g11 * pz;
    if g11 < -(f32::EPSILON as f64) || g33 < -(f32::EPSILON as f64) {
        return -1;
    } else if g11 < 0. {
        g11 = 0.
    } else if g33 < 0. {
        g33 = 0.
    }
    *g11n = g11;
    *g13n = g13;
    *g33n = g33;
    1
}
/// Compute horizontal slowness for sine of incidence angle `s`.
/// Modes 0 and 1 select qP and qSV. Invalid modes or angles return -1 without
/// writing `p`; negative squared velocity calls the C-compatible fatal error.
///
/// C: `int p_hor2DTI(Stiff2D *spar, double s, int mode, double *p)`.
pub fn p_hor2DTI(a: &Stiff2D, s: f64, mode: i32, p: &mut f64) -> i32 {
    if mode != 0 && mode != 1 {
        return -1;
    }
    let mut cc = 1. - s * s;
    if cc.abs() < (f32::EPSILON as f64) {
        cc = 0.
    }
    if cc < 0. {
        return -1;
    }
    let c = cc.sqrt();
    let ss = s * s;
    let sc = s * c;
    let g11 = a.a1111 * ss + a.a1313 * cc + 2. * a.a1113 * sc;
    let g33 = a.a3333 * cc + a.a1313 * ss + 2. * a.a3313 * sc;
    let g13 = (a.a1133 + a.a1313) * sc + a.a1113 * ss + a.a3313 * cc;
    let sqr = ((g11 + g33) * (g11 + g33) - 4. * (g11 * g33 - g13 * g13)).sqrt();
    let vp2 = g11 + g33;
    let vp2 = if mode == 0 { vp2 + sqr } else { vp2 - sqr };
    if vp2 < 0. {
        oxwave_core::err!("\n ERROR: ray initialization \n");
    }
    *p = s / (vp2 * 0.5).sqrt();
    1
}
/// Compute vertical slowness from horizontal slowness in a VTI medium.
/// Mode 0 selects P; every other mode selects SV, as in C. Detected evanescent
/// energy writes a diagnostic and returns -1 without changing the output.
///
/// C: `int p_vert2DVTI(Stiff2D *spar1, double pl, int modei, double *p_vert)`.
pub fn p_vert2DVTI(a: &Stiff2D, pl: f64, mode: i32, out: &mut f64) -> i32 {
    let p2 = pl * pl;
    let o3 = 1. / a.a3333;
    let o13 = 1. / a.a1313;
    let sqr = a.a1111 * o13 + a.a1313 * o3 - (a.a1313 + a.a1133) * (a.a1313 + a.a1133) * o3 * o13;
    let k1 = o3 + o13 - sqr * p2;
    let k2 = a.a1111 * o3 * p2 - o3;
    let k3 = p2 - o13;
    let mut sqr = k1 * k1 - 4. * k2 * k3;
    if sqr.abs() < (f32::EPSILON as f64) {
        sqr = 0.
    }
    if sqr < 0. {
        fprintf!(stderr(), " evanescent energy generated  \n");
        return -1;
    }
    let qasq = 0.5 * (k1 - sqr.sqrt());
    if qasq < 0. {
        fprintf!(stderr(), "evanescent energy generated  \n");
        return -1;
    }
    let qa = qasq.sqrt();
    let qb = (0.5 * (k1 + sqr.sqrt())).sqrt();
    *out = if mode == 0 { qa } else { qb };
    1
}
/// Rotate the stiffness tensor counterclockwise by `phi` radians.
///
/// C: `void rottens2D(Stiff2D *spar, double phi)`.
pub fn rottens2D(s: &mut Stiff2D, phi: f64) {
    let si = phi.sin();
    let co = phi.cos();
    let ss = si * si;
    let cc = co * co;
    let sc = si * co;
    let (a, c, f, l, m, n, o, p, q) = (
        s.a1111, s.a3333, s.a1133, s.a1313, s.a1113, s.a3313, s.a1212, s.a2323, s.a1223,
    );
    s.a1111 = a * cc * cc + c * ss * ss + 2. * f * cc * ss + 4. * l * cc * ss
        - 4. * m * cc * sc
        - 4. * n * ss * sc;
    s.a3333 = a * ss * ss
        + c * cc * cc
        + 2. * f * cc * ss
        + 4. * l * cc * ss
        + 4. * m * ss * sc
        + 4. * n * cc * sc;
    s.a1133 = a * ss * cc + c * ss * cc + f * (ss * ss + cc * cc)
        - 4. * l * ss * cc
        - 2. * m * (ss * sc - cc * sc)
        - 2. * n * (cc * sc - ss * sc);
    s.a1313 = a * ss * cc + c * ss * cc - 2. * f * ss * cc + l * (cc * cc + ss * ss - 2. * ss * cc)
        - 2. * m * (ss * sc - cc * sc)
        - 2. * n * (cc * sc - ss * sc);
    s.a1113 = a * cc * sc - c * ss * sc
        + f * (ss * sc - cc * sc)
        + 2. * l * (ss * sc - cc * sc)
        + m * (cc * cc - 3. * cc * ss)
        + (3. * cc * ss - ss * ss) * n;
    s.a3313 = a * ss * sc - c * cc * sc
        + f * (cc * sc - ss * sc)
        + 2. * l * (cc * sc - ss * sc)
        + m * (3. * cc * ss - ss * ss)
        + n * (cc * cc - 3. * ss * cc);
    s.a1212 = cc * o - 2. * sc * q + ss * p;
    s.a2323 = ss * o + 2. * sc * q + cc * p;
    s.a1223 = sc * o + (cc - ss) * q - sc * p;
}
/// Convert Thomsen parameters to density-normalized, tilted TI stiffnesses.
/// Returns 1 on success or -1 for rejected parameters. A rejection can leave
/// the first three computed tensor components updated, matching C.
///
/// C: `int thom2stiffTI(double vp, double vs, double eps, double delta, double gamma, double phi, Stiff2D *spar, int sign)`.
pub fn thom2stiffTI(
    vp: f64,
    vs: f64,
    eps: f64,
    delta: f64,
    gamma: f64,
    phi: f64,
    s: &mut Stiff2D,
    sign: i32,
) -> i32 {
    if vp <= 0. || vs < 0. {
        return -1;
    }
    let mut t1 = vp * vp;
    let t2 = vs * vs;
    s.a3333 = t1;
    s.a1313 = t2;
    s.a1111 = 2. * t1 * eps + t1;
    t1 = 2. * delta * t1 * (t1 - t2) + (t1 - t2) * (t1 - t2);
    if t1 < 0. {
        return -1;
    }
    if sign == 1 {
        s.a1133 = t1.sqrt() - t2
    } else if sign == -1 {
        s.a1133 = -t1.sqrt() - t2
    } else {
        return -1;
    }
    s.a1212 = 2. * t2 * gamma + t2;
    s.a1113 = 0.;
    s.a3313 = 0.;
    s.a1223 = 0.;
    s.a2323 = t2;
    rottens2D(s, phi);
    1
}
/// Convert density-normalized VTI stiffnesses to Thomsen parameters.
/// Returns -1 without writing outputs when a1111 <= 0 or a1313 < 0; otherwise
/// returns 1, including singular inputs that produce NaNs or infinities.
///
/// C: `int stiff2thomVTI(double a1111, double a3333, double a1133, double a1313, double a1212, double *vp, double *vs, double *eps, double *delta, double *gamma)`.
pub fn stiff2thomVTI(
    a: f64,
    c: f64,
    f: f64,
    l: f64,
    o: f64,
    vp: &mut f64,
    vs: &mut f64,
    eps: &mut f64,
    delta: &mut f64,
    gamma: &mut f64,
) -> i32 {
    if a <= 0. || l < 0. {
        return -1;
    }
    *vp = c.sqrt();
    *vs = l.sqrt();
    *eps = (a - c) * 0.5 / c;
    let t = (f + l) * (f + l) - (c - l) * (c - l);
    *delta = t / (c * (c - l)) * 0.5;
    *gamma = (o - l) / l * 0.5;
    1
}
/// Convert stiffnesses to the equivalent VTI Thomsen parameters.
/// Returns 1, including singular inputs that produce NaNs or infinities.
///
/// C: `int stiff2tv(Stiff2D *spar, double *alpha, double *beta, double *ev, double *dv, double *gv)`.
pub fn stiff2tv(
    s: &Stiff2D,
    alpha: &mut f64,
    beta: &mut f64,
    ev: &mut f64,
    dv: &mut f64,
    gv: &mut f64,
) -> i32 {
    *alpha = s.a3333.sqrt();
    *beta = s.a2323.sqrt();
    *ev = (s.a1111 - s.a3333) / (2. * s.a3333);
    *dv = ((s.a1133 + s.a1313) * (s.a1133 + s.a1313) - (s.a3333 - s.a1313) * (s.a3333 - s.a1313))
        / (2. * s.a3333 * (s.a3333 - s.a1313));
    *gv = (s.a1212 - s.a2323) / (2. * s.a2323);
    1
}
/// Convert generic Thomsen parameters to those of the equivalent VTI medium.
/// Returns 1, including singular inputs that produce NaNs or infinities.
///
/// C: `int thom2tv(double vp, double vs, double eps, double delta, double gamma, double *alpha, double *beta, double *ev, double *dv, double *gv)`.
pub fn thom2tv(
    vp: f64,
    vs: f64,
    eps: f64,
    delta: f64,
    gamma: f64,
    alpha: &mut f64,
    beta: &mut f64,
    ev: &mut f64,
    dv: &mut f64,
    gv: &mut f64,
) -> i32 {
    *ev = -eps / (1. + 2. * eps);
    *gv = -gamma / (1. + 2. * gamma);
    *alpha = vp * (1. + 2. * eps).sqrt();
    *beta = vs * (1. + 2. * gamma).sqrt();
    let f = vs / vp;
    let f = 1. - f * f;
    *dv = (delta - 2. * eps * (1. + eps / f)) / ((1. + 2. * eps) * (1. + 2. * eps / f));
    1
}
/// Compute phase velocity from the sine of phase angle in a VTI medium.
/// Mode 0 selects P; every other mode selects SV. A negative discriminant
/// outside the C roundoff threshold returns -1 without changing the output.
///
/// C: `int v_phase2DVTI(Stiff2D *spar1, double sangle, int mode, double *v_phase)`.
pub fn v_phase2DVTI(a: &Stiff2D, sangle: f64, mode: i32, out: &mut f64) -> i32 {
    let s2 = sangle * sangle;
    let c2 = 1. - s2;
    let mut sqr = (a.a1111 - a.a1313) * s2 - (a.a3333 - a.a1313) * c2;
    sqr *= sqr;
    sqr += 4. * (a.a1133 + a.a1313) * (a.a1133 + a.a1313) * s2 * c2;
    if sqr.abs() < (f32::EPSILON as f64) {
        sqr = 0.
    } else if sqr < 0. {
        fprintf!(stderr(), " error in v_phase \n");
        return -1;
    }
    let sqr = sqr.sqrt();
    let dummy = (a.a1111 + a.a1313) * s2 + (a.a3333 + a.a1313) * c2;
    *out = if mode == 0 {
        (0.5 * (dummy + sqr)).sqrt()
    } else {
        (0.5 * (dummy - sqr)).sqrt()
    };
    1
}

/// Compute a Graebner reflection or transmission coefficient for two TI media.
/// Modes 0 and 1 select P and SV; `rort` is 1 for reflection and 0 for
/// transmission. Returns 1 on success or -1 on a detected error. Invalid mode
/// flags can return -1 after writing the reciprocal denominator to `coeff`.
///
/// C: `int graebner2D(Stiff2D *spar1, double rho1, Stiff2D *spar2, double rho2, double pl, int modei, int modet, int rort, double *coeff)`.
#[allow(unused_assignments, unused_variables, unused_mut, unused_parens)]
pub fn graebner2D(
    spar1: &Stiff2D,
    rho1: f64,
    spar2: &Stiff2D,
    rho2: f64,
    pl: f64,
    modei: i32,
    modet: i32,
    rort: i32,
    coeff: &mut f64,
) -> i32 {
    let mut K1: f64 = 0.0;
    let mut K2: f64 = 0.0;
    let mut K3: f64 = 0.0;

    let mut qa1: f64 = 0.0;
    let mut qa2: f64 = 0.0;
    let mut qb1: f64 = 0.0;
    let mut qb2: f64 = 0.0;
    let mut sqr: f64 = 0.0;
    let mut p2: f64 = 0.0;

    let mut la1: f64 = 0.0;
    let mut la2: f64 = 0.0;
    let mut lb1: f64 = 0.0;
    let mut lb2: f64 = 0.0;

    let mut ma1: f64 = 0.0;
    let mut ma2: f64 = 0.0;
    let mut mb1: f64 = 0.0;
    let mut mb2: f64 = 0.0;

    let mut a1: f64 = 0.0;
    let mut a2: f64 = 0.0;
    let mut b1: f64 = 0.0;
    let mut b2: f64 = 0.0;

    let mut c1: f64 = 0.0;
    let mut c2: f64 = 0.0;
    let mut d1: f64 = 0.0;
    let mut d2: f64 = 0.0;

    let mut qasq: f64 = 0.0;
    let mut qbsq: f64 = 0.0;

    let mut a1111: f64 = 0.0;
    let mut a3333: f64 = 0.0;
    let mut a1133: f64 = 0.0;

    let mut a1313: f64 = 0.0;
    let mut oa3333: f64 = 0.0;
    let mut oa1313: f64 = 0.0;

    let mut m11: f64 = 0.0;
    let mut m12: f64 = 0.0;
    let mut m13: f64 = 0.0;
    let mut m14: f64 = 0.0;
    let mut m21: f64 = 0.0;
    let mut m22: f64 = 0.0;
    let mut m23: f64 = 0.0;
    let mut m24: f64 = 0.0;

    let mut m31: f64 = 0.0;
    let mut m32: f64 = 0.0;
    let mut m33: f64 = 0.0;
    let mut m34: f64 = 0.0;
    let mut m41: f64 = 0.0;
    let mut m42: f64 = 0.0;
    let mut m43: f64 = 0.0;
    let mut m44: f64 = 0.0;

    p2 = pl * pl;

    a1111 = spar1.a1111;
    a3333 = spar1.a3333;
    a1133 = spar1.a1133;
    a1313 = spar1.a1313;

    oa3333 = 1.0 / a3333;
    oa1313 = 1.0 / a1313;

    sqr = a1111 * oa1313 + a1313 * oa3333 - (a1313 + a1133) * (a1313 + a1133) * oa3333 * oa1313;
    K1 = oa3333 + oa1313 - sqr * p2;

    K2 = a1111 * oa3333 * p2 - oa3333;

    K3 = p2 - oa1313;

    sqr = K1 * K1 - 4.0 * K2 * K3;
    if sqr.abs() < (f32::EPSILON as f64) {
        sqr = 0.0;
    }

    if sqr < 0.0 {
        fprintf!(stderr(), "ERROR in slowness medium1 \n");
        return -1;
    }

    qasq = 0.5 * (K1 - sqr.sqrt());
    qa1 = qasq.sqrt();
    qbsq = 0.5 * (K1 + sqr.sqrt());
    qb1 = qbsq.sqrt();

    sqr = a1111 * p2 + a1313 * qasq + a3333 * qasq + a1313 * p2 - 2.0;
    if sqr.abs() < (f32::EPSILON as f64) {
        fprintf!(stderr(), "ERROR P eigenvector medium 1 \n");
        return -1;
    }

    sqr = 1.0 / sqr;

    la1 = a3333 * qasq + a1313 * p2 - 1.0;

    la1 = la1 * sqr;
    if la1.abs() < (f32::EPSILON as f64) {
        la1 = 0.0;
    }

    if la1 < 0.0 {
        fprintf!(stderr(), "ERROR la1 < 0\n");
        return -1;
    }

    la1 = la1.sqrt();

    ma1 = a1313 * qasq + a1111 * p2 - 1.0;
    ma1 = ma1 * sqr;
    if ma1.abs() < (f32::EPSILON as f64) {
        ma1 = 0.0;
    }

    if ma1 < 0.0 {
        fprintf!(stderr(), "ERROR ma1 <0\n");
        return -1;
    }

    ma1 = ma1.sqrt();

    sqr = a1111 * p2 + a1313 * qbsq + a3333 * qbsq + a1313 * p2 - 2.0;
    if sqr.abs() < (f32::EPSILON as f64) {
        fprintf!(stderr(), "ERROR SV eigenvector medium 1\n");
        return -1;
    }

    sqr = 1.0 / sqr;

    lb1 = a1313 * qbsq + a1111 * p2 - 1.0;

    lb1 = lb1 * sqr;

    if lb1.abs() < (f32::EPSILON as f64) {
        lb1 = 0.0;
    }

    if lb1 < 0.0 {
        fprintf!(stderr(), "ERROR lb1 <0\n");
        return -1;
    }

    lb1 = lb1.sqrt();

    mb1 = a3333 * qbsq + a1313 * p2 - 1.0;
    mb1 = mb1 * sqr;

    if mb1.abs() < (f32::EPSILON as f64) {
        mb1 = 0.0;
    }

    if mb1 < 0.0 {
        fprintf!(stderr(), "ERROR mb1 <0\n");
        return -1;
    }

    mb1 = mb1.sqrt();

    a1313 = a1313 * rho1;
    a3333 = a3333 * rho1;
    a1133 = a1133 * rho1;

    a1 = a1313 * (qa1 * la1 + pl * ma1);
    b1 = a1313 * (qb1 * mb1 - pl * lb1);
    c1 = pl * la1 * a1133 + qa1 * ma1 * a3333;
    d1 = pl * mb1 * a1133 - qb1 * lb1 * a3333;

    a1111 = spar2.a1111;
    a3333 = spar2.a3333;
    a1133 = spar2.a1133;
    a1313 = spar2.a1313;

    oa3333 = 1.0 / a3333;
    oa1313 = 1.0 / a1313;

    sqr = a1111 * oa1313 + a1313 * oa3333 - (a1313 + a1133) * (a1313 + a1133) * oa3333 * oa1313;
    K1 = oa3333 + oa1313 - sqr * p2;

    K2 = a1111 * oa3333 * p2 - oa3333;

    K3 = p2 - oa1313;

    sqr = K1 * K1 - 4.0 * K2 * K3;
    if sqr.abs() < (f32::EPSILON as f64) {
        sqr = 0.0;
    }

    if sqr < 0.0 {
        fprintf!(stderr(), "ERROR slowness medium 2 \n");
        return -1;
    }

    qasq = 0.5 * (K1 - sqr.sqrt());

    if qasq.abs() < (f32::EPSILON as f64) {
        qasq = 0.0;
    }

    if qasq < 0.0 {
        fprintf!(stderr(), "\n COMPLEX vertical P- slowness medium 2 \n");
        return -1;
    }

    qa2 = qasq.sqrt();
    qbsq = 0.5 * (K1 + sqr.sqrt());
    qb2 = qbsq.sqrt();

    sqr = a1111 * p2 + a1313 * qasq + a3333 * qasq + a1313 * p2 - 2.0;
    if sqr.abs() < (f32::EPSILON as f64) {
        fprintf!(stderr(), "ERROR P eigenvector medium 2 \n");
        return -1;
    }

    sqr = 1.0 / sqr;

    la2 = a3333 * qasq + a1313 * p2 - 1.0;

    la2 = la2 * sqr;
    if la2.abs() < (f32::EPSILON as f64) {
        la2 = 0.0;
    }

    if la2 < 0.0 {
        fprintf!(stderr(), "ERROR la2<0\n");
        return -1;
    }

    la2 = la2.sqrt();

    ma2 = a1313 * qasq + a1111 * p2 - 1.0;
    ma2 = ma2 * sqr;

    if ma2.abs() < (f32::EPSILON as f64) {
        ma2 = 0.0;
    }

    if ma2 < 0.0 {
        fprintf!(stderr(), "ERROR ma2<0\n");
        return -1;
    }

    ma2 = ma2.sqrt();

    sqr = a1111 * p2 + a1313 * qbsq + a3333 * qbsq + a1313 * p2 - 2.0;
    if sqr.abs() < (f32::EPSILON as f64) {
        fprintf!(stderr(), "ERROR SV eigenvector medium 2 \n");
        return -1;
    }

    sqr = 1.0 / sqr;

    lb2 = a1313 * qbsq + a1111 * p2 - 1.0;
    lb2 = lb2 * sqr;
    if lb2.abs() < (f32::EPSILON as f64) {
        lb2 = 0.0;
    }

    if lb2 < 0.0 {
        fprintf!(stderr(), "ERROR lb2<0\n");
        return -1;
    }

    lb2 = lb2.sqrt();

    mb2 = a3333 * qbsq + a1313 * p2 - 1.0;
    mb2 = mb2 * sqr;
    if mb2.abs() < (f32::EPSILON as f64) {
        mb2 = 0.0;
    }

    if mb2 < 0.0 {
        fprintf!(stderr(), "ERROR mb2<0\n");
        return -1;
    }

    mb2 = mb2.sqrt();

    a1313 = a1313 * rho2;
    a3333 = a3333 * rho2;
    a1133 = a1133 * rho2;

    a2 = a1313 * (qa2 * la2 + pl * ma2);
    b2 = a1313 * (qb2 * mb2 - pl * lb2);
    c2 = pl * la2 * a1133 + qa2 * ma2 * a3333;
    d2 = pl * mb2 * a1133 - qb2 * lb2 * a3333;

    m11 = la1;
    m12 = mb1;
    m13 = -la2;
    m14 = -mb2;
    m21 = c1;
    m22 = d1;
    m23 = -c2;
    m24 = -d2;
    m31 = ma1;
    m32 = -lb1;
    m33 = ma2;
    m34 = -lb2;
    m41 = a1;
    m42 = b1;
    m43 = a2;
    m44 = b2;

    sqr = (m14 * m23 * m32 * m41 - m13 * m24 * m32 * m41 - m14 * m22 * m33 * m41
        + m12 * m24 * m33 * m41
        + m13 * m22 * m34 * m41
        - m12 * m23 * m34 * m41
        - m14 * m23 * m31 * m42
        + m13 * m24 * m31 * m42
        + m14 * m21 * m33 * m42
        - m11 * m24 * m33 * m42
        - m13 * m21 * m34 * m42
        + m11 * m23 * m34 * m42
        + m14 * m22 * m31 * m43
        - m12 * m24 * m31 * m43
        - m14 * m21 * m32 * m43
        + m11 * m24 * m32 * m43
        + m12 * m21 * m34 * m43
        - m11 * m22 * m34 * m43
        - m13 * m22 * m31 * m44
        + m12 * m23 * m31 * m44
        + m13 * m21 * m32 * m44
        - m11 * m23 * m32 * m44
        - m12 * m21 * m33 * m44
        + m11 * m22 * m33 * m44);

    if sqr.abs() < (f32::EPSILON as f64) {
        fprintf!(stderr(), "ERROR Denominator == 0 \n");
        return -1;
    }

    *coeff = 1.0 / sqr;

    if modei == 0 && modet == 0 && rort == 0 {
        *coeff *= 2.0
            * (-(m14 * m21 * m32 * m41) + m11 * m24 * m32 * m41 + m12 * m21 * m34 * m41
                - m11 * m22 * m34 * m41
                + m14 * m21 * m31 * m42
                - m11 * m24 * m31 * m42
                - m12 * m21 * m31 * m44
                + m11 * m22 * m31 * m44);
    } else if modei == 0 && modet == 0 && rort == 1 {
        *coeff *= (m14 * m23 * m32 - m13 * m24 * m32 - m14 * m22 * m33
            + m12 * m24 * m33
            + m13 * m22 * m34
            - m12 * m23 * m34)
            * m41
            + m31
                * (-(m14 * m23 * m42) + m13 * m24 * m42 + m14 * m22 * m43
                    - m12 * m24 * m43
                    - m13 * m22 * m44
                    + m12 * m23 * m44)
            + m21
                * (-(m14 * m33 * m42) + m13 * m34 * m42 + m14 * m32 * m43
                    - m12 * m34 * m43
                    - m13 * m32 * m44
                    + m12 * m33 * m44)
            + m11
                * (m24 * m33 * m42 - m23 * m34 * m42 - m24 * m32 * m43
                    + m22 * m34 * m43
                    + m23 * m32 * m44
                    - m22 * m33 * m44);
    } else if modei == 0 && modet == 1 && rort == 0 {
        *coeff *= 2.0
            * (m13 * m21 * m32 * m41 - m11 * m23 * m32 * m41 - m12 * m21 * m33 * m41
                + m11 * m22 * m33 * m41
                - m13 * m21 * m31 * m42
                + m11 * m23 * m31 * m42
                + m12 * m21 * m31 * m43
                - m11 * m22 * m31 * m43);
    } else if modei == 0 && modet == 1 && rort == 1 {
        *coeff *= 2.0
            * (m14 * m21 * m33 * m41 - m11 * m24 * m33 * m41 - m13 * m21 * m34 * m41
                + m11 * m23 * m34 * m41
                - m14 * m21 * m31 * m43
                + m11 * m24 * m31 * m43
                + m13 * m21 * m31 * m44
                - m11 * m23 * m31 * m44);
    } else if modei == 1 && modet == 0 && rort == 0 {
        *coeff *= 2.0
            * (-(m14 * m22 * m32 * m41) + m12 * m24 * m32 * m41 + m14 * m22 * m31 * m42
                - m12 * m24 * m31 * m42
                + m12 * m21 * m34 * m42
                - m11 * m22 * m34 * m42
                - m12 * m21 * m32 * m44
                + m11 * m22 * m32 * m44);
    } else if modei == 1 && modet == 0 && rort == 1 {
        *coeff *= 2.0
            * (-(m14 * m22 * m33 * m42) + m12 * m24 * m33 * m42 + m13 * m22 * m34 * m42
                - m12 * m23 * m34 * m42
                + m14 * m22 * m32 * m43
                - m12 * m24 * m32 * m43
                - m13 * m22 * m32 * m44
                + m12 * m23 * m32 * m44);
    } else if modei == 1 && modet == 1 && rort == 0 {
        *coeff *= 2.0
            * (m13 * m22 * m32 * m41 - m12 * m23 * m32 * m41 - m13 * m22 * m31 * m42
                + m12 * m23 * m31 * m42
                - m12 * m21 * m33 * m42
                + m11 * m22 * m33 * m42
                + m12 * m21 * m32 * m43
                - m11 * m22 * m32 * m43);
    } else if modei == 1 && modet == 1 && rort == 1 {
        *coeff *= (-(m14 * m23 * m31) + m13 * m24 * m31 + m14 * m21 * m33
            - m11 * m24 * m33
            - m13 * m21 * m34
            + m11 * m23 * m34)
            * m42
            + m32
                * (m14 * m23 * m41 - m13 * m24 * m41 - m14 * m21 * m43
                    + m11 * m24 * m43
                    + m13 * m21 * m44
                    - m11 * m23 * m44)
            + m22
                * (m14 * m33 * m41 - m13 * m34 * m41 - m14 * m31 * m43
                    + m11 * m34 * m43
                    + m13 * m31 * m44
                    - m11 * m33 * m44)
            + m12
                * (-(m24 * m33 * m41) + m23 * m34 * m41 + m24 * m31 * m43
                    - m21 * m34 * m43
                    - m23 * m31 * m44
                    + m21 * m33 * m44);
    } else {
        return -1;
    }

    sqr = *coeff;

    return 1;
}
