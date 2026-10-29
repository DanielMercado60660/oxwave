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

//! Port of Trielas/lib/rt_ani_real.c and its private support routines in
//! Trielas/main/elaray.c. Author: Andreas Rueger, Colorado School of Mines,
//! 04/24/94 (coefficient), 02/01/94 (root finder), 03/14/94 (polarization).

#![allow(non_snake_case, unused_mut, unused_parens, clippy::too_many_arguments)]
use super::{rottensor, rotvector, solveForSlowqP, solveForSlowqSV};
use oxwave_core::cio::CFile;
use oxwave_lu::{sgeco, sgesl};

// Private support from elaray.c, preserving its polarization conventions.
fn polar(
    px: f32,
    pz: f32,
    g11: f32,
    g13: f32,
    g33: f32,
    polx: &mut f32,
    polz: &mut f32,
    mode: i32,
) {
    let (mut x, mut z);
    if g11 < f32::EPSILON {
        x = 0.0;
        z = (g33 as f64).sqrt() as f32;
    } else if g33 < f32::EPSILON {
        z = 0.0;
        x = (g11 as f64).sqrt() as f32;
    } else {
        x = ((g11 as f64).sqrt() * if g13 < 0.0 { -1.0 } else { 1.0 }) as f32;
        z = (g33 as f64).sqrt() as f32;
    }
    if ((mode == 0 || mode == 3) && px * x + pz * z < 0.0)
        || ((mode == 4 || mode == 1) && px * x < 0.0)
    {
        x = -x;
        z = -z;
    }
    *polx = x;
    *polz = z;
}

fn gvelreal(
    a: f32,
    b: f32,
    f: f32,
    l: f32,
    m: f32,
    n: f32,
    x: f32,
    z: f32,
    vx: &mut f32,
    vz: &mut f32,
    g11n: &mut f32,
    g13n: &mut f32,
    g33n: &mut f32,
) {
    let x2 = x * x;
    let z2 = z * z;
    let xz = x * z;
    let gamm11 = a * x2 + l * z2 + 2.0 * m * xz;
    let gamm33 = b * z2 + l * x2 + 2.0 * n * xz;
    let gamm13 = (f + l) * xz + m * x2 + n * z2;
    let den = 1.0 / (gamm11 + gamm33 - 2.0);
    let mut g11 = (gamm33 - 1.0) * den;
    let mut g33 = (gamm11 - 1.0) * den;
    let g13 = -gamm13 * den;
    *vx =
        a * x * g11 + (f + l) * z * g13 + n * z * g33 + m * (z * g11 + 2.0 * x * g13) + l * g33 * x;
    *vz =
        b * z * g33 + (f + l) * x * g13 + m * x * g11 + n * (x * g33 + 2.0 * z * g13) + l * g11 * z;
    if g11 < -10.0 * f32::EPSILON || g33 < -10.0 * f32::EPSILON {
        oxwave_core::err!("\n ERROR: g11 or g33 <0 in initialization \n");
    } else if g11 < 0.0 {
        g11 = 0.0;
    } else if g33 < 0.0 {
        g33 = 0.0;
    }
    *g11n = g11;
    *g13n = g13;
    *g33n = g33;
}

fn findqPqSV(
    s: f32,
    c: f32,
    pl: f32,
    a: f32,
    b: f32,
    f: f32,
    l: f32,
    m: f32,
    n: f32,
    mode: i32,
    x: &mut f32,
    z: &mut f32,
    vx: &mut f32,
    vz: &mut f32,
    g11: &mut f32,
    g13: &mut f32,
    g33: &mut f32,
    rt: f32,
    ifp: Option<&CFile>,
) -> i32 {
    let solver = if mode == 3 {
        solveForSlowqP
    } else {
        solveForSlowqSV
    };
    let mut check1 = solver(s, c, pl, a, b, f, l, m, n, x, z, rt);
    let mut check2 = false;
    if check1 == 1 {
        gvelreal(a, b, f, l, m, n, *x, *z, vx, vz, g11, g13, g33);
        let vgm = -s * *vx + c * *vz;
        check2 = (rt == 1.0 && vgm > 0.0) || (rt == -1.0 && vgm < 0.0);
        if check2 {
            if let Some(fp) = ifp {
                oxwave_core::fprintf!(fp, "\n first try in rootfinder was successfull \n");
            }
        }
    }
    if check1 != 1 || !check2 {
        if let Some(fp) = ifp {
            oxwave_core::fprintf!(fp, "\n need second try in rootfinder \n");
        }
        check1 = solver(s, c, pl, a, b, f, l, m, n, x, z, -1.0 * rt);
        if check1 != 1 {
            // The C helper dereferences a NULL stream on this failure path.
            // A missing diagnostic stream does not affect its return value.
            if let Some(fp) = ifp {
                oxwave_core::fprintf!(fp, " rootfinder not sucessful !!! \n");
            }
            return 0;
        } else {
            gvelreal(a, b, f, l, m, n, *x, *z, vx, vz, g11, g13, g33);
            let vgm = -s * *vx + c * *vz;
            if (rt == 1.0 && vgm < 0.0) || (rt == -1.0 && vgm > 0.0) {
                return 0;
            }
            check2 = true;
            if let Some(fp) = ifp {
                oxwave_core::fprintf!(fp, "\n second try in rootfinder was successfull \n");
            }
        }
    }
    let newpl = c * *x + *z * s;
    if (((newpl - pl) * (newpl - pl)) as f64) < 0.000000001 && check2 && check1 == 1 {
        1
    } else {
        0
    }
}

/// Evaluate anisotropic qP/qSV reflection or transmission, including a free surface.
/// The optional diagnostic stream is None for a C NULL stream.
///
/// C prototype: int rt_ani_real(float gz, float gx, float a1111i, float a3333i,
/// float a1133i, float a1313i, float a1113i, float a3313i, float rhoi,
/// float a1111t, float a3333t, float a1133t, float a1313t, float a1113t,
/// float a3313t, float rhot, float pxi, float pzi, float g11i, float g13i,
/// float g33i, float pxg, float pzg, float g11g, float g13g, float g33g,
/// int modei, int modeg, int rort, float *coeff, FILE *ifp).
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
    mut pxi: f32,
    mut pzi: f32,
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
    coeff: &mut f32,
    ifp: Option<&CFile>,
) -> i32 {
    let index: usize;
    let mut c1111i: f32;
    let mut c3333i: f32;
    let mut c1133i: f32;
    let mut c1313i: f32;
    let mut c1113i: f32;
    let mut c3313i: f32;
    let mut c1111t: f32;
    let mut c3333t: f32;
    let mut c1133t: f32;
    let mut c1313t: f32;
    let mut c1113t: f32;
    let mut c3313t: f32;
    let mut pxrP: f32 = 0.0;
    let mut pzrP: f32 = 0.0;
    let mut pxtP: f32 = 0.0;
    let mut pztP: f32 = 0.0;
    let mut pxrS: f32 = 0.0;
    let mut pzrS: f32 = 0.0;
    let mut pxtS: f32 = 0.0;
    let mut pztS: f32 = 0.0;
    let mut g11: f32 = 0.0;
    let mut g13: f32 = 0.0;
    let mut g33: f32 = 0.0;
    let mut vgx: f32 = 0.0;
    let mut vgz: f32 = 0.0;
    let mut A: f32;
    let mut B: f32;
    let mut C: f32;
    let mut D: f32;
    let mut H: f32;
    let mut I: f32;
    let mut den: f32;
    let mut polxi: f32 = 0.0;
    let mut polzi: f32 = 0.0;
    let mut polxrP: f32 = 0.0;
    let mut polzrP: f32 = 0.0;
    let mut polxtP: f32 = 0.0;
    let mut polztP: f32 = 0.0;
    let mut polxrS: f32 = 0.0;
    let mut polzrS: f32 = 0.0;
    let mut polxtS: f32 = 0.0;
    let mut polztS: f32 = 0.0;
    let mut a = [[0.0f32; 4]; 4];
    let mut b = [0.0f32; 4];
    let mut ipvt = [0i32; 4];
    let mut z = [0.0f32; 4];
    let mut rcond = 0.0f32;

    if ((modeg == 0 || modeg == 3) && (rort == 1 || rort == 2)) {
        index = 0;
    } else if ((modeg == 1 || modeg == 4) && (rort == 1 || rort == 2)) {
        index = 1;
    } else if ((modeg == 0 || modeg == 3) && rort == 0) {
        index = 2;
    } else if ((modeg == 1 || modeg == 4) && rort == 0) {
        index = 3;
    } else {
        oxwave_core::err!(" wrong mode in rt_ani_real");
    }

    polar(pxi, pzi, g11i, g13i, g33i, &mut polxi, &mut polzi, modei);

    rotvector(&mut pxi, &mut pzi, gz, gx);

    if (index != 0) {
        if (findqPqSV(
            gz, gx, pxi, a1111i, a3333i, a1133i, a1313i, a1113i, a3313i, 3, &mut pxrP, &mut pzrP,
            &mut vgx, &mut vgz, &mut g11, &mut g13, &mut g33, -1.0, ifp,
        ) != 1)
        {
            return 0;
        }
    } else {
        pxrP = pxg;
        pzrP = pzg;
        g11 = g11g;
        g13 = g13g;
        g33 = g33g;
    }

    polar(pxrP, pzrP, g11, g13, g33, &mut polxrP, &mut polzrP, 3);

    if (index != 1) {
        if (findqPqSV(
            gz, gx, pxi, a1111i, a3333i, a1133i, a1313i, a1113i, a3313i, 4, &mut pxrS, &mut pzrS,
            &mut vgx, &mut vgz, &mut g11, &mut g13, &mut g33, -1.0, ifp,
        ) != 1)
        {
            return 0;
        }
    } else {
        pxrS = pxg;
        pzrS = pzg;
        g11 = g11g;
        g13 = g13g;
        g33 = g33g;
    }

    polar(pxrS, pzrS, g11, g13, g33, &mut polxrS, &mut polzrS, 4);

    rotvector(&mut pxrP, &mut pzrP, gz, gx);
    rotvector(&mut pxrS, &mut pzrS, gz, gx);
    rotvector(&mut polxi, &mut polzi, gz, gx);
    rotvector(&mut polxrP, &mut polzrP, gz, gx);
    rotvector(&mut polxrS, &mut polzrS, gz, gx);

    c1111i = rhoi * a1111i;
    c3333i = rhoi * a3333i;
    c1313i = rhoi * a1313i;
    c1133i = rhoi * a1133i;
    c1113i = rhoi * a1113i;
    c3313i = rhoi * a3313i;

    rottensor(
        &mut c1111i,
        &mut c3333i,
        &mut c1133i,
        &mut c1313i,
        &mut c1113i,
        &mut c3313i,
        gz,
        gx,
    );

    if (rort == 2) {
        A = polxrP * pxrP * c1113i
            + c1313i * (polxrP * pzrP + polzrP * pxrP)
            + c3313i * polzrP * pzrP;
        B = polxrS * pxrS * c1113i
            + c1313i * (polxrS * pzrS + polzrS * pxrS)
            + c3313i * polzrS * pzrS;
        C = polxrP * pxrP * c1133i
            + c3313i * (polxrP * pzrP + polzrP * pxrP)
            + c3333i * polzrP * pzrP;
        D = polxrS * pxrS * c1133i
            + c3313i * (polxrS * pzrS + polzrS * pxrS)
            + c3333i * polzrS * pzrS;
        H = polxi * pxi * c1113i + c1313i * (polxi * pzi + polzi * pxi) + c3313i * polzi * pzi;
        I = polxi * pxi * c1133i + c3313i * (polxi * pzi + polzi * pxi) + c3333i * polzi * pzi;
        den = D * A - C * B;

        if (den.abs() < f32::EPSILON) {
            oxwave_core::warn!(" problems in rt_ani_real (free surface)");
            return -1;
        } else {
            den = (C * H - I * A) / den;
        }

        if (index == 1) {
            *coeff = den;
        } else if (A.abs() < f32::EPSILON) {
            oxwave_core::warn!(" problems in rt_ani_real (free surface)");
            return -1;
        } else if (index == 0) {
            *coeff = -(H + B * den) / A;
        } else {
            oxwave_core::warn!(" problems in rt_ani_real (free surface)");
            return -1;
        }

        return 1;
    }

    if (index != 2) {
        if (findqPqSV(
            gz, gx, pxi, a1111t, a3333t, a1133t, a1313t, a1113t, a3313t, 3, &mut pxtP, &mut pztP,
            &mut vgx, &mut vgz, &mut g11, &mut g13, &mut g33, 1.0, ifp,
        ) != 1)
        {
            return 0;
        }
    } else {
        pxtP = pxg;
        pztP = pzg;
        g11 = g11g;
        g13 = g13g;
        g33 = g33g;
    }

    polar(pxtP, pztP, g11, g13, g33, &mut polxtP, &mut polztP, 3);

    if (index != 3) {
        if (findqPqSV(
            gz, gx, pxi, a1111t, a3333t, a1133t, a1313t, a1113t, a3313t, 4, &mut pxtS, &mut pztS,
            &mut vgx, &mut vgz, &mut g11, &mut g13, &mut g33, 1.0, ifp,
        ) != 1)
        {
            return 0;
        }
    } else {
        pxtS = pxg;
        pztS = pzg;
        g11 = g11g;
        g13 = g13g;
        g33 = g33g;
    }

    polar(pxtS, pztS, g11, g13, g33, &mut polxtS, &mut polztS, 4);

    rotvector(&mut pxtP, &mut pztP, gz, gx);
    rotvector(&mut pxtS, &mut pztS, gz, gx);

    rotvector(&mut polxtP, &mut polztP, gz, gx);
    rotvector(&mut polxtS, &mut polztS, gz, gx);

    if (pxrP * polxrP + pzrP * polzrP < 0.0) {
        polxrP = -polxrP;
        polzrP = -polzrP;
    }
    if (pxtP * polxtP + pztP * polztP < 0.0) {
        polxtP = -polxtP;
        polztP = -polztP;
    }
    if (pxrS * polxrS < 0.0) {
        polxrS = -polxrS;
        polzrS = -polzrS;
    }
    if (pxtS * polxtS < 0.0) {
        polxtS = -polxtS;
        polztS = -polztS;
    }

    c1111t = rhot * a1111t;
    c3333t = rhot * a3333t;
    c1313t = rhot * a1313t;
    c1133t = rhot * a1133t;
    c1113t = rhot * a1113t;
    c3313t = rhot * a3313t;

    rottensor(
        &mut c1111t,
        &mut c3333t,
        &mut c1133t,
        &mut c1313t,
        &mut c1113t,
        &mut c3313t,
        gz,
        gx,
    );

    a[0][0] = polxrP;
    a[1][0] = polxrS;
    a[2][0] = -polxtP;
    a[3][0] = -polxtS;
    a[0][1] = polzrP;
    a[1][1] = polzrS;
    a[2][1] = -polztP;
    a[3][1] = -polztS;

    a[0][2] = c1133i * pxrP * polxrP;
    a[0][2] += c3313i * (pzrP * polxrP + pxrP * polzrP);
    a[0][2] += c3333i * pzrP * polzrP;

    a[1][2] = c1133i * pxrS * polxrS;
    a[1][2] += c3313i * (pzrS * polxrS + polzrS * pxrS);
    a[1][2] += c3333i * pzrS * polzrS;

    a[2][2] = c1133t * pxtP * polxtP;
    a[2][2] += c3313t * (pztP * polxtP + pxtP * polztP);
    a[2][2] += c3333t * pztP * polztP;
    a[2][2] = -a[2][2];

    a[3][2] = c1133t * pxtS * polxtS;
    a[3][2] += c3313t * (pztS * polxtS + pxtS * polztS);
    a[3][2] += c3333t * pztS * polztS;
    a[3][2] = -a[3][2];

    a[0][3] = c1113i * pxrP * polxrP;
    a[0][3] += c1313i * (pzrP * polxrP + pxrP * polzrP);
    a[0][3] += c3313i * pzrP * polzrP;

    a[1][3] = c1113i * pxrS * polxrS;
    a[1][3] += c1313i * (pzrS * polxrS + pxrS * polzrS);
    a[1][3] += c3313i * pzrS * polzrS;

    a[2][3] = c1113t * pxtP * polxtP;
    a[2][3] += c1313t * (pztP * polxtP + pxtP * polztP);
    a[2][3] += c3313t * pztP * polztP;
    a[2][3] = -a[2][3];

    a[3][3] = c1113t * pxtS * polxtS;
    a[3][3] += c1313t * (pztS * polxtS + pxtS * polztS);
    a[3][3] += c3313t * pztS * polztS;
    a[3][3] = -a[3][3];

    b[0] = -polxi;

    b[1] = -polzi;

    b[2] = c1133i * polxi * pxi;
    b[2] += c3313i * (polxi * pzi + polzi * pxi);
    b[2] += c3333i * polzi * pzi;
    b[2] = -b[2];

    b[3] = c1113i * polxi * pxi;
    b[3] += c1313i * (polxi * pzi + polzi * pxi);
    b[3] += c3313i * polzi * pzi;
    b[3] = -b[3];

    sgeco(
        &mut a.each_mut().map(|v| v.as_mut_slice()),
        4,
        &mut ipvt,
        &mut rcond,
        &mut z,
    );
    sgesl(&a.each_ref().map(|v| v.as_slice()), 4, &ipvt, &mut b, 0);

    if let Some(ifp) = ifp {
        oxwave_core::fprintf!(ifp, "\n Anisotropic Reflection/Transmission coeff\n");
        oxwave_core::fprintf!(ifp, " inc slowness: %f %f\n", pxi, pzi);
        oxwave_core::fprintf!(ifp, " inc polariz.: %f %f\n", polxi, polzi);
        oxwave_core::fprintf!(ifp, " Pref slowness: %f %f\n", pxrP, pzrP);
        oxwave_core::fprintf!(ifp, " Pref polariz.: %f %f\n", polxrP, polzrP);
        oxwave_core::fprintf!(ifp, " Sref slowness %f %f\n", pxrS, pzrS);
        oxwave_core::fprintf!(ifp, " Sref polariz.: %f %f\n", polxrS, polzrS);
        oxwave_core::fprintf!(ifp, " Ptrans slowness %f %f\n", pxtP, pztP);
        oxwave_core::fprintf!(ifp, " Ptrans polariz.: %f %f\n", polxtP, polztP);
        oxwave_core::fprintf!(ifp, " Strans slowness %f %f\n", pxtS, pztS);
        oxwave_core::fprintf!(ifp, " Strans polariz.: %f %f\n", polxtS, polztS);
        oxwave_core::fprintf!(ifp, " Values with respect to rotated system\n");
        oxwave_core::fprintf!(ifp, " \n *****coefficients ******* \n ");
        oxwave_core::fprintf!(ifp, " Pref=%g \n", b[0]);
        oxwave_core::fprintf!(ifp, " Sref=%g \n", b[1]);
        oxwave_core::fprintf!(ifp, " Ptrans=%g \n", b[2]);
        oxwave_core::fprintf!(ifp, " Strans=%g \n", b[3]);
    }

    *coeff = b[index];

    1
}
