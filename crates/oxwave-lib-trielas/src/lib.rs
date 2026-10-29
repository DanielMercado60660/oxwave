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

//! Ports of `Trielas/lib/findPiso.c`, `findSiso.c`, `findqSH.c`, `rottensh.c`,
//! `rottensor.c`, `rotvector.c`, `rt_SHa_real.c`, `rt_ani_real.c`,
//! `rt_iso_cmplx.c`, `rt_iso_real.c`, `solveForSlowqP.c`, and `solveForSlowqSV.c`.
//! Authors: Andreas Rueger, Colorado School of Mines, 03/14/94;
//! Tariq Alkalifah, Colorado School of Mines, 02/02/94.

mod ani;
pub use ani::rt_ani_real;

/// Rotate a two-dimensional elastic stiffness tensor in place.
///
/// C prototype: `void rottensor(float *a1111, float *a3333, float *a1133,
/// float *a1313, float *a1113, float *a3313, float si, float co)`.
#[allow(clippy::too_many_arguments)]
pub fn rottensor(
    a1111: &mut f32,
    a3333: &mut f32,
    a1133: &mut f32,
    a1313: &mut f32,
    a1113: &mut f32,
    a3313: &mut f32,
    si: f32,
    co: f32,
) {
    let ss = si * si;
    let cc = co * co;
    let sc = si * co;
    let a = *a1111;
    let c = *a3333;
    let f = *a1133;
    let l = *a1313;
    let m = *a1113;
    let n = *a3313;

    *a1111 = a * cc * cc + c * ss * ss + 2.0 * f * cc * ss + 4.0 * l * cc * ss
        - 4.0 * m * cc * sc
        - 4.0 * n * ss * sc;
    *a3333 = a * ss * ss
        + c * cc * cc
        + 2.0 * f * cc * ss
        + 4.0 * l * cc * ss
        + 4.0 * m * ss * sc
        + 4.0 * n * cc * sc;
    *a1133 = a * ss * cc + c * ss * cc + f * (ss * ss + cc * cc)
        - 4.0 * l * ss * cc
        - 2.0 * m * (ss * sc - cc * sc)
        - 2.0 * n * (cc * sc - ss * sc);
    *a1313 = a * ss * cc + c * ss * cc - 2.0 * f * ss * cc
        + l * (cc * cc + ss * ss - 2.0 * ss * cc)
        - 2.0 * m * (ss * sc - cc * sc)
        - 2.0 * n * (cc * sc - ss * sc);
    *a1113 = a * cc * sc - c * ss * sc
        + f * (ss * sc - cc * sc)
        + 2.0 * l * (ss * sc - cc * sc)
        + m * (cc * cc - 3.0 * cc * ss)
        + (3.0 * cc * ss - ss * ss) * n;
    *a3313 = a * ss * sc - c * cc * sc
        + f * (cc * sc - ss * sc)
        + 2.0 * l * (cc * sc - ss * sc)
        + m * (3.0 * cc * ss - ss * ss)
        + n * (cc * cc - 3.0 * ss * cc);
}

/// Rotate a vector. C: `void rotvector(float *x, float *z, float si, float co)`.
/// Trielas/lib/rotvector.c; Andreas Rueger, 03/14/94.
pub fn rotvector(x: &mut f32, z: &mut f32, si: f32, co: f32) {
    let o = *x;
    let p = *z;
    *x = o * co + p * si;
    *z = p * co - o * si;
}

/// Rotate an SH tensor. C: `void rottensh(float *a1212, float *a2323, float *a1223, float si, float co)`.
/// Trielas/lib/rottensh.c; Andreas Rueger, 03/14/94.
pub fn rottensh(a1212: &mut f32, a2323: &mut f32, a1223: &mut f32, si: f32, co: f32) {
    let ss = si * si;
    let cc = co * co;
    let sc = si * co;
    let o = *a1212;
    let p = *a2323;
    let q = *a1223;
    *a1212 = cc * o - 2.0 * sc * q + ss * p;
    *a2323 = ss * o + 2.0 * sc * q + cc * p;
    *a1223 = sc * o + (cc - ss) * q - sc * p;
}

/// Find an isotropic P root. C: `int findPiso(float a3333, float pl, float gx, float gz, float *px, float *pz, float *vgx, float *vgz, float *g11, float *g13, float *g33, int r)`.
/// Trielas/lib/findPiso.c; Andreas Rueger, 02/06/94.
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn findPiso(
    a: f32,
    pl: f32,
    gx: f32,
    gz: f32,
    px: &mut f32,
    pz: &mut f32,
    vgx: &mut f32,
    vgz: &mut f32,
    g11: &mut f32,
    g13: &mut f32,
    g33: &mut f32,
    r: i32,
) -> i32 {
    find_iso(a, pl, gx, gz, px, pz, vgx, vgz, g11, g13, g33, r, false)
}

/// Find an isotropic S root. C: `int findSiso(float a1313, float pl, float gx, float gz, float *px, float *pz, float *vgx, float *vgz, float *g11, float *g13, float *g33, int r)`.
/// Trielas/lib/findSiso.c; Andreas Rueger, 02/06/94.
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn findSiso(
    a: f32,
    pl: f32,
    gx: f32,
    gz: f32,
    px: &mut f32,
    pz: &mut f32,
    vgx: &mut f32,
    vgz: &mut f32,
    g11: &mut f32,
    g13: &mut f32,
    g33: &mut f32,
    r: i32,
) -> i32 {
    find_iso(a, pl, gx, gz, px, pz, vgx, vgz, g11, g13, g33, r, true)
}
#[allow(clippy::too_many_arguments)]
fn find_iso(
    a: f32,
    pl: f32,
    gx: f32,
    gz: f32,
    px: &mut f32,
    pz: &mut f32,
    vgx: &mut f32,
    vgz: &mut f32,
    g11: &mut f32,
    g13: &mut f32,
    g33: &mut f32,
    r: i32,
    shear: bool,
) -> i32 {
    let pm = if r == 1 {
        -(*pz * gx - gz * *px)
    } else if r == 0 {
        let pms = 1.0f32 / a - pl * pl;
        if pms < 0.0 {
            return -1;
        }
        (pms as f64).sqrt() as f32
    } else {
        if shear {
            oxwave_core::err!(" wrong mode in findSiso");
        } else {
            oxwave_core::err!(" wrong mode in findPiso");
        }
    };
    let x = pl * gx - pm * gz;
    let z = pm * gx + pl * gz;
    let vx = x * a;
    let vz = z * a;
    if shear {
        *g11 = vz * z;
        *g33 = vx * x;
        *g13 = -vx * z;
    } else {
        *g11 = vx * x;
        *g33 = vz * z;
        *g13 = vx * z;
    }
    *vgx = vx;
    *vgz = vz;
    *px = x;
    *pz = z;
    1
}

/// Continue SH slowness. C: `int findqSH(float si, float co, float pl, float a1212, float a1223, float a2323, float *pxnew, float *pznew, float *vgx, float *vgz, int r)`.
/// Trielas/lib/findqSH.c; Andreas Rueger, 02/01/94.
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn findqSH(
    si: f32,
    co: f32,
    pl: f32,
    a1212: f32,
    a1223: f32,
    a2323: f32,
    pxnew: &mut f32,
    pznew: &mut f32,
    vgx: &mut f32,
    vgz: &mut f32,
    r: i32,
) -> i32 {
    let n1 = -si;
    let n3 = co;
    let b1 = pl * co;
    let b3 = pl * si;
    let a = n1 * n1 * a1212 + n3 * n3 * a2323 + 2.0 * n1 * n3 * a1223;
    let b = 2.0 * (b1 * n1 * a1212 + b3 * n3 * a2323 + a1223 * (b1 * n3 + b3 * n1));
    let c = b1 * b1 * a1212 + b3 * b3 * a2323 + 2.0 * a1223 * b1 * b3 - 1.0;
    let sqr = b * b - 4.0 * a * c;
    if sqr < 0.0 || a.abs() < f32::EPSILON {
        return 0;
    }
    let sqr = (sqr as f64).sqrt() as f32;
    let eta = if r == 0 {
        (0.5 / a as f64 * (sqr - b) as f64) as f32
    } else if r == 1 {
        (0.5 / a as f64 * (-sqr - b) as f64) as f32
    } else {
        oxwave_core::err!(" ERROR in findqSH \n")
    };
    let px = b1 + eta * n1;
    let pz = b3 + eta * n3;
    let vx = a1212 * px + a1223 * pz;
    let vz = a1223 * px + a2323 * pz;
    let vgm = -si * vx + co * vz;
    if (vgm > 0.0 && r == 1) || (vgm < 0.0 && r == 0) {
        return 0;
    }
    *pxnew = px;
    *pznew = pz;
    *vgx = vx;
    *vgz = vz;
    1
}

/// SH reflection or transmission coefficient. C: `int rt_SHa_real(float c1212i, float c2323i, float c1223i, float c1212t, float c2323t, float c1223t, float pxi, float pzi, float pxt, float pzt, float pxr, float pzr, int rort, float *coeff, float gz, float gx)`.
/// Trielas/lib/rt_SHa_real.c; Andreas Rueger, 02/01/94.
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn rt_SHa_real(
    mut c1212i: f32,
    mut c2323i: f32,
    mut c1223i: f32,
    mut c1212t: f32,
    mut c2323t: f32,
    mut c1223t: f32,
    mut pxi: f32,
    mut pzi: f32,
    mut pxt: f32,
    mut pzt: f32,
    mut pxr: f32,
    mut pzr: f32,
    rort: i32,
    coeff: &mut f32,
    gz: f32,
    gx: f32,
) -> i32 {
    rotvector(&mut pxt, &mut pzt, gz, gx);
    rotvector(&mut pxi, &mut pzi, gz, gx);
    rotvector(&mut pxr, &mut pzr, gz, gx);
    rottensh(&mut c1212i, &mut c2323i, &mut c1223i, gz, gx);
    rottensh(&mut c1212t, &mut c2323t, &mut c1223t, gz, gx);
    if rort == 2 {
        let d = c1223i * pxr + c2323i * pzr;
        if d.abs() < f32::EPSILON {
            return -1;
        }
        *coeff = -(c1223i * pxi + c2323i * pzi) / d;
    } else {
        let d = c1223i * pxr + c2323i * pzr - c1223t * pxt - c2323t * pzt;
        if d.abs() < f32::EPSILON {
            return -1;
        }
        let rc = (c1223t * pxt + c2323t * pzt - c1223i * pxi - c2323i * pzi) / d;
        *coeff = if rort == 1 {
            rc
        } else if rort == 0 {
            rc + 1.0
        } else {
            oxwave_core::err!(" wrong rort in int rt_SHa_real")
        };
    }
    1
}

/// Continue qP slowness across an interface, Trielas/lib/solveForSlowqP.c.
/// Credits: Colorado School of Mines (the C file does not name an author).
/// C: `int solveForSlowqP(float s, float c, float pl, float a1111, float a3333, float a1133, float a1313, float a1113, float a3313, float *pxnew, float *pznew, float rt)`.
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn solveForSlowqP(
    s: f32,
    c: f32,
    pl: f32,
    a1111: f32,
    a3333: f32,
    a1133: f32,
    a1313: f32,
    a1113: f32,
    a3313: f32,
    pxnew: &mut f32,
    pznew: &mut f32,
    rt: f32,
) -> i32 {
    slow(
        s, c, pl, a1111, a3333, a1133, a1313, a1113, a3313, pxnew, pznew, rt, false,
    )
}
/// Continue qSV slowness across an interface, Trielas/lib/solveForSlowqSV.c.
/// Credits: Colorado School of Mines (the C file does not name an author).
/// C: `int solveForSlowqSV(float s, float c, float pl, float a1111, float a3333, float a1133, float a1313, float a1113, float a3313, float *pxnew, float *pznew, float rt)`.
#[allow(non_snake_case, clippy::too_many_arguments)]
pub fn solveForSlowqSV(
    s: f32,
    c: f32,
    pl: f32,
    a1111: f32,
    a3333: f32,
    a1133: f32,
    a1313: f32,
    a1113: f32,
    a3313: f32,
    pxnew: &mut f32,
    pznew: &mut f32,
    rt: f32,
) -> i32 {
    slow(
        s, c, pl, a1111, a3333, a1133, a1313, a1113, a3313, pxnew, pznew, rt, true,
    )
}
#[allow(clippy::too_many_arguments)]
fn slow(
    s: f32,
    c: f32,
    pl: f32,
    a: f32,
    b: f32,
    f: f32,
    l: f32,
    m: f32,
    n: f32,
    outx: &mut f32,
    outz: &mut f32,
    rt: f32,
    sv: bool,
) -> i32 {
    let vertical = ((s * s) as f64) < if sv { 0.5 } else { 0.49 };
    let base = if sv { l } else { b };
    let mut sq = (1.0f32 / base - pl * pl) as f64;
    if sq < 0.0 {
        sq = 0.0;
    }
    let (s, c, pl, rt) = (s as f64, c as f64, pl as f64, rt as f64);
    let (a, b, l, m, n) = (a as f64, b as f64, l as f64, m as f64, n as f64);
    let fl = (f + l as f32) as f64;
    // In C, 2*a1113 and 2*a3313 are evaluated as float before pxz
    // promotes the products to double.
    let twom = (2.0f32 * m as f32) as f64;
    let twon = (2.0f32 * n as f32) as f64;
    let eigen = |x: f64, z: f64, clamp: bool, zero: bool| {
        let x2 = x * x;
        let z2 = z * z;
        let xz = x * z;
        let g11 = a * x2 + l * z2 + twom * xz;
        let g33 = b * z2 + l * x2 + twon * xz;
        let g13 = fl * xz + m * x2 + n * z2;
        let mut d = (g11 + g33) * (g11 + g33) - 4.0 * (g11 * g33 - g13 * g13);
        if clamp && d < 0.0 {
            d = 0.0;
        }
        let root = if zero { 0.0 } else { d.sqrt() };
        0.5 * (g11 + g33 + if sv { -root } else { root }) - 1.0
    };
    let (mut x, mut z);
    let mut iter = 0;
    if vertical {
        let oc = 1.0 / c;
        z = (pl as f32 * s as f32) as f64 + rt * sq.sqrt() * c;
        x = (pl - z * s) * oc;
        let mut diff = eigen(x, z, true, false);
        let mut err = diff * diff;
        let mut prev = z;
        if err > 1e-11 {
            z = prev + 0.0001 / (base as f64).sqrt();
        }
        while err > 1e-11 && iter < 50 {
            x = (pl - z * s) * oc;
            let d = eigen(x, z, false, false);
            err = d * d;
            let next = (diff * z - d * prev) / (diff - d);
            diff = d;
            prev = z;
            z = next;
            iter += 1;
        }
    } else {
        let os = 1.0 / s;
        x = (pl as f32 * c as f32) as f64 - rt * sq.sqrt() * s;
        z = (pl - x * c) * os;
        // The qP C branch leaves sqr at its initial zero here; preserve it.
        let mut diff = eigen(x, z, true, !sv);
        let mut err = diff * diff;
        let mut prev = x;
        x = prev + 0.0001 / (base as f64).sqrt();
        while err > 1e-11 && iter < 50 {
            z = (pl - x * c) * os;
            let d = eigen(x, z, false, false);
            err = d * d;
            let next = (diff * x - d * prev) / (diff - d);
            diff = d;
            prev = x;
            x = next;
            iter += 1;
        }
    }
    *outx = x as f32;
    *outz = z as f32;
    if iter >= 49 {
        -1
    } else {
        1
    }
}

/// Isotropic real reflection and transmission coefficients.
/// Trielas/lib/rt_iso_real.c; Andreas Rueger, 02/01/94.
/// C: `int rt_iso_real(float vp1, float vp2, float vs1, float vs2, float rho1, float rho2, float pl, int modei, int modet, int rort, float *coeff)`.
#[allow(clippy::too_many_arguments, non_snake_case)]
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
    coeff: &mut f32,
) -> i32 {
    let si1 = pl * vp1;
    let si2 = pl * vp2;
    let sj1 = pl * vs1;
    let sj2 = pl * vs2;
    if si1 > 1.0 || si2 > 1.0 || sj1 > 1.0 || sj2 > 1.0 {
        return 0;
    }
    let cj1 = ((1.0 - sj1 * sj1) as f64).sqrt() as f32;
    let cj2 = ((1.0 - sj2 * sj2) as f64).sqrt() as f32;
    if modei == 2 && modet == 2 {
        let del = rho1 * vs1 * cj1 + rho2 * vs2 * cj2;
        if del.abs() < f32::EPSILON {
            return -1;
        }
        *coeff = if rort == 1 {
            1.0 / del * (rho1 * vs1 * cj1 - rho2 * vs2 * cj2)
        } else if rort == 0 {
            1.0 / del * 2.0 * rho1 * vs1 * cj1
        } else if rort == 2 {
            1.0
        } else {
            oxwave_core::err!("ERROR in reftransiso\n")
        };
        return 1;
    }
    let ci1 = ((1.0 - si1 * si1) as f64).sqrt() as f32;
    let ci2 = ((1.0 - si2 * si2) as f64).sqrt() as f32;
    let pl2 = pl * pl;
    let vs12 = vs1 * vs1;
    let vs22 = vs2 * vs2;
    let covi1 = ci1 / vp1;
    let covi2 = ci2 / vp2;
    let covj1 = cj1 / vs1;
    let covj2 = cj2 / vs2;
    let temp = rho2 * (1.0 - 2.0 * vs22 * pl2);
    let mut a = temp - rho1 * (1.0 - 2.0 * vs12 * pl2);
    let mut b = temp + 2.0 * rho1 * vs12 * pl2;
    let c = rho1 * (1.0 - 2.0 * vs12 * pl2) + 2.0 * rho2 * vs22 * pl2;
    let d = 2.0 * (rho2 * vs22 - rho1 * vs12);
    let e = b * covi1 + c * covi2;
    let f = b * covj1 + c * covj2;
    let g = a - d * covi1 * covj2;
    let h = a - d * covi2 * covj1;
    let dd = e * f + g * h * pl2;
    if dd.abs() < f32::EPSILON {
        return -1;
    }
    *coeff = match (modei, modet, rort) {
        (0, 0, 1) => ((b * covi1 - c * covi2) * f - (a + d * covi1 * covj2) * h * pl2) / dd,
        (0, 0, 0) => 2.0 * rho1 * ci1 * f / (vp2 * dd),
        (0, 1, 1) => -2.0 * covi1 * (a * b + c * d * covi2 * covj2) * pl * vp1 / (vs1 * dd),
        (0, 1, 0) => 2.0 * rho1 * ci1 * h * pl / (vs2 * dd),
        (1, 0, 1) => -2.0 * covj1 * (a * b + c * d * covi2 * covj2) * pl * vs1 / (vp1 * dd),
        (1, 0, 0) => -2.0 * rho1 * cj1 * g * pl / (vp2 * dd),
        (1, 1, 1) => -((b * covj1 - c * covj2) * e - (a + d * covi2 * covj1) * g * pl2) / dd,
        (1, 1, 0) => 2.0 * rho1 * cj1 * e / (vs2 * dd),
        (_, _, 2) => {
            a = 1.0 / vs12 - 2.0 * pl2;
            b = 4.0 * pl2 * ci1 * cj1;
            let t = a * a * vs1 * vp1 + b;
            if t.abs() < f32::EPSILON {
                return -1;
            }
            match (modei, modet) {
                (0, 0) => (-a * a * vs1 * vp1 + b) / t,
                (0, 1) => -4.0 * pl * vp1 * ci1 * a / t,
                (1, 1) => -(a * a * vs1 * vp1 - b) / t,
                (1, 0) => 4.0 * vs1 * cj1 * pl * a / t,
                _ => return -1,
            }
        }
        _ => return -1,
    };
    1
}

/// Isotropic complex reflection and transmission coefficients.
/// Trielas/lib/rt_iso_cmplx.c; Andreas Rueger, 02/01/94.
/// C: `int rt_iso_cmplx(float vp1, float vp2, float vs1, float vs2, float rho1, float rho2, float pl, int modei, int modet, int rort, float *coeff, float *phase)`.
#[allow(clippy::too_many_arguments, non_snake_case)]
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
    coeff: &mut f32,
    phase: &mut f32,
) -> i32 {
    use oxwave_cwp::complex::*;
    // The frozen core's cinv uses scaled division. SU complex.c instead
    // computes this unscaled reciprocal, including f32 overflow in the norm
    // and the double literal in 1.0/norm. Keep those semantics here.
    let reciprocal = |z: Complex| {
        let scale = (1.0f64 / ((z.r * z.r + z.i * z.i) as f64)) as f32;
        cmplx(z.r * scale, -z.i * scale)
    };
    // SU complex.c squares the float ratio before promotion by the double
    // literal 1.0. The frozen core squares in f64, changing the last bit of
    // some magnitudes. Preserve C rounding for both the pole guard and coeff.
    let rcabs = |z: Complex| {
        let x = z.r.abs();
        let y = z.i.abs();
        if x == 0.0 {
            y
        } else if y == 0.0 {
            x
        } else if x > y {
            let temp = y / x;
            (x as f64 * (1.0 + (temp * temp) as f64).sqrt()) as f32
        } else {
            let temp = x / y;
            (y as f64 * (1.0 + (temp * temp) as f64).sqrt()) as f32
        }
    };
    let angle = |s: f32| {
        if s > 1.0 {
            cmplx(0.0, -(((s * s - 1.0) as f64).sqrt() as f32))
        } else {
            cmplx(((1.0 - s * s) as f64).sqrt() as f32, 0.0)
        }
    };
    let ci1 = angle(pl * vp1);
    let ci2 = angle(pl * vp2);
    let cj1 = angle(pl * vs1);
    let cj2 = angle(pl * vs2);
    let result = if modei == 2 && modet == 2 {
        let t1 = crmul(cj1, rho1 * vs1);
        let t2 = crmul(cj2, rho2 * vs2);
        let del = reciprocal(cadd(t1, t2));
        if rort == 1 {
            cmul(del, csub(t1, t2))
        } else if rort == 0 {
            cmul(del, crmul(t1, 2.0))
        } else {
            oxwave_core::err!("ERROR in reftransiso\n")
        }
    } else {
        let pl2 = pl * pl;
        let vs12 = vs1 * vs1;
        let vs22 = vs2 * vs2;
        let vi1 = crmul(ci1, 1.0 / vp1);
        let vi2 = crmul(ci2, 1.0 / vp2);
        let vj1 = crmul(cj1, 1.0 / vs1);
        let vj2 = crmul(cj2, 1.0 / vs2);
        let temp = rho2 * (1.0 - 2.0 * vs22 * pl2);
        let a = temp - rho1 * (1.0 - 2.0 * vs12 * pl2);
        let b = temp + 2.0 * rho1 * vs12 * pl2;
        let c = rho1 * (1.0 - 2.0 * vs12 * pl2) + 2.0 * rho2 * vs22 * pl2;
        let d = 2.0 * (rho2 * vs22 - rho1 * vs12);
        let e = cadd(crmul(vi1, b), crmul(vi2, c));
        let f = cadd(crmul(vj1, b), crmul(vj2, c));
        let t = crmul(cmul(vi1, vj2), d);
        let g = cmplx(a - t.r, t.i);
        let t = crmul(cmul(vi2, vj1), d);
        let h = cmplx(a - t.r, t.i);
        let dd = cadd(cmul(e, f), crmul(cmul(g, h), pl2));
        if rcabs(dd) < f32::EPSILON {
            return -1;
        }
        let inv = reciprocal(dd);
        match (modei, modet, rort) {
            (0, 0, 1) => cmul(
                csub(
                    cmul(csub(crmul(vi1, b), crmul(vi2, c)), f),
                    cmul(crmul(h, pl2), cadd(cmplx(a, 0.0), crmul(cmul(vi1, vj2), d))),
                ),
                inv,
            ),
            (0, 0, 0) => cmul(crmul(cmul(f, ci1), 2.0 * rho1 / vp2), inv),
            (0, 1, 1) => cmul(
                cmul(
                    cadd(cmplx(a * b, 0.0), crmul(cmul(vi2, vj2), c * d)),
                    crmul(vi1, -2.0 * pl * vp1 / vs1),
                ),
                inv,
            ),
            (0, 1, 0) => cmul(crmul(cmul(ci1, h), 2.0 * rho1 * pl / vs2), reciprocal(dd)),
            (1, 0, 1) => cmul(
                cmul(
                    cadd(cmplx(a * b, 0.0), crmul(cmul(vi2, vj2), c * d)),
                    crmul(vj1, -2.0 * pl * vs1 / vp1),
                ),
                reciprocal(dd),
            ),
            (1, 0, 0) => cmul(crmul(cmul(cj1, g), -2.0 * rho1 * pl / vp2), reciprocal(dd)),
            (1, 1, 1) => cmul(
                cadd(
                    cmul(csub(crmul(vj2, c), crmul(vj1, b)), e),
                    cmul(crmul(g, pl2), cadd(cmplx(a, 0.0), crmul(cmul(vi2, vj1), d))),
                ),
                reciprocal(dd),
            ),
            (1, 1, 0) => cmul(crmul(cmul(e, cj1), 2.0 * rho1), reciprocal(crmul(dd, vs2))),
            _ => return -1,
        }
    };
    if result.i.abs() < f32::EPSILON {
        *coeff = result.r;
        *phase = 0.0;
    } else if result.r.abs() < f32::EPSILON && result.i > 0.0 {
        *coeff = rcabs(result);
        *phase = (std::f64::consts::PI / 2.0) as f32;
    } else if result.r.abs() < f32::EPSILON && result.i < 0.0 {
        *coeff = rcabs(result);
        *phase = (-std::f64::consts::PI / 2.0) as f32;
    } else if result.r > 0.0 && result.i > 0.0 {
        *coeff = rcabs(result);
        *phase = ((result.i / result.r) as f64).atan() as f32;
    } else if result.r < 0.0 {
        *coeff = rcabs(result);
        *phase = (std::f64::consts::PI + ((result.i / result.r) as f64).atan()) as f32;
    } else if result.r > 0.0 && result.i < 0.0 {
        *coeff = rcabs(result);
        *phase = (2.0 * std::f64::consts::PI + ((result.i / result.r) as f64).atan()) as f32;
    }
    1
}
