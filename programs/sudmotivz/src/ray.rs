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

//! Numerical helpers from su/main/dip_moveout/sudmotivz.c (Revision 1.4).
//! Credits: CWP: Tariq Alkhalifah, Craig Artley.

#![allow(non_snake_case, unused_parens, unused_assignments, unused_mut)]

use oxwave_core::prelude::*;

// Trial solution for constant velocity on the prestack migration ellipse.
pub(super) fn getX(tsg: f32, h: f32, p0: f32, v: f32, vars: &mut [f32]) {
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut a2: f32 = 0.0;
    let mut b2: f32 = 0.0;
    let mut h2: f32 = 0.0;
    let mut x0: f32 = 0.0;
    let mut ps: f32 = 0.0;
    let mut pg: f32 = 0.0;
    let mut tg: f32 = 0.0;
    let mut t0: f32 = 0.0;
    let mut sine: f32 = 0.0;
    let mut tangent: f32 = 0.0;
    h2 = (h * h);
    a2 = ((((((v * v) * tsg) * tsg) as f64) / 4f64) as f32);
    b2 = ((if (((a2 - h2) as f64) > 0f64) {
        ((a2 - h2) as f64)
    } else {
        0f64
    }) as f32);
    sine = (((0.5f64 * (p0 as f64)) * (v as f64)) as f32);
    tangent = (((sine as f64)
        / (if ((1f64 - ((sine * sine) as f64)) > 0f64) {
            (1f64 - ((sine * sine) as f64))
        } else {
            0f64
        })
        .sqrt()) as f32);
    x = (((tangent as f64) * ((a2 / (b2 + (tangent * tangent))) as f64).sqrt()) as f32);
    z = (((b2 as f64)
        * (if ((1f64 - (((x * x) / a2) as f64)) > 0f64) {
            (1f64 - (((x * x) / a2) as f64))
        } else {
            0f64
        }))
    .sqrt() as f32);
    x0 = ((h2 / a2) * x);
    tg = ((((((h - x) * (h - x)) + (z * z)) as f64).sqrt() / (v as f64)) as f32);
    t0 = (((((((x0 - x) * (x0 - x)) + (z * z)) as f64).sqrt() * 2f64) / (v as f64)) as f32);
    ps = ((x + h) / ((v * v) * (tsg - tg)));
    pg = ((x - h) / ((v * v) * tg));
    vars[(0i32) as usize] = x0;
    vars[(1i32) as usize] = ps;
    vars[(2i32) as usize] = pg;
    vars[(3i32) as usize] = tg;
    vars[(4i32) as usize] = t0;
}

// Assemble the five equations and the column-major Jacobian.
pub(super) fn getFJ(
    x: &[Vec<f32>],
    a: &[Vec<f32>],
    tau: &[Vec<f32>],
    s: &[Vec<f32>],
    nt: i32,
    np: i32,
    dt: f32,
    dp: f32,
    tsg: f32,
    h: f32,
    p0: f32,
    vars: &[f32],
    eqns: &mut [f32],
    jacob: &mut [Vec<f32>],
) {
    let mut i: i32 = 0;
    let mut j: i32 = 0;
    let mut neqns: i32 = 5i32;
    let mut x0: f32 = 0.0;
    let mut ps: f32 = 0.0;
    let mut pg: f32 = 0.0;
    let mut ts: f32 = 0.0;
    let mut tg: f32 = 0.0;
    let mut t0: f32 = 0.0;
    let mut ft: f32 = (0f64 as f32);
    let mut fp: f32 = (0f64 as f32);
    let mut ftp: f32 = 0.0;
    let mut dfdt: f32 = 0.0;
    let mut dfdp: f32 = 0.0;
    let mut sign: f32 = 0.0;
    let mut f0tp: f32 = 0.0;
    let mut df0dt: f32 = 0.0;
    let mut df0dp: f32 = 0.0;
    let mut stp: f32 = 0.0;
    let mut dsdt: f32 = 0.0;
    let mut dsdp: f32 = 0.0;
    let mut signp0: f32 = 0.0;
    x0 = vars[(0i32) as usize];
    ps = vars[(1i32) as usize];
    pg = vars[(2i32) as usize];
    tg = vars[(3i32) as usize];
    t0 = vars[(4i32) as usize];
    ts = ((if (((tsg - tg) as f64) > 0f64) {
        ((tsg - tg) as f64)
    } else {
        0f64
    }) as f32);
    i = 0i32;
    while (i < neqns) {
        eqns[(i) as usize] = (0f64 as f32);
        i += 1;
    }
    j = 0i32;
    while (j < neqns) {
        i = 0i32;
        while (i < neqns) {
            jacob[(j) as usize][(i) as usize] = (0f64 as f32);
            i += 1;
        }
        j += 1;
    }
    eqns[(0i32) as usize] = ((eqns[(0i32) as usize] as f64) + (2f64 * (h as f64))) as f32;
    eqns[(1i32) as usize] += (h - x0);
    jacob[(0i32) as usize][(1i32) as usize] =
        ((jacob[(0i32) as usize][(1i32) as usize] as f64) - 1f64) as f32;
    intder(
        x,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        ((2f64 * (ts as f64)) as f32),
        (if (ps < (0i32 as f32)) { (-ps) } else { ps }),
        &mut ftp,
        &mut dfdt,
        &mut dfdp,
    );
    sign = ((if (ps < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    eqns[(0i32) as usize] -= (sign * ftp);
    jacob[(1i32) as usize][(0i32) as usize] -= dfdp;
    jacob[(3i32) as usize][(0i32) as usize] = ((jacob[(3i32) as usize][(0i32) as usize] as f64)
        + ((2f64 * (sign as f64)) * (dfdt as f64)))
        as f32;
    intder(
        x,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        ((2f64 * (tg as f64)) as f32),
        (if (pg < (0i32 as f32)) { (-pg) } else { pg }),
        &mut ftp,
        &mut dfdt,
        &mut dfdp,
    );
    sign = ((if (pg < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    eqns[(0i32) as usize] += (sign * ftp);
    eqns[(1i32) as usize] += (sign * ftp);
    jacob[(2i32) as usize][(0i32) as usize] += dfdp;
    jacob[(2i32) as usize][(1i32) as usize] += dfdp;
    jacob[(3i32) as usize][(0i32) as usize] = ((jacob[(3i32) as usize][(0i32) as usize] as f64)
        + ((2f64 * (sign as f64)) * (dfdt as f64)))
        as f32;
    jacob[(3i32) as usize][(1i32) as usize] = ((jacob[(3i32) as usize][(1i32) as usize] as f64)
        + ((2f64 * (sign as f64)) * (dfdt as f64)))
        as f32;
    intder(
        x,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        t0,
        (if (p0 < (0i32 as f32)) { (-p0) } else { p0 }),
        &mut ftp,
        &mut dfdt,
        &mut dfdp,
    );
    sign = ((if (p0 < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    eqns[(1i32) as usize] -= (sign * ftp);
    jacob[(4i32) as usize][(1i32) as usize] -= (sign * dfdt);
    intder(
        tau,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        ((2f64 * (ts as f64)) as f32),
        (if (ps < (0i32 as f32)) { (-ps) } else { ps }),
        &mut ftp,
        &mut dfdt,
        &mut dfdp,
    );
    sign = ((if (ps < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    eqns[(2i32) as usize] -= ftp;
    jacob[(1i32) as usize][(2i32) as usize] -= (sign * dfdp);
    jacob[(3i32) as usize][(2i32) as usize] =
        ((jacob[(3i32) as usize][(2i32) as usize] as f64) + (2f64 * (dfdt as f64))) as f32;
    intder(
        tau,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        ((2f64 * (tg as f64)) as f32),
        (if (pg < (0i32 as f32)) { (-pg) } else { pg }),
        &mut ftp,
        &mut dfdt,
        &mut dfdp,
    );
    sign = ((if (pg < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    eqns[(2i32) as usize] += ftp;
    eqns[(3i32) as usize] += ftp;
    jacob[(2i32) as usize][(2i32) as usize] += (sign * dfdp);
    jacob[(2i32) as usize][(3i32) as usize] += (sign * dfdp);
    jacob[(3i32) as usize][(2i32) as usize] =
        ((jacob[(3i32) as usize][(2i32) as usize] as f64) + (2f64 * (dfdt as f64))) as f32;
    jacob[(3i32) as usize][(3i32) as usize] =
        ((jacob[(3i32) as usize][(3i32) as usize] as f64) + (2f64 * (dfdt as f64))) as f32;
    intder(
        tau,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        t0,
        (if (p0 < (0i32 as f32)) { (-p0) } else { p0 }),
        &mut ftp,
        &mut dfdt,
        &mut dfdp,
    );
    sign = ((if (p0 < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    eqns[(3i32) as usize] -= ftp;
    jacob[(4i32) as usize][(3i32) as usize] -= dfdt;
    intder(
        a,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        ((2f64 * (ts as f64)) as f32),
        (if (ps < (0i32 as f32)) { (-ps) } else { ps }),
        &mut ftp,
        &mut dfdt,
        &mut dfdp,
    );
    intder(
        a,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        t0,
        (if (p0 < (0i32 as f32)) { (-p0) } else { p0 }),
        &mut f0tp,
        &mut df0dt,
        &mut df0dp,
    );
    intder(
        s,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        ((2f64 * (ts as f64)) as f32),
        (if (ps < (0i32 as f32)) { (-ps) } else { ps }),
        &mut stp,
        &mut dsdt,
        &mut dsdp,
    );
    sign = ((if (ps < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    signp0 = ((if (p0 < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    eqns[(4i32) as usize] = ((eqns[(4i32) as usize] as f64)
        - ((stp as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).sin()))
        as f32;
    jacob[(1i32) as usize][(4i32) as usize] = ((jacob[(1i32) as usize][(4i32) as usize] as f64)
        + (((((-sign) * dsdp) as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).sin())
            - (((stp as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).cos()) * (dfdp as f64))))
        as f32;
    jacob[(3i32) as usize][(4i32) as usize] = ((jacob[(3i32) as usize][(4i32) as usize] as f64)
        + (2f64
            * (((dsdt as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).sin())
                + ((((sign * stp) as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).cos())
                    * (dfdt as f64))))) as f32;
    jacob[(4i32) as usize][(4i32) as usize] = ((jacob[(4i32) as usize][(4i32) as usize] as f64)
        + ((((2f64 * (signp0 as f64)) * (stp as f64))
            * (((sign * ftp) - (signp0 * f0tp)) as f64).cos())
            * (df0dt as f64))) as f32;
    intder(
        a,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        ((2f64 * (tg as f64)) as f32),
        (if (pg < (0i32 as f32)) { (-pg) } else { pg }),
        &mut ftp,
        &mut dfdt,
        &mut dfdp,
    );
    intder(
        s,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        ((2f64 * (tg as f64)) as f32),
        (if (pg < (0i32 as f32)) { (-pg) } else { pg }),
        &mut stp,
        &mut dsdt,
        &mut dsdp,
    );
    sign = ((if (pg < (0i32 as f32)) { (-1f64) } else { 1f64 }) as f32);
    eqns[(4i32) as usize] = ((eqns[(4i32) as usize] as f64)
        - ((stp as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).sin()))
        as f32;
    jacob[(2i32) as usize][(4i32) as usize] = ((jacob[(2i32) as usize][(4i32) as usize] as f64)
        + (((((-sign) * dsdp) as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).sin())
            - (((stp as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).cos()) * (dfdp as f64))))
        as f32;
    jacob[(3i32) as usize][(4i32) as usize] = ((jacob[(3i32) as usize][(4i32) as usize] as f64)
        - (2f64
            * (((dsdt as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).sin())
                + ((((sign * stp) as f64) * (((sign * ftp) - (signp0 * f0tp)) as f64).cos())
                    * (dfdt as f64))))) as f32;
    jacob[(4i32) as usize][(4i32) as usize] = ((jacob[(4i32) as usize][(4i32) as usize] as f64)
        + ((((2f64 * (signp0 as f64)) * (stp as f64))
            * (((sign * ftp) - (signp0 * f0tp)) as f64).cos())
            * (df0dt as f64))) as f32;
}

// Bilinear interpolation and finite differences; constant extrapolation
// outside the mesh gives zero derivatives in the extrapolated direction.
pub(super) fn intder(
    f: &[Vec<f32>],
    nx: i32,
    dx: f32,
    fx: f32,
    ny: i32,
    dy: f32,
    fy: f32,
    x: f32,
    y: f32,
    fxy: &mut f32,
    dfdx: &mut f32,
    dfdy: &mut f32,
) {
    let mut ix0: i32 = 0;
    let mut ix1: i32 = 0;
    let mut iy0: i32 = 0;
    let mut iy1: i32 = 0;
    let mut x0i: f32 = 0.0;
    let mut y0i: f32 = 0.0;
    let mut fracx: f32 = 0.0;
    let mut fracy: f32 = 0.0;
    let mut lx: f32 = 0.0;
    let mut ly: f32 = 0.0;
    let mut f00: f32 = 0.0;
    let mut f10: f32 = 0.0;
    let mut f01: f32 = 0.0;
    let mut f11: f32 = 0.0;
    lx = (fx + (((nx - 1i32) as f32) * dx));
    ly = (fy + (((ny - 1i32) as f32) * dy));
    if (x <= fx) {
        ix1 = 0i32;
        ix0 = ix1;
        fracx = (0f64 as f32);
    } else {
        if (x >= lx) {
            ix1 = (nx - 1i32);
            ix0 = ix1;
            fracx = (1f64 as f32);
        } else {
            x0i = ((x - fx) / dx);
            ix0 = (x0i as i32);
            ix1 = (ix0 + 1i32);
            fracx = (x0i - (ix0 as f32));
        }
    }
    if (y <= fy) {
        iy1 = 0i32;
        iy0 = iy1;
        fracy = (0i32 as f32);
    } else {
        if (y >= ly) {
            iy1 = (ny - 1i32);
            iy0 = iy1;
            fracy = (1f64 as f32);
        } else {
            y0i = ((y - fy) / dy);
            iy0 = (y0i as i32);
            iy1 = (iy0 + 1i32);
            fracy = (y0i - (iy0 as f32));
        }
    }
    f00 = f[(iy0) as usize][(ix0) as usize];
    f01 = f[(iy1) as usize][(ix0) as usize];
    f10 = f[(iy0) as usize][(ix1) as usize];
    f11 = f[(iy1) as usize][(ix1) as usize];
    (*fxy) = (((((((1i32 as f32) - fracx) * ((1i32 as f32) - fracy)) * f00)
        + ((((1i32 as f32) - fracx) * fracy) * f01))
        + ((fracx * ((1i32 as f32) - fracy)) * f10))
        + ((fracx * fracy) * f11));
    (*dfdx) = (((((1i32 as f32) - fracy) * (f10 - f00)) + (fracy * (f11 - f01))) / dx);
    (*dfdy) = (((((1i32 as f32) - fracx) * (f01 - f00)) + (fracx * (f11 - f10))) / dy);
}

// Trace a ray through the elastic coefficients as a function of vertical
// time. After returning to the surface, extrapolate the surface values.
pub(super) fn rayvt(
    p: f32,
    nt: i32,
    dt: f32,
    a1111: &[f32],
    a3333: &[f32],
    a1313: &[f32],
    a1133: &[f32],
    tau: &mut [f32],
    x: &mut [f32],
    a: &mut [f32],
    s: &mut [f32],
) {
    let mut it: i32 = 0;
    let mut jt: i32 = 0;
    let mut itau: i32 = 0;
    let mut turn: i32 = 0;
    let mut at: f32 = 0.0;
    let mut taut: f32 = 0.0;
    let mut xt: f32 = 0.0;
    let mut taui: f32 = 0.0;
    let mut frac: f32 = 0.0;
    let mut a1111t: f32 = 0.0;
    let mut a3333t: f32 = 0.0;
    let mut a1313t: f32 = 0.0;
    let mut a1133t: f32 = 0.0;
    let mut f1: f32 = 0.0;
    let mut f2: f32 = 0.0;
    let mut f3: f32 = 0.0;
    let mut f4: f32 = 0.0;
    let mut f5: f32 = 0.0;
    let mut f6: f32 = 0.0;
    let mut eps: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut pxz: f32 = 0.0;
    let mut px2: f32 = 0.0;
    let mut pz2: f32 = 0.0;
    let mut alpha: f32 = 0.0;
    let mut beta: f32 = 0.0;
    let mut gamma: f32 = 0.0;
    let mut det: f32 = 0.0;
    let mut rad: f32 = 0.0;
    let mut signbeta: f32 = 0.0;
    let mut q: f32 = 0.0;
    let mut vp: f32 = 0.0;
    let mut gamma11: f32 = 0.0;
    let mut gamma33: f32 = 0.0;
    let mut gamma13: f32 = 0.0;
    let mut den: f32 = 0.0;
    let mut g11: f32 = 0.0;
    let mut g33: f32 = 0.0;
    let mut g13: f32 = 0.0;
    let mut fac: f32 = 0.0;
    let mut dx: f32 = 0.0;
    let mut dz: f32 = 0.0;
    let mut dtau: f32 = 0.0;
    let mut sp: f32 = 0.0;
    a1111t = a1111[(0i32) as usize];
    a3333t = a3333[(0i32) as usize];
    a1313t = a1313[(0i32) as usize];
    a1133t = a1133[(0i32) as usize];
    f1 = (a1111t + a1313t);
    f2 = (a3333t + a1313t);
    f3 = (a1111t - a1313t);
    f4 = (a3333t - a1313t);
    f5 = (((2f64 * ((a1133t + a1313t) as f64)) * ((a1133t + a1313t) as f64)) as f32);
    f6 = (2i32 as f32);
    eps = (1.0000000000000001E-5f64 as f32);
    px2 = (p * p);
    alpha = ((f2 * f2) - (f4 * f4));
    beta = ((2i32 as f32) * (((((f1 * f2) + (f3 * f4)) - f5) * px2) - (f2 * f6)));
    gamma = ((((f6 * f6) as f64)
        - ((((2f64 * (f1 as f64)) * (f6 as f64)) - ((((f1 * f1) - (f3 * f3)) * px2) as f64))
            * (px2 as f64))) as f32);
    det = ((((beta * beta) as f64) - ((4f64 * (alpha as f64)) * (gamma as f64))) as f32);
    if (det < (0i32 as f32)) {
        return;
    }
    rad = ((det as f64).sqrt() as f32);
    if ((if (beta < (0i32 as f32)) {
        (-beta)
    } else {
        beta
    }) > eps)
    {
        signbeta = ((if (beta < (0i32 as f32)) {
            (-beta)
        } else {
            beta
        }) / beta);
    } else {
        signbeta = (1f64 as f32);
    }
    q = (((-0.5f64) * ((beta + (signbeta * rad)) as f64)) as f32);
    pz2 = (gamma / q);
    if (pz2 < (0i32 as f32)) {
        return;
    }
    sp = (((px2 + pz2) as f64).sqrt() as f32);
    vp = ((1i32 as f32) / sp);
    pz = ((pz2 as f64).sqrt() as f32);
    pxz = (p * pz);
    gamma11 = ((a1111t * px2) + (a1313t * pz2));
    gamma33 = ((a3333t * pz2) + (a1313t * px2));
    gamma13 = ((a1133t + a1313t) * pxz);
    den = ((1i32 as f32) / ((gamma11 + gamma33) - (2i32 as f32)));
    g11 = ((gamma33 - (1i32 as f32)) * den);
    g33 = ((gamma11 - (1i32 as f32)) * den);
    g13 = ((-gamma13) * den);
    fac = (((1i32 as f64) / (a3333t as f64).sqrt()) as f32);
    dx = ((((a1111t * p) * g11) + (((a1133t + a1313t) * pz) * g13)) + ((a1313t * g33) * p));
    dz = ((((a3333t * pz) * g33) + (((a1133t + a1313t) * p) * g13)) + ((a1313t * g11) * pz));
    dtau = ((if (1f64 < ((dz * fac) as f64)) {
        1f64
    } else {
        ((dz * fac) as f64)
    }) as f32);
    at = ((if (1f64 < ((pz * vp) as f64)) {
        1f64
    } else {
        ((pz * vp) as f64)
    })
    .acos() as f32);
    xt = (0f64 as f32);
    taut = xt;
    tau[(0i32) as usize] = taut;
    x[(0i32) as usize] = xt;
    a[(0i32) as usize] = at;
    s[(0i32) as usize] = sp;
    turn = 1i32;
    it = 1i32;
    while (it < nt) {
        tau[(it) as usize] = (tau[(it - 1i32) as usize] + (dt * dtau));
        if ((tau[(it) as usize] as f64) <= 0f64) {
            break;
        }
        x[(it) as usize] = (x[(it - 1i32) as usize] + (dx * dt));
        if (it < (nt - 1i32)) {
            taui = (tau[(it) as usize] / dt);
            itau = (taui as i32);
            frac = (taui - (itau as f32));
            a1111t = ((((1f64 - (frac as f64)) * (a1111[(itau) as usize] as f64))
                + ((frac * a1111[(itau + 1i32) as usize]) as f64)) as f32);
            a3333t = ((((1f64 - (frac as f64)) * (a3333[(itau) as usize] as f64))
                + ((frac * a3333[(itau + 1i32) as usize]) as f64)) as f32);
            a1313t = ((((1f64 - (frac as f64)) * (a1313[(itau) as usize] as f64))
                + ((frac * a1313[(itau + 1i32) as usize]) as f64)) as f32);
            a1133t = ((((1f64 - (frac as f64)) * (a1133[(itau) as usize] as f64))
                + ((frac * a1133[(itau + 1i32) as usize]) as f64)) as f32);
        }
        f1 = (a1111t + a1313t);
        f2 = (a3333t + a1313t);
        f3 = (a1111t - a1313t);
        f4 = (a3333t - a1313t);
        f5 = (((2f64 * ((a1133t + a1313t) as f64)) * ((a1133t + a1313t) as f64)) as f32);
        alpha = ((f2 * f2) - (f4 * f4));
        beta = ((2i32 as f32) * (((((f1 * f2) + (f3 * f4)) - f5) * px2) - (f2 * f6)));
        gamma = ((((f6 * f6) as f64)
            - ((((2f64 * (f1 as f64)) * (f6 as f64)) - ((((f1 * f1) - (f3 * f3)) * px2) as f64))
                * (px2 as f64))) as f32);
        det = ((((beta * beta) as f64) - ((4f64 * (alpha as f64)) * (gamma as f64))) as f32);
        if (det < (0i32 as f32)) {
            fprintf!(stderr(), "det=%f\n", (det as f64));
            return;
        }
        rad = ((det as f64).sqrt() as f32);
        if ((if (beta < (0i32 as f32)) {
            (-beta)
        } else {
            beta
        }) > eps)
        {
            signbeta = ((if (beta < (0i32 as f32)) {
                (-beta)
            } else {
                beta
            }) / beta);
        } else {
            signbeta = (1f64 as f32);
        }
        q = (((-0.5f64) * ((beta + (signbeta * rad)) as f64)) as f32);
        pz2 = (gamma / q);
        if (pz2 < (0i32 as f32)) {
            turn *= (-1i32);
        }
        pz2 = (if (pz2 < (0i32 as f32)) { (-pz2) } else { pz2 });
        pz = (((turn as f64) * (pz2 as f64).sqrt()) as f32);
        sp = (((px2 + pz2) as f64).sqrt() as f32);
        vp = ((1i32 as f32) / sp);
        pxz = (p * pz);
        gamma11 = ((a1111t * px2) + (a1313t * pz2));
        gamma33 = ((a3333t * pz2) + (a1313t * px2));
        gamma13 = ((a1133t + a1313t) * pxz);
        den = ((1i32 as f32) / ((gamma11 + gamma33) - (2i32 as f32)));
        g11 = ((gamma33 - (1i32 as f32)) * den);
        g33 = ((gamma11 - (1i32 as f32)) * den);
        g13 = ((-gamma13) * den);
        fac = (((1i32 as f64) / (a3333t as f64).sqrt()) as f32);
        dx = ((((a1111t * p) * g11) + (((a1133t + a1313t) * pz) * g13)) + ((a1313t * g33) * p));
        dz = ((((a3333t * pz) * g33) + (((a1133t + a1313t) * p) * g13)) + ((a1313t * g11) * pz));
        dtau = (dz * fac);
        at = ((if (1f64 < ((pz * vp) as f64)) {
            1f64
        } else {
            ((pz * vp) as f64)
        })
        .acos() as f32);
        a[(it) as usize] = at;
        s[(it) as usize] = sp;
        it += 1;
    }
    jt = it;
    while (jt < nt) {
        tau[(jt) as usize] = tau[(jt - 1i32) as usize];
        x[(jt) as usize] = x[(jt - 1i32) as usize];
        a[(jt) as usize] = a[(jt - 1i32) as usize];
        s[(jt) as usize] = s[(jt - 1i32) as usize];
        jt += 1;
    }
}
