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

//! Anisotropic Gaussian beam migration. Port of su/main/migration_inversion/sumiggbzoan.c (Revision 1.5).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: CWP: Tariq Alkhalifah, based on MIGGBZO by Dave Hale.
//! CWP: repackaged as an SU program by John Stockwell, April 2006.

#![allow(
    unused_mut,
    unused_assignments,
    unused_variables,
    non_snake_case,
    unused_parens
)]
use super::*;
pub(super) fn miggbzoaniso(
    bwh: f32,
    fmin: f32,
    fmax: f32,
    mut amin: f32,
    mut amax: f32,
    nt: i32,
    dt: f32,
    nx: i32,
    dx: f32,
    nz: i32,
    dz: f32,
    a1111: &[Vec<f32>],
    a3333: &[Vec<f32>],
    a1133: &[Vec<f32>],
    a1313: &[Vec<f32>],
    a1113: &[Vec<f32>],
    a3313: &[Vec<f32>],
    f: &[Vec<f32>],
    g: &mut [Vec<f32>],
    verbose: i32,
) {
    let mut nxb: i32 = 0;
    let mut npx: i32 = 0;
    let mut ntau: i32 = 0;
    let mut ipx: i32 = 0;
    let mut ix: i32 = 0;
    let mut ixb: i32 = 0;
    let mut ixlo: i32 = 0;
    let mut ixhi: i32 = 0;
    let mut nxw: i32 = 0;
    let mut iz: i32 = 0;
    let mut ixm: i32 = 0;
    let mut i: i32 = 0;
    let mut ft: f32 = 0.0;
    let mut fx: f32 = 0.0;
    let mut fz: f32 = 0.0;
    let mut xwh: f32 = 0.0;
    let mut dxb: f32 = 0.0;
    let mut fxb: f32 = 0.0;
    let mut xb: f32 = 0.0;
    let mut vmin: f32 = 0.0;
    let mut dpx: f32 = 0.0;
    let mut fpx: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut taupad: f32 = 0.0;
    let mut dtau: f32 = 0.0;
    let mut ftau: f32 = 0.0;
    let mut fxw: f32 = 0.0;
    let mut pxmin: f32 = 0.0;
    let mut pxmax: f32 = 0.0;
    let mut a0: f32 = 0.0;
    let mut x0: f32 = 0.0;
    let mut z0: f32 = 0.0;
    let mut bwhc: f32 = 0.0;
    let mut f1: f32 = 0.0;
    let mut f2: f32 = 0.0;
    let mut f3: f32 = 0.0;
    let mut f4: f32 = 0.0;
    let mut f5: f32 = 0.0;
    let mut f6: f32 = 0.0;
    let mut alpha: f32 = 0.0;
    let mut beta: f32 = 0.0;
    let mut gamma: f32 = 0.0;
    let mut det: f32 = 0.0;
    let mut rad: f32 = 0.0;
    let mut signbeta: f32 = 0.0;
    let mut eps: f32 = 0.0;
    let mut pz2: f32 = 0.0;
    let mut q: f32 = 0.0;
    let mut vp: f32 = 0.0;
    let mut px2: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut aa: f32 = 0.0;
    let mut bb: f32 = 0.0;
    let mut ccc: f32 = 0.0;
    let mut dd: f32 = 0.0;
    let mut ee: f32 = 0.0;
    let mut fxx: f32 = 0.0;
    let mut dfx: f32 = 0.0;
    let mut f7: f32 = 0.0;
    let mut f8: f32 = 0.0;
    let mut f9: f32 = 0.0;
    let mut f10: f32 = 0.0;
    let mut f11: f32 = 0.0;
    let mut dpz: f32 = 0.0;
    let mut gamm11: f32 = 0.0;
    let mut gamm33: f32 = 0.0;
    let mut gamm13: f32 = 0.0;
    let mut vpmin: f32 = 0.0;
    let mut vpmax: f32 = 0.0;
    let mut vp2: f32 = 0.0;
    let mut s: f32 = 0.0;
    let mut c: f32 = 0.0;
    let mut ss: f32 = 0.0;
    let mut cc: f32 = 0.0;
    /* first t, x, and z assumed to be zero */
    fz = (0_f64 as f32);
    fx = fz;
    ft = fx;
    /* convert minimum and maximum angles to radians */
    amin = (((amin as f64) * ((3.1415926535897931_f64) / 180_f64)) as f32);
    amax = (((amax as f64) * ((3.1415926535897931_f64) / 180_f64)) as f32);
    if (amin > amax) {
        let mut atemp: f32 = amin;
        amin = amax;
        amax = atemp;
    }
    /* window half-width */
    xwh = ((3_f64 * (bwh as f64)) as f32);
    /* beam center sampling */
    dxb = ((((if ((((bwh as f64) * ((2_f64 * (fmin as f64)) / (fmax as f64)).sqrt()) / (dx as f64))
        > 0_f64)
    {
        ((((bwh as f64) * ((2_f64 * (fmin as f64)) / (fmax as f64)).sqrt()) / (dx as f64))
            + 0.5_f64)
    } else {
        ((((bwh as f64) * ((2_f64 * (fmin as f64)) / (fmax as f64)).sqrt()) / (dx as f64))
            - 0.5_f64)
    }) as i32) as f32)
        * dx);
    nxb = (((1i32 as f32) + ((((nx - 1i32) as f32) * dx) / dxb)) as i32);
    fxb = (((fx as f64)
        + (0.5_f64 * (((((nx - 1i32) as f32) * dx) - (((nxb - 1i32) as f32) * dxb)) as f64)))
        as f32);
    if (verbose != 0) {
        warn!(
            "nxb=%d dxb=%g fxb=%g bwh=%g",
            nxb,
            (dxb as f64),
            (fxb as f64),
            (bwh as f64)
        );
    }
    /* minimum velocity at surface z=0 */
    ixm = 0i32;
    {
        {
            {
                ix = 1i32;
                ix
            };
            {
                vmin = a3333[(0i32) as usize][(0i32) as usize];
                vmin
            }
        };
        while (ix < nx) {
            if (a3333[(ix) as usize][(0i32) as usize] < vmin) {
                ixm = ix;
                vmin = a3333[(ix) as usize][(0i32) as usize];
            }
            ix += 1;
        }
    }
    /* phase velocity at min. and max. angles */
    c = ((amin as f64).cos() as f32);
    s = ((amin as f64).sin() as f32);
    cc = (c * c);
    ss = (s * s);
    gamm11 = (((a1111[(ixm) as usize][(0i32) as usize] * ss)
        + (a1313[(ixm) as usize][(0i32) as usize] * cc))
        + ((((2i32 as f32) * a1113[(ixm) as usize][(0i32) as usize]) * s) * c));
    gamm33 = (((a3333[(ixm) as usize][(0i32) as usize] * cc)
        + (a1313[(ixm) as usize][(0i32) as usize] * ss))
        + ((((2i32 as f32) * a3313[(ixm) as usize][(0i32) as usize]) * s) * c));
    gamm13 = (((((a1133[(ixm) as usize][(0i32) as usize]
        + a1313[(ixm) as usize][(0i32) as usize])
        * s)
        * c)
        + (a1113[(ixm) as usize][(0i32) as usize] * ss))
        + (a3313[(ixm) as usize][(0i32) as usize] * cc));
    vp2 = ((((gamm11 + gamm33) as f64)
        + ((((gamm11 + gamm33) * (gamm11 + gamm33))
            - ((4i32 as f32) * ((gamm11 * gamm33) - (gamm13 * gamm13)))) as f64)
            .sqrt()) as f32);
    vpmin = (((vp2 as f64) * 0.5_f64).sqrt() as f32);
    c = ((amax as f64).cos() as f32);
    s = ((amax as f64).sin() as f32);
    cc = (c * c);
    ss = (s * s);
    gamm11 = (((a1111[(ixm) as usize][(0i32) as usize] * ss)
        + (a1313[(ixm) as usize][(0i32) as usize] * cc))
        + ((((2i32 as f32) * a1113[(ixm) as usize][(0i32) as usize]) * s) * c));
    gamm33 = (((a3333[(ixm) as usize][(0i32) as usize] * cc)
        + (a1313[(ixm) as usize][(0i32) as usize] * ss))
        + ((((2i32 as f32) * a3313[(ixm) as usize][(0i32) as usize]) * s) * c));
    gamm13 = (((((a1133[(ixm) as usize][(0i32) as usize]
        + a1313[(ixm) as usize][(0i32) as usize])
        * s)
        * c)
        + (a1113[(ixm) as usize][(0i32) as usize] * ss))
        + (a3313[(ixm) as usize][(0i32) as usize] * cc));
    vp2 = ((((gamm11 + gamm33) as f64)
        + ((((gamm11 + gamm33) * (gamm11 + gamm33))
            - ((4i32 as f32) * ((gamm11 * gamm33) - (gamm13 * gamm13)))) as f64)
            .sqrt()) as f32);
    vpmax = (((vp2 as f64) * 0.5_f64).sqrt() as f32);
    /* beam sampling */
    pxmin = (((amin as f64).sin() / (vpmin as f64)) as f32);
    pxmax = (((amax as f64).sin() / (vpmax as f64)) as f32);
    dpx = ((((vmin as f64).sqrt() * 1_f64)
        / ((((2_f64 * (xwh as f64)) * ((fmin * fmax) as f64).sqrt()) * 0.5_f64)
            * ((vpmin + vpmax) as f64))) as f32);
    npx = (((1i32 as f32) + ((pxmax - pxmin) / dpx)) as i32);
    fpx = (((pxmin as f64) + (0.5_f64 * (((pxmax - pxmin) - (((npx - 1i32) as f32) * dpx)) as f64)))
        as f32);
    taupad = ((if ((if ((pxmin) < (0i32 as f32)) {
        (-(pxmin))
    } else {
        (pxmin)
    }) > (if ((pxmax) < (0i32 as f32)) {
        (-(pxmax))
    } else {
        (pxmax)
    })) {
        (if ((pxmin) < (0i32 as f32)) {
            (-(pxmin))
        } else {
            (pxmin)
        })
    } else {
        (if ((pxmax) < (0i32 as f32)) {
            (-(pxmax))
        } else {
            (pxmax)
        })
    }) * xwh);
    taupad = ((((if (((taupad / dt) as f64) > 0_f64) {
        (((taupad / dt) as f64) + 0.5_f64)
    } else {
        (((taupad / dt) as f64) - 0.5_f64)
    }) as i32) as f32)
        * dt);
    ftau = (ft - taupad);
    dtau = dt;
    ntau = (((nt as f64) + ((2_f64 * (taupad as f64)) / (dtau as f64))) as i32);
    fprintf!(
        stderr(),
        "npx=%d dpx=%g fpx=%g\n",
        npx,
        (dpx as f64),
        (fpx as f64)
    );
    /* zero migrated image */
    {
        ix = 0i32;
        while (ix < nx) {
            {
                iz = 0i32;
                while (iz < nz) {
                    g[(ix) as usize][(iz) as usize] = (0_f64 as f32);
                    iz += 1;
                }
            }
            ix += 1;
        }
    }
    /* loop over beam centers */
    {
        {
            {
                ixb = 0i32;
                ixb
            };
            {
                xb = fxb;
                xb
            }
        };
        while (ixb < nxb) {
            /* horizontal window */
            ix = ((if ((((xb - fx) / dx) as f64) > 0_f64) {
                ((((xb - fx) / dx) as f64) + 0.5_f64)
            } else {
                ((((xb - fx) / dx) as f64) - 0.5_f64)
            }) as i32);
            ixlo = (if ((ix
                - ((if (((xwh / dx) as f64) > 0_f64) {
                    (((xwh / dx) as f64) + 0.5_f64)
                } else {
                    (((xwh / dx) as f64) - 0.5_f64)
                }) as i32))
                > (0i32))
            {
                (ix - ((if (((xwh / dx) as f64) > 0_f64) {
                    (((xwh / dx) as f64) + 0.5_f64)
                } else {
                    (((xwh / dx) as f64) - 0.5_f64)
                }) as i32))
            } else {
                (0i32)
            });
            ixhi = (if ((ix
                + ((if (((xwh / dx) as f64) > 0_f64) {
                    (((xwh / dx) as f64) + 0.5_f64)
                } else {
                    (((xwh / dx) as f64) - 0.5_f64)
                }) as i32))
                < (nx - 1i32))
            {
                (ix + ((if (((xwh / dx) as f64) > 0_f64) {
                    (((xwh / dx) as f64) + 0.5_f64)
                } else {
                    (((xwh / dx) as f64) - 0.5_f64)
                }) as i32))
            } else {
                (nx - 1i32)
            });
            nxw = ((1i32 + ixhi) - ixlo);
            fxw = (fx + (((ixlo - ix) as f32) * dx));
            if (verbose != 0) {
                warn!("ixb/nxb = %d/%d  ix = %d\n", ixb, nxb, ix);
            }
            /* allocate space for beams */
            let mut b = vec![vec![0.0f32; ntau as usize]; npx as usize];
            /* form beams at surface */
            formBeams(
                bwh,
                dxb,
                fmin,
                nt,
                dt,
                ft,
                nxw,
                dx,
                fxw,
                &f[ixlo as usize..],
                ntau,
                dtau,
                ftau,
                npx,
                dpx,
                fpx,
                &mut b,
            );
            /* loop over beams */
            {
                {
                    {
                        ipx = 0i32;
                        ipx
                    };
                    {
                        px = fpx;
                        px
                    }
                };
                while (ipx < npx) {
                    /*calculating emergence angle from ray parameter*/
                    f1 = (a1111[(ix) as usize][(0i32) as usize]
                        + a1313[(ix) as usize][(0i32) as usize]);
                    f2 = (a3333[(ix) as usize][(0i32) as usize]
                        + a1313[(ix) as usize][(0i32) as usize]);
                    f3 = (a1111[(ix) as usize][(0i32) as usize]
                        - a1313[(ix) as usize][(0i32) as usize]);
                    f4 = (a3333[(ix) as usize][(0i32) as usize]
                        - a1313[(ix) as usize][(0i32) as usize]);
                    f5 = (((2_f64
                        * ((a1133[(ix) as usize][(0i32) as usize]
                            + a1313[(ix) as usize][(0i32) as usize])
                            as f64))
                        * ((a1133[(ix) as usize][(0i32) as usize]
                            + a1313[(ix) as usize][(0i32) as usize])
                            as f64)) as f32);
                    f6 = (2i32 as f32);
                    eps = (1.0E-4_f64 as f32);
                    px2 = (px * px);
                    alpha = ((f2 * f2) - (f4 * f4));
                    beta = ((2i32 as f32) * (((((f1 * f2) + (f3 * f4)) - f5) * px2) - (f2 * f6)));
                    gamma = ((((f6 * f6) as f64)
                        - ((((2_f64 * (f1 as f64)) * (f6 as f64))
                            - ((((f1 * f1) - (f3 * f3)) * px2) as f64))
                            * (px2 as f64))) as f32);
                    det = ((((beta * beta) as f64) - ((4_f64 * (alpha as f64)) * (gamma as f64)))
                        as f32);
                    if (det < (0i32 as f32)) {
                        {
                            ipx += 1;
                            px = (px + dpx)
                        };
                        continue;
                    }
                    rad = ((det as f64).sqrt() as f32);
                    if ((if ((beta) < (0i32 as f32)) {
                        (-(beta))
                    } else {
                        (beta)
                    }) > eps)
                    {
                        signbeta = ((if ((beta) < (0i32 as f32)) {
                            (-(beta))
                        } else {
                            (beta)
                        }) / beta);
                    } else {
                        signbeta = (1_f64 as f32);
                    }
                    q = (((-0.5_f64) * ((beta + (signbeta * rad)) as f64)) as f32);
                    pz2 = (gamma / q);
                    if (pz2 < (0i32 as f32)) {
                        {
                            ipx += 1;
                            px = (px + dpx)
                        };
                        continue;
                    }
                    pz = ((pz2 as f64).sqrt() as f32);
                    f7 = ((2i32 as f32)
                        * (a1113[(ix) as usize][(0i32) as usize]
                            + a3313[(ix) as usize][(0i32) as usize]));
                    f8 = ((2i32 as f32)
                        * (a1113[(ix) as usize][(0i32) as usize]
                            + a3313[(ix) as usize][(0i32) as usize]));
                    f9 = a3313[(ix) as usize][(0i32) as usize];
                    f10 = a1113[(ix) as usize][(0i32) as usize];
                    f11 = ((2i32 as f32)
                        * (a1133[(ix) as usize][(0i32) as usize]
                            + a1313[(ix) as usize][(0i32) as usize]));
                    aa = (((f2 * f2) - (f4 * f4)) - (f9 * f9));
                    bb = ((((((2i32 as f32) * f2) * f7) + (((2i32 as f32) * f4) * f8))
                        - (((2i32 as f32) * f11) * f9))
                        * px);
                    ccc = (((((((f7 * f7) * px2) - (((2i32 as f32) * (f6 - (f1 * px2))) * f2))
                        - ((f8 * f8) * px2))
                        + ((((2i32 as f32) * f3) * f4) * px2))
                        - ((f11 * f11) * px2))
                        - ((((2i32 as f32) * f10) * f9) * px2));
                    dd = (((((((-2i32) as f32) * (f6 - (f1 * px2))) * f7) * px)
                        - (((((2i32 as f32) * f3) * f8) * px2) * px))
                        - (((((2i32 as f32) * f11) * f10) * px) * px2));
                    ee = ((((f6 - (f1 * px2)) * (f6 - (f1 * px2))) - (((f3 * f3) * px2) * px2))
                        - (((f10 * f10) * px2) * px2));
                    {
                        i = 0i32;
                        while (i < 20i32) {
                            fxx = ((((((aa * pz2) * pz2) + ((bb * pz2) * pz)) + (ccc * pz2))
                                + (dd * pz))
                                + ee);
                            dfx = (((((((4i32 as f32) * aa) * pz) * pz2)
                                + (((3i32 as f32) * bb) * pz2))
                                + (((2i32 as f32) * ccc) * pz))
                                + dd);
                            if (dfx == (0i32 as f32)) {
                                break;
                            }
                            dpz = ((-fxx) / dfx);
                            if (((if ((dpz) < (0i32 as f32)) {
                                (-(dpz))
                            } else {
                                (dpz)
                            }) as f64)
                                < 9.9999999999999995E-7_f64)
                            {
                                break;
                            }
                            pz = (pz + dpz);
                            pz2 = (pz * pz);
                            i += 1;
                        }
                    }
                    if (pz < (0i32 as f32)) {
                        {
                            ipx += 1;
                            px = (px + dpx)
                        };
                        continue;
                    }
                    vp = (((1i32 as f64) / ((px2 + pz2) as f64).sqrt()) as f32);
                    /* sine of emergence angle; skip if out of bounds */
                    if (((px * vp) as f64) > ((amax as f64).sin() + 0.01_f64)) {
                        {
                            ipx += 1;
                            px = (px + dpx)
                        };
                        continue;
                    }
                    if (((px * vp) as f64) < ((amin as f64).sin() - 0.01_f64)) {
                        {
                            ipx += 1;
                            px = (px + dpx)
                        };
                        continue;
                    }
                    /* emergence angle and location */
                    a0 = ((-((px * vp) as f64).asin()) as f32);
                    x0 = (fx + ((ix as f32) * dx));
                    z0 = fz;
                    /* beam half-width adjusted for cosine of angle */
                    bwhc = (((bwh as f64) * (a0 as f64).cos()) as f32);
                    /* trace ray */
                    let ray = makeRay(
                        x0, z0, a0, nt, dt, ft, a1111, a3333, a1133, a1313, a1113, a3313, nx, dx,
                        fx, nz, dz, fz,
                    );
                    /* accumulate contribution of beam in migrated image */
                    accBeam(
                        ray,
                        fmin,
                        bwhc,
                        ntau,
                        dtau,
                        ftau,
                        &b[ipx as usize],
                        nx,
                        dx,
                        fx,
                        nz,
                        dz,
                        fz,
                        g,
                    );
                    /* free ray */
                    {
                        ipx += 1;
                        px = (px + dpx)
                    };
                }
            }
            /* free space for beams */
            {
                ixb += 1;
                xb = (xb + dxb)
            };
        }
    }
}
