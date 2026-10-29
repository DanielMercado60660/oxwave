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

//! Numerical helpers from par/main/ray_theory/rayt2dtopo.c.
//! Copyright (c) Colorado School of Mines, 1995. All rights reserved.
//! Author: Zhenyue Liu, Colorado School of Mines, 1995.
//! C implicit arithmetic conversions are explicit to retain float rounding.
//!
//! Safety: pointers represent the original contiguous C arrays. Callers provide
//! velocity grids of nx*nz floats, tables of nxo*nzo floats, rays of nt elements,
//! and valid scalar output pointers; npv gates optional arrays. Vec owners stay
//! alive throughout the calls. Raw pointers permit the original sum2 calls to
//! use the same scalar for an input value and an output address.

#![allow(
    non_snake_case,
    unused_parens,
    unused_variables,
    unused_assignments,
    unused_mut
)]
use super::{Geo2d, Ray, Raytr};

pub(super) unsafe fn makeRay(
    mut geo2dv: Geo2d,
    mut vel: *mut f32,
    mut vxx: *mut f32,
    mut vxz: *mut f32,
    mut vzz: *mut f32,
    mut raytr: Raytr,
    mut a0: f32,
    mut nrs: *mut i32,
    mut ray: *mut Ray,
    mut npv: i32,
    mut pv: *mut f32,
) {
    let mut nx: i32 = geo2dv.nx;
    let mut nz: i32 = geo2dv.nz;
    let mut nt: i32 = raytr.nt;
    let mut fx: f32 = geo2dv.fx;
    let mut fz: f32 = geo2dv.fz;
    let mut dt: f32 = raytr.dt;
    let mut amin: f32 = raytr.amin;
    let mut amax: f32 = raytr.amax;
    let mut it: i32 = Default::default();
    /* float tzmin,tzmax,tx,ov,vv,lx,lz,h,h2,h6; */
    let mut tzmin: f32 = Default::default();
    let mut tzmax: f32 = Default::default();
    let mut ov: f32 = Default::default();
    let mut vv: f32 = Default::default();
    let mut lx: f32 = Default::default();
    let mut lz: f32 = Default::default();
    let mut h: f32 = Default::default();
    let mut h2: f32 = Default::default();
    let mut h6: f32 = Default::default();
    let mut v: f32 = Default::default();
    let mut dvdx: f32 = Default::default();
    let mut dvdz: f32 = Default::default();
    let mut uxx: f32 = Default::default();
    let mut uxz: f32 = Default::default();
    let mut uzz: f32 = Default::default();
    let mut tzt: f32 = Default::default();
    let mut dtv: f32 = Default::default();
    let mut x: f32 = Default::default();
    let mut z: f32 = Default::default();
    let mut px: f32 = Default::default();
    let mut pz: f32 = Default::default();
    let mut p2: f32 = Default::default();
    let mut q2: f32 = Default::default();
    let mut dx: f32 = Default::default();
    let mut dz: f32 = Default::default();
    let mut dpx: f32 = Default::default();
    let mut dpz: f32 = Default::default();
    let mut dp2: f32 = Default::default();
    let mut dq2: f32 = Default::default();
    let mut xt: f32 = Default::default();
    let mut zt: f32 = Default::default();
    let mut pxt: f32 = Default::default();
    let mut pzt: f32 = Default::default();
    let mut p2t: f32 = Default::default();
    let mut q2t: f32 = Default::default();
    let mut dxt: f32 = Default::default();
    let mut dzt: f32 = Default::default();
    let mut dpxt: f32 = Default::default();
    let mut dpzt: f32 = Default::default();
    let mut dp2t: f32 = Default::default();
    let mut dq2t: f32 = Default::default();
    /* velocity and derivatives at source	*/
    x = raytr.xs;
    z = raytr.zs;
    let _ = vel2Interp(
        geo2dv,
        vel,
        vxx,
        vxz,
        vzz,
        x,
        z,
        &raw mut v,
        &raw mut dvdx,
        &raw mut dvdz,
        &raw mut uxx,
        &raw mut uxz,
        &raw mut uzz,
    );
    ov = (1.0f64 / v as f64) as f32;
    vv = v * v;
    p2 = 1.0f64 as f32;
    q2 = 0.0f64 as f32;
    *nrs = nt;
    /*	compute slowness at source	*/
    px = (ov as f64 * (a0 as f64).sin()) as f32;
    pz = (ov as f64 * (a0 as f64).cos()) as f32;
    /*	first ray step output	*/
    tzt = v * pz;
    (*ray.offset(0 as isize)).x = x;
    (*ray.offset(0 as isize)).z = z;
    (*ray.offset(0 as isize)).px = px;
    (*ray.offset(0 as isize)).pz = pz;
    (*ray.offset(0 as isize)).q2 = q2;
    (*ray.offset(0 as isize)).p2 = p2;
    (*ray.offset(0 as isize)).v = v;
    (*ray.offset(0 as isize)).dvdx = dvdx;
    (*ray.offset(0 as isize)).dvdz = dvdz;
    if (npv != 0) {
        (*ray.offset(0 as isize)).dtv = 0 as f32;
    }
    /*	compute minimum and maximum z-component of unit ray vector	*/
    tzmin = ((amax as f64).cos() - 0.01f64) as f32;
    tzmax = ((amin as f64).cos() + 0.01f64) as f32;
    /*	determine fraction steps for Rung-Kuta	*/
    h = dt;
    h2 = (h as f64 * 0.5f64) as f32;
    h6 = (h as f64 / 6.0f64) as f32;
    lx = fx + (nx - 1) as f32 * geo2dv.dx;
    lz = fz + (nz - 1) as f32 * geo2dv.dz;
    /*	loop over time steps	*/
    it = 1;
    while it < nt {
        {
            if x < fx || x > lx || z < fz || z > lz || tzt < tzmin || tzt > tzmax {
                let _ = {
                    it = ((it as i32) - (1)) as i32;
                    it
                };
                break;
            }
            /*	step 1 of 4th-order Rung-Kuta	*/
            /*	  compute k1 = f(y0)	*/
            let _ = dfrungk(
                v,
                dvdx,
                dvdz,
                uxx,
                uxz,
                uzz,
                x,
                z,
                px,
                pz,
                p2,
                q2,
                &raw mut dx,
                &raw mut dz,
                &raw mut dpx,
                &raw mut dpz,
                &raw mut dp2,
                &raw mut dq2,
            );
            /*	  compute y1 = y0+0.5*h*k1	*/
            let _ = sum2(
                h2,
                x,
                z,
                px,
                pz,
                p2,
                q2,
                dx,
                dz,
                dpx,
                dpz,
                dp2,
                dq2,
                &raw mut xt,
                &raw mut zt,
                &raw mut pxt,
                &raw mut pzt,
                &raw mut p2t,
                &raw mut q2t,
            );
            /*	velocity interpolation	*/
            let _ = vel2Interp(
                geo2dv,
                vel,
                vxx,
                vxz,
                vzz,
                xt,
                zt,
                &raw mut v,
                &raw mut dvdx,
                &raw mut dvdz,
                &raw mut uxx,
                &raw mut uxz,
                &raw mut uzz,
            );
            /*	step 2 of 4th-order Rung-Kuta	*/
            /*	  compute k2 = f(y1)	*/
            let _ = dfrungk(
                v,
                dvdx,
                dvdz,
                uxx,
                uxz,
                uzz,
                xt,
                zt,
                pxt,
                pzt,
                p2t,
                q2t,
                &raw mut dxt,
                &raw mut dzt,
                &raw mut dpxt,
                &raw mut dpzt,
                &raw mut dp2t,
                &raw mut dq2t,
            );
            /*	  compute y2 = y0+0.5*h*k2	*/
            let _ = sum2(
                h2,
                x,
                z,
                px,
                pz,
                p2,
                q2,
                dxt,
                dzt,
                dpxt,
                dpzt,
                dp2t,
                dq2t,
                &raw mut xt,
                &raw mut zt,
                &raw mut pxt,
                &raw mut pzt,
                &raw mut p2t,
                &raw mut q2t,
            );
            /*	  compute k = k1+2.0*k2	*/
            let _ = sum2(
                2.0f64 as f32,
                dx,
                dz,
                dpx,
                dpz,
                dp2,
                dq2,
                dxt,
                dzt,
                dpxt,
                dpzt,
                dp2t,
                dq2t,
                &raw mut dx,
                &raw mut dz,
                &raw mut dpx,
                &raw mut dpz,
                &raw mut dp2,
                &raw mut dq2,
            );
            /*	velocity interpolation	*/
            let _ = vel2Interp(
                geo2dv,
                vel,
                vxx,
                vxz,
                vzz,
                xt,
                zt,
                &raw mut v,
                &raw mut dvdx,
                &raw mut dvdz,
                &raw mut uxx,
                &raw mut uxz,
                &raw mut uzz,
            );
            /*	step 3 of 4th-order Rung-Kuta	*/
            /*	  compute k3 = f(y2)	*/
            let _ = dfrungk(
                v,
                dvdx,
                dvdz,
                uxx,
                uxz,
                uzz,
                xt,
                zt,
                pxt,
                pzt,
                p2t,
                q2t,
                &raw mut dxt,
                &raw mut dzt,
                &raw mut dpxt,
                &raw mut dpzt,
                &raw mut dp2t,
                &raw mut dq2t,
            );
            /*	  compute y3 = y0+h*k2	*/
            let _ = sum2(
                h,
                x,
                z,
                px,
                pz,
                p2,
                q2,
                dxt,
                dzt,
                dpxt,
                dpzt,
                dp2t,
                dq2t,
                &raw mut xt,
                &raw mut zt,
                &raw mut pxt,
                &raw mut pzt,
                &raw mut p2t,
                &raw mut q2t,
            );
            /*	  compute k = k1+2.0*k2+2.0*k3	*/
            let _ = sum2(
                2.0f64 as f32,
                dx,
                dz,
                dpx,
                dpz,
                dp2,
                dq2,
                dxt,
                dzt,
                dpxt,
                dpzt,
                dp2t,
                dq2t,
                &raw mut dx,
                &raw mut dz,
                &raw mut dpx,
                &raw mut dpz,
                &raw mut dp2,
                &raw mut dq2,
            );
            /*	velocity interpolation	*/
            let _ = vel2Interp(
                geo2dv,
                vel,
                vxx,
                vxz,
                vzz,
                xt,
                zt,
                &raw mut v,
                &raw mut dvdx,
                &raw mut dvdz,
                &raw mut uxx,
                &raw mut uxz,
                &raw mut uzz,
            );
            /*	step 4 of 4th-order Rung-Kuta	*/
            /*	  compute k4 = f(y3)	*/
            let _ = dfrungk(
                v,
                dvdx,
                dvdz,
                uxx,
                uxz,
                uzz,
                xt,
                zt,
                pxt,
                pzt,
                p2t,
                q2t,
                &raw mut dxt,
                &raw mut dzt,
                &raw mut dpxt,
                &raw mut dpzt,
                &raw mut dp2t,
                &raw mut dq2t,
            );
            /*	  compute k = k1+2.0*k2+2.0*k3+k4	*/
            let _ = sum2(
                1.0f64 as f32,
                dx,
                dz,
                dpx,
                dpz,
                dp2,
                dq2,
                dxt,
                dzt,
                dpxt,
                dpzt,
                dp2t,
                dq2t,
                &raw mut dx,
                &raw mut dz,
                &raw mut dpx,
                &raw mut dpz,
                &raw mut dp2,
                &raw mut dq2,
            );
            /*	  compute y4 = y0+h*k/6	*/
            let _ = sum2(
                h6,
                x,
                z,
                px,
                pz,
                p2,
                q2,
                dx,
                dz,
                dpx,
                dpz,
                dp2,
                dq2,
                &raw mut x,
                &raw mut z,
                &raw mut px,
                &raw mut pz,
                &raw mut p2,
                &raw mut q2,
            );
            /*	velocity interpolation	*/
            let _ = vel2Interp(
                geo2dv,
                vel,
                vxx,
                vxz,
                vzz,
                x,
                z,
                &raw mut v,
                &raw mut dvdx,
                &raw mut dvdz,
                &raw mut uxx,
                &raw mut uxz,
                &raw mut uzz,
            );
            /*	velocity itself interpolation	*/
            if (npv != 0) {
                let _ = vintp(geo2dv, pv, x, z, &raw mut dtv);
            }
            tzt = v * pz;
            /*	save ray parameters */
            (*ray.offset(it as isize)).x = x;
            (*ray.offset(it as isize)).z = z;
            (*ray.offset(it as isize)).px = px;
            (*ray.offset(it as isize)).pz = pz;
            (*ray.offset(it as isize)).p2 = p2;
            (*ray.offset(it as isize)).q2 = q2;
            (*ray.offset(it as isize)).v = v;
            (*ray.offset(it as isize)).dvdx = dvdx;
            (*ray.offset(it as isize)).dvdz = dvdz;
            if (npv != 0) {
                (*ray.offset(it as isize)).dtv = dtv / v;
            }
        }
        it += 1;
    }
    *nrs = it;
}

pub(super) unsafe fn dfrungk(
    mut v: f32,
    mut vx: f32,
    mut vz: f32,
    mut vxx: f32,
    mut vxz: f32,
    mut vzz: f32,
    mut x: f32,
    mut z: f32,
    mut px: f32,
    mut pz: f32,
    mut p2: f32,
    mut q2: f32,
    mut dx: *mut f32,
    mut dz: *mut f32,
    mut dpx: *mut f32,
    mut dpz: *mut f32,
    mut dp2: *mut f32,
    mut dq2: *mut f32,
) {
    let mut ov: f32 = Default::default();
    let mut vv: f32 = Default::default();
    let mut c: f32 = Default::default();
    let mut s: f32 = Default::default();
    let mut v11: f32 = Default::default();
    ov = (1.0f64 / v as f64) as f32;
    vv = v * v;
    *dx = vv * px;
    *dz = vv * pz;
    *dpx = -ov * vx;
    *dpz = -ov * vz;
    c = v * pz;
    s = v * px;
    v11 = ((vxx * c * c) as f64 - 2.0f64 * vxz as f64 * c as f64 * s as f64 + (vzz * s * s) as f64)
        as f32;
    *dp2 = -ov * v11 * q2;
    *dq2 = vv * p2;
}

pub(super) unsafe fn sum2(
    mut h: f32,
    mut a1: f32,
    mut a2: f32,
    mut a3: f32,
    mut a4: f32,
    mut a5: f32,
    mut a6: f32,
    mut b1: f32,
    mut b2: f32,
    mut b3: f32,
    mut b4: f32,
    mut b5: f32,
    mut b6: f32,
    mut c1: *mut f32,
    mut c2: *mut f32,
    mut c3: *mut f32,
    mut c4: *mut f32,
    mut c5: *mut f32,
    mut c6: *mut f32,
) {
    *c1 = a1 + h * b1;
    *c2 = a2 + h * b2;
    *c3 = a3 + h * b3;
    *c4 = a4 + h * b4;
    *c5 = a5 + h * b5;
    *c6 = a6 + h * b6;
}

pub(super) unsafe fn raytime2d(
    mut raytr: Raytr,
    mut geo2dv: Geo2d,
    mut v: *mut f32,
    mut geo2dt: Geo2d,
    mut t: *mut f32,
    mut npv: i32,
    mut vo: *mut f32,
    mut pv: *mut f32,
    mut tv: *mut f32,
    mut cs: *mut f32,
) {
    // Nonpositive time counts make C overrun or dereference a failed ray
    // allocation. The tested multi-ray cases terminate with SIGSEGV; use
    // an explicit signal instead of an invalid Rust pointer access.
    if raytr.nt <= 0 && raytr.na > 0 {
        super::fault(11);
    }
    let mut na: i32 = raytr.na;
    let mut nt: i32 = raytr.nt;
    let mut nxo: i32 = geo2dt.nx;
    let mut nzo: i32 = geo2dt.nz;
    let mut nx: i32 = geo2dv.nx;
    let mut nz: i32 = geo2dv.nz;
    let mut fxo: f32 = geo2dt.fx;
    let mut fzo: f32 = geo2dt.fz;
    let mut dxo: f32 = geo2dt.dx;
    let mut dzo: f32 = geo2dt.dz;
    let mut odxo: f32 = geo2dt.odx;
    let mut odzo: f32 = geo2dt.odz;
    let mut fac: f32 = raytr.fac;
    let mut dt: f32 = raytr.dt;
    let mut fa: f32 = raytr.fa;
    let mut da: f32 = raytr.da;
    let mut xs: f32 = raytr.xs;
    let mut vxx: *mut f32 = std::ptr::null_mut();
    let mut vxz: *mut f32 = std::ptr::null_mut();
    let mut vzz: *mut f32 = std::ptr::null_mut();
    let mut s: *mut f32 = std::ptr::null_mut();
    let mut zs: f32 = raytr.zs;
    let mut nrs: i32 = Default::default();
    let mut ray: *mut Ray = std::ptr::null_mut();
    let mut i: i32 = Default::default();
    let mut ia: i32 = Default::default();
    let mut ixo: i32 = Default::default();
    let mut izo: i32 = Default::default();
    let mut a: f32 = Default::default();
    let mut xo: f32 = Default::default();
    let mut zo: f32 = Default::default();
    let mut exo: f32 = Default::default();
    let mut ezo: f32 = Default::default();
    /*	variables used for extrapolation	*/
    let mut it: i32 = Default::default();
    let mut jt: i32 = Default::default();
    let mut nxf: i32 = Default::default();
    let mut nxe: i32 = Default::default();
    let mut nzf: i32 = Default::default();
    let mut nze: i32 = Default::default();
    let mut tc: f32 = Default::default();
    let mut xoc: f32 = Default::default();
    let mut zoc: f32 = Default::default();
    let mut xc: f32 = Default::default();
    let mut zc: f32 = Default::default();
    let mut r1: f32 = Default::default();
    let mut v0: f32 = Default::default();
    let mut norm2: f32 = Default::default();
    let mut terr: f32 = Default::default();
    let mut vd1: f32 = Default::default();
    let mut cuv: f32 = Default::default();
    let mut sd: f32 = Default::default();
    let mut sind: f32 = Default::default();
    let mut cosd: f32 = Default::default();
    let mut vd: f32 = Default::default();
    let mut pxd: f32 = Default::default();
    let mut pzd: f32 = Default::default();
    let mut r2: f32 = Default::default();
    let mut t1: f32 = Default::default();
    let mut t2: f32 = Default::default();
    let mut r1x: f32 = Default::default();
    let mut r2x: f32 = Default::default();
    let mut t1x: f32 = Default::default();
    let mut t2x: f32 = Default::default();
    let mut t2xz: f32 = Default::default();
    let mut p2d: f32 = Default::default();
    let mut q2d: f32 = Default::default();
    let mut gradv: [f32; 2] = Default::default();
    let mut d2t: [f32; 3] = Default::default();
    let mut dtvd: f32 = Default::default();
    let mut odt: f32 = (1.0f64 / dt as f64) as f32;
    let mut xcosd: f32 = Default::default();
    let mut tvd: *mut f32 = std::ptr::null_mut();
    let mut vxx_storage = vec![f32::default(); ((nx * nz) as usize) as usize];
    vxx = vxx_storage.as_mut_ptr();
    let mut vxz_storage = vec![f32::default(); ((nx * nz) as usize) as usize];
    vxz = vxz_storage.as_mut_ptr();
    let mut vzz_storage = vec![f32::default(); ((nx * nz) as usize) as usize];
    vzz = vzz_storage.as_mut_ptr();
    let mut tvd_storage = vec![f32::default(); (nt as usize) as usize];
    tvd = tvd_storage.as_mut_ptr();
    let mut s_storage = vec![f32::default(); ((nxo * nzo) as usize) as usize];
    s = s_storage.as_mut_ptr();
    let mut ray_storage = vec![Ray::default(); (nt as usize) as usize];
    ray = ray_storage.as_mut_ptr();
    /* compute second derivatives of velocity	*/
    let _ = dv2(geo2dv, v, vxx, vxz, vzz);
    /*	maximum range of output points	*/
    exo = fxo + (nxo - 1) as f32 * dxo;
    ezo = fzo + (nzo - 1) as f32 * dzo;
    /*	initialize traveltime and raypath length	*/
    ixo = 0;
    while ixo < nxo {
        izo = 0;
        while izo < nzo {
            {
                i = izo + ixo * nzo;
                *t.offset(i as isize) = 999999 as f32;
                *s.offset(i as isize) = 999999 as f32;
                if (npv != 0) {
                    *tv.offset(i as isize) = {
                        *cs.offset(i as isize) = 0.0f64 as f32;
                        *cs.offset(i as isize)
                    };
                }
            }
            izo += 1;
        }
        ixo += 1;
    }
    /* 	loop over shooting angles at source point 	*/
    ia = 0;
    a = fa;
    while ia < na {
        {
            /*		trace rays	*/
            let _ = makeRay(
                geo2dv,
                v,
                vxx,
                vxz,
                vzz,
                raytr,
                a,
                &raw mut nrs,
                ray,
                npv,
                pv,
            );
            /*		extropolate to grids near central rays	*/
            v0 = (*ray.offset(0 as isize)).v;
            r2 = 0.0f64 as f32;
            xc = raytr.xs;
            zc = raytr.zs;
            sd = 0 as f32;
            vd1 = v0;
            if (npv != 0) {
                dtvd = 0.0f64 as f32;
                it = 1;
                while it < nrs {
                    {
                        let _ = {
                            dtvd = ((dtvd as f64)
                                + (0.5f64
                                    * dt as f64
                                    * ((*ray.offset(it as isize)).dtv
                                        + (*ray.offset((it - 1) as isize)).dtv)
                                        as f64)) as f32;
                            dtvd
                        };
                        *tvd.offset(it as isize) = (*ray.offset(it as isize)).v * dtvd;
                    }
                    it += 1;
                }
            }
            it = 1;
            while it < nrs {
                {
                    xo = (*ray.offset(it as isize)).x;
                    zo = (*ray.offset(it as isize)).z;
                    vd = (*ray.offset(it as isize)).v;
                    sd = (sd as f64 + 0.5f64 * dt as f64 * (vd + vd1) as f64) as f32;
                    vd1 = vd;
                    /*		    	if the point is within the previous circle	*/
                    if r2 > (xc - xo) * (xc - xo) + (zc - zo) * (zc - zo) && it != nrs - 1 {
                        it += 1;
                        continue;
                    }
                    /*		    	if the point is out of output range	*/
                    if xo < fxo || xo > exo || zo < fzo || zo > ezo {
                        it += 1;
                        continue;
                    }
                    xc = xo;
                    zc = zo;
                    q2d = (*ray.offset(it as isize)).q2;
                    /*			if caustics	*/
                    if q2d as f64 == 0.0f64 {
                        r2 = 0.0f64 as f32;
                        it += 1;
                        continue;
                    }
                    vd = (*ray.offset(it as isize)).v;
                    pxd = (*ray.offset(it as isize)).px;
                    pzd = (*ray.offset(it as isize)).pz;
                    p2d = (*ray.offset(it as isize)).p2;
                    sind = vd * pxd;
                    cosd = vd * pzd;
                    *gradv.as_mut_ptr().offset(0 as isize) = (*ray.offset(it as isize)).dvdx;
                    *gradv.as_mut_ptr().offset(1 as isize) = (*ray.offset(it as isize)).dvdz;
                    /*			calculate second derivatives of traveltime*/
                    let _ = ddt(
                        p2d,
                        q2d,
                        cosd,
                        sind,
                        gradv.as_mut_ptr(),
                        vd,
                        d2t.as_mut_ptr(),
                        &raw mut cuv,
                    );
                    /*			determine radius for extrapolation	*/
                    tc = it as f32 * dt;
                    terr = tc * fac;
                    norm2 = ((*d2t.as_mut_ptr().offset(0 as isize)
                        * *d2t.as_mut_ptr().offset(0 as isize)
                        + *d2t.as_mut_ptr().offset(2 as isize)
                            * *d2t.as_mut_ptr().offset(2 as isize))
                        as f64
                        + 2.0f64
                            * *d2t.as_mut_ptr().offset(1 as isize) as f64
                            * *d2t.as_mut_ptr().offset(1 as isize) as f64)
                        .sqrt() as f32;
                    r2 = terr / norm2;
                    r1 = (r2 as f64).sqrt() as f32;
                    /*			keep ray-centered coordinate system regular */
                    if (r1 * cuv) as f64 > 0.10000000000000001f64 {
                        r1 = (0.10000000000000001f64 / cuv as f64) as f32;
                    }
                    /*			radius cannot be too large	*/
                    if r1 as f64 > 0.10000000000000001f64 * sd as f64 {
                        r1 = (0.10000000000000001f64 * sd as f64) as f32;
                    }
                    r2 = r1 * r1;
                    nxf = (((xc - r1 - fxo) * odxo) as f64 + 0.90000000000000002f64) as i32;
                    if nxf < 0 {
                        nxf = 0;
                    }
                    nxe = (((xc + r1 - fxo) * odxo) as f64 + 0.10000000000000001f64) as i32;
                    if nxe >= nxo {
                        nxe = nxo - 1;
                    }
                    /*       	fprintf(stderr,"x,z,t,r=%g %g %g %g\n",
                    xc,zc,tc,r1);*/
                    ixo = nxf;
                    while ixo <= nxe {
                        {
                            xoc = fxo - xc + ixo as f32 * dxo;
                            r2x = r2 - xoc * xoc;
                            if r2x < 0 as f32 {
                                ixo += 1;
                                continue;
                            }
                            r1x = (r2x as f64).sqrt() as f32;
                            t1x = tc + xoc * pxd;
                            t2x = tc * xoc * xoc * *d2t.as_mut_ptr().offset(0 as isize);
                            t2xz =
                                (2.0f64 * xoc as f64 * *d2t.as_mut_ptr().offset(1 as isize) as f64)
                                    as f32;
                            if (npv != 0) {
                                xcosd = pzd + xoc * *d2t.as_mut_ptr().offset(1 as isize);
                            }
                            nzf =
                                (((zc - r1x - fzo) * odzo) as f64 + 0.90000000000000002f64) as i32;
                            if nzf < 0 {
                                nzf = 0;
                            }
                            nze =
                                (((zc + r1x - fzo) * odzo) as f64 + 0.10000000000000001f64) as i32;
                            if nze >= nzo {
                                nze = nzo - 1;
                            }
                            i = ixo * nzo;
                            izo = nzf;
                            while izo <= nze {
                                {
                                    /*				    if ray path is shorter	*/
                                    if sd < *s.offset((i + izo) as isize) {
                                        zoc = fzo - zc + izo as f32 * dzo;
                                        t1 = t1x + zoc * pzd;
                                        t2 = t2x
                                            + zoc
                                                * (zoc * *d2t.as_mut_ptr().offset(2 as isize)
                                                    + t2xz);
                                        *s.offset((i + izo) as isize) = sd;
                                        *t.offset((i + izo) as isize) = t1 * t1 + tc * t2;
                                        if (npv != 0) {
                                            jt = ((t1 as f64 + 0.5f64 * t2 as f64) * odt as f64
                                                + 0.5f64)
                                                as i32;
                                            if jt < 0 {
                                                jt = 0;
                                            }
                                            if jt > nrs - 1 {
                                                jt = nrs - 1;
                                            }
                                            *tv.offset((i + izo) as isize) =
                                                *tvd.offset(jt as isize);
                                            *cs.offset((i + izo) as isize) =
                                                xcosd + *d2t.as_mut_ptr().offset(2 as isize) * zoc;
                                        }
                                    }
                                }
                                izo += 1;
                            }
                        }
                        ixo += 1;
                    }
                }
                it += 1;
            }
        }
        ia += 1;
        let _ = {
            a = ((a as f32) + (da)) as f32;
            a
        };
    }
    /*	square root of traveltimes */
    ixo = 0;
    while ixo < nxo {
        {
            i = ixo * nzo;
            izo = 0;
            while izo < nzo {
                {
                    *t.offset((i + izo) as isize) = if 0.0f64 > *t.offset((i + izo) as isize) as f64
                    {
                        0.0f64
                    } else {
                        *t.offset((i + izo) as isize) as f64
                    } as f32;
                    if *t.offset((i + izo) as isize) < 999 as f32 {
                        *t.offset((i + izo) as isize) =
                            (*t.offset((i + izo) as isize) as f64).sqrt() as f32;
                    }
                    if (npv != 0) {
                        let _ = {
                            *cs.offset((i + izo) as isize) = ((*cs.offset((i + izo) as isize)
                                as f32)
                                * (*vo.offset((i + izo) as isize)))
                                as f32;
                            *cs.offset((i + izo) as isize)
                        };
                    }
                }
                izo += 1;
            }
        }
        ixo += 1;
    }
    /*	compute traveltime near source	*/
    nxf = (((xs - fxo) / dxo) as f64 - 1.5f64) as i32;
    if nxf < 0 {
        nxf = 0;
    }
    nxe = (((xs - fxo) / dxo) as f64 + 2.5f64) as i32;
    if nxe > nxo - 1 {
        nxe = nxo - 1;
    }
    nzf = (((zs - fzo) / dzo) as f64 - 0.5f64) as i32;
    if nzf < 0 {
        nzf = 0;
    }
    nze = (((zs - fzo) / dzo) as f64 + 1.5f64) as i32;
    if nze >= nzo {
        nze = nzo - 1;
    }
    ixo = nxf;
    while ixo <= nxe {
        {
            xo = fxo + ixo as f32 * dxo;
            i = ixo * nzo;
            izo = nzf;
            while izo <= nze {
                {
                    zo = fzo + izo as f32 * dzo;
                    if *t.offset((i + izo) as isize) > 999 as f32 {
                        *t.offset((i + izo) as isize) =
                            ((((xo - xs) * (xo - xs) + (zo - zs) * (zo - zs)) as f64).sqrt()
                                / v0 as f64) as f32;
                    }
                }
                izo += 1;
            }
        }
        ixo += 1;
    }
    let _ = ();
    let _ = ();
    let _ = ();
    let _ = ();
    let _ = ();
}

pub(super) unsafe fn ddt(
    mut p: f32,
    mut q: f32,
    mut c: f32,
    mut s: f32,
    mut dv: *mut f32,
    mut v: f32,
    mut d2t: *mut f32,
    mut cuv: *mut f32,
) {
    let mut ov2: f32 = Default::default();
    let mut m11: f32 = Default::default();
    let mut m12: f32 = Default::default();
    let mut m22: f32 = Default::default();
    ov2 = (1.0f64 / (v * v) as f64) as f32;
    m11 = p / q;
    /*	calculate 2nd column of m	*/
    m12 = -(*dv.offset(0 as isize) * c - *dv.offset(1 as isize) * s) * ov2;
    m22 = -(*dv.offset(0 as isize) * s + *dv.offset(1 as isize) * c) * ov2;
    /*	compute h*m*h^T	*/
    *d2t.offset(0 as isize) = ((m11 * c * c) as f64
        + 2.0f64 * m12 as f64 * c as f64 * s as f64
        + (m22 * s * s) as f64) as f32;
    *d2t.offset(1 as isize) = (m22 - m11) * c * s + m12 * (c * c - s * s);
    *d2t.offset(2 as isize) = ((m11 * s * s) as f64 - 2.0f64 * m12 as f64 * c as f64 * s as f64
        + (m22 * c * c) as f64) as f32;
    /*	compute the curvature of raypath	*/
    *cuv = ((m12 as f64).abs() * v as f64) as f32;
}

pub(super) unsafe fn trans(
    mut nx: i32,
    mut nz: i32,
    mut nxt: i32,
    mut nx0: i32,
    mut v: *mut f32,
    mut vt: *mut f32,
) {
    let mut ix: i32 = Default::default();
    let mut iz: i32 = Default::default();
    let mut i: i32 = Default::default();
    let mut i0: i32 = Default::default();
    ix = 0;
    while ix < nxt {
        {
            i = ix * nz;
            i0 = (ix + nx0) * nz;
            iz = 0;
            while iz < nz {
                *vt.offset((i + iz) as isize) = *v.offset((i0 + iz) as isize);
                iz += 1;
            }
        }
        ix += 1;
    }
}

pub(super) unsafe fn voint(
    mut geo2dv: Geo2d,
    mut v: *mut f32,
    mut geo2dt: Geo2d,
    mut ov2: *mut f32,
    mut npv: i32,
    mut vo: *mut f32,
) {
    let mut nx: i32 = geo2dv.nx;
    let mut nz: i32 = geo2dv.nz;
    let mut nxo: i32 = geo2dt.nx;
    let mut nzo: i32 = geo2dt.nz;
    /*	int ix,iz,i,io,jx,jz,ixo,izo; */
    let mut i: i32 = Default::default();
    let mut io: i32 = Default::default();
    let mut jx: i32 = Default::default();
    let mut jz: i32 = Default::default();
    let mut ixo: i32 = Default::default();
    let mut izo: i32 = Default::default();
    let mut fx: f32 = geo2dv.fx;
    let mut fz: f32 = geo2dv.fz;
    let mut odx: f32 = geo2dv.odx;
    let mut odz: f32 = geo2dv.odz;
    let mut fxo: f32 = geo2dt.fx;
    let mut fzo: f32 = geo2dt.fz;
    let mut dxo: f32 = geo2dt.dx;
    let mut dzo: f32 = geo2dt.dz;
    let mut x: f32 = Default::default();
    let mut z: f32 = Default::default();
    let mut ax: f32 = Default::default();
    let mut az: f32 = Default::default();
    let mut sx: f32 = Default::default();
    let mut sz: f32 = Default::default();
    let mut temp: f32 = Default::default();
    ixo = 0;
    x = fxo;
    while ixo < nxo {
        {
            ax = (x - fx) * odx;
            jx = ax as i32;
            sx = ax - jx as f32;
            if jx < 0 {
                jx = 0;
            }
            if jx > nx - 2 {
                jx = nx - 2;
            }
            izo = 0;
            z = fzo;
            while izo < nzo {
                {
                    az = (z - fz) * odz;
                    jz = az as i32;
                    sz = az - jz as f32;
                    if jz < 0 {
                        jz = 0;
                    }
                    if jz > nz - 2 {
                        jz = nz - 2;
                    }
                    io = ixo * nzo + izo;
                    i = jx * nz + jz;
                    temp = ((1.0f64 - sx as f64)
                        * ((1.0f64 - sz as f64) * *v.offset(i as isize) as f64
                            + (sz * *v.offset((i + 1) as isize)) as f64)
                        + sx as f64
                            * ((1.0f64 - sz as f64) * *v.offset((i + nz) as isize) as f64
                                + (sz * *v.offset((i + nz + 1) as isize)) as f64))
                        as f32;
                    if (npv != 0) {
                        *vo.offset(io as isize) = temp;
                    }
                    *ov2.offset(io as isize) = (1.0f64 / (temp * temp) as f64) as f32;
                }
                izo += 1;
                let _ = {
                    z = ((z as f32) + (dzo)) as f32;
                    z
                };
            }
        }
        ixo += 1;
        let _ = {
            x = ((x as f32) + (dxo)) as f32;
            x
        };
    }
}

pub(super) unsafe fn dv2(
    mut geo2dv: Geo2d,
    mut v: *mut f32,
    mut vxx: *mut f32,
    mut vxz: *mut f32,
    mut vzz: *mut f32,
) {
    let mut nx: i32 = geo2dv.nx;
    let mut nz: i32 = geo2dv.nz;
    let mut ix: i32 = Default::default();
    let mut iz: i32 = Default::default();
    let mut i: i32 = Default::default();
    let mut odx: f32 = geo2dv.odx;
    let mut odz: f32 = geo2dv.odz;
    let mut odxx: f32 = Default::default();
    let mut odzz: f32 = Default::default();
    let mut odxz: f32 = Default::default();
    odxx = odx * odx;
    odzz = odz * odz;
    odxz = (0.25f64 * odx as f64 * odz as f64) as f32;
    /*	initialize	*/
    ix = 0;
    while ix < nx {
        iz = 0;
        while iz < nz {
            {
                i = ix * nz + iz;
                *vxx.offset(i as isize) = 0.0f64 as f32;
                *vxz.offset(i as isize) = 0.0f64 as f32;
                *vzz.offset(i as isize) = 0.0f64 as f32;
            }
            iz += 1;
        }
        ix += 1;
    }
    /*	finite difference 	*/
    ix = 1;
    while ix < nx - 1 {
        iz = 1;
        while iz < nz - 1 {
            {
                i = ix * nz + iz;
                *vxx.offset(i as isize) = (odxx as f64
                    * (*v.offset((i + nz) as isize) as f64 - 2.0f64 * *v.offset(i as isize) as f64
                        + *v.offset((i - nz) as isize) as f64))
                    as f32;
                *vxz.offset(i as isize) = odxz
                    * (*v.offset((i + nz + 1) as isize)
                        - *v.offset((i + nz - 1) as isize)
                        - *v.offset((i - nz + 1) as isize)
                        + *v.offset((i - nz - 1) as isize));
                *vzz.offset(i as isize) = (odzz as f64
                    * (*v.offset((i + 1) as isize) as f64 - 2.0f64 * *v.offset(i as isize) as f64
                        + *v.offset((i - 1) as isize) as f64))
                    as f32;
            }
            iz += 1;
        }
        ix += 1;
    }
}

pub(super) unsafe fn eiknl(mut geo2dt: Geo2d, mut t: *mut f32, mut ov2: *mut f32, mut tmax: f32) {
    let mut nx: i32 = geo2dt.nx;
    let mut nz: i32 = geo2dt.nz;
    let mut ix: i32 = Default::default();
    let mut iz: i32 = Default::default();
    let mut dx: f32 = geo2dt.dx;
    let mut dz: f32 = geo2dt.dz;
    /* float tx2,tz2,t0,t1,tl,tr,temp,odx2,*tt1,*tt2; */
    let mut tx2: f32 = Default::default();
    let mut tz2: f32 = Default::default();
    let mut t0: f32 = Default::default();
    let mut tl: f32 = Default::default();
    let mut tr: f32 = Default::default();
    let mut temp: f32 = Default::default();
    let mut odx2: f32 = Default::default();
    let mut tt1: *mut f32 = std::ptr::null_mut();
    let mut tt2: *mut f32 = std::ptr::null_mut();
    let mut tt1_storage = vec![f32::default(); (nx as usize) as usize];
    tt1 = tt1_storage.as_mut_ptr();
    let mut tt2_storage = vec![f32::default(); (nx as usize) as usize];
    tt2 = tt2_storage.as_mut_ptr();
    odx2 = (1.0f64 / (dx * dx) as f64) as f32;
    ix = 0;
    while ix < nx {
        *tt1.offset(ix as isize) = *t.offset((ix * nz) as isize);
        ix += 1;
    }
    iz = 1;
    while iz < nz {
        {
            ix = 0;
            while ix < nx {
                *tt2.offset(ix as isize) = *t.offset((ix * nz + iz) as isize);
                ix += 1;
            }
            ix = 1;
            while ix < nx {
                {
                    /*	if traveltime is abnormal and the upper is normal	*/
                    t0 = *tt1.offset(ix as isize);
                    if *tt2.offset(ix as isize) > tmax && t0 <= tmax {
                        tl = {
                            tr = 0.0f64 as f32;
                            tr
                        };
                        if ix > 0 {
                            tl = *tt1.offset((ix - 1) as isize);
                        }
                        if ix < nx - 1 {
                            tr = *tt1.offset((ix + 1) as isize);
                        }
                        tx2 = (tl - t0) * (tl - t0);
                        temp = (tr - t0) * (tr - t0);
                        if tx2 > temp {
                            tx2 = temp;
                        }
                        temp = (0.25f64 * (tl - tr) as f64 * (tl - tr) as f64) as f32;
                        if tx2 > temp {
                            tx2 = temp;
                        }
                        let _ = {
                            tx2 = ((tx2 as f32) * (odx2)) as f32;
                            tx2
                        };
                        tz2 = (0.5f64
                            * (*ov2.offset((ix * nz + iz - 1) as isize)
                                + *ov2.offset((ix * nz + iz) as isize))
                                as f64
                            - tx2 as f64) as f32;
                        if tz2 >= 0 as f32 {
                            *tt2.offset(ix as isize) =
                                (t0 as f64 + dz as f64 * (tz2 as f64).sqrt()) as f32;
                        }
                    }
                }
                ix += 1;
            }
            ix = 0;
            while ix < nx {
                {
                    *t.offset((ix * nz + iz) as isize) = *tt2.offset(ix as isize);
                    *tt1.offset(ix as isize) = *tt2.offset(ix as isize);
                }
                ix += 1;
            }
        }
        iz += 1;
    }
    let _ = ();
    let _ = ();
}

pub(super) unsafe fn vel2Interp(
    mut geo2dv: Geo2d,
    mut v: *mut f32,
    mut vxx: *mut f32,
    mut vxz: *mut f32,
    mut vzz: *mut f32,
    mut x: f32,
    mut z: f32,
    mut u: *mut f32,
    mut ux: *mut f32,
    mut uz: *mut f32,
    mut uxx: *mut f32,
    mut uxz: *mut f32,
    mut uzz: *mut f32,
) {
    let mut k0: i32 = Default::default();
    let mut k1: i32 = Default::default();
    let mut k2: i32 = Default::default();
    let mut k3: i32 = Default::default();
    let mut jx: i32 = Default::default();
    let mut jz: i32 = Default::default();
    let mut nx: i32 = geo2dv.nx;
    let mut nz: i32 = geo2dv.nz;
    let mut fx: f32 = geo2dv.fx;
    let mut dx: f32 = geo2dv.dx;
    let mut odx: f32 = geo2dv.odx;
    let mut fz: f32 = geo2dv.fz;
    let mut dz: f32 = geo2dv.dz;
    let mut odz: f32 = geo2dv.odz;
    let mut ax: f32 = Default::default();
    let mut az: f32 = Default::default();
    let mut sx: f32 = Default::default();
    let mut sz: f32 = Default::default();
    let mut sxx: f32 = Default::default();
    let mut szz: f32 = Default::default();
    let mut a0: f32 = Default::default();
    let mut a1: f32 = Default::default();
    let mut a2: f32 = Default::default();
    let mut a3: f32 = Default::default();
    let mut g0: f32 = Default::default();
    let mut g1: f32 = Default::default();
    let mut g2: f32 = Default::default();
    let mut g3: f32 = Default::default();
    let mut gx0: f32 = Default::default();
    let mut gx1: f32 = Default::default();
    let mut gx2: f32 = Default::default();
    let mut gx3: f32 = Default::default();
    let mut gz0: f32 = Default::default();
    let mut gz1: f32 = Default::default();
    let mut gz2: f32 = Default::default();
    let mut gz3: f32 = Default::default();
    /*	determine interpolate coefficients	*/
    ax = (x - fx) * odx;
    jx = ax as i32;
    if jx < 0 {
        jx = 0;
    }
    if jx > nx - 2 {
        jx = nx - 2;
    }
    sx = ax - jx as f32;
    az = (z - fz) * odz;
    jz = az as i32;
    if jz < 0 {
        jz = 0;
    }
    if jz > nz - 2 {
        jz = nz - 2;
    }
    sz = az - jz as f32;
    sxx = (0.5f64 * sx as f64 * (1.0f64 - sx as f64) * dx as f64 * dx as f64) as f32;
    szz = (0.5f64 * sz as f64 * (1.0f64 - sz as f64) * dz as f64 * dz as f64) as f32;
    a0 = ((1.0f64 - sx as f64) * (1.0f64 - sz as f64)) as f32;
    a1 = ((1.0f64 - sx as f64) * sz as f64) as f32;
    a2 = (sx as f64 * (1.0f64 - sz as f64)) as f32;
    a3 = sx * sz;
    /*	    set the table of indices for interpolation	*/
    k0 = nz * jx + jz;
    k1 = k0 + 1;
    k2 = k0 + nz;
    k3 = k2 + 1;
    g0 = *v.offset(k0 as isize);
    g1 = *v.offset(k1 as isize);
    g2 = *v.offset(k2 as isize);
    g3 = *v.offset(k3 as isize);
    gx0 = *vxx.offset(k0 as isize);
    gx1 = *vxx.offset(k1 as isize);
    gx2 = *vxx.offset(k2 as isize);
    gx3 = *vxx.offset(k3 as isize);
    gz0 = *vzz.offset(k0 as isize);
    gz1 = *vzz.offset(k1 as isize);
    gz2 = *vzz.offset(k2 as isize);
    gz3 = *vzz.offset(k3 as isize);
    /*	interpolation	*/
    *uxx = a0 * gx0 + a1 * gx1 + a2 * gx2 + a3 * gx3;
    *uxz = a0 * *vxz.offset(k0 as isize)
        + a1 * *vxz.offset(k1 as isize)
        + a2 * *vxz.offset(k2 as isize)
        + a3 * *vxz.offset(k3 as isize);
    *uzz = a0 * gz0 + a1 * gz1 + a2 * gz2 + a3 * gz3;
    *u = a0 * g0 + a1 * g1 + a2 * g2 + a3 * g3 - (sxx * *uxx + szz * *uzz);
    *ux = (((1.0f64 - sz as f64) * (g2 - g0 - sxx * (gx2 - gx0) - szz * (gz2 - gz0)) as f64
        + (sz * (g3 - g1 - sxx * (gx3 - gx1) - szz * (gz3 - gz1))) as f64)
        * odx as f64
        + (sx as f64 - 0.5f64) * dx as f64 * *uxx as f64) as f32;
    *uz = (((1.0f64 - sx as f64) * (g1 - g0 - sxx * (gx1 - gx0) - szz * (gz1 - gz0)) as f64
        + (sx * (g3 - g2 - sxx * (gx3 - gx2) - szz * (gz3 - gz2))) as f64)
        * odz as f64
        + (sz as f64 - 0.5f64) * dz as f64 * *uzz as f64) as f32;
}

pub(super) unsafe fn vintp(
    mut geo2dv: Geo2d,
    mut v: *mut f32,
    mut x: f32,
    mut z: f32,
    mut u: *mut f32,
) {
    let mut k0: i32 = Default::default();
    let mut k1: i32 = Default::default();
    let mut k2: i32 = Default::default();
    let mut k3: i32 = Default::default();
    let mut jx: i32 = Default::default();
    let mut jz: i32 = Default::default();
    let mut nx: i32 = geo2dv.nx;
    let mut nz: i32 = geo2dv.nz;
    /*	float fx=geo2dv.fx,dx=geo2dv.dx,odx=geo2dv.odx;
        float fz=geo2dv.fz,dz=geo2dv.dz,odz=geo2dv.odz;
    */
    let mut fx: f32 = geo2dv.fx;
    let mut odx: f32 = geo2dv.odx;
    let mut fz: f32 = geo2dv.fz;
    let mut odz: f32 = geo2dv.odz;
    let mut ax: f32 = Default::default();
    let mut az: f32 = Default::default();
    let mut sx: f32 = Default::default();
    let mut sz: f32 = Default::default();
    /*	determine interpolate coefficients	*/
    ax = (x - fx) * odx;
    jx = ax as i32;
    if jx < 0 {
        jx = 0;
    }
    if jx > nx - 2 {
        jx = nx - 2;
    }
    sx = ax - jx as f32;
    az = (z - fz) * odz;
    jz = az as i32;
    if jz < 0 {
        jz = 0;
    }
    if jz > nz - 2 {
        jz = nz - 2;
    }
    sz = az - jz as f32;
    /*	    set the table of indices for interpolation	*/
    k0 = nz * jx + jz;
    k1 = k0 + 1;
    k2 = k0 + nz;
    k3 = k2 + 1;
    /*	interpolation	*/
    *u = ((1.0f64 - sx as f64)
        * ((1.0f64 - sz as f64) * *v.offset(k0 as isize) as f64
            + (sz * *v.offset(k1 as isize)) as f64)
        + sx as f64
            * ((1.0f64 - sz as f64) * *v.offset(k2 as isize) as f64
                + (sz * *v.offset(k3 as isize)) as f64)) as f32;
}
