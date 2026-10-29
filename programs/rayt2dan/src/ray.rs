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

//! Ray integration from par/main/ray_theory/rayt2dan.c.
//! Credits: CWP: Tariq Alkhalifah.

#![allow(
    non_snake_case,
    unused_mut,
    unused_variables,
    unused_assignments,
    unused_parens
)]
use crate::vel::Vel2;

#[derive(Clone, Copy, Default)]
pub struct RayStep {
    pub t: f32,
    pub x: f32,
    pub z: f32,
    pub q1: f32,
    pub p1: f32,
    pub q2: f32,
    pub p2: f32,
    pub kmah: i32,
    pub c: f32,
    pub s: f32,
    pub v: f32,
    pub dvdx: f32,
    pub dvdz: f32,
}
#[derive(Default)]
pub struct Ray {
    pub nrs: i32,
    pub rs: Vec<RayStep>,
}

pub fn crash() -> ! {
    // The original dereferences a NULL ray when a takeoff point lies outside the grid.
    unsafe extern "C" {
        fn signal(sig: i32, handler: usize) -> usize;
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        signal(11, 0);
        raise(11);
    }
    std::process::exit(139)
}

pub fn make_ray(
    mut SV: i32,
    mut x0: f32,
    mut z0: f32,
    mut a0: f32,
    mut nt: i32,
    mut dt: f32,
    mut ft: f32,
    mut a1111xz: &[f32],
    mut a3333xz: &[f32],
    mut a1133xz: &[f32],
    mut a1313xz: &[f32],
    mut a1113xz: &[f32],
    mut a3313xz: &[f32],
    mut nx: i32,
    mut dx: f32,
    mut fx: f32,
    mut nz: i32,
    mut dz: f32,
    mut fz: f32,
    mut amax: f32,
    mut amin: f32,
) -> Ray {
    let mut it: i32 = 0;
    let mut kmah: i32 = 0;
    let mut t: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut c: f32 = 0.0;
    let mut s: f32 = 0.0;
    let mut p1: f32 = 0.0;
    let mut q1: f32 = 0.0;
    let mut p2: f32 = 0.0;
    let mut q2: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut px2: f32 = 0.0;
    let mut pz2: f32 = 0.0;
    let mut pxz: f32 = 0.0;
    let mut lx: f32 = 0.0;
    let mut lz: f32 = 0.0;
    let mut cc: f32 = 0.0;
    let mut ss: f32 = 0.0;
    let mut a1111: f32 = 0.0;
    let mut da1111dx: f32 = 0.0;
    let mut da1111dz: f32 = 0.0;
    let mut dda1111dxdx: f32 = 0.0;
    let mut dda1111dxdz: f32 = 0.0;
    let mut dda1111dzdz: f32 = 0.0;
    let mut a3333: f32 = 0.0;
    let mut da3333dx: f32 = 0.0;
    let mut da3333dz: f32 = 0.0;
    let mut dda3333dxdx: f32 = 0.0;
    let mut dda3333dxdz: f32 = 0.0;
    let mut dda3333dzdz: f32 = 0.0;
    let mut a1133: f32 = 0.0;
    let mut da1133dx: f32 = 0.0;
    let mut da1133dz: f32 = 0.0;
    let mut dda1133dxdx: f32 = 0.0;
    let mut dda1133dxdz: f32 = 0.0;
    let mut dda1133dzdz: f32 = 0.0;
    let mut a1313: f32 = 0.0;
    let mut da1313dx: f32 = 0.0;
    let mut da1313dz: f32 = 0.0;
    let mut dda1313dxdx: f32 = 0.0;
    let mut dda1313dxdz: f32 = 0.0;
    let mut dda1313dzdz: f32 = 0.0;
    let mut a1113: f32 = 0.0;
    let mut da1113dx: f32 = 0.0;
    let mut da1113dz: f32 = 0.0;
    let mut dda1113dxdx: f32 = 0.0;
    let mut dda1113dxdz: f32 = 0.0;
    let mut dda1113dzdz: f32 = 0.0;
    let mut a3313: f32 = 0.0;
    let mut da3313dx: f32 = 0.0;
    let mut da3313dz: f32 = 0.0;
    let mut dda3313dxdx: f32 = 0.0;
    let mut dda3313dxdz: f32 = 0.0;
    let mut dda3313dzdz: f32 = 0.0;
    let mut da1111dn: f32 = 0.0;
    let mut dda1111dndn: f32 = 0.0;
    let mut da3333dn: f32 = 0.0;
    let mut dda3333dndn: f32 = 0.0;
    let mut da1133dn: f32 = 0.0;
    let mut dda1133dndn: f32 = 0.0;
    let mut da1313dn: f32 = 0.0;
    let mut dda1313dndn: f32 = 0.0;
    let mut da1113dn: f32 = 0.0;
    let mut dda1113dndn: f32 = 0.0;
    let mut da3313dn: f32 = 0.0;
    let mut dda3313dndn: f32 = 0.0;
    let mut gamm11: f32 = 0.0;
    let mut gamm13: f32 = 0.0;
    let mut gamm33: f32 = 0.0;
    let mut vp2: f32 = 0.0;
    let mut vp: f32 = 0.0;
    let mut ovp: f32 = 0.0;
    let mut sqr: f32 = 0.0;
    /*Convert degrees to radians*/
    a0 = ((((a0 as f64) * 3.1415f64) / (180i32 as f64)) as f32);
    amax = ((((amax as f64) * 3.1415f64) / (180i32 as f64)) as f32);
    amin = ((((amin as f64) * 3.1415f64) / (180i32 as f64)) as f32);
    /* allocate and initialize velocities interpolator */
    let a11112 = Vel2::new(nx, dx, fx, nz, dz, fz, &a1111xz);
    let a33332 = Vel2::new(nx, dx, fx, nz, dz, fz, &a3333xz);
    let a11332 = Vel2::new(nx, dx, fx, nz, dz, fz, &a1133xz);
    let a13132 = Vel2::new(nx, dx, fx, nz, dz, fz, &a1313xz);
    let a11132 = Vel2::new(nx, dx, fx, nz, dz, fz, &a1113xz);
    let a33132 = Vel2::new(nx, dx, fx, nz, dz, fz, &a3313xz);
    /* last x and z in velocity model */
    lx = (fx + (((nx - 1i32) as f32) * dx));
    lz = (fz + (((nz - 1i32) as f32) * dz));
    /* ensure takeoff point is within model */
    if ((((x0 < fx) || (x0 > lx)) || (z0 < fz)) || (z0 > lz)) {
        crash();
    }
    // A negative allocation count returns NULL in C, then rs[0] is written.
    if nt < 0 {
        crash();
    }
    let mut rs = vec![RayStep::default(); nt.max(1) as usize];
    /* cosine and sine of takeoff angle */
    c = ((a0 as f64).cos() as f32);
    s = ((a0 as f64).sin() as f32);
    cc = (c * c);
    ss = (s * s);
    /* velocities and derivatives at takeoff point */
    (
        a1111,
        da1111dx,
        da1111dz,
        dda1111dxdx,
        dda1111dxdz,
        dda1111dzdz,
    ) = a11112.interp(x0, z0);
    da1111dn = ((da1111dx * c) - (da1111dz * s));
    dda1111dndn = (((((dda1111dxdx * cc) as f64)
        - (((2f64 * (dda1111dxdz as f64)) * (s as f64)) * (c as f64)))
        + ((dda1111dzdz * ss) as f64)) as f32);
    (
        a3333,
        da3333dx,
        da3333dz,
        dda3333dxdx,
        dda3333dxdz,
        dda3333dzdz,
    ) = a33332.interp(x0, z0);
    da3333dn = ((da3333dx * c) - (da3333dz * s));
    dda3333dndn = (((((dda3333dxdx * cc) as f64)
        - (((2f64 * (dda3333dxdz as f64)) * (s as f64)) * (c as f64)))
        + ((dda3333dzdz * ss) as f64)) as f32);
    (
        a1133,
        da1133dx,
        da1133dz,
        dda1133dxdx,
        dda1133dxdz,
        dda1133dzdz,
    ) = a11332.interp(x0, z0);
    da1133dn = ((da1133dx * c) - (da1133dz * s));
    dda1133dndn = (((((dda1133dxdx * cc) as f64)
        - (((2f64 * (dda1133dxdz as f64)) * (s as f64)) * (c as f64)))
        + ((dda1133dzdz * ss) as f64)) as f32);
    (
        a1313,
        da1313dx,
        da1313dz,
        dda1313dxdx,
        dda1313dxdz,
        dda1313dzdz,
    ) = a13132.interp(x0, z0);
    da1313dn = ((da1313dx * c) - (da1313dz * s));
    dda1313dndn = (((((dda1313dxdx * cc) as f64)
        - (((2f64 * (dda1313dxdz as f64)) * (s as f64)) * (c as f64)))
        + ((dda1313dzdz * ss) as f64)) as f32);
    (
        a1113,
        da1113dx,
        da1113dz,
        dda1113dxdx,
        dda1113dxdz,
        dda1113dzdz,
    ) = a11132.interp(x0, z0);
    da1113dn = ((da1113dx * c) - (da1113dz * s));
    dda1113dndn = (((((dda1113dxdx * cc) as f64)
        - (((2f64 * (dda1113dxdz as f64)) * (s as f64)) * (c as f64)))
        + ((dda1113dzdz * ss) as f64)) as f32);
    (
        a3313,
        da3313dx,
        da3313dz,
        dda3313dxdx,
        dda3313dxdz,
        dda3313dzdz,
    ) = a33132.interp(x0, z0);
    da3313dn = ((da3313dx * c) - (da3313dz * s));
    dda3313dndn = (((((dda3313dxdx * cc) as f64)
        - (((2f64 * (dda3313dxdz as f64)) * (s as f64)) * (c as f64)))
        + ((dda3313dzdz * ss) as f64)) as f32);
    /*computing the phase velocity for a0 angle */
    gamm11 = (((a1111 * ss) + (a1313 * cc)) + ((((2i32 as f32) * a1113) * s) * c));
    gamm33 = (((a3333 * cc) + (a1313 * ss)) + ((((2i32 as f32) * a3313) * s) * c));
    gamm13 = (((((a1133 + a1313) * s) * c) + (a1113 * ss)) + (a3313 * cc));
    sqr = (((((gamm11 + gamm33) * (gamm11 + gamm33))
        - ((4i32 as f32) * ((gamm11 * gamm33) - (gamm13 * gamm13)))) as f64)
        .sqrt() as f32);
    if (SV != 0) {
        vp2 = ((gamm11 + gamm33) - sqr);
    } else {
        vp2 = ((gamm11 + gamm33) + sqr);
    }
    vp = (((vp2 as f64) * 0.5f64).sqrt() as f32);
    ovp = ((1i32 as f32) / vp);
    px = (s * ovp);
    pz = (c * ovp);
    /* first ray step */
    rs[(0i32) as usize].t = {
        t = ft;
        t
    };
    rs[(0i32) as usize].x = {
        x = x0;
        x
    };
    rs[(0i32) as usize].z = {
        z = z0;
        z
    };
    rs[(0i32) as usize].q1 = {
        q1 = (1f64 as f32);
        q1
    };
    rs[(0i32) as usize].p1 = {
        p1 = (0f64 as f32);
        p1
    };
    rs[(0i32) as usize].q2 = {
        q2 = (0f64 as f32);
        q2
    };
    rs[(0i32) as usize].p2 = {
        p2 = (1f64 as f32);
        p2
    };
    rs[(0i32) as usize].kmah = {
        kmah = 0i32;
        kmah
    };
    rs[(0i32) as usize].c = c;
    rs[(0i32) as usize].s = s;
    rs[(0i32) as usize].v = vp;
    rs[(0i32) as usize].dvdx =
        ((((0.5f64 * (da3333dx as f64)) * (vp as f64)) / (a3333 as f64)) as f32);
    rs[(0i32) as usize].dvdz =
        ((((0.5f64 * (da3333dz as f64)) * (vp as f64)) / (a3333 as f64)) as f32);
    /* loop over time steps */
    it = 1i32;
    while (it < nt) {
        /* variables used for Runge-Kutta integration */
        let mut h: f32 = dt;
        let mut hhalf: f32 = (((dt as f64) / 2f64) as f32);
        let mut hsixth: f32 = (((dt as f64) / 6f64) as f32);
        let mut q2old: f32 = 0.0;
        let mut xt: f32 = 0.0;
        let mut zt: f32 = 0.0;
        let mut p1t: f32 = 0.0;
        let mut q1t: f32 = 0.0;
        let mut p2t: f32 = 0.0;
        let mut q2t: f32 = 0.0;
        let mut dx: f32 = 0.0;
        let mut dz: f32 = 0.0;
        let mut dp1: f32 = 0.0;
        let mut dq1: f32 = 0.0;
        let mut dp2: f32 = 0.0;
        let mut dq2: f32 = 0.0;
        let mut dxt: f32 = 0.0;
        let mut dzt: f32 = 0.0;
        let mut dp1t: f32 = 0.0;
        let mut dq1t: f32 = 0.0;
        let mut dp2t: f32 = 0.0;
        let mut dq2t: f32 = 0.0;
        let mut dxm: f32 = 0.0;
        let mut dzm: f32 = 0.0;
        let mut dp1m: f32 = 0.0;
        let mut dq1m: f32 = 0.0;
        let mut dp2m: f32 = 0.0;
        let mut dq2m: f32 = 0.0;
        let mut gamma11: f32 = 0.0;
        let mut gamma33: f32 = 0.0;
        let mut gamma13: f32 = 0.0;
        let mut g11: f32 = 0.0;
        let mut g13: f32 = 0.0;
        let mut g33: f32 = 0.0;
        let mut den: f32 = 0.0;
        let mut sxfactor: f32 = 0.0;
        let mut szfactor: f32 = 0.0;
        let mut snfact: f32 = 0.0;
        let mut dpx: f32 = 0.0;
        let mut dpz: f32 = 0.0;
        let mut pxt: f32 = 0.0;
        let mut pzt: f32 = 0.0;
        let mut dpxt: f32 = 0.0;
        let mut dpzt: f32 = 0.0;
        let mut dpxm: f32 = 0.0;
        let mut dpzm: f32 = 0.0;
        let mut dxdn: f32 = 0.0;
        let mut dzdn: f32 = 0.0;
        let mut snfactor: f32 = 0.0;
        let mut dxx: f32 = 0.0;
        let mut dzz: f32 = 0.0;
        let mut dcdp1: f32 = 0.0;
        let mut dcdp3: f32 = 0.0;
        let mut dcdp13: f32 = 0.0;
        let mut ddcdnn: f32 = 0.0;
        let mut ddcdqq: f32 = 0.0;
        let mut ddcdpn: f32 = 0.0;
        let mut dgamma11dpx: f32 = 0.0;
        let mut dgamma11dpz: f32 = 0.0;
        let mut dgamma33dpx: f32 = 0.0;
        let mut dgamma33dpz: f32 = 0.0;
        let mut dgamma13dpx: f32 = 0.0;
        let mut dgamma13dpz: f32 = 0.0;
        let mut dg11dpx: f32 = 0.0;
        let mut dg11dpz: f32 = 0.0;
        let mut dg33dpx: f32 = 0.0;
        let mut dg33dpz: f32 = 0.0;
        let mut dg13dpx: f32 = 0.0;
        let mut dg13dpz: f32 = 0.0;
        let mut ddxdpx: f32 = 0.0;
        let mut ddzdpz: f32 = 0.0;
        let mut ddxdpz: f32 = 0.0;
        let mut dgamma11dn: f32 = 0.0;
        let mut dgamma33dn: f32 = 0.0;
        let mut dgamma13dn: f32 = 0.0;
        let mut dg11dn: f32 = 0.0;
        let mut dg33dn: f32 = 0.0;
        let mut dg13dn: f32 = 0.0;
        let mut dsnfactordn: f32 = 0.0;
        let mut ddxdn: f32 = 0.0;
        let mut ddzdn: f32 = 0.0;
        /* if ray is out of bounds, break */
        if ((((((x < fx) || (x > lx)) || (z < fz)) || (z > lz))
            || ((c as f64) > ((amin as f64).cos() + 0.01f64)))
            || ((c as f64) < (((amax as f64).cos()) - 0.01f64)))
        {
            break;
        }
        /* remember old q2 */
        q2old = q2;
        /* step 1 of 4th-order Runge-Kutta */
        px2 = (px * px);
        pz2 = (pz * pz);
        pxz = (px * pz);
        /*anisotropy parameters*/
        gamma11 = (((a1111 * px2) + (a1313 * pz2)) + (((2i32 as f32) * a1113) * pxz));
        gamma33 = (((a3333 * pz2) + (a1313 * px2)) + (((2i32 as f32) * a3313) * pxz));
        gamma13 = ((((a1133 + a1313) * pxz) + (a1113 * px2)) + (a3313 * pz2));
        den = ((1i32 as f32) / ((gamma11 + gamma33) - (2i32 as f32)));
        g11 = ((gamma33 - (1i32 as f32)) * den);
        g33 = ((gamma11 - (1i32 as f32)) * den);
        g13 = ((-gamma13) * den);
        sxfactor = (((((((da1111dx * px2) * g11) + ((da3333dx * pz2) * g33))
            + ((((2i32 as f32) * (da1133dx + da1313dx)) * pxz) * g13))
            + (da1313dx * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * da3313dx) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * da1113dx) * ((pxz * g11) + (px2 * g13))));
        szfactor = (((((((da1111dz * px2) * g11) + ((da3333dz * pz2) * g33))
            + ((((2i32 as f32) * (da1133dz + da1313dz)) * pxz) * g13))
            + (da1313dz * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * da3313dz) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * da1113dz) * ((pxz * g11) + (px2 * g13))));
        snfact = ((sxfactor * c) - (szfactor * s));
        /*computing ray velocities and derivatives*/
        dx = ((((((a1111 * px) * g11) + (((a1133 + a1313) * pz) * g13)) + ((a3313 * pz) * g33))
            + (a1113 * ((pz * g11) + (((2i32 as f32) * px) * g13))))
            + ((a1313 * g33) * px));
        dz = ((((((a3333 * pz) * g33) + (((a1133 + a1313) * px) * g13)) + ((a1113 * px) * g11))
            + (a3313 * ((px * g33) + (((2i32 as f32) * pz) * g13))))
            + ((a1313 * g11) * pz));
        dgamma11dpx = ((((2i32 as f32) * a1111) * px) + (((2i32 as f32) * a1113) * pz));
        dgamma11dpz = ((((2i32 as f32) * a1313) * pz) + (((2i32 as f32) * a1113) * px));
        dgamma33dpx = ((((2i32 as f32) * a1313) * px) + (((2i32 as f32) * a3313) * pz));
        dgamma33dpz = ((((2i32 as f32) * a3333) * pz) + (((2i32 as f32) * a3313) * px));
        dgamma13dpx = (((a1133 + a1313) * pz) + (((2i32 as f32) * a1113) * px));
        dgamma13dpz = (((a1133 + a1313) * px) + (((2i32 as f32) * a3313) * pz));
        dgamma11dn = (((da1111dn * px2) + (da1313dn * pz2)) + (((2i32 as f32) * da1113dn) * pxz));
        dgamma33dn = (((da3333dn * pz2) + (da1313dn * px2)) + (((2i32 as f32) * da3313dn) * pxz));
        dgamma13dn = ((((da1133dn + da1313dn) * pxz) + (da1113dn * px2)) + (da3313dn * pz2));
        dg11dpx = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dx)))
            * den)
            * den)
            + ((dgamma33dpx - ((2i32 as f32) * dx)) * den));
        dg11dpz = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dz)))
            * den)
            * den)
            + ((dgamma33dpz - ((2i32 as f32) * dz)) * den));
        dg33dpx = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dx)))
            * den)
            * den)
            + ((dgamma11dpx - ((2i32 as f32) * dx)) * den));
        dg33dpz = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dz)))
            * den)
            * den)
            + ((dgamma11dpz - ((2i32 as f32) * dz)) * den));
        dg13dpx = ((((gamma13 * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dx))) * den)
            * den)
            - (dgamma13dpx * den));
        dg13dpz = ((((gamma13 * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dz))) * den)
            * den)
            - (dgamma13dpz * den));
        dg11dn = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact)))
            * den)
            * den)
            + ((dgamma33dn - snfact) * den));
        dg33dn = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact)))
            * den)
            * den)
            + ((dgamma11dn - snfact) * den));
        dg13dn = ((((gamma13 * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact))) * den)
            * den)
            - (dgamma13dn * den));
        ddxdpx = ((((((a1111 * px) * dg11dpx) + (((a1133 + a1313) * pz) * dg13dpx))
            + ((a3313 * pz) * dg33dpx))
            + (a1113 * ((pz * dg11dpx) + (((2i32 as f32) * px) * dg13dpx))))
            + ((a1313 * dg33dpx) * px));
        ddzdpz = ((((((a3333 * pz) * dg33dpz) + (((a1133 + a1313) * px) * dg13dpz))
            + ((a1113 * px) * dg11dpz))
            + (a3313 * ((px * dg33dpz) + (((2i32 as f32) * pz) * dg13dpz))))
            + ((a1313 * dg11dpz) * pz));
        ddxdpz = ((((((a1111 * px) * dg11dpz) + (((a1133 + a1313) * pz) * dg13dpz))
            + ((a3313 * pz) * dg33dpz))
            + (a1113 * ((pz * dg11dpz) + (((2i32 as f32) * px) * dg13dpz))))
            + ((a1313 * dg33dpz) * px));
        dsnfactordn = (((((((da1111dn * px2) * dg11dn) + ((da3333dn * pz2) * dg33dn))
            + ((((2i32 as f32) * (da1133dn + da1313dn)) * pxz) * dg13dn))
            + (da1313dn * ((px2 * dg33dn) + (pz2 * dg11dn))))
            + (((2i32 as f32) * da3313dn) * ((pz2 * dg13dn) + (pxz * dg33dn))))
            + (((2i32 as f32) * da1113dn) * ((pxz * dg11dn) + (px2 * dg13dn))));
        ddxdn = ((((((a1111 * px) * dg11dn) + (((a1133 + a1313) * pz) * dg13dn))
            + ((a3313 * pz) * dg33dn))
            + (a1113 * ((pz * dg11dn) + (((2i32 as f32) * px) * dg13dn))))
            + ((a1313 * dg33dn) * px));
        ddzdn = ((((((a3333 * pz) * dg33dn) + (((a1133 + a1313) * px) * dg13dn))
            + ((a1113 * px) * dg11dn))
            + (a3313 * ((px * dg33dn) + (((2i32 as f32) * pz) * dg13dn))))
            + ((a1313 * dg11dn) * pz));
        /*evaluating change in slowness and amplitude along ray*/
        dpx = (((-0.5f64) * (sxfactor as f64)) as f32);
        dpz = (((-0.5f64) * (szfactor as f64)) as f32);
        dcdp1 = (((((a1111 * g11) + (a1313 * g33)) + (((2i32 as f32) * a1113) * g13)) + ddxdpx)
            - (dx * dx));
        dcdp3 = (((((a3333 * g33) + (a1313 * g11)) + (((2i32 as f32) * a3313) * g13)) + ddzdpz)
            - (dz * dz));
        dcdp13 = ((((((a1133 * g13) + (a1313 * g13)) + (a1113 * g11)) + (a3313 * g33)) + ddxdpz)
            - (dx * dz));
        ddcdqq = (((((dcdp1 * cc) as f64) - (((2f64 * (dcdp13 as f64)) * (s as f64)) * (c as f64)))
            + ((dcdp3 * ss) as f64)) as f32);
        dxdn = ((((((da1111dn * px) * g11) + (((da1133dn + da1313dn) * pz) * g13))
            + ((da3313dn * pz) * g33))
            + (da1113dn * ((pz * g11) + (((2i32 as f32) * px) * g13))))
            + ((da1313dn * g33) * px));
        dzdn = ((((((da3333dn * pz) * g33) + (((da1133dn + da1313dn) * px) * g13))
            + ((da1113dn * px) * g11))
            + (da3313dn * ((px * g33) + (((2i32 as f32) * pz) * g13))))
            + ((da1313dn * g11) * pz));
        ddcdpn = (((((((((dxdn * c) - (dzdn * s)) as f64)
            - (((0.5f64 * (dx as f64)) * (sxfactor as f64)) * (cc as f64)))
            + (((0.5f64 * (((dx * szfactor) + (dz * sxfactor)) as f64)) * (s as f64))
                * (c as f64)))
            - (((0.5f64 * (dz as f64)) * (szfactor as f64)) * (ss as f64)))
            + ((ddxdn * c) as f64))
            - ((ddzdn * s) as f64)) as f32);
        snfactor = (((((((dda1111dndn * px2) * g11) + ((dda3333dndn * pz2) * g33))
            + ((((2i32 as f32) * (dda1133dndn + dda1313dndn)) * pxz) * g13))
            + (dda1313dndn * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * dda3313dndn) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * dda1113dndn) * ((pxz * g11) + (px2 * g13))));
        ddcdnn = ((((((0.5f64 * (snfactor as f64))
            - (((0.25f64 * (sxfactor as f64)) * (sxfactor as f64)) * (cc as f64)))
            + ((((0.5f64 * (sxfactor as f64)) * (szfactor as f64)) * (s as f64)) * (c as f64)))
            - (((0.25f64 * (szfactor as f64)) * (szfactor as f64)) * (ss as f64)))
            + (0.5f64 * (dsnfactordn as f64))) as f32);
        dp1 = (((-ddcdnn) * q1) - (ddcdpn * p1));
        dq1 = ((ddcdqq * p1) + (ddcdpn * q1));
        dp2 = (((-ddcdnn) * q2) - (ddcdpn * p2));
        dq2 = ((ddcdqq * p2) + (ddcdpn * q2));
        xt = (x + (hhalf * dx));
        zt = (z + (hhalf * dz));
        pxt = (px + (hhalf * dpx));
        pzt = (pz + (hhalf * dpz));
        p1t = (p1 + (hhalf * dp1));
        q1t = (q1 + (hhalf * dq1));
        p2t = (p2 + (hhalf * dp2));
        q2t = (q2 + (hhalf * dq2));
        vp = (((1i32 as f64) / (((pxt * pxt) + (pzt * pzt)) as f64).sqrt()) as f32);
        s = (pxt * vp);
        c = (pzt * vp);
        ss = (s * s);
        cc = (c * c);
        (
            a1111,
            da1111dx,
            da1111dz,
            dda1111dxdx,
            dda1111dxdz,
            dda1111dzdz,
        ) = a11112.interp(xt, zt);
        da1111dn = ((da1111dx * c) - (da1111dz * s));
        dda1111dndn = (((((dda1111dxdx * cc) as f64)
            - (((2f64 * (dda1111dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1111dzdz * ss) as f64)) as f32);
        (
            a3333,
            da3333dx,
            da3333dz,
            dda3333dxdx,
            dda3333dxdz,
            dda3333dzdz,
        ) = a33332.interp(xt, zt);
        da3333dn = ((da3333dx * c) - (da3333dz * s));
        dda3333dndn = (((((dda3333dxdx * cc) as f64)
            - (((2f64 * (dda3333dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda3333dzdz * ss) as f64)) as f32);
        (
            a1133,
            da1133dx,
            da1133dz,
            dda1133dxdx,
            dda1133dxdz,
            dda1133dzdz,
        ) = a11332.interp(xt, zt);
        da1133dn = ((da1133dx * c) - (da1133dz * s));
        dda1133dndn = (((((dda1133dxdx * cc) as f64)
            - (((2f64 * (dda1133dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1133dzdz * ss) as f64)) as f32);
        (
            a1313,
            da1313dx,
            da1313dz,
            dda1313dxdx,
            dda1313dxdz,
            dda1313dzdz,
        ) = a13132.interp(xt, zt);
        da1313dn = ((da1313dx * c) - (da1313dz * s));
        dda1313dndn = (((((dda1313dxdx * cc) as f64)
            - (((2f64 * (dda1313dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1313dzdz * ss) as f64)) as f32);
        (
            a1113,
            da1113dx,
            da1113dz,
            dda1113dxdx,
            dda1113dxdz,
            dda1113dzdz,
        ) = a11132.interp(xt, zt);
        da1113dn = ((da1113dx * c) - (da1113dz * s));
        dda1113dndn = (((((dda1113dxdx * cc) as f64)
            - (((2f64 * (dda1113dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1113dzdz * ss) as f64)) as f32);
        (
            a3313,
            da3313dx,
            da3313dz,
            dda3313dxdx,
            dda3313dxdz,
            dda3313dzdz,
        ) = a33132.interp(xt, zt);
        da3313dn = ((da3313dx * c) - (da3313dz * s));
        dda3313dndn = (((((dda3313dxdx * cc) as f64)
            - (((2f64 * (dda3313dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda3313dzdz * ss) as f64)) as f32);
        /* step 2 of 4th-order Runge-Kutta */
        px2 = (pxt * pxt);
        pz2 = (pzt * pzt);
        pxz = (pxt * pzt);
        /*anisotropy parameters*/
        gamma11 = (((a1111 * px2) + (a1313 * pz2)) + (((2i32 as f32) * a1113) * pxz));
        gamma33 = (((a3333 * pz2) + (a1313 * px2)) + (((2i32 as f32) * a3313) * pxz));
        gamma13 = ((((a1133 + a1313) * pxz) + (a1113 * px2)) + (a3313 * pz2));
        den = ((1i32 as f32) / ((gamma11 + gamma33) - (2i32 as f32)));
        g11 = ((gamma33 - (1i32 as f32)) * den);
        g33 = ((gamma11 - (1i32 as f32)) * den);
        g13 = ((-gamma13) * den);
        sxfactor = (((((((da1111dx * px2) * g11) + ((da3333dx * pz2) * g33))
            + ((((2i32 as f32) * (da1133dx + da1313dx)) * pxz) * g13))
            + (da1313dx * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * da3313dx) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * da1113dx) * ((pxz * g11) + (px2 * g13))));
        szfactor = (((((((da1111dz * px2) * g11) + ((da3333dz * pz2) * g33))
            + ((((2i32 as f32) * (da1133dz + da1313dz)) * pxz) * g13))
            + (da1313dz * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * da3313dz) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * da1113dz) * ((pxz * g11) + (px2 * g13))));
        snfact = ((sxfactor * c) - (szfactor * s));
        /*computing ray velocities and derivatives*/
        dxt = ((((((a1111 * pxt) * g11) + (((a1133 + a1313) * pzt) * g13))
            + ((a3313 * pzt) * g33))
            + (a1113 * ((pzt * g11) + (((2i32 as f32) * pxt) * g13))))
            + ((a1313 * g33) * pxt));
        dzt = ((((((a3333 * pzt) * g33) + (((a1133 + a1313) * pxt) * g13))
            + ((a1113 * pxt) * g11))
            + (a3313 * ((pxt * g33) + (((2i32 as f32) * pzt) * g13))))
            + ((a1313 * g11) * pzt));
        dpxt = (((-0.5f64) * (sxfactor as f64)) as f32);
        dpzt = (((-0.5f64) * (szfactor as f64)) as f32);
        dgamma11dpx = ((((2i32 as f32) * a1111) * pxt) + (((2i32 as f32) * a1113) * pzt));
        dgamma11dpz = ((((2i32 as f32) * a1313) * pzt) + (((2i32 as f32) * a1113) * pxt));
        dgamma33dpx = ((((2i32 as f32) * a1313) * pxt) + (((2i32 as f32) * a3313) * pzt));
        dgamma33dpz = ((((2i32 as f32) * a3333) * pzt) + (((2i32 as f32) * a3313) * pxt));
        dgamma13dpx = (((a1133 + a1313) * pzt) + (((2i32 as f32) * a1113) * pxt));
        dgamma13dpz = (((a1133 + a1313) * pxt) + (((2i32 as f32) * a3313) * pzt));
        dgamma11dn = (((da1111dn * px2) + (da1313dn * pz2)) + (((2i32 as f32) * da1113dn) * pxz));
        dgamma33dn = (((da3333dn * pz2) + (da1313dn * px2)) + (((2i32 as f32) * da3313dn) * pxz));
        dgamma13dn = ((((da1133dn + da1313dn) * pxz) + (da1113dn * px2)) + (da3313dn * pz2));
        dg11dpx = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxt)))
            * den)
            * den)
            + ((dgamma33dpx - ((2i32 as f32) * dxt)) * den));
        dg11dpz = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzt)))
            * den)
            * den)
            + ((dgamma33dpz - ((2i32 as f32) * dzt)) * den));
        dg33dpx = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxt)))
            * den)
            * den)
            + ((dgamma11dpx - ((2i32 as f32) * dxt)) * den));
        dg33dpz = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzt)))
            * den)
            * den)
            + ((dgamma11dpz - ((2i32 as f32) * dzt)) * den));
        dg13dpx = ((((gamma13 * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxt))) * den)
            * den)
            - (dgamma13dpx * den));
        dg13dpz = ((((gamma13 * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzt))) * den)
            * den)
            - (dgamma13dpz * den));
        dg11dn = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact)))
            * den)
            * den)
            + ((dgamma33dn - snfact) * den));
        dg33dn = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact)))
            * den)
            * den)
            + ((dgamma11dn - snfact) * den));
        dg13dn = ((((gamma13 * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact))) * den)
            * den)
            - (dgamma13dn * den));
        ddxdpx = ((((((a1111 * pxt) * dg11dpx) + (((a1133 + a1313) * pzt) * dg13dpx))
            + ((a3313 * pzt) * dg33dpx))
            + (a1113 * ((pzt * dg11dpx) + (((2i32 as f32) * pxt) * dg13dpx))))
            + ((a1313 * dg33dpx) * pxt));
        ddzdpz = ((((((a3333 * pzt) * dg33dpz) + (((a1133 + a1313) * pxt) * dg13dpz))
            + ((a1113 * pxt) * dg11dpz))
            + (a3313 * ((pxt * dg33dpz) + (((2i32 as f32) * pzt) * dg13dpz))))
            + ((a1313 * dg11dpz) * pzt));
        ddxdpz = ((((((a1111 * pxt) * dg11dpz) + (((a1133 + a1313) * pzt) * dg13dpz))
            + ((a3313 * pzt) * dg33dpz))
            + (a1113 * ((pzt * dg11dpz) + (((2i32 as f32) * pxt) * dg13dpz))))
            + ((a1313 * dg33dpz) * pxt));
        dsnfactordn = (((((((da1111dn * px2) * dg11dn) + ((da3333dn * pz2) * dg33dn))
            + ((((2i32 as f32) * (da1133dn + da1313dn)) * pxz) * dg13dn))
            + (da1313dn * ((px2 * dg33dn) + (pz2 * dg11dn))))
            + (((2i32 as f32) * da3313dn) * ((pz2 * dg13dn) + (pxz * dg33dn))))
            + (((2i32 as f32) * da1113dn) * ((pxz * dg11dn) + (px2 * dg13dn))));
        ddxdn = ((((((a1111 * pxt) * dg11dn) + (((a1133 + a1313) * pzt) * dg13dn))
            + ((a3313 * pzt) * dg33dn))
            + (a1113 * ((pzt * dg11dn) + (((2i32 as f32) * pxt) * dg13dn))))
            + ((a1313 * dg33dn) * pxt));
        ddzdn = ((((((a3333 * pzt) * dg33dn) + (((a1133 + a1313) * pxt) * dg13dn))
            + ((a1113 * pxt) * dg11dn))
            + (a3313 * ((pxt * dg33dn) + (((2i32 as f32) * pzt) * dg13dn))))
            + ((a1313 * dg11dn) * pzt));
        dcdp1 = (((((a1111 * g11) + (a1313 * g33)) + (((2i32 as f32) * a1113) * g13)) + ddxdpx)
            - (dxt * dxt));
        dcdp3 = (((((a3333 * g33) + (a1313 * g11)) + (((2i32 as f32) * a3313) * g13)) + ddzdpz)
            - (dzt * dzt));
        dcdp13 = ((((((a1133 * g13) + (a1313 * g13)) + (a1113 * g11)) + (a3313 * g33)) + ddxdpz)
            - (dxt * dzt));
        ddcdqq = (((((dcdp1 * cc) as f64) - (((2f64 * (dcdp13 as f64)) * (s as f64)) * (c as f64)))
            + ((dcdp3 * ss) as f64)) as f32);
        dxdn = ((((((da1111dn * pxt) * g11) + (((da1133dn + da1313dn) * pzt) * g13))
            + ((da3313dn * pzt) * g33))
            + (da1113dn * ((pzt * g11) + (((2i32 as f32) * pxt) * g13))))
            + ((da1313dn * g33) * pxt));
        dzdn = ((((((da3333dn * pzt) * g33) + (((da1133dn + da1313dn) * pxt) * g13))
            + ((da1113dn * pxt) * g11))
            + (da3313dn * ((pxt * g33) + (((2i32 as f32) * pzt) * g13))))
            + ((da1313dn * g11) * pzt));
        ddcdpn = (((((((((dxdn * c) - (dzdn * s)) as f64)
            - (((0.5f64 * (dxt as f64)) * (sxfactor as f64)) * (cc as f64)))
            + (((0.5f64 * (((dxt * szfactor) + (dzt * sxfactor)) as f64)) * (s as f64))
                * (c as f64)))
            - (((0.5f64 * (dzt as f64)) * (szfactor as f64)) * (ss as f64)))
            + ((ddxdn * c) as f64))
            - ((ddzdn * s) as f64)) as f32);
        snfactor = (((((((dda1111dndn * px2) * g11) + ((dda3333dndn * pz2) * g33))
            + ((((2i32 as f32) * (dda1133dndn + dda1313dndn)) * pxz) * g13))
            + (dda1313dndn * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * dda3313dndn) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * dda1113dndn) * ((pxz * g11) + (px2 * g13))));
        ddcdnn = ((((((0.5f64 * (snfactor as f64))
            - (((0.25f64 * (sxfactor as f64)) * (sxfactor as f64)) * (cc as f64)))
            + ((((0.5f64 * (sxfactor as f64)) * (szfactor as f64)) * (s as f64)) * (c as f64)))
            - (((0.25f64 * (szfactor as f64)) * (szfactor as f64)) * (ss as f64)))
            + (0.5f64 * (dsnfactordn as f64))) as f32);
        dp1t = (((-ddcdnn) * q1t) - (ddcdpn * p1t));
        dq1t = ((ddcdqq * p1t) + (ddcdpn * q1t));
        dp2t = (((-ddcdnn) * q2t) - (ddcdpn * p2t));
        dq2t = ((ddcdqq * p2t) + (ddcdpn * q2t));
        xt = (x + (hhalf * dxt));
        zt = (z + (hhalf * dzt));
        pxt = (px + (hhalf * dpxt));
        pzt = (pz + (hhalf * dpzt));
        p1t = (p1 + (hhalf * dp1t));
        q1t = (q1 + (hhalf * dq1t));
        p2t = (p2 + (hhalf * dp2t));
        q2t = (q2 + (hhalf * dq2t));
        vp = (((1i32 as f64) / (((pxt * pxt) + (pzt * pzt)) as f64).sqrt()) as f32);
        s = (pxt * vp);
        c = (pzt * vp);
        ss = (s * s);
        cc = (c * c);
        (
            a1111,
            da1111dx,
            da1111dz,
            dda1111dxdx,
            dda1111dxdz,
            dda1111dzdz,
        ) = a11112.interp(xt, zt);
        da1111dn = ((da1111dx * c) - (da1111dz * s));
        dda1111dndn = (((((dda1111dxdx * cc) as f64)
            - (((2f64 * (dda1111dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1111dzdz * ss) as f64)) as f32);
        (
            a3333,
            da3333dx,
            da3333dz,
            dda3333dxdx,
            dda3333dxdz,
            dda3333dzdz,
        ) = a33332.interp(xt, zt);
        da3333dn = ((da3333dx * c) - (da3333dz * s));
        dda3333dndn = (((((dda3333dxdx * cc) as f64)
            - (((2f64 * (dda3333dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda3333dzdz * ss) as f64)) as f32);
        (
            a1133,
            da1133dx,
            da1133dz,
            dda1133dxdx,
            dda1133dxdz,
            dda1133dzdz,
        ) = a11332.interp(xt, zt);
        da1133dn = ((da1133dx * c) - (da1133dz * s));
        dda1133dndn = (((((dda1133dxdx * cc) as f64)
            - (((2f64 * (dda1133dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1133dzdz * ss) as f64)) as f32);
        (
            a1313,
            da1313dx,
            da1313dz,
            dda1313dxdx,
            dda1313dxdz,
            dda1313dzdz,
        ) = a13132.interp(xt, zt);
        da1313dn = ((da1313dx * c) - (da1313dz * s));
        dda1313dndn = (((((dda1313dxdx * cc) as f64)
            - (((2f64 * (dda1313dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1313dzdz * ss) as f64)) as f32);
        (
            a1113,
            da1113dx,
            da1113dz,
            dda1113dxdx,
            dda1113dxdz,
            dda1113dzdz,
        ) = a11132.interp(xt, zt);
        da1113dn = ((da1113dx * c) - (da1113dz * s));
        dda1113dndn = (((((dda1113dxdx * cc) as f64)
            - (((2f64 * (dda1113dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1113dzdz * ss) as f64)) as f32);
        (
            a3313,
            da3313dx,
            da3313dz,
            dda3313dxdx,
            dda3313dxdz,
            dda3313dzdz,
        ) = a33132.interp(xt, zt);
        da3313dn = ((da3313dx * c) - (da3313dz * s));
        dda3313dndn = (((((dda3313dxdx * cc) as f64)
            - (((2f64 * (dda3313dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda3313dzdz * ss) as f64)) as f32);
        /* step 3 of 4th-order Runge-Kutta */
        px2 = (pxt * pxt);
        pz2 = (pzt * pzt);
        pxz = (pxt * pzt);
        /*anisotropy parameters*/
        gamma11 = (((a1111 * px2) + (a1313 * pz2)) + (((2i32 as f32) * a1113) * pxz));
        gamma33 = (((a3333 * pz2) + (a1313 * px2)) + (((2i32 as f32) * a3313) * pxz));
        gamma13 = ((((a1133 + a1313) * pxz) + (a1113 * px2)) + (a3313 * pz2));
        den = ((1i32 as f32) / ((gamma11 + gamma33) - (2i32 as f32)));
        g11 = ((gamma33 - (1i32 as f32)) * den);
        g33 = ((gamma11 - (1i32 as f32)) * den);
        g13 = ((-gamma13) * den);
        sxfactor = (((((((da1111dx * px2) * g11) + ((da3333dx * pz2) * g33))
            + ((((2i32 as f32) * (da1133dx + da1313dx)) * pxz) * g13))
            + (da1313dx * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * da3313dx) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * da1113dx) * ((pxz * g11) + (px2 * g13))));
        szfactor = (((((((da1111dz * px2) * g11) + ((da3333dz * pz2) * g33))
            + ((((2i32 as f32) * (da1133dz + da1313dz)) * pxz) * g13))
            + (da1313dz * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * da3313dz) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * da1113dz) * ((pxz * g11) + (px2 * g13))));
        snfact = ((sxfactor * c) - (szfactor * s));
        /*computing ray velocities and derivatives*/
        dxm = ((((((a1111 * pxt) * g11) + (((a1133 + a1313) * pzt) * g13))
            + ((a3313 * pzt) * g33))
            + (a1113 * ((pzt * g11) + (((2i32 as f32) * pxt) * g13))))
            + ((a1313 * g33) * pxt));
        dzm = ((((((a3333 * pzt) * g33) + (((a1133 + a1313) * pxt) * g13))
            + ((a1113 * pxt) * g11))
            + (a3313 * ((pxt * g33) + (((2i32 as f32) * pzt) * g13))))
            + ((a1313 * g11) * pzt));
        dpxm = (((-0.5f64) * (sxfactor as f64)) as f32);
        dpzm = (((-0.5f64) * (szfactor as f64)) as f32);
        dgamma11dpx = ((((2i32 as f32) * a1111) * pxt) + (((2i32 as f32) * a1113) * pzt));
        dgamma11dpz = ((((2i32 as f32) * a1313) * pzt) + (((2i32 as f32) * a1113) * pxt));
        dgamma33dpx = ((((2i32 as f32) * a1313) * pxt) + (((2i32 as f32) * a3313) * pzt));
        dgamma33dpz = ((((2i32 as f32) * a3333) * pzt) + (((2i32 as f32) * a3313) * pxt));
        dgamma13dpx = (((a1133 + a1313) * pzt) + (((2i32 as f32) * a1113) * pxt));
        dgamma13dpz = (((a1133 + a1313) * pxt) + (((2i32 as f32) * a3313) * pzt));
        dgamma11dn = (((da1111dn * px2) + (da1313dn * pz2)) + (((2i32 as f32) * da1113dn) * pxz));
        dgamma33dn = (((da3333dn * pz2) + (da1313dn * px2)) + (((2i32 as f32) * da3313dn) * pxz));
        dgamma13dn = ((((da1133dn + da1313dn) * pxz) + (da1113dn * px2)) + (da3313dn * pz2));
        dg11dpx = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxm)))
            * den)
            * den)
            + ((dgamma33dpx - ((2i32 as f32) * dxm)) * den));
        dg11dpz = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzm)))
            * den)
            * den)
            + ((dgamma33dpz - ((2i32 as f32) * dzm)) * den));
        dg33dpx = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxm)))
            * den)
            * den)
            + ((dgamma11dpx - ((2i32 as f32) * dxm)) * den));
        dg33dpz = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzm)))
            * den)
            * den)
            + ((dgamma11dpz - ((2i32 as f32) * dzm)) * den));
        dg13dpx = ((((gamma13 * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxm))) * den)
            * den)
            - (dgamma13dpx * den));
        dg13dpz = ((((gamma13 * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzm))) * den)
            * den)
            - (dgamma13dpz * den));
        dg11dn = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact)))
            * den)
            * den)
            + ((dgamma33dn - snfact) * den));
        dg33dn = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact)))
            * den)
            * den)
            + ((dgamma11dn - snfact) * den));
        dg13dn = ((((gamma13 * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact))) * den)
            * den)
            - (dgamma13dn * den));
        ddxdpx = ((((((a1111 * pxt) * dg11dpx) + (((a1133 + a1313) * pzt) * dg13dpx))
            + ((a3313 * pzt) * dg33dpx))
            + (a1113 * ((pzt * dg11dpx) + (((2i32 as f32) * pxt) * dg13dpx))))
            + ((a1313 * dg33dpx) * pxt));
        ddzdpz = ((((((a3333 * pzt) * dg33dpz) + (((a1133 + a1313) * pxt) * dg13dpz))
            + ((a1113 * pxt) * dg11dpz))
            + (a3313 * ((pxt * dg33dpz) + (((2i32 as f32) * pzt) * dg13dpz))))
            + ((a1313 * dg11dpz) * pzt));
        ddxdpz = ((((((a1111 * pxt) * dg11dpz) + (((a1133 + a1313) * pzt) * dg13dpz))
            + ((a3313 * pzt) * dg33dpz))
            + (a1113 * ((pzt * dg11dpz) + (((2i32 as f32) * pxt) * dg13dpz))))
            + ((a1313 * dg33dpz) * pxt));
        dsnfactordn = (((((((da1111dn * px2) * dg11dn) + ((da3333dn * pz2) * dg33dn))
            + ((((2i32 as f32) * (da1133dn + da1313dn)) * pxz) * dg13dn))
            + (da1313dn * ((px2 * dg33dn) + (pz2 * dg11dn))))
            + (((2i32 as f32) * da3313dn) * ((pz2 * dg13dn) + (pxz * dg33dn))))
            + (((2i32 as f32) * da1113dn) * ((pxz * dg11dn) + (px2 * dg13dn))));
        ddxdn = ((((((a1111 * pxt) * dg11dn) + (((a1133 + a1313) * pzt) * dg13dn))
            + ((a3313 * pzt) * dg33dn))
            + (a1113 * ((pzt * dg11dn) + (((2i32 as f32) * pxt) * dg13dn))))
            + ((a1313 * dg33dn) * pxt));
        ddzdn = ((((((a3333 * pzt) * dg33dn) + (((a1133 + a1313) * pxt) * dg13dn))
            + ((a1113 * pxt) * dg11dn))
            + (a3313 * ((pxt * dg33dn) + (((2i32 as f32) * pzt) * dg13dn))))
            + ((a1313 * dg11dn) * pzt));
        dcdp1 = (((((a1111 * g11) + (a1313 * g33)) + (((2i32 as f32) * a1113) * g13)) + ddxdpx)
            - (dxm * dxm));
        dcdp3 = (((((a3333 * g33) + (a1313 * g11)) + (((2i32 as f32) * a3313) * g13)) + ddzdpz)
            - (dzm * dzm));
        dcdp13 = ((((((a1133 * g13) + (a1313 * g13)) + (a1113 * g11)) + (a3313 * g33)) + ddxdpz)
            - (dxm * dzm));
        ddcdqq = (((((dcdp1 * cc) as f64) - (((2f64 * (dcdp13 as f64)) * (s as f64)) * (c as f64)))
            + ((dcdp3 * ss) as f64)) as f32);
        dxdn = ((((((da1111dn * pxt) * g11) + (((da1133dn + da1313dn) * pzt) * g13))
            + ((da3313dn * pzt) * g33))
            + (da1113dn * ((pzt * g11) + (((2i32 as f32) * pxt) * g13))))
            + ((da1313dn * g33) * pxt));
        dzdn = ((((((da3333dn * pzt) * g33) + (((da1133dn + da1313dn) * pxt) * g13))
            + ((da1113dn * pxt) * g11))
            + (da3313dn * ((pxt * g33) + (((2i32 as f32) * pzt) * g13))))
            + ((da1313dn * g11) * pzt));
        ddcdpn = (((((((((dxdn * c) - (dzdn * s)) as f64)
            - (((0.5f64 * (dxm as f64)) * (sxfactor as f64)) * (cc as f64)))
            + (((0.5f64 * (((dxm * szfactor) + (dzm * sxfactor)) as f64)) * (s as f64))
                * (c as f64)))
            - (((0.5f64 * (dzm as f64)) * (szfactor as f64)) * (ss as f64)))
            + ((ddxdn * c) as f64))
            - ((ddzdn * s) as f64)) as f32);
        snfactor = (((((((dda1111dndn * px2) * g11) + ((dda3333dndn * pz2) * g33))
            + ((((2i32 as f32) * (dda1133dndn + dda1313dndn)) * pxz) * g13))
            + (dda1313dndn * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * dda3313dndn) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * dda1113dndn) * ((pxz * g11) + (px2 * g13))));
        ddcdnn = ((((((0.5f64 * (snfactor as f64))
            - (((0.25f64 * (sxfactor as f64)) * (sxfactor as f64)) * (cc as f64)))
            + ((((0.5f64 * (sxfactor as f64)) * (szfactor as f64)) * (s as f64)) * (c as f64)))
            - (((0.25f64 * (szfactor as f64)) * (szfactor as f64)) * (ss as f64)))
            + (0.5f64 * (dsnfactordn as f64))) as f32);
        dp1m = (((-ddcdnn) * q1t) - (ddcdpn * p1t));
        dq1m = ((ddcdqq * p1t) + (ddcdpn * q1t));
        dp2m = (((-ddcdnn) * q2t) - (ddcdpn * p2t));
        dq2m = ((ddcdqq * p2t) + (ddcdpn * q2t));
        xt = (x + (hhalf * dx));
        zt = (z + (hhalf * dz));
        pxt = (px + (h * dpxm));
        pzt = (pz + (h * dpzm));
        p1t = (p1 + (h * dp1m));
        q1t = (q1 + (h * dq1m));
        p2t = (p2 + (h * dp2m));
        q2t = (q2 + (h * dq2m));
        dxm += dxt;
        dzm += dzt;
        dpxm += dpxt;
        dpzm += dpzt;
        dp1m += dp1t;
        dq1m += dq1t;
        dp2m += dp2t;
        dq2m += dq2t;
        vp = (((1i32 as f64) / (((pxt * pxt) + (pzt * pzt)) as f64).sqrt()) as f32);
        s = (pxt * vp);
        c = (pzt * vp);
        ss = (s * s);
        cc = (c * c);
        (
            a1111,
            da1111dx,
            da1111dz,
            dda1111dxdx,
            dda1111dxdz,
            dda1111dzdz,
        ) = a11112.interp(xt, zt);
        da1111dn = ((da1111dx * c) - (da1111dz * s));
        dda1111dndn = (((((dda1111dxdx * cc) as f64)
            - (((2f64 * (dda1111dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1111dzdz * ss) as f64)) as f32);
        (
            a3333,
            da3333dx,
            da3333dz,
            dda3333dxdx,
            dda3333dxdz,
            dda3333dzdz,
        ) = a33332.interp(xt, zt);
        da3333dn = ((da3333dx * c) - (da3333dz * s));
        dda3333dndn = (((((dda3333dxdx * cc) as f64)
            - (((2f64 * (dda3333dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda3333dzdz * ss) as f64)) as f32);
        (
            a1133,
            da1133dx,
            da1133dz,
            dda1133dxdx,
            dda1133dxdz,
            dda1133dzdz,
        ) = a11332.interp(xt, zt);
        da1133dn = ((da1133dx * c) - (da1133dz * s));
        dda1133dndn = (((((dda1133dxdx * cc) as f64)
            - (((2f64 * (dda1133dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1133dzdz * ss) as f64)) as f32);
        (
            a1313,
            da1313dx,
            da1313dz,
            dda1313dxdx,
            dda1313dxdz,
            dda1313dzdz,
        ) = a13132.interp(xt, zt);
        da1313dn = ((da1313dx * c) - (da1313dz * s));
        dda1313dndn = (((((dda1313dxdx * cc) as f64)
            - (((2f64 * (dda1313dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1313dzdz * ss) as f64)) as f32);
        (
            a1113,
            da1113dx,
            da1113dz,
            dda1113dxdx,
            dda1113dxdz,
            dda1113dzdz,
        ) = a11132.interp(xt, zt);
        da1113dn = ((da1113dx * c) - (da1113dz * s));
        dda1113dndn = (((((dda1113dxdx * cc) as f64)
            - (((2f64 * (dda1113dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1113dzdz * ss) as f64)) as f32);
        (
            a3313,
            da3313dx,
            da3313dz,
            dda3313dxdx,
            dda3313dxdz,
            dda3313dzdz,
        ) = a33132.interp(xt, zt);
        da3313dn = ((da3313dx * c) - (da3313dz * s));
        dda3313dndn = (((((dda3313dxdx * cc) as f64)
            - (((2f64 * (dda3313dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda3313dzdz * ss) as f64)) as f32);
        /* step 4 of 4th-order Runge-Kutta */
        px2 = (pxt * pxt);
        pz2 = (pzt * pzt);
        pxz = (pxt * pzt);
        /*anisotropy parameters*/
        gamma11 = (((a1111 * px2) + (a1313 * pz2)) + (((2i32 as f32) * a1113) * pxz));
        gamma33 = (((a3333 * pz2) + (a1313 * px2)) + (((2i32 as f32) * a3313) * pxz));
        gamma13 = ((((a1133 + a1313) * pxz) + (a1113 * px2)) + (a3313 * pz2));
        den = ((1i32 as f32) / ((gamma11 + gamma33) - (2i32 as f32)));
        g11 = ((gamma33 - (1i32 as f32)) * den);
        g33 = ((gamma11 - (1i32 as f32)) * den);
        g13 = ((-gamma13) * den);
        sxfactor = (((((((da1111dx * px2) * g11) + ((da3333dx * pz2) * g33))
            + ((((2i32 as f32) * (da1133dx + da1313dx)) * pxz) * g13))
            + (da1313dx * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * da3313dx) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * da1113dx) * ((pxz * g11) + (px2 * g13))));
        szfactor = (((((((da1111dz * px2) * g11) + ((da3333dz * pz2) * g33))
            + ((((2i32 as f32) * (da1133dz + da1313dz)) * pxz) * g13))
            + (da1313dz * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * da3313dz) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * da1113dz) * ((pxz * g11) + (px2 * g13))));
        snfact = ((sxfactor * c) - (szfactor * s));
        /*computing ray velocities and derivatives*/
        dxt = ((((((a1111 * pxt) * g11) + (((a1133 + a1313) * pzt) * g13))
            + ((a3313 * pzt) * g33))
            + (a1113 * ((pzt * g11) + (((2i32 as f32) * pxt) * g13))))
            + ((a1313 * g33) * pxt));
        dzt = ((((((a3333 * pzt) * g33) + (((a1133 + a1313) * pxt) * g13))
            + ((a1113 * pxt) * g11))
            + (a3313 * ((pxt * g33) + (((2i32 as f32) * pzt) * g13))))
            + ((a1313 * g11) * pzt));
        dpxt = (((-0.5f64) * (sxfactor as f64)) as f32);
        dpzt = (((-0.5f64) * (szfactor as f64)) as f32);
        dgamma11dpx = ((((2i32 as f32) * a1111) * pxt) + (((2i32 as f32) * a1113) * pzt));
        dgamma11dpz = ((((2i32 as f32) * a1313) * pzt) + (((2i32 as f32) * a1113) * pxt));
        dgamma33dpx = ((((2i32 as f32) * a1313) * pxt) + (((2i32 as f32) * a3313) * pzt));
        dgamma33dpz = ((((2i32 as f32) * a3333) * pzt) + (((2i32 as f32) * a3313) * pxt));
        dgamma13dpx = (((a1133 + a1313) * pzt) + (((2i32 as f32) * a1113) * pxt));
        dgamma13dpz = (((a1133 + a1313) * pxt) + (((2i32 as f32) * a3313) * pzt));
        dgamma11dn = (((da1111dn * px2) + (da1313dn * pz2)) + (((2i32 as f32) * da1113dn) * pxz));
        dgamma33dn = (((da3333dn * pz2) + (da1313dn * px2)) + (((2i32 as f32) * da3313dn) * pxz));
        dgamma13dn = ((((da1133dn + da1313dn) * pxz) + (da1113dn * px2)) + (da3313dn * pz2));
        dg11dpx = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxt)))
            * den)
            * den)
            + ((dgamma33dpx - ((2i32 as f32) * dxt)) * den));
        dg11dpz = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzt)))
            * den)
            * den)
            + ((dgamma33dpz - ((2i32 as f32) * dzt)) * den));
        dg33dpx = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxt)))
            * den)
            * den)
            + ((dgamma11dpx - ((2i32 as f32) * dxt)) * den));
        dg33dpz = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzt)))
            * den)
            * den)
            + ((dgamma11dpz - ((2i32 as f32) * dzt)) * den));
        dg13dpx = ((((gamma13 * ((dgamma11dpx + dgamma33dpx) - ((4i32 as f32) * dxt))) * den)
            * den)
            - (dgamma13dpx * den));
        dg13dpz = ((((gamma13 * ((dgamma11dpz + dgamma33dpz) - ((4i32 as f32) * dzt))) * den)
            * den)
            - (dgamma13dpz * den));
        dg11dn = (((((-(gamma33 - (1i32 as f32)))
            * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact)))
            * den)
            * den)
            + ((dgamma33dn - snfact) * den));
        dg33dn = (((((-(gamma11 - (1i32 as f32)))
            * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact)))
            * den)
            * den)
            + ((dgamma11dn - snfact) * den));
        dg13dn = ((((gamma13 * ((dgamma11dn + dgamma33dn) - ((2i32 as f32) * snfact))) * den)
            * den)
            - (dgamma13dn * den));
        ddxdpx = ((((((a1111 * pxt) * dg11dpx) + (((a1133 + a1313) * pzt) * dg13dpx))
            + ((a3313 * pzt) * dg33dpx))
            + (a1113 * ((pzt * dg11dpx) + (((2i32 as f32) * pxt) * dg13dpx))))
            + ((a1313 * dg33dpx) * pxt));
        ddzdpz = ((((((a3333 * pzt) * dg33dpz) + (((a1133 + a1313) * pxt) * dg13dpz))
            + ((a1113 * pxt) * dg11dpz))
            + (a3313 * ((pxt * dg33dpz) + (((2i32 as f32) * pzt) * dg13dpz))))
            + ((a1313 * dg11dpz) * pzt));
        ddxdpz = ((((((a1111 * pxt) * dg11dpz) + (((a1133 + a1313) * pzt) * dg13dpz))
            + ((a3313 * pzt) * dg33dpz))
            + (a1113 * ((pzt * dg11dpz) + (((2i32 as f32) * pxt) * dg13dpz))))
            + ((a1313 * dg33dpz) * pxt));
        dsnfactordn = (((((((da1111dn * px2) * dg11dn) + ((da3333dn * pz2) * dg33dn))
            + ((((2i32 as f32) * (da1133dn + da1313dn)) * pxz) * dg13dn))
            + (da1313dn * ((px2 * dg33dn) + (pz2 * dg11dn))))
            + (((2i32 as f32) * da3313dn) * ((pz2 * dg13dn) + (pxz * dg33dn))))
            + (((2i32 as f32) * da1113dn) * ((pxz * dg11dn) + (px2 * dg13dn))));
        ddxdn = ((((((a1111 * pxt) * dg11dn) + (((a1133 + a1313) * pzt) * dg13dn))
            + ((a3313 * pzt) * dg33dn))
            + (a1113 * ((pzt * dg11dn) + (((2i32 as f32) * pxt) * dg13dn))))
            + ((a1313 * dg33dn) * pxt));
        ddzdn = ((((((a3333 * pzt) * dg33dn) + (((a1133 + a1313) * pxt) * dg13dn))
            + ((a1113 * pxt) * dg11dn))
            + (a3313 * ((pxt * dg33dn) + (((2i32 as f32) * pzt) * dg13dn))))
            + ((a1313 * dg11dn) * pzt));
        dcdp1 = (((((a1111 * g11) + (a1313 * g33)) + (((2i32 as f32) * a1113) * g13)) + ddxdpx)
            - (dxt * dxt));
        dcdp3 = (((((a3333 * g33) + (a1313 * g11)) + (((2i32 as f32) * a3313) * g13)) + ddzdpz)
            - (dzt * dzt));
        dcdp13 = ((((((a1133 * g13) + (a1313 * g13)) + (a1113 * g11)) + (a3313 * g33)) + ddxdpz)
            - (dxt * dzt));
        ddcdqq = (((((dcdp1 * cc) as f64) - (((2f64 * (dcdp13 as f64)) * (s as f64)) * (c as f64)))
            + ((dcdp3 * ss) as f64)) as f32);
        dxdn = ((((((da1111dn * pxt) * g11) + (((da1133dn + da1313dn) * pzt) * g13))
            + ((da3313dn * pzt) * g33))
            + (da1113dn * ((pzt * g11) + (((2i32 as f32) * pxt) * g13))))
            + ((da1313dn * g33) * pxt));
        dzdn = ((((((da3333dn * pzt) * g33) + (((da1133dn + da1313dn) * pxt) * g13))
            + ((da1113dn * pxt) * g11))
            + (da3313dn * ((pxt * g33) + (((2i32 as f32) * pzt) * g13))))
            + ((da1313dn * g11) * pzt));
        ddcdpn = (((((((((dxdn * c) - (dzdn * s)) as f64)
            - (((0.5f64 * (dxt as f64)) * (sxfactor as f64)) * (cc as f64)))
            + (((0.5f64 * (((dxt * szfactor) + (dzt * sxfactor)) as f64)) * (s as f64))
                * (c as f64)))
            - (((0.5f64 * (dzt as f64)) * (szfactor as f64)) * (ss as f64)))
            + ((ddxdn * c) as f64))
            - ((ddzdn * s) as f64)) as f32);
        snfactor = (((((((dda1111dndn * px2) * g11) + ((dda3333dndn * pz2) * g33))
            + ((((2i32 as f32) * (dda1133dndn + dda1313dndn)) * pxz) * g13))
            + (dda1313dndn * ((px2 * g33) + (pz2 * g11))))
            + (((2i32 as f32) * dda3313dndn) * ((pz2 * g13) + (pxz * g33))))
            + (((2i32 as f32) * dda1113dndn) * ((pxz * g11) + (px2 * g13))));
        ddcdnn = ((((((0.5f64 * (snfactor as f64))
            - (((0.25f64 * (sxfactor as f64)) * (sxfactor as f64)) * (cc as f64)))
            + ((((0.5f64 * (sxfactor as f64)) * (szfactor as f64)) * (s as f64)) * (c as f64)))
            - (((0.25f64 * (szfactor as f64)) * (szfactor as f64)) * (ss as f64)))
            + (0.5f64 * (dsnfactordn as f64))) as f32);
        dp1t = (((-ddcdnn) * q1t) - (ddcdpn * p1t));
        dq1t = ((ddcdqq * p1t) + (ddcdpn * q1t));
        dp2t = (((-ddcdnn) * q2t) - (ddcdpn * p2t));
        dq2t = ((ddcdqq * p2t) + (ddcdpn * q2t));
        dxx = (((hsixth as f64) * (((dx + dxt) as f64) + (2f64 * (dxm as f64)))) as f32);
        dzz = (((hsixth as f64) * (((dz + dzt) as f64) + (2f64 * (dzm as f64)))) as f32);
        x += dxx;
        z += dzz;
        px = ((px as f64) + ((hsixth as f64) * (((dpx + dpxt) as f64) + (2f64 * (dpxm as f64)))))
            as f32;
        pz = ((pz as f64) + ((hsixth as f64) * (((dpz + dpzt) as f64) + (2f64 * (dpzm as f64)))))
            as f32;
        p1 = ((p1 as f64) + ((hsixth as f64) * (((dp1 + dp1t) as f64) + (2f64 * (dp1m as f64)))))
            as f32;
        q1 = ((q1 as f64) + ((hsixth as f64) * (((dq1 + dq1t) as f64) + (2f64 * (dq1m as f64)))))
            as f32;
        p2 = ((p2 as f64) + ((hsixth as f64) * (((dp2 + dp2t) as f64) + (2f64 * (dp2m as f64)))))
            as f32;
        q2 = ((q2 as f64) + ((hsixth as f64) * (((dq2 + dq2t) as f64) + (2f64 * (dq2m as f64)))))
            as f32;
        vp = (((1i32 as f64) / (((px * px) + (pz * pz)) as f64).sqrt()) as f32);
        s = (px * vp);
        c = (pz * vp);
        ss = (s * s);
        cc = (c * c);
        (
            a1111,
            da1111dx,
            da1111dz,
            dda1111dxdx,
            dda1111dxdz,
            dda1111dzdz,
        ) = a11112.interp(x, z);
        da1111dn = ((da1111dx * c) - (da1111dz * s));
        dda1111dndn = (((((dda1111dxdx * cc) as f64)
            - (((2f64 * (dda1111dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1111dzdz * ss) as f64)) as f32);
        (
            a3333,
            da3333dx,
            da3333dz,
            dda3333dxdx,
            dda3333dxdz,
            dda3333dzdz,
        ) = a33332.interp(x, z);
        da3333dn = ((da3333dx * c) - (da3333dz * s));
        dda3333dndn = (((((dda3333dxdx * cc) as f64)
            - (((2f64 * (dda3333dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda3333dzdz * ss) as f64)) as f32);
        (
            a1133,
            da1133dx,
            da1133dz,
            dda1133dxdx,
            dda1133dxdz,
            dda1133dzdz,
        ) = a11332.interp(x, z);
        da1133dn = ((da1133dx * c) - (da1133dz * s));
        dda1133dndn = (((((dda1133dxdx * cc) as f64)
            - (((2f64 * (dda1133dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1133dzdz * ss) as f64)) as f32);
        (
            a1313,
            da1313dx,
            da1313dz,
            dda1313dxdx,
            dda1313dxdz,
            dda1313dzdz,
        ) = a13132.interp(x, z);
        da1313dn = ((da1313dx * c) - (da1313dz * s));
        dda1313dndn = (((((dda1313dxdx * cc) as f64)
            - (((2f64 * (dda1313dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1313dzdz * ss) as f64)) as f32);
        (
            a1113,
            da1113dx,
            da1113dz,
            dda1113dxdx,
            dda1113dxdz,
            dda1113dzdz,
        ) = a11132.interp(x, z);
        da1113dn = ((da1113dx * c) - (da1113dz * s));
        dda1113dndn = (((((dda1113dxdx * cc) as f64)
            - (((2f64 * (dda1113dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda1113dzdz * ss) as f64)) as f32);
        (
            a3313,
            da3313dx,
            da3313dz,
            dda3313dxdx,
            dda3313dxdz,
            dda3313dzdz,
        ) = a33132.interp(x, z);
        da3313dn = ((da3313dx * c) - (da3313dz * s));
        dda3313dndn = (((((dda3313dxdx * cc) as f64)
            - (((2f64 * (dda3313dxdz as f64)) * (s as f64)) * (c as f64)))
            + ((dda3313dzdz * ss) as f64)) as f32);
        /* update time */
        t += dt;
        /* update kmah index */
        if ((((q2 as f64) <= 0f64) && ((q2old as f64) > 0f64))
            || (((q2 as f64) >= 0f64) && ((q2old as f64) < 0f64)))
        {
            kmah += 1;
        }
        /* save ray parameters */
        rs[(it) as usize].t = t;
        rs[(it) as usize].x = x;
        rs[(it) as usize].z = z;
        rs[(it) as usize].c = c;
        rs[(it) as usize].s = s;
        rs[(it) as usize].q1 = q1;
        rs[(it) as usize].p1 = p1;
        rs[(it) as usize].q2 = q2;
        rs[(it) as usize].p2 = p2;
        rs[(it) as usize].kmah = kmah;
        rs[(it) as usize].v = vp;
        rs[(it) as usize].dvdx =
            ((((0.5f64 * (da3333dx as f64)) * (vp as f64)) / (a3333 as f64)) as f32);
        rs[(it) as usize].dvdz =
            ((((0.5f64 * (da3333dz as f64)) * (vp as f64)) / (a3333 as f64)) as f32);
        it += 1;
    }
    /* free velocity interpolator */
    return Ray { nrs: it, rs };
}
