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

//! Numerical helpers from su/main/synthetics_waveforms_testpatterns/sufctanismod.c
//! (Revision 1.8). Author: Tong Fei, CWP, Colorado School of Mines, Dec 1993.

#![allow(unused_mut, unused_variables, unused_assignments, unused_parens)]

use super::Grid;
use oxwave_core::prelude::*;
use std::f64::consts::PI;

pub(super) fn anis_solver2(
    mut it: i32,
    mut u: &mut Grid,
    mut v: &mut Grid,
    mut w: &mut Grid,
    mut e11: &mut Grid,
    mut e33: &mut Grid,
    mut e23: &mut Grid,
    mut e12: &mut Grid,
    mut e13: &mut Grid,
    mut aa: &Grid,
    mut cc: &Grid,
    mut ff: &Grid,
    mut ll: &Grid,
    mut nn: &Grid,
    mut rho: &Grid,
    mut xzsource: &Grid,
    mut fx: &Vec<f32>,
    mut fy: &Vec<f32>,
    mut fz: &Vec<f32>,
    mut dx: f32,
    mut dz: f32,
    mut dt: f32,
    mut nx: i32,
    mut nz: i32,
    mut dofct: i32,
    mut isurf: i32,
    mut eta0: f32,
    mut eta: f32,
    mut deta0dx: f32,
    mut deta0dz: f32,
    mut detadx: f32,
    mut detadz: f32,
    mut mbx1: i32,
    mut mbx2: i32,
    mut mbz1: i32,
    mut mbz2: i32,
) {
    let mut ix: i32 = 0i32;
    let mut iz: i32 = 0i32;
    let mut dtdx: f32 = (dt / dx);
    let mut dtdz: f32 = (dt / dz);
    let mut rdxdz: f32 = (((1.0f64 / (dx as f64)) / (dz as f64)) as f32);
    let mut u1: Grid = Grid::default();
    let mut f0: Grid = Grid::default();
    let mut f1: Grid = Grid::default();
    u1 = Grid::new(nz, nx);
    f0 = Grid::new(nz, nx);
    f1 = Grid::new(nz, nx);
    // Evolve u, then w, then v, then the five strain components. The same
    // temporary u1 is reused between components, just as in the C solver.
    difference(
        mbx1, mbx2, mbz1, mbz2, it, &mut u1, &u, &xzsource, &fx, dt, rdxdz,
    );
    difference_2x(mbx1, mbx2, mbz1, mbz2, (-1i32), dtdx, &mut u1, &e11, &aa);
    difference_2x(mbx1, mbx2, mbz1, mbz2, (-1i32), dtdx, &mut u1, &e33, &ff);
    difference_2z(
        mbx1,
        mbx2,
        mbz1,
        mbz2,
        (-1i32),
        ((2.0f64 * (dtdz as f64)) as f32),
        &mut u1,
        &e13,
        &ll,
    );
    if (dofct != 0) {
        fct2d1o(
            &u, &mut u1, nx, nz, eta0, eta, deta0dx, deta0dz, detadx, detadz, dx, dz, mbx1, mbx2,
            mbz1, mbz2, &mut f0, &mut f1,
        );
    }
    absorb1(
        it, nx, nz, dx, dz, dt, &u, &mut u1, &cc, &rho, isurf, mbx1, mbx2, mbz1, mbz2, 1i32, 1i32,
    );
    absorb1(
        it, nx, nz, dx, dz, dt, &u, &mut u1, &ll, &rho, isurf, mbx1, mbx2, mbz1, mbz2, 1i32, 1i32,
    );
    ix = mbx1;
    while (ix < mbx2) {
        {
            iz = mbz1;
            while (iz < mbz2) {
                {
                    u[(ix) as usize][(iz) as usize] = u1[(ix) as usize][(iz) as usize];
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    difference(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        it,
        &mut u1,
        &w,
        &xzsource,
        &fz,
        dt,
        rdxdz,
    );
    difference_2x(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        0i32,
        ((2.0f64 * (dtdx as f64)) as f32),
        &mut u1,
        &e13,
        &ll,
    );
    difference_2z(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        0i32,
        dtdx,
        &mut u1,
        &e11,
        &ff,
    );
    difference_2z(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        0i32,
        dtdz,
        &mut u1,
        &e33,
        &cc,
    );
    if (dofct != 0) {
        fct2d1o(
            &w,
            &mut u1,
            nx,
            nz,
            eta0,
            eta,
            deta0dx,
            deta0dz,
            detadx,
            detadz,
            dx,
            dz,
            mbx1,
            (mbx2 - 1i32),
            mbz1,
            (mbz2 - 1i32),
            &mut f0,
            &mut f1,
        );
    }
    absorb1(
        it,
        nx,
        nz,
        dx,
        dz,
        dt,
        &w,
        &mut u1,
        &cc,
        &rho,
        isurf,
        mbx1,
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        0i32,
        1i32,
    );
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = mbz1;
            while (iz < (mbz2 - 1i32)) {
                {
                    w[(ix) as usize][(iz) as usize] = u1[(ix) as usize][(iz) as usize];
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    difference(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        it,
        &mut u1,
        &v,
        &xzsource,
        &fy,
        dt,
        rdxdz,
    );
    difference_2x(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        0i32,
        ((2.0f64 * (dtdx as f64)) as f32),
        &mut u1,
        &e12,
        &nn,
    );
    difference_2z(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        0i32,
        ((2.0f64 * (dtdz as f64)) as f32),
        &mut u1,
        &e23,
        &ll,
    );
    if (dofct != 0) {
        fct2d1o(
            &v,
            &mut u1,
            nx,
            nz,
            eta0,
            eta,
            deta0dx,
            deta0dz,
            detadx,
            detadz,
            dx,
            dz,
            mbx1,
            (mbx2 - 1i32),
            mbz1,
            (mbz2 - 1i32),
            &mut f0,
            &mut f1,
        );
    }
    absorb1(
        it,
        nx,
        nz,
        dx,
        dz,
        dt,
        &v,
        &mut u1,
        &ll,
        &rho,
        isurf,
        mbx1,
        (mbx2 - 1i32),
        mbz1,
        (mbz2 - 1i32),
        0i32,
        1i32,
    );
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = mbz1;
            while (iz < (mbz2 - 1i32)) {
                {
                    v[(ix) as usize][(iz) as usize] = u1[(ix) as usize][(iz) as usize];
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    difference(
        (mbx1 - 1i32),
        mbx2,
        (mbz1 - 1i32),
        (mbz2 + 1i32),
        it,
        &mut u1,
        &e11,
        &xzsource,
        &fz,
        (0.0f64 as f32),
        (0.0f64 as f32),
    );
    difference_2x(
        (mbx1 - 1i32),
        mbx2,
        (mbz1 - 1i32),
        (mbz2 + 1i32),
        0i32,
        dtdx,
        &mut u1,
        &u,
        &rho,
    );
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = mbz1;
            while (iz < mbz2) {
                {
                    e11[(ix) as usize][(iz) as usize] = u1[(ix) as usize][(iz) as usize];
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    difference(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        mbz2,
        it,
        &mut u1,
        &e33,
        &xzsource,
        &fz,
        (0.0f64 as f32),
        (0.0f64 as f32),
    );
    difference_2z(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        mbz2,
        (-1i32),
        dtdx,
        &mut u1,
        &w,
        &rho,
    );
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = (mbz1 + 1i32);
            while (iz < (mbz2 - 1i32)) {
                {
                    e33[(ix) as usize][(iz) as usize] = u1[(ix) as usize][(iz) as usize];
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    difference(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        mbz2,
        it,
        &mut u1,
        &e23,
        &xzsource,
        &fz,
        (0.0f64 as f32),
        (0.0f64 as f32),
    );
    difference_2z(
        (mbx1 - 1i32),
        mbx2,
        mbz1,
        mbz2,
        (-1i32),
        ((0.5f64 * (dtdx as f64)) as f32),
        &mut u1,
        &v,
        &rho,
    );
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = (mbz1 + 1i32);
            while (iz < (mbz2 - 1i32)) {
                {
                    e23[(ix) as usize][(iz) as usize] = u1[(ix) as usize][(iz) as usize];
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    difference(
        mbx1,
        mbx2,
        (mbz1 - 1i32),
        mbz2,
        it,
        &mut u1,
        &e12,
        &xzsource,
        &fz,
        (0.0f64 as f32),
        (0.0f64 as f32),
    );
    difference_2x(
        mbx1,
        mbx2,
        (mbz1 - 1i32),
        mbz2,
        (-1i32),
        ((0.5f64 * (dtdx as f64)) as f32),
        &mut u1,
        &v,
        &rho,
    );
    absorb1(
        it,
        nx,
        nz,
        dx,
        dz,
        dt,
        &e12,
        &mut u1,
        &ll,
        &rho,
        isurf,
        mbx1,
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        1i32,
        0i32,
    );
    ix = mbx1;
    while (ix < mbx2) {
        {
            iz = mbz1;
            while (iz < (mbz2 - 1i32)) {
                {
                    e12[(ix) as usize][(iz) as usize] = u1[(ix) as usize][(iz) as usize];
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    difference(
        mbx1,
        mbx2,
        (mbz1 - 1i32),
        mbz2,
        it,
        &mut u1,
        &e13,
        &xzsource,
        &fz,
        (0.0f64 as f32),
        (0.0f64 as f32),
    );
    difference_2x(
        mbx1,
        mbx2,
        (mbz1 - 1i32),
        mbz2,
        (-1i32),
        ((0.5f64 * (dtdx as f64)) as f32),
        &mut u1,
        &w,
        &rho,
    );
    difference_2z(
        mbx1,
        mbx2,
        (mbz1 - 1i32),
        mbz2,
        0i32,
        ((0.5f64 * (dtdx as f64)) as f32),
        &mut u1,
        &u,
        &rho,
    );
    absorb1(
        it,
        nx,
        nz,
        dx,
        dz,
        dt,
        &e13,
        &mut u1,
        &cc,
        &rho,
        isurf,
        mbx1,
        mbx2,
        mbz1,
        (mbz2 - 1i32),
        1i32,
        0i32,
    );
    ix = mbx1;
    while (ix < mbx2) {
        {
            iz = mbz1;
            while (iz < (mbz2 - 1i32)) {
                {
                    e13[(ix) as usize][(iz) as usize] = u1[(ix) as usize][(iz) as usize];
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    ();
    ();
    ();
}

pub(super) fn absorb1(
    mut it: i32,
    mut nx: i32,
    mut nz: i32,
    mut dx: f32,
    mut dz: f32,
    mut dt: f32,
    mut u: &Grid,
    mut u1: &mut Grid,
    mut cc: &Grid,
    mut rho: &Grid,
    mut isurf: i32,
    mut mbx1: i32,
    mut mbx2: i32,
    mut mbz1: i32,
    mut mbz2: i32,
    mut side: i32,
    mut tb: i32,
) {
    let mut i: i32 = 0i32;
    let mut j: i32 = 0i32;
    let mut k: f32 = 0.0f32;
    let mut dtdz: f32 = (dt / dz);
    let mut dtdx: f32 = (dt / dx);
    let mut vell: Vec<f32> = Vec::new();
    let mut velr: Vec<f32> = Vec::new();
    let mut velt: Vec<f32> = Vec::new();
    let mut velb: Vec<f32> = Vec::new();
    // Match the minimal C fix: absorb1's side stencil reads both adjacent
    // rows (sufctanismod.c:1101,1103), which do not exist for nx=1.
    if side != 0 && mbx2 - mbx1 < 2 {
        err!("absorb1: side boundaries require at least two x grid points");
    }
    vell = vec![0.0f32; (nz) as usize];
    velr = vec![0.0f32; (nz) as usize];
    velt = vec![0.0f32; (nx) as usize];
    velb = vec![0.0f32; (nx) as usize];
    boundary_vel(
        &cc, &rho, &mut vell, &mut velr, &mut velt, &mut velb, nx, nz, mbx1, mbz1, mbx2, mbz2,
    );
    if (tb != 0) {
        i = mbx1;
        while (i < mbx2) {
            {
                k = ((-dtdz) * velt[(i) as usize]);
                if ((k > (1i32 as f32)) || (k < ((-1i32) as f32))) {
                    warn!("WARNING:Big k-values.");
                }
                if (isurf == 1i32) {
                    u1[(i) as usize][(mbz1) as usize] = (u[(i) as usize][(mbz1) as usize]
                        - (k * (u[(i) as usize][(mbz1 + 1i32) as usize]
                            - u[(i) as usize][(mbz1) as usize])));
                } else {
                    if (isurf == 2i32) {
                        u1[(i) as usize][(mbz1) as usize] =
                            u1[(i) as usize][(mbz1 + 1i32) as usize];
                    } else {
                        if (isurf == 3i32) {
                            u1[(i) as usize][(mbz1) as usize] = (0.0f64 as f32);
                        }
                    }
                }
                k = (dtdz * velb[(i) as usize]);
                u1[(i) as usize][(mbz2 - 1i32) as usize] = (u[(i) as usize]
                    [(mbz2 - 1i32) as usize]
                    - (k * (u[(i) as usize][(mbz2 - 1i32) as usize]
                        - u[(i) as usize][(mbz2 - 2i32) as usize])));
            }
            i += 1;
        }
    }
    if (side != 0) {
        j = mbz1;
        while (j < mbz2) {
            {
                k = (dtdx * velr[(j) as usize]);
                u1[(mbx2 - 1i32) as usize][(j) as usize] = (u[(mbx2 - 1i32) as usize]
                    [(j) as usize]
                    - (k * (u[(mbx2 - 1i32) as usize][(j) as usize]
                        - u[(mbx2 - 2i32) as usize][(j) as usize])));
                k = ((-dtdx) * vell[(j) as usize]);
                u1[(mbx1) as usize][(j) as usize] = (u[(mbx1) as usize][(j) as usize]
                    - (k * (u[(mbx1 + 1i32) as usize][(j) as usize]
                        - u[(mbx1) as usize][(j) as usize])));
            }
            j += 1;
        }
    }
    ();
    ();
    ();
    ();
}

pub(super) fn boundary_vel(
    mut cc: &Grid,
    mut rho: &Grid,
    mut vell: &mut Vec<f32>,
    mut velr: &mut Vec<f32>,
    mut velt: &mut Vec<f32>,
    mut velb: &mut Vec<f32>,
    mut nx: i32,
    mut nz: i32,
    mut mbx1: i32,
    mut mbz1: i32,
    mut mbx2: i32,
    mut mbz2: i32,
) {
    let mut i: i32 = 0i32;
    let mut k: i32 = 0i32;
    k = 0i32;
    while (k < nz) {
        {
            vell[(k) as usize] = (((cc[(mbx1) as usize][(k) as usize]
                * rho[(mbx1) as usize][(k) as usize]) as f64)
                .sqrt() as f32);
            velr[(k) as usize] = (((cc[(mbx2 - 1i32) as usize][(k) as usize]
                * rho[(mbx2 - 1i32) as usize][(k) as usize])
                as f64)
                .sqrt() as f32);
        }
        k += 1;
    }
    i = 0i32;
    while (i < nx) {
        {
            velt[(i) as usize] = (((cc[(i) as usize][(mbz1) as usize]
                * rho[(i) as usize][(mbz1) as usize]) as f64)
                .sqrt() as f32);
            velb[(i) as usize] = (((cc[(i) as usize][(mbz2 - 1i32) as usize]
                * rho[(i) as usize][(mbz2 - 1i32) as usize])
                as f64)
                .sqrt() as f32);
        }
        i += 1;
    }
}

pub(super) fn fct2d1o(
    mut u: &Grid,
    mut u1: &mut Grid,
    mut nx: i32,
    mut nz: i32,
    mut eta0: f32,
    mut eta: f32,
    mut deta0dx: f32,
    mut deta0dz: f32,
    mut detadx: f32,
    mut detadz: f32,
    mut dx: f32,
    mut dz: f32,
    mut mbx1: i32,
    mut mbx2: i32,
    mut mbz1: i32,
    mut mbz2: i32,
    mut f0: &mut Grid,
    mut f1: &mut Grid,
) {
    let mut ix: i32 = 0i32;
    let mut iz: i32 = 0i32;
    let mut sf: f32 = 0.0f32;
    let mut eta00: f32 = 0.0f32;
    let mut eta11: f32 = 0.0f32;
    let mut min: f32 = 0.0f32;
    let mut deta0dxdx: f32 = (deta0dx * dx);
    let mut deta0dzdz: f32 = (deta0dz * dz);
    let mut detadxdx: f32 = (detadx * dx);
    let mut detadzdz: f32 = (detadz * dz);
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = mbz1;
            while (iz < (mbz2 - 1i32)) {
                {
                    eta00 = (((eta0 as f64)
                        * ((1.0f64 + ((deta0dxdx * ((ix - 1i32) as f32)) as f64))
                            + ((deta0dzdz * ((iz - 1i32) as f32)) as f64)))
                        as f32);
                    eta11 = (((eta as f64)
                        * ((1.0f64 + ((detadxdx * ((ix - 1i32) as f32)) as f64))
                            + ((detadzdz * ((iz - 1i32) as f32)) as f64)))
                        as f32);
                    f0[(ix) as usize][(iz) as usize] = (eta00
                        * (u[(ix + 1i32) as usize][(iz) as usize]
                            - u[(ix) as usize][(iz) as usize]));
                    f1[(ix) as usize][(iz) as usize] = (eta11
                        * (u1[(ix + 1i32) as usize][(iz) as usize]
                            - u1[(ix) as usize][(iz) as usize]));
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            eta00 = (((eta0 as f64)
                * ((1.0f64 + ((deta0dxdx * ((ix - 1i32) as f32)) as f64))
                    + ((deta0dzdz * ((mbz2 - 1i32) as f32)) as f64))) as f32);
            eta11 = (((eta as f64)
                * ((1.0f64 + ((detadxdx * ((ix - 1i32) as f32)) as f64))
                    + ((detadzdz * ((mbz2 - 1i32) as f32)) as f64))) as f32);
            f0[(ix) as usize][(mbz2 - 1i32) as usize] = (eta00
                * (u[(ix + 1i32) as usize][(mbz2 - 1i32) as usize]
                    - u[(ix) as usize][(mbz2 - 1i32) as usize]));
            f1[(ix) as usize][(mbz2 - 1i32) as usize] = (eta11
                * (u1[(ix + 1i32) as usize][(mbz2 - 1i32) as usize]
                    - u1[(ix) as usize][(mbz2 - 1i32) as usize]));
        }
        ix += 1;
    }
    ix = (mbx1 + 1i32);
    while (ix < (mbx2 - 1i32)) {
        {
            iz = (mbz1 + 1i32);
            while (iz < (mbz2 - 1i32)) {
                {
                    u1[(ix) as usize][(iz) as usize] = (u1[(ix) as usize][(iz) as usize]
                        + (f0[(ix) as usize][(iz) as usize]
                            - f0[(ix - 1i32) as usize][(iz) as usize]));
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = mbz1;
            while (iz < (mbz2 - 1i32)) {
                {
                    f0[(ix) as usize][(iz) as usize] = (u1[(ix + 1i32) as usize][(iz) as usize]
                        - u1[(ix) as usize][(iz) as usize]);
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            f0[(ix) as usize][(mbz2 - 1i32) as usize] = (u1[(ix + 1i32) as usize]
                [(mbz2 - 1i32) as usize]
                - u1[(ix) as usize][(mbz2 - 1i32) as usize]);
        }
        ix += 1;
    }
    ix = (mbx1 + 1i32);
    while (ix < (mbx2 - 2i32)) {
        {
            iz = mbz1;
            while (iz < mbz2) {
                {
                    if (f1[(ix) as usize][(iz) as usize] > (0i32 as f32)) {
                        sf = (1.0f64 as f32);
                    } else {
                        sf = ((-1.0f64) as f32);
                    }
                    if ((sf * f0[(ix - 1i32) as usize][(iz) as usize])
                        > (sf * f1[(ix) as usize][(iz) as usize]))
                    {
                        min = (sf * f1[(ix) as usize][(iz) as usize]);
                    } else {
                        min = (sf * f0[(ix - 1i32) as usize][(iz) as usize]);
                    }
                    if (min > (sf * f0[(ix + 1i32) as usize][(iz) as usize])) {
                        min = (sf * f0[(ix + 1i32) as usize][(iz) as usize]);
                    }
                    if ((min as f64) > 0.0f64) {
                        f1[(ix) as usize][(iz) as usize] = (sf * min);
                    } else {
                        f1[(ix) as usize][(iz) as usize] = (0.0f64 as f32);
                    }
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    iz = mbz1;
    while (iz < mbz2) {
        {
            if ((f1[(mbx1) as usize][(iz) as usize] as f64) > 0.0f64) {
                sf = (1.0f64 as f32);
            } else {
                sf = ((-1.0f64) as f32);
            }
            if ((sf * f1[(mbx1) as usize][(iz) as usize])
                > (sf * f0[(mbx1 + 1i32) as usize][(iz) as usize]))
            {
                min = (sf * f0[(mbx1 + 1i32) as usize][(iz) as usize]);
            } else {
                min = (sf * f1[(mbx1) as usize][(iz) as usize]);
            }
            if ((min as f64) > 0.0f64) {
                f1[(mbx1) as usize][(iz) as usize] = (sf * min);
            } else {
                f1[(mbx1) as usize][(iz) as usize] = (0.0f64 as f32);
            }
            if ((f1[(mbx2 - 2i32) as usize][(iz) as usize] as f64) > 0.0f64) {
                sf = (1.0f64 as f32);
            } else {
                sf = ((-1.0f64) as f32);
            }
            if ((sf * f0[(mbx2 - 3i32) as usize][(iz) as usize])
                > (sf * f1[(mbx2 - 2i32) as usize][(iz) as usize]))
            {
                min = (sf * f1[(mbx2 - 2i32) as usize][(iz) as usize]);
            } else {
                min = (sf * f0[(mbx2 - 3i32) as usize][(iz) as usize]);
            }
            if ((min as f64) > 0.0f64) {
                f1[(mbx2 - 2i32) as usize][(iz) as usize] = (sf * min);
            } else {
                f1[(mbx2 - 2i32) as usize][(iz) as usize] = (0.0f64 as f32);
            }
        }
        iz += 1;
    }
    ix = (mbx1 + 1i32);
    while (ix < (mbx2 - 1i32)) {
        {
            iz = (mbz1 + 1i32);
            while (iz < (mbz2 - 1i32)) {
                {
                    u1[(ix) as usize][(iz) as usize] = (u1[(ix) as usize][(iz) as usize]
                        - (f1[(ix) as usize][(iz) as usize]
                            - f1[(ix - 1i32) as usize][(iz) as usize]));
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = mbz1;
            while (iz < (mbz2 - 1i32)) {
                {
                    eta00 = (((eta0 as f64)
                        * ((1.0f64 + ((deta0dxdx * ((ix - 1i32) as f32)) as f64))
                            + ((deta0dzdz * ((iz - 1i32) as f32)) as f64)))
                        as f32);
                    eta11 = (((eta as f64)
                        * ((1.0f64 + ((detadxdx * ((ix - 1i32) as f32)) as f64))
                            + ((detadzdz * ((iz - 1i32) as f32)) as f64)))
                        as f32);
                    f0[(ix) as usize][(iz) as usize] = (eta00
                        * (u[(ix) as usize][(iz + 1i32) as usize]
                            - u[(ix) as usize][(iz) as usize]));
                    f1[(ix) as usize][(iz) as usize] = (eta11
                        * (u1[(ix) as usize][(iz + 1i32) as usize]
                            - u1[(ix) as usize][(iz) as usize]));
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    iz = mbz1;
    while (iz < (mbz2 - 1i32)) {
        {
            eta00 = (((eta0 as f64)
                * ((1.0f64 + ((deta0dxdx * ((mbx2 - 1i32) as f32)) as f64))
                    + ((deta0dzdz * ((mbz2 - 1i32) as f32)) as f64))) as f32);
            eta11 = (((eta as f64)
                * ((1.0f64 + ((detadxdx * ((mbx2 - 1i32) as f32)) as f64))
                    + ((detadzdz * ((mbz2 - 1i32) as f32)) as f64))) as f32);
            f0[(mbx2 - 1i32) as usize][(iz) as usize] = (eta00
                * (u[(mbx2 - 1i32) as usize][(iz + 1i32) as usize]
                    - u[(mbx2 - 1i32) as usize][(iz) as usize]));
            f1[(mbx2 - 1i32) as usize][(iz) as usize] = (eta11
                * (u1[(mbx2 - 1i32) as usize][(iz + 1i32) as usize]
                    - u1[(mbx2 - 1i32) as usize][(iz) as usize]));
        }
        iz += 1;
    }
    ix = (mbx1 + 1i32);
    while (ix < (mbx2 - 1i32)) {
        {
            iz = (mbz1 + 1i32);
            while (iz < (mbz2 - 1i32)) {
                {
                    u1[(ix) as usize][(iz) as usize] = (u1[(ix) as usize][(iz) as usize]
                        + (f0[(ix) as usize][(iz) as usize]
                            - f0[(ix) as usize][(iz - 1i32) as usize]));
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    ix = mbx1;
    while (ix < (mbx2 - 1i32)) {
        {
            iz = mbz1;
            while (iz < (mbz2 - 1i32)) {
                {
                    f0[(ix) as usize][(iz) as usize] = (u1[(ix) as usize][(iz + 1i32) as usize]
                        - u1[(ix) as usize][(iz) as usize]);
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    // Preserve the original x-boundary index used as this z-loop's start.
    iz = mbx1;
    while (iz < (mbz2 - 1i32)) {
        {
            f0[(mbx2 - 1i32) as usize][(iz) as usize] = (u1[(mbx2 - 1i32) as usize]
                [(iz + 1i32) as usize]
                - u1[(mbx2 - 1i32) as usize][(iz) as usize]);
        }
        iz += 1;
    }
    ix = mbx1;
    while (ix < mbx2) {
        {
            iz = (mbz1 + 1i32);
            while (iz < (mbz2 - 2i32)) {
                {
                    if (f1[(ix) as usize][(iz) as usize] > (0i32 as f32)) {
                        sf = (1.0f64 as f32);
                    } else {
                        sf = ((-1.0f64) as f32);
                    }
                    if ((sf * f0[(ix) as usize][(iz - 1i32) as usize])
                        > (sf * f1[(ix) as usize][(iz) as usize]))
                    {
                        min = (sf * f1[(ix) as usize][(iz) as usize]);
                    } else {
                        min = (sf * f0[(ix) as usize][(iz - 1i32) as usize]);
                    }
                    if (min > (sf * f0[(ix) as usize][(iz + 1i32) as usize])) {
                        min = (sf * f0[(ix) as usize][(iz + 1i32) as usize]);
                    }
                    if ((min as f64) > 0.0f64) {
                        f1[(ix) as usize][(iz) as usize] = (sf * min);
                    } else {
                        f1[(ix) as usize][(iz) as usize] = (0.0f64 as f32);
                    }
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    ix = mbx1;
    while (ix < mbx2) {
        {
            if ((f1[(ix) as usize][(mbz1) as usize] as f64) > 0.0f64) {
                sf = (1.0f64 as f32);
            } else {
                sf = ((-1.0f64) as f32);
            }
            if ((sf * f1[(ix) as usize][(mbz1) as usize])
                > (sf * f0[(ix) as usize][(mbz1 + 1i32) as usize]))
            {
                min = (sf * f0[(ix) as usize][(mbz1 + 1i32) as usize]);
            } else {
                min = (sf * f1[(ix) as usize][(mbz1) as usize]);
            }
            if ((min as f64) > 0.0f64) {
                f1[(ix) as usize][(mbz1) as usize] = (sf * min);
            } else {
                f1[(ix) as usize][(mbz1) as usize] = (0.0f64 as f32);
            }
            if ((f1[(ix) as usize][(mbz2 - 2i32) as usize] as f64) > 0.0f64) {
                sf = (1.0f64 as f32);
            } else {
                sf = ((-1.0f64) as f32);
            }
            if ((sf * f0[(ix) as usize][(mbz2 - 3i32) as usize])
                > (sf * f1[(ix) as usize][(mbz2 - 2i32) as usize]))
            {
                min = (sf * f1[(ix) as usize][(mbz2 - 2i32) as usize]);
            } else {
                min = (sf * f0[(ix) as usize][(mbz2 - 3i32) as usize]);
            }
            if ((min as f64) > 0.0f64) {
                f1[(ix) as usize][(mbz2 - 2i32) as usize] = (sf * min);
            } else {
                f1[(ix) as usize][(mbz2 - 2i32) as usize] = (0.0f64 as f32);
            }
        }
        ix += 1;
    }
    ix = (mbx1 + 1i32);
    while (ix < (mbx2 - 1i32)) {
        {
            iz = (mbz1 + 1i32);
            while (iz < (mbz2 - 1i32)) {
                {
                    u1[(ix) as usize][(iz) as usize] = (u1[(ix) as usize][(iz) as usize]
                        - (f1[(ix) as usize][(iz) as usize]
                            - f1[(ix) as usize][(iz - 1i32) as usize]));
                }
                iz += 1;
            }
        }
        ix += 1;
    }
}

pub(super) fn tforce_ricker(mut nt: i32, mut tforce: &mut Vec<f32>, mut dt: f32, mut fpeak: f32) {
    let mut it: i32 = 0i32;
    let mut t1: f32 = 0.0f32;
    let mut t0: f32 = 0.0f32;
    t0 = ((1.0f64 / (fpeak as f64)) as f32);
    it = 0i32;
    while (it < nt) {
        {
            t1 = ((it as f32) * dt);
            tforce[(it) as usize] = ((((((((-PI) * PI) * (fpeak as f64)) * (fpeak as f64))
                * ((t1 - t0) as f64))
                * ((t1 - t0) as f64))
                .exp()
                * (1.0f64
                    - ((((((2.0f64 * PI) * PI) * (fpeak as f64)) * (fpeak as f64))
                        * ((t1 - t0) as f64))
                        * ((t1 - t0) as f64)))) as f32);
        }
        it += 1;
    }
}

pub(super) fn tforce_akb(mut nt: i32, mut tforce: &mut Vec<f32>, mut dt: f32, mut fpeak: f32) {
    let mut it: i32 = 0i32;
    let mut t1: f32 = 0.0f32;
    let mut t0: f32 = ((1.8f64 / (fpeak as f64)) as f32);
    it = 0i32;
    while (it < nt) {
        {
            t1 = ((it as f32) * dt);
            tforce[(it) as usize] = ((((-60.0f64) * ((t1 - t0) as f64))
                * (((((-2.0f64) * (fpeak as f64)) * (fpeak as f64)) * ((t1 - t0) as f64))
                    * ((t1 - t0) as f64))
                    .exp()) as f32);
        }
        it += 1;
    }
}

pub(super) fn tforce_spike(mut nt: i32, mut tforce: &mut Vec<f32>, mut dt: f32, mut fpeak: f32) {
    let mut it: i32 = 0i32;
    it = 0i32;
    while (it < nt) {
        {
            tforce[(it) as usize] = (0.0f64 as f32);
        }
        it += 1;
    }
    tforce[(1i32) as usize] = (1.0f64 as f32);
}

pub(super) fn tforce_unit(mut nt: i32, mut tforce: &mut Vec<f32>, mut dt: f32, mut fpeak: f32) {
    let mut it: i32 = 0i32;
    let mut t1: f32 = 0.0f32;
    it = 0i32;
    while (it < nt) {
        {
            t1 = ((it as f32) * dt);
            tforce[(it) as usize] = (1.0f64 as f32);
        }
        it += 1;
    }
}

pub(super) fn locate_source(
    mut nx: i32,
    mut nz: i32,
    mut sx: i32,
    mut sz: i32,
    mut xzsource: &mut Grid,
    mut source: i32,
) {
    let mut ix: i32 = 0i32;
    let mut iz: i32 = 0i32;
    ix = 0i32;
    while (ix < nx) {
        {
            iz = 0i32;
            while (iz < nz) {
                {
                    xzsource[(ix) as usize][(iz) as usize] = (0.0f64 as f32);
                }
                iz += 1;
            }
        }
        ix += 1;
    }
    if (source == 1i32) {
        set_source(xzsource, sx, sz);
    }
    if (source == 2i32) {
        ix = 0i32;
        while (ix < (nx / 14i32)) {
            {
                set_source(xzsource, ix, nz / 3i32);
            }
            ix += 1;
        }
        ix = (nx / 14i32);
        while (ix < ((nx * 17i32) / 70i32)) {
            {
                set_source(
                    xzsource,
                    ix,
                    (((nz / 3i32) as f64)
                        + ((PI / 4.0f64).tan() * ((ix as f64) - ((nx as f64) / 14.0f64))))
                        as i32,
                );
            }
            ix += 1;
        }
        ix = ((17i32 * nx) / 70i32);
        while (ix < ((80i32 * nx) / 80i32)) {
            {
                set_source(xzsource, ix, (22i32 * nz) / 30i32);
            }
            ix += 1;
        }
    }
    if (source == 3i32) {
        ix = 0i32;
        while (ix < (nx / 14i32)) {
            {
                set_source(xzsource, ix, nz / 3i32);
            }
            ix += 1;
        }
        ix = (nx / 14i32);
        while (ix < ((nx * 34i32) / 70i32)) {
            {
                set_source(
                    xzsource,
                    ix,
                    (((nz / 3i32) as f64)
                        + ((PI / 18.0f64).tan() * ((ix as f64) - ((nx as f64) / 14.0f64))))
                        as i32,
                );
            }
            ix += 1;
        }
        ix = ((17i32 * nx) / 70i32);
        while (ix < ((80i32 * nx) / 80i32)) {
            {
                set_source(xzsource, ix, (22i32 * nz) / 30i32);
            }
            ix += 1;
        }
    }
}

fn set_source(xzsource: &mut Grid, sx: i32, sz: i32) {
    let Some(index) = xzsource.c_index(sx, sz) else {
        err!(
            "Source sx=%i sz=%i is outside the allocated grid. Exit",
            sx,
            sz
        );
    };
    xzsource.values[index] = 1.0;
}

pub(super) fn moving_bc(
    mut it: i32,
    mut nx: i32,
    mut nz: i32,
    mut sx: i32,
    mut sz: i32,
    mut dx: f32,
    mut dz: f32,
    mut impulse: i32,
    mut movebc: i32,
    mut t: &Vec<f32>,
    mut vmax: f32,
    mut mbx1: &mut i32,
    mut mbz1: &mut i32,
    mut mbx2: &mut i32,
    mut mbz2: &mut i32,
) {
    let mut temp: f32 = 0.0f32;
    let mut dzdx: f32 = (dz / dx);
    vmax = (((vmax as f64) * 1.5f64) as f32);
    temp = ((vmax * t[(it) as usize]) / dz);
    if (movebc != 0) {
        if (impulse != 0) {
            if ((0i32 as f64) > ((((sz as f32) - temp) as f64) - (0.04f64 * (nz as f64)))) {
                *mbz1 = 0i32;
            } else {
                *mbz1 = (((((sz as f32) - temp) as f64) - (0.04f64 * (nz as f64))) as i32);
            }
            if ((nz as f64) > ((((sz as f32) + temp) as f64) + (0.04f64 * (nz as f64)))) {
                *mbz2 = (((((sz as f32) + temp) as f64) + (0.04f64 * (nz as f64))) as i32);
            } else {
                *mbz2 = nz;
            }
            if ((0i32 as f64) > ((((sx as f32) - (temp * dzdx)) as f64) - (0.04f64 * (nx as f64))))
            {
                *mbx1 = 0i32;
            } else {
                *mbx1 = (((((sx as f32) - (temp * dzdx)) as f64) - (0.04f64 * (nx as f64))) as i32);
            }
            if ((nx as f64) > ((((sx as f32) + (temp * dzdx)) as f64) + (0.04f64 * (nx as f64)))) {
                *mbx2 = (((((sx as f32) + (temp * dzdx)) as f64) + (0.04f64 * (nx as f64))) as i32);
            } else {
                *mbx2 = nx;
            }
        } else {
            *mbz1 = 0i32;
            if ((nz as f64) > ((temp as f64) + (0.04f64 * (nz as f64)))) {
                *mbz2 = (((temp as f64) + (0.04f64 * (nz as f64))) as i32);
            } else {
                *mbz2 = nz;
            }
            *mbx1 = 0i32;
            *mbx2 = nx;
        }
    } else {
        *mbx1 = 0i32;
        *mbx2 = nx;
        *mbz1 = 0i32;
        *mbz2 = nz;
    }
}

pub(super) fn moving_fctbc(
    mut mbx1: i32,
    mut mbz1: i32,
    mut mbx2: i32,
    mut mbz2: i32,
    mut nxcc1: i32,
    mut nzcc1: i32,
    mut nxcc2: i32,
    mut nzcc2: i32,
    mut fctxbeg: &mut i32,
    mut fctzbeg: &mut i32,
    mut fctxend: &mut i32,
    mut fctzend: &mut i32,
) {
    if (mbx1 > nxcc1) {
        *fctxbeg = mbx1;
    } else {
        *fctxbeg = nxcc1;
    }
    if (mbx2 < nxcc2) {
        *fctxend = mbx2;
    } else {
        *fctxend = nxcc2;
    }
    if (mbz1 > nzcc1) {
        *fctzbeg = mbz1;
    } else {
        *fctzbeg = nzcc1;
    }
    if (mbz2 < nzcc2) {
        *fctzend = mbz2;
    } else {
        *fctzend = nzcc2;
    }
}

pub(super) fn difference(
    mut mbx1: i32,
    mut mbx2: i32,
    mut mbz1: i32,
    mut mbz2: i32,
    mut it: i32,
    mut u1: &mut Grid,
    mut u: &Grid,
    mut xzsource: &Grid,
    mut f: &Vec<f32>,
    mut dt: f32,
    mut rdxdz: f32,
) {
    let mut i: i32 = 0i32;
    let mut k: i32 = 0i32;
    i = (mbx1 + 1i32);
    while (i < (mbx2 - 1i32)) {
        {
            k = (mbz1 + 1i32);
            while (k < (mbz2 - 1i32)) {
                {
                    u1[(i) as usize][(k) as usize] = (u[(i) as usize][(k) as usize]
                        + (((f[(it) as usize] * xzsource[(i) as usize][(k) as usize]) * dt)
                            * rdxdz));
                }
                k += 1;
            }
        }
        i += 1;
    }
}

pub(super) fn difference_2x(
    mut mbx1: i32,
    mut mbx2: i32,
    mut mbz1: i32,
    mut mbz2: i32,
    mut shift: i32,
    mut dtdx: f32,
    mut u: &mut Grid,
    mut e: &Grid,
    mut c: &Grid,
) {
    let mut i: i32 = 0i32;
    let mut k: i32 = 0i32;
    i = (mbx1 + 1i32);
    while (i < (mbx2 - 1i32)) {
        {
            k = (mbz1 + 1i32);
            while (k < (mbz2 - 1i32)) {
                {
                    u[(i) as usize][(k) as usize] = (u[(i) as usize][(k) as usize]
                        + ((dtdx * c[(i) as usize][(k) as usize])
                            * (e[((i + 1i32) + shift) as usize][(k) as usize]
                                - e[(i + shift) as usize][(k) as usize])));
                }
                k += 1;
            }
        }
        i += 1;
    }
}

pub(super) fn difference_2z(
    mut mbx1: i32,
    mut mbx2: i32,
    mut mbz1: i32,
    mut mbz2: i32,
    mut shift: i32,
    mut dtdz: f32,
    mut u: &mut Grid,
    mut e: &Grid,
    mut c: &Grid,
) {
    let mut i: i32 = 0i32;
    let mut k: i32 = 0i32;
    i = (mbx1 + 1i32);
    while (i < (mbx2 - 1i32)) {
        {
            k = (mbz1 + 1i32);
            while (k < (mbz2 - 1i32)) {
                {
                    u[(i) as usize][(k) as usize] = (u[(i) as usize][(k) as usize]
                        + ((dtdz * c[(i) as usize][(k) as usize])
                            * (e[(i) as usize][((k + 1i32) + shift) as usize]
                                - e[(i) as usize][(k + shift) as usize])));
                }
                k += 1;
            }
        }
        i += 1;
    }
}
