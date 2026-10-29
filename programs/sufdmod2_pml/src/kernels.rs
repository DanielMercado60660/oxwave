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

// Copyright (c) Colorado School of Mines, 2011. All rights reserved.

//! Local finite-difference and PML routines from sufdmod2_pml.c (Revision 1.9).
//! Credits: Dave Hale; Michael Holzrichter, TAMU, 1998 (PML);
//! Greg Wimpey, CWP/WesternGeco, 2006 (density corrections).
#![allow(unused_variables, unused_assignments, unused_parens)]
use super::Grid;
use oxwave_core::prelude::*;
pub(super) fn star1(
    nx: i32,
    dx: f32,
    nz: i32,
    dz: f32,
    dt: f32,
    dvv: &Grid,
    od: &Grid,
    has_od: bool,
    s: &Grid,
    pm: &Grid,
    p: &Grid,
    pp: &mut Grid,
) {
    let mut ix: i32 = 0;
    let mut iz: i32 = 0;
    let mut xscale1: f32 = 0.0;
    let mut zscale1: f32 = 0.0;
    let mut xscale2: f32 = 0.0;
    let mut zscale2: f32 = 0.0;
    // determine constants
    xscale1 = ((dt * dt) / (dx * dx));
    zscale1 = ((dt * dt) / (dz * dz));
    xscale2 = ((0.25 * (xscale1 as f64)) as f32);
    zscale2 = ((0.25 * (zscale1 as f64)) as f32);
    // do the finite-difference star
    ix = 1;
    while (ix < (nx - 1)) {
        iz = 1;
        while (iz < (nz - 1)) {
            pp[(ix, iz)] = (((((2.0 * (p[(ix, iz)] as f64)) - (pm[(ix, iz)] as f64))
                + ((dvv[(ix, iz)] as f64)
                    * (((od[(ix, iz)] as f64)
                        * (((xscale1 as f64)
                            * (((p[((ix + 1), iz)] + p[((ix - 1), iz)]) as f64)
                                - (2.0 * (p[(ix, iz)] as f64))))
                            + ((zscale1 as f64)
                                * (((p[(ix, (iz + 1))] + p[(ix, (iz - 1))]) as f64)
                                    - (2.0 * (p[(ix, iz)] as f64))))))
                        + (((xscale2
                            * ((od[((ix + 1), iz)] - od[((ix - 1), iz)])
                                * (p[((ix + 1), iz)] - p[((ix - 1), iz)])))
                            + (zscale2
                                * ((od[(ix, (iz + 1))] - od[(ix, (iz - 1))])
                                    * (p[(ix, (iz + 1))] - p[(ix, (iz - 1))]))))
                            as f64))))
                + (s[(ix, iz)] as f64)) as f32);
            iz += 1;
        }
        ix += 1;
    }
}

pub(super) fn star2(
    nx: i32,
    dx: f32,
    nz: i32,
    dz: f32,
    dt: f32,
    dvv: &Grid,
    od: &Grid,
    has_od: bool,
    s: &Grid,
    pm: &Grid,
    p: &Grid,
    pp: &mut Grid,
) {
    let mut ix: i32 = 0;
    let mut iz: i32 = 0;
    let mut scale1: f32 = 0.0;
    let mut scale2: f32 = 0.0;
    if (dx != dz) {
        warn!("ASSERT FAILED: dx != dz in star2");
    }
    // determine constants
    scale1 = ((dt * dt) / (dx * dx));
    scale2 = ((0.25 * (scale1 as f64)) as f32);
    // do the finite-difference star
    ix = 1;
    while (ix < (nx - 1)) {
        iz = 1;
        while (iz < (nz - 1)) {
            pp[(ix, iz)] = (((((2.0 * (p[(ix, iz)] as f64)) - (pm[(ix, iz)] as f64))
                + ((dvv[(ix, iz)] as f64)
                    * (((od[(ix, iz)] as f64)
                        * ((scale1 as f64)
                            * (((((p[((ix + 1), iz)] + p[((ix - 1), iz)]) + p[(ix, (iz + 1))])
                                + p[(ix, (iz - 1))]) as f64)
                                - (4.0 * (p[(ix, iz)] as f64)))))
                        + ((scale2
                            * (((od[((ix + 1), iz)] - od[((ix - 1), iz)])
                                * (p[((ix + 1), iz)] - p[((ix - 1), iz)]))
                                + ((od[(ix, (iz + 1))] - od[(ix, (iz - 1))])
                                    * (p[(ix, (iz + 1))] - p[(ix, (iz - 1))]))))
                            as f64))))
                + (s[(ix, iz)] as f64)) as f32);
            iz += 1;
        }
        ix += 1;
    }
}

pub(super) fn star3(
    nx: i32,
    dx: f32,
    nz: i32,
    dz: f32,
    dt: f32,
    dvv: &Grid,
    od: &Grid,
    has_od: bool,
    s: &Grid,
    pm: &Grid,
    p: &Grid,
    pp: &mut Grid,
) {
    let mut ix: i32 = 0;
    let mut iz: i32 = 0;
    let mut xscale: f32 = 0.0;
    let mut zscale: f32 = 0.0;
    if has_od {
        warn!("ASSERT FAILED: od !=  NULL in star3");
    }
    // determine constants
    xscale = ((dt * dt) / (dx * dx));
    zscale = ((dt * dt) / (dz * dz));
    // do the finite-difference star
    ix = 1;
    while (ix < (nx - 1)) {
        iz = 1;
        while (iz < (nz - 1)) {
            pp[(ix, iz)] = (((((2.0 * (p[(ix, iz)] as f64)) - (pm[(ix, iz)] as f64))
                + ((dvv[(ix, iz)] as f64)
                    * (((xscale as f64)
                        * (((p[((ix + 1), iz)] + p[((ix - 1), iz)]) as f64)
                            - (2.0 * (p[(ix, iz)] as f64))))
                        + ((zscale as f64)
                            * (((p[(ix, (iz + 1))] + p[(ix, (iz - 1))]) as f64)
                                - (2.0 * (p[(ix, iz)] as f64)))))))
                + (s[(ix, iz)] as f64)) as f32);
            iz += 1;
        }
        ix += 1;
    }
}

pub(super) fn star4(
    nx: i32,
    dx: f32,
    nz: i32,
    dz: f32,
    dt: f32,
    dvv: &Grid,
    od: &Grid,
    has_od: bool,
    s: &Grid,
    pm: &Grid,
    p: &Grid,
    pp: &mut Grid,
) {
    let mut ix: i32 = 0;
    let mut iz: i32 = 0;
    let mut scale: f32 = 0.0;
    // determine constants
    if has_od {
        warn!("ASSERT FAILED: od !=  NULL in star4");
    }
    if (dz != dx) {
        warn!("ASSERT FAILED: dz !=  dx in star4");
    }
    scale = ((dt * dt) / (dx * dz));
    // do the finite-difference star
    ix = 1;
    while (ix < (nx - 1)) {
        iz = 1;
        while (iz < (nz - 1)) {
            pp[(ix, iz)] = (((((2.0 * (p[(ix, iz)] as f64)) - (pm[(ix, iz)] as f64))
                + (((scale * dvv[(ix, iz)]) as f64)
                    * (((((p[((ix + 1), iz)] + p[((ix - 1), iz)]) + p[(ix, (iz + 1))])
                        + p[(ix, (iz - 1))]) as f64)
                        - (4.0 * (p[(ix, iz)] as f64)))))
                + (s[(ix, iz)] as f64)) as f32);
            iz += 1;
        }
        ix += 1;
    }
}

pub(super) fn absorb(
    nx: i32,
    dx: f32,
    nz: i32,
    dz: f32,
    dt: f32,
    dvv: &Grid,
    od: &Grid,
    has_od: bool,
    pm: &Grid,
    p: &Grid,
    pp: &mut Grid,
    abs: &[i32; 4],
) {
    let mut ix: i32 = 0;
    let mut iz: i32 = 0;
    let mut ov: f32 = 0.0;
    let mut ovs: f32 = 0.0;
    let mut cosa: f32 = 0.0;
    let mut beta: f32 = 0.0;
    let mut gamma: f32 = 0.0;
    let mut dpdx: f32 = 0.0;
    let mut dpdz: f32 = 0.0;
    let mut dpdt: f32 = 0.0;
    let mut dpdxs: f32 = 0.0;
    let mut dpdzs: f32 = 0.0;
    let mut dpdts: f32 = 0.0;
    // solve for upper boundary
    iz = 1;
    ix = 0;
    while (ix < nx) {
        if (abs[0 as usize] != 0) {
            if has_od {
                ovs = ((1.0 / ((od[(ix, iz)] * dvv[(ix, iz)]) as f64)) as f32);
            } else {
                ovs = ((1.0 / (dvv[(ix, iz)] as f64)) as f32);
            }
            ov = ((ovs as f64).sqrt() as f32);
            if (ix == 0) {
                dpdx = ((p[(1, iz)] - p[(0, iz)]) / dx);
            } else {
                if (ix == (nx - 1)) {
                    dpdx = ((p[((nx - 1), iz)] - p[((nx - 2), iz)]) / dx);
                } else {
                    dpdx = ((((p[((ix + 1), iz)] - p[((ix - 1), iz)]) as f64) / (2.0 * (dx as f64)))
                        as f32);
                }
            }
            dpdt = ((((pp[(ix, iz)] - pm[(ix, iz)]) as f64) / (2.0 * (dt as f64))) as f32);
            dpdxs = (dpdx * dpdx);
            dpdts = (dpdt * dpdt);
            if ((ovs * dpdts) > dpdxs) {
                cosa = ((1.0 - ((dpdxs / (ovs * dpdts)) as f64)).sqrt() as f32);
            } else {
                cosa = (0.0 as f32);
            }
            beta = (((ov * dz) / dt) * cosa);
            gamma = (((1.0 - (beta as f64)) / (1.0 + (beta as f64))) as f32);
            pp[(ix, (iz - 1))] = ((gamma * (pp[(ix, iz)] - p[(ix, (iz - 1))])) + p[(ix, iz)]);
        } else {
            pp[(ix, (iz - 1))] = (0.0 as f32);
        }
        ix += 1;
    }
    // extrapolate along left boundary
    ix = 1;
    iz = 0;
    while (iz < nz) {
        if (abs[1 as usize] != 0) {
            if has_od {
                ovs = ((1.0 / ((od[(ix, iz)] * dvv[(ix, iz)]) as f64)) as f32);
            } else {
                ovs = ((1.0 / (dvv[(ix, iz)] as f64)) as f32);
            }
            ov = ((ovs as f64).sqrt() as f32);
            if (iz == 0) {
                dpdz = ((p[(ix, 1)] - p[(ix, 0)]) / dz);
            } else {
                if (iz == (nz - 1)) {
                    dpdz = ((p[(ix, (nz - 1))] - p[(ix, (nz - 2))]) / dz);
                } else {
                    dpdz = ((((p[(ix, (iz + 1))] - p[(ix, (iz - 1))]) as f64) / (2.0 * (dz as f64)))
                        as f32);
                }
            }
            dpdt = ((((pp[(ix, iz)] - pm[(ix, iz)]) as f64) / (2.0 * (dt as f64))) as f32);
            dpdzs = (dpdz * dpdz);
            dpdts = (dpdt * dpdt);
            if ((ovs * dpdts) > dpdzs) {
                cosa = ((1.0 - ((dpdzs / (ovs * dpdts)) as f64)).sqrt() as f32);
            } else {
                cosa = (0.0 as f32);
            }
            beta = (((ov * dx) / dt) * cosa);
            gamma = (((1.0 - (beta as f64)) / (1.0 + (beta as f64))) as f32);
            pp[((ix - 1), iz)] = ((gamma * (pp[(ix, iz)] - p[((ix - 1), iz)])) + p[(ix, iz)]);
        } else {
            pp[((ix - 1), iz)] = (0.0 as f32);
        }
        iz += 1;
    }
    // extrapolate along lower boundary
    iz = (nz - 2);
    ix = 0;
    while (ix < nx) {
        if (abs[2 as usize] != 0) {
            if has_od {
                ovs = ((1.0 / ((od[(ix, iz)] * dvv[(ix, iz)]) as f64)) as f32);
            } else {
                ovs = ((1.0 / (dvv[(ix, iz)] as f64)) as f32);
            }
            ov = ((ovs as f64).sqrt() as f32);
            if (ix == 0) {
                dpdx = ((p[(1, iz)] - p[(0, iz)]) / dx);
            } else {
                if (ix == (nx - 1)) {
                    dpdx = ((p[((nx - 1), iz)] - p[((nx - 2), iz)]) / dx);
                } else {
                    dpdx = ((((p[((ix + 1), iz)] - p[((ix - 1), iz)]) as f64) / (2.0 * (dx as f64)))
                        as f32);
                }
            }
            dpdt = ((((pp[(ix, iz)] - pm[(ix, iz)]) as f64) / (2.0 * (dt as f64))) as f32);
            dpdxs = (dpdx * dpdx);
            dpdts = (dpdt * dpdt);
            if ((ovs * dpdts) > dpdxs) {
                cosa = ((1.0 - ((dpdxs / (ovs * dpdts)) as f64)).sqrt() as f32);
            } else {
                cosa = (0.0 as f32);
            }
            beta = (((ov * dz) / dt) * cosa);
            gamma = (((1.0 - (beta as f64)) / (1.0 + (beta as f64))) as f32);
            pp[(ix, (iz + 1))] = ((gamma * (pp[(ix, iz)] - p[(ix, (iz + 1))])) + p[(ix, iz)]);
        } else {
            pp[(ix, (iz + 1))] = (0.0 as f32);
        }
        ix += 1;
    }
    // extrapolate along right boundary
    ix = (nx - 2);
    iz = 0;
    while (iz < nz) {
        if (abs[3 as usize] != 0) {
            if has_od {
                ovs = ((1.0 / ((od[(ix, iz)] * dvv[(ix, iz)]) as f64)) as f32);
            } else {
                ovs = ((1.0 / (dvv[(ix, iz)] as f64)) as f32);
            }
            ov = ((ovs as f64).sqrt() as f32);
            if (iz == 0) {
                dpdz = ((p[(ix, 1)] - p[(ix, 0)]) / dz);
            } else {
                if (iz == (nz - 1)) {
                    dpdz = ((p[(ix, (nz - 1))] - p[(ix, (nz - 2))]) / dz);
                } else {
                    dpdz = ((((p[(ix, (iz + 1))] - p[(ix, (iz - 1))]) as f64) / (2.0 * (dz as f64)))
                        as f32);
                }
            }
            dpdt = ((((pp[(ix, iz)] - pm[(ix, iz)]) as f64) / (2.0 * (dt as f64))) as f32);
            dpdzs = (dpdz * dpdz);
            dpdts = (dpdt * dpdt);
            if ((ovs * dpdts) > dpdzs) {
                cosa = ((1.0 - ((dpdzs / (ovs * dpdts)) as f64)).sqrt() as f32);
            } else {
                cosa = (0.0 as f32);
            }
            beta = (((ov * dx) / dt) * cosa);
            gamma = (((1.0 - (beta as f64)) / (1.0 + (beta as f64))) as f32);
            pp[((ix + 1), iz)] = ((gamma * (pp[(ix, iz)] - p[((ix + 1), iz)])) + p[(ix, iz)]);
        } else {
            pp[((ix + 1), iz)] = (0.0 as f32);
        }
        iz += 1;
    }
}

pub(super) struct Pml {
    pub cax_b: Grid,
    pub cax_r: Grid,
    pub cbx_b: Grid,
    pub cbx_r: Grid,
    pub caz_b: Grid,
    pub caz_r: Grid,
    pub cbz_b: Grid,
    pub cbz_r: Grid,
    pub dax_b: Grid,
    pub dax_r: Grid,
    pub dbx_b: Grid,
    pub dbx_r: Grid,
    pub daz_b: Grid,
    pub daz_r: Grid,
    pub dbz_b: Grid,
    pub dbz_r: Grid,
    pub ux_b: Grid,
    pub ux_r: Grid,
    pub uz_b: Grid,
    pub uz_r: Grid,
    pub v_b: Grid,
    pub v_r: Grid,
    pub w_b: Grid,
    pub w_r: Grid,
    pub pml_max: f32,
    pub pml_thickness: i32,
}
impl Pml {
    pub(super) fn new(pml_max: f32, pml_thickness: i32) -> Self {
        Self {
            pml_max,
            pml_thickness,
            cax_b: Grid::new(0, 0),
            cax_r: Grid::new(0, 0),
            cbx_b: Grid::new(0, 0),
            cbx_r: Grid::new(0, 0),
            caz_b: Grid::new(0, 0),
            caz_r: Grid::new(0, 0),
            cbz_b: Grid::new(0, 0),
            cbz_r: Grid::new(0, 0),
            dax_b: Grid::new(0, 0),
            dax_r: Grid::new(0, 0),
            dbx_b: Grid::new(0, 0),
            dbx_r: Grid::new(0, 0),
            daz_b: Grid::new(0, 0),
            daz_r: Grid::new(0, 0),
            dbz_b: Grid::new(0, 0),
            dbz_r: Grid::new(0, 0),
            ux_b: Grid::new(0, 0),
            ux_r: Grid::new(0, 0),
            uz_b: Grid::new(0, 0),
            uz_r: Grid::new(0, 0),
            v_b: Grid::new(0, 0),
            v_r: Grid::new(0, 0),
            w_b: Grid::new(0, 0),
            w_r: Grid::new(0, 0),
        }
    }
    pub(super) fn pml_init(
        &mut self,
        nx: i32,
        nz: i32,
        dx: f32,
        dz: f32,
        dt: f32,
        dvv: &Grid,
        od: &Grid,
        has_od: bool,
        verbose: i32,
    ) {
        let mut dvv_0: f32 = 0.0;
        let mut dvv_1: f32 = 0.0;
        let mut dvv_2: f32 = 0.0;
        let mut dvv_3: f32 = 0.0;
        let mut sigma_ex: f32 = 0.0;
        let mut sigma_ez: f32 = 0.0;
        let mut sigma_mx: f32 = 0.0;
        let mut sigma_mz: f32 = 0.0;
        let mut ix: i32 = 0;
        let mut iz: i32 = 0;
        // Allocate arrays for pad on right
        self.cax_r = Grid::new(nz, (self.pml_thickness + 2));
        self.cbx_r = Grid::new(nz, (self.pml_thickness + 2));
        self.caz_r = Grid::new(nz, (self.pml_thickness + 3));
        self.cbz_r = Grid::new(nz, (self.pml_thickness + 3));
        self.dax_r = Grid::new(nz, (self.pml_thickness + 2));
        self.dbx_r = Grid::new(nz, (self.pml_thickness + 2));
        self.daz_r = Grid::new(nz, (self.pml_thickness + 2));
        self.dbz_r = Grid::new(nz, (self.pml_thickness + 2));
        self.ux_r = Grid::new(nz, (self.pml_thickness + 2));
        self.uz_r = Grid::new(nz, (self.pml_thickness + 2));
        self.v_r = Grid::new(nz, (self.pml_thickness + 2));
        self.w_r = Grid::new(nz, (self.pml_thickness + 3));
        // Zero out arrays for pad on right
        ix = 0;
        while (ix < (self.pml_thickness + 2)) {
            iz = 0;
            while (iz < nz) {
                self.uz_r[(ix, iz)] = (0.0 as f32);
                self.ux_r[(ix, iz)] = self.uz_r[(ix, iz)];
                self.w_r[(ix, iz)] = (0.0 as f32);
                self.v_r[(ix, iz)] = self.w_r[(ix, iz)];
                self.cbx_r[(ix, iz)] = (0.0 as f32);
                self.cax_r[(ix, iz)] = self.cbx_r[(ix, iz)];
                self.cbz_r[(ix, iz)] = (0.0 as f32);
                self.caz_r[(ix, iz)] = self.cbz_r[(ix, iz)];
                self.dbx_r[(ix, iz)] = (0.0 as f32);
                self.dax_r[(ix, iz)] = self.dbx_r[(ix, iz)];
                self.dbz_r[(ix, iz)] = (0.0 as f32);
                self.daz_r[(ix, iz)] = self.dbz_r[(ix, iz)];
                iz += 1;
            }
            ix += 1;
        }
        // Zero out extra bit on right pad
        ix = (self.pml_thickness + 2);
        iz = 0;
        while (iz < nz) {
            self.cbz_r[(ix, iz)] = (0.0 as f32);
            self.caz_r[(ix, iz)] = self.cbz_r[(ix, iz)];
            self.w_r[(ix, iz)] = (0.0 as f32);
            iz += 1;
        }
        // Allocate arrays for pad on bottom
        self.cax_b = Grid::new((self.pml_thickness + 3), (nx + self.pml_thickness));
        self.cbx_b = Grid::new((self.pml_thickness + 3), (nx + self.pml_thickness));
        self.caz_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        self.cbz_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        self.dax_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        self.dbx_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        self.daz_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        self.dbz_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        self.ux_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        self.uz_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        self.v_b = Grid::new((self.pml_thickness + 3), (nx + self.pml_thickness));
        self.w_b = Grid::new((self.pml_thickness + 2), (nx + self.pml_thickness));
        // Zero out arrays for pad on bottom
        ix = 0;
        while (ix < (nx + self.pml_thickness)) {
            iz = 0;
            while (iz < (self.pml_thickness + 2)) {
                self.uz_b[(ix, iz)] = (0.0 as f32);
                self.ux_b[(ix, iz)] = self.uz_b[(ix, iz)];
                self.w_b[(ix, iz)] = (0.0 as f32);
                self.v_b[(ix, iz)] = self.w_b[(ix, iz)];
                self.cbx_b[(ix, iz)] = (0.0 as f32);
                self.cax_b[(ix, iz)] = self.cbx_b[(ix, iz)];
                self.cbz_b[(ix, iz)] = (0.0 as f32);
                self.caz_b[(ix, iz)] = self.cbz_b[(ix, iz)];
                self.dbx_b[(ix, iz)] = (0.0 as f32);
                self.dax_b[(ix, iz)] = self.dbx_b[(ix, iz)];
                self.dbz_b[(ix, iz)] = (0.0 as f32);
                self.daz_b[(ix, iz)] = self.dbz_b[(ix, iz)];
                iz += 1;
            }
            ix += 1;
        }
        // Zero out extra bit on bottom pad
        ix = 0;
        iz = (self.pml_thickness + 2);
        while (ix < (nx + self.pml_thickness)) {
            self.cbx_b[(ix, iz)] = (0.0 as f32);
            self.cax_b[(ix, iz)] = self.cbx_b[(ix, iz)];
            self.v_b[(ix, iz)] = (0.0 as f32);
            ix += 1;
        }
        // Initialize cax & cbx arrays
        ix = 0;
        while (ix < (self.pml_thickness + 2)) {
            iz = 0;
            while (iz < nz) {
                sigma_ez = (0.0 as f32);
                self.cax_r[(ix, iz)] =
                    (((2.0 - ((sigma_ez * dt) as f64)) / (2.0 + ((sigma_ez * dt) as f64))) as f32);
                self.cbx_r[(ix, iz)] = ((((2.0 * (dt as f64)) / (dz as f64))
                    / (2.0 + ((sigma_ez * dt) as f64)))
                    as f32);
                iz += 1;
            }
            ix += 1;
        }
        ix = 0;
        while (ix < (nx + self.pml_thickness)) {
            iz = 0;
            while (iz < (self.pml_thickness + 3)) {
                if ((iz == 0) || (iz == (self.pml_thickness + 2))) {
                    sigma_ez = (0.0 as f32);
                } else {
                    sigma_ez = ((((self.pml_max as f64) * 0.5)
                        * (1.0
                            - ((((2 as f64) * std::f64::consts::PI) * ((iz as f64) - 0.5))
                                / ((self.pml_thickness + 1) as f64))
                                .cos())) as f32);
                }
                self.cax_b[(ix, iz)] =
                    (((2.0 - ((sigma_ez * dt) as f64)) / (2.0 + ((sigma_ez * dt) as f64))) as f32);
                self.cbx_b[(ix, iz)] = ((((2.0 * (dt as f64)) / (dz as f64))
                    / (2.0 + ((sigma_ez * dt) as f64)))
                    as f32);
                iz += 1;
            }
            ix += 1;
        }
        // Initialize caz & cbz arrays
        ix = 0;
        while (ix < (self.pml_thickness + 3)) {
            iz = 0;
            while (iz < nz) {
                if ((ix == 0) || (ix == (self.pml_thickness + 2))) {
                    sigma_ex = (0.0 as f32);
                } else {
                    sigma_ex = ((((self.pml_max as f64) * 0.5)
                        * (1.0
                            - ((((2 as f64) * std::f64::consts::PI) * ((ix as f64) - 0.5))
                                / ((self.pml_thickness + 1) as f64))
                                .cos())) as f32);
                }
                self.caz_r[(ix, iz)] =
                    (((2.0 - ((sigma_ex * dt) as f64)) / (2.0 + ((sigma_ex * dt) as f64))) as f32);
                self.cbz_r[(ix, iz)] = ((((2.0 * (dt as f64)) / (dz as f64))
                    / (2.0 + ((sigma_ex * dt) as f64)))
                    as f32);
                iz += 1;
            }
            ix += 1;
        }
        ix = 0;
        while (ix < (nx + self.pml_thickness)) {
            iz = 0;
            while (iz < (self.pml_thickness + 2)) {
                if (ix == 0) {
                    sigma_ex = ((((self.pml_max as f64) * 0.5)
                        * (1.0
                            - ((((2 as f64) * std::f64::consts::PI) * 0.5)
                                / ((self.pml_thickness + 1) as f64))
                                .cos())) as f32);
                } else {
                    if (ix < nx) {
                        sigma_ex = (0.0 as f32);
                    } else {
                        sigma_ex = ((((self.pml_max as f64) * 0.5)
                            * (1.0
                                - ((((2 as f64) * std::f64::consts::PI)
                                    * (((ix - nx) as f64) + 0.5))
                                    / ((self.pml_thickness + 1) as f64))
                                    .cos())) as f32);
                    }
                }
                self.caz_b[(ix, iz)] =
                    (((2.0 - ((sigma_ex * dt) as f64)) / (2.0 + ((sigma_ex * dt) as f64))) as f32);
                self.cbz_b[(ix, iz)] = ((((2.0 * (dt as f64)) / (dz as f64))
                    / (2.0 + ((sigma_ex * dt) as f64)))
                    as f32);
                iz += 1;
            }
            ix += 1;
        }
        // Initialize right pad's dax dbx, daz, & dbz arrays
        ix = 0;
        while (ix < (self.pml_thickness + 2)) {
            iz = 0;
            while (iz < nz) {
                // Determine sigma_mx and sigma_mz
                if ((ix == 0) || (ix == (self.pml_thickness + 1))) {
                    sigma_mx = (0.0 as f32);
                } else {
                    sigma_mx = ((((self.pml_max as f64) * 0.5)
                        * (1.0
                            - ((((2 as f64) * std::f64::consts::PI) * (ix as f64))
                                / ((self.pml_thickness + 1) as f64))
                                .cos())) as f32);
                }
                sigma_mz = (0.0 as f32);
                // Determine velocity, interpolate
                if has_od {
                    if (ix == 0) {
                        dvv_0 = (((dvv[((nx - 1), iz)] * od[((nx - 1), iz)]) as f64).sqrt() as f32);
                    } else {
                        if (ix == (self.pml_thickness + 1)) {
                            dvv_0 = (((dvv[(0, iz)] * od[(0, iz)]) as f64).sqrt() as f32);
                        } else {
                            dvv_0 =
                                (((dvv[((nx - 1), iz)] * od[((nx - 1), iz)]) as f64).sqrt() as f32);
                            dvv_1 = (((dvv[(0, iz)] * od[(0, iz)]) as f64).sqrt() as f32);
                            dvv_0 = (((ix as f32) * dvv_1)
                                + ((((1 + self.pml_thickness) - ix) as f32) * dvv_0));
                            dvv_0 = (dvv_0 / ((self.pml_thickness + 1) as f32));
                        }
                    }
                } else {
                    if (ix == 0) {
                        dvv_0 = ((dvv[((nx - 1), iz)] as f64).sqrt() as f32);
                    } else {
                        if (ix == (self.pml_thickness + 1)) {
                            dvv_0 = ((dvv[(0, iz)] as f64).sqrt() as f32);
                        } else {
                            dvv_0 = ((dvv[((nx - 1), iz)] as f64).sqrt() as f32);
                            dvv_1 = ((dvv[(0, iz)] as f64).sqrt() as f32);
                            dvv_0 = (((ix as f32) * dvv_1)
                                + ((((1 + self.pml_thickness) - ix) as f32) * dvv_0));
                            dvv_0 = (dvv_0 / ((self.pml_thickness + 1) as f32));
                        }
                    }
                }
                dvv_0 = (dvv_0 * dvv_0);
                self.dax_r[(ix, iz)] =
                    (((2.0 - ((sigma_mx * dt) as f64)) / (2.0 + ((sigma_mx * dt) as f64))) as f32);
                self.dbx_r[(ix, iz)] = (((((2.0 * (dt as f64)) * (dvv_0 as f64)) / (dx as f64))
                    / (2.0 + ((sigma_mx * dt) as f64)))
                    as f32);
                self.daz_r[(ix, iz)] =
                    (((2.0 - ((sigma_mz * dt) as f64)) / (2.0 + ((sigma_mz * dt) as f64))) as f32);
                self.dbz_r[(ix, iz)] = (((((2.0 * (dt as f64)) * (dvv_0 as f64)) / (dz as f64))
                    / (2.0 + ((sigma_mz * dt) as f64)))
                    as f32);
                iz += 1;
            }
            ix += 1;
        }
        // Initialize bottom pad's dax, dbx, daz, & dbz arrays
        ix = 0;
        while (ix < (nx + self.pml_thickness)) {
            iz = 0;
            while (iz < (self.pml_thickness + 2)) {
                // Determine sigma_mx and sigma_mz
                if (ix < nx) {
                    sigma_mx = (0.0 as f32);
                } else {
                    sigma_mx = ((((self.pml_max as f64) * 0.5)
                        * (1.0
                            - ((((2 as f64) * std::f64::consts::PI) * (((ix - nx) + 1) as f64))
                                / ((self.pml_thickness + 1) as f64))
                                .cos())) as f32);
                }
                if ((iz == 0) || (iz == (self.pml_thickness + 1))) {
                    sigma_mz = (0.0 as f32);
                } else {
                    sigma_mz = ((((self.pml_max as f64) * 0.5)
                        * (1.0
                            - ((((2 as f64) * std::f64::consts::PI) * (iz as f64))
                                / ((self.pml_thickness + 1) as f64))
                                .cos())) as f32);
                }
                // Determine velocity, interpolate
                if has_od {
                    if (ix < nx) {
                        if (iz == 0) {
                            dvv_0 =
                                (((dvv[(ix, (nz - 1))] * od[(ix, (nz - 1))]) as f64).sqrt() as f32);
                        } else {
                            if (iz == (self.pml_thickness + 1)) {
                                dvv_0 = (((dvv[(ix, 0)] * od[(ix, 0)]) as f64).sqrt() as f32);
                            } else {
                                dvv_0 = (((dvv[(ix, (nz - 1))] * od[(ix, (nz - 1))]) as f64).sqrt()
                                    as f32);
                                dvv_1 = (((dvv[(ix, 0)] * od[(ix, 0)]) as f64).sqrt() as f32);
                                dvv_0 = (((iz as f32) * dvv_1)
                                    + ((((1 + self.pml_thickness) - iz) as f32) * dvv_0));
                                dvv_0 = (dvv_0 / ((self.pml_thickness + 1) as f32));
                            }
                        }
                    } else {
                        if (iz == 0) {
                            dvv_0 = (((dvv[((nx - 1), (nz - 1))] * od[((nx - 1), (nz - 1))]) as f64)
                                .sqrt() as f32);
                            dvv_1 =
                                (((dvv[(0, (nz - 1))] * od[(0, (nz - 1))]) as f64).sqrt() as f32);
                        } else {
                            if (iz == (self.pml_thickness + 1)) {
                                dvv_0 = (((dvv[((nx - 1), 0)] * od[((nx - 1), 0)]) as f64).sqrt()
                                    as f32);
                                dvv_1 = (((dvv[(0, 0)] * od[(0, 0)]) as f64).sqrt() as f32);
                            } else {
                                dvv_2 = (((dvv[((nx - 1), (nz - 1))] * od[((nx - 1), (nz - 1))])
                                    as f64)
                                    .sqrt() as f32);
                                dvv_3 = (((dvv[((nx - 1), 0)] * od[((nx - 1), 0)]) as f64).sqrt()
                                    as f32);
                                dvv_0 = (((iz as f32) * dvv_3)
                                    + ((((1 + self.pml_thickness) - iz) as f32) * dvv_2));
                                dvv_0 = (dvv_0 / ((self.pml_thickness + 1) as f32));
                                dvv_2 = (((dvv[(0, (nz - 1))] * od[(0, (nz - 1))]) as f64).sqrt()
                                    as f32);
                                dvv_3 = (((dvv[(0, 0)] * od[(0, 0)]) as f64).sqrt() as f32);
                                dvv_1 = (((iz as f32) * dvv_3)
                                    + ((((1 + self.pml_thickness) - iz) as f32) * dvv_2));
                                dvv_1 = (dvv_1 / ((self.pml_thickness + 1) as f32));
                            }
                        }
                        dvv_0 = (((((ix - nx) + 1) as f32) * dvv_1)
                            + ((((nx + self.pml_thickness) - ix) as f32) * dvv_0));
                        dvv_0 = (dvv_0 / ((self.pml_thickness + 1) as f32));
                    }
                } else {
                    if (ix < nx) {
                        if (iz == 0) {
                            dvv_0 = ((dvv[(ix, (nz - 1))] as f64).sqrt() as f32);
                        } else {
                            if (iz == (self.pml_thickness + 1)) {
                                dvv_0 = ((dvv[(ix, 0)] as f64).sqrt() as f32);
                            } else {
                                dvv_0 = ((dvv[(ix, (nz - 1))] as f64).sqrt() as f32);
                                dvv_1 = ((dvv[(ix, 0)] as f64).sqrt() as f32);
                                dvv_0 = (((iz as f32) * dvv_1)
                                    + ((((1 + self.pml_thickness) - iz) as f32) * dvv_0));
                                dvv_0 = (dvv_0 / ((self.pml_thickness + 1) as f32));
                            }
                        }
                    } else {
                        if (iz == 0) {
                            dvv_0 = ((dvv[((nx - 1), (nz - 1))] as f64).sqrt() as f32);
                            dvv_1 = ((dvv[(0, (nz - 1))] as f64).sqrt() as f32);
                        } else {
                            if (iz == (self.pml_thickness + 1)) {
                                dvv_0 = ((dvv[((nx - 1), 0)] as f64).sqrt() as f32);
                                dvv_1 = ((dvv[(0, 0)] as f64).sqrt() as f32);
                            } else {
                                dvv_2 = ((dvv[((nx - 1), (nz - 1))] as f64).sqrt() as f32);
                                dvv_3 = ((dvv[((nx - 1), 0)] as f64).sqrt() as f32);
                                dvv_0 = (((iz as f32) * dvv_3)
                                    + ((((1 + self.pml_thickness) - iz) as f32) * dvv_2));
                                dvv_0 = (dvv_0 / ((self.pml_thickness + 1) as f32));
                                dvv_2 = ((dvv[(0, (nz - 1))] as f64).sqrt() as f32);
                                dvv_3 = ((dvv[(0, 0)] as f64).sqrt() as f32);
                                dvv_1 = (((iz as f32) * dvv_3)
                                    + ((((1 + self.pml_thickness) - iz) as f32) * dvv_2));
                                dvv_1 = (dvv_1 / ((self.pml_thickness + 1) as f32));
                            }
                        }
                        dvv_0 = (((((ix - nx) + 1) as f32) * dvv_1)
                            + ((((nx + self.pml_thickness) - ix) as f32) * dvv_0));
                        dvv_0 = (dvv_0 / ((self.pml_thickness + 1) as f32));
                    }
                }
                dvv_0 = (dvv_0 * dvv_0);
                self.dax_b[(ix, iz)] =
                    (((2.0 - ((sigma_mx * dt) as f64)) / (2.0 + ((sigma_mx * dt) as f64))) as f32);
                self.dbx_b[(ix, iz)] = (((((2.0 * (dt as f64)) * (dvv_0 as f64)) / (dx as f64))
                    / (2.0 + ((sigma_mx * dt) as f64)))
                    as f32);
                self.daz_b[(ix, iz)] =
                    (((2.0 - ((sigma_mz * dt) as f64)) / (2.0 + ((sigma_mz * dt) as f64))) as f32);
                self.dbz_b[(ix, iz)] = (((((2.0 * (dt as f64)) * (dvv_0 as f64)) / (dz as f64))
                    / (2.0 + ((sigma_mz * dt) as f64)))
                    as f32);
                iz += 1;
            }
            ix += 1;
        }
    }
    pub(super) fn pml_absorb(
        &mut self,
        nx: i32,
        dx: f32,
        nz: i32,
        dz: f32,
        dt: f32,
        dvv: &Grid,
        od: &Grid,
        has_od: bool,
        pm: &Grid,
        p: &Grid,
        pp: &mut Grid,
        abs: &[i32; 4],
    ) {
        let mut ix: i32 = 0;
        let mut iz: i32 = 0;
        let mut jx: i32 = 0;
        let mut jz: i32 = 0;
        // Calculate v for bottom pad above and below main domain
        ix = 0;
        jz = (self.pml_thickness + 2);
        while (ix < nx) {
            self.v_b[(ix, 0)] = ((((self.cax_b[(ix, 0)] * self.v_b[(ix, 0)]) as f64)
                + ((self.cbx_b[(ix, 0)] as f64)
                    * (((self.ux_b[(ix, 0)] + self.uz_b[(ix, 0)]) as f64)
                        - (if (abs[2 as usize] != 0) {
                            (p[(ix, (nz - 2))] as f64)
                        } else {
                            0.0
                        })))) as f32);
            self.v_b[(ix, jz)] = ((((self.cax_b[(ix, jz)] * self.v_b[(ix, jz)]) as f64)
                + ((self.cbx_b[(ix, jz)] as f64)
                    * (((if (abs[0 as usize] != 0) {
                        (p[(ix, 1)] as f64)
                    } else {
                        0.0
                    }) - (self.ux_b[(ix, (jz - 1))] as f64))
                        - (self.uz_b[(ix, (jz - 1))] as f64))))
                as f32);
            ix += 1;
        }
        // Calculate v for bottom pad above and below right pad
        ix = nx;
        jx = 1;
        jz = (self.pml_thickness + 2);
        while (ix < (nx + self.pml_thickness)) {
            self.v_b[(ix, 0)] = ((self.cax_b[(ix, 0)] * self.v_b[(ix, 0)])
                + (self.cbx_b[(ix, 0)]
                    * (((self.ux_b[(ix, 0)] + self.uz_b[(ix, 0)]) - self.ux_r[(jx, (nz - 2))])
                        - self.uz_r[(jx, (nz - 2))])));
            self.v_b[(ix, jz)] = ((self.cax_b[(ix, jz)] * self.v_b[(ix, jz)])
                + (self.cbx_b[(ix, jz)]
                    * (((self.ux_r[(jx, 1)] + self.uz_r[(jx, 1)]) - self.ux_b[(ix, (jz - 1))])
                        - self.uz_b[(ix, (jz - 1))])));
            ix += 1;
            jx += 1;
        }
        // Calculate v for main part of bottom pad
        ix = 0;
        while (ix < (nx + self.pml_thickness)) {
            iz = 1;
            while (iz < (self.pml_thickness + 2)) {
                self.v_b[(ix, iz)] = ((self.cax_b[(ix, iz)] * self.v_b[(ix, iz)])
                    + (self.cbx_b[(ix, iz)]
                        * (((self.ux_b[(ix, iz)] + self.uz_b[(ix, iz)])
                            - self.ux_b[(ix, (iz - 1))])
                            - self.uz_b[(ix, (iz - 1))])));
                iz += 1;
            }
            ix += 1;
        }
        // Calculate w for left edge of bottom pad
        iz = 0;
        ix = ((nx + self.pml_thickness) - 1);
        while (iz < (self.pml_thickness + 2)) {
            self.w_b[(0, iz)] = ((self.caz_b[(0, iz)] * self.w_b[(0, iz)])
                + (self.cbz_b[(0, iz)]
                    * (((self.ux_b[(ix, iz)] + self.uz_b[(ix, iz)]) - self.ux_b[(0, iz)])
                        - self.uz_b[(0, iz)])));
            iz += 1;
        }
        // Calculate w for main part of bottom pad
        ix = 1;
        while (ix < (nx + self.pml_thickness)) {
            iz = 0;
            while (iz < (self.pml_thickness + 2)) {
                self.w_b[(ix, iz)] = ((self.caz_b[(ix, iz)] * self.w_b[(ix, iz)])
                    + (self.cbz_b[(ix, iz)]
                        * (((self.ux_b[((ix - 1), iz)] + self.uz_b[((ix - 1), iz)])
                            - self.ux_b[(ix, iz)])
                            - self.uz_b[(ix, iz)])));
                iz += 1;
            }
            ix += 1;
        }
        // Calculate v along top and bottom edge of right pad
        ix = 0;
        jx = (nx - 1);
        jz = self.pml_thickness;
        while (ix < (self.pml_thickness + 2)) {
            if (jx == (nx + self.pml_thickness)) {
                jx = 0;
            }
            self.v_r[(ix, 0)] = ((self.cax_r[(ix, 0)] * self.v_r[(ix, 0)])
                + (self.cbx_r[(ix, 0)]
                    * (((self.ux_r[(ix, 0)] + self.uz_r[(ix, 0)]) - self.ux_b[(jx, jz)])
                        - self.uz_b[(jx, jz)])));
            self.v_r[(ix, (nz - 1))] = ((self.cax_r[(ix, (nz - 1))] * self.v_r[(ix, (nz - 1))])
                + (self.cbx_r[(ix, (nz - 1))]
                    * (((self.ux_b[(jx, 0)] + self.uz_b[(jx, 0)]) - self.ux_r[(ix, (nz - 2))])
                        - self.uz_r[(ix, (nz - 2))])));
            ix += 1;
            jx += 1;
        }
        // Calculate v in rest of right pad
        ix = 0;
        while (ix < (self.pml_thickness + 2)) {
            iz = 1;
            while (iz < (nz - 1)) {
                self.v_r[(ix, iz)] = ((self.cax_r[(ix, iz)] * self.v_r[(ix, iz)])
                    + (self.cbx_r[(ix, iz)]
                        * (((self.ux_r[(ix, iz)] + self.uz_r[(ix, iz)])
                            - self.ux_r[(ix, (iz - 1))])
                            - self.uz_r[(ix, (iz - 1))])));
                iz += 1;
            }
            ix += 1;
        }
        // Calculate w along left and right sides of right pad
        iz = 0;
        jx = (self.pml_thickness + 2);
        while (iz < nz) {
            self.w_r[(0, iz)] = ((((self.caz_r[(0, iz)] * self.w_r[(0, iz)]) as f64)
                + ((self.cbz_r[(0, iz)] as f64)
                    * (((if (abs[3 as usize] != 0) {
                        (p[((nx - 2), iz)] as f64)
                    } else {
                        0.0
                    }) - (self.ux_r[(0, iz)] as f64))
                        - (self.uz_r[(0, iz)] as f64)))) as f32);
            self.w_r[(jx, iz)] = ((((self.caz_r[(jx, iz)] * self.w_r[(jx, iz)]) as f64)
                + ((self.cbz_r[(jx, iz)] as f64)
                    * (((self.ux_r[((jx - 1), iz)] + self.uz_r[((jx - 1), iz)]) as f64)
                        - (if (abs[1 as usize] != 0) {
                            (p[(1, iz)] as f64)
                        } else {
                            0.0
                        })))) as f32);
            iz += 1;
        }
        // Calculate w in main part of right pad
        ix = 1;
        while (ix < (self.pml_thickness + 2)) {
            iz = 0;
            while (iz < nz) {
                self.w_r[(ix, iz)] = ((self.caz_r[(ix, iz)] * self.w_r[(ix, iz)])
                    + (self.cbz_r[(ix, iz)]
                        * (((self.ux_r[((ix - 1), iz)] + self.uz_r[((ix - 1), iz)])
                            - self.ux_r[(ix, iz)])
                            - self.uz_r[(ix, iz)])));
                iz += 1;
            }
            ix += 1;
        }
        // Calculate ux and uz in bottom pad
        ix = 0;
        while (ix < ((nx + self.pml_thickness) - 1)) {
            iz = 0;
            while (iz < (self.pml_thickness + 2)) {
                self.ux_b[(ix, iz)] = ((self.dax_b[(ix, iz)] * self.ux_b[(ix, iz)])
                    + (self.dbx_b[(ix, iz)] * (self.w_b[(ix, iz)] - self.w_b[((ix + 1), iz)])));
                iz += 1;
            }
            ix += 1;
        }
        ix = ((nx + self.pml_thickness) - 1);
        iz = 0;
        while (iz < (self.pml_thickness + 2)) {
            self.ux_b[(ix, iz)] = ((self.dax_b[(ix, iz)] * self.ux_b[(ix, iz)])
                + (self.dbx_b[(ix, iz)] * (self.w_b[(ix, iz)] - self.w_b[(0, iz)])));
            iz += 1;
        }
        ix = 0;
        while (ix < (nx + self.pml_thickness)) {
            iz = 0;
            while (iz < (self.pml_thickness + 2)) {
                self.uz_b[(ix, iz)] = ((self.daz_b[(ix, iz)] * self.uz_b[(ix, iz)])
                    + (self.dbz_b[(ix, iz)] * (self.v_b[(ix, (iz + 1))] - self.v_b[(ix, iz)])));
                iz += 1;
            }
            ix += 1;
        }
        // Calculate ux and uz in right pad
        ix = 0;
        while (ix < (self.pml_thickness + 2)) {
            iz = 0;
            while (iz < nz) {
                self.ux_r[(ix, iz)] = ((self.dax_r[(ix, iz)] * self.ux_r[(ix, iz)])
                    + (self.dbx_r[(ix, iz)] * (self.w_r[(ix, iz)] - self.w_r[((ix + 1), iz)])));
                iz += 1;
            }
            ix += 1;
        }
        ix = 0;
        while (ix < (self.pml_thickness + 2)) {
            iz = 0;
            while (iz < (nz - 1)) {
                self.uz_r[(ix, iz)] = ((self.daz_r[(ix, iz)] * self.uz_r[(ix, iz)])
                    + (self.dbz_r[(ix, iz)] * (self.v_r[(ix, (iz + 1))] - self.v_r[(ix, iz)])));
                iz += 1;
            }
            ix += 1;
        }
        ix = 0;
        iz = (nz - 1);
        jx = (nx - 1);
        while (ix < (self.pml_thickness + 2)) {
            if (jx == (nx + self.pml_thickness)) {
                jx = 0;
            }
            self.uz_r[(ix, 0)] = self.uz_b[(jx, (self.pml_thickness + 1))];
            self.uz_r[(ix, iz)] = self.uz_b[(jx, 0)];
            ix += 1;
            jx += 1;
        }
        // Update top and bottom edge of main grid with new field values
        ix = 0;
        jz = (self.pml_thickness + 1);
        while (ix < nx) {
            if (abs[0 as usize] != 0) {
                pp[(ix, 0)] = (self.ux_b[(ix, jz)] + self.uz_b[(ix, jz)]);
            }
            if (abs[2 as usize] != 0) {
                pp[(ix, (nz - 1))] = (self.ux_b[(ix, 0)] + self.uz_b[(ix, 0)]);
            }
            ix += 1;
        }
        // Update left and right edges of main grid with new field values
        iz = 1;
        jx = (self.pml_thickness + 1);
        while (iz < (nz - 1)) {
            if (abs[1 as usize] != 0) {
                pp[(0, iz)] = (self.ux_r[(jx, iz)] + self.uz_r[(jx, iz)]);
            }
            if (abs[3 as usize] != 0) {
                pp[((nx - 1), iz)] = (self.ux_r[(0, iz)] + self.uz_r[(0, iz)]);
            }
            iz += 1;
        }
    }
}
