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

//! Piecewise cubic Hermite interpolation from par/main/ray_theory/rayt2dan.c.
//! Credits: CWP: Dave Hale (1992).

#![allow(unused_parens, unused_mut, unused_assignments)]

fn h00(x: f32) -> f32 {
    return ((((((2f64 * (x as f64)) * (x as f64)) * (x as f64))
        - ((3f64 * (x as f64)) * (x as f64)))
        + 1f64) as f32);
}

fn h01(x: f32) -> f32 {
    return ((((6f64 * (x as f64)) * (x as f64)) - (6f64 * (x as f64))) as f32);
}

fn h02(x: f32) -> f32 {
    return (((12f64 * (x as f64)) - 6f64) as f32);
}

fn h10(x: f32) -> f32 {
    return ((((((-2f64) * (x as f64)) * (x as f64)) * (x as f64))
        + ((3f64 * (x as f64)) * (x as f64))) as f32);
}

fn h11(x: f32) -> f32 {
    return (((((-6f64) * (x as f64)) * (x as f64)) + (6f64 * (x as f64))) as f32);
}

fn h12(x: f32) -> f32 {
    return ((((-12f64) * (x as f64)) + 6f64) as f32);
}

fn k00(x: f32) -> f32 {
    return ((((((x * x) * x) as f64) - ((2f64 * (x as f64)) * (x as f64))) + (x as f64)) as f32);
}

fn k01(x: f32) -> f32 {
    return (((((3f64 * (x as f64)) * (x as f64)) - (4f64 * (x as f64))) + 1f64) as f32);
}

fn k02(x: f32) -> f32 {
    return (((6f64 * (x as f64)) - 4f64) as f32);
}

fn k10(x: f32) -> f32 {
    return (((x * x) * x) - (x * x));
}

fn k11(x: f32) -> f32 {
    return ((((3f64 * (x as f64)) * (x as f64)) - (2f64 * (x as f64))) as f32);
}

fn k12(x: f32) -> f32 {
    return (((6f64 * (x as f64)) - 2f64) as f32);
}

fn build_tables() -> [[[f32; 4]; 101]; 3] {
    let mut tbl = [[[0.0; 4]; 101]; 3];
    {
        let mut i: i32 = 0;
        let mut x: f32 = 0.0;
        /* tabulate interpolator for 0th derivative */
        i = 0i32;
        while (i < 101i32) {
            x = ((((i as f32) as f64) / ((101i32 as f64) - 1f64)) as f32);
            tbl[(0i32) as usize][(i) as usize][(0i32) as usize] =
                (((-0.5f64) * (k00(x) as f64)) as f32);
            tbl[(0i32) as usize][(i) as usize][(1i32) as usize] =
                (((h00(x) as f64) - (0.5f64 * (k10(x) as f64))) as f32);
            tbl[(0i32) as usize][(i) as usize][(2i32) as usize] =
                (((h10(x) as f64) + (0.5f64 * (k00(x) as f64))) as f32);
            tbl[(0i32) as usize][(i) as usize][(3i32) as usize] =
                ((0.5f64 * (k10(x) as f64)) as f32);
            tbl[(1i32) as usize][(i) as usize][(0i32) as usize] =
                (((-0.5f64) * (k01(x) as f64)) as f32);
            tbl[(1i32) as usize][(i) as usize][(1i32) as usize] =
                (((h01(x) as f64) - (0.5f64 * (k11(x) as f64))) as f32);
            tbl[(1i32) as usize][(i) as usize][(2i32) as usize] =
                (((h11(x) as f64) + (0.5f64 * (k01(x) as f64))) as f32);
            tbl[(1i32) as usize][(i) as usize][(3i32) as usize] =
                ((0.5f64 * (k11(x) as f64)) as f32);
            tbl[(2i32) as usize][(i) as usize][(0i32) as usize] =
                (((-0.5f64) * (k02(x) as f64)) as f32);
            tbl[(2i32) as usize][(i) as usize][(1i32) as usize] =
                (((h02(x) as f64) - (0.5f64 * (k12(x) as f64))) as f32);
            tbl[(2i32) as usize][(i) as usize][(2i32) as usize] =
                (((h12(x) as f64) + (0.5f64 * (k02(x) as f64))) as f32);
            tbl[(2i32) as usize][(i) as usize][(3i32) as usize] =
                ((0.5f64 * (k12(x) as f64)) as f32);
            i += 1;
        }
        /* remember that tables have been built */
    }

    tbl
}

pub struct Vel2<'a> {
    nx: i32,
    nz: i32,
    nxm: i32,
    nzm: i32,
    xs: f32,
    xb: f32,
    zs: f32,
    zb: f32,
    sx: [f32; 3],
    sz: [f32; 3],
    vu: &'a [f32],
    tbl: &'static [[[f32; 4]; 101]; 3],
}
impl<'a> Vel2<'a> {
    pub fn new(nx: i32, dx: f32, fx: f32, nz: i32, dz: f32, fz: f32, v: &'a [f32]) -> Self {
        let xs = (1.0f64 / dx as f64) as f32;
        let xb = 4.0f32 - fx * xs;
        let zs = (1.0f64 / dz as f64) as f32;
        let zb = 4.0f32 - fz * zs;
        static TABLE: std::sync::OnceLock<[[[f32; 4]; 101]; 3]> = std::sync::OnceLock::new();
        Self {
            nx,
            nz,
            nxm: nx - 4,
            nzm: nz - 4,
            xs,
            xb,
            zs,
            zb,
            sx: [1.0, xs, xs * xs],
            sz: [1.0, zs, zs * zs],
            vu: v,
            tbl: TABLE.get_or_init(build_tables),
        }
    }
    pub fn interp(&self, x: f32, z: f32) -> (f32, f32, f32, f32, f32, f32) {
        let ax = self.xb + x * self.xs;
        let mut jx = ax as i32;
        let bx = ax - jx as f32;
        let lx = (if bx >= 0.0 {
            (bx * 100.0f32) as f64 + 0.5
        } else {
            (bx as f64 + 1.0) * 100.0 - 0.5
        }) as usize;
        let mx = -5 + jx;
        let az = self.zb + z * self.zs;
        let mut jz = az as i32;
        let bz = az - jz as f32;
        let lz = (if bz >= 0.0 {
            (bz * 100.0f32) as f64 + 0.5
        } else {
            (bz as f64 + 1.0) * 100.0 - 0.5
        }) as usize;
        let mz = -5 + jz;
        let kx = [0usize, 1, 0, 2, 1, 0];
        let kz = [0usize, 0, 1, 0, 1, 2];
        let mut vd = [0.0f32; 6];
        if mx >= 0 && mx <= self.nxm && mz >= 0 && mz <= self.nzm {
            let mut i = 0;
            while i < 6 {
                let px = &self.tbl[kx[i]][lx];
                let pz = &self.tbl[kz[i]][lz];
                vd[i] = self.sx[kx[i]]
                    * self.sz[kz[i]]
                    * (self.vu[((mx + 0) * self.nz + mz + 0) as usize] * px[0] * pz[0]
                        + self.vu[((mx + 0) * self.nz + mz + 1) as usize] * px[0] * pz[1]
                        + self.vu[((mx + 0) * self.nz + mz + 2) as usize] * px[0] * pz[2]
                        + self.vu[((mx + 0) * self.nz + mz + 3) as usize] * px[0] * pz[3]
                        + self.vu[((mx + 1) * self.nz + mz + 0) as usize] * px[1] * pz[0]
                        + self.vu[((mx + 1) * self.nz + mz + 1) as usize] * px[1] * pz[1]
                        + self.vu[((mx + 1) * self.nz + mz + 2) as usize] * px[1] * pz[2]
                        + self.vu[((mx + 1) * self.nz + mz + 3) as usize] * px[1] * pz[3]
                        + self.vu[((mx + 2) * self.nz + mz + 0) as usize] * px[2] * pz[0]
                        + self.vu[((mx + 2) * self.nz + mz + 1) as usize] * px[2] * pz[1]
                        + self.vu[((mx + 2) * self.nz + mz + 2) as usize] * px[2] * pz[2]
                        + self.vu[((mx + 2) * self.nz + mz + 3) as usize] * px[2] * pz[3]
                        + self.vu[((mx + 3) * self.nz + mz + 0) as usize] * px[3] * pz[0]
                        + self.vu[((mx + 3) * self.nz + mz + 1) as usize] * px[3] * pz[1]
                        + self.vu[((mx + 3) * self.nz + mz + 2) as usize] * px[3] * pz[2]
                        + self.vu[((mx + 3) * self.nz + mz + 3) as usize] * px[3] * pz[3]);
                i += 1;
            }
        } else {
            let mut i = 0;
            while i < 6 {
                let px = &self.tbl[kx[i]][lx];
                let pz = &self.tbl[kz[i]][lz];
                jx = 0;
                let mut jmx = mx;
                let mut sum = 0.0f32;
                while jx < 4 {
                    let mut mmx = jmx;
                    if mmx < 0 {
                        mmx = 0;
                    } else if mmx >= self.nx {
                        mmx = self.nx - 1;
                    }
                    jz = 0;
                    let mut jmz = mz;
                    while jz < 4 {
                        let mut mmz = jmz;
                        if mmz < 0 {
                            mmz = 0;
                        } else if mmz >= self.nz {
                            mmz = self.nz - 1;
                        }
                        sum += self.vu[(mmx * self.nz + mmz) as usize]
                            * px[jx as usize]
                            * pz[jz as usize];
                        jz += 1;
                        jmz += 1;
                    }
                    jx += 1;
                    jmx += 1;
                }
                vd[i] = self.sx[kx[i]] * self.sz[kz[i]] * sum;
                i += 1;
            }
        }
        (vd[0], vd[1], vd[2], vd[3], vd[4], vd[5])
    }
}
