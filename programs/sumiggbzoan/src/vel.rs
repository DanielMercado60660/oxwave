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

#![allow(unused_mut, unused_assignments, unused_parens, non_snake_case)]
fn h00(x: f32) -> f32 {
    return ((((((2_f64 * (x as f64)) * (x as f64)) * (x as f64))
        - ((3_f64 * (x as f64)) * (x as f64)))
        + 1_f64) as f32);
}
fn h01(x: f32) -> f32 {
    return ((((6_f64 * (x as f64)) * (x as f64)) - (6_f64 * (x as f64))) as f32);
}
fn h02(x: f32) -> f32 {
    return (((12_f64 * (x as f64)) - 6_f64) as f32);
}
fn h10(x: f32) -> f32 {
    return ((((((-2_f64) * (x as f64)) * (x as f64)) * (x as f64))
        + ((3_f64 * (x as f64)) * (x as f64))) as f32);
}
fn h11(x: f32) -> f32 {
    return (((((-6_f64) * (x as f64)) * (x as f64)) + (6_f64 * (x as f64))) as f32);
}
fn h12(x: f32) -> f32 {
    return ((((-12_f64) * (x as f64)) + 6_f64) as f32);
}
fn k00(x: f32) -> f32 {
    return ((((((x * x) * x) as f64) - ((2_f64 * (x as f64)) * (x as f64))) + (x as f64)) as f32);
}
fn k01(x: f32) -> f32 {
    return (((((3_f64 * (x as f64)) * (x as f64)) - (4_f64 * (x as f64))) + 1_f64) as f32);
}
fn k02(x: f32) -> f32 {
    return (((6_f64 * (x as f64)) - 4_f64) as f32);
}
fn k10(x: f32) -> f32 {
    return (((x * x) * x) - (x * x));
}
fn k11(x: f32) -> f32 {
    return ((((3_f64 * (x as f64)) * (x as f64)) - (2_f64 * (x as f64))) as f32);
}
fn k12(x: f32) -> f32 {
    return (((6_f64 * (x as f64)) - 2_f64) as f32);
}
fn buildTables() -> [[[f32; 4]; 101]; 3] {
    let mut tbl = [[[0.0; 4]; 101]; 3];
    let mut i: i32 = 0;
    let mut x: f32 = 0.0;
    /* tabulate interpolator for 0th derivative */
    {
        i = 0i32;
        while (i < 101i32) {
            x = ((((i as f32) as f64) / ((101i32 as f64) - 1_f64)) as f32);
            tbl[(0i32) as usize][(i) as usize][(0i32) as usize] =
                (((-0.5_f64) * (k00(x) as f64)) as f32);
            tbl[(0i32) as usize][(i) as usize][(1i32) as usize] =
                (((h00(x) as f64) - (0.5_f64 * (k10(x) as f64))) as f32);
            tbl[(0i32) as usize][(i) as usize][(2i32) as usize] =
                (((h10(x) as f64) + (0.5_f64 * (k00(x) as f64))) as f32);
            tbl[(0i32) as usize][(i) as usize][(3i32) as usize] =
                ((0.5_f64 * (k10(x) as f64)) as f32);
            tbl[(1i32) as usize][(i) as usize][(0i32) as usize] =
                (((-0.5_f64) * (k01(x) as f64)) as f32);
            tbl[(1i32) as usize][(i) as usize][(1i32) as usize] =
                (((h01(x) as f64) - (0.5_f64 * (k11(x) as f64))) as f32);
            tbl[(1i32) as usize][(i) as usize][(2i32) as usize] =
                (((h11(x) as f64) + (0.5_f64 * (k01(x) as f64))) as f32);
            tbl[(1i32) as usize][(i) as usize][(3i32) as usize] =
                ((0.5_f64 * (k11(x) as f64)) as f32);
            tbl[(2i32) as usize][(i) as usize][(0i32) as usize] =
                (((-0.5_f64) * (k02(x) as f64)) as f32);
            tbl[(2i32) as usize][(i) as usize][(1i32) as usize] =
                (((h02(x) as f64) - (0.5_f64 * (k12(x) as f64))) as f32);
            tbl[(2i32) as usize][(i) as usize][(2i32) as usize] =
                (((h12(x) as f64) + (0.5_f64 * (k02(x) as f64))) as f32);
            tbl[(2i32) as usize][(i) as usize][(3i32) as usize] =
                ((0.5_f64 * (k12(x) as f64)) as f32);
            i += 1;
        }
    }
    /* remember that tables have been built */

    tbl
}

/// Piecewise cubic Hermite interpolation and constant extrapolation, as in C.
pub(super) struct Vel2<'a> {
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
    vu: &'a [Vec<f32>],
}
impl<'a> Vel2<'a> {
    pub(super) fn new(
        nx: i32,
        dx: f32,
        fx: f32,
        nz: i32,
        dz: f32,
        fz: f32,
        vu: &'a [Vec<f32>],
    ) -> Self {
        let xs = (1.0 / dx as f64) as f32;
        let zs = (1.0 / dz as f64) as f32;
        Self {
            nx,
            nz,
            nxm: nx - 4,
            nzm: nz - 4,
            xs,
            xb: 4.0 - fx * xs,
            zs,
            zb: 4.0 - fz * zs,
            sx: [1.0, xs, xs * xs],
            sz: [1.0, zs, zs * zs],
            vu,
        }
    }
    pub(super) fn interp(&self, x: f32, z: f32) -> (f32, f32, f32, f32, f32, f32) {
        static TABLE: std::sync::OnceLock<[[[f32; 4]; 101]; 3]> = std::sync::OnceLock::new();
        let tbl = TABLE.get_or_init(buildTables);
        let ax = self.xb + x * self.xs;
        let jx = ax as i32;
        let bx = ax - jx as f32;
        let lx = if bx >= 0.0 {
            (bx * 100.0) as f64 + 0.5
        } else {
            (bx as f64 + 1.0) * 100.0 - 0.5
        } as usize;
        let mx = -5 + jx;
        let az = self.zb + z * self.zs;
        let jz = az as i32;
        let bz = az - jz as f32;
        let lz = if bz >= 0.0 {
            (bz * 100.0) as f64 + 0.5
        } else {
            (bz as f64 + 1.0) * 100.0 - 0.5
        } as usize;
        let mz = -5 + jz;
        let kx = [0, 1, 0, 2, 1, 0];
        let kz = [0, 0, 1, 0, 1, 2];
        let mut vd = [0.0f32; 6];
        for i in 0..6 {
            let px = &tbl[kx[i]][lx];
            let pz = &tbl[kz[i]][lz];
            let mut sum = 0.0;
            if mx >= 0 && mx <= self.nxm && mz >= 0 && mz <= self.nzm {
                // The C fast path is a left-associated sum, with no initial zero.
                sum = self.vu[mx as usize][mz as usize] * px[0] * pz[0];
                for j in 1..16 {
                    let a = j / 4;
                    let b = j % 4;
                    sum = sum + self.vu[mx as usize + a][mz as usize + b] * px[a] * pz[b];
                }
            } else {
                for jx in 0..4 {
                    let mmx = (mx + jx).clamp(0, self.nx - 1);
                    for jz in 0..4 {
                        let mmz = (mz + jz).clamp(0, self.nz - 1);
                        sum +=
                            self.vu[mmx as usize][mmz as usize] * px[jx as usize] * pz[jz as usize];
                    }
                }
            }
            vd[i] = self.sx[kx[i]] * self.sz[kz[i]] * sum;
        }
        (vd[0], vd[1], vd[2], vd[3], vd[4], vd[5])
    }
}
