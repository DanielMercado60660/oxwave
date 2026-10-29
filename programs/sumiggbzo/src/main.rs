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

//! SUMIGGBZO - Migration via Gaussian beams of zero-offset SU data.
//! Port of su/main/migration_inversion/sumiggbzo.c (Revision 1.15).
//! Credits: CWP: Dave Hale (algorithm), Jack K. Cohen, and John Stockwell
//! (reformatting for SU).

#![allow(
    non_snake_case,
    unused_parens,
    unused_mut,
    dead_code,
    clippy::too_many_arguments
)]
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_pfa_fft::*;
use oxwave_sinc_interp::ints8c;
use std::f64::consts::PI;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGGBZO - MIGration via Gaussian Beams of Zero-Offset SU data\t",
    "\t\t\t\t\t\t\t\t\t",
    " sumiggbzo <infile >outfile vfile=  nz= [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " vfile=                 name of file containing v(z,x)\t\t\t",
    " nz=                    number of depth samples\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt=from header\t\ttime sampling interval\t\t\t",
    " dx=from header(d2) or 1.0\tspatial sampling interval \t\t",
    " dz=1.0                 depth sampling interval\t\t\t",
    " fmin=0.025/dt          minimum frequency\t\t\t\t",
    " fmax=10*fmin           maximum frequency\t\t\t\t",
    " amin=-amax             minimum emergence angle; must be > -90 degrees\t",
    " amax=60                maximum emergence angle; must be < 90 degrees\t",
    " bwh=0.5*vavg/fmin      beam half-width; vavg denotes average velocity\t",
    " verbose=0\t\t =0 silent; =1 chatty\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: spatial units of v(z,x) must be the same as those of dx.\t",
    " v(z,x) is represented numerically in C-style binary floats v[x][z],\t",
    " where the depth direction is the fast direction in the data. Such\t",
    " models can be created with unif2 or makevel.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "(In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  \t\t",
    " denotes a v(z,x) array, the opposite of what Matlab and Fortran\t",
    " programmers may expect.)\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t\t",
    " In the event of a \"Segmentation Violation\" try reducing the value of",
    " the \"bwh\" parameter. Run program with verbose=1 do see the default\t",
    " value.\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let mut vfile = String::new();
    let mut nz = 0i32;
    let mut dz = 0.0f32;
    mustgetparstring("vfile", &mut vfile);
    mustgetparint("nz", &mut nz);
    mustgetparfloat("dz", &mut dz);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut dt = 0.0f32;
    let mut dx = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        if tr.dt != 0 {
            dt = (tr.dt as f64 / 1000000.0) as f32;
        } else {
            dt = 0.004;
            if verbose != 0 {
                warn!("tr.dt not set, assuming dt=0.004");
            }
        }
    }
    if getparfloat("dx", &mut dx) == 0 {
        if tr.d2 != 0.0 {
            dx = tr.d2;
        } else {
            dx = 1.0;
            if verbose != 0 {
                warn!("tr.d2 not set, assuming d2=1.0");
            }
        }
    }
    if getparfloat("dz", &mut dz) == 0 {
        dz = 1.0;
    }
    let mut fmin = 0.0f32;
    let mut fmax = 0.0f32;
    let mut amax = 0.0f32;
    let mut amin = 0.0f32;
    if getparfloat("fmin", &mut fmin) == 0 {
        fmin = (0.025 / dt as f64) as f32;
    }
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (10.0 * fmin as f64) as f32;
    }
    if getparfloat("amax", &mut amax) == 0 {
        amax = 60.0;
    }
    if getparfloat("amin", &mut amin) == 0 {
        amin = -amax;
    }
    let tracefp = etmpfile();
    let mut nx = 0usize;
    loop {
        nx += 1;
        efwrite(&tr.data[..nt], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    let nz = nz as usize;
    let mut fflat = vec![0.0f32; nt * nx];
    let mut vflat = vec![0.0f32; nz * nx];
    let mut g = vec![vec![0.0f32; nz]; nx];
    erewind(&tracefp);
    efread(&mut fflat, &tracefp);
    efclose(&tracefp);
    let vfp = efopen(&vfile, "r");
    if efread(&mut vflat, &vfp) != nz * nx {
        err!("error reading vfile=%s!\n", vfile);
    }
    let mut vavg = 0.0f32;
    for ix in 0..nx {
        for iz in 0..nz {
            vflat[ix * nz + iz] = (vflat[ix * nz + iz] as f64 * 0.5) as f32;
            vavg += vflat[ix * nz + iz];
        }
    }
    vavg /= (nx * nz) as f32;
    let mut bwh = 0.0f32;
    if getparfloat("bwh", &mut bwh) == 0 {
        bwh = vavg / fmin;
    }
    checkpars();
    if verbose != 0 {
        warn!("bhw=%f", bwh);
    }
    let f: Vec<Vec<f32>> = fflat.chunks(nt).map(|x| x.to_vec()).collect();
    let v: Vec<Vec<f32>> = vflat.chunks(nz).map(|x| x.to_vec()).collect();
    miggbzo(
        bwh, fmin, fmax, amin, amax, nt, dt, nx, dx, nz, dz, &f, &v, &mut g, verbose,
    );
    tr.ns = nz as u16;
    tr.trid = 130; // TRID_DEPTH
    tr.d1 = dz;
    tr.d2 = dx;
    for ix in 0..nx {
        tr.tracl = (ix + 1) as i32;
        tr.tracr = (ix + 1) as i32;
        tr.data[..nz].copy_from_slice(&g[ix]);
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

fn miggbzo(
    bwh: f32,
    fmin: f32,
    fmax: f32,
    mut amin: f32,
    mut amax: f32,
    nt: usize,
    dt: f32,
    nx: usize,
    dx: f32,
    nz: usize,
    dz: f32,
    f: &[Vec<f32>],
    v: &[Vec<f32>],
    g: &mut [Vec<f32>],
    verbose: i32,
) {
    let (ft, fx, fz) = (0.0f32, 0.0f32, 0.0f32);
    amin = (amin as f64 * (PI / 180.0)) as f32;
    amax = (amax as f64 * (PI / 180.0)) as f32;
    if amin > amax {
        std::mem::swap(&mut amin, &mut amax);
    }
    let xwh = (3.0 * bwh as f64) as f32;
    let dxb =
        (nint((bwh as f64 * (2.0 * fmin as f64 / fmax as f64).sqrt()) / dx as f64) as f32) * dx;
    let nxb = (1.0f32 + (nx as i32 - 1) as f32 * dx / dxb) as i32;
    let fxb =
        (fx as f64 + 0.5 * ((nx as i32 - 1) as f32 * dx - (nxb - 1) as f32 * dxb) as f64) as f32;
    if verbose != 0 {
        warn!("nxb=%d dxb=%f fxb=%f ", nxb, dxb, fxb);
    }
    let mut vmin = v[0][0];
    for ix in 1..nx {
        if v[ix][0] < vmin {
            vmin = v[ix][0];
        }
    }
    if verbose != 0 {
        warn!("vmin=%f fmin=%f fmax=%f", vmin, fmin, fmax);
    }
    let pxmin = ((amin as f64).sin() / vmin as f64) as f32;
    let pxmax = ((amax as f64).sin() / vmin as f64) as f32;
    let dpx = (1.0 / (2.0 * xwh as f64 * ((fmin * fmax) as f64).sqrt())) as f32;
    let npx = (1.0f32 + (pxmax - pxmin) / dpx) as i32;
    let fpx = (pxmin as f64 + 0.5 * (pxmax - pxmin - (npx - 1) as f32 * dpx) as f64) as f32;
    let mut taupad = (if pxmin.abs() > pxmax.abs() {
        pxmin.abs()
    } else {
        pxmax.abs()
    }) * xwh;
    taupad = nint((taupad / dt) as f64) as f32 * dt;
    let ftau = ft - taupad;
    let dtau = dt;
    let ntau = (nt as f64 + 2.0 * taupad as f64 / dtau as f64) as usize;
    if verbose != 0 {
        warn!("npx=%d dpx=%g fpx=%g", npx, dpx, fpx);
        warn!("pxmin=%f pxmax=%f", pxmin, pxmax);
        warn!("amin=%f amax=%f", amin, amax);
    }
    for row in g.iter_mut() {
        row.fill(0.0);
    }
    let mut xb = fxb;
    for ixb in 0..nxb {
        let ix = nint(((xb - fx) / dx) as f64);
        let ixlo = (ix - nint((xwh / dx) as f64)).max(0);
        let ixhi = (ix + nint((xwh / dx) as f64)).min(nx as i32 - 1);
        let nxw = 1 + ixhi - ixlo;
        let fxw = fx + (ixlo - ix) as f32 * dx;
        if verbose != 0 {
            warn!("ixb/nxb = %d/%d  ix = %d", ixb, nxb, ix);
        }
        let mut b = vec![vec![0.0f32; ntau]; npx as usize];
        formBeams(
            bwh,
            dxb,
            fmin,
            nt,
            dt,
            ft,
            nxw as usize,
            dx,
            fxw,
            &f[ixlo as usize..],
            ntau,
            dtau,
            ftau,
            npx as usize,
            dpx,
            fpx,
            &mut b,
        );
        let mut px = fpx;
        for ipx in 0..npx as usize {
            if !(((px * v[ix as usize][0]) as f64) > (amax as f64).sin() + 0.01)
                && !(((px * v[ix as usize][0]) as f64) < (amin as f64).sin() - 0.01)
            {
                let a0 = -((px * v[ix as usize][0]) as f64).asin() as f32;
                let x0 = fx + ix as f32 * dx;
                let bwhc = (bwh as f64 * (a0 as f64).cos()) as f32;
                let mut ray = makeRay(x0, fz, a0, nt, dt, ft, nx, dx, fx, nz, dz, fz, v);
                accBeam(
                    &mut ray, fmin, bwhc, ntau, dtau, ftau, &b[ipx], nx, dx, fx, nz, dz, fz, g,
                );
            }
            px += dpx;
        }
        xb += dxb;
    }
}

#[derive(Clone, Copy, Default)]
struct RayStep {
    t: f32,
    a: f32,
    x: f32,
    z: f32,
    q1: f32,
    p1: f32,
    q2: f32,
    p2: f32,
    kmah: i32,
    c: f32,
    s: f32,
    v: f32,
    dvdx: f32,
    dvdz: f32,
}
struct Ray {
    rs: Vec<RayStep>,
    ic: usize,
    circles: Vec<Circle>,
}
#[derive(Clone, Copy)]
struct Circle {
    irsf: usize,
    irsl: usize,
    x: f32,
    z: f32,
    r: f32,
}

struct Vel2<'a> {
    nx: usize,
    nz: usize,
    xs: f32,
    xb: f32,
    zs: f32,
    zb: f32,
    sx: [f32; 3],
    sz: [f32; 3],
    vu: &'a [Vec<f32>],
    tbl: [[[f32; 4]; 101]; 3],
}
impl<'a> Vel2<'a> {
    fn new(nx: usize, dx: f32, fx: f32, nz: usize, dz: f32, fz: f32, v: &'a [Vec<f32>]) -> Self {
        let xs = (1.0 / dx as f64) as f32;
        let zs = (1.0 / dz as f64) as f32;
        let mut tbl = [[[0.0f32; 4]; 101]; 3];
        buildTables(&mut tbl);
        Self {
            nx,
            nz,
            xs,
            xb: 4.0 - fx * xs,
            zs,
            zb: 4.0 - fz * zs,
            sx: [1.0, xs, xs * xs],
            sz: [1.0, zs, zs * zs],
            vu: v,
            tbl,
        }
    }
    fn interp(&self, x: f32, z: f32) -> [f32; 6] {
        let ax = self.xb + x * self.xs;
        let jx = ax as i32;
        let bx = ax - jx as f32;
        let lx = (if bx >= 0.0 {
            ((bx * 100.0) as f64 + 0.5) as i32
        } else {
            ((bx as f64 + 1.0) * 100.0 - 0.5) as i32
        }) as usize;
        let mx = -5 + jx;
        let az = self.zb + z * self.zs;
        let jz = az as i32;
        let bz = az - jz as f32;
        let lz = (if bz >= 0.0 {
            ((bz * 100.0) as f64 + 0.5) as i32
        } else {
            ((bz as f64 + 1.0) * 100.0 - 0.5) as i32
        }) as usize;
        let mz = -5 + jz;
        let kx = [0, 1, 0, 2, 1, 0];
        let kz = [0, 0, 1, 0, 1, 2];
        let mut vd = [0.0f32; 6];
        for i in 0..6 {
            let px = &self.tbl[kx[i]][lx];
            let pz = &self.tbl[kz[i]][lz];
            let interior =
                mx >= 0 && mx <= self.nx as i32 - 4 && mz >= 0 && mz <= self.nz as i32 - 4;
            if interior {
                let mx = mx as usize;
                let mz = mz as usize;
                let vu = self.vu;
                vd[i] = self.sx[kx[i]]
                    * self.sz[kz[i]]
                    * (vu[mx][mz] * px[0] * pz[0]
                        + vu[mx][mz + 1] * px[0] * pz[1]
                        + vu[mx][mz + 2] * px[0] * pz[2]
                        + vu[mx][mz + 3] * px[0] * pz[3]
                        + vu[mx + 1][mz] * px[1] * pz[0]
                        + vu[mx + 1][mz + 1] * px[1] * pz[1]
                        + vu[mx + 1][mz + 2] * px[1] * pz[2]
                        + vu[mx + 1][mz + 3] * px[1] * pz[3]
                        + vu[mx + 2][mz] * px[2] * pz[0]
                        + vu[mx + 2][mz + 1] * px[2] * pz[1]
                        + vu[mx + 2][mz + 2] * px[2] * pz[2]
                        + vu[mx + 2][mz + 3] * px[2] * pz[3]
                        + vu[mx + 3][mz] * px[3] * pz[0]
                        + vu[mx + 3][mz + 1] * px[3] * pz[1]
                        + vu[mx + 3][mz + 2] * px[3] * pz[2]
                        + vu[mx + 3][mz + 3] * px[3] * pz[3]);
            } else {
                // Constant extrapolation at the edges of the velocity model.
                let mut sum = 0.0f32;
                for jx in 0..4 {
                    let mmx = (mx + jx as i32).clamp(0, self.nx as i32 - 1) as usize;
                    for jz in 0..4 {
                        let mmz = (mz + jz as i32).clamp(0, self.nz as i32 - 1) as usize;
                        sum += self.vu[mmx][mmz] * px[jx] * pz[jz];
                    }
                }
                vd[i] = self.sx[kx[i]] * self.sz[kz[i]] * sum;
            }
        }
        vd
    }
}
fn h00(x: f32) -> f32 {
    ((((((2.0f64 * (x as f64)) * (x as f64)) * (x as f64)) - ((3.0f64 * (x as f64)) * (x as f64)))
        + 1.0f64) as f32)
}
fn h01(x: f32) -> f32 {
    ((((6.0f64 * (x as f64)) * (x as f64)) - (6.0f64 * (x as f64))) as f32)
}
fn h02(x: f32) -> f32 {
    (((12.0f64 * (x as f64)) - 6.0f64) as f32)
}
fn h10(x: f32) -> f32 {
    ((((((-2.0f64) * (x as f64)) * (x as f64)) * (x as f64)) + ((3.0f64 * (x as f64)) * (x as f64)))
        as f32)
}
fn h11(x: f32) -> f32 {
    (((((-6.0f64) * (x as f64)) * (x as f64)) + (6.0f64 * (x as f64))) as f32)
}
fn h12(x: f32) -> f32 {
    ((((-12.0f64) * (x as f64)) + 6.0f64) as f32)
}
fn k00(x: f32) -> f32 {
    ((((((x * x) * x) as f64) - ((2.0f64 * (x as f64)) * (x as f64))) + (x as f64)) as f32)
}
fn k01(x: f32) -> f32 {
    (((((3.0f64 * (x as f64)) * (x as f64)) - (4.0f64 * (x as f64))) + 1.0f64) as f32)
}
fn k02(x: f32) -> f32 {
    (((6.0f64 * (x as f64)) - 4.0f64) as f32)
}
fn k10(x: f32) -> f32 {
    (((x * x) * x) - (x * x))
}
fn k11(x: f32) -> f32 {
    ((((3.0f64 * (x as f64)) * (x as f64)) - (2.0f64 * (x as f64))) as f32)
}
fn k12(x: f32) -> f32 {
    (((6.0f64 * (x as f64)) - 2.0f64) as f32)
}
fn buildTables(tbl: &mut [[[f32; 4]; 101]; 3]) {
    for i in 0..101 {
        let x = (i as f64 / 100.0) as f32;
        tbl[0][i][0] = (((-0.5f64) * (k00(x) as f64)) as f32);
        tbl[0][i][1] = (((h00(x) as f64) - (0.5f64 * (k10(x) as f64))) as f32);
        tbl[0][i][2] = (((h10(x) as f64) + (0.5f64 * (k00(x) as f64))) as f32);
        tbl[0][i][3] = ((0.5f64 * (k10(x) as f64)) as f32);
        tbl[1][i][0] = (((-0.5f64) * (k01(x) as f64)) as f32);
        tbl[1][i][1] = (((h01(x) as f64) - (0.5f64 * (k11(x) as f64))) as f32);
        tbl[1][i][2] = (((h11(x) as f64) + (0.5f64 * (k01(x) as f64))) as f32);
        tbl[1][i][3] = ((0.5f64 * (k11(x) as f64)) as f32);
        tbl[2][i][0] = (((-0.5f64) * (k02(x) as f64)) as f32);
        tbl[2][i][1] = (((h02(x) as f64) - (0.5f64 * (k12(x) as f64))) as f32);
        tbl[2][i][2] = (((h12(x) as f64) + (0.5f64 * (k02(x) as f64))) as f32);
        tbl[2][i][3] = ((0.5f64 * (k12(x) as f64)) as f32);
    }
}

fn makeRay(
    x0: f32,
    z0: f32,
    a0: f32,
    nt: usize,
    dt: f32,
    ft: f32,
    nx: usize,
    dx: f32,
    fx: f32,
    nz: usize,
    dz: f32,
    fz: f32,
    vxz: &[Vec<f32>],
) -> Ray {
    let vel2 = Vel2::new(nx, dx, fx, nz, dz, fz, vxz);
    let lx = fx + (nx as i32 - 1) as f32 * dx;
    let lz = fz + (nz as i32 - 1) as f32 * dz;
    let mut c = (a0 as f64).cos() as f32;
    let mut s = (a0 as f64).sin() as f32;
    let mut cc = c * c;
    let mut ss = s * s;
    let [mut v, mut dvdx, mut dvdz, mut ddvdxdx, mut ddvdxdz, mut ddvdzdz] = vel2.interp(x0, z0);
    let mut ddvdndn = (((((ddvdxdx * cc) as f64)
        - (((2.0f64 * (ddvdxdz as f64)) * (s as f64)) * (c as f64)))
        + ((ddvdzdz * ss) as f64)) as f32);

    let mut vv = v * v;
    let mut t = ft;
    let mut a = a0;
    let mut x = x0;
    let mut z = z0;
    let mut q1 = 1.0f32;
    let mut p1 = 0.0f32;
    let mut q2 = 0.0f32;
    let mut p2 = 1.0f32;
    let mut kmah = 0i32;
    let mut rs = Vec::with_capacity(nt);
    rs.push(RayStep {
        t,
        a,
        x,
        z,
        q1,
        p1,
        q2,
        p2,
        kmah,
        c,
        s,
        v,
        dvdx,
        dvdz,
    });
    for _it in 1..nt {
        let h = dt;
        let hhalf = (dt as f64 / 2.0) as f32;
        let hsixth = (dt as f64 / 6.0) as f32;
        if x < fx || x > lx || z < fz || z > lz {
            break;
        }
        let q2old = q2;
        let mut xt: f32;
        let mut zt: f32;
        let mut at: f32;
        let mut p1t: f32;
        let mut q1t: f32;
        let mut p2t: f32;
        let mut q2t: f32;
        let mut dx: f32;
        let mut dz: f32;
        let mut da: f32;
        let mut dp1: f32;
        let mut dq1: f32;
        let mut dp2: f32;
        let mut dq2: f32;
        let mut dxt: f32;
        let mut dzt: f32;
        let mut dat: f32;
        let mut dp1t: f32;
        let mut dq1t: f32;
        let mut dp2t: f32;
        let mut dq2t: f32;
        let mut dxm: f32;
        let mut dzm: f32;
        let mut dam: f32;
        let mut dp1m: f32;
        let mut dq1m: f32;
        let mut dp2m: f32;
        let mut dq2m: f32;
        dx = (v * s);
        dz = (v * c);
        da = ((dvdz * s) - (dvdx * c));
        dp1 = (((-ddvdndn) * q1) / v);
        dq1 = (vv * p1);
        dp2 = (((-ddvdndn) * q2) / v);
        dq2 = (vv * p2);
        xt = (x + (hhalf * dx));
        zt = (z + (hhalf * dz));
        at = (a + (hhalf * da));
        p1t = (p1 + (hhalf * dp1));
        q1t = (q1 + (hhalf * dq1));
        p2t = (p2 + (hhalf * dp2));
        q2t = (q2 + (hhalf * dq2));
        c = ((at as f64).cos() as f32);
        s = ((at as f64).sin() as f32);
        cc = (c * c);
        ss = (s * s);
        [v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz] = vel2.interp(xt, zt);
        ddvdndn = (((((ddvdxdx * cc) as f64)
            - (((2.0f64 * (ddvdxdz as f64)) * (s as f64)) * (c as f64)))
            + ((ddvdzdz * ss) as f64)) as f32);
        vv = (v * v);
        dxt = (v * s);
        dzt = (v * c);
        dat = ((dvdz * s) - (dvdx * c));
        dp1t = (((-ddvdndn) * q1t) / v);
        dq1t = (vv * p1t);
        dp2t = (((-ddvdndn) * q2t) / v);
        dq2t = (vv * p2t);
        xt = (x + (hhalf * dxt));
        zt = (z + (hhalf * dzt));
        at = (a + (hhalf * dat));
        p1t = (p1 + (hhalf * dp1t));
        q1t = (q1 + (hhalf * dq1t));
        p2t = (p2 + (hhalf * dp2t));
        q2t = (q2 + (hhalf * dq2t));
        c = ((at as f64).cos() as f32);
        s = ((at as f64).sin() as f32);
        cc = (c * c);
        ss = (s * s);
        [v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz] = vel2.interp(xt, zt);
        ddvdndn = (((((ddvdxdx * cc) as f64)
            - (((2.0f64 * (ddvdxdz as f64)) * (s as f64)) * (c as f64)))
            + ((ddvdzdz * ss) as f64)) as f32);
        vv = (v * v);
        dxm = (v * s);
        dzm = (v * c);
        dam = ((dvdz * s) - (dvdx * c));
        dp1m = (((-ddvdndn) * q1t) / v);
        dq1m = (vv * p1t);
        dp2m = (((-ddvdndn) * q2t) / v);
        dq2m = (vv * p2t);
        xt = (x + (h * dxm));
        zt = (z + (h * dzm));
        at = (a + (h * dam));
        p1t = (p1 + (h * dp1m));
        q1t = (q1 + (h * dq1m));
        p2t = (p2 + (h * dp2m));
        q2t = (q2 + (h * dq2m));
        dxm = (dxm + dxt);
        dzm = (dzm + dzt);
        dam = (dam + dat);
        dp1m = (dp1m + dp1t);
        dq1m = (dq1m + dq1t);
        dp2m = (dp2m + dp2t);
        dq2m = (dq2m + dq2t);
        c = ((at as f64).cos() as f32);
        s = ((at as f64).sin() as f32);
        cc = (c * c);
        ss = (s * s);
        [v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz] = vel2.interp(xt, zt);
        ddvdndn = (((((ddvdxdx * cc) as f64)
            - (((2.0f64 * (ddvdxdz as f64)) * (s as f64)) * (c as f64)))
            + ((ddvdzdz * ss) as f64)) as f32);
        vv = (v * v);
        dxt = (v * s);
        dzt = (v * c);
        dat = ((dvdz * s) - (dvdx * c));
        dp1t = (((-ddvdndn) * q1t) / v);
        dq1t = (vv * p1t);
        dp2t = (((-ddvdndn) * q2t) / v);
        dq2t = (vv * p2t);
        x = (((x as f64) + ((hsixth as f64) * (((dx + dxt) as f64) + (2.0f64 * (dxm as f64)))))
            as f32);
        z = (((z as f64) + ((hsixth as f64) * (((dz + dzt) as f64) + (2.0f64 * (dzm as f64)))))
            as f32);
        a = (((a as f64) + ((hsixth as f64) * (((da + dat) as f64) + (2.0f64 * (dam as f64)))))
            as f32);
        p1 = (((p1 as f64) + ((hsixth as f64) * (((dp1 + dp1t) as f64) + (2.0f64 * (dp1m as f64)))))
            as f32);
        q1 = (((q1 as f64) + ((hsixth as f64) * (((dq1 + dq1t) as f64) + (2.0f64 * (dq1m as f64)))))
            as f32);
        p2 = (((p2 as f64) + ((hsixth as f64) * (((dp2 + dp2t) as f64) + (2.0f64 * (dp2m as f64)))))
            as f32);
        q2 = (((q2 as f64) + ((hsixth as f64) * (((dq2 + dq2t) as f64) + (2.0f64 * (dq2m as f64)))))
            as f32);
        c = ((a as f64).cos() as f32);
        s = ((a as f64).sin() as f32);
        cc = (c * c);
        ss = (s * s);
        [v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz] = vel2.interp(x, z);
        ddvdndn = (((((ddvdxdx * cc) as f64)
            - (((2.0f64 * (ddvdxdz as f64)) * (s as f64)) * (c as f64)))
            + ((ddvdzdz * ss) as f64)) as f32);
        vv = (v * v);

        if (q2 <= 0.0 && q2old > 0.0) || (q2 >= 0.0 && q2old < 0.0) {
            kmah += 1;
        }
        t += dt;
        rs.push(RayStep {
            t,
            a,
            x,
            z,
            q1,
            p1,
            q2,
            p2,
            kmah,
            c,
            s,
            v,
            dvdx,
            dvdz,
        });
    }
    Ray {
        rs,
        ic: 0,
        circles: Vec::new(),
    }
}
fn makeCircles(nc: usize, rs: &[RayStep]) -> Vec<Circle> {
    let nrs = rs.len();
    let nrsc = 1 + (nrs - 1) / nc;
    let mut circles = Vec::with_capacity(nc);
    for ic in 0..nc {
        let irsf = (ic * nrsc).min(nrs - 1);
        let irsl = (ic * nrsc + nrsc - 1).min(nrs - 1);
        let mut xmin = rs[irsf].x;
        let mut xmax = xmin;
        let mut zmin = rs[irsf].z;
        let mut zmax = zmin;
        for irs in irsf + 1..=irsl {
            if rs[irs].x < xmin {
                xmin = rs[irs].x;
            }
            if rs[irs].x > xmax {
                xmax = rs[irs].x;
            }
            if rs[irs].z < zmin {
                zmin = rs[irs].z;
            }
            if rs[irs].z > zmax {
                zmax = rs[irs].z;
            }
        }
        let x = (0.5 * (xmin + xmax) as f64) as f32;
        let z = (0.5 * (zmin + zmax) as f64) as f32;
        let r = (((x - xmin) * (x - xmin) + (z - zmin) * (z - zmin)) as f64).sqrt() as f32;
        circles.push(Circle {
            irsf,
            irsl,
            x,
            z,
            r,
        });
    }
    circles
}
fn nearestRayStep(ray: &mut Ray, x: f32, z: f32) -> usize {
    if ray.circles.is_empty() {
        let nc = (ray.rs.len() as f32 as f64).sqrt() as usize;
        ray.circles = makeCircles(nc, &ray.rs);
        ray.ic = 0;
    }
    let nc = ray.circles.len() as i32;
    let c = &ray.circles;
    let mut ic = ray.ic;
    let mut dx = x - c[ic].x;
    let mut dz = z - c[ic].z;
    let mut dmin = (2.0 * (((dx * dx + dz * dz) as f64).sqrt() + c[ic].r as f64)) as f32;
    let mut dsmin = dmin * dmin;
    let mut irsmin = 0;
    let mut jc = ic as i32;
    let mut js = 0i32;
    for _kc in 0..nc {
        dx = x - c[jc as usize].x;
        dz = z - c[jc as usize].z;
        let ds = dx * dx + dz * dz;
        let rdmin = c[jc as usize].r + dmin;
        if ds <= rdmin * rdmin {
            let mut update = false;
            for irs in c[jc as usize].irsf..=c[jc as usize].irsl {
                dx = x - ray.rs[irs].x;
                dz = z - ray.rs[irs].z;
                let ds = dx * dx + dz * dz;
                if ds < dsmin {
                    dsmin = ds;
                    irsmin = irs;
                    update = true;
                }
            }
            if update {
                dmin = (dsmin as f64).sqrt() as f32;
                ic = jc as usize;
            }
        }
        js = if js > 0 { -js - 1 } else { -js + 1 };
        jc += js;
        if jc < 0 || jc >= nc {
            js = if js > 0 { -js - 1 } else { -js + 1 };
            jc += js;
        }
    }
    ray.ic = ic;
    irsmin
}

fn formBeams(
    bwh: f32,
    dxb: f32,
    fmin: f32,
    nt: usize,
    dt: f32,
    ft: f32,
    nx: usize,
    dx: f32,
    fx: f32,
    f: &[Vec<f32>],
    ntau: usize,
    _dtau: f32,
    ftau: f32,
    npx: usize,
    dpx: f32,
    fpx: f32,
    g: &mut [Vec<f32>],
) {
    let wmin = (2.0 * PI * fmin as f64) as f32;
    let pxmax = if dpx < 0.0 {
        fpx
    } else {
        fpx + (npx as i32 - 1) as f32 * dpx
    };
    let xmax = if dx < 0.0 {
        fx
    } else {
        fx + (nx as i32 - 1) as f32 * dx
    };
    let ntpad = ((pxmax * xmax).abs() / dt) as i32;
    let ntfft = npfar((nt as i32 + ntpad).max(ntau as i32)) as usize;
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fftscl = (1.0 / ntfft as f64) as f32;
    let mut fpad = vec![0.0f32; ntfft];
    let mut gpad = vec![0.0f32; ntfft];
    let mut cf = vec![vec![Complex::default(); nw]; nx];
    let mut cg = vec![vec![Complex::default(); nw]; npx];
    let mut cfx = vec![Complex::default(); nx];
    let mut cgpx = vec![Complex::default(); npx];
    for ix in 0..nx {
        fpad[..nt].copy_from_slice(&f[ix][..nt]);
        fpad[nt..].fill(0.0);
        pfarc(1, &fpad, &mut cf[ix]);
    }
    let mut w = 0.0f32;
    for iw in 0..nw {
        let scale = (-0.5 * w as f64 / (wmin * bwh * bwh) as f64) as f32;
        let amp = ((fftscl * dpx * w) as f64 / (2.0 * PI)
            * (w as f64 / (PI * wmin as f64)).sqrt()
            * dxb as f64
            / bwh as f64) as f32;
        let phase = w * (ft - ftau);
        let a = (amp as f64 * (phase as f64).cos()) as f32;
        let b = (amp as f64 * (phase as f64).sin()) as f32;
        let mut x = fx;
        for ix in 0..nx {
            let es = ((scale * x * x) as f64).exp() as f32;
            let as_ = a * es;
            let bs = b * es;
            let cfr = cf[ix][iw].r;
            let cfi = cf[ix][iw].i;
            cfx[ix].r = as_ * cfr - bs * cfi;
            cfx[ix].i = bs * cfr + as_ * cfi;
            x += dx;
        }
        xtop(w, nx, dx, fx, &cfx, npx, dpx, fpx, &mut cgpx);
        for ipx in 0..npx {
            cg[ipx][iw] = cgpx[ipx];
        }
        w += dw;
    }
    for ipx in 0..npx {
        pfacr(-1, &cg[ipx], &mut gpad);
        g[ipx][..ntau].copy_from_slice(&gpad[..ntau]);
    }
}

fn xtop(
    w: f32,
    nx: usize,
    dx: f32,
    fx: f32,
    g: &[Complex],
    np: usize,
    dp: f32,
    fp: f32,
    h: &mut [Complex],
) {
    let lwrap = 8usize;
    let nxfft = npfa(((nx + lwrap) * 2) as i32) as usize;
    let nk = nxfft;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let fk = (-PI / dx as f64) as f32;
    let ek = (PI / dx as f64) as f32;
    let fka = fk - lwrap as f32 * dk;
    let nka = lwrap + nk + lwrap;
    let mut gka = vec![Complex::default(); nka];
    let mut hp = vec![Complex::default(); np];
    let mut kp = vec![0.0f32; np];
    for ix in 0..nx {
        gka[lwrap + ix].r = dx * g[ix].r;
        gka[lwrap + ix].i = dx * g[ix].i;
    }
    for ix in (1..nx).step_by(2) {
        gka[lwrap + ix].r = -gka[lwrap + ix].r;
        gka[lwrap + ix].i = -gka[lwrap + ix].i;
    }
    pfacc(-1, &mut gka[lwrap..lwrap + nk]);
    for ik in 0..lwrap {
        gka[ik] = gka[ik + nk];
    }
    for ik in lwrap + nk..lwrap + nk + lwrap {
        gka[ik] = gka[ik - nk];
    }
    let xshift = (0.5 * (nx as i32 - 1) as f64 * dx as f64) as f32;
    let mut k = fka;
    for ik in 0..nka {
        let phase = k * xshift;
        let c = (phase as f64).cos() as f32;
        let s = (phase as f64).sin() as f32;
        let temp = gka[ik].r * c - gka[ik].i * s;
        gka[ik].i = gka[ik].r * s + gka[ik].i * c;
        gka[ik].r = temp;
        k += dk;
    }
    let mut p = fp;
    for ip in 0..np {
        kp[ip] = w * p;
        if kp[ip] < fk && kp[ip] < ek {
            kp[ip] = (fk as f64 - 1000.0 * dk as f64) as f32;
        } else if kp[ip] > fk && kp[ip] > ek {
            kp[ip] = (ek as f64 + 1000.0 * dk as f64) as f32;
        }
        p += dp;
    }
    ints8c(
        dk,
        fka,
        &gka,
        Complex::default(),
        Complex::default(),
        &kp,
        &mut hp,
    );
    let xshift = (-fx as f64 - 0.5 * (nx as i32 - 1) as f64 * dx as f64) as f32;
    for ip in 0..np {
        let phase = kp[ip] * xshift;
        let c = (phase as f64).cos() as f32;
        let s = (phase as f64).sin() as f32;
        h[ip].r = hp[ip].r * c - hp[ip].i * s;
        h[ip].i = hp[ip].r * s + hp[ip].i * c;
    }
}

struct BeamData {
    ntr: usize,
    dtr: f32,
    ftr: f32,
    nti: usize,
    dti: f32,
    fti: f32,
    cf: Vec<Vec<Complex>>,
}
fn beamData(wmin: f32, nt: usize, dt: f32, ft: f32, f: &[f32]) -> BeamData {
    let ntpad = 25usize;
    let ntfft = npfaro((nt + ntpad) as i32, (2 * (nt + ntpad)) as i32) as usize;
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let iwnyq = nw - 1;
    let ntrfft = npfao((4 * ntfft) as i32, (5 * ntfft) as i32) as usize;
    let dtr = dt * ntfft as f32 / ntrfft as f32;
    let ftr = ft;
    let ntr = (1.0f32 + (nt + ntpad - 1) as f32 * dt / dtr) as usize;
    let nti = 6usize;
    let dti = (-5.0 / (wmin * (nti - 1) as f32) as f64) as f32;
    let fti = 0.0f32;
    let mut cf = vec![vec![Complex::default(); ntr]; nti];
    let mut fa = vec![0.0f32; ntfft];
    let mut ca = vec![Complex::default(); nw];
    let mut cb = vec![Complex::default(); ntrfft];
    fa[..nt].copy_from_slice(&f[..nt]);
    pfarc(1, &fa, &mut ca);
    for iw in 1..iwnyq {
        ca[iw].r = (ca[iw].r as f64 * 2.0) as f32;
        ca[iw].i = (ca[iw].i as f64 * 2.0) as f32;
    }
    let mut ti = fti;
    for iti in 0..nti {
        let mut w = 0.0f32;
        for iw in 0..nw {
            let scale = ((w * ti) as f64).exp() as f32;
            cb[iw].r = ca[iw].r * scale;
            cb[iw].i = ca[iw].i * scale;
            w += dw;
        }
        cb[nw..].fill(Complex::default());
        pfacc(-1, &mut cb);
        let scale = (1.0 / ntfft as f64) as f32;
        for itr in 0..ntr {
            cf[iti][itr].r = scale * cb[itr].r;
            cf[iti][itr].i = scale * cb[itr].i;
        }
        ti += dti;
    }
    BeamData {
        ntr,
        dtr,
        ftr,
        nti,
        dti,
        fti,
        cf,
    }
}

#[derive(Clone, Copy, Default)]
struct Cell {
    live: i32,
    dead: i32,
    tr: f32,
    ti: f32,
    ar: f32,
    ai: f32,
}
struct Cells<'a> {
    lx: i32,
    mx: i32,
    nx: usize,
    dx: f32,
    fx: f32,
    lz: i32,
    mz: i32,
    nz: usize,
    dz: f32,
    fz: f32,
    live: i32,
    dead: i32,
    wmin: f32,
    lmin: f32,
    cell: Vec<Vec<Cell>>,
    ray: &'a mut Ray,
    bd: BeamData,
    g: &'a mut [Vec<f32>],
}
fn accBeam(
    ray: &mut Ray,
    fmin: f32,
    lmin: f32,
    nt: usize,
    dt: f32,
    ft: f32,
    f: &[f32],
    nx: usize,
    dx: f32,
    fx: f32,
    nz: usize,
    dz: f32,
    fz: f32,
    g: &mut [Vec<f32>],
) {
    let wmin = (2.0 * PI * fmin as f64) as f32;
    // franuni is not provided by the frozen libraries. The C random values
    // are solely visitation tags for newly allocated (uninitialized) cells.
    // Fresh zeroed cells with nonzero tags preserve that traversal without
    // reimplementing the unavailable library generator.
    let live = 1;
    let dead = 1;
    let lx = 8;
    let lz = 8;
    let mx = 2 + (nx as i32 - 1) / lx;
    let mz = 2 + (nz as i32 - 1) / lz;
    let bd = beamData(wmin, nt, dt, ft, f);
    let jx = nint(((ray.rs[0].x - fx) / dx / lx as f32) as f64);
    let jz = nint(((ray.rs[0].z - fz) / dz / lz as f32) as f64);
    let mut cells = Cells {
        lx,
        mx,
        nx,
        dx,
        fx,
        lz,
        mz,
        nz,
        dz,
        fz,
        live,
        dead,
        wmin,
        lmin,
        cell: vec![vec![Cell::default(); mz as usize]; mx as usize],
        ray,
        bd,
        g,
    };
    setCell(&mut cells, jx, jz);
    accCell(&mut cells, jx, jz);
}
fn setCell(cells: &mut Cells, jx: i32, jz: i32) {
    if jx < 0 || jx >= cells.mx || jz < 0 || jz >= cells.mz {
        return;
    }
    if cells.cell[jx as usize][jz as usize].live == cells.live {
        return;
    }
    cells.cell[jx as usize][jz as usize].live = cells.live;
    if cellTimeAmp(cells, jx, jz) {
        setCell(cells, jx + 1, jz);
        setCell(cells, jx - 1, jz);
        setCell(cells, jx, jz + 1);
        setCell(cells, jx, jz - 1);
    }
}
fn accCell(cells: &mut Cells, jx: i32, jz: i32) {
    if jx < 0 || jx >= cells.mx - 1 || jz < 0 || jz >= cells.mz - 1 {
        return;
    }
    if cells.cell[jx as usize][jz as usize].dead == cells.dead {
        return;
    }
    cells.cell[jx as usize][jz as usize].dead = cells.dead;
    if cells.cell[jx as usize][jz as usize].live != cells.live {
        return;
    }
    if cells.cell[jx as usize + 1][jz as usize].live == cells.live
        && cells.cell[jx as usize][jz as usize + 1].live == cells.live
        && cells.cell[jx as usize + 1][jz as usize + 1].live == cells.live
    {
        cellBeam(cells, jx as usize, jz as usize);
    }
    accCell(cells, jx + 1, jz);
    accCell(cells, jx - 1, jz);
    accCell(cells, jx, jz + 1);
    accCell(cells, jx, jz - 1);
}

fn cellTimeAmp(cells: &mut Cells, jx: i32, jz: i32) -> bool {
    let tmax = cells.ray.rs[cells.ray.rs.len() - 1].t;
    let xc = cells.fx + (jx * cells.lx) as f32 * cells.dx;
    let zc = cells.fz + (jz * cells.lz) as f32 * cells.dz;
    let irs = nearestRayStep(cells.ray, xc, zc);
    let rs = cells.ray.rs[irs];
    let t0 = rs.t;
    let x0 = rs.x;
    let z0 = rs.z;
    let wmin = cells.wmin;
    let lmin = cells.lmin;
    let kmah = rs.kmah;
    let v = rs.v;
    let q1 = rs.q1;
    let p1 = rs.p1;
    let q2 = rs.q2;
    let p2 = rs.p2;
    let e = ((wmin * lmin) * lmin);
    let es = (e * e);
    let scale = ((v * v) / ((q2 * q2) + ((es * q1) * q1)));
    let vvmr = (((p2 * q2) + ((es * p1) * q1)) * scale);
    let vvmi = ((-e) * scale);
    let c = rs.c;
    let s = rs.s;
    let ov = ((1.0f64 / (v as f64)) as f32);
    let px = (s * ov);
    let pz = (c * ov);
    let pzpz = (pz * pz);
    let pxpx = (px * px);
    let pxpz = (px * pz);
    let dvdx = rs.dvdx;
    let dvdz = rs.dvdz;
    let dvds = ((s * dvdx) + (c * dvdz));
    let dvdn = ((c * dvdx) - (s * dvdz));
    let wxxr = (((((pzpz * vvmr) as f64) - ((2.0f64 * (pxpz as f64)) * (dvdn as f64)))
        - ((pxpx * dvds) as f64)) as f32);
    let wxzr = (((pxpx - pzpz) * dvdn) - (pxpz * (vvmr + dvds)));
    let wzzr = (((((pxpx * vvmr) as f64) + ((2.0f64 * (pxpz as f64)) * (dvdn as f64)))
        - ((pzpz * dvds) as f64)) as f32);
    let wxxi = (pzpz * vvmi);
    let wxzi = ((-pxpz) * vvmi);
    let wzzi = (pxpx * vvmi);
    let x = (xc - x0);
    let z = (zc - z0);
    let tr = ((((t0 as f64)
        + (((px as f64) + ((0.5f64 * (wxxr as f64)) * (x as f64))) * (x as f64)))
        + ((((pz + (wxzr * x)) as f64) + ((0.5f64 * (wzzr as f64)) * (z as f64))) * (z as f64)))
        as f32);
    let ti = (((((0.5f64 * (wxxi as f64)) * (x as f64)) * (x as f64))
        + ((((wxzi * x) as f64) + ((0.5f64 * (wzzi as f64)) * (z as f64))) * (z as f64)))
        as f32);
    let scale = (((es / ((q2 * q2) + ((es * q1) * q1))) as f64).powf(0.25f64) as f32);
    let phase = ((0.5f64
        * ((q2 as f64).atan2(((q1 * e) as f64))
            + ((2.0f64 * PI) * (((kmah + 1i32) / 2i32) as f64)))) as f32);
    let ar = (((scale as f64) * (phase as f64).cos()) as f32);
    let ai = (((scale as f64) * (phase as f64).sin()) as f32);

    let cell = &mut cells.cell[jx as usize][jz as usize];
    cell.tr = tr;
    cell.ti = ti;
    cell.ar = ar;
    cell.ai = ai;
    wmin * ti > -5.0 && tr <= tmax
}

fn cellBeam(cells: &mut Cells, jx: usize, jz: usize) {
    let lx = cells.lx as usize;
    let lz = cells.lz as usize;
    let nx = cells.nx;
    let nz = cells.nz;
    let cell = &cells.cell;
    let bd = &cells.bd;
    let ntr = bd.ntr as i32;
    let nti = bd.nti as i32;
    let dtr = bd.dtr;
    let ftr = bd.ftr;
    let dti = bd.dti;
    let fti = bd.fti;
    let cf = &bd.cf;
    let odtr = (1.0 / dtr as f64) as f32;
    let odti = (1.0 / dti as f64) as f32;
    let t00r = cell[jx][jz].tr;
    let t01r = cell[jx][jz + 1].tr;
    let t10r = cell[jx + 1][jz].tr;
    let t11r = cell[jx + 1][jz + 1].tr;
    let t00i = cell[jx][jz].ti;
    let t01i = cell[jx][jz + 1].ti;
    let t10i = cell[jx + 1][jz].ti;
    let t11i = cell[jx + 1][jz + 1].ti;
    let a00r = cell[jx][jz].ar;
    let a01r = cell[jx][jz + 1].ar;
    let a10r = cell[jx + 1][jz].ar;
    let a11r = cell[jx + 1][jz + 1].ar;
    let a00i = cell[jx][jz].ai;
    let a01i = cell[jx][jz + 1].ai;
    let a10i = cell[jx + 1][jz].ai;
    let a11i = cell[jx + 1][jz + 1].ai;
    let kxlo = jx * lx;
    let kxhi = (kxlo + lx).min(nx);
    let kzlo = jz * lz;
    let kzhi = (kzlo + lz).min(nz);
    let xdelta = (1.0 / lx as f64) as f32;
    let zdelta = (1.0 / lz as f64) as f32;
    let dtx0r = (t10r - t00r) * xdelta;
    let dtx1r = (t11r - t01r) * xdelta;
    let dtx0i = (t10i - t00i) * xdelta;
    let dtx1i = (t11i - t01i) * xdelta;
    let dax0r = (a10r - a00r) * xdelta;
    let dax1r = (a11r - a01r) * xdelta;
    let dax0i = (a10i - a00i) * xdelta;
    let dax1i = (a11i - a01i) * xdelta;
    let mut tx0r = t00r;
    let mut tx1r = t01r;
    let mut tx0i = t00i;
    let mut tx1i = t01i;
    let mut ax0r = a00r;
    let mut ax1r = a01r;
    let mut ax0i = a00i;
    let mut ax1i = a01i;
    for kx in kxlo..kxhi {
        let dtxzr = (tx1r - tx0r) * zdelta;
        let dtxzi = (tx1i - tx0i) * zdelta;
        let daxzr = (ax1r - ax0r) * zdelta;
        let daxzi = (ax1i - ax0i) * zdelta;
        let mut txzr = tx0r;
        let mut txzi = tx0i;
        let mut axzr = ax0r;
        let mut axzi = ax0i;
        for kz in kzlo..kzhi {
            let tin = (txzi - fti) * odti;
            let iti = tin as i32;
            if iti < 0 || iti >= nti - 1 {
                continue;
            }
            let cf0 = &cf[iti as usize];
            let cf1 = &cf[iti as usize + 1];
            let tifrac = tin - iti as f32;
            let mtifrac = (1.0 - tifrac as f64) as f32;
            let trn = (txzr - ftr) * odtr;
            let itr = trn as i32;
            if itr < 0 || itr >= ntr - 1 {
                continue;
            }
            let trfrac = trn - itr as f32;
            let mtrfrac = (1.0 - trfrac as f64) as f32;
            let itr = itr as usize;
            let cf0r = mtrfrac * cf0[itr].r + trfrac * cf0[itr + 1].r;
            let cf1r = mtrfrac * cf1[itr].r + trfrac * cf1[itr + 1].r;
            let cfr = mtifrac * cf0r + tifrac * cf1r;
            let cf0i = mtrfrac * cf0[itr].i + trfrac * cf0[itr + 1].i;
            let cf1i = mtrfrac * cf1[itr].i + trfrac * cf1[itr + 1].i;
            let cfi = mtifrac * cf0i + tifrac * cf1i;
            cells.g[kx][kz] += axzr * cfr - axzi * cfi;
            txzr += dtxzr;
            txzi += dtxzi;
            axzr += daxzr;
            axzi += daxzi;
        }
        tx0r += dtx0r;
        tx1r += dtx1r;
        tx0i += dtx0i;
        tx1i += dtx1i;
        ax0r += dax0r;
        ax1r += dax1r;
        ax0i += dax0i;
        ax1i += dax1i;
    }
}
