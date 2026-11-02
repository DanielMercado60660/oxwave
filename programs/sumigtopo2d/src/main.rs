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

//! SUMIGTOPO2D - Kirchhoff depth migration from a variable recording surface.
//! Port of su/main/migration_inversion/sumigtopo2d.c (Revision 1.7).
//! Credits: Zhenyue Liu, 03/01/95, Colorado School of Mines.
//! Trino Salinas, 07/01/96, Colorado School of Mines: nonflat recording surface.
//! Surface subroutines adapted from Dave Hale, 09/17/91.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::{
    complex::Complex,
    cwp::{abs, max, min, PI},
};
use oxwave_f_cubicspline::csplin;
use oxwave_inverse_interp::yxtoxy;
use oxwave_pfa_fft::{npfaro, pfacr, pfarc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUMIGTOPO2D - Kirchhoff Depth Migration of 2D postack/prestack data\t",
    "\t     from the (variable topography) recording surface\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "    sumigtopo2d  infile=  outfile=  [parameters] \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " infile=stdin\t\tfile for input seismic traces\t\t\t",
    " outfile=stdout\tfile for common offset migration output  \t",
    " ttfile\t\tfile for input traveltime tables\t\t",
    "   The following 9 parameters describe traveltime tables:\t\t",
    " fzt \t\t\tfirst depth sample in traveltime table\t\t",
    " nzt \t\t\tnumber of depth samples in traveltime table\t",
    " dzt\t\t\tdepth interval in traveltime table\t\t",
    " fxt\t\t\tfirst lateral sample in traveltime table\t",
    " nxt\t\t\tnumber of lateral samples in traveltime table\t",
    " dxt\t\t\tlateral interval in traveltime table\t\t",
    " fs \t\t\tx-coordinate of first source\t\t\t",
    " ns \t\t\tnumber of sources\t\t\t\t",
    " ds \t\t\tx-coordinate increment of sources\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " fxi                   x-coordinate of the first input trace           ",
    " dxi                   horizontal spacing of input data                ",
    " nxi                   number of input trace locations in surface      ",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt= or from header (dt) \ttime sampling interval of input data\t",
    " ft= or from header (ft) \tfirst time sample of input data\t\t",
    " dxm= or from header (d2) \tsampling interval of midpoints \t\t",
    " surf=\"0,0;99999,0\"  Recording surface \"x1,z1;x2,z2;x3,z3;...\"     ",
    " fzo=fzt               z-coordinate of first point in output trace \t",
    " dzo=0.2*dzt\t\tvertical spacing of output trace \t\t",
    " nzo=5*(nzt-1)+1 \tnumber of points in output trace\t\t",
    " fxo=fxt               x-coordinate of first output trace \t\t",
    " dxo=0.5*dxt\t\thorizontal spacing of output trace \t\t",
    " nxo=2*(nxt-1)+1  \tnumber of output traces \t\t\t",
    " off0=0               \tfirst offest in output \t\t\t\t",
    " doff=99999\t\toffset increment in output \t\t\t",
    " noff=1       \t\tnumber of offsets in output \t\t\t",
    " fmax=0.25/dt\t\tfrequency-highcut for input traces\t\t",
    " offmax=99999\t\tmaximum absolute offset allowed in migration \t",
    " aperx=nxt*dxt/2  \tmigration lateral aperature \t\t\t",
    " angmax=60\t\tmigration angle aperature from vertical \t",
    " v0=1500(m/s)\t\treference velocity value at surface\t\t",
    " dvz=0.0  \t\treference velocity vertical gradient\t\t",
    " ls=1\t                flag for line source\t\t\t\t",
    " jpfile=stderr\t\tjob print file name \t\t\t\t",
    " mtr=100  \t\tprint verbal information at every mtr traces\t",
    " ntr=100000\t\tmaximum number of input traces to be migrated\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " 1. Traveltime tables were generated by program rayt2dtopo (or any \t",
    "    other one that considers topography )on relatively coarse grids,\t",
    "    with dimension ns*nxt*nzt. In the migration process, traveltimes\t",
    "    are interpolated into shot/gephone positions and output grids.\t",
    " 2. Input seismic traces must be SU format and can be any type of \t",
    "    gathers (common shot, common offset, common CDP, and so on).\t",
    " 3. Migrated traces are output in CDP gathers if velocity analysis\t",
    "    is required, with dimension nxo*noff*nzo.  \t\t\t",
    " 4. If the offset value of an input trace is not in the offset array \t",
    "    of output, the nearest one in the array is chosen. \t\t",
    " 5. Amplitudes are computed using the reference velocity profile, v(z),",
    "    specified by the parameters v0= and dvz=.\t\t\t\t",
    " 6. Input traces must specify source and receiver positions via the header",
    "    fields tr.sx and tr.gx. Offset is computed automatically.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// C's alloc2float/ealloc2float stores all rows in one contiguous block.
// In particular, the traveltime-table fread below spans every row.
#[derive(Clone)]
struct Table {
    data: Vec<f32>,
    n1: usize,
}

impl std::ops::Index<usize> for Table {
    type Output = [f32];

    fn index(&self, i2: usize) -> &Self::Output {
        &self.data[i2 * self.n1..(i2 + 1) * self.n1]
    }
}

impl std::ops::IndexMut<usize> for Table {
    fn index_mut(&mut self, i2: usize) -> &mut Self::Output {
        &mut self.data[i2 * self.n1..(i2 + 1) * self.n1]
    }
}

#[derive(Clone, Copy, Default)]
struct SurfaceSegment {
    x: f32,
    z: f32,
    s: f32,
    c: f32,
}

struct Surface {
    ss: Vec<SurfaceSegment>,
}

fn decode_surfaces() -> Vec<(Vec<f32>, Vec<f32>)> {
    let mut nr = countparname("surf");
    if nr == 0 {
        nr = 1;
    }
    let mut surfaces = Vec::new();
    for ir in 0..nr {
        let mut s = String::new();
        if getnparstring(ir + 1, "surf", &mut s) == 0 {
            s = "0,0;99999,0".to_string();
        }
        match decode_surface(&s) {
            Some(surface) => surfaces.push(surface),
            None => err!("Surface number %d specified incorrectly!\n", ir + 1),
        }
    }
    surfaces
}

fn decode_surface(string: &str) -> Option<(Vec<f32>, Vec<f32>)> {
    let string = string.split('\0').next().unwrap_or("");
    let mut string = string.as_bytes().to_vec();
    string.push(0);
    // C first counts tokens using strtok, which replaces only the delimiter
    // ending each token. Delimiters skipped before a token remain in the buffer.
    let mut first = None;
    let mut nxz = 0;
    let mut s = 0;
    while string[s] != 0 {
        while matches!(string[s], b',' | b';') {
            s += 1;
        }
        if string[s] == 0 {
            break;
        }
        first.get_or_insert(s);
        nxz += 1;
        while string[s] != 0 && !matches!(string[s], b',' | b';') {
            s += 1;
        }
        if string[s] != 0 {
            string[s] = 0;
            s += 1;
        }
    }
    if nxz % 2 != 0 {
        return None;
    }
    let nxz = nxz / 2;
    if nxz < 2 {
        return None;
    }
    let mut x = vec![0.0; nxz];
    let mut z = vec![0.0; nxz];
    let mut s = first.unwrap();
    for ixz in 0..nxz {
        // Conversion advances by strlen(s)+1, not by the strtok token starts.
        // For example, "-10,,5;30,,7" becomes the flat surface (-10,0),(30,0).
        let end = s + string[s..].iter().position(|&c| c == 0).unwrap();
        x[ixz] = atof(std::str::from_utf8(&string[s..end]).unwrap()) as f32;
        s = end + 1;
        let end = s + string[s..].iter().position(|&c| c == 0).unwrap();
        z[ixz] = atof(std::str::from_utf8(&string[s..end]).unwrap()) as f32;
        s = end + 1;
    }
    Some((x, z))
}

fn makesurf(dsmax: f32, surfaces: &[(Vec<f32>, Vec<f32>)]) -> Vec<Surface> {
    let mut rr = Vec::new();
    for (xu, zu) in surfaces {
        let nu = xu.len();
        let mut u = vec![0.0; nu];
        for iu in 0..nu {
            u[iu] = iu as f32;
        }
        let mut xud = vec![[0.0; 4]; nu];
        csplin(&u, xu, &mut xud);
        let mut zud = vec![[0.0; 4]; nu];
        csplin(&u, zu, &mut zud);

        let nuu = 20 * nu;
        let duu = (u[nu - 1] - u[0]) / (nuu - 1) as f32;
        let mut s = vec![0.0f32; nuu];
        let mut xlast = xu[0];
        let mut zlast = zu[0];
        let mut uu = duu;
        for iuu in 1..nuu {
            let mut x = [0.0f32];
            let mut z = [0.0f32];
            intcub(0, &u, &xud, &[uu], &mut x);
            intcub(0, &u, &zud, &[uu], &mut z);
            let dx = x[0] - xlast;
            let dz = z[0] - zlast;
            s[iuu] = (s[iuu - 1] as f64 + ((dx * dx + dz * dz) as f64).sqrt()) as f32;
            xlast = x[0];
            zlast = z[0];
            uu += duu;
        }
        let ns = (1.0 + s[nuu - 1] / dsmax) as i32;
        let ds = s[nuu - 1] / ns as f32;
        let fs = (0.5 * ds as f64) as f32;
        let mut us = vec![0.0; ns as usize];
        yxtoxy(duu, 0.0, &s, ds, fs, 0.0, (nu - 1) as f32, &mut us);
        let mut ss = vec![SurfaceSegment::default(); ns as usize];
        for is in 0..ns as usize {
            let mut rsx = [0.0];
            let mut rsz = [0.0];
            let mut rsxd = [0.0];
            let mut rszd = [0.0];
            intcub(0, &u, &xud, &[us[is]], &mut rsx);
            intcub(0, &u, &zud, &[us[is]], &mut rsz);
            intcub(1, &u, &xud, &[us[is]], &mut rsxd);
            intcub(1, &u, &zud, &[us[is]], &mut rszd);
            ss[is].x = rsx[0];
            ss[is].z = rsz[0];
            ss[is].c =
                (rsxd[0] as f64 / ((rsxd[0] * rsxd[0] + rszd[0] * rszd[0]) as f64).sqrt()) as f32;
            ss[is].s =
                (-rszd[0] as f64 / ((rsxd[0] * rsxd[0] + rszd[0] * rszd[0]) as f64).sqrt()) as f32;
        }
        rr.push(Surface { ss });
    }
    rr
}

#[allow(clippy::too_many_arguments)]
fn zcoor_surfaces(
    fx: f32,
    dx: f32,
    nx: i32,
    fxs: f32,
    dxs: f32,
    nxs: i32,
    srf: &[Surface],
    szif: &mut [f32],
    sz: &mut [f32],
    nangl: &mut [f32],
) {
    let ss = &srf[0].ss;
    let ns = ss.len();
    let mut ixi = 0;
    for ix in 0..nx as usize {
        let x = fx + ix as f32 * dx;
        for ik in ixi..ns - 1 {
            if ss[ik].x <= x && ss[ik + 1].x >= x {
                let az = ss[ik].z - ss[ik + 1].z;
                let ax0 = ss[ik + 1].x - x;
                let ax = ss[ik + 1].x - ss[ik].x;
                szif[ix] = ax0 * az / ax + ss[ik + 1].z;
                let az = ss[ik].s - ss[ik + 1].s;
                let temp = ax0 * az / ax + ss[ik + 1].s;
                nangl[ix] = (temp as f64).asin() as f32;
                ixi = ik;
            }
        }
    }
    ixi = 0;
    for ix in 0..nxs as usize {
        let x = fxs + ix as f32 * dxs;
        for ik in ixi..ns - 1 {
            if ss[ik].x <= x && ss[ik + 1].x >= x {
                let az = ss[ik].z - ss[ik + 1].z;
                let ax0 = ss[ik + 1].x - x;
                let ax = ss[ik + 1].x - ss[ik].x;
                sz[ix] = ax0 * az / ax + ss[ik + 1].z;
                ixi = ik;
            }
        }
    }
    let nx = nx as usize;
    let nxs = nxs as usize;
    // C lines 728-733 extrapolate using three points even for shorter arrays.
    // Retain the interpolated values when extrapolation lacks enough points.
    // This bounds fix leaves the four round-1 surface cases byte-identical to SU.
    if nx >= 3 {
        szif[0] = 2.0 * szif[1] - szif[2];
        szif[nx - 1] = 2.0 * szif[nx - 2] - szif[nx - 3];
    }
    if nxs >= 3 {
        sz[0] = 2.0 * sz[1] - sz[2];
    }
    if nx >= 3 {
        nangl[0] = 2.0 * nangl[1] - nangl[2];
    }
    if nxs >= 3 {
        sz[nxs - 1] = 2.0 * sz[nxs - 2] - sz[nxs - 3];
    }
    if nx >= 3 {
        nangl[nx - 1] = 2.0 * nangl[nx - 2] - nangl[nx - 3];
    }
}

#[allow(clippy::too_many_arguments)]
fn resit(nx: i32, fx: f32, dx: f32, nz: i32, nr: i32, dr: f32, tb: &Table, t: &mut Table, x0: f32) {
    for ix in 0..nx as usize {
        let xi = fx + ix as f32 * dx - x0;
        let ar = xi.abs() / dr;
        let mut jr = ar as i32;
        let sr = ar - jr as f32;
        let sr0 = (1.0 - sr as f64) as f32;
        if jr > nr - 2 {
            jr = nr - 2;
        }
        for iz in 0..nz as usize {
            t[ix][iz] -= sr0 * tb[jr as usize][iz] + sr * tb[jr as usize + 1][iz];
        }
    }
}

fn sum2(nx: i32, nz: i32, a1: f32, a2: f32, t1: &Table, t2: &Table, t: &mut Table) {
    for ix in 0..nx as usize {
        for iz in 0..nz as usize {
            t[ix][iz] = a1 * t1[ix][iz] + a2 * t2[ix][iz];
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn timeb(
    nr: i32,
    nz: i32,
    dr: f32,
    dz: f32,
    fz: f32,
    z0: f32,
    a: f32,
    v0: f32,
    t: &mut Table,
    p: &mut Table,
    sig: &mut Table,
    ang: &mut Table,
) {
    if a == 0.0 {
        let mut r = 0.0f32;
        for ir in 0..nr as usize {
            let mut z = fz - z0;
            for iz in 0..nz as usize {
                let mut rou = ((r * r + z * z) as f64).sqrt() as f32;
                if rou < dz {
                    rou = dz;
                }
                t[ir][iz] = rou / v0;
                p[ir][iz] = r / (rou * v0);
                sig[ir][iz] = v0 * rou;
                ang[ir][iz] = ((r / rou) as f64).asin() as f32;
                z += dz;
            }
            r += dr;
        }
    } else {
        let oa = (1.0 / a as f64) as f32;
        let zc = v0 * oa;
        let mut r = 0.0f32;
        for ir in 0..nr as usize {
            let mut z = fz + zc - z0;
            for iz in 0..nz as usize {
                // The C initial rou=sqrt(r*r+z*z) is overwritten on this branch.
                let v = v0 + a * (z - zc + z0);
                if ir == 0 {
                    // log(...)*oa is evaluated in double before assignment.
                    t[ir][iz] = (((v / v0) as f64).ln() * oa as f64) as f32;
                    p[ir][iz] = 0.0;
                    ang[ir][iz] = 0.0;
                    sig[ir][iz] = (0.5
                        * (z as f64 + 0.1 * dz as f64 - zc as f64 + z0 as f64)
                        * (v0 + v) as f64) as f32;
                } else {
                    let rc = ((r * r + z * z - zc * zc) as f64 / (2.0 * r as f64)) as f32;
                    let rou = ((zc * zc + rc * rc) as f64).sqrt() as f32;
                    t[ir][iz] = ((((v * (rou + rc)) / (v0 * (rou + rc - r))) as f64).ln()
                        * oa as f64) as f32;
                    p[ir][iz] = (((rou * rou - rc * rc) as f64).sqrt() / (rou * v0) as f64) as f32;
                    let mut temp = v0 * p[ir][iz];
                    if temp > 1.0 {
                        temp = 1.0;
                    }
                    ang[ir][iz] = (temp as f64).asin() as f32;
                    sig[ir][iz] = a * rou * r;
                }
                z += dz;
            }
            r += dr;
        }
    }
}

fn table(n1: i32, n2: i32) -> Table {
    Table {
        data: vec![0.0; n1 as usize * n2 as usize],
        n1: n1 as usize,
    }
}

// The C does not check outfile/jpfile fopen failures before using FILE*.
// Preserve its SIGSEGV without calling C or introducing a Rust null dereference
// (which would be undefined behavior). The assembly explicitly faults on the
// unmapped zero page; it does not access any Rust allocation.
fn null_file() -> ! {
    #[cfg(target_arch = "aarch64")]
    unsafe {
        std::arch::asm!("str xzr, [{address}]", address = in(reg) 0usize);
    }
    #[cfg(target_arch = "x86_64")]
    unsafe {
        std::arch::asm!("mov qword ptr [{address}], 0", address = in(reg) 0usize);
    }
    exit(139)
}

#[allow(clippy::too_many_arguments)]
fn mig2d(
    trace: &mut [f32],
    nt: i32,
    ft: f32,
    dt: f32,
    sx: f32,
    gx: f32,
    mig: &mut Table,
    aperx: f32,
    nx: i32,
    fx: f32,
    dx: f32,
    nz: f32,
    fz: f32,
    dz: f32,
    ls: i32,
    mtmax: i32,
    dxm: f32,
    fmax: f32,
    angmax: f32,
    tbs: &Table,
    tbg: &Table,
    pbs: &Table,
    pbg: &Table,
    angbs: &Table,
    angbg: &Table,
    nr: i32,
    tsum: &Table,
    nzt: i32,
    fzt: f32,
    dzt: f32,
    nxt: i32,
    fxt: f32,
    dxt: f32,
    nangls: f32,
    nanglg: f32,
    sigbs: &Table,
    sigbg: &Table,
    filter: &mut Option<Filter>,
) {
    // Match the fixed C: depth interpolation needs two table samples.
    if nzt < 2 {
        err!("nzt must be at least 2 for depth interpolation");
    }
    let mut tmt = table(nzt, nxt);
    let mut ampt = table(nzt, nxt);
    let mut ampti = table(nzt, nxt);
    let mut tm = vec![0.0f32; nzt as usize];
    let mut tzt = vec![0.0f32; nzt as usize];
    let mut amp = vec![0.0f32; nzt as usize];
    let mut ampi = vec![0.0f32; nzt as usize];
    let mut zpt = vec![0.0f32; nxt as usize];
    let mut trf = vec![0.0f32; (nt + 2 * mtmax) as usize];
    let odt = (1.0 / dt as f64) as f32;
    let z0 = (fz - fzt) / dzt;
    let rdz = dz / dzt;
    let pmin = (1.0 / (2.0 * dxm as f64 * fmax as f64)) as f32;
    filt(trace, nt, dt, fmax, ls, mtmax, &mut trf, filter);
    let xm = (0.5 * (sx + gx) as f64) as f32;
    let rxz = if angmax == 90.0 {
        0.0
    } else {
        (1.0 / (angmax as f64 * PI / 180.0).tan()) as f32
    };
    let mut nxtf = ((xm - aperx - fxt) / dxt) as i32;
    if nxtf < 0 {
        nxtf = 0;
    }
    let mut nxte = ((xm + aperx - fxt) / dxt + 1.0) as i32;
    if nxte >= nxt {
        nxte = nxt - 1;
    }

    /* compute amplitudes and filter length */
    for ix in nxtf..=nxte {
        let x = fxt + ix as f32 * dxt;
        let dis = if xm >= x { xm - x } else { x - xm };
        let mut izt0 = (((dis - dxt) * rxz - fzt) / dzt - 1.0) as i32;
        if izt0 < 0 {
            izt0 = 0;
        }
        if izt0 >= nzt {
            izt0 = nzt - 1;
        }
        let ar = if sx >= x {
            (sx - x) / dx
        } else {
            (x - sx) / dx
        };
        let mut jrs = ar as i32;
        if jrs > nr - 2 {
            jrs = nr - 2;
        }
        let srs = ar - jrs as f32;
        let srs0 = (1.0 - srs as f64) as f32;
        let ar = if gx >= x {
            (gx - x) / dx
        } else {
            (x - gx) / dx
        };
        let mut jrg = ar as i32;
        if jrg > nr - 2 {
            jrg = nr - 2;
        }
        let srg = ar - jrg as f32;
        let srg0 = (1.0 - srg as f64) as f32;
        let sigp = if (sx - x) * (gx - x) > 0.0 {
            1.0f32
        } else {
            -1.0
        };
        let sigs = if sx - x < 0.0 { 1.0f32 } else { -1.0 };
        let sigg = if gx - x < 0.0 { 1.0f32 } else { -1.0 };
        let ix = ix as usize;
        let jrs = jrs as usize;
        let jrg = jrg as usize;
        zpt[ix] = fzt + (nzt - 1) as f32 * dzt;
        for iz in izt0 as usize..nzt as usize {
            let mut angs = srs0 * angbs[jrs][iz] + srs * angbs[jrs + 1][iz];
            angs = abs(sigs * angs - nangls);
            let mut angg = srg0 * angbg[jrg][iz] + srg * angbg[jrg + 1][iz];
            angg = abs(sigg * angg - nanglg);
            let sigmas = srs0 * sigbs[jrs][iz] + srs * sigbs[jrs + 1][iz];
            let sigmag = srg0 * sigbg[jrg][iz] + srg * sigbg[jrg + 1][iz];
            let mut ampd =
                ((angs as f64).cos() * sigmag as f64 + (angg as f64).cos() * sigmas as f64) as f32;
            ampd = (ampd as f64
                / ((nangls as f64).cos()
                    * (nanglg as f64).cos()
                    * ((sigmas * sigmag) as f64).sqrt())) as f32;
            /* Filter of 90 degrees to the operator (Peels,1988) */
            if abs(angs) as f64 >= PI / 2.0 || abs(angg) as f64 >= PI / 2.0 {
                ampd = 0.0;
            }
            if ampd < 0.0 {
                ampd = -ampd;
            }
            ampt[ix][iz] = ampd;
            let mut pd = srs0 * pbs[jrs][iz]
                + srs * pbs[jrs + 1][iz]
                + sigp * (srg0 * pbg[jrg][iz] + srg * pbg[jrg + 1][iz]);
            if pd < 0.0 {
                pd = -pd;
            }
            let mut temp = pd * dxm * odt;
            if temp < 1.0 {
                temp = 1.0;
            }
            if temp > mtmax as f32 {
                temp = mtmax as f32;
            }
            ampti[ix][iz] = ampd / (temp * temp);
            tmt[ix][iz] = temp;
            if pd < pmin && zpt[ix] as f64 > fzt as f64 + (nzt as f64 - 1.1) * dzt as f64 {
                zpt[ix] = fzt + iz as f32 * dzt;
            }
        }
    }
    let mut nxf = (((xm - aperx - fx) / dx) as f64 + 0.5) as i32;
    if nxf < 0 {
        nxf = 0;
    }
    let mut nxe = (((xm + aperx - fx) / dx) as f64 + 0.5) as i32;
    if nxe >= nx {
        nxe = nx - 1;
    }

    /* interpolate amplitudes and filter length along lateral */
    for ix in nxf..=nxe {
        let x = fx + ix as f32 * dx;
        let dis = if xm >= x { xm - x } else { x - xm };
        let mut izt0 = ((dis * rxz - fzt) / dzt) as i32;
        if izt0 < 0 {
            izt0 = 0;
        }
        if izt0 >= nzt {
            izt0 = nzt - 1;
        }
        let mut iz0 = ((dis * rxz - fz) / dz) as i32;
        if iz0 < 0 {
            iz0 = 0;
        }
        if iz0 as f32 >= nz {
            iz0 = (nz - 1.0) as i32;
        }
        let mut ax = (x - fxt) / dxt;
        let mut jx = ax as i32;
        ax -= jx as f32;
        if ax as f64 <= 0.01 {
            ax = 0.0;
        }
        if ax as f64 >= 0.99 {
            ax = 1.0;
        }
        let ax0 = (1.0 - ax as f64) as f32;
        if jx > nxte - 1 {
            jx = nxte - 1;
        }
        if jx < nxtf {
            jx = nxtf;
        }
        let ar = if sx >= x {
            (sx - x) / dx
        } else {
            (x - sx) / dx
        };
        let mut jrs = ar as i32;
        if jrs > nr - 2 {
            jrs = nr - 2;
        }
        let srs = ar - jrs as f32;
        let srs0 = (1.0 - srs as f64) as f32;
        let ar = if gx >= x {
            (gx - x) / dx
        } else {
            (x - gx) / dx
        };
        let mut jrg = ar as i32;
        if jrg > nr - 2 {
            jrg = nr - 2;
        }
        let srg = ar - jrg as f32;
        let srg0 = (1.0 - srg as f64) as f32;
        let jx = jx as usize;
        let jrs = jrs as usize;
        let jrg = jrg as usize;
        let ix = ix as usize;
        for iz in izt0 as usize..nzt as usize {
            tzt[iz] = ax0 * tsum[jx][iz]
                + ax * tsum[jx + 1][iz]
                + srs0 * tbs[jrs][iz]
                + srs * tbs[jrs + 1][iz]
                + srg0 * tbg[jrg][iz]
                + srg * tbg[jrg + 1][iz];
            amp[iz] = ax0 * ampt[jx][iz] + ax * ampt[jx + 1][iz];
            ampi[iz] = ax0 * ampti[jx][iz] + ax * ampti[jx + 1][iz];
            tm[iz] = ax0 * tmt[jx][iz] + ax * tmt[jx + 1][iz];
        }
        let mut nzp = (((ax0 * zpt[jx] + ax * zpt[jx + 1] - fz) / dz) as f64 + 1.5) as i32;
        if nzp < iz0 {
            nzp = iz0;
        }
        if nzp as f32 > nz {
            nzp = nz as i32;
        }
        /* interpolate along depth if operator aliasing */
        for iz in iz0..nzp {
            let az = z0 + iz as f32 * rdz;
            let mut jz = az as i32;
            if jz >= nzt - 1 {
                jz = nzt - 2;
            }
            let sz = az - jz as f32;
            let sz0 = (1.0 - sz as f64) as f32;
            let jz = jz as usize;
            let td = sz0 * tzt[jz] + sz * tzt[jz + 1];
            let at = (td - ft) * odt + mtmax as f32;
            let jt = at as i32;
            if jt > mtmax && jt < nt + mtmax - 1 {
                let ampd = sz0 * ampi[jz] + sz * ampi[jz + 1];
                // The leading 0.5 promotes each following multiplication to
                // double only after the float products have been evaluated.
                let mt = (0.5 + (sz0 * tm[jz]) as f64 + (sz * tm[jz + 1]) as f64) as i32;
                let res = at - jt as f32;
                let res0 = (1.0 - res as f64) as f32;
                let temp = (res0 as f64
                    * (-trf[(jt - mt) as usize] as f64 + 2.0 * trf[jt as usize] as f64
                        - trf[(jt + mt) as usize] as f64)
                    + res as f64
                        * (-trf[(jt - mt + 1) as usize] as f64
                            + 2.0 * trf[(jt + 1) as usize] as f64
                            - trf[(jt + mt + 1) as usize] as f64))
                    * ampd as f64;
                let temp = temp as f32;
                mig[ix][iz as usize] += temp;
            }
        }
        /* interpolate along depth if not operator aliasing */
        let mut iz = nzp;
        while (iz as f32) < nz {
            let az = z0 + iz as f32 * rdz;
            let mut jz = az as i32;
            if jz >= nzt - 1 {
                jz = nzt - 2;
            }
            let sz = az - jz as f32;
            let sz0 = (1.0 - sz as f64) as f32;
            let jz = jz as usize;
            let td = sz0 * tzt[jz] + sz * tzt[jz + 1];
            let at = (td - ft) * odt;
            let jt = at as i32;
            if jt > 0 && jt < nt - 1 {
                let ampd = sz0 * amp[jz] + sz * amp[jz + 1];
                let res = at - jt as f32;
                let res0 = (1.0 - res as f64) as f32;
                let temp = (res0 * trace[jt as usize] + res * trace[jt as usize + 1]) * ampd;
                mig[ix][iz as usize] += temp;
            }
            iz += 1;
        }
    }
}

struct Filter {
    nfft: i32,
    itaper: i32,
    nw: i32,
    nwf: i32,
    taper: Vec<f32>,
    amp: Vec<f32>,
    ampi: Vec<f32>,
}

#[allow(clippy::too_many_arguments)]
fn filt(
    trace: &mut [f32],
    nt: i32,
    dt: f32,
    mut fmax: f32,
    ls: i32,
    m: i32,
    trf: &mut [f32],
    state: &mut Option<Filter>,
) {
    fmax = (fmax as f64 * (2.0 * PI)) as f32;
    let ftaper = (0.1 * fmax as f64) as f32;
    let const2 = (0.5 * 2.0f64.sqrt()) as f32;
    if state.is_none() {
        let nfft = npfaro(nt + 2 * m, 2 * (nt + 2 * m));
        if nfft >= SU_NFLTS as i32 || nfft >= 720720 {
            err!("Padded nt=%d -- too big", nfft);
        }
        let nw = nfft / 2 + 1;
        let dw = (2.0 * PI / (nfft as f32 * dt) as f64) as f32;
        let itaper = (0.5 + (ftaper / dw) as f64) as i32;
        let mut taper = vec![0.0f32; (2 * itaper + 1) as usize];
        for iw in -itaper..=itaper {
            let temp = (iw as f32 as f64 / (1.0 + itaper as f64)) as f32;
            taper[(iw + itaper) as usize] = (1.0 - temp) * (1.0 - temp) * (temp + 2.0) / 4.0;
        }
        let mut nwf = (0.5 + (fmax / dw) as f64) as i32;
        if nwf > nw - itaper - 1 {
            nwf = nw - itaper - 1;
        }
        let mut amp = vec![0.0f32; (nwf + itaper + 1) as usize];
        let mut ampi = vec![0.0f32; (nwf + itaper + 1) as usize];
        for iw in 1..=nwf + itaper {
            amp[iw as usize] = (((dw * iw as f32) as f64).sqrt() / nfft as f64) as f32;
            ampi[iw as usize] = (0.5 / (1.0 - ((iw as f32 * dw * dt) as f64).cos())) as f32;
        }
        *state = Some(Filter {
            nfft,
            itaper,
            nw,
            nwf,
            taper,
            amp,
            ampi,
        });
    }
    let f = state.as_ref().unwrap();
    let mut rt = vec![0.0f32; f.nfft as usize];
    let mut ct = vec![Complex::default(); f.nw as usize];
    rt[..nt as usize].copy_from_slice(&trace[..nt as usize]);
    pfarc(1, &rt, &mut ct);
    // C filt() can start before ct[0] for a wide taper. The minimal C fix
    // skips those invalid writes while preserving the taper's original index.
    for iw in max(0, f.nwf - f.itaper)..=f.nwf + f.itaper {
        let itemp = iw - (f.nwf - f.itaper);
        ct[iw as usize].r = f.taper[itemp as usize] * ct[iw as usize].r;
        ct[iw as usize].i = f.taper[itemp as usize] * ct[iw as usize].i;
    }
    for iw in f.nwf + f.itaper + 1..f.nw {
        ct[iw as usize].r = 0.0;
        ct[iw as usize].i = 0.0;
    }
    if ls == 0 {
        for iw in 0..=(f.nwf + f.itaper) as usize {
            let temp = -ct[iw].i * f.amp[iw];
            ct[iw].i = ct[iw].r * f.amp[iw];
            ct[iw].r = temp;
        }
    } else {
        for iw in 0..=(f.nwf + f.itaper) as usize {
            let temp = (ct[iw].r - ct[iw].i) * f.amp[iw] * const2;
            ct[iw].i = (ct[iw].r + ct[iw].i) * f.amp[iw] * const2;
            ct[iw].r = temp;
        }
    }
    pfacr(-1, &ct, &mut rt);
    for it in 0..nt as usize {
        trace[it] = rt[it];
    }
    for iw in 0..=(f.nwf + f.itaper) as usize {
        ct[iw].i = ct[iw].i * f.ampi[iw];
        ct[iw].r = ct[iw].r * f.ampi[iw];
    }
    pfacr(-1, &ct, &mut rt);
    for it in 0..m as usize {
        trf[it] = rt[(f.nfft - m) as usize + it];
    }
    for it in 0..(nt + m) as usize {
        trf[it + m as usize] = rt[it];
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* open input and output files */
    let mut infile = "stdin".to_string();
    let infp = if getparstring("infile", &mut infile) == 0 {
        stdin()
    } else {
        fopen(&infile, "r").unwrap_or_else(|| err!("cannot open infile=%s\n", infile))
    };
    let mut outfile = "stdout".to_string();
    let outfp = if getparstring("outfile", &mut outfile) == 0 {
        Some(stdout())
    } else {
        fopen(&outfile, "w")
    };
    efseeko(&infp, 0, SEEK_CUR);
    // C passes a null FILE* on an unsuccessful outfile fopen.
    let outfp = outfp.unwrap_or_else(|| null_file());
    efseeko(&outfp, 0, SEEK_CUR);
    let mut ttfile = String::new();
    if getparstring("ttfile", &mut ttfile) == 0 {
        err!("must specify ttfile!\n");
    }
    let ttfp = fopen(&ttfile, "r").unwrap_or_else(|| err!("cannot open ttfile=%s\n", ttfile));
    let mut jpfile = String::new();
    let jpfp = if getparstring("jpfile", &mut jpfile) == 0 {
        Some(stderr())
    } else {
        fopen(&jpfile, "w")
    };

    /* get information from the first header */
    if fgettr(&infp, &mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    if (dt as f64) < 0.0000001 {
        err!("dt must be positive!\n");
    }
    let mut ft = 0.0f32;
    if getparfloat("ft", &mut ft) == 0 {
        ft = (tr.delrt as f32 as f64 / 1000.0) as f32;
    }
    let mut dxm = 0.0f32;
    if getparfloat("dxm", &mut dxm) == 0 {
        dxm = tr.d2;
    }
    if (dxm as f64) < 0.0000001 {
        err!("dxm must be positive!\n");
    }

    /* get traveltime table parameters */
    macro_rules! required_int {
        ($name:ident) => {
            let mut $name = 0i32;
            if getparint(stringify!($name), &mut $name) == 0 {
                err!(concat!("must specify ", stringify!($name), "!\n"));
            }
        };
    }
    macro_rules! required_float {
        ($name:ident) => {
            let mut $name = 0.0f32;
            if getparfloat(stringify!($name), &mut $name) == 0 {
                err!(concat!("must specify ", stringify!($name), "!\n"));
            }
        };
    }
    required_int!(nxt);
    required_float!(fxt);
    required_float!(dxt);
    required_int!(nzt);
    required_float!(fzt);
    required_float!(dzt);
    required_int!(ns);
    required_float!(fs);
    required_float!(ds);
    let ext = fxt + (nxt - 1) as f32 * dxt;
    let ezt = fzt + (nzt - 1) as f32 * dzt;
    let es = fs + (ns - 1) as f32 * ds;
    required_int!(nxi);
    required_float!(fxi);
    required_float!(dxi);

    /* optional parameters */
    macro_rules! optional_int {
        ($name:ident, $default:expr) => {
            let mut $name = 0i32;
            if getparint(stringify!($name), &mut $name) == 0 {
                $name = $default;
            }
        };
    }
    macro_rules! optional_float {
        ($name:ident, $default:expr) => {
            let mut $name = 0.0f32;
            if getparfloat(stringify!($name), &mut $name) == 0 {
                $name = $default;
            }
        };
    }
    optional_int!(nxo, (nxt - 1) * 2 + 1);
    optional_float!(fxo, fxt);
    optional_float!(dxo, (dxt as f64 * 0.5) as f32);
    optional_int!(nzo, (nzt - 1) * 5 + 1);
    optional_float!(fzo, fzt);
    optional_float!(dzo, (dzt as f64 * 0.2) as f32);
    let exo = fxo + (nxo - 1) as f32 * dxo;
    let ezo = fzo + (nzo - 1) as f32 * dzo;
    if fxt > fxo || ext < exo || fzt > fzo || ezt < ezo {
        err!(" migration output range is out of traveltime table!\n");
    }
    optional_float!(v0, 1500.0);
    optional_float!(dvz, 0.0);
    optional_float!(angmax, 60.0);
    if (angmax as f64) < 0.00001 {
        err!("angmax must be positive!\n");
    }
    let mut mtmax =
        ((2.0 * dxm) as f64 * (angmax as f64 * PI / 180.0).sin() / (v0 * dt) as f64) as i32;
    if mtmax < 1 {
        mtmax = 1;
    }
    optional_float!(aperx, (0.5 * nxt as f64 * dxt as f64) as f32);
    optional_float!(offmax, 3000.0);
    optional_float!(fmax, (0.25 / dt as f64) as f32);
    optional_int!(noff, 1);
    optional_float!(off0, 0.0);
    optional_float!(doff, 99999.0);
    optional_int!(ls, 1);
    optional_int!(ntr, 100000);
    optional_int!(mtr, 100);
    checkpars();

    /* Estimation of the recording surface */
    let surfaces = decode_surfaces();
    let srf = makesurf((dxm as f64 * 0.025) as f32, &surfaces);
    let jpfp = jpfp.unwrap_or_else(|| null_file());

    fprintf!(&jpfp, "\n");
    fprintf!(&jpfp, " Migration parameters\n");
    fprintf!(&jpfp, " ================\n");
    fprintf!(&jpfp, " infile=%s \n", infile);
    fprintf!(&jpfp, " outfile=%s \n", outfile);
    fprintf!(&jpfp, " ttfile=%s \n", ttfile);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nzt=%d fzt=%g dzt=%g\n", nzt, fzt, dzt);
    fprintf!(&jpfp, " nxt=%d fxt=%g dxt=%g\n", nxt, fxt, dxt);
    fprintf!(&jpfp, " ns=%d fs=%g ds=%g\n", ns, fs, ds);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nxi=%d fxi=%g dxi=%g\n", nxi, fxi, dxi);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nzo=%d fzo=%g dzo=%g\n", nzo, fzo, dzo);
    fprintf!(&jpfp, " nxo=%d fxo=%g dxo=%g\n", nxo, fxo, dxo);
    fprintf!(&jpfp, " \n");

    let mut szif = vec![0.0; nxi as usize];
    let mut sz = vec![0.0; ns as usize];
    let mut nangl = vec![0.0; nxi as usize];
    zcoor_surfaces(
        fxi, dxi, nxi, fs, ds, ns, &srf, &mut szif, &mut sz, &mut nangl,
    );

    /* compute reference traveltime and slowness */
    let mut rmax = max(es - fxt, ext - fs);
    rmax = min(rmax, (0.5 * offmax as f64 + aperx as f64) as f32);
    let nr = 2 + (rmax / dxo) as i32;
    let mut tb = vec![table(nzt, nr); ns as usize];
    let mut pb = vec![table(nzt, nr); ns as usize];
    let mut sigb = vec![table(nzt, nr); ns as usize];
    let mut angb = vec![table(nzt, nr); ns as usize];
    for is in 0..ns as usize {
        let v00 = v0 + sz[is] * dvz;
        timeb(
            nr,
            nzt,
            dxo,
            dzt,
            fzt,
            sz[is],
            dvz,
            v00,
            &mut tb[is],
            &mut pb[is],
            &mut sigb[is],
            &mut angb[is],
        );
    }
    fprintf!(&jpfp, " nt=%d ft=%g dt=%g \n", nt, ft, dt);
    fprintf!(&jpfp, " dxm=%g fmax=%g\n", dxm, fmax);
    fprintf!(&jpfp, " noff=%d off0=%g doff=%g\n", noff, off0, doff);
    fprintf!(&jpfp, " v0=%g dvz=%g \n", v0, dvz);
    fprintf!(
        &jpfp,
        " aperx=%g offmax=%g angmax=%g\n",
        aperx,
        offmax,
        angmax
    );
    fprintf!(&jpfp, " ntr=%d mtr=%d ls=%d\n", ntr, mtr, ls);
    fprintf!(&jpfp, " ================\n");
    fflush(&jpfp);

    let mut mig = vec![table(nzo, nxo); noff as usize];
    let mut ttab = vec![table(nzt, nxt); ns as usize];
    let mut tt = table(nzt, nxt);
    let mut tsum = table(nzt, nxt);
    let mut tbs = table(nzt, nr);
    let mut tbg = table(nzt, nr);
    let mut pbs = table(nzt, nr);
    let mut pbg = table(nzt, nr);
    let mut sigbs = table(nzt, nr);
    let mut sigbg = table(nzt, nr);
    let mut angbs = table(nzt, nr);
    let mut angbg = table(nzt, nr);

    fprintf!(&jpfp, " input traveltime tables \n");
    for is in 0..ns as usize {
        let nseek = nxt as i64 * nzt as i64 * is as i64;
        efseeko(&ttfp, nseek * 4, SEEK_SET);
        fread(&mut ttab[is].data, &ttfp);
        let s = fs + is as f32 * ds;
        resit(nxt, fxt, dxt, nzt, nr, dxo, &tb[is], &mut ttab[is], s);
    }
    fprintf!(&jpfp, " start migration ... \n");
    fprintf!(&jpfp, " \n");
    fflush(&jpfp);

    let mut jtr = 1i32;
    let mut ktr = 0i32;
    let mut filter = None;
    loop {
        let sx = tr.sx as f32;
        let gx = tr.gx as f32;
        let mut io = (((gx - sx - off0) / doff) as f64 + 0.5) as i32;
        if io < 0 {
            io = 0;
        }
        if io >= noff {
            io = noff - 1;
        }
        if min(sx, gx) >= fs && max(sx, gx) <= es && max(gx - sx, sx - gx) <= offmax {
            let a_s = (sx - fs) / ds;
            let mut is = a_s as i32;
            if is == ns - 1 {
                is = ns - 2;
            }
            let mut res = a_s - is as f32;
            if (res as f64) <= 0.01 {
                res = 0.0;
            }
            if (res as f64) >= 0.99 {
                res = 1.0;
            }
            let is = is as usize;
            sum2(
                nxt,
                nzt,
                1.0 - res,
                res,
                &ttab[is],
                &ttab[is + 1],
                &mut tsum,
            );
            sum2(nr, nzt, 1.0 - res, res, &pb[is], &pb[is + 1], &mut pbs);
            sum2(nr, nzt, 1.0 - res, res, &tb[is], &tb[is + 1], &mut tbs);
            sum2(
                nr,
                nzt,
                1.0 - res,
                res,
                &sigb[is],
                &sigb[is + 1],
                &mut sigbs,
            );
            sum2(
                nr,
                nzt,
                1.0 - res,
                res,
                &angb[is],
                &angb[is + 1],
                &mut angbs,
            );

            let a_s = (gx - fs) / ds;
            let mut is = a_s as i32;
            if is == ns - 1 {
                is = ns - 2;
            }
            let mut res = a_s - is as f32;
            if (res as f64) <= 0.01 {
                res = 0.0;
            }
            if (res as f64) >= 0.99 {
                res = 1.0;
            }
            let is = is as usize;
            sum2(nxt, nzt, 1.0 - res, res, &ttab[is], &ttab[is + 1], &mut tt);
            // sum2(...,1,1,tt,tsum,tsum) permits aliasing in the original.
            for ix in 0..nxt as usize {
                for iz in 0..nzt as usize {
                    tsum[ix][iz] = tt[ix][iz] + tsum[ix][iz];
                }
            }
            sum2(nr, nzt, 1.0 - res, res, &pb[is], &pb[is + 1], &mut pbg);
            sum2(nr, nzt, 1.0 - res, res, &tb[is], &tb[is + 1], &mut tbg);
            sum2(
                nr,
                nzt,
                1.0 - res,
                res,
                &sigb[is],
                &sigb[is + 1],
                &mut sigbg,
            );
            sum2(
                nr,
                nzt,
                1.0 - res,
                res,
                &angb[is],
                &angb[is + 1],
                &mut angbg,
            );

            let a_s = (sx - fxi) / dxi;
            let mut is = a_s as i32;
            if is == nxi - 1 {
                is = nxi - 2;
            }
            let mut res = a_s - is as f32;
            if (res as f64) <= 0.01 {
                res = 0.0;
            }
            if (res as f64) >= 0.99 {
                res = 1.0;
            }
            let nangls = (1.0 - res) * nangl[is as usize] + res * nangl[is as usize + 1];
            let a_s = (gx - fxi) / dxi;
            let mut is = a_s as i32;
            if is == nxi - 1 {
                is = nxi - 2;
            }
            let mut res = a_s - is as f32;
            if (res as f64) <= 0.01 {
                res = 0.0;
            }
            if (res as f64) >= 0.99 {
                res = 1.0;
            }
            let nanglg = (1.0 - res) * nangl[is as usize] + res * nangl[is as usize + 1];

            mig2d(
                &mut tr.data,
                nt,
                ft,
                dt,
                sx,
                gx,
                &mut mig[io as usize],
                aperx,
                nxo,
                fxo,
                dxo,
                nzo as f32,
                fzo,
                dzo,
                ls,
                mtmax,
                dxm,
                fmax,
                angmax,
                &tbs,
                &tbg,
                &pbs,
                &pbg,
                &angbs,
                &angbg,
                nr,
                &tsum,
                nzt,
                fzt,
                dzt,
                nxt,
                fxt,
                dxt,
                nangls,
                nanglg,
                &sigbs,
                &sigbg,
                &mut filter,
            );
            ktr += 1;
            // ARM64 C remainder by zero returns its dividend.
            if (jtr - 1).checked_rem(mtr).unwrap_or(jtr - 1) == 0 {
                fprintf!(&jpfp, " migrated trace %d\n", jtr);
                fflush(&jpfp);
            }
        }
        jtr += 1;
        // Preserve the read before the ntr test (and the do/while first iteration).
        if fgettr(&infp, &mut tr) == 0 || jtr >= ntr {
            break;
        }
    }
    fprintf!(&jpfp, " migrated %d traces in total\n", ktr);
    let mut tro = Trace::new();
    tro.ns = nzo as u16;
    tro.d1 = dzo;
    tro.dt = (1000 * ((1000.0 * dt) as f64 + 0.5) as i32) as u16;
    tro.delrt = 0;
    tro.f1 = fzo;
    tro.f2 = fxo;
    tro.d2 = dxo;
    tro.trid = 200;
    let scal = (4.0 / PI.sqrt() * dxm as f64 / v0 as f64) as f32;
    for ixo in 0..nxo as usize {
        for io in 0..noff as usize {
            tro.data[..nzo as usize].copy_from_slice(&mig[io][ixo]);
            tro.offset = (off0 + io as f32 * doff) as i32;
            tro.tracr = 1 + ixo as i32;
            tro.tracl = 1 + io as i32 + noff * ixo as i32;
            tro.cdp = (fxo + ixo as f32 * dxo) as i32;
            tro.cdpt = 1 + io as i32;
            for izo in 0..nzo as usize {
                tro.data[izo] *= scal;
            }
            fputtr(&outfp, &mut tro);
        }
    }
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " output done\n");
    fflush(&jpfp);
    efclose(&jpfp);
    efclose(&outfp);
    exit(CWP_Exit());
}
