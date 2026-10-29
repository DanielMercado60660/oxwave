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

//! SUINVZCO3D - Seismic inversion of common offset data with V(Z) in 3D.
//! Port of su/main/migration_inversion/suinvzco3d.c (Revision 1.4).
//! Credits: Meng Xu, Center for Wave Phenomena, 07/8/95.

#![allow(unused_parens)]

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{PI, abs, nint};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUINVZCO3D - Seismic INVersion of Common Offset data with V(Z) velocity",
    "             function in 3D\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "     suinvzco3d <infile >outfile [optional parameters] \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " vfile                  file containing velocity array v[nz]\t\t",
    " nz=                    number of z samples (1st dimension) in velocity",
    " nxm=\t\t\tnumber of midpoints of input traces\t\t",
    " ny=\t\t\tnumber of lines \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt= or from header (dt) \ttime sampling interval of input data\t",
    " offs= or from header (offset) \tsource-receiver offset\t \t",
    " dxm= or from header (d2) \tsampling interval of midpoints \t\t",
    " fxm=0                  first midpoint in input trace\t\t\t",
    " nxd=5\t\t\tskipped number of midpoints (see note)\t\t",
    " dx=50.0                x sampling interval of velocity\t\t",
    " fx=0.0                 first x sample of velocity\t\t\t",
    " dz=50.0                z sampling interval of velocity\t\t",
    " nxb=nx/2\t\tband centered at midpoints (see note)\t\t",
    " fxo=0.0                x-coordinate of first output trace \t\t",
    " dxo=15.0\t\thorizontal spacing of output trace \t\t",
    " nxo=101                number of output traces \t\t\t",
    " fyo=0.0\t\ty-coordinate of first output trace\t\t",
    " dyo=15.0\t\ty-coordinate spacing of output trace\t\t",
    " nyo=101\t\tnumber of output traces in y-direction\t\t",
    " fzo=0.0                z-coordinate of first point in output trace \t",
    " dzo=15.0               vertical spacing of output trace \t\t",
    " nzo=101                number of points in output trace\t\t",
    " fmax=0.25/dt\t\tMaximum frequency set for operator antialiasing ",
    " ang=180\t\tMaximum dip angle allowed in the image\t\t",
    " verbose=1              =1 to print some useful information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " This algorithm is based on formula (50) in Geophysics Vol. 51, \t",
    " 1552-1558, by Cohen, J., Hagin, F., and Bleistein, N.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Traveltime and amplitude are calculated by ray tracing.\t\t",
    " Interpolation is used to calculate traveltime and amplitude.\t\t",
    " For each midpoint, traveltime and amplitude are calculated in the \t",
    " horizontal range of (xm-nxb*dx, xm+nxb*dx). Velocity is changed by \t",
    " linear interpolation in two upper trianglar corners whose width is \t",
    " nxc*dx and height is nzc*dz.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Eikonal equation will fail to solve if there is a polar turned ray.\t",
    " In this case, the program shows the related geometric information. \t",
    "\t\t\t\t\t\t\t\t\t",
    " Offsets are signed - may be positive or negative. \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// C alloc2/alloc3 arrays have contiguous sample storage. In particular, the
// sigma stencil at j=-1 addresses the preceding row, and iz+1 at the last
// depth addresses the next column. Preserve those addresses rather than
// giving each row a separate Vec. The tested SU runs yield zeros for
// uninitialized grid entries and tail reads. Zero these explicitly and
// provide one guard row for the otherwise undefined tail accesses.
struct Grid {
    data: Vec<f32>,
    ny: i32,
    nx: i32,
    nz: i32,
}
impl Grid {
    fn new(ny: i32, nx: i32, nz: i32) -> Self {
        Self {
            data: vec![
                0.0;
                (ny.max(0) as usize * nx.max(0) as usize * nz.max(0) as usize)
                    + nz.max(1) as usize
                    + 1
            ],
            ny,
            nx,
            nz,
        }
    }
}
impl std::ops::Index<(i32, i32, i32)> for Grid {
    type Output = f32;
    fn index(&self, (y, x, z): (i32, i32, i32)) -> &f32 {
        let row = y * self.nx + x;
        if y < 0 || y >= self.ny || row < 0 || row >= self.ny * self.nx {
            crash();
        }
        let pos = row * self.nz + z;
        self.data.get(pos as usize).unwrap_or_else(|| crash())
    }
}
impl std::ops::IndexMut<(i32, i32, i32)> for Grid {
    fn index_mut(&mut self, (y, x, z): (i32, i32, i32)) -> &mut f32 {
        &mut self.data[((y * self.nx + x) * self.nz + z) as usize]
    }
}
impl std::ops::Index<(i32, i32)> for Grid {
    type Output = f32;
    fn index(&self, (x, z): (i32, i32)) -> &f32 {
        if x < 0 || x >= self.nx {
            crash();
        }
        self.data
            .get((x * self.nz + z) as usize)
            .unwrap_or_else(|| crash())
    }
}
impl std::ops::IndexMut<(i32, i32)> for Grid {
    fn index_mut(&mut self, (x, z): (i32, i32)) -> &mut f32 {
        &mut self.data[(x * self.nz + z) as usize]
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    let temp = tr.dt as f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (((temp as f64) * 0.000001f64) as f32);
    }
    if (dt as f64) < 0.0000001 {
        err!("dt must be positive!\n");
    }
    let mut offs = 0.0f32;
    if getparfloat("offs", &mut offs) == 0 {
        offs = tr.offset as f32;
    }
    let mut dxm = 0.0f32;
    if getparfloat("dxm", &mut dxm) == 0 {
        dxm = tr.d2;
    }
    if (dxm as f64) < 0.0000001 {
        err!("dxm must be positive!\n");
    }
    let (mut ny, mut nz, mut nxm) = (0i32, 0i32, 0i32);
    if getparint("ny", &mut ny) == 0 {
        err!("must specify ny!\n");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!\n");
    }
    if getparint("nxm", &mut nxm) == 0 {
        err!("must specify nxm!\n");
    }
    let mut vfile = String::new();
    if getparstring("vfile", &mut vfile) == 0 {
        err!("must specify vfile!\n");
    }
    macro_rules! fp {
        ($name:ident,$default:expr) => {
            let mut $name = 0.0f32;
            if getparfloat(stringify!($name), &mut $name) == 0 {
                $name = $default;
            }
        };
    }
    macro_rules! ip {
        ($name:ident,$default:expr) => {
            let mut $name = 0i32;
            if getparint(stringify!($name), &mut $name) == 0 {
                $name = $default;
            }
        };
    }
    fp!(fxm, 0.0);
    let mut dx = 0.0f32;
    if getparfloat("dxm", &mut dx) == 0 {
        dx = 50.0;
    }
    fp!(dy, 50.0);
    fp!(dz, 50.0);
    ip!(nxb, nxm / 2);
    ip!(nxc, 0);
    ip!(na, 100);
    nxc = nxc.min(nxb);
    ip!(nzc, 0);
    nzc = nzc.min(nz);
    fp!(fxo, 0.0);
    if fxo < fxm {
        err!("the first output trace is out of velocity model.\n");
    }
    fp!(dxo, 15.0);
    ip!(nxo, 101);
    fp!(fyo, 0.0);
    fp!(dyo, 15.0);
    ip!(nyo, 101);
    fp!(fzo, 0.0);
    fp!(dzo, 15.0);
    ip!(nzo, 101);
    ip!(verbose, 1);
    fp!(fmax, ((0.25f64 / (dt as f64)) as f32));
    fp!(ang, 180.0);
    // The C parses these, but its aliasing/dip filters are commented out.
    let _ = (nxc, nzc, fmax, ang);
    checkpars();
    let ey = (((ny - 1i32) as f32) * dy);
    let ez = (((nz - 1i32) as f32) * dz);
    let odt = ((1.0f64 / (dt as f64)) as f32);
    let mut outtrace = Grid::new(nyo, nxo, nzo);
    let mut ts = Grid::new(nyo, nxo, nzo);
    let mut amp = Grid::new(nyo, nxo, nzo);
    let mut p3 = Grid::new(nyo, nxo, nzo);
    let mut v = vec![0.0f32; nz.max(0) as usize + 1];
    let mut vnew = vec![0.0f32; nzo.max(0) as usize + 1];
    // fopen is deliberately unchecked by the C program.
    let Some(vfp) = fopen(&vfile, "r") else {
        crash();
    };
    if fread(&mut v[..nz.max(0) as usize], &vfp) != nz as usize {
        err!("cannot read %d velocities from file %s", nz, vfile);
    }
    for i in 0..nzo {
        let mut temp = (((i as f32) * dzo) / dz);
        let mut nsamp = temp as i32;
        temp = (temp - (nsamp as f32));
        if nsamp > nz - 1 {
            nsamp = nz - 1;
        }
        vnew[i as usize] = ((((1i32 as f32) - temp) * v[(nsamp as usize)])
            + (temp * v[((nsamp + 1i32) as usize)]));
    }
    if abs(offs) > (((2i32 * nxb) as f32) * dx) {
        err!("\t band NXB is too small!");
    }
    if ((nxm as f32) * dxm) < ((nxo as f32) * dxo) || ((ny as f32) * dy) < ((nyo as f32) * dyo) {
        err!("\t not enough range in input  \t\tfor the imaging");
    }
    amplitude(
        nxo, nyo, nzo, dxo, dyo, dzo, na, &vnew, offs, &mut ts, &mut amp, &mut p3,
    );
    fprintf!(stderr(), "\nTables generated...\n");
    for iym in 0..ny {
        for ixm in 0..nxm {
            let xm = (fxm + ((ixm as f32) * dxm));
            let xs = (((xm as f64) - (0.5f64 * (offs as f64))) as f32);
            let xg = (xs + offs);
            let ys = ((iym as f32) * dy);
            let yg = ((iym as f32) * dy);
            let mut yo = fyo;
            for iyo in 0..nyo {
                let mut xo = fxo;
                for ixo in 0..nxo {
                    let i = ((((((xo - xm) * (xo - xm)) + ((yo - ys) * (yo - ys))) as f64).sqrt()
                        / (dx as f64)) as i32);
                    if i <= nxb {
                        let mut zo = fzo + dzo;
                        for izo in 1..nzo {
                            if !(zo - ez > 0.0
                                || yo - ey > 0.0
                                || abs(yo - ys) / dyo > (nyo - 2) as f32
                                || abs(xo - xs) / dxo > (nxo - 2) as f32)
                            {
                                let xi = abs(((xo - xs) / dxo));
                                let ix = xi as i32;
                                let yi = (abs((yo - ys)) / dyo);
                                let iy = yi as i32;
                                let zi = (zo / dzo);
                                let iz = zi as i32;
                                let sx = (xi - (ix as f32));
                                let sy = (yi - (iy as f32));
                                let sz = (zi - (iz as f32));
                                let mut tsd = (((1.0f64 - (sy as f64))
                                    * (((1.0f64 - (sz as f64))
                                        * (((1.0f64 - (sx as f64)) * (ts[(iy, ix, iz)] as f64))
                                            + ((sx * ts[(iy, (ix + 1i32), iz)]) as f64)))
                                        + ((sz as f64)
                                            * (((1.0f64 - (sx as f64))
                                                * (ts[(iy, ix, (iz + 1i32))] as f64))
                                                + ((sx * ts[(iy, (ix + 1i32), (iz + 1i32))])
                                                    as f64)))))
                                    as f32);
                                tsd = (((tsd as f64)
                                    + ((sy as f64)
                                        * (((1.0f64 - (sz as f64))
                                            * (((1.0f64 - (sx as f64))
                                                * (ts[((iy + 1i32), ix, iz)] as f64))
                                                + ((sx * ts[((iy + 1i32), (ix + 1i32), iz)])
                                                    as f64)))
                                            + ((sz as f64)
                                                * (((1.0f64 - (sx as f64))
                                                    * (ts[((iy + 1i32), ix, (iz + 1i32))]
                                                        as f64))
                                                    + ((sx
                                                        * ts[(
                                                            (iy + 1i32),
                                                            (ix + 1i32),
                                                            (iz + 1i32),
                                                        )])
                                                        as f64))))))
                                    as f32);
                                // C uses the output-grid index in the original
                                // velocity array. Finer output spacing can put
                                // iz beyond that allocation even inside the
                                // model. Use zero for these undefined tail
                                // reads; sins/sing only feed a disabled filter
                                // and cannot affect the output.
                                let velocity = v.get(iz as usize).copied().unwrap_or(0.0);
                                let mut sins = ((p3[(iy, ix, iz)] * p3[(iy, ix, iz)])
                                    * (velocity * velocity));
                                if sins >= 1.0 {
                                    sins = ((1.0f64 / (sins as f64)) as f32);
                                    sins = ((1.0f64 - (sins as f64)).sqrt() as f32);
                                }
                                let xi = (abs((xo - xg)) / dxo);
                                let ix = xi as i32;
                                let yi = (abs((yo - yg)) / dyo);
                                let iy = yi as i32;
                                let zi = (zo / dzo);
                                let iz = zi as i32;
                                let sx = (xi - (ix as f32));
                                let sy = (yi - (iy as f32));
                                let sz = (zi - (iz as f32));
                                let tgd = ((((1.0f64 - (sy as f64))
                                    * (((1.0f64 - (sz as f64))
                                        * (((1.0f64 - (sx as f64)) * (ts[(iy, ix, iz)] as f64))
                                            + ((sx * ts[(iy, (ix + 1i32), iz)]) as f64)))
                                        + ((sz as f64)
                                            * (((1.0f64 - (sx as f64))
                                                * (ts[(iy, ix, (iz + 1i32))] as f64))
                                                + ((sx * ts[(iy, (ix + 1i32), (iz + 1i32))])
                                                    as f64)))))
                                    + ((sy as f64)
                                        * (((1.0f64 - (sz as f64))
                                            * (((1.0f64 - (sx as f64))
                                                * (ts[((iy + 1i32), ix, iz)] as f64))
                                                + ((sx * ts[((iy + 1i32), (ix + 1i32), iz)])
                                                    as f64)))
                                            + ((sz as f64)
                                                * (((1.0f64 - (sx as f64))
                                                    * (ts[((iy + 1i32), ix, (iz + 1i32))]
                                                        as f64))
                                                    + ((sx
                                                        * ts[(
                                                            (iy + 1i32),
                                                            (ix + 1i32),
                                                            (iz + 1i32),
                                                        )])
                                                        as f64))))))
                                    as f32);
                                let velocity = v.get(iz as usize).copied().unwrap_or(0.0);
                                let mut sing = ((p3[(iy, ix, iz)] * p3[(iy, ix, iz)])
                                    * (velocity * velocity));
                                if sing >= 1.0 {
                                    sing = ((1.0f64 / (sing as f64)) as f32);
                                    sing = ((1.0f64 - (sing as f64)).sqrt() as f32);
                                }
                                let _ = (sins, sing);
                                let ampd = if xo <= xs {
                                    amp[(iy, ix, iz)]
                                } else {
                                    let ix = ((abs((xo - xs)) / dxo) as i32);
                                    amp[(iy, ix, iz)]
                                };
                                let mut samp = ((tsd + tgd) * odt);
                                let nsamp = samp as i32;
                                if nsamp <= nt - 2 {
                                    samp = (samp - (nsamp as f32));
                                    outtrace[(iyo, ixo, izo)] = (((outtrace[(iyo, ixo, izo)]
                                        as f64)
                                        + ((ampd as f64)
                                            * (((samp * tr.data[((nsamp + 1i32) as usize)])
                                                as f64)
                                                + ((1.0f64 - (samp as f64))
                                                    * (tr.data[(nsamp as usize)] as f64)))))
                                        as f32);
                                }
                            }
                            zo += dzo;
                        }
                    }
                    xo += dxo;
                }
                yo += dyo;
            }
            if ixm < nxm - 1 || iym < ny - 1 {
                gettr(&mut tr);
            }
        }
        if verbose != 0 {
            fprintf!(stderr(), "\tfinish line%d\n", iym + 1);
        }
    }
    let temp = ((((4i32 as f64) / (PI).sqrt()) * (dxm as f64)) as f32);
    for iyo in 0..nyo {
        for ixo in 0..nxo {
            for izo in 0..nzo {
                outtrace[(iyo, ixo, izo)] *= temp;
            }
            tro.offset = offs as i32;
            tro.tracl = 1 + ixo * nyo + iyo;
            tro.tracr = tro.tracl;
            tro.ns = nzo as u16;
            tro.d1 = dzo;
            tro.f1 = fzo;
            tro.f2 = fxo;
            tro.d2 = dxo;
            tro.cdp = ixo;
            tro.trid = 200;
            for izo in 0..nzo {
                tro.data[izo as usize] = outtrace[(iyo, ixo, izo)];
            }
            puttr(&mut tro);
        }
    }
    exit(CWP_Exit());
}

// The original dereferences invalid pointers on these paths. Preserve the
// shell-visible status (128 + SIGSEGV) without calling C or introducing Rust
// undefined behavior. Core has no API to terminate by signal; a parent that
// inspects wait status directly can distinguish this exit from SIGSEGV.
fn crash() -> ! {
    exit(139)
}

// Calculation of a 3x3 determinant.
fn determ(a: [[f32; 3]; 3]) -> f32 {
    let mut det = ((a[0][0] * a[1][1]) * a[2][2]);
    det = (det + ((a[1][0] * a[2][1]) * a[0][2]));
    det = (det + ((a[2][0] * a[1][2]) * a[0][1]));
    det = (det - ((a[0][2] * a[1][1]) * a[2][0]));
    det = (det - ((a[1][2] * a[2][1]) * a[0][0]));
    det = (det - ((a[2][2] * a[1][0]) * a[0][1]));
    det
}

// Ray equation solution in a V(Z) medium in polar coordinates.
#[allow(non_snake_case)]
fn SingleRay(
    nz: i32,
    a: f32,
    v: &[f32],
    dz: f32,
    x: &mut Grid,
    p3: &mut Grid,
    sigma: &mut Grid,
    t: &mut Grid,
    dzda2: &mut Grid,
    ray: i32,
) {
    let (mut sig, mut tt, mut temp) = (0.0f32, 0.0f32, 0.0f32);
    for i in 0..nz {
        sigma[(ray, i)] = 0.0;
        p3[(ray, i)] = 0.0;
        x[(ray, i)] = 0.0;
        t[(ray, i)] = 0.0;
        dzda2[(ray, i)] = 0.0;
    }
    for i in 0..nz {
        if (((1.0f64 / (v[(i as usize)] as f64)) / (v[(i as usize)] as f64))
            - ((((a as f64).sin() * (a as f64).sin()) / (v[0] as f64)) / (v[0] as f64)))
            > 0.0
        {
            p3[(ray, i)] = ((abs((((1.0f64 / (v[(i as usize)] as f64))
                / (v[(i as usize)] as f64))
                - ((((a as f64).sin() * (a as f64).sin()) / (v[0] as f64)) / (v[0] as f64)))))
            .sqrt() as f32);
            sigma[(ray, i)] = sig;
            sig = (sig + (dz / p3[(ray, i)]));
            x[(ray, i)] = ((((sigma[(ray, i)] as f64) * (a as f64).sin()) / (v[0] as f64)) as f32);
            t[(ray, i)] = tt;
            tt = (tt + (((dz / v[(i as usize)]) / v[(i as usize)]) / p3[(ray, i)]));
            dzda2[(ray, i)] = temp;
            temp = (((temp as f64)
                - ((((((dz / p3[(ray, i)]) as f64) * (a as f64).sin()) / (v[0] as f64))
                    * (p3[(ray, 0i32)] as f64))
                    / (p3[(ray, i)] as f64))) as f32);
            if (p3[(ray, i)] as f64) <= 0.000001 {
                fprintf!(stderr(), "\nTurning point: z=%f\n", ((i as f32) * dz));
                continue;
            }
        } else {
            fprintf!(
                stderr(),
                "\nTurning point: ang=%f,z=%f",
                a,
                ((i as f32) * dz)
            );
            for ii in i..nz {
                sigma[(ray, ii)] = 0.0;
                x[(ray, ii)] = -1.0;
                t[(ray, ii)] = 0.0;
                dzda2[(ray, ii)] = 0.0;
            }
            break;
        }
    }
}

#[allow(clippy::too_many_arguments)]
// Interpolate from ray coordinates to the 2D (x,z) grid.
fn para2d(
    na: i32,
    nx: i32,
    nz: i32,
    dx: f32,
    dz: f32,
    v: &[f32],
    t2d: &mut Grid,
    sigma2d: &mut Grid,
    p2d: &mut Grid,
    dzda2d: &mut Grid,
    a2d: &mut Grid,
) {
    let da = (((3.14159265f64 / (na as f64)) / 2.0f64) as f32);
    let mut p3 = Grid::new(1, na, nz);
    let mut sigma = Grid::new(1, na, nz);
    let mut t = Grid::new(1, na, nz);
    let mut x = Grid::new(1, na, nz);
    let mut dzda2 = Grid::new(1, na, nz);
    let mut a2 = Grid::new(1, na, nz);
    let mut a = 0.0f32;
    for i in 0..na {
        SingleRay(
            nz, a, v, dz, &mut x, &mut p3, &mut sigma, &mut t, &mut dzda2, i,
        );
        for k in 0..nz {
            a2[(i, k)] = a;
        }
        a += da;
    }
    for k in 1..nz {
        for i in 0..nx {
            let x0 = (dx * (i as f32));
            let mut iflag = 0;
            for ii in 0..na {
                if x[(ii, k)] > x0 && x[(0, k)] <= x0 {
                    iflag = ii;
                    break;
                }
            }
            let sx = if iflag == 0 {
                iflag = na - 1;
                0.0f32
            } else {
                ((x[(iflag, k)] - x0) / (x[(iflag, k)] - x[((iflag - 1i32), k)]))
            };
            t2d[(i, k)] = ((((1.0f64 - (sx as f64)) * (t[(iflag, k)] as f64))
                + ((sx * t[((iflag - 1i32), k)]) as f64)) as f32);
            sigma2d[(i, k)] = ((((1.0f64 - (sx as f64)) * (sigma[(iflag, k)] as f64))
                + ((sx * sigma[((iflag - 1i32), k)]) as f64))
                as f32);
            p2d[(i, k)] = ((((1.0f64 - (sx as f64)) * (p3[(iflag, k)] as f64))
                + ((sx * p3[((iflag - 1i32), k)]) as f64)) as f32);
            a2d[(i, k)] = ((((1.0f64 - (sx as f64)) * (a2[(iflag, k)] as f64))
                + ((sx * a2[((iflag - 1i32), k)]) as f64)) as f32);
            dzda2d[(i, k)] = ((((1.0f64 - (sx as f64)) * (dzda2[(iflag, k)] as f64))
                + ((sx * dzda2[((iflag - 1i32), k)]) as f64)) as f32);
            if t2d[(i, k)] == 0.0 {
                fprintf!(stderr(), "z=%f\n", ((k as f32) * dz));
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
// Interpolate ray parameters from 2D to 3D.
fn para3d(
    nxo: i32,
    nyo: i32,
    nzo: i32,
    dxo: f32,
    dyo: f32,
    t3d: &mut Grid,
    sigma3d: &mut Grid,
    p3d: &mut Grid,
    dzda23d: &mut Grid,
    a23d: &mut Grid,
    t2d: &Grid,
    sigma2d: &Grid,
    p2d: &Grid,
    dzda2d: &Grid,
    a2d: &Grid,
) {
    for iy in 0..nyo {
        for ix in 0..nxo {
            for iz in 0..nzo {
                let mut temp1 = (((((((ix as f32) * dxo) * ((ix as f32) * dxo))
                    + (((iy as f32) * dyo) * ((iy as f32) * dyo)))
                    as f64)
                    .sqrt()
                    / (dxo as f64)) as f32);
                let ix2 = temp1 as i32;
                temp1 = (temp1 - (ix2 as f32));
                t3d[(iy, ix, iz)] = ((((1.0f64 - (temp1 as f64)) * (t2d[(ix2, iz)] as f64))
                    + ((temp1 * t2d[((ix2 + 1i32), iz)]) as f64))
                    as f32);
                p3d[(iy, ix, iz)] = ((((1.0f64 - (temp1 as f64)) * (p2d[(ix2, iz)] as f64))
                    + ((temp1 * p2d[((ix2 + 1i32), iz)]) as f64))
                    as f32);
                sigma3d[(iy, ix, iz)] = ((((1.0f64 - (temp1 as f64)) * (sigma2d[(ix2, iz)] as f64))
                    + ((temp1 * sigma2d[((ix2 + 1i32), iz)]) as f64))
                    as f32);
                a23d[(iy, ix, iz)] = ((((1.0f64 - (temp1 as f64)) * (a2d[(ix2, iz)] as f64))
                    + ((temp1 * a2d[((ix2 + 1i32), iz)]) as f64))
                    as f32);
                dzda23d[(iy, ix, iz)] = ((((1.0f64 - (temp1 as f64)) * (dzda2d[(ix2, iz)] as f64))
                    + ((temp1 * dzda2d[((ix2 + 1i32), iz)]) as f64))
                    as f32);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
// Calculate the Jacobian and Beylkin determinant for the amplitude weight.
fn amplitude(
    nxo: i32,
    nyo: i32,
    nzo: i32,
    dxo: f32,
    dyo: f32,
    dzo: f32,
    na: i32,
    v: &[f32],
    off: f32,
    t3d: &mut Grid,
    amp: &mut Grid,
    p3d: &mut Grid,
) {
    let nxr = nint(
        (((((((nxo as f32) * dxo) * (nxo as f32)) * dxo)
            + ((((nyo as f32) * dyo) * (nyo as f32)) * dyo)) as f64)
            .sqrt()
            / (dxo as f64)),
    ) + 1;
    let mut p2d = Grid::new(1, nxr, nzo);
    let mut t2d = Grid::new(1, nxr, nzo);
    let mut sigma2d = Grid::new(1, nxr, nzo);
    let mut dzda2d = Grid::new(1, nxr, nzo);
    let mut a2d = Grid::new(1, nxr, nzo);
    let mut jacob = Grid::new(nyo, nxo, nzo);
    let mut asag = Grid::new(nyo, nxo, nzo);
    let mut h = Grid::new(nyo, nxo, nzo);
    let mut sigma3d = Grid::new(nyo, nxo, nzo);
    let mut dzda23d = Grid::new(nyo, nxo, nzo);
    let mut a23d = Grid::new(nyo, nxo, nzo);
    para2d(
        na,
        nxr,
        nzo,
        dxo,
        dzo,
        v,
        &mut t2d,
        &mut sigma2d,
        &mut p2d,
        &mut dzda2d,
        &mut a2d,
    );
    para3d(
        nxo,
        nyo,
        nzo,
        dxo,
        dyo,
        t3d,
        &mut sigma3d,
        p3d,
        &mut dzda23d,
        &mut a23d,
        &t2d,
        &sigma2d,
        &p2d,
        &dzda2d,
        &a2d,
    );
    let xi = (off / dxo);
    let ix = xi as i32;
    let sx = (xi - (ix as f32));
    let _ = sx;
    let mut det = [[0.0f32; 3]; 3];
    for i in 0..nyo {
        for j in 0..nxo {
            for k in 1..nzo {
                let a1 = if j == 0 {
                    ((3.14159265f64 / 2.0f64) as f32)
                } else {
                    ((((((i as f32) * dyo) / (j as f32)) * dxo) as f64).atan() as f32)
                };
                let a2 = a23d[(i, j, k)];
                det[0][0] = ((((a2 as f64).sin() * (a1 as f64).cos()) / (v[0] as f64)) as f32);
                det[0][1] = ((((a1 as f64).sin() * (a2 as f64).sin()) / (v[0] as f64)) as f32);
                det[0][2] = p3d[(i, j, k)];
                det[1][0] = ((((((-sigma3d[(i, j, k)]) as f64) * (a1 as f64).sin())
                    * (a2 as f64).sin())
                    / (v[0] as f64)) as f32);
                det[1][1] = (((((sigma3d[(i, j, k)] as f64) * (a1 as f64).cos())
                    * (a2 as f64).sin())
                    / (v[0] as f64)) as f32);
                det[1][2] = 0.0;
                det[2][0] = (((((sigma3d[(i, j, k)] as f64) * (a1 as f64).cos())
                    * (a2 as f64).cos())
                    / (v[0] as f64)) as f32);
                det[2][1] = (((((sigma3d[(i, j, k)] as f64) * (a1 as f64).sin())
                    * (a2 as f64).cos())
                    / (v[0] as f64)) as f32);
                det[2][2] = dzda23d[(i, j, k)];
                let temp = (((4.0f64 * PI) * ((v[0] * abs(determ(det))) as f64).sqrt()) as f32);
                if temp > 0.0 {
                    jacob[(i, j, k)] =
                        (((100.0f64 * ((a2 as f64).sin()).sqrt()) / (temp as f64)) as f32);
                }
            }
        }
    }
    for i in 0..nyo {
        for j in 0..nxo {
            for k in 0..nzo {
                asag[(i, j, k)] = (jacob[(i, j, k)] * jacob[(i, abs((j - ix)), k)]);
            }
        }
    }
    for i in 1..nyo {
        for j in 0..nxo {
            for k in 1..nzo {
                if t3d[(i, j, k)] == 0.0 {
                    amp[(i, j, k)] = 0.0;
                } else {
                    let a1s = if j == 0 {
                        ((PI / 2.0f64) as f32)
                    } else {
                        ((((((i as f32) * dyo) / (j as f32)) / dxo) as f64).atan() as f32)
                    };
                    let a1g = if abs((((j as f32) * dxo) - off)) == 0.0 {
                        ((PI / (2i32 as f64)) as f32)
                    } else {
                        (((((i as f32) * dyo) / abs((((j as f32) * dxo) - off))) as f64).atan()
                            as f32)
                    };
                    let a2s = a23d[(i, j, k)];
                    let a2g = a23d[(i, abs(j - ix), k)];
                    let (mut jl, mut jr) = if j == nxo - 1 {
                        (j - 1, j)
                    } else {
                        (j - 1, j + 1)
                    };
                    let mut temp1 = (((((sigma3d[(i, jl, k)] - sigma3d[(i, jr, k)]) / dxo) as f64)
                        / 2.0f64) as f32);
                    if abs(j - ix) == 0 {
                        jl = abs(j - ix);
                        jr = jl + 1;
                    } else {
                        jl = abs(j - ix) - 1;
                        jr = jl + 1;
                    }
                    let mut temp2 = ((sigma3d[(i, jl, k)] - sigma3d[(i, jr, k)]) / dxo);
                    let mut ps = [0.0f32; 3];
                    let mut pg = [0.0f32; 3];
                    ps[0] = ((((1.0f64 / (v[0] as f64)) * (a1s as f64).cos()) * (a2s as f64).sin())
                        as f32);
                    ps[1] = ((((1.0f64 / (v[0] as f64)) * (a1s as f64).sin()) * (a2s as f64).sin())
                        as f32);
                    ps[2] = p3d[(i, j, k)];
                    pg[0] = ((((1.0f64 / (v[0] as f64)) * (a1g as f64).cos()) * (a2g as f64).sin())
                        as f32);
                    pg[1] = ((((1.0f64 / (v[0] as f64)) * (a1g as f64).sin()) * (a2g as f64).sin())
                        as f32);
                    pg[2] = p3d[(i, abs(j - ix), k)];
                    det[0][0] = (ps[0] + pg[0]);
                    det[0][1] = (ps[1] + pg[1]);
                    det[0][2] = (ps[2] + pg[2]);
                    det[1][0] = ((((-(1.0f64 + ((ps[0] * temp1) as f64)))
                        / (sigma3d[(i, j, k)] as f64))
                        - ((1.0f64 + ((pg[0] * temp2) as f64))
                            / (sigma3d[(i, abs((j - ix)), k)] as f64)))
                        as f32);
                    det[1][1] = (((-(ps[1] * temp1)) / sigma3d[(i, j, k)])
                        - ((pg[1] * temp2) / sigma3d[(i, abs((j - ix)), k)]));
                    det[1][2] = (((((((ps[0] as f64) * (1.0f64 + ((ps[0] * temp1) as f64)))
                        / (sigma3d[(i, j, k)] as f64))
                        + (((ps[1] * (ps[1] * temp1)) / sigma3d[(i, j, k)]) as f64))
                        / (ps[2] as f64))
                        + (((((pg[0] as f64) * (1.0f64 + ((pg[0] * temp2) as f64)))
                            / (sigma3d[(i, abs((j - ix)), k)] as f64))
                            + (((pg[1] * (pg[1] * temp2)) / sigma3d[(i, abs((j - ix)), k)])
                                as f64))
                            / (pg[2] as f64))) as f32);
                    let (il, ir) = if i == nyo - 1 {
                        (i - 1, i)
                    } else {
                        (i - 1, i + 1)
                    };
                    temp1 = (((((sigma3d[(il, j, k)] - sigma3d[(ir, j, k)]) / dyo) as f64) / 2.0f64)
                        as f32);
                    temp2 = (((((sigma3d[(il, abs((j - ix)), k)] - sigma3d[(ir, abs((j - ix)), k)])
                        as f64)
                        / 2.0f64)
                        / (dyo as f64)) as f32);
                    det[2][0] = (((-(ps[0] * temp1)) / sigma3d[(i, j, k)])
                        - ((pg[0] * temp2) / sigma3d[(i, abs((j - ix)), k)]));
                    det[2][1] = ((((-(1.0f64 + ((ps[1] * temp1) as f64)))
                        / (sigma3d[(i, j, k)] as f64))
                        - ((1.0f64 + ((pg[1] * temp2) as f64))
                            / (sigma3d[(i, abs((j - ix)), k)] as f64)))
                        as f32);
                    det[2][2] = (((((((ps[0] * (ps[0] * temp1)) / sigma3d[(i, j, k)]) as f64)
                        + (((ps[1] as f64) * (1.0f64 + ((ps[1] * temp1) as f64)))
                            / (sigma3d[(i, j, k)] as f64)))
                        / (ps[2] as f64))
                        + (((((pg[0] * (pg[0] * temp2)) / sigma3d[(i, abs((j - ix)), k)]) as f64)
                            + (((pg[1] as f64) * (1.0f64 + ((pg[1] * temp2) as f64)))
                                / (sigma3d[(i, abs((j - ix)), k)] as f64)))
                            / (pg[2] as f64))) as f32);
                    h[(i, j, k)] = determ(det);
                    let r = ((((((ps[0] + pg[0]) * (ps[0] + pg[0]))
                        + ((ps[1] + pg[1]) * (ps[1] + pg[1])))
                        + ((ps[2] + pg[2]) * (ps[2] + pg[2]))) as f64)
                        .sqrt() as f32);
                    amp[(i, j, k)] = ((h[(i, j, k)] / asag[(i, j, k)]) / r);
                }
            }
        }
    }
    for i in 1..nxo {
        for k in 1..nzo {
            amp[(0, i, k)] = amp[(1, i, k)];
        }
    }
}
