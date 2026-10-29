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

//! SUFDMOD2_PML - second-order acoustic modeling with PML boundaries.
//! Port of su/main/synthetics_waveforms_testpatterns/sufdmod2_pml.c, Revision 1.9.
//! Credits: Dave Hale; John Stockwell, 1993; Craig Artley, 1993;
//! Michael Holzrichter, TAMU, 1998 (PML); Greg Wimpey,
//! CWP/WesternGeco, 2006 (density-variation corrections).

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::cwp::{abs, max, min, nint};
use oxwave_f_cubicspline::cmonot;
use std::ops::{Index, IndexMut};
mod kernels;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUFDMOD2_PML - Finite-Difference MODeling (2nd order) for acoustic wave",
    "    equation with PML absorbing boundary conditions.\t\t\t",
    " Caveat: experimental PML absorbing boundary condition version,\t",
    "may be buggy!\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sufdmod2_pml <vfile >wfile nx= nz= tmax= xs= zs= [optional parameters]",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " <vfile\t\tfile containing velocity[nx][nz]\t\t",
    " >wfile\t\tfile containing waves[nx][nz] for time steps\t",
    " nx=\t\t\tnumber of x samples (2nd dimension)\t\t",
    " nz=\t\t\tnumber of z samples (1st dimension)\t\t",
    " xs=\t\t\tx coordinates of source\t\t\t\t",
    " zs=\t\t\tz coordinates of source\t\t\t\t",
    " tmax=\t\t\tmaximum time\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nt=1+tmax/dt\t\tnumber of time samples (dt determined for stability)",
    " mt=1\t\t\tnumber of time steps (dt) per output time step\t",
    " \t\t\t\t\t\t\t\t\t",
    " dx=1.0\t\tx sampling interval\t\t\t\t",
    " fx=0.0\t\tfirst x sample\t\t\t\t\t",
    " dz=1.0\t\tz sampling interval\t\t\t\t",
    " fz=0.0\t\tfirst z sample\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " fmax = vmin/(10.0*h)\tmaximum frequency in source wavelet\t\t",
    " fpeak=0.5*fmax\tpeak frequency in ricker wavelet\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dfile=\t\tinput file containing density[nx][nz]\t\t",
    " vsx=\t\t\tx coordinate of vertical line of seismograms\t",
    " hsz=\t\t\tz coordinate of horizontal line of seismograms\t",
    " vsfile=\t\toutput file for vertical line of seismograms[nz][nt]",
    " hsfile=\t\toutput file for horizontal line of seismograms[nx][nt]",
    " ssfile=\t\toutput file for source point seismograms[nt]\t",
    " verbose=0\t\t=1 for diagnostic messages, =2 for more\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " abs=1,1,1,1\t\tAbsorbing boundary conditions on top,left,bottom,right",
    " \t\t\tsides of the model. \t\t\t\t",
    " \t\t=0,1,1,1 for free surface condition on the top\t\t",
    "                                                                       ",
    " ...PML parameters....                                                 ",
    " pml_max=1000.0        PML absorption parameter                        ",
    " pml_thick=0           half-thickness of pml layer (0 = do not use PML)",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This program uses the traditional explicit second order differencing\t",
    " method. \t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Two different absorbing boundary condition schemes are available. The ",
    " first is a traditional absorbing boundary condition scheme created by ",
    " Hale, 1990. The second is based on the perfectly matched layer (PML)\t",
    " method of Berenger, 1995.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

/// Contiguous alloc2float storage; the first dimension is the fast axis.
struct Grid {
    n1: i32,
    data: Vec<f32>,
}
impl Grid {
    fn new(n1: i32, n2: i32) -> Self {
        Self {
            n1,
            data: vec![0.0; (n1 as usize).wrapping_mul(n2 as usize)],
        }
    }
}
impl Index<(i32, i32)> for Grid {
    type Output = f32;
    fn index(&self, (i, j): (i32, i32)) -> &f32 {
        &self.data[(i * self.n1 + j) as usize]
    }
}
impl IndexMut<(i32, i32)> for Grid {
    fn index_mut(&mut self, (i, j): (i32, i32)) -> &mut f32 {
        &mut self.data[(i * self.n1 + j) as usize]
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut nx, mut nz, mut tmax) = (0, 0, 0.0f32);
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx!");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!");
    }
    if getparfloat("tmax", &mut tmax) == 0 {
        err!("must specify tmax!");
    }
    let nxs = countparval("xs");
    let nzs = countparval("zs");
    if nxs != nzs {
        err!("number of xs = %d must equal number of zs = %d", nxs, nzs);
    }
    let ns = nxs;
    if ns == 0 {
        err!("must specify xs and zs!");
    }
    let mut xs = vec![0.0f32; ns as usize];
    let mut zs = vec![0.0f32; ns as usize];
    let mut ixs = vec![0; ns as usize];
    let mut izs = vec![0; ns as usize];
    getparfloat("xs", &mut xs);
    getparfloat("zs", &mut zs);
    let mut absorb_flags = [1i32; 4];
    let nabs = countparval("abs");
    if nabs == 4 {
        getparint("abs", &mut absorb_flags[..]);
    } else if nabs != 0 {
        warn!("Number of abs %d, using abs=1,1,1,1", nabs);
    }
    let mut nt = 0;
    let mut mt = 1;
    let mut dx = 1.0f32;
    let mut fx = 0.0f32;
    let mut dz = 1.0f32;
    let mut fz = 0.0f32;
    getparint("nt", &mut nt);
    getparint("mt", &mut mt);
    getparfloat("dx", &mut dx);
    getparfloat("fx", &mut fx);
    getparfloat("dz", &mut dz);
    getparfloat("fz", &mut fz);
    let mut pml_max = 1000.0f32;
    let mut pml_thick = 0i32;
    getparfloat("pml_max", &mut pml_max);
    getparint("pml_thick", &mut pml_thick);
    let pml_thickness = 2 * pml_thick;
    for is in 0..ns as usize {
        ixs[is] = nint(((xs[is] - fx) / dx) as f64);
        // The original uses dx for this z index as well.
        izs[is] = nint(((zs[is] - fz) / dx) as f64);
    }
    let mut hsz = 0.0f32;
    getparfloat("hsz", &mut hsz);
    let hs1 = nint(((hsz - fz) / dz) as f64);
    let mut vsx = 0.0f32;
    getparfloat("vsx", &mut vsx);
    let vs2 = nint(((vsx - fx) / dx) as f64);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut dfile = String::new();
    let mut hsfile = String::new();
    let mut vsfile = String::new();
    let mut ssfile = String::new();
    getparstring("dfile", &mut dfile);
    getparstring("hsfile", &mut hsfile);
    getparstring("vsfile", &mut vsfile);
    getparstring("ssfile", &mut ssfile);
    // There is no checkpars() call in the C program.
    let mut s = Grid::new(nz, nx);
    let mut dvv = Grid::new(nz, nx);
    let mut od = Grid::new(nz, nx);
    let mut pm = Grid::new(nz, nx);
    let mut p = Grid::new(nz, nx);
    let mut pp = Grid::new(nz, nx);
    // The original velocity fread is unchecked.
    fread(&mut dvv.data, &stdin());
    let mut vmin = dvv[(0, 0)];
    let mut vmax = dvv[(0, 0)];
    for ix in 0..nx {
        for iz in 0..nz {
            vmin = min(vmin, dvv[(ix, iz)]);
            vmax = max(vmax, dvv[(ix, iz)]);
        }
    }
    let h = min(abs(dx), abs(dz));
    let dt = (h as f64 / (2.0 * vmax as f64)) as f32;
    let mut fmax = (vmin as f64 / (10.0 * h as f64)) as f32;
    getparfloat("fmax", &mut fmax);
    let mut fpeak = (0.5 * fmax as f64) as f32;
    getparfloat("fpeak", &mut fpeak);
    if nt == 0 {
        nt = (1.0f32 + tmax / dt) as i32;
    }

    // Open outputs in the same order as C, before reading the density file.
    let hseisfp = if !hsfile.is_empty() {
        Some(fopen(&hsfile, "w").unwrap_or_else(|| err!("cannot open hsfile=%s", hsfile)))
    } else {
        None
    };
    let mut hs = if hseisfp.is_some() {
        Some(Grid::new(nt, nx))
    } else {
        None
    };
    let vseisfp = if !vsfile.is_empty() {
        Some(fopen(&vsfile, "w").unwrap_or_else(|| err!("cannot open vsfile=%s", vsfile)))
    } else {
        None
    };
    let mut vs = if vseisfp.is_some() {
        Some(Grid::new(nt, nz))
    } else {
        None
    };
    let sseisfp = if !ssfile.is_empty() {
        Some(fopen(&ssfile, "w").unwrap_or_else(|| err!("cannot open ssfile=%s", ssfile)))
    } else {
        None
    };
    let mut ss = if sseisfp.is_some() {
        Some(Grid::new(nt, ns))
    } else {
        None
    };
    let mut dmin = 0.0f32;
    let mut dmax = 0.0f32;
    if !dfile.is_empty() {
        let densityfp = fopen(&dfile, "r").unwrap_or_else(|| err!("cannot open dfile=%s", dfile));
        if fread(&mut od.data, &densityfp) != (nx * nz) as usize {
            err!("error reading dfile=%s", dfile);
        }
        fclose(&densityfp);
        dmin = od[(0, 0)];
        dmax = od[(0, 0)];
        for ix in 0..nx {
            for iz in 0..nz {
                dmin = min(dmin, od[(ix, iz)]);
                dmax = max(dmax, od[(ix, iz)]);
            }
        }
    }
    if dfile.is_empty() || dmin == dmax {
        for ix in 0..nx {
            for iz in 0..nz {
                od[(ix, iz)] = 1.0;
            }
        }
        dmin = 1.0;
        dmax = 1.0;
    }
    for ix in 0..nx {
        for iz in 0..nz {
            dvv[(ix, iz)] = od[(ix, iz)] * dvv[(ix, iz)] * dvv[(ix, iz)];
            od[(ix, iz)] = (1.0 / od[(ix, iz)] as f64) as f32;
            pp[(ix, iz)] = 0.0;
            p[(ix, iz)] = 0.0;
            pm[(ix, iz)] = 0.0;
        }
    }
    let has_od = dmin != dmax;
    if verbose != 0 {
        warn!("nx = %d", nx);
        warn!("dx = %g", dx);
        warn!("nz = %d", nz);
        warn!("dz = %g", dz);
        warn!("nt = %d", nt);
        warn!("dt = %g", dt);
        warn!("tmax = %g", tmax);
        warn!("fmax = %g", fmax);
        warn!("fpeak = %g", fpeak);
        warn!("vmin = %g", vmin);
        warn!("vmax = %g", vmax);
        warn!("mt = %d", mt);
        warn!("pml_max = %g", pml_max);
        warn!("pml_half = %d", pml_thick);
        warn!("pml_thickness = %d", pml_thickness);
        if !has_od {
            warn!("constant density");
        } else {
            warn!("dfile=%s", dfile);
            warn!("dmin = %g", dmin);
            warn!("dmax = %g", dmax);
        }
    }
    let mut pml = kernels::Pml::new(pml_max, pml_thickness);
    if pml_thick != 0 {
        pml.pml_init(nx, nz, dx, dz, dt, &dvv, &od, has_od, verbose);
    }
    let mut cubetr = Trace::new();
    let tdelay = 0.0f32;
    let mut tracr = 0i32;
    let mut t = 0.0f32;
    let mut extended_source = None;
    for it in 0..nt {
        if verbose > 1 {
            warn!("it=%d  t=%g", it, t);
        }
        if ns == 1 {
            ptsrc(
                xs[0], zs[0], nx, dx, fx, nz, dz, fz, dt, t, fmax, fpeak, tdelay, &mut s,
            );
        } else {
            exsrc(
                ns,
                &xs,
                &zs,
                nx,
                dx,
                fx,
                nz,
                dz,
                fz,
                dt,
                t,
                fmax,
                &mut s,
                &mut extended_source,
            );
        }
        tstep2(
            nx,
            dx,
            nz,
            dz,
            dt,
            &dvv,
            &od,
            has_od,
            &s,
            &pm,
            &p,
            &mut pp,
            &absorb_flags,
            pml_thick,
            &mut pml,
        );
        // On arm64 C integer remainder by zero yields the dividend.
        if it.checked_rem(mt).unwrap_or(it) == 0 {
            cubetr.sx = xs[0] as i32;
            cubetr.sdepth = zs[0] as i32;
            cubetr.trid = 130;
            cubetr.ns = nz as u16;
            cubetr.d1 = dz;
            cubetr.d2 = dx;
            cubetr.delrt = (-1000.0 * tdelay as f64) as i16;
            let mut tracl = 0;
            for ix in 0..nx {
                tracl += 1;
                tracr += 1;
                cubetr.offset = (ix as f32 * dx - xs[0]) as i32;
                cubetr.tracl = tracl;
                cubetr.tracr = tracr;
                for iz in 0..nz {
                    cubetr.data[iz as usize] = pp[(ix, iz)];
                }
                fputtr(&stdout(), &mut cubetr);
            }
        }
        if let Some(hs) = hs.as_mut() {
            for ix in 0..nx {
                hs[(ix, it)] = pp[(ix, hs1)];
            }
        }
        if let Some(vs) = vs.as_mut() {
            for iz in 0..nz {
                vs[(iz, it)] = pp[(vs2, iz)];
            }
        }
        if let Some(ss) = ss.as_mut() {
            for is in 0..ns as usize {
                ss[(is as i32, it)] = pp[(ixs[is], izs[is])];
            }
        }
        std::mem::swap(&mut pm, &mut p);
        std::mem::swap(&mut p, &mut pp);
        t += dt;
    }
    if let Some(hs) = hs {
        let fp = hseisfp.unwrap();
        let mut horiztr = Trace::new();
        horiztr.sx = xs[0] as i32;
        horiztr.sdepth = zs[0] as i32;
        horiztr.trid = 1;
        horiztr.ns = nt as u16;
        horiztr.dt = (1000000.0f32 * dt) as i32 as u16;
        horiztr.d2 = dx;
        horiztr.delrt = (-1000.0 * tdelay as f64) as i16;
        for ix in 0..nx {
            horiztr.offset = (ix as f32 * dx - xs[0]) as i32;
            horiztr.tracl = ix + 1;
            horiztr.tracr = ix + 1;
            for it in 0..nt {
                horiztr.data[it as usize] = hs[(ix, it)];
            }
            fputtr(&fp, &mut horiztr);
        }
        fclose(&fp);
    }
    if let Some(vs) = vs {
        let fp = vseisfp.unwrap();
        let mut verttr = Trace::new();
        verttr.trid = 1;
        verttr.ns = nt as u16;
        verttr.sx = xs[0] as i32;
        verttr.sdepth = zs[0] as i32;
        verttr.dt = (1000000.0f32 * dt) as i32 as u16;
        verttr.d2 = dx;
        verttr.delrt = (-1000.0 * tdelay as f64) as i16;
        for iz in 0..nz {
            verttr.offset = (iz as f32 * dz - zs[0]) as i32;
            verttr.tracl = iz + 1;
            verttr.tracr = iz + 1;
            for it in 0..nt {
                verttr.data[it as usize] = vs[(iz, it)];
            }
            fputtr(&fp, &mut verttr);
        }
        fclose(&fp);
    }
    if let Some(ss) = ss {
        let fp = sseisfp.unwrap();
        let mut srctr = Trace::new();
        srctr.trid = 1;
        srctr.ns = nt as u16;
        srctr.dt = (1000000.0f32 * dt) as i32 as u16;
        srctr.d2 = dx;
        srctr.delrt = (-1000.0 * tdelay as f64) as i16;
        for is in 0..ns {
            srctr.sx = xs[is as usize] as i32;
            srctr.sdepth = zs[is as usize] as i32;
            srctr.tracl = is + 1;
            srctr.tracr = is + 1;
            for it in 0..nt {
                srctr.data[it as usize] = ss[(is, it)];
            }
            fputtr(&fp, &mut srctr);
        }
        fclose(&fp);
    }
    exit(CWP_Exit());
}

#[allow(clippy::too_many_arguments, unused_assignments)]
fn ptsrc(
    xs: f32,
    zs: f32,
    nx: i32,
    dx: f32,
    fx: f32,
    nz: i32,
    dz: f32,
    fz: f32,
    dt: f32,
    t: f32,
    fmax: f32,
    fpeak: f32,
    mut tdelay: f32,
    s: &mut Grid,
) {
    for ix in 0..nx {
        for iz in 0..nz {
            s[(ix, iz)] = (0.0 * dt as f64 * fmax as f64) as f32;
        }
    }
    tdelay = (1.0 / fpeak as f64) as f32;
    if t as f64 > 2.0 * tdelay as f64 {
        return;
    }
    let ts = ricker(t - tdelay, fpeak);
    let xsn = (xs - fx) / dx;
    let zsn = (zs - fz) / dz;
    let ixs = nint(xsn as f64);
    let izs = nint(zsn as f64);
    for ix in max(0, ixs - 3)..=min(nx - 1, ixs + 3) {
        for iz in max(0, izs - 3)..=min(nz - 1, izs + 3) {
            let xn = ix as f32 - xsn;
            let zn = iz as f32 - zsn;
            s[(ix, iz)] = (ts as f64 * ((-xn * xn - zn * zn) as f64).exp()) as f32;
        }
    }
}

fn ricker(t: f32, fpeak: f32) -> f32 {
    let x = (std::f64::consts::PI * fpeak as f64 * t as f64) as f32;
    let xx = x * x;
    ((-xx as f64).exp() * (1.0 - 2.0 * xx as f64)) as f32
}

struct ExtendedSource {
    vs: Vec<f32>,
    xsd: Vec<[f32; 4]>,
    zsd: Vec<[f32; 4]>,
}

/// Update source pressure for an extended source (Dave Hale, 03/01/90).
#[allow(clippy::too_many_arguments)]
fn exsrc(
    ns: i32,
    xs: &[f32],
    zs: &[f32],
    nx: i32,
    dx: f32,
    fx: f32,
    nz: i32,
    dz: f32,
    fz: f32,
    dt: f32,
    t: f32,
    fmax: f32,
    s: &mut Grid,
    made: &mut Option<ExtendedSource>,
) {
    if made.is_none() {
        let mut vs = vec![0.0; ns as usize];
        let mut xsd = vec![[0.0; 4]; ns as usize];
        let mut zsd = vec![[0.0; 4]; ns as usize];
        for is in 0..ns {
            vs[is as usize] = is as f32;
        }
        cmonot(&vs, xs, &mut xsd);
        cmonot(&vs, zs, &mut zsd);
        *made = Some(ExtendedSource { vs, xsd, zsd });
    }
    let ExtendedSource { vs, xsd, zsd } = made.as_ref().unwrap();
    for ix in 0..nx {
        for iz in 0..nz {
            s[(ix, iz)] = (0.0 * dt as f64) as f32;
        }
    }
    let sigma = (0.25 / fmax as f64) as f32;
    let tbias = (3.0 * sigma as f64) as f32;
    let ascale = (-0.5f64.exp() / sigma as f64) as f32;
    let tscale = (0.5 / (sigma * sigma) as f64) as f32;
    if t as f64 > 2.0 * tbias as f64 {
        return;
    }
    let ts = ((ascale * (t - tbias)) as f64 * ((-tscale * (t - tbias) * (t - tbias)) as f64).exp())
        as f32;
    let mut v = vs[0];
    let mut distprev = 0.0f32;
    let mut dv = 1.0f32;
    while dv != 0.0 {
        let mut xv = [0.0f32];
        let mut zv = [0.0f32];
        let mut dxdv = [0.0f32];
        let mut dzdv = [0.0f32];
        intcub(0, vs, xsd, &[v], &mut xv);
        intcub(0, vs, zsd, &[v], &mut zv);
        intcub(1, vs, xsd, &[v], &mut dxdv);
        intcub(1, vs, zsd, &[v], &mut dzdv);
        let (xv, zv, dxdv, dzdv) = (xv[0], zv[0], dxdv[0], dzdv[0]);
        if dxdv == 0.0 {
            dv = dz / abs(dzdv);
        } else if dzdv == 0.0 {
            dv = dx / abs(dxdv);
        } else {
            dv = min(dz / abs(dzdv), dx / abs(dxdv));
        }
        if v + dv > vs[ns as usize - 1] {
            dv = vs[ns as usize - 1] - v;
        }
        let dist = (dv as f64 * ((dzdv * dzdv + dxdv * dxdv) as f64).sqrt()
            / ((dx * dx + dz * dz) as f64).sqrt()) as f32;
        let amp = ((dist + distprev) as f64 / 2.0) as f32;
        let xvn = (xv - fx) / dx;
        let zvn = (zv - fz) / dz;
        let ixv = nint(xvn as f64);
        let izv = nint(zvn as f64);
        for ix in max(0, ixv - 3)..=min(nx - 1, ixv + 3) {
            for iz in max(0, izv - 3)..=min(nz - 1, izv + 3) {
                let xn = ix as f32 - xvn;
                let zn = iz as f32 - zvn;
                s[(ix, iz)] = (s[(ix, iz)] as f64
                    + (ts * amp) as f64 * ((-xn * xn - zn * zn) as f64).exp())
                    as f32;
            }
        }
        distprev = dist;
        v += dv;
    }
}

#[allow(clippy::too_many_arguments)]
fn tstep2(
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
    abs: &[i32; 4],
    pml_thick: i32,
    pml: &mut kernels::Pml,
) {
    if has_od && dx != dz {
        kernels::star1(nx, dx, nz, dz, dt, dvv, od, has_od, s, pm, p, pp);
    } else if has_od && dx == dz {
        kernels::star2(nx, dx, nz, dz, dt, dvv, od, has_od, s, pm, p, pp);
    } else if !has_od && dx != dz {
        kernels::star3(nx, dx, nz, dz, dt, dvv, od, has_od, s, pm, p, pp);
    } else {
        kernels::star4(nx, dx, nz, dz, dt, dvv, od, has_od, s, pm, p, pp);
    }
    if pml_thick == 0 {
        kernels::absorb(nx, dx, nz, dz, dt, dvv, od, has_od, pm, p, pp, abs);
    } else {
        pml.pml_absorb(nx, dx, nz, dz, dt, dvv, od, has_od, pm, p, pp, abs);
    }
}
