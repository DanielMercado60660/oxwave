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

//! SUFDMOD2 - Finite-Difference MODeling (2nd order) for acoustic wave equation.
//! Port of su/main/synthetics_waveforms_testpatterns/sufdmod2.c, Revision 1.29.
//! Authors: CWP: Dave Hale; modified for SU by John Stockwell, 1993.
//! U Houston: plane wave and monochromatic source options, Chris Liner, 2010.
//! Modifications: Tony Kocurko, Memorial University in Newfoundland and Labrador:
//! source location files and corrected source z sample calculation.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::cwp::{abs, max, min, nint, PI};
use oxwave_f_cubicspline::cmonot;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUFDMOD2 - Finite-Difference MODeling (2nd order) for acoustic wave equation",
    " \t\t\t\t\t\t\t\t\t",
    " sufdmod2 <vfile >wfile nx= nz= tmax= xs= zs= [optional parameters]\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " <vfile\t\tfile containing velocity[nx][nz]\t\t",
    " >wfile\t\tfile containing waves[nx][nz] for time steps\t",
    " nx=\t\t\tnumber of x samples (2nd dimension)\t\t",
    " nz=\t\t\tnumber of z samples (1st dimension)\t\t",
    " xs=\t\t\tx coordinates of source, or, alternatively, the name",
    "\t\t\tof a file that contains the x- and z-coordinates,",
    "\t\t\twith the number of pairs as the first record and",
    "\t\t\tthe actual pairs of (x,z) locations following.  ",
    " zs=\t\t\tz coordinates of source\t\t\t\t",
    " tmax=\t\t\tmaximum time\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " sstrength=1.0\t\tstrength of source\t\t\t\t",
    " pw=0\t\t\tuse point or extended source geometry parameters",
    " \t\t\t=1  use horizontal plane wave source \t\t",
    " pwt=20\t\tamp taper on ends of line src (in grid points)  ",
    " mono=0\t\tuse ricker wavelet as source function \t\t",
    " \t\t\t=1  use single frequency src (freq=2*fpeak)\t",
    " nt=1+tmax/dt\t\tnumber of time samples (dt determined for stability)",
    " mt=1\t\t\tnumber of time steps (dt) per output time step\t",
    " \t\t\t\t\t\t\t\t\t",
    " dx=1.0\t\tx sampling interval\t\t\t\t",
    " fx=0.0\t\tfirst x sample\t\t\t\t\t",
    " dz=1.0\t\tz sampling interval\t\t\t\t",
    " fz=0.0\t\tfirst z sample\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " fmax = vmin/(10.0*h)\tmaximum frequency in source \t\t\t",
    " fpeak=0.5*fmax\tpeak frequency in ricker wavelet\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dfile=\t\tinput file containing density[nx][nz]\t\t",
    " vsx=\t\t\tx coordinate of vertical line of seismograms\t",
    " hsz=\t\t\tz coordinate of horizontal line of seismograms\t",
    " vsfile=\t\toutput file for vertical line of seismograms[nz][nt]",
    " hsfile=\t\toutput file for horizontal line of seismograms[nx][nt]",
    " ssfile=\t\toutput file for source point seismograms[nt]\t",
    " verbose=0\t\t=1 for diagnostic messages, =2 for more\t\t",
    " abs=1,1,1,1\t\tabsorbing boundary conditions on top,left,bottom,right",
    "\t\t\tsides of the model. \t\t\t\t",
    " \t\t\t=0,1,1,1 for free surface condition on the top\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " This program uses the traditional explicit second order differencing\t",
    " method. \t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut nx, mut nz) = (0i32, 0i32);
    let mut tmax = 0.0f32;
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx!\n");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!\n");
    }
    if getparfloat("tmax", &mut tmax) == 0 {
        err!("must specify tmax!\n");
    }
    // An xs argument naming a readable file takes precedence over numeric xs/zs.
    let mut spfile = String::new();
    getparstring("xs", &mut spfile);
    let srpntfp = fopen(&spfile, "r");
    let (nxs, nzs) = if let Some(ref fp) = srpntfp {
        let Some(n) = fscanf_i32(fp) else {
            fclose(fp);
            err!("error reading nxs from spfile=%s\n", spfile);
        };
        (n, n)
    } else {
        (countparval("xs"), countparval("zs"))
    };
    if nxs != nzs {
        err!("number of xs = %d must equal number of zs = %d\n", nxs, nzs);
    }
    let mut ns = nxs;
    let mut pw = 0i32;
    getparint("pw", &mut pw);
    if pw == 1 {
        ns = nx;
    }
    if ns == 0 {
        err!("must specify xs and zs!\n");
    }
    let mut xs = vec![0.0f32; ns as usize];
    let mut zs = vec![0.0f32; ns as usize];
    let mut ixs = vec![0i32; ns as usize];
    let mut izs = vec![0i32; ns as usize];
    if let Some(ref fp) = srpntfp {
        for it in 0..ns as usize {
            let x = fscanf_f32(fp);
            let z = if x.is_some() { fscanf_f32(fp) } else { None };
            if x.is_none() || z.is_none() {
                fclose(fp);
                err!(
                    "error reading xs[%d], zs[%d] from spfile=%s\n",
                    it as i32,
                    it as i32,
                    spfile
                );
            }
            xs[it] = x.unwrap();
            zs[it] = z.unwrap();
        }
        fclose(fp);
    } else {
        getparfloat("xs", &mut xs);
        getparfloat("zs", &mut zs);
    }
    let mut sstrength = 1.0f32;
    getparfloat("sstrength", &mut sstrength);
    let ztmp = if pw == 1 { zs[0] } else { 0.0 };
    let nabs = countparval("abs");
    let mut bounds = [1i32; 4];
    if nabs == 4 {
        getparint("abs", &mut bounds[..]);
    } else if nabs != 0 {
        warn!("Number of abs %d, using abs=1,1,1,1", nabs);
    }
    let (mut nt, mut mt) = (0i32, 1i32);
    getparint("nt", &mut nt);
    getparint("mt", &mut mt);
    let (mut dx, mut fx, mut dz, mut fz) = (1.0f32, 0.0f32, 1.0f32, 0.0f32);
    getparfloat("dx", &mut dx);
    getparfloat("fx", &mut fx);
    getparfloat("dz", &mut dz);
    getparfloat("fz", &mut fz);
    let (mut pwt, mut mono) = (20i32, 0i32);
    getparint("pwt", &mut pwt);
    getparint("mono", &mut mono);
    for is in 0..ns as usize {
        ixs[is] = nint(((xs[is] - fx) / dx) as f64);
        izs[is] = nint(((zs[is] - fz) / dz) as f64);
    }
    let mut hsz = 0.0f32;
    getparfloat("hsz", &mut hsz);
    let hs1 = nint(((hsz - fz) / dz) as f64);
    let mut vsx = 0.0f32;
    getparfloat("vsx", &mut vsx);
    let vs2 = nint(((vsx - fx) / dx) as f64);
    let mut verbose = 0i32;
    getparint("verbose", &mut verbose);
    let (mut dfile, mut hsfile, mut vsfile, mut ssfile) =
        (String::new(), String::new(), String::new(), String::new());
    getparstring("dfile", &mut dfile);
    getparstring("hsfile", &mut hsfile);
    getparstring("vsfile", &mut vsfile);
    getparstring("ssfile", &mut ssfile);
    // The C program never calls checkpars().
    let (nx, nz) = (nx as usize, nz as usize);
    let mut s = vec![vec![0.0f32; nz]; nx];
    // fread addresses dvv[0] as one contiguous allocation.
    let mut velocities = vec![0.0f32; nx * nz];
    fread(&mut velocities, &stdin());
    let mut dvv: Vec<Vec<f32>> = velocities.chunks(nz).map(|r| r.to_vec()).collect();
    let mut od = vec![vec![0.0f32; nz]; nx];
    let mut pm = vec![vec![0.0f32; nz]; nx];
    let mut p = vec![vec![0.0f32; nz]; nx];
    let mut pp = vec![vec![0.0f32; nz]; nx];
    let (mut vmin, mut vmax) = (dvv[0][0], dvv[0][0]);
    for ix in 0..nx {
        for iz in 0..nz {
            vmin = min(vmin, dvv[ix][iz]);
            if verbose == 1 && dvv[ix][iz] == 0.0 {
                warn!("v=0 at (ix,iz)=(%i,%i)", ix as i32, iz as i32);
            }
            vmax = max(vmax, dvv[ix][iz]);
        }
    }
    let h = min(abs(dx), abs(dz));
    let dt = (h as f64 / (2.0 * vmax as f64)) as f32;
    warn!("stable dt=%g", dt);
    let mut fmax = (vmin as f64 / (10.0 * h as f64)) as f32;
    getparfloat("fmax", &mut fmax);
    let mut fpeak = (0.5 * fmax as f64) as f32;
    getparfloat("fpeak", &mut fpeak);
    if nt == 0 {
        nt = (1.0f32 + tmax / dt) as i32;
    }
    let hseisfp = if hsfile.is_empty() {
        None
    } else {
        Some(fopen(&hsfile, "w").unwrap_or_else(|| err!("cannot open hsfile=%s\n", hsfile)))
    };
    let mut hs = hseisfp
        .as_ref()
        .map(|_| vec![vec![0.0f32; nt as usize]; nx]);
    let vseisfp = if vsfile.is_empty() {
        None
    } else {
        Some(fopen(&vsfile, "w").unwrap_or_else(|| err!("cannot open vsfile=%s\n", vsfile)))
    };
    let mut vs = vseisfp
        .as_ref()
        .map(|_| vec![vec![0.0f32; nt as usize]; nz]);
    let sseisfp = if ssfile.is_empty() {
        None
    } else {
        Some(fopen(&ssfile, "w").unwrap_or_else(|| err!("cannot open ssfile=%s\n", ssfile)))
    };
    let mut ss = sseisfp
        .as_ref()
        .map(|_| vec![vec![0.0f32; nt as usize]; ns as usize]);
    let (mut dmin, mut dmax) = (0.0f32, 0.0f32);
    if !dfile.is_empty() {
        let fp = fopen(&dfile, "r").unwrap_or_else(|| err!("cannot open dfile=%s\n", dfile));
        let mut densities = vec![0.0f32; nx * nz];
        if fread(&mut densities, &fp) != nx * nz {
            err!("error reading dfile=%s\n", dfile);
        }
        fclose(&fp);
        for ix in 0..nx {
            for iz in 0..nz {
                od[ix][iz] = densities[ix * nz + iz];
            }
        }
        dmin = od[0][0];
        dmax = od[0][0];
        for ix in 0..nx {
            for iz in 0..nz {
                dmin = min(dmin, od[ix][iz]);
                dmax = max(dmax, od[ix][iz]);
            }
        }
    }
    if dfile.is_empty() || dmin == dmax {
        for ix in 0..nx {
            for iz in 0..nz {
                od[ix][iz] = 1.0;
            }
        }
        dmin = 1.0;
        dmax = 1.0;
    }
    for ix in 0..nx {
        for iz in 0..nz {
            dvv[ix][iz] = od[ix][iz] * dvv[ix][iz] * dvv[ix][iz];
            od[ix][iz] = (1.0 / od[ix][iz] as f64) as f32;
            pp[ix][iz] = 0.0;
            p[ix][iz] = 0.0;
            pm[ix][iz] = 0.0;
        }
    }
    let od = if dmin == dmax { None } else { Some(od) };
    if verbose != 0 {
        fprintf!(stderr(), "nx = %d\n", nx as i32);
        fprintf!(stderr(), "dx = %g\n", dx);
        fprintf!(stderr(), "nz = %d\n", nz as i32);
        fprintf!(stderr(), "dz = %g\n", dz);
        fprintf!(stderr(), "nt = %d\n", nt);
        fprintf!(stderr(), "dt = %g\n", dt);
        fprintf!(stderr(), "tmax = %g\n", tmax);
        fprintf!(stderr(), "fmax = %g\n", fmax);
        fprintf!(stderr(), "fpeak = %g\n", fpeak);
        fprintf!(stderr(), "vmin = %g\n", vmin);
        fprintf!(stderr(), "vmax = %g\n", vmax);
        fprintf!(stderr(), "mt = %d\n", mt);
        if dmin == dmax {
            fprintf!(stderr(), "constant density\n");
        } else {
            fprintf!(stderr(), "dfile=%s\n", dfile);
            fprintf!(stderr(), "dmin = %g\n", dmin);
            fprintf!(stderr(), "dmax = %g\n", dmax);
        }
    }
    if pw == 1 {
        for ix in 0..nx {
            xs[ix] = ix as f32 * dx;
            ixs[ix] = ix as i32;
            zs[ix] = ztmp;
            izs[ix] = (ztmp / dz) as i32;
        }
    }
    let tdelay = 0.0f32; // ptsrc receives this by value; the C caller keeps zero.
    let mut cubetr = Trace::new();
    let mut tracr = 0i32;
    let mut extended = None;
    let mut t = 0.0f32;
    for it in 0..nt.max(0) as usize {
        if verbose > 1 {
            fprintf!(stderr(), "it=%d  t=%g\n", it as i32, t);
        }
        if ns == 1 {
            ptsrc(
                sstrength, xs[0], zs[0], nx, dx, fx, nz, dz, fz, dt, t, fmax, fpeak, mono, &mut s,
            );
        } else {
            if extended.is_none() {
                extended = Some(Extended::new(&xs, &zs));
            }
            extended
                .as_ref()
                .unwrap()
                .exsrc(nx, dx, fx, nz, dz, fz, dt, t, fmax, pwt, mono, &mut s);
        }
        tstep2(
            nx,
            dx,
            nz,
            dz,
            dt,
            &dvv,
            od.as_deref(),
            &s,
            &pm,
            &p,
            &mut pp,
            &bounds,
        );
        // ARM sdiv by zero yields quotient zero; msub leaves the dividend as remainder.
        if (it as i32).checked_rem(mt).unwrap_or(it as i32) == 0 {
            cubetr.sx = xs[0] as i32;
            cubetr.sdepth = zs[0] as i32;
            cubetr.trid = 30;
            cubetr.ns = nz as u16;
            cubetr.d1 = dz;
            cubetr.d2 = dx;
            cubetr.delrt = (-1000.0 * tdelay as f64) as i16;
            for ix in 0..nx {
                tracr = tracr.wrapping_add(1);
                cubetr.offset = (ix as f32 * dx - xs[0]) as i32;
                cubetr.gx = (ix as f32 * dx) as i32;
                cubetr.tracl = (ix + 1) as i32;
                cubetr.tracr = tracr;
                cubetr.data[..nz].copy_from_slice(&pp[ix]);
                fputtr(&stdout(), &mut cubetr);
            }
        }
        if let Some(ref mut hs) = hs {
            for ix in 0..nx {
                hs[ix][it] = pp[ix][hs1 as usize];
            }
        }
        if let Some(ref mut vs) = vs {
            for iz in 0..nz {
                vs[iz][it] = pp[vs2 as usize][iz];
            }
        }
        if let Some(ref mut ss) = ss {
            for is in 0..ns as usize {
                ss[is][it] = pp[ixs[is] as usize][izs[is] as usize];
            }
        }
        std::mem::swap(&mut pm, &mut p);
        std::mem::swap(&mut p, &mut pp);
        t += dt;
    }
    if let Some(hs) = hs {
        let mut tr = Trace::new();
        tr.sx = xs[0] as i32;
        tr.sdepth = zs[0] as i32;
        tr.trid = 1;
        tr.ns = nt as u16;
        tr.dt = (1000000.0f32 * dt) as u32 as u16;
        tr.d2 = dx;
        tr.delrt = (-1000.0 * tdelay as f64) as i16;
        for ix in 0..nx {
            tr.offset = (ix as f32 * dx - xs[0]) as i32;
            tr.gx = (ix as f32 * dx) as i32;
            tr.tracl = (ix + 1) as i32;
            tr.tracr = (ix + 1) as i32;
            tr.data[..nt as usize].copy_from_slice(&hs[ix]);
            fputtr(hseisfp.as_ref().unwrap(), &mut tr);
        }
        fclose(hseisfp.as_ref().unwrap());
    }
    if let Some(vs) = vs {
        let mut tr = Trace::new();
        tr.trid = 1;
        tr.ns = nt as u16;
        tr.sx = xs[0] as i32;
        tr.sdepth = zs[0] as i32;
        tr.dt = (1000000.0f32 * dt) as u32 as u16;
        tr.d2 = dx;
        tr.delrt = (-1000.0 * tdelay as f64) as i16;
        for iz in 0..nz {
            tr.offset = (iz as f32 * dz - zs[0]) as i32;
            tr.tracl = (iz + 1) as i32;
            tr.tracr = (iz + 1) as i32;
            tr.data[..nt as usize].copy_from_slice(&vs[iz]);
            fputtr(vseisfp.as_ref().unwrap(), &mut tr);
        }
        fclose(vseisfp.as_ref().unwrap());
    }
    if let Some(ss) = ss {
        let mut tr = Trace::new();
        tr.trid = 1;
        tr.ns = nt as u16;
        tr.dt = (1000000.0f32 * dt) as u32 as u16;
        tr.d2 = dx;
        tr.delrt = (-1000.0 * tdelay as f64) as i16;
        for is in 0..ns as usize {
            tr.sx = xs[is] as i32;
            tr.sdepth = zs[is] as i32;
            tr.tracl = (is + 1) as i32;
            tr.tracr = (is + 1) as i32;
            tr.data[..nt as usize].copy_from_slice(&ss[is]);
            fputtr(sseisfp.as_ref().unwrap(), &mut tr);
        }
        fclose(sseisfp.as_ref().unwrap());
    }
    exit(CWP_Exit());
}

type Grid = Vec<Vec<f32>>;

#[allow(clippy::too_many_arguments)]
fn tstep2(
    nx: usize,
    dx: f32,
    nz: usize,
    dz: f32,
    dt: f32,
    dvv: &[Vec<f32>],
    od: Option<&[Vec<f32>]>,
    s: &[Vec<f32>],
    pm: &[Vec<f32>],
    p: &[Vec<f32>],
    pp: &mut Grid,
    bounds: &[i32; 4],
) {
    if let Some(od) = od {
        if dx != dz {
            star1(nx, dx, nz, dz, dt, dvv, od, s, pm, p, pp);
        } else {
            star2(nx, dx, nz, dz, dt, dvv, od, s, pm, p, pp);
        }
    } else if dx != dz {
        star3(nx, dx, nz, dz, dt, dvv, s, pm, p, pp);
    } else {
        star4(nx, dx, nz, dz, dt, dvv, s, pm, p, pp);
    }
    absorb(nx, dx, nz, dz, dt, dvv, od, pm, p, pp, bounds);
}

// Convolve with finite-difference star for variable density and dx!=dz.
#[allow(clippy::too_many_arguments)]
fn star1(
    nx: usize,
    dx: f32,
    nz: usize,
    dz: f32,
    dt: f32,
    dvv: &[Vec<f32>],
    od: &[Vec<f32>],
    s: &[Vec<f32>],
    pm: &[Vec<f32>],
    p: &[Vec<f32>],
    pp: &mut Grid,
) {
    let xscale1 = (dt * dt) / (dx * dx);
    let zscale1 = (dt * dt) / (dz * dz);
    let xscale2 = (0.25 * xscale1 as f64) as f32;
    let zscale2 = (0.25 * zscale1 as f64) as f32;
    for ix in 1..nx - 1 {
        for iz in 1..nz - 1 {
            pp[ix][iz] = (2.0 * p[ix][iz] as f64 - pm[ix][iz] as f64
                + dvv[ix][iz] as f64
                    * (od[ix][iz] as f64
                        * (xscale1 as f64
                            * ((p[ix + 1][iz] + p[ix - 1][iz]) as f64 - 2.0 * p[ix][iz] as f64)
                            + zscale1 as f64
                                * ((p[ix][iz + 1] + p[ix][iz - 1]) as f64
                                    - 2.0 * p[ix][iz] as f64))
                        + (xscale2
                            * ((od[ix + 1][iz] - od[ix - 1][iz]) * (p[ix + 1][iz] - p[ix - 1][iz]))
                            + zscale2
                                * ((od[ix][iz + 1] - od[ix][iz - 1])
                                    * (p[ix][iz + 1] - p[ix][iz - 1])))
                            as f64)
                + s[ix][iz] as f64) as f32;
        }
    }
}

// Convolve with finite-difference star for variable density and dx==dz.
#[allow(clippy::too_many_arguments)]
fn star2(
    nx: usize,
    dx: f32,
    nz: usize,
    dz: f32,
    dt: f32,
    dvv: &[Vec<f32>],
    od: &[Vec<f32>],
    s: &[Vec<f32>],
    pm: &[Vec<f32>],
    p: &[Vec<f32>],
    pp: &mut Grid,
) {
    if dx != dz {
        fprintf!(stderr(), "ASSERT FAILED: dx != dz in star2\n");
    }
    let scale1 = (dt * dt) / (dx * dx);
    let scale2 = (0.25 * scale1 as f64) as f32;
    for ix in 1..nx - 1 {
        for iz in 1..nz - 1 {
            pp[ix][iz] = (2.0 * p[ix][iz] as f64 - pm[ix][iz] as f64
                + dvv[ix][iz] as f64
                    * (od[ix][iz] as f64
                        * (scale1 as f64
                            * ((p[ix + 1][iz] + p[ix - 1][iz] + p[ix][iz + 1] + p[ix][iz - 1])
                                as f64
                                - 4.0 * p[ix][iz] as f64))
                        + (scale2
                            * ((od[ix + 1][iz] - od[ix - 1][iz]) * (p[ix + 1][iz] - p[ix - 1][iz])
                                + (od[ix][iz + 1] - od[ix][iz - 1])
                                    * (p[ix][iz + 1] - p[ix][iz - 1])))
                            as f64)
                + s[ix][iz] as f64) as f32;
        }
    }
}

// Convolve with finite-difference star for density==1.0 and dx!=dz.
#[allow(clippy::too_many_arguments)]
fn star3(
    nx: usize,
    dx: f32,
    nz: usize,
    dz: f32,
    dt: f32,
    dvv: &[Vec<f32>],
    s: &[Vec<f32>],
    pm: &[Vec<f32>],
    p: &[Vec<f32>],
    pp: &mut Grid,
) {
    let xscale = (dt * dt) / (dx * dx);
    let zscale = (dt * dt) / (dz * dz);
    for ix in 1..nx - 1 {
        for iz in 1..nz - 1 {
            pp[ix][iz] = (2.0 * p[ix][iz] as f64 - pm[ix][iz] as f64
                + dvv[ix][iz] as f64
                    * (xscale as f64
                        * ((p[ix + 1][iz] + p[ix - 1][iz]) as f64 - 2.0 * p[ix][iz] as f64)
                        + zscale as f64
                            * ((p[ix][iz + 1] + p[ix][iz - 1]) as f64 - 2.0 * p[ix][iz] as f64))
                + s[ix][iz] as f64) as f32;
        }
    }
}

// Convolve with finite-difference star for density==1.0 and dx==dz.
#[allow(clippy::too_many_arguments)]
fn star4(
    nx: usize,
    dx: f32,
    nz: usize,
    dz: f32,
    dt: f32,
    dvv: &[Vec<f32>],
    s: &[Vec<f32>],
    pm: &[Vec<f32>],
    p: &[Vec<f32>],
    pp: &mut Grid,
) {
    if dz != dx {
        fprintf!(stderr(), "ASSERT FAILED: dz !=  dx in star4\n");
    }
    let scale = (dt * dt) / (dx * dz);
    for ix in 1..nx - 1 {
        for iz in 1..nz - 1 {
            pp[ix][iz] = (2.0 * p[ix][iz] as f64 - pm[ix][iz] as f64
                + (scale * dvv[ix][iz]) as f64
                    * ((p[ix + 1][iz] + p[ix - 1][iz] + p[ix][iz + 1] + p[ix][iz - 1]) as f64
                        - 4.0 * p[ix][iz] as f64)
                + s[ix][iz] as f64) as f32;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn absorb(
    nx: usize,
    dx: f32,
    nz: usize,
    dz: f32,
    dt: f32,
    dvv: &[Vec<f32>],
    od: Option<&[Vec<f32>]>,
    pm: &[Vec<f32>],
    p: &[Vec<f32>],
    pp: &mut Grid,
    bounds: &[i32; 4],
) {
    // Solve for upper boundary.
    let iz = 1;
    for ix in 0..nx {
        if bounds[0] != 0 {
            let ovs = if let Some(od) = od {
                (1.0 / (od[ix][iz] * dvv[ix][iz]) as f64) as f32
            } else {
                (1.0 / dvv[ix][iz] as f64) as f32
            };
            let ov = (ovs as f64).sqrt() as f32;
            let dpdx = if ix == 0 {
                (p[1][iz] - p[0][iz]) / dx
            } else if ix == nx - 1 {
                (p[nx - 1][iz] - p[nx - 2][iz]) / dx
            } else {
                ((p[ix + 1][iz] - p[ix - 1][iz]) as f64 / (2.0 * dx as f64)) as f32
            };
            let dpdt = ((pp[ix][iz] - pm[ix][iz]) as f64 / (2.0 * dt as f64)) as f32;
            let dpdxs = dpdx * dpdx;
            let dpdts = dpdt * dpdt;
            let cosa = if ovs * dpdts > dpdxs {
                (1.0 - (dpdxs / (ovs * dpdts)) as f64).sqrt() as f32
            } else {
                0.0
            };
            let beta = ov * dz / dt * cosa;
            let gamma = ((1.0 - beta as f64) / (1.0 + beta as f64)) as f32;
            pp[ix][iz - 1] = gamma * (pp[ix][iz] - p[ix][iz - 1]) + p[ix][iz];
        } else {
            pp[ix][iz - 1] = 0.0;
        }
    }
    // Extrapolate along left boundary.
    let ix = 1;
    for iz in 0..nz {
        if bounds[1] != 0 {
            let ovs = if let Some(od) = od {
                (1.0 / (od[ix][iz] * dvv[ix][iz]) as f64) as f32
            } else {
                (1.0 / dvv[ix][iz] as f64) as f32
            };
            let ov = (ovs as f64).sqrt() as f32;
            let dpdz = if iz == 0 {
                (p[ix][1] - p[ix][0]) / dz
            } else if iz == nz - 1 {
                (p[ix][nz - 1] - p[ix][nz - 2]) / dz
            } else {
                ((p[ix][iz + 1] - p[ix][iz - 1]) as f64 / (2.0 * dz as f64)) as f32
            };
            let dpdt = ((pp[ix][iz] - pm[ix][iz]) as f64 / (2.0 * dt as f64)) as f32;
            let dpdzs = dpdz * dpdz;
            let dpdts = dpdt * dpdt;
            let cosa = if ovs * dpdts > dpdzs {
                (1.0 - (dpdzs / (ovs * dpdts)) as f64).sqrt() as f32
            } else {
                0.0
            };
            let beta = ov * dx / dt * cosa;
            let gamma = ((1.0 - beta as f64) / (1.0 + beta as f64)) as f32;
            pp[ix - 1][iz] = gamma * (pp[ix][iz] - p[ix - 1][iz]) + p[ix][iz];
        } else {
            pp[ix - 1][iz] = 0.0;
        }
    }
    // Extrapolate along lower boundary.
    // Keep the C signed index: with nz=1 and this boundary disabled,
    // iz=-1 is never dereferenced, and iz+1 selects the only sample.
    let iz = nz as i32 - 2;
    for ix in 0..nx {
        if bounds[2] != 0 {
            let iz = iz as usize;
            let ovs = if let Some(od) = od {
                (1.0 / (od[ix][iz] * dvv[ix][iz]) as f64) as f32
            } else {
                (1.0 / dvv[ix][iz] as f64) as f32
            };
            let ov = (ovs as f64).sqrt() as f32;
            let dpdx = if ix == 0 {
                (p[1][iz] - p[0][iz]) / dx
            } else if ix == nx - 1 {
                (p[nx - 1][iz] - p[nx - 2][iz]) / dx
            } else {
                ((p[ix + 1][iz] - p[ix - 1][iz]) as f64 / (2.0 * dx as f64)) as f32
            };
            let dpdt = ((pp[ix][iz] - pm[ix][iz]) as f64 / (2.0 * dt as f64)) as f32;
            let dpdxs = dpdx * dpdx;
            let dpdts = dpdt * dpdt;
            let cosa = if ovs * dpdts > dpdxs {
                (1.0 - (dpdxs / (ovs * dpdts)) as f64).sqrt() as f32
            } else {
                0.0
            };
            let beta = ov * dz / dt * cosa;
            let gamma = ((1.0 - beta as f64) / (1.0 + beta as f64)) as f32;
            pp[ix][iz + 1] = gamma * (pp[ix][iz] - p[ix][iz + 1]) + p[ix][iz];
        } else {
            pp[ix][(iz + 1) as usize] = 0.0;
        }
    }
    // Extrapolate along right boundary.
    // Likewise, nx=1 is valid here when the right boundary is disabled.
    let ix = nx as i32 - 2;
    for iz in 0..nz {
        if bounds[3] != 0 {
            let ix = ix as usize;
            let ovs = if let Some(od) = od {
                (1.0 / (od[ix][iz] * dvv[ix][iz]) as f64) as f32
            } else {
                (1.0 / dvv[ix][iz] as f64) as f32
            };
            let ov = (ovs as f64).sqrt() as f32;
            let dpdz = if iz == 0 {
                (p[ix][1] - p[ix][0]) / dz
            } else if iz == nz - 1 {
                (p[ix][nz - 1] - p[ix][nz - 2]) / dz
            } else {
                ((p[ix][iz + 1] - p[ix][iz - 1]) as f64 / (2.0 * dz as f64)) as f32
            };
            let dpdt = ((pp[ix][iz] - pm[ix][iz]) as f64 / (2.0 * dt as f64)) as f32;
            let dpdzs = dpdz * dpdz;
            let dpdts = dpdt * dpdt;
            let cosa = if ovs * dpdts > dpdzs {
                (1.0 - (dpdzs / (ovs * dpdts)) as f64).sqrt() as f32
            } else {
                0.0
            };
            let beta = ov * dx / dt * cosa;
            let gamma = ((1.0 - beta as f64) / (1.0 + beta as f64)) as f32;
            pp[ix + 1][iz] = gamma * (pp[ix][iz] - p[ix + 1][iz]) + p[ix][iz];
        } else {
            pp[(ix + 1) as usize][iz] = 0.0;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn ptsrc(
    sstrength: f32,
    xs: f32,
    zs: f32,
    nx: usize,
    dx: f32,
    fx: f32,
    nz: usize,
    dz: f32,
    fz: f32,
    dt: f32,
    t: f32,
    fmax: f32,
    fpeak: f32,
    mono: i32,
    s: &mut Grid,
) {
    for ix in 0..nx {
        for iz in 0..nz {
            s[ix][iz] = (0.0 * dt as f64 * fmax as f64) as f32;
        }
    }
    let tdelay = (1.0 / fpeak as f64) as f32;
    if t as f64 > 2.0 * tdelay as f64 && mono == 0 {
        return;
    }
    let ts = ricker(t - tdelay, fpeak, mono);
    let xsn = (xs - fx) / dx;
    let zsn = (zs - fz) / dz;
    let ixs = nint(xsn as f64);
    let izs = nint(zsn as f64);
    for ix in max(0, ixs - 3)..=min(nx as i32 - 1, ixs + 3) {
        for iz in max(0, izs - 3)..=min(nz as i32 - 1, izs + 3) {
            let xn = ix as f32 - xsn;
            let zn = iz as f32 - zsn;
            s[ix as usize][iz as usize] =
                ((sstrength * ts) as f64 * ((-xn * xn - zn * zn) as f64).exp()) as f32;
        }
    }
}

fn ricker(t: f32, fpeak: f32, mono: i32) -> f32 {
    let x = (PI * fpeak as f64 * t as f64) as f32;
    let xx = x * x;
    if mono == 0 {
        ((-xx as f64).exp() * (1.0 - 2.0 * xx as f64)) as f32
    } else {
        ((4.0f32 * x) as f64).sin() as f32
    }
}

struct Extended {
    vs: Vec<f32>,
    xsd: Vec<[f32; 4]>,
    zsd: Vec<[f32; 4]>,
}
impl Extended {
    fn new(xs: &[f32], zs: &[f32]) -> Self {
        // Make the spline coefficients once, as in exsrc's static made block.
        let ns = xs.len();
        let mut vs = vec![0.0f32; ns];
        let mut xsd = vec![[0.0f32; 4]; ns];
        let mut zsd = vec![[0.0f32; 4]; ns];
        for is in 0..ns {
            vs[is] = is as f32;
        }
        cmonot(&vs, xs, &mut xsd);
        cmonot(&vs, zs, &mut zsd);
        Self { vs, xsd, zsd }
    }

    #[allow(clippy::too_many_arguments)]
    fn exsrc(
        &self,
        nx: usize,
        dx: f32,
        fx: f32,
        nz: usize,
        dz: f32,
        fz: f32,
        dt: f32,
        t: f32,
        fmax: f32,
        pwt: i32,
        mono: i32,
        s: &mut Grid,
    ) {
        for ix in 0..nx {
            for iz in 0..nz {
                s[ix][iz] = (0.0 * dt as f64) as f32;
            }
        }
        let fpeak = if mono == 0 {
            (0.5 * fmax as f64) as f32
        } else {
            fmax
        };
        let tdelay = (1.0 / fpeak as f64) as f32;
        if t as f64 > 2.0 * tdelay as f64 && mono == 0 {
            return;
        }
        let ts = ricker(t - tdelay, fpeak, mono);
        let mut v = self.vs[0];
        let mut distprev = 0.0f32;
        loop {
            let mut xv = [0.0];
            let mut zv = [0.0];
            let mut dxdv = [0.0];
            let mut dzdv = [0.0];
            intcub(0, &self.vs, &self.xsd, &[v], &mut xv);
            intcub(0, &self.vs, &self.zsd, &[v], &mut zv);
            intcub(1, &self.vs, &self.xsd, &[v], &mut dxdv);
            intcub(1, &self.vs, &self.zsd, &[v], &mut dzdv);
            let (dxdv, dzdv) = (dxdv[0], dzdv[0]);
            let mut dv = if dxdv == 0.0 {
                dz / abs(dzdv)
            } else if dzdv == 0.0 {
                dx / abs(dxdv)
            } else {
                min(dz / abs(dzdv), dx / abs(dxdv))
            };
            let end = self.vs[self.vs.len() - 1];
            if v + dv > end {
                dv = end - v;
            }
            let dist = (dv as f64 * ((dzdv * dzdv + dxdv * dxdv) as f64).sqrt()
                / ((dx * dx + dz * dz) as f64).sqrt()) as f32;
            let amp = ((dist + distprev) as f64 / 2.0) as f32;
            let xvn = (xv[0] - fx) / dx;
            let zvn = (zv[0] - fz) / dz;
            let izv = nint(zvn as f64);
            let pio2 = (3.1415926 / 2.0) as f32;
            let opwt = (1.0 / pwt as f32 as f64) as f32;
            for ix in 0..nx {
                for iz in max(0, izv - 3)..=min(nz as i32 - 1, izv + 3) {
                    let xn = ix as f32 - xvn;
                    let zn = iz as f32 - zvn;
                    let mut a = 1.0f32;
                    if (ix as i32) < pwt {
                        a = ((pio2 * (ix as f32 * opwt)) as f64).sin() as f32;
                    }
                    if (ix as i32) > nx as i32 - pwt {
                        a = (((ix as i32 - (nx as i32 - pwt)) as f32 * pio2 * opwt) as f64).cos()
                            as f32;
                    }
                    s[ix][iz as usize] = (s[ix][iz as usize] as f64
                        + (a * ts * amp) as f64 * ((-xn * xn - zn * zn) as f64).exp())
                        as f32;
                }
            }
            if dv == 0.0 {
                break;
            }
            distprev = dist;
            v += dv;
        }
    }
}
