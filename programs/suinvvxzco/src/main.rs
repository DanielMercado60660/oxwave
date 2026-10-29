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

//! SUINVXZCO - Seismic inversion of common offset data with V(X,Z).
//! Port of su/main/migration_inversion/suinvvxzco.c (Revision 1.6).
//! Author: Zhenyue Liu, 08/28/93, Colorado School of Mines.
//! Reference: Bleistein, N., Cohen, J. K., and Hagin, F., 1987,
//! Two-and-one-half dimensional Born inversion with an arbitrary reference,
//! Geophysics Vol. 52, no.1, 26-36.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, max, min, nint, PI};
use oxwave_f_upweik::*;
use oxwave_polar_grid::*;
use oxwave_tridiag::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUINVXZCO - Seismic INVersion of Common Offset data for a smooth \t",
    "             velocity function V(X,Z) plus a slowness perturbation vp(x,z)",
    " \t\t\t\t\t\t\t\t\t",
    "     suinvvxzco <infile >outfile [optional parameters] \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " vfile                  file containing velocity array v[nx][nz]\t",
    " nx=                    number of x samples (2nd dimension) in velocity",
    " nz=                    number of z samples (1st dimension) in velocity",
    " nxm=\t\t\tnumber of midpoints of input traces\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt= or from header (dt) \ttime sampling interval of input data\t",
    " offs= or from header (offset) \tsource-receiver offset\t \t",
    " dxm= or from header (d2) \tsampling interval of midpoints \t\t",
    " fxm=0\t\tfirst midpoint in input trace\t\t\t\t",
    " nxd=5\t\tskipped number of midpoints (see note)\t\t\t",
    " dx=50.0\tx sampling interval of velocity\t\t\t\t",
    " fx=0.0\tfirst x sample of velocity\t\t\t\t",
    " dz=50.0\tz sampling interval of velocity\t\t\t\t",
    " nxb=nx/2\tband centered at midpoints (see note)\t\t\t",
    " nxc=0\t\thozizontal range in which velocity is changed\t\t",
    " nzc=0\t\tvertical range in which velocity is changed\t\t",
    " fxo=0.0\tx-coordinate of first output trace \t\t\t",
    " dxo=15.0\thorizontal spacing of output trace \t\t\t",
    " nxo=101\tnumber of output traces \t\t\t\t",
    " fzo=0.0\tz-coordinate of first point in output trace \t\t",
    " dzo=15.0\tvertical spacing of output trace \t\t\t",
    " nzo=101\tnumber of points in output trace\t\t\t",
    " fmax=0.25/dt\tMaximum frequency set for operator antialiasing\t\t",
    " ang=180\tMaximum dip angle allowed in the image\t\t\t",
    " ls=0\t\t=1 for line source; =0 for point source\t\t\t",
    " pert=0\t=1 calculate time correction from v_p[nx][nz]\t\t",
    " vpfile\tfile containing slowness perturbation array v_p[nx][nz]\t",
    " verbose=1              =1 to print some useful information\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Traveltime and amplitude are calculated by finite difference which\t",
    " is done only in one of every NXD midpoints; in the skipped midpoint, \t",
    " interpolation is used to calculate traveltime and amplitude.\t\t",
    " For each midpoint, traveltime and amplitude are calculated in the \t",
    " horizontal range of (xm-nxb*dx, xm+nxb*dx). Velocity is changed by \t",
    " constant extropolation in two upper trianglar corners whose width is \t",
    " nxc*dx and height is nzc*dz.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Eikonal equation will fail to solve if there is a polar turned ray.\t",
    " In this case, the program shows the related geometric information. \t",
    " There are three ways to remove the turned ray: smoothing velocity, \t",
    " reducing nxb, and increaing nxc and nzc (if the turned ray occurs  \t",
    " in the shallow areas). To prevent traveltime distortion from a over\t",
    " smoothed velocity, traveltime is corrected based on the slowness \t",
    " perturbation.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Offsets are signed - may be positive or negative. \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

#[allow(unused_assignments)] // Retain C's coordinate updates in the interpolation loop.
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get information from the first header.
    if gettr(&mut tr) == 0 {
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

    let (mut nx, mut nz, mut nxm) = (0i32, 0i32, 0i32);
    let mut vfile = String::new();
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx!\n");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!\n");
    }
    if getparint("nxm", &mut nxm) == 0 {
        err!("must specify nxm!\n");
    }
    if getparstring("vfile", &mut vfile) == 0 {
        err!("must specify vfile!\n");
    }

    // Get optional parameters in C statement order.
    let mut fxm = 0.0f32;
    if getparfloat("fxm", &mut fxm) == 0 {
        fxm = 0.0;
    }
    let mut nxd = 0i32;
    if getparint("nxd", &mut nxd) == 0 {
        nxd = 5;
    }
    let mut dx = 0.0f32;
    if getparfloat("dx", &mut dx) == 0 {
        dx = 50.0;
    }
    let mut fx = 0.0f32;
    if getparfloat("fx", &mut fx) == 0 {
        fx = 0.0;
    }
    let mut dz = 0.0f32;
    if getparfloat("dz", &mut dz) == 0 {
        dz = 50.0;
    }
    let mut nxb = 0i32;
    if getparint("nxb", &mut nxb) == 0 {
        nxb = nx / 2;
    }
    let mut nxc = 0i32;
    if getparint("nxc", &mut nxc) == 0 {
        nxc = 0;
    }
    nxc = min(nxc, nxb);
    let mut nzc = 0i32;
    if getparint("nzc", &mut nzc) == 0 {
        nzc = 0;
    }
    nzc = min(nzc, nz);
    let mut fxo = 0.0f32;
    if getparfloat("fxo", &mut fxo) == 0 {
        fxo = 0.0;
    }
    if fxo < fx {
        err!("the first output trace is out of velocity model.\n");
    }
    let mut dxo = 0.0f32;
    if getparfloat("dxo", &mut dxo) == 0 {
        dxo = 15.0;
    }
    let mut nxo = 0i32;
    if getparint("nxo", &mut nxo) == 0 {
        nxo = 101;
    }
    let mut fzo = 0.0f32;
    if getparfloat("fzo", &mut fzo) == 0 {
        fzo = 0.0;
    }
    let mut dzo = 0.0f32;
    if getparfloat("dzo", &mut dzo) == 0 {
        dzo = 15.0;
    }
    let mut nzo = 0i32;
    if getparint("nzo", &mut nzo) == 0 {
        nzo = 101;
    }
    let mut ls = 0i32;
    if getparint("ls", &mut ls) == 0 {
        ls = 0;
    }
    let mut pert = 0i32;
    if getparint("pert", &mut pert) == 0 {
        pert = 0;
    }
    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 1;
    }
    let mut fmax = 0.0f32;
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (0.25 / dt as f64) as f32;
    }
    let mut ang = 0.0f32;
    if getparfloat("ang", &mut ang) == 0 {
        ang = 180.0;
    }
    let cosa = (ang as f64 * PI / 180.0).cos() as f32;

    let ex = fx + (nx - 1) as f32 * dx;
    let ez = (nz - 1) as f32 * dz;
    for ixm in 0..nxm {
        let xs = ((fxm + ixm as f32 * dxm) as f64 - 0.5 * offs as f64) as f32;
        let xg = xs + offs;
        if fx > xs || ex < xs || fx > xg || ex < xg {
            err!("shot or receiver lie outside of specified (x,z) grid\n");
        }
    }
    let odx = (1.0 / dx as f64) as f32;
    let odz = (1.0 / dz as f64) as f32;
    let odt = (1.0 / dt as f64) as f32;
    let smax = (0.5 / (dxm * fmax) as f64) as f32;
    // fread uses the entire contiguous C allocation, with depth fastest.
    let mut vold = vec![0.0f32; (nx * nz) as usize];
    let mut vpold = if pert != 0 {
        vec![0.0f32; (nx * nz) as usize]
    } else {
        vec![]
    };
    let mut outtrace = vec![vec![0.0f32; nzo as usize]; nxo as usize];

    // Read velocities and slowness perturbation (raw native floats).
    let fp = fopen(&vfile, "r").unwrap_or_else(|| failed_model_open());
    if fread(&mut vold, &fp) != (nx * nz) as usize {
        err!("cannot read %d velocities from file %s", nx * nz, vfile);
    }
    if pert != 0 {
        let mut vpfile = String::new();
        if getparstring("vpfile", &mut vpfile) == 0 {
            err!("must specify vpfile!\n");
        }
        let fp = fopen(&vpfile, "r").unwrap_or_else(|| failed_model_open());
        if fread(&mut vpold, &fp) != (nx * nz) as usize {
            err!(
                "cannot read %d slowness perturbation from file %s",
                nx * nz,
                vpfile
            );
        }
        for ix in 0..nx as usize {
            for iz in 0..nz as usize {
                vpold[ix * nz as usize + iz] /= vold[ix * nz as usize + iz];
            }
        }
    }
    checkpars();
    for ixo in 0..nxo as usize {
        for izo in 0..nzo as usize {
            outtrace[ixo][izo] = 0.0;
        }
    }
    let mut nx1 = 1 + 2 * nxb;
    let mut ts = vec![vec![0.0f32; nz as usize]; nx1 as usize];
    let mut angles = ts.clone(); // C's variable `as` is a Rust keyword.
    let mut sgs = ts.clone();
    let mut bas = ts.clone();
    let mut tg = ts.clone();
    let mut ag = ts.clone();
    let mut bag = ts.clone();
    let mut sgg = ts.clone();
    let mut amp = ts.clone();
    let mut v = ts.clone();
    let mut vp = if pert != 0 { ts.clone() } else { vec![] };
    let mut ts1 = if nxd > 1 { ts.clone() } else { vec![] };
    let mut tg1 = ts1.clone();
    let mut amp1 = ts1.clone();
    if abs(offs) > (2 * nxb) as f32 * dx {
        err!("\t band NXB is too small!");
    }
    let mut nxd1 = nxd;
    let nz0 = max(0.0f32, fzo / dz) as usize;
    let mut state = RayError::default();
    let mut ixm = 0i32;
    while ixm < nxm {
        let mut xm = fxm + ixm as f32 * dxm;
        let mut xs = (xm as f64 - 0.5 * offs as f64) as f32;
        let mut xg = xs + offs;
        let mut fx1 = xm - nxb as f32 * dx;
        let mut ex1 = min(ex + (nxd1 - 1) as f32 * dxm, xm + nxb as f32 * dx);
        nx1 = 1 + nint(((ex1 - fx1) / dx) as f64);
        let nx0 = ((fx1 - fx) / dx) as i32;
        let mut temp = (fx1 - fx) / dx - nx0 as f32;
        // Transpose velocity such that the first row is at fx1.
        for ix in 0..nx1 as usize {
            if (ix as i32) < -nx0 {
                for iz in 0..nz as usize {
                    v[ix][iz] = vold[iz];
                }
            } else if ix as i32 + nx0 > nx - 2 {
                for iz in 0..nz as usize {
                    v[ix][iz] = vold[(nx - 1) as usize * nz as usize + iz];
                }
            } else {
                for iz in 0..nz as usize {
                    v[ix][iz] = (vold[(ix as i32 + nx0) as usize * nz as usize + iz] as f64
                        * (1.0 - temp as f64)
                        + (temp * vold[(ix as i32 + nx0 + 1) as usize * nz as usize + iz]) as f64)
                        as f32;
                }
            }
        }
        if pert != 0 {
            for ix in 0..nx1 as usize {
                if (ix as i32) < -nx0 {
                    for iz in 0..nz as usize {
                        vp[ix][iz] = vpold[iz];
                    }
                } else if ix as i32 + nx0 > nx - 2 {
                    for iz in 0..nz as usize {
                        vp[ix][iz] = vpold[(nx - 1) as usize * nz as usize + iz];
                    }
                } else {
                    for iz in 0..nz as usize {
                        vp[ix][iz] = (vpold[(ix as i32 + nx0) as usize * nz as usize + iz] as f64
                            * (1.0 - temp as f64)
                            + (temp * vpold[(ix as i32 + nx0 + 1) as usize * nz as usize + iz])
                                as f64) as f32;
                    }
                }
            }
        }
        // Change velocity within triangular upper corners to avoid polar-turned rays.
        for iz in 0..nzc {
            for ix in iz * nxc / nzc..nxc {
                v[(nxc - ix - 1) as usize][iz as usize] =
                    v[(nxc - iz * nxc / nzc) as usize][iz as usize];
                v[(nx1 - nxc + ix) as usize][iz as usize] =
                    v[(nx1 - nxc + iz * nxc / nzc - 1) as usize][iz as usize];
            }
        }
        if ixm == 0 || nxd1 == 1 {
            // No interpolation for midpoints.
            eiktam_pert(
                xs,
                0.0,
                nz,
                dz,
                0.0,
                nx1,
                dx,
                fx1,
                &v,
                &vp,
                pert,
                &mut ts,
                &mut angles,
                &mut sgs,
                &mut bas,
                &mut state,
            );
            if state.ierr != 0 {
                fprintf!(
                    stderr(),
                    "\tEkonal equation fails to solve at x=%g, z=%g\n\tfrom shot xs=%g .\n",
                    state.x_err + xs,
                    state.z_err,
                    xs
                );
                break;
            }
            eiktam_pert(
                xg, 0.0, nz, dz, 0.0, nx1, dx, fx1, &v, &vp, pert, &mut tg, &mut ag, &mut sgg,
                &mut bag, &mut state,
            );
            if state.ierr != 0 {
                fprintf!(
                    stderr(),
                    "\tEkonal equation fails to solve at x=%g, z=%g\n\tfrom receiver xg=%g .\n",
                    state.x_err + xg,
                    state.z_err,
                    xg
                );
                break;
            }
            temp = (xs - fx) / dx;
            let ixd = temp as usize;
            temp -= ixd as f32;
            let vs0 =
                (1.0f32 - temp) * vold[ixd * nz as usize] + temp * vold[(ixd + 1) * nz as usize];
            temp = (xg - fx) / dx;
            let ixd = temp as usize;
            temp -= ixd as f32;
            let vg0 =
                (1.0f32 - temp) * vold[ixd * nz as usize] + temp * vold[(ixd + 1) * nz as usize];
            // Calculate amplitude.
            for ix in 0..nx1 as usize {
                amp[ix][nz0] = 0.0;
                for iz in nz0 + 1..nz as usize {
                    if abs((bas[ix][iz] as f64).sin() / vs0 as f64
                        + (bag[ix][iz] as f64).sin() / vg0 as f64)
                        >= smax as f64
                    {
                        amp[ix][iz] = 0.0;
                    } else if (((angles[ix][iz] + ag[ix][iz]) as f64) * 0.5).cos() < cosa as f64 {
                        amp[ix][iz] = 0.0;
                    } else {
                        temp = ((vs0 * sgg[ix][iz] / (vg0 * sgs[ix][iz])) as f64).sqrt() as f32;
                        amp[ix][iz] = (((bas[ix][iz] as f64).cos() * temp as f64 / vs0 as f64
                            + (bag[ix][iz] as f64).cos() / (temp * vg0) as f64)
                            * abs(((angles[ix][iz] - ag[ix][iz]) as f64 * 0.5).cos()))
                            as f32;
                        if ls == 0 {
                            amp[ix][iz] = (amp[ix][iz] as f64
                                * ((sgs[ix][iz] + sgg[ix][iz]) as f64).sqrt())
                                as f32;
                        }
                    }
                }
            }
            // Loop over output points.
            let mut xo = fxo;
            for ixo in 0..nxo as usize {
                let mut zo = fzo + dzo;
                for izo in 1..nzo as usize {
                    if xo - fx1 < 0.0 || xo - ex1 >= 0.0 || zo - ez >= 0.0 {
                        zo += dzo;
                        continue;
                    }
                    let xi = (xo - fx1) * odx;
                    let ix = xi as usize;
                    let zi = zo * odz;
                    let iz = zi as usize;
                    let sx = xi - ix as f32;
                    let sz = zi - iz as f32;
                    let tsd = ((1.0 - sz as f64)
                        * ((1.0 - sx as f64) * ts[ix][iz] as f64 + (sx * ts[ix + 1][iz]) as f64)
                        + sz as f64
                            * ((1.0 - sx as f64) * ts[ix][iz + 1] as f64
                                + (sx * ts[ix + 1][iz + 1]) as f64))
                        as f32;
                    let tgd = ((1.0 - sz as f64)
                        * ((1.0 - sx as f64) * tg[ix][iz] as f64 + (sx * tg[ix + 1][iz]) as f64)
                        + sz as f64
                            * ((1.0 - sx as f64) * tg[ix][iz + 1] as f64
                                + (sx * tg[ix + 1][iz + 1]) as f64))
                        as f32;
                    let mut samp = (tsd + tgd) * odt;
                    let nsamp = samp as i32;
                    if nsamp > nt - 2 {
                        zo += dzo;
                        continue;
                    }
                    samp -= nsamp as f32;
                    let ampd = ((1.0 - sz as f64)
                        * ((1.0 - sx as f64) * amp[ix][iz] as f64 + (sx * amp[ix + 1][iz]) as f64)
                        + sz as f64
                            * ((1.0 - sx as f64) * amp[ix][iz + 1] as f64
                                + (sx * amp[ix + 1][iz + 1]) as f64))
                        as f32;
                    outtrace[ixo][izo] = (outtrace[ixo][izo] as f64
                        + ampd as f64
                            * ((samp * tr.data[(nsamp + 1) as usize]) as f64
                                + (1.0 - samp as f64) * tr.data[nsamp as usize] as f64))
                        as f32;
                    zo += dzo;
                }
                xo += dxo;
            }
            if ixm < nxm - 1 {
                gettr(&mut tr);
            }
        } else {
            // Linear interpolation for midpoints.
            eiktam_pert(
                xs,
                0.0,
                nz,
                dz,
                0.0,
                nx1,
                dx,
                fx1,
                &v,
                &vp,
                pert,
                &mut ts1,
                &mut angles,
                &mut sgs,
                &mut bas,
                &mut state,
            );
            if state.ierr != 0 {
                fprintf!(
                    stderr(),
                    "\tEkonal equation fails to solve at x=%g, z=%g\n\tfrom shot xs=%g .\n",
                    state.x_err + xs,
                    state.z_err,
                    xs
                );
                break;
            }
            eiktam_pert(
                xg, 0.0, nz, dz, 0.0, nx1, dx, fx1, &v, &vp, pert, &mut tg1, &mut ag, &mut sgg,
                &mut bag, &mut state,
            );
            if state.ierr != 0 {
                fprintf!(
                    stderr(),
                    "\tEkonal equation fails to solve at x=%g, z=%g\n\tfrom receiver xg=%g .\n",
                    state.x_err + xg,
                    state.z_err,
                    xg
                );
                break;
            }
            temp = (xs - fx) / dx;
            let ixd = temp as usize;
            temp -= ixd as f32;
            let vs0 =
                (1.0f32 - temp) * vold[ixd * nz as usize] + temp * vold[(ixd + 1) * nz as usize];
            temp = (xg - fx) / dx;
            let ixd = temp as usize;
            temp -= ixd as f32;
            let vg0 =
                (1.0f32 - temp) * vold[ixd * nz as usize] + temp * vold[(ixd + 1) * nz as usize];
            // Calculate amplitude.
            for ix in 0..nx1 as usize {
                amp1[ix][nz0] = 0.0;
                for iz in nz0 + 1..nz as usize {
                    if abs((bas[ix][iz] as f64).sin() / vs0 as f64
                        + (bag[ix][iz] as f64).sin() / vg0 as f64)
                        >= smax as f64
                    {
                        amp1[ix][iz] = 0.0;
                    } else if (((angles[ix][iz] + ag[ix][iz]) as f64) * 0.5).cos() < cosa as f64 {
                        amp1[ix][iz] = 0.0;
                    } else {
                        temp = ((vs0 * sgg[ix][iz] / (vg0 * sgs[ix][iz])) as f64).sqrt() as f32;
                        amp1[ix][iz] = (((bas[ix][iz] as f64).cos() * temp as f64 / vs0 as f64
                            + (bag[ix][iz] as f64).cos() / (temp * vg0) as f64)
                            * abs(((angles[ix][iz] - ag[ix][iz]) as f64 * 0.5).cos()))
                            as f32;
                        if ls == 0 {
                            amp1[ix][iz] = (amp1[ix][iz] as f64
                                * ((sgs[ix][iz] + sgg[ix][iz]) as f64).sqrt())
                                as f32;
                        }
                    }
                }
            }
            // Interpolate quantities between two midpoints.
            xm -= nxd1 as f32 * dxm;
            for i in 1..=nxd1 {
                xm += dxm;
                xs = (xm as f64 - 0.5 * offs as f64) as f32;
                xg = xs + offs;
                fx1 = xm - nxb as f32 * dx;
                ex1 = min(ex + (nxd1 - 1) as f32 * dxm, xm + nxb as f32 * dx);
                nx1 = 1 + nint(((ex1 - fx1) / dx) as f64);
                temp = (nxd1 - i) as f32;
                let temp1 = (1.0 / (nxd1 - i + 1) as f64) as f32;
                for ix in 0..nx1 as usize {
                    for iz in nz0..nz as usize {
                        ts[ix][iz] = (temp * ts[ix][iz] + ts1[ix][iz]) * temp1;
                        tg[ix][iz] = (temp * tg[ix][iz] + tg1[ix][iz]) * temp1;
                        amp[ix][iz] = (temp * amp[ix][iz] + amp1[ix][iz]) * temp1;
                    }
                }
                ex1 = min(ex, xm + nxb as f32 * dx);
                // Loop over output points.
                let mut xo = fxo;
                for ixo in 0..nxo as usize {
                    let mut zo = fzo + dzo;
                    for izo in 1..nzo as usize {
                        if xo - fx1 < 0.0 || xo - ex1 >= 0.0 || zo - ez >= 0.0 {
                            zo += dzo;
                            continue;
                        }
                        let xi = (xo - fx1) * odx;
                        let ix = xi as usize;
                        let zi = zo * odz;
                        let iz = zi as usize;
                        let sx = xi - ix as f32;
                        let sz = zi - iz as f32;
                        let tsd = ((1.0 - sz as f64)
                            * ((1.0 - sx as f64) * ts[ix][iz] as f64
                                + (sx * ts[ix + 1][iz]) as f64)
                            + sz as f64
                                * ((1.0 - sx as f64) * ts[ix][iz + 1] as f64
                                    + (sx * ts[ix + 1][iz + 1]) as f64))
                            as f32;
                        let tgd = ((1.0 - sz as f64)
                            * ((1.0 - sx as f64) * tg[ix][iz] as f64
                                + (sx * tg[ix + 1][iz]) as f64)
                            + sz as f64
                                * ((1.0 - sx as f64) * tg[ix][iz + 1] as f64
                                    + (sx * tg[ix + 1][iz + 1]) as f64))
                            as f32;
                        let mut samp = (tsd + tgd) * odt;
                        let nsamp = samp as i32;
                        if nsamp > nt - 2 {
                            zo += dzo;
                            continue;
                        }
                        samp -= nsamp as f32;
                        let ampd = ((1.0 - sz as f64)
                            * ((1.0 - sx as f64) * amp[ix][iz] as f64
                                + (sx * amp[ix + 1][iz]) as f64)
                            + sz as f64
                                * ((1.0 - sx as f64) * amp[ix][iz + 1] as f64
                                    + (sx * amp[ix + 1][iz + 1]) as f64))
                            as f32;
                        outtrace[ixo][izo] = (outtrace[ixo][izo] as f64
                            + ampd as f64
                                * ((samp * tr.data[(nsamp + 1) as usize]) as f64
                                    + (1.0 - samp as f64) * tr.data[nsamp as usize] as f64))
                            as f32;
                        zo += dzo;
                    }
                    xo += dxo;
                }
                if ixm - nxd1 + i < nxm - 1 {
                    gettr(&mut tr);
                }
            }
        }
        if ixm < nxm - 1 && ixm > nxm - 1 - nxd1 {
            nxd1 = nxm - 1 - ixm;
        }
        if verbose != 0 {
            fprintf!(stderr(), "\tfinish midpoint %i\n", ixm + 1);
        }
        if state.ierr != 0 {
            break;
        }
        ixm += nxd1;
    }
    let temp = (4.0 / PI.sqrt() * dxm as f64) as f32;
    for ixo in 0..nxo as usize {
        for izo in 0..nzo as usize {
            outtrace[ixo][izo] *= temp;
        }
        tro.offset = offs as i32;
        tro.tracl = 1 + ixo as i32;
        tro.tracr = 1 + ixo as i32;
        tro.ns = nzo as u16;
        tro.d1 = dzo;
        tro.ns = nzo as u16;
        tro.f1 = fzo;
        tro.f2 = fxo;
        tro.d2 = dxo;
        tro.cdp = ixo as i32;
        tro.trid = 200;
        tro.data[..nzo as usize].copy_from_slice(&outtrace[ixo]);
        puttr(&mut tro);
    }
    exit(CWP_Exit());
}

// C passes a null stream to fread after an unsuccessful fopen and dies with
// SIGSEGV. Preserve the pipeline status without linking to C's raise(). Exact
// signal termination (and the shell's diagnostic) needs a core signal helper.
fn failed_model_open() -> ! {
    exit(139)
}

#[derive(Default)]
struct RayError {
    ierr: i32,
    x_err: f32,
    z_err: f32,
    r_err: f32,
    ia_err: i32,
}

// Difference equation extrapolation of delta t in polar coordinates.
// Author: Zhenyue Liu, Colorado School of Mines, 07/8/92.
#[allow(clippy::too_many_arguments)]
fn delta_t(
    na: usize,
    da: f32,
    r: f32,
    dr: f32,
    uc: &[f32],
    wc: &[f32],
    sc: &[f32],
    bc: &[f32],
    un: &[f32],
    wn: &[f32],
    sn: &mut [f32],
    bn: &[f32],
) {
    let mut d = vec![0.0f32; na - 2];
    let mut b = vec![0.0f32; na - 2];
    let mut c = vec![0.0f32; na - 2];
    let mut e = vec![0.0f32; na - 2];
    let r1 = r + dr;
    for i in 0..na - 2 {
        d[i] = ((uc[i + 1] + un[i + 1]) as f64 / (2.0 * dr as f64)) as f32;
        e[i] = ((wn[i + 1] / (r1 * r1) + wc[i + 1] / (r * r)) as f64 / (8.0 * da as f64)) as f32;
        b[i] = ((bc[i + 1] + bn[i + 1]) as f64 * 0.5 - ((sc[i + 2] - sc[i]) * e[i]) as f64
            + (d[i] * sc[i + 1]) as f64) as f32;
        c[i] = -e[i];
    }
    d[0] += c[0];
    d[na - 3] += e[na - 3];
    tripp(&mut d, &mut e, &mut c, &mut b);
    for i in 0..na - 2 {
        sn[i + 1] = b[i];
    }
    sn[0] = sn[1];
    sn[na - 1] = sn[na - 2];
}

// Compute traveltimes, propagation angle, sigma and beta, including perturbation.
// Revisor: Zhenyue Liu, Colorado School of Mines, 7/8/92.
#[allow(clippy::too_many_arguments)]
fn eiktam_pert(
    xs: f32,
    zs: f32,
    nz: i32,
    dz: f32,
    mut fz: f32,
    nx: i32,
    dx: f32,
    mut fx: f32,
    vel: &[Vec<f32>],
    vel_p: &[Vec<f32>],
    pert: i32,
    time: &mut [Vec<f32>],
    angle: &mut [Vec<f32>],
    sig: &mut [Vec<f32>],
    bet: &mut [Vec<f32>],
    state: &mut RayError,
) {
    // Shift coordinates so source is at (x=0,z=0).
    fx -= xs;
    fz -= zs;
    let ex = fx + (nx - 1) as f32 * dx;
    let ez = fz + (nz - 1) as f32 * dz;
    let mut rmaxs = fx * fx + fz * fz;
    rmaxs = max(rmaxs, fx * fx + ez * ez);
    rmaxs = max(rmaxs, ex * ex + ez * ez);
    rmaxs = max(rmaxs, ex * ex + fz * fz);
    let rmax = (rmaxs as f64).sqrt() as f32;
    let mut dr = min(abs(dx), abs(dz));
    let nr = (1 + nint((rmax / dr) as f64)) as usize;
    dr = rmax / (nr - 1) as f32;
    let fr = 0.0f32;
    let (fa, ea): (f32, f32) = if fx == 0.0 && fz == 0.0 {
        (0.0, (PI / 2.0) as f32)
    } else if fx < 0.0 && fz == 0.0 {
        ((-PI / 2.0) as f32, (PI / 2.0) as f32)
    } else if fx == 0.0 && fz < 0.0 {
        (0.0, PI as f32)
    } else {
        (-PI as f32, PI as f32)
    };
    let mut da = dr / rmax;
    let mut na = (1 + nint(((ea - fa) / da) as f64)) as usize;
    da = (ea - fa) / (na - 1) as f32;
    if fa as f64 == -PI && ea as f64 == PI {
        na -= 1;
    }

    let mut s = vec![vec![0.0f32; nz as usize]; nx as usize];
    let mut time_p = s.clone();
    let mut sp = vec![vec![0.0f32; na]; nr];
    let mut tp = sp.clone();
    let mut up = sp.clone();
    let mut wp = sp.clone();
    let mut ap = sp.clone();
    for ix in 0..nx as usize {
        for iz in 0..nz as usize {
            s[ix][iz] = (1.0 / vel[ix][iz] as f64) as f32;
        }
    }
    recttopolar(
        dz,
        fz,
        dx,
        fx,
        &s.iter().map(|v| &v[..]).collect::<Vec<_>>(),
        da,
        fa,
        dr,
        fr,
        &mut sp.iter_mut().map(|v| &mut v[..]).collect::<Vec<_>>(),
    );
    let mut ss = 0.0f32;
    for ir in 0..2 {
        for ia in 0..na {
            ss += sp[ir][ia];
        }
    }
    ss /= (2 * na) as f32;
    let mut r = 0.0f32;
    for ir in 0..2 {
        for ia in 0..na {
            up[ir][ia] = ss;
            wp[ir][ia] = 0.0;
            tp[ir][ia] = r * ss;
        }
        r += dr;
    }
    r = dr;
    for ir in 1..nr - 1 {
        let (up0, up1) = up.split_at_mut(ir + 1);
        let (wp0, wp1) = wp.split_at_mut(ir + 1);
        let (tp0, tp1) = tp.split_at_mut(ir + 1);
        eikpex(
            da,
            r,
            dr,
            &sp[ir],
            &up0[ir],
            &wp0[ir],
            &tp0[ir],
            &sp[ir + 1],
            &mut up1[0],
            &mut wp1[0],
            &mut tp1[0],
        );
        if state.ierr != 0 {
            state.x_err =
                (state.r_err as f64 * ((state.ia_err as f32 * da + fa) as f64).sin()) as f32;
            state.z_err =
                (state.r_err as f64 * ((state.ia_err as f32 * da + fa) as f64).cos()) as f32;
            return;
        }
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        fr,
        &tp.iter().map(|v| &v[..]).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        &mut time[..nx as usize]
            .iter_mut()
            .map(|v| &mut v[..])
            .collect::<Vec<_>>(),
    );
    let mut a = fa;
    for ia in 0..na {
        ap[0][ia] = a;
        a += da;
    }
    r = fr + dr;
    for ir in 1..nr {
        a = fa;
        for ia in 0..na {
            ap[ir][ia] = (a as f64 + ((wp[ir][ia] / (sp[ir][ia] * r)) as f64).asin()) as f32;
            a += da;
        }
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        fr,
        &ap.iter().map(|v| &v[..]).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        &mut angle[..nx as usize]
            .iter_mut()
            .map(|v| &mut v[..])
            .collect::<Vec<_>>(),
    );
    r = 0.0;
    for ir in 0..2 {
        for ia in 0..na {
            tp[ir][ia] = r / ss;
        }
        r += dr;
    }
    r = dr;
    for ir in 1..nr - 1 {
        let (tp0, tp1) = tp.split_at_mut(ir + 1);
        ray_theoretic_sigma(
            da,
            r,
            dr,
            &up[ir],
            &wp[ir],
            &tp0[ir],
            &up[ir + 1],
            &wp[ir + 1],
            &mut tp1[0],
        );
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        fr,
        &tp.iter().map(|v| &v[..]).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        &mut sig[..nx as usize]
            .iter_mut()
            .map(|v| &mut v[..])
            .collect::<Vec<_>>(),
    );
    for ir in 0..2 {
        a = fa;
        for ia in 0..na {
            tp[ir][ia] = a;
            a += da;
        }
    }
    r = dr;
    for ir in 1..nr - 1 {
        let (tp0, tp1) = tp.split_at_mut(ir + 1);
        ray_theoretic_beta(
            da,
            r,
            dr,
            &up[ir],
            &wp[ir],
            &tp0[ir],
            &up[ir + 1],
            &wp[ir + 1],
            &mut tp1[0],
        );
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        fr,
        &tp.iter().map(|v| &v[..]).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        &mut bet[..nx as usize]
            .iter_mut()
            .map(|v| &mut v[..])
            .collect::<Vec<_>>(),
    );
    if pert != 0 {
        recttopolar(
            dz,
            fz,
            dx,
            fx,
            &vel_p[..nx as usize]
                .iter()
                .map(|v| &v[..])
                .collect::<Vec<_>>(),
            da,
            fa,
            dr,
            fr,
            &mut sp.iter_mut().map(|v| &mut v[..]).collect::<Vec<_>>(),
        );
        r = 0.0;
        for ir in 0..2 {
            for ia in 0..na {
                tp[ir][ia] = r / ss * sp[ir][ia];
            }
            r += dr;
        }
        r = dr;
        for ir in 1..nr - 1 {
            let (tp0, tp1) = tp.split_at_mut(ir + 1);
            delta_t(
                na,
                da,
                r,
                dr,
                &up[ir],
                &wp[ir],
                &tp0[ir],
                &sp[ir],
                &up[ir + 1],
                &wp[ir + 1],
                &mut tp1[0],
                &sp[ir + 1],
            );
            r += dr;
        }
        polartorect(
            da,
            fa,
            dr,
            fr,
            &tp.iter().map(|v| &v[..]).collect::<Vec<_>>(),
            dz,
            fz,
            dx,
            fx,
            &mut time_p.iter_mut().map(|v| &mut v[..]).collect::<Vec<_>>(),
        );
        for ix in 0..nx as usize {
            for iz in 0..nz as usize {
                time[ix][iz] += time_p[ix][iz];
            }
        }
    }
}
