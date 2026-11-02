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

//! SUSYNVXZCS - common-shot synthetic seismograms in V(X,Z) media.
//! Port of su/main/synthetics_waveforms_testpatterns/susynvxzcs.c (Revision 1.17).
//! Credits (from the C source):
//! Zhenyue Liu, 07/20/92, Center for Wave Phenomena;
//! many subroutines borrowed from Dave Hale's program SUSYNLV.
//! Trino Salinas, 07/30/96, fixed moving-spread geometry.
//! Chris Liner, 12/10/08, added cable option and fldr header word.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, max, min, nint, PI};
use oxwave_f_convolution::convolve_cwp;
use oxwave_f_modeling::*;
use oxwave_f_upweik::{eikpex, ray_theoretic_beta, ray_theoretic_sigma};
use oxwave_finite_difference::mkhdiff;
use oxwave_polar_grid::{polartorect, recttopolar};
use oxwave_tridiag::tripp;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUSYNVXZCS - SYNthetic seismograms of common shot in V(X,Z) media via\t",
    " \t\tKirchhoff-style modeling\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " susynvxzcs<vfile >outfile  nx= nz= [optional parameters]\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " <vfile        file containing velocities v[nx][nz]\t\t\t",
    " >outfile      file containing seismograms of common ofset\t\t",
    " nx=           number of x samples (2nd dimension) in velocity ",
    " nz=           number of z samples (1st dimension) in velocity ",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nt=501        \tnumber of time samples\t\t\t\t",
    " dt=0.004      \ttime sampling interval (sec)\t\t\t",
    " ft=0.0        \tfirst time (sec)\t\t\t\t",
    " fpeak=0.2/dt\t\tpeak frequency of symmetric Ricker wavelet (Hz)\t",
    " nxg=\t\t\tnumber of receivers of input traces\t\t",
    " dxg=15\t\treceiver sampling interval (m)\t\t\t",
    " fxg=0.0\t\tfirst receiver (m)\t\t\t\t",
    " nxd=5         \tskipped number of receivers\t\t\t",
    " nxs=1\t\t\tnumber of offsets\t\t\t\t",
    " dxs=50\t\tshot sampling interval (m)\t\t\t",
    " fxs=0.0\t\tfirst shot (m)\t\t\t\t",
    " dx=50         \tx sampling interval (m)\t\t\t\t",
    " fx=0.         \tfirst x sample (m)\t\t\t\t",
    " dz=50         \tz sampling interval (m)\t\t\t\t",
    " nxb=nx/2    \tband width centered at midpoint (see note)\t",
    " nxc=0         hozizontal range in which velocity is changed\t",
    " nzc=0         vertical range in which velocity is changed\t",
    " pert=0        =1 calculate time correction from v_p[nx][nz]\t",
    " vpfile        file containing slowness perturbation array v_p[nx][nz]\t",
    " ref=\"1:1,2;4,2\"\treflector(s):  \"amplitude:x1,z1;x2,z2;x3,z3;...\"",
    " smooth=0\t\t=1 for smooth (piecewise cubic spline) reflectors",
    " ls=0\t\t\t=1 for line source; =0 for point source\t\t",
    " tmin=10.0*dt\t\tminimum time of interest (sec)\t\t\t",
    " ndpfz=5\t\tnumber of diffractors per Fresnel zone\t\t",
    " cable=1\t\troll reciever spread with shot\t\t\t",
    " \t\t\t=0 static reciever spread\t\t\t",
    " verbose=0\t\t=1 to print some useful information\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This algorithm is based on formula (58) in Geo. Pros. 34, 686-703,\t",
    " by N. Bleistein.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Traveltime and amplitude are calculated by finite difference which\t",
    " is done only in one of every NXD receivers; in skipped receivers, \t",
    " interpolation is used to calculate traveltime and amplitude.\t\t",
    " For each receiver, traveltime and amplitude are calculated in the \t",
    " horizontal range of (xg-nxb*dx, xg+nxb*dx). Velocity is changed by \t",
    " constant extropolation in two upper trianglar corners whose width is \t",
    " nxc*dx and height is nzc*dz.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Eikonal equation will fail to solve if there is a polar turned ray.\t",
    " In this case, the program shows the related geometric information. \t",
    " There are three ways to remove the turned rays: smoothing velocity, \t",
    " reducing nxb, and increaing nxc and nzc (if the turned ray occurs  \t",
    " in shallow areas). To prevent traveltime distortion from an over-\t",
    " smoothed velocity, traveltime is corrected based on the slowness \t",
    " perturbation.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " More than one ref (reflector) may be specified.\t\t\t",
    " Note that reflectors are encoded as quoted strings, with an optional\t",
    " reflector amplitude: preceding the x,z coordinates of each reflector.\t",
    " Default amplitude is 1.0 if amplitude: part of the string is omitted.\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    /* Get required parameters. */
    let mut nx = 0i32;
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx!\n");
    }
    let mut nz = 0i32;
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!\n");
    }

    /* Get optional parameters, retaining C evaluation order and precision. */
    let mut nt = 0i32;
    if getparint("nt", &mut nt) == 0 {
        nt = 501;
    }
    CHECK_NT("nt", nt);
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = 0.004;
    }
    let mut ft = 0.0f32;
    if getparfloat("ft", &mut ft) == 0 {
        ft = 0.0;
    }
    let mut fpeak = 0.0f32;
    if getparfloat("fpeak", &mut fpeak) == 0 {
        fpeak = (0.2 / dt as f64) as f32;
    }
    let mut nxg = 0i32;
    if getparint("nxg", &mut nxg) == 0 {
        nxg = 101;
    }
    let mut dxg = 0.0f32;
    if getparfloat("dxg", &mut dxg) == 0 {
        dxg = 15.0;
    }
    let mut fxg = 0.0f32;
    if getparfloat("fxg", &mut fxg) == 0 {
        fxg = 0.0;
    }
    let mut nxd = 0i32;
    if getparint("nxd", &mut nxd) == 0 {
        nxd = 5;
    }
    let mut nxs = 0i32;
    if getparint("nxs", &mut nxs) == 0 {
        nxs = 1;
    }
    let mut dxs = 0.0f32;
    if getparfloat("dxs", &mut dxs) == 0 {
        dxs = 50.0;
    }
    let mut fxs = 0.0f32;
    if getparfloat("fxs", &mut fxs) == 0 {
        fxs = 0.0;
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
    nxc = if nxc < nxb { nxc } else { nxb };
    let mut nzc = 0i32;
    if getparint("nzc", &mut nzc) == 0 {
        nzc = 0;
    }
    nzc = if nzc < nz { nzc } else { nz };
    let mut ls = 0i32;
    if getparint("ls", &mut ls) == 0 {
        ls = 0;
    }
    let mut pert = 0i32;
    if getparint("pert", &mut pert) == 0 {
        pert = 0;
    }
    let mut tmin = 0.0f32;
    if getparfloat("tmin", &mut tmin) == 0 {
        tmin = (10.0 * dt as f64) as f32;
    }
    let mut ndpfz = 0i32;
    if getparint("ndpfz", &mut ndpfz) == 0 {
        ndpfz = 5;
    }
    let mut smooth = 0i32;
    if getparint("smooth", &mut smooth) == 0 {
        smooth = 0;
    }
    let mut cable = 0i32;
    if getparint("cable", &mut cable) == 0 {
        cable = 1;
    }
    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* Check the ranges of shots and receivers. */
    let ex = fx + (nx - 1) as f32 * dx;
    let ez = (nz - 1) as f32 * dz;
    for ixs in 0..nxs {
        let xs = fxs + ixs as f32 * dxs;
        if fx > xs || ex < xs {
            err!("shot %i lies outside of specified (x,z) grid\n", ixs);
        }
    }
    for ixs in 0..2 {
        for ixg in 0..nxg {
            let xg = if cable == 1 {
                fxg + ixg as f32 * dxg + (ixs * (nxs - 1)) as f32 * dxs
            } else {
                fxg + ixg as f32 * dxg
            };
            if fx > xg || ex < xg {
                err!("receiver %i lies outside of specified (x,z) grid\n", ixs);
            }
        }
    }

    let mut reflectors = decodeReflectors();
    if smooth == 0 {
        reflectors = breakReflectors(&reflectors);
    }
    let nr = reflectors.len() as i32;

    /* Allocate space; read velocities and slowness perturbation. */
    let mut vold = grid(nz, nx);
    let mut aold = grid(nz, nx);
    let mut told = grid(nz, nx);
    let mut bas = grid(nz, nx);
    let mut sgold = grid(nz, nx);
    let mut vpold = if pert != 0 { grid(nz, nx) } else { Vec::new() };
    let mut input = vec![0.0f32; (nx * nz) as usize];
    if fread(&mut input, &stdin()) != (nx * nz) as usize {
        // C erroneously formats stdin's FILE* as a string. On this build
        // it prints an empty filename; no invalid pointer is passed here.
        err!("cannot read %d velocities from file %s", nx * nz, "");
    }
    for ix in 0..nx as usize {
        vold[ix].copy_from_slice(&input[ix * nz as usize..(ix + 1) * nz as usize]);
    }
    if pert != 0 {
        let mut vpfile = String::new();
        if getparstring("vpfile", &mut vpfile) == 0 {
            err!("must specify vpfile!\n");
        }
        let vpfp = fopen(&vpfile, "r");
        if vpfp.as_ref().map(|fp| fread(&mut input, fp)) != Some((nx * nz) as usize) {
            err!(
                "cannot read %d slowness perturbation from file %s",
                nx * nz,
                vpfile
            );
        }
        for ix in 0..nx as usize {
            for iz in 0..nz as usize {
                vpold[ix][iz] = input[ix * nz as usize + iz] / vold[ix][iz];
            }
        }
    }

    /* Determine maximum reflector segment length, make reflectors/wavelet. */
    tmin = max(tmin, max(ft, dt));
    let dsmax = ((vold[0][0] / (2 * ndpfz) as f32) as f64 * ((tmin / fpeak) as f64).sqrt()) as f32;
    let r = makeref(dsmax, &reflectors);
    let mut ns = 0;
    for ir in 0..nr as usize {
        ns += r[ir].ns;
    }
    let w = makericker(fpeak, dt);
    if verbose != 0 {
        warn!("\nSUSYNVXZCS:");
        warn!("Total number of small reflecting segments is %d.\n", ns);
    }

    /* Set constant segy trace header parameters. */
    let mut tr = Trace::new();
    tr.trid = 1;
    tr.counit = 1;
    tr.ns = nt as u16;
    tr.dt = (1.0e6 * dt as f64) as i32 as u16;
    tr.delrt = (1.0e3 * ft as f64) as i32 as i16;
    tr.d2 = dxg;

    let nx1 = 1 + 2 * nxb;
    let mut ts = grid(nz, nx1);
    let mut angles = grid(nz, nx1); // C calls this array "as".
    let mut sgs = grid(nz, nx1);
    let mut tg = grid(nz, nx1);
    let mut ag = grid(nz, nx1);
    let mut sgg = grid(nz, nx1);
    let mut v = grid(nz, nx1);
    let mut bag = grid(nz, nx1);
    let mut vp = if pert != 0 { grid(nz, nx1) } else { Vec::new() };
    let mut tg1 = if nxd > 1 { grid(nz, nx1) } else { Vec::new() };
    let mut ag1 = if nxd > 1 { grid(nz, nx1) } else { Vec::new() };
    let mut sgg1 = if nxd > 1 { grid(nz, nx1) } else { Vec::new() };

    /* Change velocities in two triangular upper corners. */
    for iz in 0..nzc {
        // arm64 C integer division returns zero for a zero divisor.
        let corner = (iz * nxc).checked_div(nzc).unwrap_or(0);
        for ix in corner..nxc {
            vold[(nxc - ix - 1) as usize][iz as usize] = vold[(nxc - corner) as usize][iz as usize];
            vold[(nx1 - nxc + ix) as usize][iz as usize] =
                vold[(nx1 - nxc + corner - 1) as usize][iz as usize];
        }
    }

    let mut error = EikError::default();
    let mut hd = None;
    let mut tracl = 0;
    for ixs in 0..nxs {
        let xs = fxs + ixs as f32 * dxs;
        eiktam_pert(
            xs, 0.0, nz, dz, 0.0, nx, dx, fx, &vold, &vpold, pert, &mut told, &mut aold,
            &mut sgold, &mut bas, &mut error,
        );
        if error.ierr != 0 {
            err!(
                "\tEikonal equation fails to solve at x=%g, z=%g\n\tfrom shot xs=%g .\n",
                error.x_err + xs,
                error.z_err,
                xs
            );
        }
        let mut temp = (xs - fx) / dx;
        let ixd = temp as i32;
        temp -= ixd as f32;
        // C uses integer 1 in (1-temp), so this interpolation stays f32.
        let vs0 = (1.0f32 - temp) * vold[ixd as usize][0] + temp * vold[(ixd + 1) as usize][0];
        let mut nxd1 = nxd;
        let mut ixg = 0;
        while ixg < nxg {
            let mut xg = if cable == 1 {
                fxg + ixs as f32 * dxs + ixg as f32 * dxg
            } else {
                fxg + ixg as f32 * dxg
            };
            let mut fx1 = xg - nxb as f32 * dx;
            let mut nx0 = ((fx1 - fx) / dx) as i32;
            temp = (fx1 - fx) / dx - nx0 as f32;
            for ix in 0..nx1 {
                for iz in 0..nz as usize {
                    if ix < -nx0 {
                        v[ix as usize][iz] = vold[0][iz];
                    } else if ix + nx0 > nx - 2 {
                        v[ix as usize][iz] = vold[(nx - 1) as usize][iz];
                    } else {
                        let j = (ix + nx0) as usize;
                        v[ix as usize][iz] = (vold[j][iz] as f64 * (1.0 - temp as f64)
                            + (temp * vold[j + 1][iz]) as f64)
                            as f32;
                        ts[ix as usize][iz] = (told[j][iz] as f64 * (1.0 - temp as f64)
                            + (temp * told[j + 1][iz]) as f64)
                            as f32;
                        angles[ix as usize][iz] = (aold[j][iz] as f64 * (1.0 - temp as f64)
                            + (temp * aold[j + 1][iz]) as f64)
                            as f32;
                        sgs[ix as usize][iz] = (sgold[j][iz] as f64 * (1.0 - temp as f64)
                            + (temp * sgold[j + 1][iz]) as f64)
                            as f32;
                    }
                }
            }
            if pert != 0 {
                for ix in 0..nx1 {
                    for iz in 0..nz as usize {
                        vp[ix as usize][iz] = if ix < -nx0 {
                            vpold[0][iz]
                        } else if ix + nx0 > nx - 2 {
                            vpold[(nx - 1) as usize][iz]
                        } else {
                            let j = (ix + nx0) as usize;
                            (vpold[j][iz] as f64 * (1.0 - temp as f64)
                                + (temp * vpold[j + 1][iz]) as f64)
                                as f32
                        };
                    }
                }
            }

            if ixg == 0 || nxd1 == 1 {
                /* No interpolation. */
                eiktam_pert(
                    xg, 0.0, nz, dz, 0.0, nx1, dx, fx1, &v, &vp, pert, &mut tg, &mut ag, &mut sgg,
                    &mut bag, &mut error,
                );
                if error.ierr != 0 {
                    warn!(
                        "\tEikonal equation fails to solve at x=%g, z=%g\n\tfrom receiver xg=%g .",
                        error.x_err + xg,
                        error.z_err,
                        xg
                    );
                    break;
                }
                let vg0 = (1.0f32 - temp) * vold[(nx0 + nxb) as usize][0]
                    + temp * vold[(nx0 + nxb + 1) as usize][0];
                let ex1 = min(ex, xg + nxb as f32 * dx);
                makeone(
                    &ts,
                    &angles,
                    &sgs,
                    &tg,
                    &ag,
                    &sgg,
                    ex1,
                    ez,
                    dx,
                    dz,
                    fx1,
                    vs0,
                    vg0,
                    ls,
                    &w,
                    nr,
                    &r,
                    nt,
                    dt,
                    ft,
                    &mut tr.data,
                    &mut hd,
                );
                tracl += 1;
                tr.tracl = tracl;
                tr.tracr = tracl;
                tr.fldr = ixs + 1;
                tr.tracf = ixg + 1;
                tr.sx = nint(xs as f64);
                tr.gx = nint(xg as f64);
                puttr(&mut tr);
            } else {
                /* Linear interpolation. */
                eiktam_pert(
                    xg, 0.0, nz, dz, 0.0, nx1, dx, fx1, &v, &vp, pert, &mut tg1, &mut ag1,
                    &mut sgg1, &mut bag, &mut error,
                );
                if error.ierr != 0 {
                    warn!(
                        "\tEikonal equation fails to solve at x=%g, z=%g\n\tfrom receiver xg=%g .",
                        error.x_err + xg,
                        error.z_err,
                        xg
                    );
                    break;
                }
                let vg0 = (1.0f32 - temp) * vold[(nx0 + nxb) as usize][0]
                    + temp * vold[(nx0 + nxb + 1) as usize][0];
                xg -= nxd1 as f32 * dxg;
                for i in 1..=nxd1 {
                    xg += dxg;
                    fx1 = xg - nxb as f32 * dx;
                    temp = (nxd1 - i) as f32;
                    let temp1 = (1.0 / (nxd1 - i + 1) as f64) as f32;
                    for ix in 0..nx1 as usize {
                        for iz in 0..nz as usize {
                            if i == nxd1 {
                                tg[ix][iz] = tg1[ix][iz];
                                ag[ix][iz] = ag1[ix][iz];
                                sgg[ix][iz] = sgg1[ix][iz];
                            } else {
                                tg[ix][iz] = (temp * tg[ix][iz] + tg1[ix][iz]) * temp1;
                                ag[ix][iz] = (temp * ag[ix][iz] + ag1[ix][iz]) * temp1;
                                sgg[ix][iz] = (temp * sgg[ix][iz] + sgg1[ix][iz]) * temp1;
                            }
                        }
                    }
                    nx0 = ((fx1 - fx) / dx) as i32;
                    temp = (fx1 - fx) / dx - nx0 as f32;
                    for ix in max(0, -nx0)..min(nx1, nx - 1 - nx0) {
                        for iz in 0..nz as usize {
                            let j = (ix + nx0) as usize;
                            ts[ix as usize][iz] = (told[j][iz] as f64 * (1.0 - temp as f64)
                                + (temp * told[j + 1][iz]) as f64)
                                as f32;
                            sgs[ix as usize][iz] = (sgold[j][iz] as f64 * (1.0 - temp as f64)
                                + (temp * sgold[j + 1][iz]) as f64)
                                as f32;
                            angles[ix as usize][iz] = (aold[j][iz] as f64 * (1.0 - temp as f64)
                                + (temp * aold[j + 1][iz]) as f64)
                                as f32;
                        }
                    }
                    let ex1 = min(ex, xg + nxb as f32 * dx);
                    makeone(
                        &ts,
                        &angles,
                        &sgs,
                        &tg,
                        &ag,
                        &sgg,
                        ex1,
                        ez,
                        dx,
                        dz,
                        fx1,
                        vs0,
                        vg0,
                        ls,
                        &w,
                        nr,
                        &r,
                        nt,
                        dt,
                        ft,
                        &mut tr.data,
                        &mut hd,
                    );
                    tracl += 1;
                    tr.tracl = tracl;
                    tr.tracr = tracl;
                    tr.fldr = ixs + 1;
                    tr.tracf = ixg + 1;
                    tr.sx = nint(xs as f64);
                    tr.gx = nint(xg as f64);
                    puttr(&mut tr);
                }
            }
            if ixg < nxg - 1 && ixg > nxg - 1 - nxd1 {
                nxd1 = nxg - 1 - ixg;
            }
            ixg += nxd1;
        }
        warn!("\t finish shot %f", xs);
    }
    exit(CWP_Exit());
}

fn grid(nz: i32, nx: i32) -> Vec<Vec<f32>> {
    vec![vec![0.0; nz as usize]; nx as usize]
}

/// Program-local makeone from susynvxzcs.c: one synthetic seismogram.
#[allow(clippy::too_many_arguments)]
fn makeone(
    ts: &[Vec<f32>],
    angles: &[Vec<f32>],
    sgs: &[Vec<f32>],
    tg: &[Vec<f32>],
    ag: &[Vec<f32>],
    sgg: &[Vec<f32>],
    ex: f32,
    ez: f32,
    dx: f32,
    dz: f32,
    fx: f32,
    vs0: f32,
    vg0: f32,
    ls: i32,
    w: &Wavelet,
    nr: i32,
    r: &[Reflector],
    nt: i32,
    dt: f32,
    ft: f32,
    trace: &mut [f32],
    hd: &mut Option<[f32; 41]>,
) {
    const LHD: i32 = 20;
    const NHD: i32 = 1 + 2 * LHD;
    /* Make the half-derivative filter only on the first call. */
    if hd.is_none() {
        let mut filter = [0.0; NHD as usize];
        mkhdiff(dt, LHD, &mut filter);
        *hd = Some(filter);
    }
    for it in 0..nt {
        trace[it as usize] = 0.0;
    }
    for ir in 0..nr as usize {
        let ar = r[ir].a;
        let ns = r[ir].ns;
        let ds = r[ir].ds;
        let rs = &r[ir].rs;
        for is in 0..ns as usize {
            let xd = rs[is].x;
            let zd = rs[is].z;
            let cd = rs[is].c;
            let sd = rs[is].s;
            if xd < fx || xd >= ex || zd >= ez {
                continue;
            }
            let xi = (xd - fx) / dx;
            let ix = xi as i32;
            let zi = zd / dz;
            let iz = zi as i32;
            let sx = xi - ix as f32;
            let sz = zi - iz as f32;
            let ix = ix as usize;
            let iz = iz as usize;
            // Repeat the original expression for each field, including its
            // f32 products and the promotions from unsuffixed 1.0 literals.
            macro_rules! bilinear {
                ($field:expr) => {
                    ((1.0 - sz as f64)
                        * ((1.0 - sx as f64) * $field[ix][iz] as f64
                            + (sx * $field[ix + 1][iz]) as f64)
                        + sz as f64
                            * ((1.0 - sx as f64) * $field[ix][iz + 1] as f64
                                + (sx * $field[ix + 1][iz + 1]) as f64)) as f32
                };
            }
            let tsd = bilinear!(ts);
            let asd = bilinear!(angles);
            let sgsd = bilinear!(sgs);
            let tgd = bilinear!(tg);
            let agd = bilinear!(ag);
            let sggd = bilinear!(sgg);
            let ci = (cd as f64 * (asd as f64).cos() + sd as f64 * (asd as f64).sin()) as f32;
            let cr = (cd as f64 * (agd as f64).cos() + sd as f64 * (agd as f64).sin()) as f32;
            let time = tsd + tgd;
            let mut amp = if ls != 0 {
                ((vs0 * vg0 / (sgsd * sggd)) as f64).sqrt() as f32
            } else {
                ((vs0 * vg0 / (sgsd * sggd * (sgsd + sggd))) as f64).sqrt() as f32
            };
            amp *= abs(ci + cr) * ar * ds;
            // modeling.c:520-540 truncates negative sample positions toward
            // zero, which can select a row before the sinc table. Match the
            // minimally bounds-checked C before calling the shared routine.
            let tn = (time - ft) / dt;
            let jt = tn as i32;
            let jsinc = ((tn - jt as f32) * 100.0) as i32;
            if jt + 4 >= 0 && jt - 3 < nt && !(0..101).contains(&jsinc) {
                err!("addsinc: sinc index %d out of bounds", jsinc);
            }
            addsinc(time, amp, dt, ft, &mut trace[..nt as usize]);
        }
    }
    let mut temp = vec![0.0; nt as usize];
    convolve_cwp(
        -LHD,
        hd.as_ref().unwrap(),
        0,
        &trace[..nt as usize],
        0,
        &mut temp,
    );
    convolve_cwp(w.iw, &w.wv, 0, &temp, 0, &mut trace[..nt as usize]);
}

// Legacy globals in the C program. The supplied C upweik.c calls err() on
// excessive curvature rather than setting these program-owned globals.
#[derive(Default)]
struct EikError {
    ierr: i32,
    x_err: f32,
    z_err: f32,
    r_err: f32,
    ia_err: i32,
}

/// Program-local eiktam_pert from susynvxzcs.c.
/// Revisor: Zhenyue Liu, Colorado School of Mines, 7/8/92.
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
    error: &mut EikError,
) {
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
    let nr = 1 + nint((rmax / dr) as f64);
    dr = rmax / (nr - 1) as f32;
    let fr = 0.0f32;
    let (fa, ea) = if fx == 0.0 && fz == 0.0 {
        (0.0f32, (PI / 2.0) as f32)
    } else if fx < 0.0 && fz == 0.0 {
        ((-PI / 2.0) as f32, (PI / 2.0) as f32)
    } else if fx == 0.0 && fz < 0.0 {
        (0.0f32, PI as f32)
    } else {
        (-PI as f32, PI as f32)
    };
    let mut da = dr / rmax;
    let mut na = 1 + nint(((ea - fa) / da) as f64);
    da = (ea - fa) / (na - 1) as f32;
    // C compares the promoted floats to the double PI constant verbatim.
    if fa as f64 == -PI && ea as f64 == PI {
        na -= 1;
    }
    let mut time_p = grid(nz, nx);
    let mut s = grid(nz, nx);
    let mut sp = grid(na, nr);
    let mut tp = grid(na, nr);
    let mut up = grid(na, nr);
    let mut wp = grid(na, nr);
    let mut ap = grid(na, nr);
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
        &s.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        da,
        fa,
        dr,
        fr,
        &mut sp.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
    );
    let mut ss = 0.0f32;
    for ir in 0..2 {
        for ia in 0..na as usize {
            ss += sp[ir][ia];
        }
    }
    ss /= (2 * na) as f32;
    let mut r = 0.0f32;
    for ir in 0..2 {
        for ia in 0..na as usize {
            up[ir][ia] = ss;
            wp[ir][ia] = 0.0;
            tp[ir][ia] = r * ss;
        }
        r += dr;
    }
    r = dr;
    for ir in 1..(nr - 1) as usize {
        let (uc, un) = up.split_at_mut(ir + 1);
        let (wc, wn) = wp.split_at_mut(ir + 1);
        let (tc, tn) = tp.split_at_mut(ir + 1);
        eikpex(
            da,
            r,
            dr,
            &sp[ir],
            &uc[ir],
            &wc[ir],
            &tc[ir],
            &sp[ir + 1],
            &mut un[0],
            &mut wn[0],
            &mut tn[0],
        );
        if error.ierr != 0 {
            error.x_err =
                (error.r_err as f64 * ((error.ia_err as f32 * da + fa) as f64).sin()) as f32;
            error.z_err =
                (error.r_err as f64 * ((error.ia_err as f32 * da + fa) as f64).cos()) as f32;
            return;
        }
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        fr,
        &tp.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        &mut time.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
    );

    /* Propagation angles. */
    let mut a = fa;
    for ia in 0..na as usize {
        ap[0][ia] = a;
        a += da;
    }
    r = fr + dr;
    for ir in 1..nr as usize {
        a = fa;
        for ia in 0..na as usize {
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
        &ap.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        &mut angle.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
    );

    /* Sigmas. */
    r = 0.0;
    for ir in 0..2 {
        for ia in 0..na as usize {
            tp[ir][ia] = r / ss;
        }
        r += dr;
    }
    r = dr;
    for ir in 1..(nr - 1) as usize {
        let (sc, sn) = tp.split_at_mut(ir + 1);
        ray_theoretic_sigma(
            da,
            r,
            dr,
            &up[ir],
            &wp[ir],
            &sc[ir],
            &up[ir + 1],
            &wp[ir + 1],
            &mut sn[0],
        );
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        fr,
        &tp.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        &mut sig.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
    );

    /* Betas. */
    for ir in 0..2 {
        a = fa;
        for ia in 0..na as usize {
            tp[ir][ia] = a;
            a += da;
        }
    }
    r = dr;
    for ir in 1..(nr - 1) as usize {
        let (bc, bn) = tp.split_at_mut(ir + 1);
        ray_theoretic_beta(
            da,
            r,
            dr,
            &up[ir],
            &wp[ir],
            &bc[ir],
            &up[ir + 1],
            &wp[ir + 1],
            &mut bn[0],
        );
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        fr,
        &tp.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        &mut bet.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
    );

    if pert != 0 {
        recttopolar(
            dz,
            fz,
            dx,
            fx,
            &vel_p.iter().map(Vec::as_slice).collect::<Vec<_>>(),
            da,
            fa,
            dr,
            fr,
            &mut sp.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
        );
        r = 0.0;
        for ir in 0..2 {
            for ia in 0..na as usize {
                tp[ir][ia] = r / ss * sp[ir][ia];
            }
            r += dr;
        }
        r = dr;
        for ir in 1..(nr - 1) as usize {
            let (sc, sn) = tp.split_at_mut(ir + 1);
            delta_t(
                na as usize,
                da,
                r,
                dr,
                &up[ir],
                &wp[ir],
                &sc[ir],
                &sp[ir],
                &up[ir + 1],
                &wp[ir + 1],
                &mut sn[0],
                &sp[ir + 1],
            );
            r += dr;
        }
        polartorect(
            da,
            fa,
            dr,
            fr,
            &tp.iter().map(Vec::as_slice).collect::<Vec<_>>(),
            dz,
            fz,
            dx,
            fx,
            &mut time_p.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
        );
        for ix in 0..nx as usize {
            for iz in 0..nz as usize {
                time[ix][iz] += time_p[ix][iz];
            }
        }
    }
}

/// Crank-Nicolson extrapolation of delta t in polar coordinates.
/// Program-local function from susynvxzcs.c; tripp is provided by its crate.
/// Author: Zhenyue Liu, Colorado School of Mines, 07/08/92.
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
