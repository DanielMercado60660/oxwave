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

//! RAYT2DAN - P- and SV-wave raytracing in 2D anisotropic media.
//! Port of par/main/ray_theory/rayt2dan.c. Credits: Debashish Sarkar.

#![allow(
    non_snake_case,
    unused_mut,
    unused_variables,
    unused_assignments,
    unused_parens
)]
use oxwave_core::prelude::*;
mod ray;
mod vel;
use ray::{make_ray, Ray};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " RAYT2DAN -- P- and SV-wave raytracing in 2D anisotropic media\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " rayt2dan > ttime parameterfiles= nt= nx= nz= [optional parameters]\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " VP0file=\t\t name of file containing VP0(x,z)\t\t",
    " nt=                    number of time samples\tfor each ray\t\t",
    " nx=                    number of samples (x) for the parameter fields\t",
    " nz=                    number of samples (z) for the parameter fields\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " SV=0\t\t\t for P-waves, SV=1 for Shear waves\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Parameters defining the velocity field\t\t\t\t",
    " dt=0.008                 time sampling interval\t\t\t",
    " fx=0                     first lateral sample (x) in parameter field  ",
    " fz=0                     first lateral sample (z) in parameter field  ",
    " dx=100.0                 sample spacing (x) for the parameter fields\t",
    " dz=100.0                 sample spacing (z) for the parameter fields\t",
    " \t\t\t\t\t\t\t\t\t",
    " Parameters defining the takeoff angle of a ray at a source position\t",
    " fa=-60                 first take-off angle of rays (degrees)         ",
    " na=61                  number of rays\t\t\t\t\t",
    " da=2                   increment of take-off angle\t\t\t",
    " amin=0                 minimum emergence angle; must be > -90 degrees\t",
    " amax=90                maximum emergence angle; must be < 90 degrees\t",
    " \t\t\t\t\t\t\t\t\t",
    " Parameters defining the output traveltime table\t\t\t",
    " fxo=fx                 first lateral sample in traveltime table\t",
    " nxo=nx                 number of later samples in traveltime table\t",
    " dxo=dx                 lateral interval in traveltime table\t\t",
    " fzo=fz                 first depth sample in traveltime table\t\t",
    " nzo=nz                 number of depth samples in traveltime table\t",
    " dzo=dz                 depth interval in traveltime table\t\t",
    " fac=0.01               factor to determine the radius of extrap. \t",
    " \t\t\t\t\t\t\t\t\t",
    " Parameters defining the source positions\t\t\t        ",
    " fsx=fx                 x-coordinate of first source                   ",
    " nsx=1                  number of sources                              ",
    " dsx=2*dxo              x-coordinate increment of sources              ",
    " aperx=0.5*nx*dx        ray tracing aperature in x-direction           ",
    "\t\t\t\t\t\t\t\t\t",
    " Files for general anisotropic parameters confined to a vertical plane:",
    " a1111file\t\tname of file containing a1111(x,z)\t\t",
    " a1133file          \tname of file containing a1133(x,z)\t\t",
    " VS0file\t       \tname of file containing VS0(x,z)\t\t",
    " a1113file          \tname of file containing a1113(x,z)\t\t",
    " a3313file          \tname of file containing a3313(x,z)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For transversely isotropic media Thomsen's parameters could be used:\t",
    " deltafile\t\tname of file containing delta(x,z)\t\t",
    " epsilonfile\t\tname of file containing epsilon(x,z)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " if anisotropy parameters are not given the program considers\t\t",
    " isotropic media.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut nx: i32 = 0;
    let mut nz: i32 = 0;
    let mut nt: i32 = 0;
    let mut ix: i32 = 0;
    let mut iz: i32 = 0;
    let mut ia1111: i32 = 1i32;
    let mut ia1133: i32 = 1i32;
    let mut iVS0: i32 = 1i32;
    let mut ia1113: i32 = 1i32;
    let mut ia3313: i32 = 1i32;
    let mut it: i32 = 0;
    let mut nxf: i32 = 0;
    let mut nxe: i32 = 0;
    let mut nzf: i32 = 0;
    let mut nze: i32 = 0;
    let mut ixo: i32 = 0;
    let mut izo: i32 = 0;
    let mut i: i32 = 0;
    let mut NRS: i32 = 0;
    let mut SV: i32 = 0;
    let mut idelta: i32 = 1i32;
    let mut iepsilon: i32 = 1i32;
    let mut cntt: i32 = 0;
    let mut isx: i32 = 0;
    let mut ia: i32 = 0;
    let mut nxo: i32 = 0;
    let mut nzo: i32 = 0;
    let mut nsx: i32 = 0;
    let mut na: i32 = 0;
    let mut sx: f32 = 0.0;
    let mut dx: f32 = 0.0;
    let mut sz: f32 = 0.0;
    let mut dz: f32 = 0.0;
    let mut dt: f32 = 0.0;
    let mut a: f32 = 0.0;
    let mut fxo: f32 = 0.0;
    let mut dxo: f32 = 0.0;
    let mut fzo: f32 = 0.0;
    let mut dzo: f32 = 0.0;
    let mut fx: f32 = 0.0;
    let mut fz: f32 = 0.0;
    let mut ft: f32 = 0.0;
    let mut v0: f32 = 0.0;
    let mut vd: f32 = 0.0;
    let mut vd1: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut exo: f32 = 0.0;
    let mut r2: f32 = 0.0;
    let mut r1: f32 = 0.0;
    let mut r1x: f32 = 0.0;
    let mut t1x: f32 = 0.0;
    let mut t2x: f32 = 0.0;
    let mut t2xz: f32 = 0.0;
    let mut xoc: f32 = 0.0;
    let mut zoc: f32 = 0.0;
    let mut xc: f32 = 0.0;
    let mut zc: f32 = 0.0;
    let mut terr: f32 = 0.0;
    let mut tc: f32 = 0.0;
    let mut norm2: f32 = 0.0;
    let mut curv: f32 = 0.0;
    let mut x1: f32 = 0.0;
    let mut c1: f32 = 0.0;
    let mut s1: f32 = 0.0;
    let mut z1: f32 = 0.0;
    let mut Det: f32 = 0.0;
    let mut fac: f32 = 0.0;
    let mut r2x: f32 = 0.0;
    let mut v1: f32 = 0.0;
    let mut t1: f32 = 0.0;
    let mut t2: f32 = 0.0;
    let mut a1: f32 = 0.0;
    let mut aperx: f32 = 0.0;
    let mut ezo: f32 = 0.0;
    let mut xo: f32 = 0.0;
    let mut zo: f32 = 0.0;
    let mut sd: f32 = 0.0;
    let mut fsx: f32 = 0.0;
    let mut dsx: f32 = 0.0;
    let mut amin: f32 = 0.0;
    let mut amax: f32 = 0.0;
    let mut da: f32 = 0.0;
    let mut fa: f32 = 0.0;

    /* get required parameters */
    let mut VP0file = String::new();
    if getparstring("VP0file", &mut VP0file) == 0 {
        err!("Must specify VP0file!");
    }
    if getparint("nt", &mut nt) == 0 {
        err!("Must specify nt!");
    }
    if getparint("nx", &mut nx) == 0 {
        err!("Must specify nx!");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("Must specify nz!");
    }

    /* get optional parameters */
    if getparint("SV", &mut SV) == 0 {
        SV = 0;
    }
    if getparfloat("fx", &mut fx) == 0 {
        fx = 0.0;
    }
    if getparfloat("fz", &mut fz) == 0 {
        fz = 0.0;
    }
    if getparfloat("dx", &mut dx) == 0 {
        dx = 100.0;
    }
    if getparfloat("dz", &mut dz) == 0 {
        dz = 100.0;
    }
    if getparfloat("dt", &mut dt) == 0 {
        dt = 0.008;
    }
    if getparfloat("fac", &mut fac) == 0 {
        fac = 0.01;
    }
    if getparfloat("amin", &mut amin) == 0 {
        amin = 0.0;
    }
    if getparfloat("amax", &mut amax) == 0 {
        amax = 90.0;
    }
    if getparfloat("fa", &mut fa) == 0 {
        fa = -60.0;
    }
    if getparint("na", &mut na) == 0 {
        na = 61;
    }
    if getparfloat("da", &mut da) == 0 {
        da = 2.0;
    }
    if getparfloat("fxo", &mut fxo) == 0 {
        fxo = fx;
    }
    if getparint("nxo", &mut nxo) == 0 {
        nxo = nx;
    }
    if getparfloat("dxo", &mut dxo) == 0 {
        dxo = dx;
    }
    if getparfloat("fzo", &mut fzo) == 0 {
        fzo = fz;
    }
    if getparint("nzo", &mut nzo) == 0 {
        nzo = nz;
    }
    if getparfloat("dzo", &mut dzo) == 0 {
        dzo = dz;
    }
    if getparfloat("fsx", &mut fsx) == 0 {
        fsx = fx;
    }
    if getparint("nsx", &mut nsx) == 0 {
        nsx = 1;
    }
    if getparfloat("dsx", &mut dsx) == 0 {
        dsx = 2.0f32 * dxo;
    }
    if getparfloat("aperx", &mut aperx) == 0 {
        aperx = (0.5f64 * nx as f64 * dx as f64) as f32;
    }

    /* output traveltime file */
    let mut tfile = String::new();
    let tfp = if getparstring("tfile", &mut tfile) == 0 {
        Some(stdout())
    } else {
        fopen(&tfile, "w")
    };

    /* ensure takeoff point is within model */
    if fsx < fx
        || fsx > fx + dx * (nx - 1) as f32
        || fsx + dsx * ((nsx - 1) as f32) < fx
        || fsx + dsx * ((nsx - 1) as f32) > fx + dx * (nx - 1) as f32
    {
        err!("Shot point locations exceed model dimensions");
    }

    /* allocate workspace */
    let ng = (nz * nx) as usize;
    let mut VP0 = vec![0.0f32; ng];
    let mut VS0 = vec![0.0f32; ng];
    let mut a1111 = vec![0.0f32; ng];
    let mut a3333 = vec![0.0f32; ng];
    let mut a1133 = vec![0.0f32; ng];
    let mut a1313 = vec![0.0f32; ng];
    let mut a1113 = vec![0.0f32; ng];
    let mut a3313 = vec![0.0f32; ng];
    let mut delta = vec![0.0f32; ng];
    let mut epsilon = vec![0.0f32; ng];
    // The C extrapolator reads [nt] on its final step. If malloc has
    // padding there, its fresh bytes are zero. At an exact allocation
    // boundary the read aliases the following allocation (see below).
    unsafe extern "C" {
        fn malloc_good_size(size: usize) -> usize;
    }
    let nt_bytes = nt.max(0) as usize * std::mem::size_of::<f32>();
    let adjacent = nt_bytes != 0 && unsafe { malloc_good_size(nt_bytes) } == nt_bytes;
    let mut tt = vec![0.0f32; nt.max(0) as usize + 1];
    let mut xx = vec![0.0f32; nt.max(0) as usize + 1];
    let mut zz = vec![0.0f32; nt.max(0) as usize + 1];
    let mut cc = vec![0.0f32; nt.max(0) as usize + 1];
    let mut ss = vec![0.0f32; nt.max(0) as usize + 1];
    let mut vv = vec![0.0f32; nt.max(0) as usize + 1];
    let mut dvdx1 = vec![0.0f32; nt.max(0) as usize + 1];
    let mut dvdz1 = vec![0.0f32; nt.max(0) as usize + 1];
    let mut t = vec![0.0f32; (nxo * nzo).max(0) as usize];
    let mut s = vec![0.0f32; (nxo * nzo).max(0) as usize];
    let mut N = [[0.0f32; 2]; 2];
    let mut X = [[0.0f32; 2]; 2];
    let mut Y = [[0.0f32; 2]; 2];
    let mut ray = Ray::default();
    let mut a1111file = String::new();
    ia1111 = getparstring("a1111file", &mut a1111file);
    let mut a1133file = String::new();
    ia1133 = getparstring("a1133file", &mut a1133file);
    let mut VS0file = String::new();
    iVS0 = getparstring("VS0file", &mut VS0file);
    let mut a1113file = String::new();
    ia1113 = getparstring("a1113file", &mut a1113file);
    let mut a3313file = String::new();
    ia3313 = getparstring("a3313file", &mut a3313file);
    let mut deltafile = String::new();
    idelta = getparstring("deltafile", &mut deltafile);
    let mut epsilonfile = String::new();
    iepsilon = getparstring("epsilonfile", &mut epsilonfile);

    checkpars();
    /* read required velocities */
    read_grid(&VP0file, "VP0file", &mut VP0);
    ix = 0;
    while ix < nx {
        iz = 0;
        while iz < nz {
            let k = (ix * nz + iz) as usize;
            a3333[k] = VP0[k] * VP0[k];
            iz += 1;
        }
        ix += 1;
    }
    if ia1111 != 0 {
        read_grid(&a1111file, "a1111file", &mut a1111);
    } else {
        ix = 0;
        while ix < nx {
            iz = 0;
            while iz < nz {
                let k = (ix * nz + iz) as usize;
                a1111[k] = a3333[k];
                iz += 1;
            }
            ix += 1;
        }
    }
    if ia1133 != 0 {
        read_grid(&a1133file, "a1133file", &mut a1133);
    } else {
        ix = 0;
        while ix < nx {
            iz = 0;
            while iz < nz {
                let k = (ix * nz + iz) as usize;
                a1133[k] = (0.4f64 * a3333[k] as f64) as f32;
                iz += 1;
            }
            ix += 1;
        }
    }
    if iVS0 != 0 {
        read_grid(&VS0file, "VS0file", &mut VS0);
    } else {
        ix = 0;
        while ix < nx {
            iz = 0;
            while iz < nz {
                let k = (ix * nz + iz) as usize;
                VS0[k] = (0.5477f64 * VP0[k] as f64) as f32;
                iz += 1;
            }
            ix += 1;
        }
    }
    /* compute a1313 */
    ix = 0;
    while ix < nx {
        iz = 0;
        while iz < nz {
            let k = (ix * nz + iz) as usize;
            a1313[k] = VS0[k] * VS0[k];
            iz += 1;
        }
        ix += 1;
    }
    if ia1113 != 0 {
        read_grid(&a1113file, "a1113file", &mut a1113);
    } else {
        ix = 0;
        while ix < nx {
            iz = 0;
            while iz < nz {
                let k = (ix * nz + iz) as usize;
                a1113[k] = 0.0;
                iz += 1;
            }
            ix += 1;
        }
    }
    if ia3313 != 0 {
        read_grid(&a3313file, "a3313file", &mut a3313);
    } else {
        ix = 0;
        while ix < nx {
            iz = 0;
            while iz < nz {
                let k = (ix * nz + iz) as usize;
                a3313[k] = 0.0;
                iz += 1;
            }
            ix += 1;
        }
    }

    /* if specified read Thomsen's parameters */
    if idelta == 1 || iepsilon == 1 {
        if idelta != 0 {
            read_grid(&deltafile, "deltafile", &mut delta);
        }
        if iepsilon != 0 {
            read_grid(&epsilonfile, "epsilonfile", &mut epsilon);
        }
        ix = 0;
        while ix < nx {
            iz = 0;
            while iz < nz {
                let k = (ix * nz + iz) as usize;
                a1111[k] = a3333[k] + 2.0f32 * epsilon[k] * a3333[k];
                // sqrt returns double; subtract before narrowing, as in C.
                a1133[k] = (((2.0f32 * delta[k] * a3333[k] * (a3333[k] - a1313[k])
                    + (a3333[k] - a1313[k]) * (a3333[k] - a1313[k]))
                    as f64)
                    .sqrt()
                    - a1313[k] as f64) as f32;
                iz += 1;
            }
            ix += 1;
        }
    }

    /* iterate over sources	*/
    isx = 0i32;
    while (isx < nsx) {
        sx = (fsx + ((isx as f32) * dsx));
        /* initialize time and raypath length*/
        ixo = 0i32;
        while (ixo < nxo) {
            izo = 0i32;
            while (izo < nzo) {
                i = (izo + (ixo * nzo));
                t[(i) as usize] = (99999999i32 as f32);
                s[(i) as usize] = (99999999i32 as f32);
                izo += 1;
            }
            ixo += 1;
        }
        /* iterate over takeoff angle */
        ia = 0i32;
        while (ia < na) {
            a = (fa + ((ia as f32) * da));
            printf!("sx = %f angle = %f\n", (sx as f64), (a as f64));
            /* Assume 0s as initial time */
            ft = (0i32 as f32);
            sz = (0f64 as f32);
            ray = make_ray(
                SV, sx, sz, a, nt, dt, ft, &a1111, &a3333, &a1133, &a1313, &a1113, &a3313, nx, dx,
                fx, nz, dz, fz, amax, amin,
            );
            /* store the ray */
            cntt = 0i32;
            while (cntt < nt) {
                tt[(cntt) as usize] = ray.rs[(cntt) as usize].t;
                xx[(cntt) as usize] = ray.rs[(cntt) as usize].x;
                zz[(cntt) as usize] = ray.rs[(cntt) as usize].z;
                cc[(cntt) as usize] = ray.rs[(cntt) as usize].c;
                ss[(cntt) as usize] = ray.rs[(cntt) as usize].s;
                vv[(cntt) as usize] = ray.rs[(cntt) as usize].v;
                dvdx1[(cntt) as usize] = ray.rs[(cntt) as usize].dvdx;
                dvdz1[(cntt) as usize] = ray.rs[(cntt) as usize].dvdz;
                cntt += 1;
            }
            NRS = ray.nrs;
            /* free ray */
            /* Extrapolate traveltimes around the central ray
            (Similar to Zhenyue Liu's code, rayt2d)*/
            v0 = vv[(0i32) as usize];
            r2 = (0f64 as f32);
            xc = sx;
            zc = sz;
            sd = (0i32 as f32);
            vd1 = v0;
            /* trace the auxillary ray */
            a1 = (((a as f64) + 0.1f64) as f32);
            ray = make_ray(
                SV, sx, sz, a1, NRS, dt, ft, &a1111, &a3333, &a1133, &a1313, &a1113, &a3313, nx,
                dx, fx, nz, dz, fz, amax, amin,
            );
            it = 1i32;
            while (it < NRS) {
                if adjacent && it == nt - 1 {
                    // C allocation order: tt,xx,zz,t,s,cc,ss,vv,dvdx1,dvdz1.
                    let end = nt as usize;
                    xx[end] = zz[0];
                    zz[end] = t.first().copied().unwrap_or(0.0);
                    cc[end] = ss[0];
                    ss[end] = vv[0];
                    vv[end] = dvdx1[0];
                }

                xo = xx[(it) as usize];
                zo = zz[(it) as usize];
                vd = vv[(it) as usize];
                sd = (((sd as f64) + ((0.5f64 * (dt as f64)) * ((vd + vd1) as f64))) as f32);
                vd1 = vd;
                curv = ((((-((dvdx1[(it) as usize] * cc[(it) as usize])
                    - (dvdz1[(it) as usize] * ss[(it) as usize])))
                    as f64)
                    .abs()
                    / (vv[(it) as usize] as f64)) as f32);
                /* Check if the point is in the previous circle */
                if ((r2 > (((xc - xo) * (xc - xo)) + ((zc - zo) * (zc - zo))))
                    && (it != (NRS - 1i32)))
                {
                    it += 1;
                    continue;
                }
                exo = (fxo + (((nxo - 1i32) as f32) * dxo));
                ezo = (fzo + (((nzo - 1i32) as f32) * dzo));
                /* Check if the point is out of output range */
                if ((((xo < fxo) || (xo > exo)) || (zo < fzo)) || (zo > ezo)) {
                    it += 1;
                    continue;
                }
                xc = xo;
                zc = zo;
                px = (ss[(it) as usize] / vd1);
                pz = (cc[(it) as usize] / vd1);
                /* set parameters of the first auxilliary ray */
                x1 = ray.rs[(it) as usize].x;
                z1 = ray.rs[(it) as usize].z;
                c1 = ray.rs[(it) as usize].c;
                s1 = ray.rs[(it) as usize].s;
                v1 = ray.rs[(it) as usize].v;
                /* Compute the second derivatives of traveltimes */
                /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
                /* compute Xij and Yij */
                X[(0i32) as usize][(0i32) as usize] = ((((x1 - xc) as f64) / 0.0017f64) as f32);
                X[(1i32) as usize][(0i32) as usize] = ((((z1 - zc) as f64) / 0.0017f64) as f32);
                Y[(0i32) as usize][(0i32) as usize] =
                    (((((s1 / v1) - px) as f64) / 0.0017f64) as f32);
                Y[(1i32) as usize][(0i32) as usize] =
                    (((((c1 / v1) - pz) as f64) / 0.0017f64) as f32);
                X[(0i32) as usize][(1i32) as usize] =
                    ((xx[(it + 1i32) as usize] - xx[(it) as usize]) / dt);
                X[(1i32) as usize][(1i32) as usize] =
                    ((zz[(it + 1i32) as usize] - zz[(it) as usize]) / dt);
                Y[(0i32) as usize][(1i32) as usize] = (((ss[(it + 1i32) as usize]
                    / vv[(it + 1i32) as usize])
                    - (ss[(it) as usize] / vv[(it) as usize]))
                    / dt);
                Y[(1i32) as usize][(1i32) as usize] = (((cc[(it + 1i32) as usize]
                    / vv[(it + 1i32) as usize])
                    - (cc[(it) as usize] / vv[(it) as usize]))
                    / dt);
                /* compute Nij */
                Det = ((X[(0i32) as usize][(0i32) as usize] * X[(1i32) as usize][(1i32) as usize])
                    - (X[(0i32) as usize][(1i32) as usize] * X[(1i32) as usize][(0i32) as usize]));
                N[(0i32) as usize][(0i32) as usize] = (((Y[(0i32) as usize][(0i32) as usize]
                    * X[(1i32) as usize][(1i32) as usize])
                    - (Y[(0i32) as usize][(1i32) as usize] * X[(1i32) as usize][(0i32) as usize]))
                    / Det);
                N[(0i32) as usize][(1i32) as usize] = ((((-Y[(0i32) as usize][(0i32) as usize])
                    * X[(0i32) as usize][(1i32) as usize])
                    + (Y[(0i32) as usize][(1i32) as usize] * X[(0i32) as usize][(0i32) as usize]))
                    / Det);
                N[(1i32) as usize][(0i32) as usize] = (((Y[(1i32) as usize][(0i32) as usize]
                    * X[(1i32) as usize][(1i32) as usize])
                    - (Y[(1i32) as usize][(1i32) as usize] * X[(1i32) as usize][(0i32) as usize]))
                    / Det);
                N[(1i32) as usize][(1i32) as usize] = ((((-Y[(1i32) as usize][(0i32) as usize])
                    * X[(0i32) as usize][(1i32) as usize])
                    + (Y[(1i32) as usize][(1i32) as usize] * X[(0i32) as usize][(0i32) as usize]))
                    / Det);
                /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
                /* determine the radius of extrapolation */
                tc = ((it as f32) * dt);
                terr = (tc * fac);
                norm2 = (((((N[(0i32) as usize][(0i32) as usize]
                    * N[(0i32) as usize][(0i32) as usize])
                    + (N[(1i32) as usize][(1i32) as usize] * N[(1i32) as usize][(1i32) as usize]))
                    + (((2i32 as f32) * N[(0i32) as usize][(1i32) as usize])
                        * N[(0i32) as usize][(1i32) as usize])) as f64)
                    .sqrt() as f32);
                r2 = (terr / norm2);
                r1 = ((r2 as f64).sqrt() as f32);
                /* Make sure that the coordinates are regular */
                if (((r1 * curv) as f64) > 0.1f64) {
                    r1 = ((0.1f64 / (curv as f64)) as f32);
                }
                /* Ensure that the radius is not large */
                if ((r1 as f64) > (0.1f64 * (sd as f64))) {
                    r1 = ((0.1f64 * (sd as f64)) as f32);
                }
                r2 = (r1 * r1);
                nxf = ((((((xc - r1) - fxo) / dxo) as f64) + 0.9f64) as i32);
                if (nxf < 0i32) {
                    nxf = 0i32;
                }
                nxe = ((((((xc + r1) - fxo) / dxo) as f64) + 0.1f64) as i32);
                if (nxe >= nxo) {
                    nxe = (nxo - 1i32);
                }
                ixo = nxf;
                while (ixo <= nxe) {
                    xoc = ((fxo - xc) + ((ixo as f32) * dxo));
                    r2x = (r2 - (xoc * xoc));
                    if (r2x < (0i32 as f32)) {
                        ixo += 1;
                        continue;
                    }
                    r1x = ((r2x as f64).sqrt() as f32);
                    t1x = (tc + (xoc * px));
                    t2x = ((xoc * xoc) * N[(0i32) as usize][(0i32) as usize]);
                    t2xz = (((2i32 as f32) * xoc) * N[(0i32) as usize][(1i32) as usize]);
                    nzf = ((((((zc - r1x) - fzo) / dzo) as f64) + 0.9f64) as i32);
                    if (nzf < 0i32) {
                        nzf = 0i32;
                    }
                    nze = ((((((zc + r1x) - fzo) / dzo) as f64) + 0.1f64) as i32);
                    if (nze >= nzo) {
                        nze = (nzo - 1i32);
                    }
                    izo = nzf;
                    while (izo <= nze) {
                        if (sd < s[((nzo * ixo) + izo) as usize]) {
                            zoc = ((fzo - zc) + ((izo as f32) * dzo));
                            t1 = (t1x + (zoc * pz));
                            t2 = (t2x
                                + (zoc * ((zoc * N[(1i32) as usize][(1i32) as usize]) + t2xz)));
                            s[((nzo * ixo) + izo) as usize] = sd;
                            t[((nzo * ixo) + izo) as usize] =
                                (((t1 as f64) + (0.5f64 * (t2 as f64))) as f32);
                        }

                        izo += 1;
                    }
                    ixo += 1;
                }
                it += 1;
            }
            /* free ray */
            ia += 1;
        }
        fwrite(
            &t[..(nxo * nzo).max(0) as usize],
            tfp.as_ref().unwrap_or_else(|| ray::crash()),
        );
        isx += 1;
    }

    exit(CWP_Exit());
}

/* Read the native-float model grids, retaining the program's error messages. */
fn read_grid(file: &str, field: &str, grid: &mut [f32]) {
    let fp = fopen(file, "r").unwrap_or_else(|| match field {
        "VP0file" => err!("error opening VP0file=%s!\n", file),
        "a1111file" => err!("error opening a1111file=%s!\n", file),
        "a1133file" => err!("error opening a1133file=%s!\n", file),
        "VS0file" => err!("error opening VS0file=%s!\n", file),
        "a1113file" => err!("error opening a1113file=%s!\n", file),
        "a3313file" => err!("error opening a3313file=%s!\n", file),
        "deltafile" => err!("error opening deltafile=%s!\n", file),
        "epsilonfile" => err!("error opening epsilonfile=%s!\n", file),
        _ => unreachable!(),
    });
    if fread(grid, &fp) != grid.len() {
        match field {
            "VP0file" => err!("error reading VP0file=%s!\n", file),
            "a1111file" => err!("error reading a1111file=%s!\n", file),
            "a1133file" => err!("error reading a1133file=%s!\n", file),
            "VS0file" => err!("error reading VS0file=%s!\n", file),
            "a1113file" => err!("error reading a1113file=%s!\n", file),
            "a3313file" => err!("error reading a3313file=%s!\n", file),
            "deltafile" => err!("error reading deltafile=%s!\n", file),
            "epsilonfile" => err!("error reading epsilonfile=%s!\n", file),
            _ => unreachable!(),
        }
    }
}
