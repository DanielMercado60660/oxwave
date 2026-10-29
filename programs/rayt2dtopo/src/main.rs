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

//! RAYT2DTOPO - traveltime tables calculated by 2D paraxial ray tracing.
//! Port of par/main/ray_theory/rayt2dtopo.c, Revision 1.2.
//! Copyright (c) Colorado School of Mines, 1995. All rights reserved.
//! Author: Zhenyue Liu, Colorado School of Mines, 10/11/94.
//! Trino Salinas, 01/01/96: nonflat reference surfaces.
//! Subroutines from Dave Hale's modeling library were adapted in this code
//! to define topography using cubic splines (Dave Hale, CSM, 09/17/91).

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;
mod numerics;
mod surface;
use numerics::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " RAYT2DTOPO - traveltime Tables calculated by 2D paraxial RAY tracing\t",
    "\t\t\t\t\t\t\t\t\t",
    "     rayt2d vfile= tfile= [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " vfile=stdin\t\tfile containning velocity v[nx][nz]\t\t",
    " tfile=stdout\t\tfile containning traveltime tables\t\t",
    "\t\t\tt[nxs][nxo][nzo]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters\t\t\t\t\t\t\t",
    " dt=0.008\t\ttime sample interval in ray tracing\t\t",
    " nt=401\t\tnumber of time samples in ray tracing\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " fz=0\t\t\tfirst depth sample in velocity\t\t\t",
    " nz=101\t\tnumber of depth samples in velocity\t\t",
    " dz=100\t\tdepth interval in velocity\t\t\t",
    " fx=0\t\t\tfirst lateral sample in velocity\t\t",
    " nx=101\t\tnumber of lateral samples in velocity\t\t",
    " dx=100\t\tlateral interval in velocity\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " fzo=fz\t\tfirst depth sample in traveltime table\t\t",
    " nzo=nz\t\tnumber of depth samples in traveltime table\t",
    " dzo=dz\t\tdepth interval in traveltime table\t\t",
    " fxo=fx\t\tfirst lateral sample in traveltime table\t",
    " nxo=nx\t\tnumber of lateral samples in traveltime table\t",
    " dxo=dx\t\tlateral interval in traveltime table\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " surf=\"0,0;99999,0\"  Recording surface \"x1,z1;x2,z2;x3,z3;...\"\t",
    " fxs=fx\t\tx-coordinate of first source\t\t\t",
    " nxs=1\t\t\tnumber of sources\t\t\t\t",
    " dxs=2*dxo\t\tx-coordinate increment of sources\t\t",
    " aperx=0.5*nx*dx  \tray tracing aperature in x-direction\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " fa=-60\t\tfirst take-off angle of rays (degrees)\t\t",
    " na=61\t\t\tnumber of rays  \t\t\t\t",
    " da=2\t\t\tincrement of take-off angle  \t\t\t",
    " amin=0\t\tminimum emergence angle \t\t\t",
    " amax=90\t\tmaximum emergence angle \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " fac=0.01\t\tfactor to determine radius for extrapolation\t",
    " ek=1\t\t\tflag of implementing eikonal in shadow zones \t",
    " ms=10\t\t\tprint verbal information at every ms sources\t",
    " restart=n\t\tjob is restarted (=y yes; =n no)\t\t",
    " npv=0\t\t\tflag of computing quantities for velocity analysis",
    " if npv>0 specify the following three files\t\t\t\t",
    " pvfile=pvfile\t\tinput file of velocity variation pv[nxo][nzo]\t",
    " tvfile=tvfile\t\toutput file of traveltime variation tables  \t",
    "\t\t\ttv[nxs][nxo][nzo]\t\t\t\t",
    " csfile=csfile\t\toutput file of cosine tables cs[nxs][nxo][nzo]\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " 1. Each traveltime table is calculated by paraxial ray tracing; then \t",
    "    traveltimes in shadow zones are compensated by solving eikonal\t",
    "    equation.\t\t\t\t\t\t\t\t",
    " 2. Input velocity is uniformly sampled and smooth one preferred.\t",
    " 3. Traveltime table and source ranges must be within velocity model.\t",
    " 4. Ray tracing aperature can be chosen as sum of migration aperature\t",
    "    plus half of maximum offset.\t\t\t\t\t",
    " 5. Memory requirement for this program is about\t\t\t",
    "      [nx*nz+4*mx*nz+3*nxo*nzo+npa*(nx*nz+mx*nz+3*nxo*nzo)]*4 bytes\t",
    "    where mx = min(nx,2*(1+aperx/dx)).\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

#[derive(Clone, Copy, Default)]
struct Geo2d {
    nx: i32,
    fx: f32,
    dx: f32,
    odx: f32,
    nz: i32,
    fz: f32,
    dz: f32,
    odz: f32,
}
#[derive(Clone, Copy, Default)]
struct Raytr {
    xs: f32,
    zs: f32,
    nt: i32,
    dt: f32,
    na: i32,
    fa: f32,
    da: f32,
    amin: f32,
    amax: f32,
    fac: f32,
}
#[derive(Clone, Copy, Default)]
#[allow(dead_code)]
struct Ray {
    t: f32,
    px: f32,
    pz: f32,
    x: f32,
    z: f32,
    q2: f32,
    p2: f32,
    v: f32,
    dvdx: f32,
    dvdz: f32,
    dtv: f32,
}

// The C program uses unchecked fopen and faults when a null stream is used.
fn stream(fp: Option<CFile>) -> CFile {
    fp.unwrap_or_else(|| fault(11))
}
fn fault(sig: i32) -> ! {
    unsafe extern "C" {
        fn signal(sig: i32, handler: usize) -> usize;
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        signal(sig, 0);
        raise(sig);
    }
    exit(128 + sig)
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut vfile = String::from("stdin");
    let vfp = if getparstring("vfile", &mut vfile) == 0 {
        Some(stdin())
    } else {
        fopen(&vfile, "r")
    };
    let (mut nx, mut nz) = (101i32, 101i32);
    getparint("nx", &mut nx);
    getparint("nz", &mut nz);
    if nx < 3 || nz < 3 {
        err!("nx and nz must be not less than 3!\n");
    }
    let (mut fx, mut fz, mut dx, mut dz) = (0.0f32, 0.0f32, 100.0f32, 100.0f32);
    getparfloat("fx", &mut fx);
    getparfloat("fz", &mut fz);
    getparfloat("dx", &mut dx);
    getparfloat("dz", &mut dz);
    let ex = fx + (nx - 1) as f32 * dx;
    let ez = fz + (nz - 1) as f32 * dz;
    let mut geo2dv = Geo2d {
        nx,
        fx,
        dx,
        odx: (1.0 / dx as f64) as f32,
        nz,
        fz,
        dz,
        odz: (1.0 / dz as f64) as f32,
    };
    let (mut nxo, mut nzo) = (nx, nz);
    getparint("nxo", &mut nxo);
    getparint("nzo", &mut nzo);
    if nxo < 3 || nzo < 3 {
        err!("nxo and nzo must be not less than 3!\n");
    }
    let (mut fxo, mut fzo, mut dxo, mut dzo) = (fx, fz, dx, dz);
    getparfloat("fxo", &mut fxo);
    getparfloat("fzo", &mut fzo);
    getparfloat("dxo", &mut dxo);
    getparfloat("dzo", &mut dzo);
    let exo = fxo + (nxo - 1) as f32 * dxo;
    let ezo = fzo + (nzo - 1) as f32 * dzo;
    let geo2dt = Geo2d {
        nx: nxo,
        fx: fxo,
        dx: dxo,
        odx: (1.0 / dxo as f64) as f32,
        nz: nzo,
        fz: fzo,
        dz: dzo,
        odz: (1.0 / dzo as f64) as f32,
    };
    let mut nxs = 1i32;
    getparint("nxs", &mut nxs);
    let mut fxs = fx;
    getparfloat("fxs", &mut fxs);
    let mut dxs = 2.0f32 * dxo;
    getparfloat("dxs", &mut dxs);
    let exs = fxs + (nxs - 1) as f32 * dxs;
    let mut aperx = (0.5 * (ex - fx) as f64) as f32;
    getparfloat("aperx", &mut aperx);
    if nxs > 1 {
        aperx += dxs;
    }
    let mut mx = (2.0 * aperx as f64 / dx as f64 + 2.99) as i32;
    if mx > nx {
        mx = nx;
    }
    let mut nt = 401i32;
    getparint("nt", &mut nt);
    let mut dt = 0.008f32;
    getparfloat("dt", &mut dt);
    let tmax = (nt - 1) as f32 * dt;
    let (mut fa, mut da, mut na) = (-60.0f32, 2.0f32, 61i32);
    getparfloat("fa", &mut fa);
    getparfloat("da", &mut da);
    getparint("na", &mut na);
    let (mut amin, mut amax) = (0.0f32, 90.0f32);
    getparfloat("amin", &mut amin);
    getparfloat("amax", &mut amax);
    if amax > 180.0 || amin < 0.0 {
        err!("amin and amx must be within 0 to 180 degrees!\n");
    }
    fa = (fa as f64 * (PI / 180.0)) as f32;
    da = (da as f64 * (PI / 180.0)) as f32;
    let ea = fa + (na - 1) as f32 * da;
    amin = (amin as f64 * (PI / 180.0)) as f32;
    amax = (amax as f64 * (PI / 180.0)) as f32;
    let mut fac = 0.01f32;
    getparfloat("fac", &mut fac);
    let mut raytr = Raytr {
        nt,
        dt,
        da,
        fac,
        amin,
        amax,
        ..Default::default()
    };
    let mut ek = 1i32;
    getparint("ek", &mut ek);
    let mut restart = String::from("n");
    getparstring("restart", &mut restart);
    let mut ms = 10i32;
    getparint("ms", &mut ms);
    let mut npv = 0i32;
    getparint("npv", &mut npv);
    let nsize = nzo * nxo;
    let mut ixs0 = 0i32;
    let mut tfile = String::from("stdout");
    let tfp = if getparstring("tfile", &mut tfile) == 0 {
        stdout()
    } else {
        let fp = if let Some(old) = fopen(&tfile, "r") {
            fclose(&old);
            stream(fopen(&tfile, "r+"))
        } else {
            stream(fopen(&tfile, "w"))
        };
        let fp = if restart.starts_with('y') {
            fseeko(&fp, 0, SEEK_END);
            ixs0 = (ftello(&fp) / 4 / (nsize as i64)) as i32;
            fp
        } else {
            fclose(&fp);
            stream(fopen(&tfile, "w"))
        };
        fseeko(&fp, ixs0 as i64 * nsize as i64 * 4, SEEK_SET);
        fp
    };
    let mut jpfile = String::new();
    let jpfp = if getparstring("jpfile", &mut jpfile) == 0 {
        stderr()
    } else {
        stream(fopen(&jpfile, "w"))
    };
    let (mut pvfile, mut tvfile, mut csfile) = (
        String::from("pvfile"),
        String::from("tvfile"),
        String::from("csfile"),
    );
    let (mut pvfp, mut tvfp, mut csfp) = (None, None, None);
    if npv != 0 {
        getparstring("pvfile", &mut pvfile);
        pvfp = fopen(&pvfile, "r");
        getparstring("tvfile", &mut tvfile);
        tvfp = fopen(&tvfile, "w");
        getparstring("csfile", &mut csfile);
        csfp = fopen(&csfile, "w");
    }
    fprintf!(&jpfp, "\n");
    fprintf!(&jpfp, " RAYT2D parameters\n");
    fprintf!(&jpfp, " ================\n");
    fprintf!(&jpfp, " vfile=%s \n", vfile);
    fprintf!(&jpfp, " tfile=%s \n", tfile);
    fprintf!(&jpfp, " one-way time: dt=%g nt=%d \n", dt, nt);
    fprintf!(&jpfp, " nz=%d fz=%g dz=%g\n", nz, fz, dz);
    fprintf!(&jpfp, " nx=%d fx=%g dx=%g\n", nx, fx, dx);
    fprintf!(&jpfp, " nzo=%d fzo=%g dzo=%g\n", nzo, fzo, dzo);
    fprintf!(&jpfp, " nxo=%d fxo=%g dxo=%g\n", nxo, fxo, dxo);
    fprintf!(&jpfp, " nxs=%d fxs=%g dxs=%g\n", nxs, fxs, dxs);
    fprintf!(&jpfp, " mx=%d aperx=%g \n", mx, aperx);
    fprintf!(
        &jpfp,
        " na=%d fa=%g da=%g\n",
        na,
        (fa * 180.0f32) as f64 / PI,
        (da * 180.0f32) as f64 / PI
    );
    fprintf!(
        &jpfp,
        " amin=%g amax=%g \n",
        (amin * 180.0f32) as f64 / PI,
        (amax * 180.0f32) as f64 / PI
    );
    fprintf!(
        &jpfp,
        " ek=%d ms=%d npv=%d restart=%s\n",
        ek,
        ms,
        npv,
        restart
    );
    if npv != 0 {
        fprintf!(
            &jpfp,
            " pvfile=%s tvfile=%s csfile=%s\n",
            pvfile,
            tvfile,
            csfile
        );
    }
    fflush(&jpfp);
    if fx > fxs || ex < exs || fz > 0.0 {
        err!("source lies outside of specified x grid  \n");
    }
    if fx > fxo || ex < exo || fz > fzo || ez < ezo {
        err!("output lies outside of specified x grid  \n");
    }
    let mut v = vec![0.0f32; (nx * nz) as usize];
    let mut ov2 = vec![0.0f32; nsize as usize];
    let mut t = vec![0.0f32; nsize as usize];
    let mut pv = if npv != 0 {
        vec![0.0f32; (nx * nz) as usize]
    } else {
        Vec::new()
    };
    let mut tv = if npv != 0 {
        vec![0.0f32; nsize as usize]
    } else {
        Vec::new()
    };
    let mut vo = if npv != 0 {
        vec![0.0f32; nsize as usize]
    } else {
        Vec::new()
    };
    let mut cs = if npv != 0 {
        vec![0.0f32; nsize as usize]
    } else {
        Vec::new()
    };
    let vfp = stream(vfp);
    if fread(&mut v, &vfp) != v.len() {
        // The C format mistakenly passes FILE* as %s (undefined); an empty
        // string is printed here.
        err!("cannot read %d velocities from file %s", nx * nz, "");
    }
    if npv != 0 {
        let fp = stream(pvfp);
        if fread(&mut pv, &fp) != pv.len() {
            err!("cannot read %d values from file %s", nx * nz, "");
        }
    }
    fprintf!(&jpfp, " finish velocity input\n");
    unsafe {
        voint(
            geo2dv,
            v.as_mut_ptr(),
            geo2dt,
            ov2.as_mut_ptr(),
            npv,
            vo.as_mut_ptr(),
        );
    }
    fprintf!(&jpfp, " begin traveltime calculation\n");
    fflush(&jpfp);
    let srf = surface::decode_surfaces();
    let srf = surface::makesurf((dxs as f64 * 0.025) as f32, &srf);
    // The unchecked szi allocation fails first; the following ealloc1float
    // for nangl reports this error for a negative source count.
    if nxs < 0 {
        err!(
            "%s: malloc failed (%s)",
            "ealloc.c",
            "Cannot allocate memory"
        );
    }
    let (szi, nangl) = surface::zcoor_topog(fxs, dxs, nxs, &srf);
    let mut ixs = ixs0;
    let mut xs = fxs + ixs0 as f32 * dxs;
    while ixs < nxs {
        let mut temp = fxo;
        if xs < temp {
            temp = xs;
        }
        if xs - aperx > temp {
            temp = xs - aperx;
        }
        let nx0 = ((temp - fx) / dx) as i32;
        let fxt = fx + nx0 as f32 * dx;
        temp = exo;
        if xs > temp {
            temp = xs;
        }
        if xs + aperx < temp {
            temp = xs + aperx;
        }
        let nxt = 1 + (((temp - fx) / dx) as f64 + 0.99) as i32 - nx0;
        let ext = fxt + (nxt - 1) as f32 * dx;
        let mut vt = vec![0.0f32; (nxt * nz) as usize];
        let mut pvt = if npv != 0 {
            vec![0.0f32; (nxt * nz) as usize]
        } else {
            Vec::new()
        };
        unsafe {
            trans(nx, nz, nxt, nx0, v.as_mut_ptr(), vt.as_mut_ptr());
        }
        if npv != 0 {
            unsafe {
                trans(nx, nz, nxt, nx0, pv.as_mut_ptr(), pvt.as_mut_ptr());
            }
        }
        let mut fat = fa + nangl[ixs as usize];
        let mut eat = ea + nangl[ixs as usize];
        let mut nat = na;
        if (fat as f64) < -(PI / 2.0) {
            fat = (-PI / 2.0) as f32;
        }
        if eat as f64 > PI / 2.0 {
            eat = (PI / 2.0) as f32;
        }
        if xs == fxt && fat < 0.0 {
            fat = 0.0;
            nat = ((eat / da) as f64 + 1.5) as i32;
        } else if xs == ext && eat > 0.0 {
            nat = (1.5 - (fa / da) as f64) as i32;
        }
        raytr.xs = xs;
        raytr.zs = szi[ixs as usize];
        raytr.na = nat;
        raytr.fa = fat;
        geo2dv.nx = nxt;
        geo2dv.fx = fxt;
        unsafe {
            raytime2d(
                raytr,
                geo2dv,
                vt.as_mut_ptr(),
                geo2dt,
                t.as_mut_ptr(),
                npv,
                vo.as_mut_ptr(),
                pvt.as_mut_ptr(),
                tv.as_mut_ptr(),
                cs.as_mut_ptr(),
            );
        }
        if ek != 0 {
            unsafe {
                eiknl(geo2dt, t.as_mut_ptr(), ov2.as_mut_ptr(), tmax);
            }
        }
        fwrite(&t, &tfp);
        if npv != 0 {
            fwrite(&tv, &stream(tvfp.clone()));
            fwrite(&cs, &stream(csfp.clone()));
        }
        // arm64 integer division/remainder by zero yields zero.
        if ixs.checked_rem(ms).unwrap_or(0) == 0 {
            fprintf!(&jpfp, " traveltime computed at source xs=%g\n", xs);
            fflush(&jpfp);
        }
        ixs += 1;
        xs += dxs;
    }
    fprintf!(&jpfp, " finish program rayt2d\n\n");
    fclose(&tfp);
    fclose(&vfp);
    fclose(&jpfp);
    exit(CWP_Exit());
}
