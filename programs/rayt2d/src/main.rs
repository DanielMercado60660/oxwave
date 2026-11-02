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

//! RAYT2D - traveltime tables calculated by 2D paraxial ray tracing.
//! Port of par/main/ray_theory/rayt2d.c (Revision 1.19).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: Zhenyue Liu, 10/11/94, Colorado School of Mines.
//! Trino Salinas, 01/01/96: nonflat reference surfaces.
//! Subroutines from Dave Hale's modeling library were adapted in this code
//! to define topography using cubic splines (makeref: Dave Hale, 09/17/91).

#![allow(non_snake_case)]
use oxwave_core::prelude::*;
use oxwave_cubic_interp::*;
use oxwave_f_cubicspline::*;
use oxwave_inverse_interp::*;
use std::f64::consts::PI;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " RAYT2D - traveltime Tables calculated by 2D paraxial RAY tracing\t",
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
    "      [nx*nz+4*mx*nz+3*nxo*nzo+na*(nx*nz+mx*nz+3*nxo*nzo)]*4 bytes\t",
    "    where mx = min(nx,2*(1+aperx/dx)).\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: spatial units of v(z,x) must be the same as those of dx. \t",
    " v(z,x) is represented numerically in C-style binary floats v[xn][zn],\t",
    " where the depth direction is the fast direction in the data. Such     ",
    " models can be created with unif2 or makevel.    \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

const INITIAL_T: f32 = 999999.0;
#[derive(Clone, Copy, Default)]
struct Ray {
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
#[derive(Clone, Copy)]
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
struct SurfaceSegment {
    x: f32,
    z: f32,
    s: f32,
    c: f32,
}
struct Surface {
    ss: Vec<SurfaceSegment>,
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut vfile = "stdin".to_string();
    let vfp = if getparstring("vfile", &mut vfile) == 0 {
        Some(stdin())
    } else {
        fopen(&vfile, "r")
    };
    let mut nx = 101;
    getparint("nx", &mut nx);
    let mut nz = 101;
    getparint("nz", &mut nz);
    if nx < 3 || nz < 3 {
        err!("nx and nz must be not less than 3!\n");
    }
    let mut fx = 0.0f32;
    getparfloat("fx", &mut fx);
    let mut fz = 0.0f32;
    getparfloat("fz", &mut fz);
    let mut dx = 100.0f32;
    getparfloat("dx", &mut dx);
    let mut dz = 100.0f32;
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
    let mut nxo = nx;
    getparint("nxo", &mut nxo);
    let mut nzo = nz;
    getparint("nzo", &mut nzo);
    if nxo < 3 || nzo < 3 {
        err!("nxo and nzo must be not less than 3!\n");
    }
    let mut fxo = fx;
    getparfloat("fxo", &mut fxo);
    let mut fzo = fz;
    getparfloat("fzo", &mut fzo);
    let mut dxo = dx;
    getparfloat("dxo", &mut dxo);
    let mut dzo = dz;
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
    let mut nxs = 1;
    getparint("nxs", &mut nxs);
    let mut fxs = fx;
    getparfloat("fxs", &mut fxs);
    let mut dxs = 2.0 * dxo;
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
    let mut nt = 401;
    getparint("nt", &mut nt);
    let mut dt = 0.008f32;
    getparfloat("dt", &mut dt);
    let tmax = (nt - 1) as f32 * dt;
    let mut fa = -60.0f32;
    getparfloat("fa", &mut fa);
    let mut da = 2.0f32;
    getparfloat("da", &mut da);
    let mut na = 61;
    getparint("na", &mut na);
    let mut amin = 0.0f32;
    getparfloat("amin", &mut amin);
    let mut amax = 90.0f32;
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
        ..Raytr::default()
    };
    let mut ek = 1;
    getparint("ek", &mut ek);
    let mut restart = "n".to_string();
    getparstring("restart", &mut restart);
    let mut ms = 10;
    getparint("ms", &mut ms);
    let mut npv = 0;
    getparint("npv", &mut npv);
    let mut ixs0 = 0;
    let mut isize = 0;
    let nsize = nzo * nxo;
    let mut tfile = "stdout".to_string();
    let tfp = if getparstring("tfile", &mut tfile) == 0 {
        stdout()
    } else {
        let fp = if let Some(fp) = fopen(&tfile, "r") {
            fclose(&fp);
            fopen(&tfile, "r+")
        } else {
            fopen(&tfile, "w")
        };
        let fp = if restart.starts_with('y') {
            let fp = fp.unwrap_or_else(|| null_file());
            efseeko(&fp, 0, SEEK_END);
            isize = (eftello(&fp) / 4 / (nsize as i64)) as i32;
            ixs0 = isize;
            fp
        } else {
            // rayt2d.c:291 can close NULL after an unchecked fopen failure.
            // Match the minimal C fix: skip that close, then use efopen's error.
            if let Some(fp) = fp {
                efclose(&fp);
            }
            efopen(&tfile, "w")
        };
        efseeko(&fp, (isize * nsize) as i64 * 4, SEEK_SET);
        fp
    };
    let mut jpfile = String::new();
    let jpfp = if getparstring("jpfile", &mut jpfile) == 0 {
        Some(stderr())
    } else {
        fopen(&jpfile, "w")
    };
    let mut pvfile = String::new();
    let mut tvfile = String::new();
    let mut csfile = String::new();
    let (pvfp, tvfp, csfp) = if npv != 0 {
        if getparstring("pvfile", &mut pvfile) == 0 {
            pvfile = "pvfile".to_string();
        }
        let pvfp = fopen(&pvfile, "r");
        if getparstring("tvfile", &mut tvfile) == 0 {
            tvfile = "tvfile".to_string();
        }
        let tvfp = fopen(&tvfile, "w");
        if getparstring("csfile", &mut csfile) == 0 {
            csfile = "csfile".to_string();
        }
        let csfp = fopen(&csfile, "w");
        (pvfp, tvfp, csfp)
    } else {
        (None, None, None)
    };
    checkpars();
    let jpfp = jpfp.as_ref().unwrap_or_else(|| null_file());
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
        (fa * 180.0) as f64 / PI,
        (da * 180.0) as f64 / PI
    );
    fprintf!(
        &jpfp,
        " amin=%g amax=%g \n",
        (amin * 180.0) as f64 / PI,
        (amax * 180.0) as f64 / PI
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
        warn!("This condition must NOT be satisfied: fx>fxs || ex<exs || fz>0");
        warn!("fx=%g fxs=%g ex=%g exs=%g fz=%g", fx, fxo, ex, exo, fz);
        err!("source lies outside of specified x grid  \n");
    }
    if fx > fxo || ex < exo || fz > fzo || ez < ezo {
        warn!("This condition must NOT be satisfied: fx>fxo || ex<exo || fz>fzo || ez<ezo");
        warn!(
            "fx=%g fxo=%g ex=%g exo=%g fz=%g fzo=%g,ez=%g ezo=%g",
            fx, fxo, ex, exo, fz, fzo, ez, ezo
        );
        err!("output lies outside of specified x grid  \n");
    }
    let mut v = vec![0.0f32; (nx * nz) as usize];
    let mut ov2 = vec![0.0f32; nsize as usize];
    let mut t = vec![0.0f32; nsize as usize];
    let mut pv = if npv != 0 {
        vec![0.0f32; (nx * nz) as usize]
    } else {
        vec![]
    };
    let mut tv = if npv != 0 {
        vec![0.0f32; nsize as usize]
    } else {
        vec![]
    };
    let mut vo = tv.clone();
    let mut cs = tv.clone();
    let vfp = vfp.as_ref().unwrap_or_else(|| null_file());
    if fread(&mut v, &vfp) != (nx * nz) as usize {
        // C formats FILE* as a string; its first byte is NUL on exhausted input.
        err!("cannot read %d velocities from file %s", nx * nz, "");
    }
    if npv != 0 {
        let fp = pvfp.as_ref().unwrap_or_else(|| null_file());
        if fread(&mut pv, fp) != (nx * nz) as usize {
            err!("cannot read %d values from file %s", nx * nz, "");
        }
    }
    fprintf!(&jpfp, " finish velocity input\n");
    voint(geo2dv, &v, geo2dt, &mut ov2, npv, &mut vo);
    fprintf!(&jpfp, " begin traveltime calculation\n");
    fflush(&jpfp);
    let surfaces = decodeSurfaces();
    let srf = makesurf((dxs as f64 * 0.025) as f32, &surfaces);
    let (szi, nangl) = zcoorTopog(fxs, dxs, nxs, &srf[0]);
    let mut xs = fxs + ixs0 as f32 * dxs;
    for ixs in ixs0..nxs {
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
            vec![]
        };
        trans(nx, nz, nxt, nx0, &v, &mut vt);
        if npv != 0 {
            trans(nx, nz, nxt, nx0, &pv, &mut pvt);
        }
        let mut fat = fa + nangl[ixs as usize];
        let mut eat = ea + nangl[ixs as usize];
        let mut nat = na;
        if (fat as f64) < -(PI / 2.0) {
            fat = (-PI / 2.0) as f32;
        }
        if (eat as f64) > PI / 2.0 {
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
        raytime2d(
            raytr, geo2dv, &vt, geo2dt, &mut t, npv, &vo, &pvt, &mut tv, &mut cs,
        );
        if ek != 0 {
            eiknl(geo2dt, &mut t, &ov2, tmax);
        }
        fwrite(&t, &tfp);
        if npv != 0 {
            fwrite(&tv, tvfp.as_ref().unwrap_or_else(|| null_file()));
            fwrite(&cs, csfp.as_ref().unwrap_or_else(|| null_file()));
        }
        // ARM64 C remainder with a zero divisor yields its dividend.
        if ixs.checked_rem(ms).unwrap_or(ixs) == 0 {
            fprintf!(&jpfp, " traveltime computed at source xs=%g\n", xs);
            fflush(&jpfp);
        }
        xs += dxs;
    }
    fprintf!(&jpfp, " finish program rayt2d\n\n");
    fclose(&tfp);
    fclose(&vfp);
    fclose(&jpfp);
    exit(CWP_Exit());
}

fn null_file() -> ! {
    unsafe extern "C" {
        fn raise(sig: i32) -> i32;
        fn signal(sig: i32, handler: usize) -> usize;
    }
    unsafe {
        // Rust installs a stack-overflow SIGSEGV handler; C uses SIG_DFL here.
        signal(11, 0);
        raise(11);
    }
    exit(139)
}

fn decodeSurfaces() -> Vec<(Vec<f32>, Vec<f32>)> {
    let mut nr = countparname("surf");
    if nr == 0 {
        nr = 1;
    }
    let mut result = Vec::new();
    for ir in 0..nr {
        let mut s = String::new();
        if getnparstring(ir + 1, "surf", &mut s) == 0 {
            s = "0,0;99999,0".to_string();
        }
        if let Some(surface) = decodeSurface(&s) {
            result.push(surface);
        } else {
            err!("Surface number %d specified incorrectly!\n", ir + 1);
        }
    }
    result
}

fn decodeSurface(string: &str) -> Option<(Vec<f32>, Vec<f32>)> {
    let mut t = string.as_bytes().to_vec();
    t.push(0);
    let mut cursor = 0;
    let mut first = 0;
    let mut nxz = 0;
    // strtok replaces each token's terminating delimiter, while leaving skipped
    // delimiters intact. Keep that detail: C walks the mutated buffer below.
    loop {
        while t[cursor] == b',' || t[cursor] == b';' {
            cursor += 1;
        }
        if t[cursor] == 0 {
            break;
        }
        if nxz == 0 {
            first = cursor;
        }
        nxz += 1;
        while t[cursor] != 0 && t[cursor] != b',' && t[cursor] != b';' {
            cursor += 1;
        }
        if t[cursor] == 0 {
            break;
        }
        t[cursor] = 0;
        cursor += 1;
    }
    if nxz % 2 != 0 {
        return None;
    }
    nxz /= 2;
    if nxz < 2 {
        return None;
    }
    let mut x = vec![0.0; nxz];
    let mut z = vec![0.0; nxz];
    let mut s = first;
    for ixz in 0..nxz {
        let end = s + t[s..].iter().position(|&c| c == 0).unwrap();
        x[ixz] = atof(&String::from_utf8_lossy(&t[s..end])) as f32;
        s = end + 1;
        let end = s + t[s..].iter().position(|&c| c == 0).unwrap();
        z[ixz] = atof(&String::from_utf8_lossy(&t[s..end])) as f32;
        s = end + 1;
    }
    Some((x, z))
}

fn makesurf(dsmax: f32, surfaces: &[(Vec<f32>, Vec<f32>)]) -> Vec<Surface> {
    let mut rr = Vec::new();
    for (xu, zu) in surfaces {
        let nu = xu.len();
        let mut u = vec![0.0f32; nu];
        for iu in 0..nu {
            u[iu] = iu as f32;
        }
        let mut xud = vec![[0.0f32; 4]; nu];
        csplin(&u, xu, &mut xud);
        let mut zud = vec![[0.0f32; 4]; nu];
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
        if ns < 0 {
            // ealloc1float converts the negative C count to size_t.
            err!("ealloc.c: malloc failed (%s)", "Cannot allocate memory");
        }
        let mut us = vec![0.0f32; ns.max(0) as usize];
        yxtoxy(duu, 0.0, &s, ds, fs, 0.0, (nu - 1) as f32, &mut us);
        let mut ss = vec![SurfaceSegment::default(); ns.max(0) as usize];
        for is in 0..ss.len() {
            let mut rsx = [0.0f32];
            let mut rsz = [0.0f32];
            let mut rsxd = [0.0f32];
            let mut rszd = [0.0f32];
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

fn zcoorTopog(fxs: f32, dxs: f32, nxs: i32, srf: &Surface) -> (Vec<f32>, Vec<f32>) {
    if nxs < 0 {
        err!("ealloc.c: malloc failed (%s)", "Cannot allocate memory");
    }
    let ss = &srf.ss;
    // The C unconditionally accesses samples [-2..2], even for nxs=1 or 2.
    // Pad both sides explicitly; fresh SU allocations have zeroed spare slots.
    let mut sz = vec![0.0f32; nxs.max(3) as usize + 4];
    let mut nangl = sz.clone();
    let mut ixi = 0;
    for ix in 0..nxs {
        let x = fxs + ix as f32 * dxs;
        for ik in ixi..ss.len().saturating_sub(1) {
            if ss[ik].x <= x && ss[ik + 1].x >= x {
                let az = ss[ik].z - ss[ik + 1].z;
                let ax0 = ss[ik + 1].x - x;
                let ax = ss[ik + 1].x - ss[ik].x;
                sz[ix as usize + 2] = ax0 * az / ax + ss[ik + 1].z;
                let az = ss[ik].s - ss[ik + 1].s;
                let temp = ax0 * az / ax + ss[ik + 1].s;
                nangl[ix as usize + 2] = (temp as f64).asin() as f32;
                ixi = ik;
            }
        }
    }
    sz[2] = 2.0 * sz[3] - sz[4];
    nangl[2] = 2.0 * nangl[3] - nangl[4];
    if nxs >= 1 {
        let i = nxs as usize;
        sz[i + 1] = 2.0 * sz[i] - sz[i - 1];
        nangl[i + 1] = 2.0 * nangl[i] - nangl[i - 1];
    }
    for is in 0..nxs as usize {
        warn!("sz[%d]=%g,nangl[%d]=%g", is, sz[is + 2], is, nangl[is + 2]);
    }
    (
        sz[2..2 + nxs.max(0) as usize].to_vec(),
        nangl[2..2 + nxs.max(0) as usize].to_vec(),
    )
}

fn trans(_nx: i32, nz: i32, nxt: i32, nx0: i32, v: &[f32], vt: &mut [f32]) {
    for ix in 0..nxt {
        let i = (ix * nz) as usize;
        let i0 = ((ix + nx0) * nz) as usize;
        for iz in 0..nz as usize {
            vt[i + iz] = v[i0 + iz];
        }
    }
}

fn voint(g: Geo2d, v: &[f32], o: Geo2d, ov2: &mut [f32], npv: i32, vo: &mut [f32]) {
    let mut x = o.fx;
    for ixo in 0..o.nx {
        let ax = (x - g.fx) * g.odx;
        let mut jx = ax as i32;
        let sx = ax - jx as f32;
        if jx < 0 {
            jx = 0;
        }
        if jx > g.nx - 2 {
            jx = g.nx - 2;
        }
        let mut z = o.fz;
        for izo in 0..o.nz {
            let az = (z - g.fz) * g.odz;
            let mut jz = az as i32;
            let sz = az - jz as f32;
            if jz < 0 {
                jz = 0;
            }
            if jz > g.nz - 2 {
                jz = g.nz - 2;
            }
            let io = (ixo * o.nz + izo) as usize;
            let i = (jx * g.nz + jz) as usize;
            let nz = g.nz as usize;
            let temp = ((1.0 - sx as f64)
                * ((1.0 - sz as f64) * v[i] as f64 + (sz * v[i + 1]) as f64)
                + sx as f64 * ((1.0 - sz as f64) * v[i + nz] as f64 + (sz * v[i + nz + 1]) as f64))
                as f32;
            if npv != 0 {
                vo[io] = temp;
            }
            ov2[io] = (1.0 / (temp * temp) as f64) as f32;
            z += o.dz;
        }
        x += o.dx;
    }
}

fn dv2(g: Geo2d, v: &[f32], vxx: &mut [f32], vxz: &mut [f32], vzz: &mut [f32]) {
    let odxx = g.odx * g.odx;
    let odzz = g.odz * g.odz;
    let odxz = (0.25 * g.odx as f64 * g.odz as f64) as f32;
    for ix in 0..g.nx {
        for iz in 0..g.nz {
            let i = (ix * g.nz + iz) as usize;
            vxx[i] = 0.0;
            vxz[i] = 0.0;
            vzz[i] = 0.0;
        }
    }
    let nz = g.nz as usize;
    for ix in 1..g.nx - 1 {
        for iz in 1..g.nz - 1 {
            let i = (ix * g.nz + iz) as usize;
            vxx[i] =
                (odxx as f64 * (v[i + nz] as f64 - 2.0 * v[i] as f64 + v[i - nz] as f64)) as f32;
            vxz[i] = odxz * (v[i + nz + 1] - v[i + nz - 1] - v[i - nz + 1] + v[i - nz - 1]);
            vzz[i] = (odzz as f64 * (v[i + 1] as f64 - 2.0 * v[i] as f64 + v[i - 1] as f64)) as f32;
        }
    }
}

fn eiknl(g: Geo2d, t: &mut [f32], ov2: &[f32], tmax: f32) {
    let nx = g.nx as usize;
    let nz = g.nz as usize;
    let mut tt1 = vec![0.0f32; nx];
    let mut tt2 = tt1.clone();
    let odx2 = (1.0 / (g.dx * g.dx) as f64) as f32;
    for ix in 0..nx {
        tt1[ix] = t[ix * nz];
    }
    for iz in 1..nz {
        for ix in 0..nx {
            tt2[ix] = t[ix * nz + iz];
        }
        for ix in 1..nx {
            let t0 = tt1[ix];
            if tt2[ix] > tmax && t0 <= tmax {
                let mut tl = 0.0f32;
                let mut tr = 0.0f32;
                if ix > 0 {
                    tl = tt1[ix - 1];
                }
                if ix < nx - 1 {
                    tr = tt1[ix + 1];
                }
                let mut tx2 = (tl - t0) * (tl - t0);
                let temp = (tr - t0) * (tr - t0);
                if tx2 > temp {
                    tx2 = temp;
                }
                let temp = (0.25 * (tl - tr) as f64 * (tl - tr) as f64) as f32;
                if tx2 > temp {
                    tx2 = temp;
                }
                tx2 *= odx2;
                let tz2 =
                    (0.5 * (ov2[ix * nz + iz - 1] + ov2[ix * nz + iz]) as f64 - tx2 as f64) as f32;
                if tz2 >= 0.0 {
                    tt2[ix] = (t0 as f64 + g.dz as f64 * (tz2 as f64).sqrt()) as f32;
                }
            }
        }
        for ix in 0..nx {
            t[ix * nz + iz] = tt2[ix];
            tt1[ix] = tt2[ix];
        }
    }
}

// Return velocity and its first and second derivatives, in the C argument order.
fn vel2Interp(
    g: Geo2d,
    v: &[f32],
    vxx: &[f32],
    vxz: &[f32],
    vzz: &[f32],
    x: f32,
    z: f32,
) -> [f32; 6] {
    let ax = (x - g.fx) * g.odx;
    let mut jx = ax as i32;
    if jx < 0 {
        jx = 0;
    }
    if jx > g.nx - 2 {
        jx = g.nx - 2;
    }
    let sx = ax - jx as f32;
    let az = (z - g.fz) * g.odz;
    let mut jz = az as i32;
    if jz < 0 {
        jz = 0;
    }
    if jz > g.nz - 2 {
        jz = g.nz - 2;
    }
    let sz = az - jz as f32;
    let sxx = (0.5 * sx as f64 * (1.0 - sx as f64) * g.dx as f64 * g.dx as f64) as f32;
    let szz = (0.5 * sz as f64 * (1.0 - sz as f64) * g.dz as f64 * g.dz as f64) as f32;
    let a0 = ((1.0 - sx as f64) * (1.0 - sz as f64)) as f32;
    let a1 = ((1.0 - sx as f64) * sz as f64) as f32;
    let a2 = (sx as f64 * (1.0 - sz as f64)) as f32;
    let a3 = sx * sz;
    let k0 = (g.nz * jx + jz) as usize;
    let k1 = k0 + 1;
    let k2 = k0 + g.nz as usize;
    let k3 = k2 + 1;
    let g0 = v[k0];
    let g1 = v[k1];
    let g2 = v[k2];
    let g3 = v[k3];
    let gx0 = vxx[k0];
    let gx1 = vxx[k1];
    let gx2 = vxx[k2];
    let gx3 = vxx[k3];
    let gz0 = vzz[k0];
    let gz1 = vzz[k1];
    let gz2 = vzz[k2];
    let gz3 = vzz[k3];
    let uxx = a0 * gx0 + a1 * gx1 + a2 * gx2 + a3 * gx3;
    let uxz = a0 * vxz[k0] + a1 * vxz[k1] + a2 * vxz[k2] + a3 * vxz[k3];
    let uzz = a0 * gz0 + a1 * gz1 + a2 * gz2 + a3 * gz3;
    let u = a0 * g0 + a1 * g1 + a2 * g2 + a3 * g3 - (sxx * uxx + szz * uzz);
    let ux = (((1.0 - sz as f64) * (g2 - g0 - sxx * (gx2 - gx0) - szz * (gz2 - gz0)) as f64
        + (sz * (g3 - g1 - sxx * (gx3 - gx1) - szz * (gz3 - gz1))) as f64)
        * g.odx as f64
        + (sx as f64 - 0.5) * g.dx as f64 * uxx as f64) as f32;
    let uz = (((1.0 - sx as f64) * (g1 - g0 - sxx * (gx1 - gx0) - szz * (gz1 - gz0)) as f64
        + (sx * (g3 - g2 - sxx * (gx3 - gx2) - szz * (gz3 - gz2))) as f64)
        * g.odz as f64
        + (sz as f64 - 0.5) * g.dz as f64 * uzz as f64) as f32;
    [u, ux, uz, uxx, uxz, uzz]
}

fn vintp(g: Geo2d, v: &[f32], x: f32, z: f32) -> f32 {
    let ax = (x - g.fx) * g.odx;
    let mut jx = ax as i32;
    if jx < 0 {
        jx = 0;
    }
    if jx > g.nx - 2 {
        jx = g.nx - 2;
    }
    let sx = ax - jx as f32;
    let az = (z - g.fz) * g.odz;
    let mut jz = az as i32;
    if jz < 0 {
        jz = 0;
    }
    if jz > g.nz - 2 {
        jz = g.nz - 2;
    }
    let sz = az - jz as f32;
    let k0 = (g.nz * jx + jz) as usize;
    let k1 = k0 + 1;
    let k2 = k0 + g.nz as usize;
    let k3 = k2 + 1;
    ((1.0 - sx as f64) * ((1.0 - sz as f64) * v[k0] as f64 + (sz * v[k1]) as f64)
        + sx as f64 * ((1.0 - sz as f64) * v[k2] as f64 + (sz * v[k3]) as f64)) as f32
}

// The six-component states follow x,z,px,pz,p2,q2 throughout Runge-Kutta.
fn dfrungk(v: f32, vx: f32, vz: f32, vxx: f32, vxz: f32, vzz: f32, y: [f32; 6]) -> [f32; 6] {
    let [x, z, px, pz, p2, q2] = y;
    let ov = (0.0 * z as f64 + 0.0 * x as f64 + 1.0 / v as f64) as f32;
    let vv = v * v;
    let dx = vv * px;
    let dz = vv * pz;
    let dpx = -ov * vx;
    let dpz = -ov * vz;
    let c = v * pz;
    let s = v * px;
    let v11 = ((vxx * c * c) as f64 - 2.0 * vxz as f64 * c as f64 * s as f64 + (vzz * s * s) as f64)
        as f32;
    let dp2 = -ov * v11 * q2;
    let dq2 = vv * p2;
    [dx, dz, dpx, dpz, dp2, dq2]
}

fn sum2(h: f32, a: [f32; 6], b: [f32; 6]) -> [f32; 6] {
    [
        a[0] + h * b[0],
        a[1] + h * b[1],
        a[2] + h * b[2],
        a[3] + h * b[3],
        a[4] + h * b[4],
        a[5] + h * b[5],
    ]
}

fn makeRay(
    g: Geo2d,
    vel: &[f32],
    vxx: &[f32],
    vxz: &[f32],
    vzz: &[f32],
    r: Raytr,
    a0: f32,
    ray: &mut [Ray],
    npv: i32,
    pv: &[f32],
) -> i32 {
    if r.nt <= 0 {
        // C's zero/negative alloc1 result is NULL, and ray[0] is still written.
        null_file();
    }
    let mut y = [r.xs, r.zs, 0.0, 0.0, 1.0, 0.0];
    let [mut v, mut dvdx, mut dvdz, mut uxx, mut uxz, mut uzz] =
        vel2Interp(g, vel, vxx, vxz, vzz, y[0], y[1]);
    let ov = (1.0 / v as f64) as f32;
    y[2] = (ov as f64 * (a0 as f64).sin()) as f32;
    y[3] = (ov as f64 * (a0 as f64).cos()) as f32;
    let mut tzt = v * y[3];
    ray[0] = Ray {
        x: y[0],
        z: y[1],
        px: y[2],
        pz: y[3],
        p2: y[4],
        q2: y[5],
        v,
        dvdx,
        dvdz,
        dtv: 0.0,
    };
    let tzmin = ((r.amax as f64).cos() - 0.01) as f32;
    let tzmax = ((r.amin as f64).cos() + 0.01) as f32;
    let h = r.dt;
    let h2 = (h as f64 * 0.5) as f32;
    let h6 = (h as f64 / 6.0) as f32;
    let lx = g.fx + (g.nx - 1) as f32 * g.dx;
    let lz = g.fz + (g.nz - 1) as f32 * g.dz;
    let mut it = 1;
    while it < r.nt {
        if y[0] < g.fx || y[0] > lx || y[1] < g.fz || y[1] > lz || tzt < tzmin || tzt > tzmax {
            it -= 1;
            break;
        }
        // k1=f(y0); y1=y0+0.5*h*k1
        let mut d = dfrungk(v, dvdx, dvdz, uxx, uxz, uzz, y);
        let mut yt = sum2(h2, y, d);
        [v, dvdx, dvdz, uxx, uxz, uzz] = vel2Interp(g, vel, vxx, vxz, vzz, yt[0], yt[1]);
        // k2=f(y1); y2=y0+0.5*h*k2; k=k1+2*k2
        let mut dt = dfrungk(v, dvdx, dvdz, uxx, uxz, uzz, yt);
        yt = sum2(h2, y, dt);
        d = sum2(2.0, d, dt);
        [v, dvdx, dvdz, uxx, uxz, uzz] = vel2Interp(g, vel, vxx, vxz, vzz, yt[0], yt[1]);
        // k3=f(y2); y3=y0+h*k3; k=k1+2*k2+2*k3
        dt = dfrungk(v, dvdx, dvdz, uxx, uxz, uzz, yt);
        yt = sum2(h, y, dt);
        d = sum2(2.0, d, dt);
        [v, dvdx, dvdz, uxx, uxz, uzz] = vel2Interp(g, vel, vxx, vxz, vzz, yt[0], yt[1]);
        // k4=f(y3); y4=y0+h*(k1+2*k2+2*k3+k4)/6
        dt = dfrungk(v, dvdx, dvdz, uxx, uxz, uzz, yt);
        d = sum2(1.0, d, dt);
        y = sum2(h6, y, d);
        [v, dvdx, dvdz, uxx, uxz, uzz] = vel2Interp(g, vel, vxx, vxz, vzz, y[0], y[1]);
        let dtv = if npv != 0 {
            vintp(g, pv, y[0], y[1])
        } else {
            0.0
        };
        tzt = v * y[3];
        ray[it as usize] = Ray {
            x: y[0],
            z: y[1],
            px: y[2],
            pz: y[3],
            p2: y[4],
            q2: y[5],
            v,
            dvdx,
            dvdz,
            dtv: if npv != 0 { dtv / v } else { 0.0 },
        };
        it += 1;
    }
    it
}

fn ddt(p: f32, q: f32, c: f32, s: f32, dv: [f32; 2], v: f32) -> ([f32; 3], f32) {
    let ov2 = (1.0 / (v * v) as f64) as f32;
    let m11 = p / q;
    let m12 = -(dv[0] * c - dv[1] * s) * ov2;
    let m22 = -(dv[0] * s + dv[1] * c) * ov2;
    let d0 = ((m11 * c * c) as f64 + 2.0 * m12 as f64 * c as f64 * s as f64 + (m22 * s * s) as f64)
        as f32;
    let d1 = (m22 - m11) * c * s + m12 * (c * c - s * s);
    let d2 = ((m11 * s * s) as f64 - 2.0 * m12 as f64 * c as f64 * s as f64 + (m22 * c * c) as f64)
        as f32;
    let cuv = ((m12 as f64).abs() * v as f64) as f32;
    ([d0, d1, d2], cuv)
}

fn raytime2d(
    r: Raytr,
    g: Geo2d,
    v: &[f32],
    o: Geo2d,
    t: &mut [f32],
    npv: i32,
    vo: &[f32],
    pv: &[f32],
    tv: &mut [f32],
    cs: &mut [f32],
) {
    let nxo = o.nx;
    let nzo = o.nz;
    let mut vxx = vec![0.0f32; (g.nx * g.nz) as usize];
    let mut vxz = vxx.clone();
    let mut vzz = vxx.clone();
    let mut tvd = vec![0.0f32; r.nt.max(1) as usize];
    let mut s = vec![0.0f32; (nxo * nzo) as usize];
    let mut ray = vec![Ray::default(); r.nt.max(1) as usize];
    dv2(g, v, &mut vxx, &mut vxz, &mut vzz);
    let exo = o.fx + (nxo - 1) as f32 * o.dx;
    let ezo = o.fz + (nzo - 1) as f32 * o.dz;
    for ixo in 0..nxo {
        for izo in 0..nzo {
            let i = (izo + ixo * nzo) as usize;
            t[i] = INITIAL_T;
            s[i] = INITIAL_T;
            if npv != 0 {
                tv[i] = 0.0;
                cs[i] = 0.0;
            }
        }
    }
    let odt = (1.0 / r.dt as f64) as f32;
    let mut v0 = 0.0f32;
    let mut a = r.fa;
    for _ia in 0..r.na {
        let nrs = makeRay(g, v, &vxx, &vxz, &vzz, r, a, &mut ray, npv, pv);
        v0 = ray[0].v;
        let mut r2 = 0.0f32;
        let mut xc = r.xs;
        let mut zc = r.zs;
        let mut sd = 0.0f32;
        let mut vd1 = v0;
        if npv != 0 {
            let mut dtvd = 0.0f32;
            for it in 1..nrs as usize {
                dtvd = (dtvd as f64 + 0.5 * r.dt as f64 * (ray[it].dtv + ray[it - 1].dtv) as f64)
                    as f32;
                tvd[it] = ray[it].v * dtvd;
            }
        }
        for it in 1..nrs {
            let rr = ray[it as usize];
            let xo = rr.x;
            let zo = rr.z;
            let vd = rr.v;
            sd = (sd as f64 + 0.5 * r.dt as f64 * (vd + vd1) as f64) as f32;
            vd1 = vd;
            if r2 > (xc - xo) * (xc - xo) + (zc - zo) * (zc - zo) && it != nrs - 1 {
                continue;
            }
            if xo < o.fx || xo > exo || zo < o.fz || zo > ezo {
                continue;
            }
            xc = xo;
            zc = zo;
            let q2d = rr.q2;
            if q2d == 0.0 {
                r2 = 0.0;
                continue;
            }
            let vd = rr.v;
            let pxd = rr.px;
            let pzd = rr.pz;
            let p2d = rr.p2;
            let sind = vd * pxd;
            let cosd = vd * pzd;
            let (d2t, cuv) = ddt(p2d, q2d, cosd, sind, [rr.dvdx, rr.dvdz], vd);
            let tc = it as f32 * r.dt;
            let terr = tc * r.fac;
            let norm2 = (((d2t[0] * d2t[0] + d2t[2] * d2t[2]) as f64
                + 2.0 * d2t[1] as f64 * d2t[1] as f64)
                .sqrt()) as f32;
            r2 = terr / norm2;
            let mut r1 = (r2 as f64).sqrt() as f32;
            if (r1 * cuv) as f64 > 0.1 {
                r1 = (0.1 / cuv as f64) as f32;
            }
            if r1 as f64 > 0.1 * sd as f64 {
                r1 = (0.1 * sd as f64) as f32;
            }
            r2 = r1 * r1;
            let mut nxf = (((xc - r1 - o.fx) * o.odx) as f64 + 0.9) as i32;
            if nxf < 0 {
                nxf = 0;
            }
            let mut nxe = (((xc + r1 - o.fx) * o.odx) as f64 + 0.1) as i32;
            if nxe >= nxo {
                nxe = nxo - 1;
            }
            for ixo in nxf..=nxe {
                let xoc = o.fx - xc + ixo as f32 * o.dx;
                let r2x = r2 - xoc * xoc;
                if r2x < 0.0 {
                    continue;
                }
                let r1x = (r2x as f64).sqrt() as f32;
                let t1x = tc + xoc * pxd;
                let t2x = tc * xoc * xoc * d2t[0];
                let t2xz = (2.0 * xoc as f64 * d2t[1] as f64) as f32;
                let xcosd = if npv != 0 { pzd + xoc * d2t[1] } else { 0.0 };
                let mut nzf = (((zc - r1x - o.fz) * o.odz) as f64 + 0.9) as i32;
                if nzf < 0 {
                    nzf = 0;
                }
                let mut nze = (((zc + r1x - o.fz) * o.odz) as f64 + 0.1) as i32;
                if nze >= nzo {
                    nze = nzo - 1;
                }
                let i = (ixo * nzo) as usize;
                for izo in nzf..=nze {
                    let j = i + izo as usize;
                    if sd < s[j] {
                        let zoc = o.fz - zc + izo as f32 * o.dz;
                        let t1 = t1x + zoc * pzd;
                        let t2 = t2x + zoc * (zoc * d2t[2] + t2xz);
                        s[j] = sd;
                        t[j] = t1 * t1 + tc * t2;
                        if npv != 0 {
                            let mut jt = ((t1 as f64 + 0.5 * t2 as f64) * odt as f64 + 0.5) as i32;
                            if jt < 0 {
                                jt = 0;
                            }
                            if jt > nrs - 1 {
                                jt = nrs - 1;
                            }
                            tv[j] = tvd[jt as usize];
                            cs[j] = xcosd + d2t[2] * zoc;
                        }
                    }
                }
            }
        }
        a += r.da;
    }
    for ixo in 0..nxo {
        let i = (ixo * nzo) as usize;
        for izo in 0..nzo as usize {
            // MAX is a conditional expression, including its NaN behavior.
            t[i + izo] = if 0.0 > t[i + izo] { 0.0 } else { t[i + izo] };
            if t[i + izo] < 999.0 {
                t[i + izo] = (t[i + izo] as f64).sqrt() as f32;
            }
            if npv != 0 {
                cs[i + izo] *= vo[i + izo];
            }
        }
    }
    let mut nxf = (((r.xs - o.fx) / o.dx) as f64 - 1.5) as i32;
    if nxf < 0 {
        nxf = 0;
    }
    let mut nxe = (((r.xs - o.fx) / o.dx) as f64 + 2.5) as i32;
    if nxe > nxo - 1 {
        nxe = nxo - 1;
    }
    let mut nzf = (((r.zs - o.fz) / o.dz) as f64 - 0.5) as i32;
    if nzf < 0 {
        nzf = 0;
    }
    let mut nze = (((r.zs - o.fz) / o.dz) as f64 + 1.5) as i32;
    if nze >= nzo {
        nze = nzo - 1;
    }
    for ixo in nxf..=nxe {
        let xo = o.fx + ixo as f32 * o.dx;
        let i = (ixo * nzo) as usize;
        for izo in nzf..=nze {
            let zo = o.fz + izo as f32 * o.dz;
            let j = i + izo as usize;
            if t[j] > 999.0 {
                // The sqrt/v0 expression is double in C, with one final rounding.
                t[j] = ((((xo - r.xs) * (xo - r.xs) + (zo - r.zs) * (zo - r.zs)) as f64).sqrt()
                    / v0 as f64) as f32;
            }
        }
    }
}
