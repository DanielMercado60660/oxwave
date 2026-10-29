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

//! TRISEIS - Gaussian beam synthetic seismograms for a sloth model.
//! Port of tri/main/triseis.c (Test Release 1.1, 2011/11/21).
//! Author: Dave Hale, Colorado School of Mines, 02/09/91.
//! Modified: Andreas Rueger, Colorado School of Mines, 08/18/93:
//! 2.5-D amplitudes, ref/transm, timewindow, lscale, aperture, beam width.

#![allow(non_snake_case, unused_assignments, clippy::too_many_arguments)]
use oxwave_core::prelude::*;
use oxwave_cubic_interp::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::PI;
use oxwave_f_cubicspline::*;
use oxwave_inverse_interp::*;
use oxwave_lib_tri::*;
use oxwave_pfa_fft::*;
use std::ffi::c_void;
use std::ptr::null_mut;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " TRISEIS - Gaussian beam synthetic seismograms for a sloth model\t",
    "\t\t\t\t\t\t\t\t\t",
    "  triseis <modelfile >seisfile xs= zs= xg= zg= [optional parameters]\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " xs=            x coordinates of source surface\t\t\t",
    " zs=            z coordinates of source surface\t\t\t",
    " xg=            x coordinates of receiver surface\t\t\t",
    " zg=            z coordinates of receiver surface\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " ns=1           number of sources uniformly distributed along s surface",
    " ds=            increment between source locations (see notes below)\t",
    " fs=0.0         first source location (relative to start of s surface)\t",
    " ng=101         number of receivers uniformly distributed along g surface",
    " dg=            increment between receiver locations (see notes below)\t",
    " fg=0.0         first receiver location (relative to start of g surface)",
    " dgds=0.0       change in receiver location with respect to source location",
    " krecord=1      integer index of receiver surface (see notes below)\t",
    " kreflect=-1    integer index of reflecting surface (see notes below)\t",
    " prim           =1, only single-reflected rays are considered \t\t",
    "                =0, only direct hits are considered  \t\t\t",
    " bw=0           beamwidth at peak frequency\t\t\t\t",
    " nt=251         number of time samples\t\t\t\t\t",
    " dt=0.004       time sampling interval\t\t\t\t\t",
    " ft=0.0         first time sample\t\t\t\t\t",
    " nangle=101     number of ray takeoff angles\t\t\t\t",
    " fangle=-45     first ray takeoff angle (in degrees)\t\t\t",
    " langle=45      last ray takeoff angle (in degrees)\t\t\t",
    " reftrans=0     =1 complex refl/transm. coefficients considered \t",
    " atten=0        =1 add noncausal attenuation \t\t\t\t",
    "                =2 add causal attenuation \t\t\t\t",
    " lscale=        if defined restricts range of extrapolation\t\t",
    " fpeak=0.1/dt   peak frequency of ricker wavelet \t\t\t",
    " aperture=      maximum angle of receiver aperture \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " NOTES:\t\t\t\t\t\t\t\t",
    " Only rays that terminate with index krecord will contribute to the\t",
    " synthetic seismograms at the receiver (xg,zg) locations.  The\t\t",
    " source and receiver locations are determined by cubic spline\t\t",
    " interpolation of the specified (xs,zs) and (xg,zg) coordinates.\t",
    " The default source location increment (ds) is determined to span\t",
    " the source surface defined by (xs,zs).  Likewise for dg.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// Attribute layouts from Triangles/sloth.h; the library owns the raw graph.
#[repr(C)]
struct FaceAttributes {
    s00: f32,
    dsdx: f32,
    dsdz: f32,
    dens: f32,
    qfac: f32,
}
#[repr(C)]
struct EdgeAttributes {
    k: i32,
}
#[repr(C)]
struct EdgeUseAttributes {
    tx: f32,
    tz: f32,
    c: f32,
}
#[derive(Clone, Copy, Default)]
struct RayStep {
    sigma: f32,
    x: f32,
    z: f32,
    px: f32,
    pz: f32,
    t: f32,
    q1: f32,
    p1: f32,
    q2: f32,
    p2: f32,
    atten: f32,
    ampli: f32,
    ampliphase: f32,
    kmah: i32,
    nref: i32,
    eu: *mut EdgeUse,
    f: *mut Face,
}
#[derive(Clone, Copy, Default)]
struct RayEnd {
    rs: RayStep,
    sb: f32,
    se: f32,
    dsdxe: f32,
    dsdze: f32,
    kend: i32,
    dangle: f32,
}
extern "C" {
    static mut __stdinp: *mut c_void;
    fn __error() -> *mut i32;
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    // readModel takes a libc FILE*, whereas the core uses its own CFile.
    let m = unsafe { readModel(__stdinp) };
    let nxs = countparval("xs");
    let nzs = countparval("zs");
    if nxs == 0 || nzs == 0 {
        err!("must specify both xs and zs");
    }
    if nxs != nzs {
        err!("number of xs must equal number of zs");
    }
    let mut xsu = Vec::<f32>::new();
    let mut zsu = Vec::<f32>::new();
    getparfloat("xs", &mut xsu);
    getparfloat("zs", &mut zsu);
    let mut ns = 1i32;
    getparint("ns", &mut ns);
    let mut ds = 0.0f32;
    getparfloat("ds", &mut ds);
    let mut fs = 0.0f32;
    getparfloat("fs", &mut fs);
    let nxg = countparval("xg");
    let nzg = countparval("zg");
    if nxg == 0 || nzg == 0 {
        err!("must specify both xg and zg");
    }
    if nxg != nzg {
        err!("number of xg must equal number of zg");
    }
    let mut xgu = Vec::<f32>::new();
    let mut zgu = Vec::<f32>::new();
    getparfloat("xg", &mut xgu);
    getparfloat("zg", &mut zgu);
    let mut ng = 101i32;
    getparint("ng", &mut ng);
    let mut dg = 0.0f32;
    getparfloat("dg", &mut dg);
    let mut fg = 0.0f32;
    getparfloat("fg", &mut fg);
    let mut dgds = 0.0f32;
    getparfloat("dgds", &mut dgds);
    let mut krecord = 1i32;
    getparint("krecord", &mut krecord);
    let mut prim = i32::MAX;
    getparint("prim", &mut prim);
    let mut kreflect = -1i32;
    getparint("kreflect", &mut kreflect);
    let mut bw = 0.0f32;
    getparfloat("bw", &mut bw);
    let mut nt = 251i32;
    getparint("nt", &mut nt);
    let mut dt = 0.004f32;
    getparfloat("dt", &mut dt);
    let mut ft = 0.0f32;
    getparfloat("ft", &mut ft);
    let mut fpeak = (0.1 / dt as f64) as f32;
    getparfloat("fpeak", &mut fpeak);
    if 0.5 / (dt as f64) < 2.0 * fpeak as f64 {
        fprintf!(stderr(), "WARNING: ALIASING POSSIBLE\n");
        fprintf!(stderr(), "decrease dt or reduce fpeak\n");
    }
    let mut nangle = 101i32;
    getparint("nangle", &mut nangle);
    let mut fangle = -45.0f32;
    getparfloat("fangle", &mut fangle);
    let mut langle = 45.0f32;
    getparfloat("langle", &mut langle);
    let mut reftrans = 0i32;
    getparint("reftrans", &mut reftrans);
    let mut natten = 0i32;
    getparint("atten", &mut natten);
    let mut aperture = f32::MAX;
    getparfloat("aperture", &mut aperture);
    let mut lscale = f32::MAX;
    getparfloat("lscale", &mut lscale);
    checkpars();
    if aperture != f32::MAX {
        aperture = (aperture as f64 * PI / 180.0).asin() as f32;
    }
    fangle = (fangle as f64 * (PI / 180.0)) as f32;
    langle = (langle as f64 * (PI / 180.0)) as f32;
    let dangle = if nangle > 1 {
        (langle - fangle) / (nangle - 1) as f32
    } else {
        1.0
    };
    // The corresponding C size_t allocations cannot succeed for negative
    // dimensions. Report their allocator failure before converting to usize.
    if ns < 0 || ng < 0 || nt < 0 || nangle < 0 {
        unsafe {
            *__error() = 12;
        }
        syserr!("%s: malloc failed", "ealloc.c");
    }
    let mut xs = vec![0.0; ns as usize];
    let mut zs = vec![0.0; ns as usize];
    let mut xg = vec![0.0; ng as usize];
    let mut zg = vec![0.0; ng as usize];
    let zeros = vec![0.0f32; nt as usize];
    let mut re = vec![RayEnd::default(); nangle as usize];
    makexzs(&xsu, &zsu, ds, fs, &mut xs, &mut zs);
    let mut s = fs;
    for is in 0..ns as usize {
        fprintf!(stderr(), "shot %3d of %3d at x=%g\n", is as i32, ns, xs[is]);
        if !unsafe { insideModel(m, xs[is], zs[is]) } {
            for _ in 0..ng {
                if efwrite(&zeros, &stdout()) != nt as usize {
                    err!("error writing output traces!\n");
                }
            }
            s += ds;
            continue;
        }
        unsafe {
            shootRays(
                m, xs[is], zs[is], dangle, fangle, krecord, kreflect, &mut re,
            );
        }
        makexzs(&xgu, &zgu, dg, fg + dgds * (s - fs), &mut xg, &mut zg);
        makesyn(
            fpeak, prim, &re, lscale, krecord, &xg, &zg, bw, nt, dt, ft, natten, reftrans, aperture,
        );
        s += ds;
    }
    exit(CWP_Exit());
}

fn makexzs(xu: &[f32], zu: &[f32], mut ds: f32, fs: f32, xs: &mut [f32], zs: &mut [f32]) {
    let nu = xu.len();
    let ns = xs.len();
    let mut xud = vec![[0.0; 4]; nu];
    let mut zud = vec![[0.0; 4]; nu];
    let mut u = vec![0.0; nu];
    for iu in 0..nu {
        u[iu] = iu as f32;
    }
    csplin(&u, xu, &mut xud);
    csplin(&u, zu, &mut zud);
    let nuu = 10 * nu;
    let duu = (u[nu - 1] - u[0]) / (nuu - 1) as f32;
    // C reads s[-1] (heap garbage) on the first iteration; large ray
    // allocations can change it and mislocate every receiver. Deliberately
    // use zero, the intended value, for stable geometry.
    let mut s = vec![0.0f32; nuu];
    let mut xlast = xu[0];
    let mut zlast = zu[0];
    let mut uu = 0.0;
    for iuu in 0..nuu {
        let mut x = [0.0];
        let mut z = [0.0];
        intcub(0, &u, &xud, &[uu], &mut x);
        intcub(0, &u, &zud, &[uu], &mut z);
        let dx = x[0] - xlast;
        let dz = z[0] - zlast;
        s[iuu] = ((if iuu == 0 { 0.0 } else { s[iuu - 1] }) as f64
            + ((dx * dx + dz * dz) as f64).sqrt()) as f32;
        xlast = x[0];
        zlast = z[0];
        uu += duu;
    }
    let mut us = vec![0.0; ns];
    if ds == 0.0 {
        ds = (s[nuu - 1] - s[0]) / if ns > 1 { (ns - 1) as f32 } else { 1.0 };
    }
    if ds < 0.0 {
        yxtoxy(
            duu,
            0.0,
            &s,
            -ds,
            fs + (ns as i32 - 1) as f32 * ds,
            0.0,
            (nu - 1) as f32,
            &mut us,
        );
        us.reverse();
    } else {
        yxtoxy(duu, 0.0, &s, ds, fs, 0.0, (nu - 1) as f32, &mut us);
    }
    intcub(0, &u, &xud, &us, xs);
    intcub(0, &u, &zud, &us, zs);
}

fn makesyn(
    fpeak: f32,
    prim: i32,
    re: &[RayEnd],
    lscale: f32,
    krecord: i32,
    xg: &[f32],
    zg: &[f32],
    bw: f32,
    nt: i32,
    dt: f32,
    ft: f32,
    natten: i32,
    reftrans: i32,
    aperture: f32,
) {
    let ntfft = npfaro(nt, 2 * nt);
    let nw = (ntfft / 2 + 1) as usize;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let wpeak = (2.0 * PI * fpeak as f64) as f32;
    let wref = 1000.0f32;
    let delay = 0.0f32;
    let sqlscale = lscale * lscale;
    let c = 0.99f32;
    let eps = 0.0001f32;
    let mut syn = vec![0.0f32; ntfft as usize];
    let mut cwave = vec![cmplx(0.0, 0.0); nw];
    let mut csyn = vec![vec![cmplx(0.0, 0.0); nw]; xg.len()];
    let mut lnvec = vec![0.0f32; nw];
    if natten == 2 {
        let temp = (1.0 / PI) as f32;
        lnvec[0] = -1.0;
        let mut w = dw;
        for iw in 1..nw {
            lnvec[iw] = (((w / wref) as f64).ln() * temp as f64) as f32;
            w += dw;
        }
    }
    for end in re {
        if end.kend != krecord {
            continue;
        }
        let r = end.rs;
        let nref = r.nref;
        if nref != prim && prim != i32::MAX {
            continue;
        }
        let px = r.px;
        let v = (1.0 / (end.se as f64).sqrt()) as f32;
        if (px * v).abs() > aperture {
            continue;
        }
        let sigma = (1.0 / (r.sigma as f64).sqrt()) as f32;
        let (x, z, t, pz) = (r.x, r.z, r.t, r.pz);
        let (q1, p1, q2, p2) = (r.q1, r.p1, r.q2, r.p2);
        let (kmah, atten, ampli, ampliphase) = (r.kmah, r.atten, r.ampli, r.ampliphase);
        let dangle = end.dangle;
        let (dsdx, dsdz) = (end.dsdxe, end.dsdze);
        let v0 = (1.0 / (end.sb as f64).sqrt()) as f32;
        let v2 = v * v;
        let v4 = v2 * v2;
        let dvds = (-0.5 * v4 as f64 * (dsdx * px + dsdz * pz) as f64) as f32;
        let dvdn = (-0.5 * v4 as f64 * (dsdx * pz - dsdz * px) as f64) as f32;
        let scale = ((sigma * dangle) as f64 * (v as f64 / ((v0 * 2.0) as f64 * PI)).sqrt()) as f32;
        let (eps1, eps2);
        if bw != 0.0 {
            eps1 = 0.0;
            eps2 = (-0.5 * wpeak as f64 * bw as f64 * bw as f64) as f32;
        } else {
            let temp = c * (p2 / q2).abs() + eps;
            let e1 = -q2 / q1 * (p2 * p1 / (q2 * q1) + temp * temp);
            eps1 = e1 / (p1 * p1 / (q1 * q1) + temp * temp);
            eps2 = -temp / (p1 * p1 + q1 * q1 * temp * temp);
        }
        let ceps = cmplx(eps1, eps2);
        let mut cp = crmul(ceps, p1);
        cp.r += p2;
        let mut cq = crmul(ceps, q1);
        cq.r += q2;
        let cpoq = cdiv(cp, cq);
        let (cpoqr, cpoqi) = (cpoq.r, cpoq.i);
        if cpoqi < 0.0 {
            continue;
        }
        let cqr = if nref % 2 != 0 { cneg(cq) } else { cq };
        let mut phase = (-PI * ((kmah + 1) / 2) as f64) as f32;
        let mut camp = cwp_csqrt(cdiv(ceps, cqr));
        if reftrans != 0 {
            phase += ampliphase;
            camp = crmul(camp, ampli);
        }
        phase = (phase as f64 + PI / 4.0) as f32;
        camp = crmul(camp, scale);
        camp = cmul(
            camp,
            cmplx((phase as f64).cos() as f32, (phase as f64).sin() as f32),
        );
        let campr = camp.r as f64;
        let campi = camp.i as f64;
        if natten == 2 {
            for iw in 0..nw {
                lnvec[iw] = atten * lnvec[iw];
            }
        }
        for ig in 0..xg.len() {
            let dx = xg[ig] - x;
            let dz = zg[ig] - z;
            let disq = dx * dx + dz * dz;
            let n = v * (pz * dx - px * dz);
            let nn = n * n;
            if lscale != f32::MAX && disq > sqlscale {
                continue;
            }
            let dels = v * (px * dx + pz * dz) * (1.0 - dvdn / v * n);
            let delt = ((dels / v) as f64
                - 0.5 * dvds as f64 * dels as f64 * dels as f64 / v2 as f64)
                as f32;
            let ta = ((t + delt) as f64 + 0.5 * cpoqr as f64 * nn as f64 - ft as f64) as f32;
            if ta < 0.0 || ta > (nt - 1) as f32 * dt {
                continue;
            }
            let cosd = ((dw * ta) as f64).cos();
            let sind = ((dw * ta) as f64).sin();
            let expd = if natten == 0 {
                (-0.5 * dw as f64 * cpoqi as f64 * nn as f64).exp()
            } else {
                (-0.5 * dw as f64 * (cpoqi * nn + atten) as f64).exp()
            };
            let mut cosw = 1.0f64;
            let mut sinw = 0.0f64;
            let mut expw = 1.0f64;
            if natten == 0 || natten == 1 {
                for iw in 0..nw {
                    csyn[ig][iw].r =
                        (csyn[ig][iw].r as f64 + expw * (campr * cosw - campi * sinw)) as f32;
                    csyn[ig][iw].i =
                        (csyn[ig][iw].i as f64 + expw * (campr * sinw + campi * cosw)) as f32;
                    let tempd = cosw * cosd - sinw * sind;
                    sinw = cosw * sind + sinw * cosd;
                    cosw = tempd;
                    expw *= expd;
                    if expw < 0.01 {
                        break;
                    }
                }
            } else {
                let mut w = 0.0;
                for iw in 0..nw {
                    let tempd = (w * (ta - lnvec[iw])) as f64;
                    cosw = tempd.cos();
                    sinw = tempd.sin();
                    csyn[ig][iw].r =
                        (csyn[ig][iw].r as f64 + expw * (campr * cosw - campi * sinw)) as f32;
                    csyn[ig][iw].i =
                        (csyn[ig][iw].i as f64 + expw * (campr * sinw + campi * cosw)) as f32;
                    expw *= expd;
                    if expw < 0.01 {
                        break;
                    }
                    w += dw;
                }
            }
        }
    }
    let mut w = 0.0;
    for iw in 0..nw {
        cwave[iw] = cricker(w, wpeak, delay);
        cwave[iw] = crmul(cwave[iw], (w as f64).sqrt() as f32);
        w += dw;
    }
    for ig in 0..xg.len() {
        for iw in 0..nw {
            csyn[ig][iw] = cmul(csyn[ig][iw], cwave[iw]);
        }
        pfacr(-1, &csyn[ig], &mut syn);
        fwrite(&syn[..nt as usize], &stdout());
    }
}

fn cricker(w: f32, wpeak: f32, delay: f32) -> Complex {
    crmul(
        cwp_cexp(cmplx(-((w / wpeak) as f64).powf(2.0) as f32, delay * w)),
        (4.0 * w as f64 * w as f64 * PI.sqrt() / (wpeak as f64).powf(3.0)) as f32,
    )
}

unsafe fn insideModel(m: *mut Model, x: f32, z: f32) -> bool {
    (*m).ymin <= x && x <= (*m).ymax && (*m).xmin <= z && z <= (*m).xmax
}

unsafe fn shootRays(
    m: *mut Model,
    mut xs: f32,
    mut zs: f32,
    dangle: f32,
    fangle: f32,
    kstop: i32,
    kreflect: i32,
    re: &mut [RayEnd],
) {
    let mut tris = insideTriInModel(m, null_mut(), zs, xs);
    let dec = checkIfSourceOnEdge(tris, zs, xs);
    if dec != 0 {
        let zmax = (*m).xmax;
        let mut eps = (*m).eps;
        if dec == 2 {
            eps *= 5.0;
        }
        zs = if zs < zmax - eps { zs + eps } else { zs - eps };
        if checkIfSourceOnEdge(tris, zs, xs) != 0 {
            let xmax = (*m).ymax;
            xs = if xs < xmax - eps { xs + eps } else { xs - eps };
        }
        tris = insideTriInModel(m, null_mut(), zs, xs);
    }
    let fa = &*((*tris).fa as *mut FaceAttributes);
    let sb = fa.s00 + fa.dsdx * xs + fa.dsdz * zs;
    let mut angle = fangle;
    for end in re {
        let mut rs = RayStep {
            sigma: 0.0,
            x: xs,
            z: zs,
            px: ((angle as f64).sin() * (sb as f64).sqrt()) as f32,
            pz: ((angle as f64).cos() * (sb as f64).sqrt()) as f32,
            q1: 1.0,
            p1: 0.0,
            q2: 0.0,
            p2: 1.0,
            t: 0.0,
            kmah: 0,
            nref: 0,
            eu: null_mut(),
            f: tris,
            ampli: 1.0,
            ampliphase: 0.0,
            atten: 0.0,
        };
        let stopped;
        loop {
            rs = traceRayInTri(rs);
            let ea = (*(*rs.eu).e).ea as *mut EdgeAttributes;
            if !ea.is_null() && (*ea).k == kstop {
                stopped = true;
                break;
            } else if !ea.is_null() && (*ea).k == kreflect {
                rs = reflectRayFromEdge(rs);
            } else if let Some(next) = traceRayAcrossEdge(rs) {
                rs = next;
            } else {
                stopped = false;
                break;
            }
        }
        let fa = &*((*rs.f).fa as *mut FaceAttributes);
        end.rs = rs;
        if rs.ampliphase < -3.0 && kstop == 1 && rs.nref == 0 {
            fprintf!(stderr(), "phase=%f    ampli=%f\n", rs.ampliphase, rs.ampli);
        }
        end.sb = sb;
        end.se = fa.s00 + fa.dsdx * rs.x + fa.dsdz * rs.z;
        end.dsdxe = fa.dsdx;
        end.dsdze = fa.dsdz;
        end.kend = if stopped { kstop } else { -1 };
        end.dangle = dangle;
        angle += dangle;
    }
}

unsafe fn traceRayInTri(mut rs: RayStep) -> RayStep {
    let (mut x, mut z, mut px, mut pz) = (rs.x, rs.z, rs.px, rs.pz);
    let fa = &*((*rs.f).fa as *mut FaceAttributes);
    let (s00, dsdx, dsdz, qfac) = (fa.s00, fa.dsdx, fa.dsdz, fa.qfac);
    let mut eusmall: *mut EdgeUse = null_mut();
    let mut euzero: *mut EdgeUse = null_mut();
    let mut small = f32::MAX;
    let mut eut = (*rs.f).eu;
    loop {
        let va = &*(*(*eut).vu).v;
        let vb = &*(*(*(*eut).euCW).vu).v;
        let (xa, za, xb, zb) = (va.y, va.x, vb.y, vb.x);
        let dx = xb - xa;
        let dz = zb - za;
        let b = dx * pz - dz * px;
        let c = if eut == rs.eu {
            0.0
        } else {
            dx * (z - za) - dz * (x - xa)
        };
        let (dsigma1, dsigma2);
        if dsdx == 0.0 && dsdz == 0.0 {
            dsigma1 = if b != 0.0 { -c / b } else { f32::MAX };
            dsigma2 = f32::MAX;
        } else {
            let a = (0.25 * (dx * dsdz - dz * dsdx) as f64) as f32;
            let ds = ((b * b) as f64 - 4.0 * a as f64 * c as f64) as f32;
            if ds < 0.0 {
                dsigma1 = f32::MAX;
                dsigma2 = f32::MAX;
            } else {
                let q = (-0.5
                    * if b < 0.0 {
                        b as f64 - (ds as f64).sqrt()
                    } else {
                        b as f64 + (ds as f64).sqrt()
                    }) as f32;
                dsigma1 = if a != 0.0 { q / a } else { f32::MAX };
                dsigma2 = if q != 0.0 { c / q } else { f32::MAX };
            }
        }
        if 0.0 < dsigma1 && dsigma1 < small {
            small = dsigma1;
            eusmall = eut;
        }
        if 0.0 < dsigma2 && dsigma2 < small {
            small = dsigma2;
            eusmall = eut;
        }
        if dsigma1 == 0.0 {
            euzero = eut;
        }
        if dsigma2 == 0.0 {
            euzero = eut;
        }
        eut = (*eut).euCW;
        if eut == (*rs.f).eu {
            break;
        }
    }
    let dsigma;
    if !eusmall.is_null() {
        dsigma = small;
        rs.eu = eusmall;
    } else {
        dsigma = 0.0;
        rs.eu = euzero;
    }
    evaluateDynamic(
        dsigma,
        px,
        pz,
        dsdx,
        dsdz,
        &mut rs.q1,
        &mut rs.p1,
        &mut rs.q2,
        &mut rs.p2,
        &mut rs.kmah,
    );
    rs.sigma += dsigma;
    let ddt = (dsigma as f64
        * ((s00 + dsdx * x + dsdz * z) as f64
            + dsigma as f64
                * (0.5 * (dsdx * px + dsdz * pz) as f64
                    + dsigma as f64 * (0.0833333 * (dsdx * dsdx + dsdz * dsdz) as f64))))
        as f32;
    rs.t += ddt;
    x = (x as f64 + dsigma as f64 * (px as f64 + 0.25 * dsdx as f64 * dsigma as f64)) as f32;
    z = (z as f64 + dsigma as f64 * (pz as f64 + 0.25 * dsdz as f64 * dsigma as f64)) as f32;
    px = (px as f64 + 0.5 * dsdx as f64 * dsigma as f64) as f32;
    pz = (pz as f64 + 0.5 * dsdz as f64 * dsigma as f64) as f32;
    let va = &*(*(*rs.eu).vu).v;
    let vb = &*(*(*(*rs.eu).euCW).vu).v;
    let (xa, za) = (va.y, va.x);
    let dx = vb.y - xa;
    let dz = vb.x - za;
    let frac = if dx.abs() > dz.abs() {
        (x - xa) / dx
    } else {
        (z - za) / dz
    };
    if (frac as f64) < 0.0001 {
        x = (xa as f64 + 0.0001 * dx as f64) as f32;
        z = (za as f64 + 0.0001 * dz as f64) as f32;
    } else if (frac as f64) > 0.9999 {
        x = (xa as f64 + 0.9999 * dx as f64) as f32;
        z = (za as f64 + 0.9999 * dz as f64) as f32;
    }
    if qfac != f32::MAX {
        rs.atten += ddt / qfac;
    } else {
        rs.atten = 0.0;
    }
    rs.x = x;
    rs.z = z;
    rs.px = px;
    rs.pz = pz;
    rs
}

unsafe fn traceRayAcrossEdge(mut rs: RayStep) -> Option<RayStep> {
    let (x, z, mut px, mut pz) = (rs.x, rs.z, rs.px, rs.pz);
    let eu = rs.eu;
    if (*(*eu).euMate).f.is_null() {
        return None;
    }
    let fa = &*((*rs.f).fa as *mut FaceAttributes);
    let (ds1dx, ds1dz) = (fa.dsdx, fa.dsdz);
    let s1 = fa.s00 + ds1dx * x + ds1dz * z;
    let mut dens1 = fa.dens;
    let eum = (*eu).euMate;
    let f = (*eum).f;
    let fa = &*((*f).fa as *mut FaceAttributes);
    let (ds2dx, ds2dz) = (fa.dsdx, fa.dsdz);
    let s2 = fa.s00 + ds2dx * x + ds2dz * z;
    let mut dens2 = fa.dens;
    if dens1 == f32::MAX {
        dens1 = 1.0;
    }
    if dens2 == f32::MAX {
        dens2 = 1.0;
    }
    if s1 != s2 || ds1dx != ds2dx || ds1dz != ds2dz || dens1 != dens2 {
        let va = &*(*(*eu).vu).v;
        let vb = &*(*(*eum).vu).v;
        let dx = vb.y - va.y;
        let dz = vb.x - va.x;
        let frac = if dx.abs() > dz.abs() {
            (x - va.y) / dx
        } else {
            (z - va.x) / dz
        };
        let eua = (*eu).eua as *mut EdgeUseAttributes;
        let euma = (*eum).eua as *mut EdgeUseAttributes;
        let (mut gx, mut gz);
        if !eua.is_null() && !euma.is_null() {
            gx = ((frac * (*euma).tx) as f64 - (1.0 - frac as f64) * (*eua).tx as f64) as f32;
            gz = ((frac * (*euma).tz) as f64 - (1.0 - frac as f64) * (*eua).tz as f64) as f32;
        } else {
            gx = -dx;
            gz = -dz;
        }
        let scale = (1.0 / ((gx * gx + gz * gz) as f64).sqrt()) as f32;
        gx *= scale;
        gz *= scale;
        let hx = -gz;
        let h_z = gx;
        let (px1, pz1) = (px, pz);
        let px1r = px * h_z - pz * hx;
        let pz1r = px * hx + pz * h_z;
        let px2r = px1r;
        let pz2rs = s2 - px2r * px2r;
        if pz2rs <= 0.0 {
            return None;
        }
        if (pz1r * pz1r) as f64 <= 0.008 * s1 as f64 {
            return None;
        }
        let pz2r = (pz2rs as f64).sqrt() as f32;
        px = px2r * h_z + pz2r * hx;
        pz = pz2r * h_z - px2r * hx;
        let c1ov1 = pz1r;
        let c2ov2 = pz2r;
        let oc2s = s2 / pz2rs;
        let cterm;
        if !eua.is_null() && !euma.is_null() {
            let g = ((frac * (*euma).c) as f64 - (1.0 - frac as f64) * (*eua).c as f64) as f32;
            cterm = g * oc2s * (c1ov1 - c2ov2);
        } else {
            cterm = 0.0;
        }
        let scale = (pz2r as f64 * (s1 as f64).sqrt() / (pz1r as f64 * (s2 as f64).sqrt())) as f32;
        rs.q1 = rs.q1 * scale;
        rs.p1 = rs.p1 / scale + cterm * rs.q1;
        rs.q2 = rs.q2 * scale;
        rs.p2 = rs.p2 / scale + cterm * rs.q2;
        let fac1 = (-0.5 / (s1 * s1) as f64) as f32;
        let fac2 = (-0.5 / (s2 * s2) as f64) as f32;
        let dv1dl = fac1 * (px1 * ds1dx + pz1 * ds1dz);
        let dv2dl = fac2 * (px * ds2dx + pz * ds2dz);
        let dv1dm = fac1 * (pz1 * ds1dx - px1 * ds1dz);
        let dv2dm = fac2 * (pz * ds2dx - px * ds2dz);
        let iterm = ((-px1r * oc2s) as f64
            * (2.0 * (dv1dm * c1ov1 - dv2dm * c2ov2) as f64 + (px1r * (dv1dl - dv2dl)) as f64))
            as f32;
        rs.p1 += iterm * rs.q1;
        rs.p2 += iterm * rs.q2;
        if s1 != s2 || dens1 != dens2 {
            let taper = (px1r * px1r) / s2;
            let taper = if taper as f64 > 0.97 {
                ((taper as f64 - 0.97) * 5.0 * PI).cos()
            } else {
                1.0
            } as f32;
            let coeff1 = dens2 / dens1;
            let coeff2 = pz2r / pz1r;
            let coeff = (1.0 + ((coeff1 - coeff2) / (coeff1 + coeff2)) as f64) as f32;
            rs.ampli = (rs.ampli as f64
                * (coeff as f64 * ((pz2r / pz1r) as f64).sqrt() * taper as f64))
                as f32;
        }
    }
    rs.px = px;
    rs.pz = pz;
    rs.eu = eum;
    rs.f = f;
    Some(rs)
}

unsafe fn reflectRayFromEdge(mut rs: RayStep) -> RayStep {
    let (x, z, mut px, mut pz) = (rs.x, rs.z, rs.px, rs.pz);
    let eu = rs.eu;
    let fa = &*((*rs.f).fa as *mut FaceAttributes);
    let (dsdx, dsdz) = (fa.dsdx, fa.dsdz);
    let s = fa.s00 + dsdx * x + dsdz * z;
    let mut dens1 = fa.dens;
    let eum = (*eu).euMate;
    let fa = &*((*(*eum).f).fa as *mut FaceAttributes);
    let s2 = fa.s00 + fa.dsdx * x + fa.dsdz * z;
    let mut dens2 = fa.dens;
    if dens1 == f32::MAX {
        dens1 = 1.0;
    }
    if dens2 == f32::MAX {
        dens2 = 1.0;
    }
    let va = &*(*(*eu).vu).v;
    let vb = &*(*(*eum).vu).v;
    let dx = vb.y - va.y;
    let dz = vb.x - va.x;
    let frac = if dx.abs() > dz.abs() {
        (x - va.y) / dx
    } else {
        (z - va.x) / dz
    };
    let eua = &*((*eu).eua as *mut EdgeUseAttributes);
    let euma = &*((*eum).eua as *mut EdgeUseAttributes);
    let mut gx = ((frac * euma.tx) as f64 - (1.0 - frac as f64) * eua.tx as f64) as f32;
    let mut gz = ((frac * euma.tz) as f64 - (1.0 - frac as f64) * eua.tz as f64) as f32;
    let scale = (1.0 / ((gx * gx + gz * gz) as f64).sqrt()) as f32;
    gx *= scale;
    gz *= scale;
    let hx = -gz;
    let h_z = gx;
    let (px1, pz1) = (px, pz);
    let pxr = px * h_z - pz * hx;
    let mut pzr = px * hx + pz * h_z;
    pzr = -pzr;
    px = pxr * h_z + pzr * hx;
    pz = pzr * h_z - pxr * hx;
    let c1ov1 = -pzr;
    let c2ov2 = -pzr;
    let oc2s = s / (pzr * pzr);
    let g = ((frac * euma.c) as f64 - (1.0 - frac as f64) * eua.c as f64) as f32;
    let cterm = g * oc2s * (c1ov1 + c2ov2);
    let scale = -c2ov2 / c1ov1;
    rs.q1 = rs.q1 * scale;
    rs.p1 = rs.p1 / scale + cterm * rs.q1;
    rs.q2 = rs.q2 * scale;
    rs.p2 = rs.p2 / scale + cterm * rs.q2;
    if dsdx != 0.0 || dsdz != 0.0 {
        let scale = (-0.5 / (s * s) as f64) as f32;
        let dv1dl = scale * (px1 * dsdx + pz1 * dsdz);
        let dv2dl = scale * (px * dsdx + pz * dsdz);
        let dv1dm = scale * (pz1 * dsdx - px1 * dsdz);
        let dv2dm = scale * (pz * dsdx - px * dsdz);
        let iterm = ((-pxr * oc2s) as f64
            * (2.0 * (dv1dm * c1ov1 + dv2dm * c2ov2) as f64 + (pxr * (dv1dl - dv2dl)) as f64))
            as f32;
        rs.p1 += iterm * rs.q1;
        rs.p2 += iterm * rs.q2;
    }
    rs.nref += 1;
    if s2 - pxr * pxr >= 0.0 {
        let temp1 = ((-pzr * dens2 / dens1) as f64 / (s as f64).sqrt()) as f32;
        let temp2 = ((s2 / s - (pxr * pxr / s)) as f64).sqrt() as f32;
        rs.ampli *= (temp1 - temp2) / (temp1 + temp2);
    } else {
        rs.ampliphase = (rs.ampliphase as f64
            - 2.0
                * (-((pxr * pxr - s2) as f64).sqrt() * dens1 as f64 / (dens2 * pzr) as f64).atan())
            as f32;
    }
    rs.px = px;
    rs.pz = pz;
    rs
}

fn evaluateDynamic(
    sigma: f32,
    px: f32,
    pz: f32,
    dsdx: f32,
    dsdz: f32,
    q1: &mut f32,
    p1: &mut f32,
    q2: &mut f32,
    p2: &mut f32,
    kmah: &mut i32,
) {
    let (q1i, p1i, q2i, p2i) = (*q1 as f64, *p1 as f64, *q2 as f64, *p2 as f64);
    let (q1o, p1o, q2o, p2o);
    if pz * dsdx == px * dsdz {
        q1o = q1i + p1i * sigma as f64;
        q2o = q2i + p2i * sigma as f64;
        p2o = p2i;
        p1o = p1i;
    } else {
        let s0 = (px * px + pz * pz) as f64;
        let s1 = (px * dsdx + pz * dsdz) as f64;
        let s2 = 0.25 * (dsdx * dsdx + dsdz * dsdz) as f64;
        let sigma = sigma as f64;
        let ss = (s2 * sigma + s1) * sigma + s0;
        let scale = 1.0 / (s0 * ss).sqrt();
        let a = s0 + 0.5 * s1 * sigma;
        let mut b = (0.25 * s1 * s1 / s0 - s2) * sigma;
        let c = s0 + s1 * sigma;
        let d = 0.5 * s1 + 2.0 * b;
        b *= sigma;
        let e = (0.5 * s1 + s2 * sigma) / ss;
        q1o = scale * (a * (q1i + p1i * sigma) + b * q1i);
        q2o = scale * (a * (q2i + p2i * sigma) + b * q2i);
        p1o = scale * (c * p1i + d * q1i) - e * q1o;
        p2o = scale * (c * p2i + d * q2i) - e * q2o;
    }
    if q2i * q2o >= 0.0 && p2i * p2o < 0.0 && q2i * p2i < 0.0 {
        *kmah += 2;
    } else if q2o == 0.0 || q2i * q2o < 0.0 || (q2i == 0.0 && p2i * p2o < 0.0 && q2o * p2o > 0.0) {
        *kmah += 1;
    }
    *q1 = q1o as f32;
    *p1 = p1o as f32;
    *q2 = q2o as f32;
    *p2 = p2o as f32;
}

unsafe fn checkIfSourceOnEdge(tris: *mut Face, zs: f32, xs: f32) -> i32 {
    let eu = (*tris).eu;
    let eps = (*(*tris).m).eps;
    let v1 = &*(*(*eu).vu).v;
    let v2 = &*(*(*(*eu).euCW).vu).v;
    let v3 = &*(*(*(*eu).euCCW).vu).v;
    let (x1, z1, x2, z2, x3, z3) = (v1.y, v1.x, v2.y, v2.x, v3.y, v3.x);
    if (xs - x1) * (xs - x1) + (zs - z1) * (zs - z1) < eps
        || (xs - x2) * (xs - x2) + (zs - z2) * (zs - z2) < eps
        || (xs - x3) * (xs - x3) + (zs - z3) * (zs - z3) < eps
    {
        return 2;
    }
    let (m12, m13, m23);
    if ((x1 - x2).abs() as f64) < 0.1 * eps as f64 {
        if ((xs - x1).abs() as f64) < 0.1 * eps as f64 {
            return 1;
        } else {
            m12 = (0.01 * f32::MAX as f64) as f32;
        }
    } else {
        m12 = (z1 - z2) / (x1 - x2);
    }
    if ((x1 - x3).abs() as f64) < 0.1 * eps as f64 {
        if ((xs - x1).abs() as f64) < 0.1 * eps as f64 {
            return 1;
        } else {
            m13 = (0.01 * f32::MAX as f64) as f32;
        }
    } else {
        m13 = (z1 - z3) / (x1 - x3);
    }
    if ((x2 - x3).abs() as f64) < 0.1 * eps as f64 {
        if ((xs - x2).abs() as f64) < 0.1 * eps as f64 {
            return 1;
        } else {
            m23 = (0.01 * f32::MAX as f64) as f32;
        }
    } else {
        m23 = (z2 - z3) / (x2 - x3);
    }
    let b12 = z1 - m12 * x1;
    let b13 = z1 - m13 * x1;
    let b23 = z2 - m23 * x2;
    if ((zs - m12 * xs - b12).abs() as f64) < 0.1 * eps as f64
        || ((zs - m13 * xs - b13).abs() as f64) < 0.1 * eps as f64
        || ((zs - m23 * xs - b23).abs() as f64) < 0.1 * eps as f64
    {
        return 1;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_receiver_positions_with_130_rayends() {
        // triseis generates its rayends internally. Exercise the main-loop
        // order with a real >94-record ray fan before computing receivers;
        // the C s[-1] read makes a large-input golden an unreliable oracle.
        let m = makeModel(0.0, -1.0, 2.0, 5.0);
        let mut fa = FaceAttributes {
            s00: 1.0,
            dsdx: 0.0,
            dsdz: 0.0,
            dens: 1.0,
            qfac: f32::MAX,
        };
        let fa_ptr = (&mut fa as *mut FaceAttributes).cast();
        let mut re = vec![RayEnd::default(); 130];
        let fangle = (-65.0 * PI / 180.0) as f32;
        let langle = (65.0 * PI / 180.0) as f32;
        let dangle = (langle - fangle) / (re.len() - 1) as f32;
        unsafe {
            let mut f = (*m).f;
            loop {
                (*f).fa = fa_ptr;
                f = (*f).fNext;
                if f == (*m).f {
                    break;
                }
            }
            shootRays(m, 1.5, 0.1, dangle, fangle, 1, -1, &mut re);
        }

        let mut xg = [0.0; 7];
        let mut zg = [0.0; 7];
        makexzs(&[0.0, 4.0], &[0.0, 0.0], 0.0, 0.0, &mut xg, &mut zg);

        // killModel frees the topology, but does not own the face attributes.
        unsafe { killModel(m) };
        assert!(re.iter().all(|end| end.rs.t > 0.0));
        for ig in 0..7 {
            let expected = 4.0 * ig as f32 / 6.0;
            assert!(
                (xg[ig] - expected).abs() <= 1e-6,
                "receiver {ig}: x={}, expected {expected}",
                xg[ig]
            );
            assert_eq!(zg[ig], 0.0);
        }
    }
}
