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

//! GBBEAM - Gaussian beam synthetic seismograms for a sloth model.
//! Port of tri/main/gbbeam.c (Test Release 1.1, 2011/11/21).
//! Author: Dave Hale, Colorado School of Mines, 02/09/91.
//! Modified: Andreas Rueger, Colorado School of Mines, 08/18/93:
//! 2.5-D amplitudes, reflection/transmission losses, attenuation,
//! timewindow, lscale, aperture, beam width, etc.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::{abs, PI};
use oxwave_f_cubicspline::csplin;
use oxwave_inverse_interp::yxtoxy;
use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::{npfaro, pfacr};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " GBBEAM - Gaussian beam synthetic seismograms for a sloth model \t",
    "\t\t\t\t\t\t\t\t\t",
    " gbbeam <rayends >syntraces xg= zg= [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " xg=              x coordinates of receiver surface\t\t\t",
    " zg=              z coordinates of receiver surface\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " ng=101           number of receivers (uniform distributed along surface)",
    " krecord=1        integer index of receiver surface (see notes below)\t",
    " ang=0.0          array of angles corresponding to amplitudes in amp\t",
    " amp=1.0          array of amplitudes corresponding to angles in ang\t",
    " bw=0             beamwidth at peak frequency\t\t\t\t",
    " nt=251           number of time samples\t\t\t\t",
    " dt=0.004         time sampling interval\t\t\t\t",
    " ft=0.0           first time sample\t\t\t\t\t",
    " reftrans=0       =1 complex refl/transm. coefficients considered\t",
    " prim             =1, only single-reflected rays are considered\t",
    "                  =0, only direct hits are considered\t\t\t",
    " atten=0          =1 add noncausal attenuation\t\t\t\t",
    "                  =2 add causal attenuation\t\t\t\t",
    " lscale=          if defined restricts range of extrapolation\t\t",
    " aperture=        maximum angle of receiver aperture\t\t\t",
    " fpeak=0.1/dt     peak frequency of ricker wavelet\t\t\t",
    " infofile         ASCII-file to store useful information\t\t",
    " NOTES:\t\t\t\t\t\t\t\t",
    " Only rays that terminate with index krecord will contribute to the\t",
    " synthetic seismograms at the receiver (xg,zg) locations.  The\t\t",
    " receiver locations are determined by cubic spline interpolation\t",
    " of the specified (xg,zg) coordinates.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

/// Native layout from Triangles/sloth.h; the input is not an SU trace stream.
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct RayEnd {
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
    kmah: i32,
    nref: i32,
    sb: f32,
    se: f32,
    dsdxe: f32,
    dsdze: f32,
    ab: f32,
    kend: i32,
    ampli: f32,
    ampliphase: f32,
    atten: f32,
    dangle: f32,
}

const _: () = assert!(std::mem::size_of::<RayEnd>() == 88);
// SAFETY: repr(C), exclusively four-byte floats and integers, without padding;
// all bit patterns are valid, just as for the C fread of RayEnd.
unsafe impl oxwave_core::cio::Pod for RayEnd {}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    let nxg = countparval("xg");
    let nzg = countparval("zg");
    if nxg == 0 || nzg == 0 {
        err!("must specify both xg and zg");
    }
    if nxg != nzg {
        err!("number of xg must equal number of zg");
    }
    let nxz = nxg as usize;
    let mut ng = 101i32;
    getparint("ng", &mut ng);
    let mut xg = vec![0.0f32; nxz.max(ng.max(0) as usize)];
    let mut zg = vec![0.0f32; nxz.max(ng.max(0) as usize)];
    getparfloat("xg", &mut xg);
    getparfloat("zg", &mut zg);

    let mut infofile = String::new();
    let infofp = if getparstring("infofile", &mut infofile) != 0 {
        Some(efopen(&infofile, "w"))
    } else {
        None
    };
    let mut lscale = f32::MAX;
    let mut bw = 0.0f32;
    let mut krecord = 1i32;
    let mut nt = 251i32;
    let mut dt = 0.004f32;
    let mut ft = 0.0f32;
    let mut reftrans = 0i32;
    let mut natten = 0i32;
    let mut aperture = f32::MAX;
    getparfloat("lscale", &mut lscale);
    getparfloat("bw", &mut bw);
    getparint("krecord", &mut krecord);
    getparint("nt", &mut nt);
    if nt < 0 {
        // C converts nt to size_t in efwrite and reads beyond syn. Reject
        // that undefined case before FFT allocation or output slicing.
        err!("nt must be nonnegative");
    }
    getparfloat("dt", &mut dt);
    getparfloat("ft", &mut ft);
    getparint("reftrans", &mut reftrans);
    getparint("atten", &mut natten);
    getparfloat("aperture", &mut aperture);
    if aperture != f32::MAX {
        aperture = (aperture as f64 * PI / 180.0).asin() as f32;
    }
    let mut fpeak = 0.0f32;
    if getparfloat("fpeak", &mut fpeak) == 0 {
        fpeak = (0.1 / dt as f64) as f32;
    }
    let mut prim = i32::MAX;
    getparint("prim", &mut prim);

    if natten != 0 && natten != 1 && natten != 2 {
        err!("unknown attenuation specified:  atten=%d\n", natten);
    }
    if 0.5 / (dt as f64) < 2.0 * fpeak as f64 {
        fprintf!(stderr(), "WARNING: ALIASING POSSIBLE\n");
        fprintf!(stderr(), "decrease dt or reduce fpeak\n");
    }

    let mut nang = countparval("ang");
    if nang == 0 {
        nang = 1;
    }
    let mut namp = countparval("amp");
    if namp == 0 {
        namp = 1;
    }
    if nang != namp {
        err!("number of angs must equal number of amps");
    }
    let mut ang = vec![0.0f32; nang as usize];
    let mut amp = vec![0.0f32; namp as usize];
    if getparfloat("ang", &mut ang) == 0 {
        ang[0] = 0.0;
    }
    if getparfloat("amp", &mut amp) == 0 {
        amp[0] = 1.0;
    }
    checkpars();

    for iang in 0..nang as usize {
        if iang < nang as usize - 1 && ang[iang] > ang[iang + 1] {
            err!("angs must increase monotonically");
        }
        ang[iang] = (ang[iang] as f64 * (PI / 180.0)) as f32;
    }
    if let Some(fp) = &infofp {
        fprintf!(fp, "THIS FILE CONTAINS BEAM INFORMATION\n\n");
    }

    let re = read_rayends(&stdin());
    if ng < 0 {
        // C's conversion of negative ns to size_t fails ealloc1float in
        // makexzs. Preserve that diagnostic instead of a Vec capacity panic.
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    makexzs(nxz, &mut xg, &mut zg, ng as usize);
    makesyn(
        prim,
        reftrans,
        &re,
        &ang,
        &amp,
        lscale,
        krecord,
        ng as usize,
        &xg,
        &zg,
        fpeak,
        aperture,
        bw,
        nt,
        dt,
        ft,
        natten,
        infofp.as_ref(),
    );
    exit(CWP_Exit());
}

/// Read whole native RayEnd records; a final partial record is discarded.
fn read_rayends(fp: &CFile) -> Vec<RayEnd> {
    let mut re = Vec::with_capacity(200);
    let mut ray = [RayEnd::default()];
    while efread(&mut ray, fp) == 1 {
        re.push(ray[0]);
        if re.len() == re.capacity() {
            re.reserve_exact(200);
        }
    }
    re
}

/// Interpolate (x,z) coordinates uniformly sampled in distance s.
fn makexzs(nu: usize, xs: &mut [f32], zs: &mut [f32], ns: usize) {
    let mut xud = vec![[0.0f32; 4]; nu];
    let mut zud = vec![[0.0f32; 4]; nu];
    let mut u = vec![0.0f32; nu];
    for iu in 0..nu {
        u[iu] = iu as f32;
    }
    csplin(&u, &xs[..nu], &mut xud);
    csplin(&u, &zs[..nu], &mut zud);
    let nuu = 20 * nu;
    let duu = (u[nu - 1] - u[0]) / (nuu - 1) as f32;
    // C reads s[-1] (heap garbage) on the first iteration, so its receiver
    // positions depend on ray input size/order and open files. At >=94
    // records (>8 KB) even the stock demo can differ. Using zero as the
    // predecessor preserves the intended geometry for every input.
    let mut storage = vec![0.0f32; nuu + 1];
    let mut xlast = xs[0];
    let mut zlast = zs[0];
    let mut uu = 0.0f32;
    for iuu in 0..nuu {
        let mut x = [0.0f32];
        let mut z = [0.0f32];
        intcub(0, &u, &xud, &[uu], &mut x);
        intcub(0, &u, &zud, &[uu], &mut z);
        let dx = x[0] - xlast;
        let dz = z[0] - zlast;
        // storage[0] defines C's invalid s[-1] read as zero when iuu == 0.
        storage[iuu + 1] = (storage[iuu] as f64 + ((dx * dx + dz * dz) as f64).sqrt()) as f32;
        xlast = x[0];
        zlast = z[0];
        uu += duu;
    }
    let s = &storage[1..];
    let mut us = vec![0.0f32; ns];
    let ds = (s[nuu - 1] - s[0]) / (if ns > 1 { ns - 1 } else { 1 }) as f32;
    yxtoxy(duu, 0.0, s, ds, 0.0, 0.0, (nu - 1) as f32, &mut us);
    intcub(0, &u, &xud, &us, &mut xs[..ns]);
    intcub(0, &u, &zud, &us, &mut zs[..ns]);
}

/// Make and write synthetic seismograms.
#[allow(clippy::too_many_arguments)]
fn makesyn(
    prim: i32,
    reftrans: i32,
    re: &[RayEnd],
    ang: &[f32],
    amp: &[f32],
    lscale: f32,
    krecord: i32,
    ng: usize,
    xg: &[f32],
    zg: &[f32],
    fpeak: f32,
    aperture: f32,
    bw: f32,
    nt: i32,
    dt: f32,
    ft: f32,
    natten: i32,
    infofp: Option<&CFile>,
) {
    let ntfft = npfaro(nt, nt.wrapping_mul(2)) as usize;
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let wpeak = (2.0 * PI * fpeak as f64) as f32;
    let delay = 0.0f32;
    let wref = 1000.0f32;
    let sqlscale = lscale * lscale;
    let c = 0.99f32;
    let eps = 0.0001f32;

    let mut syn = vec![0.0f32; ntfft];
    let mut cwave = vec![Complex::default(); nw];
    let mut csyn = vec![vec![Complex::default(); nw]; ng];
    let mut lnvec = vec![0.0f32; nw];
    let mut attenvec = vec![0.0f32; nw];
    let mut count = vec![0i32; ng];

    if natten == 2 {
        let temp = (1.0 / PI) as f32;
        lnvec[0] = -1.0;
        let mut w = dw;
        for iw in 1..nw {
            lnvec[iw] = (((w / wref) as f64).ln() * temp as f64) as f32;
            w += dw;
        }
    }
    if let Some(fp) = infofp {
        fprintf!(fp, "***** %d rays in rayendfile.\n***** %s attenuation is simulated.\n***** Geometrical spreading is two-and-one-half dimensional.\n",
            re.len(), if natten == 0 { "No" } else if natten == 1 { "Noncausal" } else { "Causal" });
    }

    for (ire, ray) in re.iter().enumerate() {
        let kend = ray.kend;
        if kend != krecord {
            continue;
        }
        let nref = ray.nref;
        if nref != prim && prim != i32::MAX {
            continue;
        }
        let px = ray.px;
        let v = (1.0 / (ray.se as f64).sqrt()) as f32;
        if abs(px * v) > aperture {
            if let Some(fp) = infofp {
                fprintf!(fp, "Ray %d is out of the aperature\n", ire);
            }
            continue;
        }
        let sigma = (1.0 / (ray.sigma as f64).sqrt()) as f32;
        let x = ray.x;
        let z = ray.z;
        let t = ray.t;
        let pz = ray.pz;
        let q1 = ray.q1;
        let p1 = ray.p1;
        let q2 = ray.q2;
        let p2 = ray.p2;
        let kmah = ray.kmah;
        let atten = ray.atten;
        let ampli = ray.ampli;
        let ampliphase = ray.ampliphase;
        let mut dangle = ray.dangle;
        let dsdx = ray.dsdxe;
        let dsdz = ray.dsdze;
        let v0 = (1.0 / (ray.sb as f64).sqrt()) as f32;
        let a0 = ray.ab;
        let v2 = v * v;
        let v4 = v2 * v2;
        let dvds = (-0.5 * v4 as f64 * (dsdx * px + dsdz * pz) as f64) as f32;
        let dvdn = (-0.5 * v4 as f64 * (dsdx * pz - dsdz * px) as f64) as f32;
        if let Some(fp) = infofp {
            fprintf!(
                fp,
                "\n***Ray emerges at (x=%g,z=%g) on interface %d.\n",
                x,
                z,
                kend
            );
            fprintf!(
                fp,
                "   Emergence angle=%g    time=%g\n",
                ((v * px) as f64).asin() * 180.0 / PI,
                t
            );
            fprintf!(fp, "   v=%g    dsdx=%g    dsdz=%g\n", v, dsdx, dsdz);
            fprintf!(fp, "   dvds=%g    dvdn=%g\n", dvds, dvdn);
            fprintf!(fp, "   ampli=%g    ampliphase=%g\n", ampli, ampliphase);
        }
        if dangle == 0.0 {
            dangle = 1.0;
        }
        let scale =
            ((sigma * dangle) as f64 * (v as f64 / ((v0 * 2.0f32) as f64 * PI)).sqrt()) as f32;
        let mut ampa = [0.0f32];
        intlin(ang, amp, amp[0], amp[amp.len() - 1], &[a0], &mut ampa);

        let (eps1, eps2);
        if bw != 0.0 {
            eps1 = 0.0f32;
            eps2 = (-0.5 * wpeak as f64 * bw as f64 * bw as f64) as f32;
        } else {
            let temp = c * abs(p2 / q2) + eps;
            let mut e1 = -q2 / q1 * (p2 * p1 / (q2 * q1) + temp * temp);
            e1 = e1 / (p1 * p1 / (q1 * q1) + temp * temp);
            eps1 = e1;
            eps2 = -temp / (p1 * p1 + q1 * q1 * temp * temp);
        }
        let ceps = cmplx(eps1, eps2);
        let mut cp = crmul(ceps, p1);
        cp.r += p2;
        let mut cq = crmul(ceps, q1);
        cq.r += q2;
        let cpoq = cdiv(cp, cq);
        let cpoqr = cpoq.r;
        let cpoqi = cpoq.i;
        if cpoqi < 0.0 {
            if let Some(fp) = infofp {
                fprintf!(fp, "cpoqi negative, skipping ray!");
            }
            continue;
        }
        let cqr = if nref % 2 != 0 { cneg(cq) } else { cq };
        let mut phase = (-PI * (kmah.wrapping_add(1) / 2) as f64) as f32;
        let mut camp = cwp_csqrt(cdiv(ceps, cqr));
        if let Some(fp) = infofp {
            let temp = ((2.0f32 * ((eps1 * q1 + q2) * (eps1 * q1 + q2) + eps2 * eps2 * q1 * q1)
                / (-wpeak * eps2)) as f64)
                .sqrt() as f32;
            fprintf!(fp, "   Gaussian beam information:\n");
            fprintf!(fp, "   Re(M)=%g \t Im(M)=%g\n", cpoqr, cpoqi);
            fprintf!(
                fp,
                "   const beam phase=%g degrees\n",
                180.0 / PI * ((camp.i / (camp.r + eps)) as f64).atan()
            );
            fprintf!(fp, "   epsilon1=%g \t epsilon2=%g\n", eps1, eps2);
            fprintf!(fp, "   ray-end beam width=%g\n", temp);
        }
        if reftrans != 0 {
            phase += ampliphase;
            camp = crmul(camp, ampli);
        }
        phase = (phase as f64 + PI / 4.0) as f32;
        camp = crmul(camp, scale * ampa[0]);
        camp = cmul(
            camp,
            cmplx((phase as f64).cos() as f32, (phase as f64).sin() as f32),
        );
        let campr = camp.r as f64;
        let campi = camp.i as f64;
        if natten == 2 {
            for iw in 0..nw {
                attenvec[iw] = atten * lnvec[iw];
            }
        }
        for ig in 0..ng {
            let dx = xg[ig] - x;
            let dz = zg[ig] - z;
            let disq = dx * dx + dz * dz;
            let n = v * (pz * dx - px * dz);
            let nn = n * n;
            if lscale != f32::MAX && disq > sqlscale {
                continue;
            }
            count[ig] = count[ig].wrapping_add(1);
            let dels = ((v * (px * dx + pz * dz)) as f64 * (1.0 - (dvdn / v * n) as f64)) as f32;
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
                let mut w = 0.0f32;
                for iw in 0..nw {
                    let tempd = (w * (ta - attenvec[iw])) as f64;
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
    let mut w = 0.0f32;
    for iw in 0..nw {
        cwave[iw] = cricker(w, wpeak, delay);
        cwave[iw] = crmul(cwave[iw], (w as f64).sqrt() as f32);
        w += dw;
    }
    for ig in 0..ng {
        for iw in 0..nw {
            csyn[ig][iw] = cmul(csyn[ig][iw], cwave[iw]);
        }
        pfacr(-1, &csyn[ig], &mut syn);
        if efwrite(&syn[..nt as usize], &stdout()) != nt as usize {
            err!("error writing output traces!\n");
        }
    }
    if let Some(fp) = infofp {
        fprintf!(
            fp,
            "\nBeam contribution to each receiver *************************\n"
        );
        fprintf!(fp, "Receiver No\t x_coord. \t z_coord. \t No. of beams\n");
        for ig in 0..ng {
            fprintf!(
                fp,
                "%d \t\t %g \t\t %g \t\t %d\n",
                ig,
                xg[ig],
                zg[ig],
                count[ig]
            );
        }
    }
}

/// Fourier transform of the Ricker wavelet (Dave Hale, 02/28/91).
fn cricker(w: f32, wpeak: f32, delay: f32) -> Complex {
    crmul(
        cwp_cexp(cmplx(-((w / wpeak) as f64).powf(2.0) as f32, delay * w)),
        (4.0 * w as f64 * w as f64 * PI.sqrt() / (wpeak as f64).powf(3.0)) as f32,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn receiver_positions_are_uniform_after_large_rayend_input() {
        // Audit M1 (audits/claude-wave3-tri-progs.md): the C s[-1] read
        // corrupts this geometry with 130 demo rays.
        // Also cover both growth steps of the ray buffer. These inputs
        // exceed the safe size for C goldens, so pin the intended result.
        for nre in [130, 201, 401] {
            let fp = etmpfile();
            for ire in 0..nre {
                let ray = RayEnd {
                    sigma: 1.0,
                    x: 4.0 * ire as f32 / (nre - 1) as f32,
                    pz: 1.0,
                    t: 0.2,
                    q1: 1.0,
                    q2: 0.2,
                    p2: 1.0,
                    sb: 1.0,
                    se: 1.0,
                    kend: 1,
                    ampli: 1.0,
                    dangle: 0.1,
                    ..RayEnd::default()
                };
                assert_eq!(efwrite(&[ray], &fp), 1);
            }
            erewind(&fp);
            let re = read_rayends(&fp);
            assert_eq!(re.len(), nre);
            assert!(std::mem::size_of_val(re.as_slice()) > 8192);

            let mut xg = [0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0];
            let mut zg = [0.0; 7];
            makexzs(2, &mut xg, &mut zg, 7);
            for ig in 0..7 {
                let expected = 4.0 * ig as f32 / 6.0;
                assert!(
                    (xg[ig] - expected).abs() < 1e-6,
                    "{nre} rays, receiver {ig}: x={} instead of {expected}",
                    xg[ig]
                );
                assert_eq!(zg[ig], 0.0);
            }
            // Keep the ray allocation live while interpolating, as in main.
            assert_eq!(re[nre - 1].x, 4.0);
            efclose(&fp);
        }
    }
}
