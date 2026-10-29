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

//! SUKDSYN2D - Kirchhoff depth synthesis of 2D seismic data.
//! Port of su/main/synthetics_waveforms_testpatterns/sukdsyn2d.c (Revision 1.12).
//! Credits: CWP: Zhenyue Liu, 07/24/95, Colorado School of Mines.

use oxwave_core::prelude::*;
use oxwave_cwp::{
    complex::Complex,
    cwp::{abs, nint, PI},
};
use oxwave_f_convolution::convolve_cwp;
use oxwave_f_modeling::{makericker, Wavelet};
use oxwave_finite_difference::mkhdiff;
use oxwave_pfa_fft::{npfaro, pfacr, pfarc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUKDSYN2D - Kirchhoff Depth SYNthesis of 2D seismic data from a\t",
    "             migrated seismic section\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "   sukdsyn2d  infile  outfile [parameters] \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters: \t\t\t\t\t\t\t",
    " infile=stdin\t\tinput migrated section\t\t\t\t",
    " outfile=stdout\tfile for output seismic traces  \t\t",
    " ttfile\t\tfile for input traveltime tables\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The following 9 parameters describe traveltime tables:\t\t",
    " fzt= \t\t\tfirst depth sample \t\t\t\t",
    " nzt= \t\t\tnumber of depth samples \t\t\t",
    " dzt=\t\t\tdepth interval \t\t\t\t\t",
    " fxt=\t\t\tfirst lateral sample \t\t\t\t",
    " nxt=\t\t\tnumber of lateral samples \t\t\t",
    " dxt=\t\t\tlateral interval \t\t\t\t",
    " fs=\t\t\tx-coordinate of first source\t\t\t",
    " ns=\t\t\tnumber of sources\t\t\t\t",
    " ds=\t\t\tx-coordinate increment of sources\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The following 6 parameters describe the migration section:\t\t",
    " fz=                   first z-coordinate in migrated section \t\t",
    " dz=     \t\tvertical spacing of migrated section \t\t",
    " nz=           \tnumber of depth points in migrated section\t",
    " fx=                   first x-coordinate of migrated section \t\t",
    " dx=     \t\thorizontal spacing of migrated section \t\t",
    " nx=           \tnumber of lateral points in migrated section  \t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nt=501        \tnumber of time samples\t\t\t\t",
    " dt=0.004      \ttime sampling interval (sec)\t\t\t",
    " ft=0.0        \tfirst time (sec)\t\t\t\t",
    " nxo=1                 number of source-receiver offsets\t\t",
    " dxo=25                offset sampling interval  \t\t\t",
    " fxo=0.0               first offset  \t\t\t\t\t",
    " nxs=101               number of shotpoints  \t\t\t\t",
    " dxs=25                shotpoint sampling interval  \t\t\t",
    " fxs=0.0               first shotpoint \t\t\t\t",
    " fmax=1/(4*dt)         maximum frequency in migration section (Hz)\t",
    " aperx=nxt*dxt/2  \tmodeling lateral aperature \t\t\t",
    " angmax=60\t\tmodeling angle aperature from vertical\t\t",
    " v0=1500(m/s)\t\treference velocity value at surface\t\t",
    " dvz=0.0  \t\treference velocity vertical gradient\t\t",
    " ls=1\t                flag for line source\t\t\t\t",
    " jpfile=stderr\t\tjob print file name \t\t\t\t",
    " mtr=100  \t\tprint verbal information at every mtr traces\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This program takes a migrated seismic section and a set of travel time",
    " tables generated using rayt2d for a specific background velocity model",
    " and generates synthetic seismic data in the form of common shot gathers.",
    " (Common offset gathers may be generated by using nxo=1.) (Demigration.)",
    "\t\t\t\t\t\t\t\t\t",
    " This program is a tool which may be used for the migration residual\t",
    " statics estimation technique of Tjan, Audebert, and Larner 1994.\t",
    "\t\t\t\t\t\t\t\t\t",
    "1. The traveltime tables are generated by program rayt2d (or other ones)",
    "   on relatively coarse grids, with dimension ns*nxt*nzt. In the\t",
    "   modeling process, traveltimes are interpolated into shot/geophone \t",
    "   positions and migration section grids. \t\t\t\t",
    "2. The input migration section must be an array of binary floats (no SU",
    "   headers).\t\t\t\t\t\t\t\t",
    "3. The synthesized traces are output in common-shot gathers in SU format.",
    "4. The memory requirement for this program is about\t\t\t",
    "    \t(ns*nxt*nzt+nx*nz+4*nr*nzt+3*nxt*nzt)*4 bytes \t\t\t",
    "    where nr = 1+min(nxt*dxt,0.5*offmax+aperx)/dxo. \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut infile = String::from("stdin");
    let mut outfile = String::from("stdout");
    let mut ttfile = String::new();
    let mut jpfile = String::new();
    let infp = if getparstring("infile", &mut infile) == 0 {
        stdin()
    } else {
        fopen(&infile, "r").unwrap_or_else(|| err!("cannot open infile=%s\n", infile))
    };
    let outfp = if getparstring("outfile", &mut outfile) == 0 {
        Some(stdout())
    } else {
        // C leaves a failed fopen unchecked until the stream is used.
        fopen(&outfile, "w")
    };
    efseeko(&infp, 0, SEEK_CUR);
    warn!("Got A");
    let outfp = outfp.unwrap_or_else(|| unsafe { fault() });
    efseeko(&outfp, 0, SEEK_CUR);
    if getparstring("ttfile", &mut ttfile) == 0 {
        err!("must specify ttfile!\n");
    }
    let ttfp = fopen(&ttfile, "r").unwrap_or_else(|| err!("cannot open ttfile=%s\n", ttfile));
    let jpfp = if getparstring("jpfile", &mut jpfile) == 0 {
        Some(stderr())
    } else {
        fopen(&jpfile, "w")
    };

    let mut nt = 501;
    getparint("nt", &mut nt);
    let mut dt = 0.004f32;
    getparfloat("dt", &mut dt);
    let mut ft = 0.0f32;
    getparfloat("ft", &mut ft);
    let mut fmax = (1.0 / (4.0 * dt) as f64) as f32;
    getparfloat("fmax", &mut fmax);
    let fpeak = (0.2 / dt as f64) as f32;
    let mut nxs = 101;
    getparint("nxs", &mut nxs);
    let mut dxs = 15.0f32;
    getparfloat("dxs", &mut dxs);
    let mut fxs = 0.0f32;
    getparfloat("fxs", &mut fxs);
    let mut nxo = 1;
    getparint("nxo", &mut nxo);
    let mut dxo = 50.0f32;
    getparfloat("dxo", &mut dxo);
    let mut fxo = 0.0f32;
    getparfloat("fxo", &mut fxo);

    // Retain the original required-parameter order and diagnostic strings.
    let mut nxt = 0;
    if getparint("nxt", &mut nxt) == 0 {
        err!("must specify nxt!\n");
    }
    let mut fxt = 0.0f32;
    if getparfloat("fxt", &mut fxt) == 0 {
        err!("must specify fxt!\n");
    }
    let mut dxt = 0.0f32;
    if getparfloat("dxt", &mut dxt) == 0 {
        err!("must specify dxt!\n");
    }
    let mut nzt = 0;
    if getparint("nzt", &mut nzt) == 0 {
        err!("must specify nzt!\n");
    }
    let mut fzt = 0.0f32;
    if getparfloat("fzt", &mut fzt) == 0 {
        err!("must specify fzt!\n");
    }
    let mut dzt = 0.0f32;
    if getparfloat("dzt", &mut dzt) == 0 {
        err!("must specify dzt!\n");
    }
    let mut ns = 0;
    if getparint("ns", &mut ns) == 0 {
        err!("must specify ns!\n");
    }
    let mut fs = 0.0f32;
    if getparfloat("fs", &mut fs) == 0 {
        err!("must specify fs!\n");
    }
    let mut ds = 0.0f32;
    if getparfloat("ds", &mut ds) == 0 {
        err!("must specify ds!\n");
    }
    let ext = fxt + (nxt - 1) as f32 * dxt;
    let ezt = fzt + (nzt - 1) as f32 * dzt;
    let es = fs + (ns - 1) as f32 * ds;

    for ixs in 0..nxs {
        let xs = fxs + ixs as f32 * dxs;
        for ixo in 0..nxo {
            let xg = xs + fxo + ixo as f32 * dxo;
            if fs > xs || es < xs || fs > xg || es < xg {
                err!("shot or receiver lie outside of specified (x,z) grid\n");
            }
        }
    }

    let mut nx = 0;
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx!\n");
    }
    let mut fx = 0.0f32;
    if getparfloat("fx", &mut fx) == 0 {
        err!("must specify fx!\n");
    }
    let mut dx = 0.0f32;
    if getparfloat("dx", &mut dx) == 0 {
        err!("must specify dx!\n");
    }
    let mut nz = 0;
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!\n");
    }
    let mut fz = 0.0f32;
    if getparfloat("fz", &mut fz) == 0 {
        err!("must specify fz!\n");
    }
    let mut dz = 0.0f32;
    if getparfloat("dz", &mut dz) == 0 {
        err!("must specify dz!\n");
    }
    let ex = fx + (nx - 1) as f32 * dx;
    let ez = fz + (nz - 1) as f32 * dz;
    if fxt > fx || ext < ex || fzt > fz || ezt < ez {
        err!(" migration section is out of traveltime table!\n");
    }

    let mut v0 = 1500.0f32;
    getparfloat("v0", &mut v0);
    let mut dvz = 0.0f32;
    getparfloat("dvz", &mut dvz);
    let mut angmax = 60.0f32;
    getparfloat("angmax", &mut angmax);
    if (angmax as f64) < 0.00001 {
        err!("angmax must be positive!\n");
    }
    let mut mzmax = (dx as f64 * (angmax as f64 * PI / 180.0).sin() / dz as f64) as i32;
    if mzmax < 1 {
        mzmax = 1;
    }
    let mut aperx = (0.5 * nxt as f64 * dxt as f64) as f32;
    getparfloat("aperx", &mut aperx);
    let mut ls = 1;
    getparint("ls", &mut ls);
    let mut mtr = 100;
    getparint("mtr", &mut mtr);

    // A null job stream faults at the first fprintf, after parameter checks.
    let jpfp = jpfp.unwrap_or_else(|| unsafe { fault() });
    fprintf!(&jpfp, "\n");
    fprintf!(&jpfp, " Modeling parameters\n");
    fprintf!(&jpfp, " ================\n");
    fprintf!(&jpfp, " infile=%s \n", infile);
    fprintf!(&jpfp, " outfile=%s \n", outfile);
    fprintf!(&jpfp, " ttfile=%s \n", ttfile);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nzt=%d fzt=%g dzt=%g\n", nzt, fzt, dzt);
    fprintf!(&jpfp, " nxt=%d fxt=%g dxt=%g\n", nxt, fxt, dxt);
    fprintf!(&jpfp, " ns=%d fs=%g ds=%g\n", ns, fs, ds);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nz=%d fz=%g dz=%g\n", nz, fz, dz);
    fprintf!(&jpfp, " nx=%d fx=%g dx=%g\n", nx, fx, dx);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nxs=%d fxs=%g dxs=%g\n", nxs, fxs, dxs);
    fprintf!(&jpfp, " nxo=%d fxo=%g dxo=%g\n", nxo, fxo, dxo);
    fprintf!(&jpfp, " \n");

    let lastoff = abs(fxo + (nxo - 1) as f32 * dxo);
    let offmax = if abs(fxo) > lastoff {
        abs(fxo)
    } else {
        lastoff
    };
    let mut rmax = if es - fxt > ext - fs {
        es - fxt
    } else {
        ext - fs
    };
    let limit = 0.5 * offmax as f64 + aperx as f64;
    rmax = if (rmax as f64) < limit {
        rmax as f64
    } else {
        limit
    } as f32;
    let nr = 2 + (rmax / dx) as i32;
    let mut tb = table(nr, nzt);
    let mut pb = table(nr, nzt);
    let mut sigb = table(nr, nzt);
    let mut cosb = table(nr, nzt);
    timeb(
        nr, nzt, dx, dzt, fzt, dvz, v0, &mut tb, &mut pb, &mut sigb, &mut cosb,
    );
    fprintf!(&jpfp, " nt=%d ft=%g dt=%g fpeak=%g \n", nt, ft, dt, fpeak);
    fprintf!(&jpfp, " v0=%g dvz=%g \n", v0, dvz);
    fprintf!(
        &jpfp,
        " aperx=%g angmax=%g offmax=%g\n",
        aperx,
        angmax,
        offmax
    );
    fprintf!(&jpfp, " mtr=%d ls=%d\n", mtr, ls);
    fprintf!(&jpfp, " ================\n");
    fflush(&jpfp);

    let mut mig = table(nx, nz);
    let mut migi = table(nx, nz + 2 * mzmax);
    let mut ttab = vec![table(nxt, nzt); ns as usize];
    let mut tt = table(nxt, nzt);
    let mut tsum = table(nxt, nzt);

    let w = makericker(fpeak, dt);
    let mut tr = Trace::new();
    tr.trid = 1;
    tr.counit = 1;
    tr.ns = nt as u16;
    tr.dt = (1.0e6 * dt as f64) as u32 as u16;
    tr.delrt = (1.0e3 * ft as f64) as i32 as i16;

    fprintf!(&jpfp, " read traveltime tables \n");
    fflush(&jpfp);
    for is in 0..ns {
        let nseek = nxt as i64 * nzt as i64 * is as i64;
        efseeko(&ttfp, nseek * 4, SEEK_SET);
        fread(&mut ttab[is as usize], &ttfp);
        let s = fs + is as f32 * ds;
        resit(nxt, fxt, dxt, nzt, nr, dx, &tb, &mut ttab[is as usize], s);
    }
    fprintf!(&jpfp, " read migration section \n");
    fflush(&jpfp);
    fread(&mut mig, &infp);
    integ(&mig, nz, dz, nx, mzmax, &mut migi);
    fprintf!(&jpfp, " start synthesis ... \n");
    fprintf!(&jpfp, " \n");
    fflush(&jpfp);

    let mut jtr = 0;
    let mut xs = fxs;
    let mut tracl = 0;
    let mut hd = [0.0f32; 41];
    let mut madehd = false;
    for ixs in 0..nxs {
        let mut xo = fxo;
        for ixo in 0..nxo {
            let xg = xs + xo;
            tracl += 1;
            tr.tracl = tracl;
            tr.tracr = tracl;
            tr.fldr = 1 + ixs;
            tr.tracf = 1 + ixo;
            tr.offset = nint(xo as f64);
            tr.sx = nint(xs as f64);
            tr.gx = nint(xg as f64);

            let a_s = (xs - fs) / ds;
            let mut is = a_s as i32;
            if is == ns - 1 {
                is = ns - 2;
            }
            let mut res = a_s - is as f32;
            if (res as f64) <= 0.01 {
                res = 0.0;
            }
            if (res as f64) >= 0.99 {
                res = 1.0;
            }
            sum2(
                nxt,
                nzt,
                1.0 - res,
                res,
                &ttab[is as usize],
                &ttab[(is + 1) as usize],
                &mut tsum,
            );

            let a_s = (xg - fs) / ds;
            let mut is = a_s as i32;
            if is == ns - 1 {
                is = ns - 2;
            }
            let mut res = a_s - is as f32;
            if (res as f64) <= 0.01 {
                res = 0.0;
            }
            if (res as f64) >= 0.99 {
                res = 1.0;
            }
            sum2(
                nxt,
                nzt,
                1.0 - res,
                res,
                &ttab[is as usize],
                &ttab[(is + 1) as usize],
                &mut tt,
            );
            // sum2(nxt,nzt,1,1,tt,tsum,tsum): output aliases its second input.
            for ix in 0..nxt {
                for iz in 0..nzt {
                    let k = (ix * nzt + iz) as usize;
                    tsum[k] = 1.0 * tt[k] + 1.0 * tsum[k];
                }
            }
            fflush(&jpfp);
            maketrace(
                &mut tr.data,
                nt,
                ft,
                dt,
                xs,
                xg,
                &mig,
                &migi,
                aperx,
                nx,
                fx,
                dx,
                nz,
                fz,
                dz,
                mzmax,
                ls,
                angmax,
                v0,
                fmax,
                &w,
                &tb,
                &pb,
                &sigb,
                &cosb,
                nr,
                &tsum,
                nxt,
                fxt,
                dxt,
                nzt,
                fzt,
                dzt,
                &mut hd,
                &mut madehd,
            );
            fputtr(&outfp, &mut tr);
            jtr += 1;
            // arm64 SDIV with a zero divisor yields zero, so MSUB leaves jtr-1.
            let rem = if mtr == 0 { jtr - 1 } else { (jtr - 1) % mtr };
            if rem == 0 {
                fprintf!(&jpfp, " generated trace %d\n", jtr);
                fflush(&jpfp);
            }
            xo += dxo;
        }
        xs += dxs;
    }
    fprintf!(&jpfp, " generated %d traces in total\n", jtr);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " output done\n");
    fflush(&jpfp);
    efclose(&jpfp);
    efclose(&outfp);
    exit(CWP_Exit());
}

// Preserve the C null-FILE fault for unchecked fopen calls.
unsafe fn fault() -> ! {
    unsafe extern "C" {
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        raise(11);
    }
    exit(139)
}

fn table(n2: i32, n1: i32) -> Vec<f32> {
    vec![0.0; (n2 as usize).wrapping_mul(n1 as usize)]
}

fn integ(mig: &[f32], nz: i32, dz: f32, nx: i32, m: i32, migi: &mut [f32]) {
    let nfft = npfaro(nz + m, 2 * (nz + m));
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 {
        err!("Padded nt=%d -- too big", nfft);
    }
    let nw = nfft / 2 + 1;
    let dw = (2.0 * PI / (nfft as f32 * dz) as f64) as f32;
    let mut amp = vec![0.0f32; nw as usize];
    for iw in 1..nw {
        amp[iw as usize] =
            (0.5 / (nfft as f64 * (1.0 - ((iw as f32 * dw * dz) as f64).cos()))) as f32;
    }
    amp[0] = amp[1];
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nw as usize];
    for ix in 0..nx {
        rt[..nz as usize].copy_from_slice(&mig[(ix * nz) as usize..((ix + 1) * nz) as usize]);
        rt[nz as usize..].fill(0.0);
        pfarc(1, &rt, &mut ct);
        for iw in 0..nw as usize {
            ct[iw].i *= amp[iw];
            ct[iw].r *= amp[iw];
        }
        pfacr(-1, &ct, &mut rt);
        for iz in 0..m {
            migi[(ix * (nz + 2 * m) + iz) as usize] = rt[(nfft - m + iz) as usize];
        }
        for iz in 0..nz + m {
            migi[(ix * (nz + 2 * m) + iz + m) as usize] = rt[iz as usize];
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn resit(nx: i32, fx: f32, dx: f32, nz: i32, nr: i32, dr: f32, tb: &[f32], t: &mut [f32], x0: f32) {
    for ix in 0..nx {
        let xi = fx + ix as f32 * dx - x0;
        let ar = xi.abs() / dr;
        let mut jr = ar as i32;
        let sr = ar - jr as f32;
        let sr0 = (1.0 - sr as f64) as f32;
        if jr > nr - 2 {
            jr = nr - 2;
        }
        for iz in 0..nz {
            t[(ix * nz + iz) as usize] -=
                sr0 * tb[(jr * nz + iz) as usize] + sr * tb[((jr + 1) * nz + iz) as usize];
        }
    }
}

fn sum2(nx: i32, nz: i32, a1: f32, a2: f32, t1: &[f32], t2: &[f32], t: &mut [f32]) {
    for ix in 0..nx {
        for iz in 0..nz {
            let k = (ix * nz + iz) as usize;
            t[k] = a1 * t1[k] + a2 * t2[k];
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn timeb(
    nr: i32,
    nz: i32,
    dr: f32,
    dz: f32,
    fz: f32,
    a: f32,
    v0: f32,
    t: &mut [f32],
    p: &mut [f32],
    sig: &mut [f32],
    cs: &mut [f32],
) {
    if a == 0.0 {
        let mut r = 0.0f32;
        for ir in 0..nr {
            let mut z = fz;
            for iz in 0..nz {
                let k = (ir * nz + iz) as usize;
                let mut rou = ((r * r + z * z) as f64).sqrt() as f32;
                t[k] = rou / v0;
                if rou < dz {
                    rou = dz;
                }
                t[k] = rou / v0;
                p[k] = r / (rou * v0);
                sig[k] = v0 * rou;
                cs[k] = z / rou;
                z += dz;
            }
            r += dr;
        }
    } else {
        let oa = (1.0 / a as f64) as f32;
        let zc = v0 * oa;
        let mut r = 0.0f32;
        for ir in 0..nr {
            let mut z = fz + zc;
            for iz in 0..nz {
                let k = (ir * nz + iz) as usize;
                let _rou = ((r * r + z * z) as f64).sqrt() as f32;
                let v = v0 + a * (z - zc);
                if ir == 0 {
                    t[k] = (((v / v0) as f64).ln() * oa as f64) as f32;
                    p[k] = 0.0;
                    sig[k] = (0.5 * (z - zc) as f64 * (v0 + v) as f64) as f32;
                    cs[k] = 1.0;
                } else {
                    let rc = ((r * r + z * z - zc * zc) as f64 / (2.0 * r as f64)) as f32;
                    let rou = ((zc * zc + rc * rc) as f64).sqrt() as f32;
                    t[k] = ((((v * (rou + rc)) / (v0 * (rou + rc - r))) as f64).ln() * oa as f64)
                        as f32;
                    p[k] = (((rou * rou - rc * rc) as f64).sqrt() / (rou * v0) as f64) as f32;
                    cs[k] = (rc - r) / rou;
                    sig[k] = a * rou * r;
                }
                z += dz;
            }
            r += dr;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn maketrace(
    trace: &mut [f32],
    nt: i32,
    ft: f32,
    dt: f32,
    xs: f32,
    xg: f32,
    mig: &[f32],
    migi: &[f32],
    aperx: f32,
    nx: i32,
    fx: f32,
    dx: f32,
    nz: i32,
    fz: f32,
    dz: f32,
    mzmax: i32,
    ls: i32,
    angmax: f32,
    v0: f32,
    fmax: f32,
    w: &Wavelet,
    tb: &[f32],
    pb: &[f32],
    sigb: &[f32],
    cosb: &[f32],
    nr: i32,
    tsum: &[f32],
    nxt: i32,
    fxt: f32,
    dxt: f32,
    nzt: i32,
    fzt: f32,
    dzt: f32,
    hd: &mut [f32; 41],
    madehd: &mut bool,
) {
    let odt = (1.0 / dt as f64) as f32;
    if !*madehd {
        mkhdiff(dt, 20, hd);
        *madehd = true;
    }
    for jt in 0..nt {
        trace[jt as usize] = 0.0;
    }
    let mut zmt = table(nxt, nzt);
    let mut ampt = table(nxt, nzt);
    let mut ampti = table(nxt, nzt);
    let mut amp = vec![0.0f32; nzt as usize];
    let mut ampi = vec![0.0f32; nzt as usize];
    let mut zm = vec![0.0f32; nzt as usize];
    let mut tzt = vec![0.0f32; nzt as usize];
    let mut zpt = vec![0.0f32; nxt as usize];
    let mut work1 = vec![0.0f32; nt as usize];
    let z0 = (fz - fzt) / dzt;
    let pmin = (1.0 / (2.0 * dx as f64 * fmax as f64)) as f32;
    let rdz = dz / dzt;
    let xm = (0.5 * (xs + xg) as f64) as f32;
    let rxz = if angmax == 90.0 {
        0.0
    } else {
        (1.0 / (angmax as f64 * PI / 180.0).tan()) as f32
    };
    let mut nxtf = ((xm - aperx - fxt) / dxt) as i32;
    if nxtf < 0 {
        nxtf = 0;
    }
    let mut nxte = (((xm + aperx - fxt) / dxt) as f64 + 1.0) as i32;
    if nxte >= nxt {
        nxte = nxt - 1;
    }

    // Compute amplitudes.
    for ix in nxtf..=nxte {
        let x = fxt + ix as f32 * dxt;
        let dis = if xm >= x { xm - x } else { x - xm };
        let mut izt0 = (((dis - dxt) * rxz - fzt) / dzt - 1.0) as i32;
        if izt0 < 0 {
            izt0 = 0;
        }
        if izt0 >= nzt {
            izt0 = nzt - 1;
        }
        let ar = if xs >= x {
            (xs - x) / dx
        } else {
            (x - xs) / dx
        };
        let mut jrs = ar as i32;
        if jrs > nr - 2 {
            jrs = nr - 2;
        }
        let srs = ar - jrs as f32;
        let srs0 = (1.0 - srs as f64) as f32;
        let ar = if xg >= x {
            (xg - x) / dx
        } else {
            (x - xg) / dx
        };
        let mut jrg = ar as i32;
        if jrg > nr - 2 {
            jrg = nr - 2;
        }
        let srg = ar - jrg as f32;
        let srg0 = (1.0 - srg as f64) as f32;
        let sigp = if (xs - x) * (xg - x) > 0.0 {
            1.0f32
        } else {
            -1.0f32
        };
        zpt[ix as usize] = fzt + (nzt - 1) as f32 * dzt;

        for iz in izt0..nzt {
            let rs = (jrs * nzt + iz) as usize;
            let rg = (jrg * nzt + iz) as usize;
            let k = (ix * nzt + iz) as usize;
            let sigs = srs0 * sigb[rs] + srs * sigb[rs + nzt as usize];
            let sigg = srg0 * sigb[rg] + srg * sigb[rg + nzt as usize];
            let coss = srs0 * cosb[rs] + srs * cosb[rs + nzt as usize];
            let cosg = srg0 * cosb[rg] + srg * cosb[rg + nzt as usize];
            let mut ampd = v0 * dx * (coss + cosg) * dz;
            if ampd < 0.0 {
                ampd = -ampd;
            }
            if ls != 0 {
                ampt[k] = (ampd as f64 / ((v0 * sigs * sigg) as f64).sqrt()) as f32;
            } else {
                ampt[k] = (ampd as f64 / ((sigs * sigg * (sigs + sigg)) as f64).sqrt()) as f32;
            }
            let mut pd = srs0 * pb[rs]
                + srs * pb[rs + nzt as usize]
                + sigp * (srg0 * pb[rg] + srg * pb[rg + nzt as usize]);
            if pd < 0.0 {
                pd = -pd;
            }
            if pd < 0.0 {
                pd = -pd;
            }
            let mut temp = (0.5 * pd as f64 * v0 as f64 * dx as f64 / dz as f64) as f32;
            if temp < 1.0 {
                temp = 1.0;
            }
            if temp > mzmax as f32 {
                temp = mzmax as f32;
            }
            ampti[k] = ampt[k] / (temp * temp);
            zmt[k] = temp;
            if pd < pmin && zpt[ix as usize] as f64 > fzt as f64 + (nzt as f64 - 1.1) * dzt as f64 {
                zpt[ix as usize] = fzt + iz as f32 * dzt;
            }
        }
    }
    let mut nxf = (((xm - aperx - fx) / dx) as f64 + 0.5) as i32;
    if nxf < 0 {
        nxf = 0;
    }
    let mut nxe = (((xm + aperx - fx) / dx) as f64 + 0.5) as i32;
    if nxe >= nx {
        nxe = nx - 1;
    }

    // Interpolate amplitudes.
    for ix in nxf..=nxe {
        let x = fx + ix as f32 * dx;
        let dis = if xm >= x { xm - x } else { x - xm };
        let mut izt0 = ((dis * rxz - fzt) / dzt - 1.0) as i32;
        if izt0 < 0 {
            izt0 = 0;
        }
        if izt0 >= nzt {
            izt0 = nzt - 1;
        }
        let mut iz0 = ((dis * rxz - fz) / dz) as i32;
        if iz0 < 0 {
            iz0 = 0;
        }
        if iz0 >= nz {
            iz0 = nz - 1;
        }
        let mut ax = (x - fxt) / dxt;
        let mut jx = ax as i32;
        ax -= jx as f32;
        if ax as f64 <= 0.01 {
            ax = 0.0;
        }
        if ax as f64 >= 0.99 {
            ax = 1.0;
        }
        let ax0 = (1.0 - ax as f64) as f32;
        if jx > nxte - 1 {
            jx = nxte - 1;
        }
        if jx < nxtf {
            jx = nxtf;
        }
        let ar = if xs >= x {
            (xs - x) / dx
        } else {
            (x - xs) / dx
        };
        let mut jrs = ar as i32;
        if jrs >= nr - 1 {
            jrs = nr - 2;
        }
        let srs = ar - jrs as f32;
        let srs0 = (1.0 - srs as f64) as f32;
        let ar = if xg >= x {
            (xg - x) / dx
        } else {
            (x - xg) / dx
        };
        let mut jrg = ar as i32;
        if jrg >= nr - 1 {
            jrg = nr - 2;
        }
        let srg = ar - jrg as f32;
        let srg0 = (1.0 - srg as f64) as f32;
        for iz in izt0..nzt {
            let k = (jx * nzt + iz) as usize;
            let rs = (jrs * nzt + iz) as usize;
            let rg = (jrg * nzt + iz) as usize;
            tzt[iz as usize] = ax0 * tsum[k]
                + ax * tsum[k + nzt as usize]
                + srs0 * tb[rs]
                + srs * tb[rs + nzt as usize]
                + srg0 * tb[rg]
                + srg * tb[rg + nzt as usize];
            amp[iz as usize] = ax0 * ampt[k] + ax * ampt[k + nzt as usize];
            ampi[iz as usize] = ax0 * ampti[k] + ax * ampti[k + nzt as usize];
            zm[iz as usize] = ax0 * zmt[k] + ax * zmt[k + nzt as usize];
        }
        let mut nzp = (((ax0 * zpt[jx as usize] + ax * zpt[(jx + 1) as usize] - fz) / dz) as f64
            + 1.5) as i32;
        if nzp < iz0 {
            nzp = iz0;
        }
        if nzp > nz {
            nzp = nz;
        }

        // Interpolate along depth if operator aliasing.
        for iz in iz0..nzp {
            let az = z0 + iz as f32 * rdz;
            let mut jz = az as i32;
            if jz >= nzt - 1 {
                jz = nzt - 2;
            }
            let sz = az - jz as f32;
            let sz0 = (1.0 - sz as f64) as f32;
            let td = sz0 * tzt[jz as usize] + sz * tzt[(jz + 1) as usize];
            let at = (td - ft) * odt;
            let jt = at as i32;
            if jt > 0 && jt < nt - 1 {
                let ampd = sz0 * ampi[jz as usize] + sz * ampi[(jz + 1) as usize];
                let mz = (0.5
                    + (sz0 * zm[jz as usize]) as f64
                    + (sz * zm[(jz + 1) as usize]) as f64) as i32;
                let res = at - jt as f32;
                let iz1 = iz + mzmax;
                let base = ix * (nz + 2 * mzmax);
                let temp = ((-migi[(base + iz1 - mz) as usize] as f64
                    + 2.0 * migi[(base + iz1) as usize] as f64
                    - migi[(base + iz1 + mz) as usize] as f64)
                    * ampd as f64) as f32;
                trace[jt as usize] =
                    (trace[jt as usize] as f64 + (1.0 - res as f64) * temp as f64) as f32;
                trace[(jt + 1) as usize] += res * temp;
            }
        }
        // Interpolate along depth if not operator aliasing.
        for iz in nzp..nz {
            let az = z0 + iz as f32 * rdz;
            let mut jz = az as i32;
            if jz >= nzt - 1 {
                jz = nzt - 2;
            }
            let sz = az - jz as f32;
            let sz0 = (1.0 - sz as f64) as f32;
            let td = sz0 * tzt[jz as usize] + sz * tzt[(jz + 1) as usize];
            let at = (td - ft) * odt;
            let jt = at as i32;
            if jt > 0 && jt < nt - 1 {
                let ampd = sz0 * amp[jz as usize] + sz * amp[(jz + 1) as usize];
                let res = at - jt as f32;
                let temp = mig[(ix * nz + iz) as usize] * ampd;
                trace[jt as usize] =
                    (trace[jt as usize] as f64 + (1.0 - res as f64) * temp as f64) as f32;
                trace[(jt + 1) as usize] += res * temp;
            }
        }
    }

    // Apply half-derivative filter to trace.
    convolve_cwp(-20, hd, 0, &trace[..nt as usize], 0, &mut work1);

    // Convolve wavelet with trace.
    convolve_cwp(
        w.iw,
        &w.wv[..w.lw as usize],
        0,
        &work1,
        0,
        &mut trace[..nt as usize],
    );
}
