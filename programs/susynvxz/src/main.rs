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

//! SUSYNVXZ - common-offset synthetic seismograms in V(X,Z) media.
//! Port of su/main/synthetics_waveforms_testpatterns/susynvxz.c (Revision 1.22).
//! Credits: CWP: Zhenyue Liu, 07/20/92.
//! Many subroutines borrowed from Dave Hale's program: SUSYNLV.
//! Original source copyright (c) Colorado School of Mines, 2011.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, max, min, nint};
use oxwave_f_convolution::convolve_cwp;
use oxwave_f_modeling::{
    addsinc, breakReflectors, decodeReflectors, makeref, makericker, Reflector, Wavelet,
};
use oxwave_f_upweik::eiktam;
use oxwave_finite_difference::mkhdiff;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUSYNVXZ - SYNthetic seismograms of common offset V(X,Z) media via\t",
    " \t\tKirchhoff-style modeling\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " susynvxz >outfile [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " <vfile\t\tfile containing velocities v[nx][nz]\t\t",
    " nx=\t\t\tnumber of x samples (2nd dimension)\t\t",
    " nz=\t\t\tnumber of z samples (1st dimension)\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nxb=nx\t\tband centered at midpoint\t\t\t",
    " nxd=1\t\t\tskipped number of midponits\t\t\t",
    " dx=100\t\tx sampling interval (m)\t\t\t\t",
    " fx=0.0\t\tfirst x sample\t\t\t\t\t",
    " dz=100\t\tz sampling interval (m)\t\t\t\t",
    " nt=101\t\tnumber of time samples\t\t\t\t",
    " dt=0.04\t\ttime sampling interval (sec)\t\t\t",
    " ft=0.0\t\tfirst time (sec)\t\t\t\t",
    " nxo=1\t\t \tnumber of offsets\t\t\t\t",
    " dxo=50\t\toffset sampling interval (m)\t\t\t",
    " fxo=0.0\t\tfirst offset (m)\t\t\t\t",
    " nxm=101\t\tnumber of midpoints\t\t\t\t",
    " dxm=50\t\tmidpoint sampling interval (m)\t\t\t",
    " fxm=0.0\t\tfirst midpoint (m)\t\t\t\t",
    " fpeak=0.2/dt\t\tpeak frequency of symmetric Ricker wavelet (Hz)\t",
    " ref=\"1:1,2;4,2\"\treflector(s):  \"amplitude:x1,z1;x2,z2;x3,z3;...\"",
    " smooth=0\t\t=1 for smooth (piecewise cubic spline) reflectors",
    " ls=0\t\t\t=1 for line source; default is point source\t",
    " tmin=10.0*dt\t\tminimum time of interest (sec)\t\t\t",
    " ndpfz=5\t\tnumber of diffractors per Fresnel zone\t\t",
    " verbose=0\t\t=1 to print some useful information\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This algorithm is based on formula (58) in Geo. Pros. 34, 686-703,\t",
    " by N. Bleistein.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Offsets are signed - may be positive or negative.\t\t\t",
    " Traveltime and amplitude are calculated by finite differences which\t",
    " is done only in part of midpoints; in the skiped midpoint, interpolation",
    " is used to calculate traveltime and amplitude.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " More than one ref (reflector) may be specified.\t\t\t",
    " Note that reflectors are encoded as quoted strings, with an optional\t",
    " reflector amplitude: preceding the x,z coordinates of each reflector.\t",
    " Default amplitude is 1.0 if amplitude: part of the string is omitted.\t",
    "\t\t\t\t\t\t\t\t\t",
];

const LHD: i32 = 20;
const NHD: usize = 1 + 2 * LHD as usize;

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    /* get required parameters */
    let mut nx = 0;
    let mut nz = 0;
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx!\n");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!\n");
    }

    /* get optional parameters */
    let mut nxb: i32 = nx;
    if getparint("nxb", &mut nxb) == 0 {
        nxb = nx;
    }
    let mut nxd: i32 = 1;
    if getparint("nxd", &mut nxd) == 0 {
        nxd = 1;
    }
    let mut dx: f32 = 100.0;
    if getparfloat("dx", &mut dx) == 0 {
        dx = 100.0;
    }
    let mut fx: f32 = 0.0;
    if getparfloat("fx", &mut fx) == 0 {
        fx = 0.0;
    }
    let mut dz: f32 = 100.0;
    if getparfloat("dz", &mut dz) == 0 {
        dz = 100.0;
    }
    let mut nt: i32 = 101;
    if getparint("nt", &mut nt) == 0 {
        nt = 101;
    }
    CHECK_NT("nt", nt);
    let mut dt: f32 = 0.04;
    if getparfloat("dt", &mut dt) == 0 {
        dt = 0.04;
    }
    let mut ft: f32 = 0.0;
    if getparfloat("ft", &mut ft) == 0 {
        ft = 0.0;
    }
    let mut nxo: i32 = 1;
    if getparint("nxo", &mut nxo) == 0 {
        nxo = 1;
    }
    let mut dxo: f32 = 50.0;
    if getparfloat("dxo", &mut dxo) == 0 {
        dxo = 50.0;
    }
    let mut fxo: f32 = 0.0;
    if getparfloat("fxo", &mut fxo) == 0 {
        fxo = 0.0;
    }
    let mut nxm: i32 = 101;
    if getparint("nxm", &mut nxm) == 0 {
        nxm = 101;
    }
    let mut dxm: f32 = 50.0;
    if getparfloat("dxm", &mut dxm) == 0 {
        dxm = 50.0;
    }
    let mut fxm: f32 = 0.0;
    if getparfloat("fxm", &mut fxm) == 0 {
        fxm = 0.0;
    }
    let mut fpeak: f32 = (0.2 / dt as f64) as f32;
    if getparfloat("fpeak", &mut fpeak) == 0 {
        fpeak = (0.2 / dt as f64) as f32;
    }
    let mut ls: i32 = 0;
    if getparint("ls", &mut ls) == 0 {
        ls = 0;
    }
    let mut tmin: f32 = (10.0 * dt as f64) as f32;
    if getparfloat("tmin", &mut tmin) == 0 {
        tmin = (10.0 * dt as f64) as f32;
    }
    let mut ndpfz: i32 = 5;
    if getparint("ndpfz", &mut ndpfz) == 0 {
        ndpfz = 5;
    }
    let mut smooth: i32 = 0;
    if getparint("smooth", &mut smooth) == 0 {
        smooth = 0;
    }
    let mut verbose: i32 = 0;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* check the ranges of shots and receivers */
    let ex = fx + (nx - 1) as f32 * dx;
    let ez = (nz - 1) as f32 * dz;
    for ixm in 0..nxm {
        for ixo in 0..nxo {
            let xs =
                ((fxm + ixm as f32 * dxm) as f64 - 0.5 * (fxo + ixo as f32 * dxo) as f64) as f32;
            let xg = xs + fxo + ixo as f32 * dxo;
            if fx > xs || ex < xs || fx > xg || ex < xg {
                err!("shot or receiver lie outside of specified (x,z) grid\n");
            }
        }
    }

    let mut reflectors = decodeReflectors();
    checkpars();
    if smooth == 0 {
        reflectors = breakReflectors(&reflectors);
    }
    let nr = reflectors.len() as i32;

    /* read the contiguous, native-endian velocity grid */
    let vfp = stdin();
    let mut vold_flat = vec![0.0f32; nx.wrapping_mul(nz) as usize];
    if fread(&mut vold_flat, &vfp) != nx.wrapping_mul(nz) as usize {
        // C passes its FILE * to %s (undefined behavior). Empty stdin was
        // observed to print an empty string; partial reads print unstable
        // FILE-internal bytes. This fallback cannot reproduce those bytes.
        err!(
            "cannot read %d velocities from file %s",
            nx.wrapping_mul(nz),
            ""
        );
    }
    let mut vold = vec![vec![0.0f32; nz as usize]; nx as usize];
    for ix in 0..nx as usize {
        for iz in 0..nz as usize {
            vold[ix][iz] = vold_flat[ix * nz as usize + iz];
        }
    }

    /* determine maximum reflector segment length and make reflectors */
    tmin = max(tmin, max(ft, dt));
    let dsmax = ((vold[0][0] / (2 * ndpfz) as f32) as f64 * ((tmin / fpeak) as f64).sqrt()) as f32;
    let r = makeref(dsmax, &reflectors);
    let mut ns = 0;
    for ir in 0..nr as usize {
        ns += r[ir].ns;
    }
    let w = makericker(fpeak, dt);
    if verbose != 0 {
        warn!("\nSUSYNVXZ:");
        warn!("Total number of small reflecting segments is %d.\n", ns);
    }

    let mut tr = Trace::new();
    tr.trid = 1;
    tr.counit = 1;
    tr.ns = nt as u16;
    tr.dt = (1.0e6 * dt as f64) as i32 as u16;
    tr.delrt = (1.0e3 * ft as f64) as i32 as i16;

    /* allocate space; each row has the same z stride as the C arrays */
    let mut nx1 = 1 + 2 * nxb;
    let mut ts = vec![vec![0.0f32; nz as usize]; nx1 as usize];
    let mut a_s = ts.clone();
    let mut sgs = ts.clone();
    let mut tg = ts.clone();
    let mut ag = ts.clone();
    let mut sgg = ts.clone();
    let mut v = ts.clone();
    let mut bas = ts.clone();
    let mut bag = ts.clone();
    let mut ts1 = Vec::new();
    let mut as1 = Vec::new();
    let mut sgs1 = Vec::new();
    let mut tg1 = Vec::new();
    let mut ag1 = Vec::new();
    let mut sgg1 = Vec::new();
    if nxd > 1 {
        ts1 = ts.clone();
        as1 = ts.clone();
        sgs1 = ts.clone();
        tg1 = ts.clone();
        ag1 = ts.clone();
        sgg1 = ts.clone();
    }
    // C makeone's static filter/cache is owned here and shared across traces.
    let mut hd = [0.0f32; NHD];
    let mut madehd = false;

    /* loop over offsets and midpoints */
    let mut tracl = 0;
    for ixo in 0..nxo {
        let xo = fxo + ixo as f32 * dxo;
        if abs(xo) > nxb as f32 * dx {
            err!("\t band NXB is too small!\n");
        }
        let mut nxd1 = nxd;
        let mut ixm = 0;
        while ixm < nxm {
            let mut xm = fxm + ixm as f32 * dxm;
            let mut xs = (xm as f64 - 0.5 * xo as f64) as f32;
            let mut xg = xs + xo;
            let mut fx1 = xm - nxb as f32 * dx;
            let mut ex1 = min(ex + (nxd1 - 1) as f32 * dxm, xm + nxb as f32 * dx);
            nx1 = (1.0f32 + (ex1 - fx1) / dx) as i32;
            let nx0 = ((fx1 - fx) / dx) as i32;
            let mut temp = (fx1 - fx) / dx - nx0 as f32;

            /* transpose velocity such that the first row is at fx1 */
            for ix in 0..nx1 {
                for iz in 0..nz as usize {
                    v[ix as usize][iz] = if ix < -nx0 {
                        vold[0][iz]
                    } else if ix + nx0 > nx - 2 {
                        vold[(nx - 1) as usize][iz]
                    } else {
                        (vold[(ix + nx0) as usize][iz] as f64 * (1.0 - temp as f64)
                            + (temp * vold[(ix + nx0 + 1) as usize][iz]) as f64)
                            as f32
                    };
                }
            }

            if ixm == 0 || nxd1 == 1 {
                /* no interpolation */
                eiktam(
                    xs,
                    0.0,
                    dz,
                    0.0,
                    dx,
                    fx1,
                    &v[..nx1 as usize]
                        .iter()
                        .map(Vec::as_slice)
                        .collect::<Vec<_>>(),
                    &mut ts[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut a_s[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut sgs[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut bas[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                );
                eiktam(
                    xg,
                    0.0,
                    dz,
                    0.0,
                    dx,
                    fx1,
                    &v[..nx1 as usize]
                        .iter()
                        .map(Vec::as_slice)
                        .collect::<Vec<_>>(),
                    &mut tg[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut ag[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut sgg[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut bag[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                );
                let mut ixd = nint(((xs - fx) / dx) as f64);
                let vs0 = vold[ixd as usize][0];
                ixd = nint(((xg - fx) / dx) as f64);
                let vg0 = vold[ixd as usize][0];
                ex1 = min(ex, xm + nxb as f32 * dx);
                makeone(
                    &ts,
                    &a_s,
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
                    &mut madehd,
                );
                tracl += 1;
                tr.tracl = tracl;
                tr.tracr = tracl;
                tr.cdp = 1 + ixm;
                tr.cdpt = 1 + ixo;
                tr.offset = nint(xo as f64);
                tr.sx = nint(xs as f64);
                tr.gx = nint(xg as f64);
                puttr(&mut tr);
            } else {
                /* linear interpolation */
                eiktam(
                    xs,
                    0.0,
                    dz,
                    0.0,
                    dx,
                    fx1,
                    &v[..nx1 as usize]
                        .iter()
                        .map(Vec::as_slice)
                        .collect::<Vec<_>>(),
                    &mut ts1[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut as1[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut sgs1[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut bas[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                );
                eiktam(
                    xg,
                    0.0,
                    dz,
                    0.0,
                    dx,
                    fx1,
                    &v[..nx1 as usize]
                        .iter()
                        .map(Vec::as_slice)
                        .collect::<Vec<_>>(),
                    &mut tg1[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut ag1[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut sgg1[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                    &mut bag[..nx1 as usize]
                        .iter_mut()
                        .map(Vec::as_mut_slice)
                        .collect::<Vec<_>>(),
                );
                let mut ixd = nint(((xs - fx) / dx) as f64);
                let vs0 = vold[ixd as usize][0];
                ixd = nint(((xg - fx) / dx) as f64);
                let vg0 = vold[ixd as usize][0];
                xm -= nxd1 as f32 * dxm;
                for i in 1..=nxd1 {
                    xm += dxm;
                    xs = (xm as f64 - 0.5 * xo as f64) as f32;
                    xg = xs + xo;
                    fx1 = xm - nxb as f32 * dx;
                    ex1 = min(ex + (nxd1 - 1) as f32 * dxm, xm + nxb as f32 * dx);
                    nx1 = (1.0f32 + (ex1 - fx1) / dx) as i32;
                    temp = (nxd1 - i) as f32;
                    let temp1 = (1.0 / (nxd1 - i + 1) as f64) as f32;
                    for ix in 0..nx1 as usize {
                        for iz in 0..nz as usize {
                            if i == nxd1 {
                                ts[ix][iz] = ts1[ix][iz];
                                tg[ix][iz] = tg1[ix][iz];
                                sgs[ix][iz] = sgs1[ix][iz];
                                ag[ix][iz] = ag1[ix][iz];
                                sgg[ix][iz] = sgg1[ix][iz];
                                a_s[ix][iz] = as1[ix][iz];
                            } else {
                                ts[ix][iz] = (temp * ts[ix][iz] + ts1[ix][iz]) * temp1;
                                tg[ix][iz] = (temp * tg[ix][iz] + tg1[ix][iz]) * temp1;
                                a_s[ix][iz] = (temp * a_s[ix][iz] + as1[ix][iz]) * temp1;
                                sgs[ix][iz] = (temp * sgs[ix][iz] + sgs1[ix][iz]) * temp1;
                                ag[ix][iz] = (temp * ag[ix][iz] + ag1[ix][iz]) * temp1;
                                sgg[ix][iz] = (temp * sgg[ix][iz] + sgg1[ix][iz]) * temp1;
                            }
                        }
                    }
                    ex1 = min(ex, xm + nxb as f32 * dx);
                    makeone(
                        &ts,
                        &a_s,
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
                        &mut madehd,
                    );
                    tracl += 1;
                    tr.tracl = tracl;
                    tr.tracr = tracl;
                    tr.cdp = 1 + ixm - nxd1 + i;
                    tr.cdpt = 1 + ixo;
                    tr.offset = nint(xo as f64);
                    tr.d2 = dxm;
                    tr.f2 = fxm;
                    tr.sx = nint(xs as f64);
                    tr.gx = nint(xg as f64);
                    puttr(&mut tr);
                }
            }
            /* set skip parameter */
            if ixm < nxm - 1 && ixm > nxm - 1 - nxd1 {
                nxd1 = nxm - 1 - ixm;
            }
            ixm += nxd1;
        }
        warn!("\t finish offset %f\n", xo);
    }
    exit(CWP_Exit());
}

/* Make one synthetic seismogram (the program's static C helper). */
#[allow(clippy::too_many_arguments)]
fn makeone(
    ts: &[Vec<f32>],
    a_s: &[Vec<f32>],
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
    hd: &mut [f32; NHD],
    madehd: &mut bool,
) {
    let trace = &mut trace[..nt as usize];
    if !*madehd {
        mkhdiff(dt, LHD, hd);
        *madehd = true;
    }
    for it in 0..nt as usize {
        trace[it] = 0.0;
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

            /* bilinear interpolation, retaining C's mixed f32/f64 products */
            macro_rules! bilinear {
                ($a:expr) => {{
                    let a = $a;
                    let ix = ix as usize;
                    let iz = iz as usize;
                    ((1.0 - sz as f64)
                        * ((1.0 - sx as f64) * a[ix][iz] as f64 + (sx * a[ix + 1][iz]) as f64)
                        + sz as f64
                            * ((1.0 - sx as f64) * a[ix][iz + 1] as f64
                                + (sx * a[ix + 1][iz + 1]) as f64)) as f32
                }};
            }
            let tsd = bilinear!(ts);
            let asd = bilinear!(a_s);
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
            addsinc(time, amp, dt, ft, trace);
        }
    }
    let mut temp = vec![0.0f32; nt as usize];
    convolve_cwp(-LHD, hd, 0, trace, 0, &mut temp);
    convolve_cwp(w.iw, &w.wv[..w.lw as usize], 0, &temp, 0, trace);
}
