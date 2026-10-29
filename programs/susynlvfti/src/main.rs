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

//! SUSYNLVFTI - synthetic seismograms in a factorized transversely isotropic medium.
//! Port of su/main/synthetics_waveforms_testpatterns/susynlvfti.c, Revision 1.8.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: Dave Hale, 09/17/91, Colorado School of Mines.
//! Upgrade to Transversely Isotropic medium: T. Alkhallfah, 4/15/93.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{max, nint, PI};
use oxwave_f_convolution::convolve_cwp;
use oxwave_f_modeling::{
    addsinc, breakReflectors, decodeReflectors, makeref, makericker, raylv2, Reflector, Wavelet,
};
use oxwave_finite_difference::mkhdiff;

mod ray;
use ray::raylvt;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSYNLVFTI - SYNthetic seismograms for Linear Velocity function in a  ",
    "              Factorized Transversely Isotropic medium\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " susynlvfti >outfile [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nt=101\t\tnumber of time samples\t\t\t\t",
    " dt=0.04\t\ttime sampling interval (sec)\t\t\t",
    " ft=0.0\t\tfirst time (sec)\t\t\t\t",
    " kilounits=1            input length units are km or kilo-feet         ",
    "                        =0 for m or ft                                 ",
    "                        Note: Output (sx,gx,offset) are always m or ft ",
    " nxo=1\t\t\tnumber of source-receiver offsets\t\t",
    " dxo=0.05\t\toffset sampling interval (kilounits)\t\t",
    " fxo=0.0\t\tfirst offset (kilounits, see notes below)\t",
    " xo=fxo,fxo+dxo,...    array of offsets (use only for non-uniform offsets)",
    " nxm=101\t\tnumber of midpoints (see notes below)\t\t",
    " dxm=0.05\t\tmidpoint sampling interval (kilounits)\t\t",
    " fxm=0.0\t\tfirst midpoint (kilounits)\t\t\t",
    " nxs=101\t\tnumber of shotpoints (see notes below)\t\t",
    " dxs=0.05\t\tshotpoint sampling interval (kilounits)\t\t",
    " fxs=0.0\t\tfirst shotpoint (kilounits)\t\t\t",
    " x0=0.0\t\tdistance x at which v00 is specified\t\t",
    " z0=0.0\t\tdepth z at which v00 is specified\t\t",
    " v00=2.0\t\tvelocity at x0,z0 (kilounits/sec)\t\t",
    " dvdx=0.0\t\tderivative of velocity with distance x (dv/dx)\t",
    " dvdz=0.0\t\tderivative of velocity with depth z (dv/dz)\t",
    " fpeak=0.2/dt\t\tpeak frequency of symmetric Ricker wavelet (Hz)\t",
    " ref=1:1,2;4,2\t\treflector(s):  \"amplitude:x1,z1;x2,z2;x3,z3;...\"",
    " smooth=0\t\t=1 for smooth (piecewise cubic spline) reflectors",
    " er=0\t\t\t=1 for exploding reflector amplitudes\t\t",
    " ls=0\t\t\t=1 for line source; default is point source\t",
    " ob=0\t\t\t=1 to include obliquity factors\t\t\t",
    " tmin=10.0*dt\t\tminimum time of interest (sec)\t\t\t",
    " ndpfz=5\t\tnumber of diffractors per Fresnel zone\t\t",
    " verbose=1\t\t=1 to print some useful information\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " For transversely isotropic media:\t\t\t\t\t",
    " angxs=0.0\t\tangle of symmetry axis with the vertical (degrees)",
    " define the media using either\t\t\t\t\t\t",
    " a=1.0\t\tcorresponding to the ratio of elastic coef.(c1111/c3333)",
    " f=0.4\t\tcorresponding to the ratio of elastic coef. (c1133/c3333)",
    " l=0.3\t\tcorresponding to the ratio of elastic coef. (c1313/c3333)",
    " Alternately use Tompson's parameters:\t\t\t\t",
    " delta=0\tThomsen's 1986 defined parameter\t\t\t",
    " epsilon=0\tThomsen's 1986 defined parameter\t\t\t",
    " ntries=40\tnumber of iterations in Snell's law and offset searches ",
    " epsx=.001\tlateral offset tolerance\t\t\t\t",
    " epst=.0001\treflection time tolerance\t\t\t\t",
    " nitmax=12\tmax number of iterations in travel time integrations\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Offsets are signed - may be positive or negative.  Receiver locations\t",
    " are computed by adding the signed offset to the source location.\t",
    " \t\t\t\t\t\t\t\t\t",
    " Specify either midpoint sampling or shotpoint sampling, but not both.\t",
    " If neither is specified, the default is the midpoint sampling above.\t",
    " \t\t\t\t\t\t\t\t\t",
    " More than one ref (reflector) may be specified.  When obliquity factors",
    " are included, then only the left side of each reflector (as the x,z\t",
    " reflector coordinates are traversed) is reflecting.  For example, if x",
    " coordinates increase, then the top side of a reflector is reflecting.\t",
    " Note that reflectors are encoded as quoted strings, with an optional\t",
    " reflector amplitude: preceding the x,z coordinates of each reflector.\t",
    " Default amplitude is 1.0 if amplitude: part of the string is omitted.\t",
    " \t\t\t\t\t\t\t\t\t",
    " Concerning the choice of delta and epsilon. The difference between delta",
    " and epsilon should not exceed one. A possible break down of the program",
    " is the result. This is caused primarly by the break down in the two point",
    " ray-tracing. Also keep the values of delta and epsilon between 2 and -2.",
];

#[allow(unused_assignments)] // Preserve the C's repeated nxm parameter read.
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut nt = 101;
    let mut dt = 0.04f32;
    let mut ft = 0.0f32;
    let mut kilounits = 1;
    if getparint("nt", &mut nt) == 0 {
        nt = 101;
    }
    if getparfloat("dt", &mut dt) == 0 {
        dt = 0.04;
    }
    if getparfloat("ft", &mut ft) == 0 {
        ft = 0.0;
    }
    if getparint("kilounits", &mut kilounits) == 0 {
        kilounits = 1;
    }
    let mut nxo = countparval("xo");
    // C leaves these uninitialized for explicit xo. The reference arm64
    // executable consistently leaves bit pattern 1 in both slots, including
    // when dxo/fxo arguments are present. Preserve its headers and zeroff test.
    let mut dxo = f32::from_bits(1);
    let mut fxo = f32::from_bits(1);
    let mut xo;
    if nxo != 0 {
        xo = vec![0.0f32; nxo as usize];
        getparfloat("xo", &mut xo);
    } else {
        if getparint("nxo", &mut nxo) == 0 {
            nxo = 1;
        }
        if getparfloat("dxo", &mut dxo) == 0 {
            dxo = 0.05;
        }
        if getparfloat("fxo", &mut fxo) == 0 {
            fxo = 0.0;
        }
        xo = vec![0.0f32; nxo as usize];
        for ixo in 0..nxo {
            xo[ixo as usize] = fxo + ixo as f32 * dxo;
        }
    }
    let (mut nxs, mut nxm) = (0, 0);
    let (mut dxs, mut fxs, mut dxm, mut fxm) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let shots = getparint("nxs", &mut nxs) != 0
        || getparfloat("dxs", &mut dxs) != 0
        || getparfloat("fxs", &mut fxs) != 0;
    let midpoints = getparint("nxm", &mut nxm) != 0
        || getparfloat("dxm", &mut dxm) != 0
        || getparfloat("fxm", &mut fxm) != 0;
    if shots && midpoints {
        err!("cannot specify both shot and midpoint sampling!\n");
    }
    let (nxsm, dxsm);
    if shots {
        if getparint("nxs", &mut nxs) == 0 {
            nxs = 101;
        }
        if getparfloat("dxs", &mut dxs) == 0 {
            dxs = 0.05;
        }
        if getparfloat("fxs", &mut fxs) == 0 {
            fxs = 0.0;
        }
        nxsm = nxs;
        dxsm = dxs;
    } else {
        if getparint("nxm", &mut nxm) == 0 {
            nxm = 101;
        }
        if getparfloat("dxm", &mut dxm) == 0 {
            dxm = 0.05;
        }
        if getparfloat("fxm", &mut fxm) == 0 {
            fxm = 0.0;
        }
        nxsm = nxm;
        dxsm = dxm;
    }
    if getparint("nxm", &mut nxm) == 0 {
        nxm = 101;
    }
    if getparfloat("dxm", &mut dxm) == 0 {
        dxm = 0.05;
    }
    if getparfloat("fxm", &mut fxm) == 0 {
        fxm = 0.0;
    }
    let (mut x0, mut z0, mut v00, mut dvdx, mut dvdz) = (0.0f32, 0.0f32, 2.0f32, 0.0f32, 0.0f32);
    if getparfloat("x0", &mut x0) == 0 {
        x0 = 0.0;
    }
    if getparfloat("z0", &mut z0) == 0 {
        z0 = 0.0;
    }
    if getparfloat("v00", &mut v00) == 0 {
        v00 = 2.0;
    }
    if getparfloat("dvdx", &mut dvdx) == 0 {
        dvdx = 0.0;
    }
    if getparfloat("dvdz", &mut dvdz) == 0 {
        dvdz = 0.0;
    }
    let mut fpeak = 0.0f32;
    if getparfloat("fpeak", &mut fpeak) == 0 {
        fpeak = (0.2 / dt as f64) as f32;
    }
    let (mut ls, mut er, mut ob) = (0, 0, 0);
    if getparint("ls", &mut ls) == 0 {
        ls = 0;
    }
    if getparint("er", &mut er) == 0 {
        er = 0;
    }
    if getparint("ob", &mut ob) == 0 {
        ob = 0;
    }
    let mut tmin = 0.0f32;
    if getparfloat("tmin", &mut tmin) == 0 {
        tmin = (10.0 * dt as f64) as f32;
    }
    let (mut ndpfz, mut smooth, mut verbose) = (5, 0, 1);
    if getparint("ndpfz", &mut ndpfz) == 0 {
        ndpfz = 5;
    }
    if getparint("smooth", &mut smooth) == 0 {
        smooth = 0;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 1;
    }
    let mut input = decodeReflectors();
    if smooth == 0 {
        input = breakReflectors(&input);
    }

    /* Optional parameters relating to anisotropy. */
    let (mut angxs, mut a, mut f, mut l, mut delta, mut epsilon) =
        (0.0f32, 1.0f32, 0.4f32, 0.3f32, 0.0f32, 0.0f32);
    if getparfloat("angxs", &mut angxs) == 0 {
        angxs = 0.0;
    }
    if getparfloat("a", &mut a) == 0 {
        a = 1.0;
    }
    if getparfloat("f", &mut f) == 0 {
        f = 0.4;
    }
    if getparfloat("l", &mut l) == 0 {
        l = 0.3;
    }
    if getparfloat("delta", &mut delta) == 0 {
        delta = 0.0;
    }
    if getparfloat("epsilon", &mut epsilon) == 0 {
        epsilon = 0.0;
    }
    let (mut ntries, mut nitmax) = (20, 12);
    let (mut epsx, mut epst) = (0.001f32, 0.0001f32);
    if getparint("ntries", &mut ntries) == 0 {
        ntries = 20;
    }
    if getparfloat("epsx", &mut epsx) == 0 {
        epsx = 0.001;
    }
    if getparfloat("epst", &mut epst) == 0 {
        epst = 0.0001;
    }
    if getparint("nitmax", &mut nitmax) == 0 {
        nitmax = 12;
    }
    let mut trans = a != 1.0 || f + 2.0f32 * l != 1.0;
    if delta != 0.0 || epsilon != 0.0 {
        a = 1.0f32 + 2.0f32 * epsilon;
        // sqrt() returns double, so the subtraction rounds only on assignment.
        f = (((2.0f32 * delta * (1.0f32 - l) + (1.0f32 - l) * (1.0f32 - l)) as f64).sqrt()
            - l as f64) as f32;
        trans = true;
    }
    let zeroff = nxo == 1 && fxo == 0.0;
    v00 -= dvdx * x0 + dvdz * z0;
    let torad = (PI / 180.0) as f32;
    let mut vmin = f32::MAX;
    let mut tminr = f32::MAX;
    for reflector in &input {
        for ixz in 0..reflector.x.len() {
            let x = reflector.x[ixz];
            let z = reflector.z[ixz];
            let v = v00 + dvdx * x + dvdz * z;
            if v < vmin {
                vmin = v;
            }
            let t = (2.0 * z as f64 / v as f64) as f32;
            if t < tminr {
                tminr = t;
            }
        }
    }
    tmin = max(tmin, max(ft, dt));
    let dsmax = ((vmin / (2 * ndpfz) as f32) as f64 * ((tmin / fpeak) as f64).sqrt()) as f32;
    let r = makeref(dsmax, &input);
    let mut ns = 0i32;
    for reflector in &r {
        ns += reflector.ns;
    }
    let w = makericker(fpeak, dt);
    if verbose != 0 {
        fprintf!(stderr(), "\nSYNLVXZ:\n");
        fprintf!(stderr(), "Minimum possible reflection time (assuming sources\nand receivers are at the surface Z=0) is %g s.\nYou may want to adjust the \"minimum time of \ninterest\" parameter.\n", tminr);
        fprintf!(
            stderr(),
            "Total number of small reflecting\nsegments is %d.\n",
            ns
        );
        fprintf!(stderr(), "\n");
    }
    tr.trid = 1;
    tr.counit = 1;
    tr.ns = nt as u16;
    tr.dt = (1.0e6 * dt as f64) as i32 as u16;
    tr.delrt = (1.0e3 * ft as f64) as i32 as i16;
    let mut tracl = 0;
    let mut hd = [0.0f32; NHD];
    let mut madehd = false;
    for ixsm in 0..nxsm {
        for ixo in 0..nxo {
            let xs = if shots {
                fxs + ixsm as f32 * dxs
            } else {
                ((fxm + ixsm as f32 * dxm) as f64 - 0.5 * xo[ixo as usize] as f64) as f32
            };
            let zs = 0.0f32;
            let xg = xs + xo[ixo as usize];
            let zg = 0.0f32;
            tracl += 1;
            tr.tracl = tracl;
            tr.tracr = tracl;
            if shots {
                tr.fldr = 1 + ixsm;
                tr.tracf = 1 + ixo;
                tr.d2 = dxo;
                tr.f2 = fxo;
            } else {
                tr.cdp = 1 + ixsm;
                tr.cdpt = 1 + ixo;
                tr.d2 = dxm;
                tr.f2 = fxm;
            }
            let offset = if dxsm > 0.0 {
                xo[ixo as usize]
            } else {
                -xo[ixo as usize]
            };
            if kilounits == 1 {
                tr.offset = nint(1000.0 * offset as f64);
                tr.sx = nint(1000.0 * xs as f64);
                tr.gx = nint(1000.0 * xg as f64);
            } else {
                tr.offset = nint(offset as f64);
                tr.sx = nint(xs as f64);
                tr.gx = nint(xg as f64);
            }
            makeone(
                v00,
                dvdx,
                dvdz,
                ls,
                er,
                ob,
                &w,
                trans,
                nitmax,
                epst,
                zeroff,
                ntries,
                epsx,
                angxs * torad,
                xs,
                zs,
                xg,
                zg,
                a,
                f,
                l,
                &r,
                nt,
                dt,
                ft,
                &mut tr.data,
                &mut hd,
                &mut madehd,
            );
            puttr(&mut tr);
        }
    }
    exit(CWP_Exit());
}

const LHD: i32 = 20;
const NHD: usize = 1 + 2 * LHD as usize;

/// Make one synthetic seismogram for v(x,z) = v00 + dvdx*x + dvdz*z.
#[allow(clippy::too_many_arguments)]
fn makeone(
    v00: f32,
    dvdx: f32,
    dvdz: f32,
    ls: i32,
    er: i32,
    ob: i32,
    w: &Wavelet,
    trans: bool,
    nitmax: i32,
    epst: f32,
    zeroff: bool,
    ntries: i32,
    epsx: f32,
    angxs: f32,
    xs: f32,
    zs: f32,
    xg: f32,
    zg: f32,
    a: f32,
    f: f32,
    l: f32,
    r: &[Reflector],
    nt: i32,
    dt: f32,
    ft: f32,
    trace: &mut [f32],
    hd: &mut [f32; NHD],
    madehd: &mut bool,
) {
    if !*madehd {
        mkhdiff(dt, LHD, hd);
        *madehd = true;
    }
    for it in 0..nt {
        trace[it as usize] = 0.0;
    }
    let vs = v00 + dvdx * xs + dvdz * zs;
    let vg = v00 + dvdx * xg + dvdz * zg;
    for reflector in r {
        let ar = reflector.a;
        let ns = reflector.ns;
        let ds = reflector.ds;
        let rs = &reflector.rs;
        for is in 0..ns as usize {
            let xd = rs[is].x;
            let zd = rs[is].z;
            let cd = rs[is].c;
            let sd = rs[is].s;
            let vd = v00 + dvdx * xd + dvdz * zd;
            let (mut cs, mut ss, mut ts, mut qs) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
            let (mut cg, mut sg, mut tg, mut qg) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
            if !zeroff {
                if trans {
                    raylvt(
                        v00,
                        dvdx,
                        dvdz,
                        xs,
                        zs,
                        xd,
                        zd,
                        a,
                        f,
                        l,
                        ntries,
                        nitmax as f32,
                        epst,
                        epsx,
                        angxs,
                        &mut cs,
                        &mut ss,
                        &mut ts,
                        &mut qs,
                    );
                    raylvt(
                        v00,
                        dvdx,
                        dvdz,
                        xg,
                        zg,
                        xd,
                        zd,
                        a,
                        f,
                        l,
                        ntries,
                        nitmax as f32,
                        epst,
                        epsx,
                        angxs,
                        &mut cg,
                        &mut sg,
                        &mut tg,
                        &mut qg,
                    );
                } else {
                    (cs, ss, ts, qs) = raylv2(v00, dvdx, dvdz, xs, zs, xd, zd);
                    (cg, sg, tg, qg) = raylv2(v00, dvdx, dvdz, xg, zg, xd, zd);
                }
            } else {
                if trans {
                    raylvt(
                        v00,
                        dvdx,
                        dvdz,
                        xs,
                        zs,
                        xd,
                        zd,
                        a,
                        f,
                        l,
                        ntries,
                        nitmax as f32,
                        epst,
                        epsx,
                        angxs,
                        &mut cs,
                        &mut ss,
                        &mut ts,
                        &mut qs,
                    );
                } else {
                    (cs, ss, ts, qs) = raylv2(v00, dvdx, dvdz, xs, zs, xd, zd);
                }
                cg = cs;
                sg = ss;
                tg = ts;
                qg = qs;
            }
            let (ci, cr) = if ob != 0 {
                (cd * cs + sd * ss, cd * cg + sd * sg)
            } else {
                (1.0f32, 1.0f32)
            };
            if ci < 0.0 || cr < 0.0 {
                continue;
            }
            let time = ts + tg;
            let mut amp = if er != 0 {
                ((vg * vd / qg) as f64).sqrt() as f32
            } else if ls != 0 {
                ((vs * vd * vd * vg / (qs * qg)) as f64).sqrt() as f32
            } else {
                ((vs * vd * vd * vg / (qs * qg * (qs + qg))) as f64).sqrt() as f32
            };
            amp *= (ci + cr) * ar * ds;
            // In C, a saturated sample index wraps jt +/- 3/4 into the
            // unchecked interior branch of addsinc, causing SIGBUS (also for
            // nt=0). Reproduce that failure without an invalid memory access.
            // In particular, a constant-velocity TI ray with ntries<=0 has
            // vcp=0 and an infinite travel time. NaN instead converts to 0.
            let jt = ((time - ft) / dt) as i32;
            if jt > i32::MAX - 4 || jt < i32::MIN + 3 {
                reference_sigbus();
            }
            addsinc(time, amp, dt, ft, &mut trace[..nt as usize]);
        }
    }
    let mut temp = vec![0.0f32; nt as usize];
    convolve_cwp(-LHD, hd, 0, &trace[..nt as usize], 0, &mut temp);
    convolve_cwp(w.iw, &w.wv, 0, &temp, 0, &mut trace[..nt as usize]);
}

/// Terminate with the reference's signal, retaining its unflushed stdout.
fn reference_sigbus() -> ! {
    unsafe extern "C" {
        fn signal(sig: i32, handler: usize) -> usize;
        fn raise(sig: i32) -> i32;
    }
    #[cfg(target_vendor = "apple")]
    const SIGBUS: i32 = 10;
    #[cfg(not(target_vendor = "apple"))]
    const SIGBUS: i32 = 7;
    unsafe {
        signal(SIGBUS, 0);
        raise(SIGBUS);
    }
    unreachable!("SIGBUS has the default disposition");
}
