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

//! SUSYNLVCW - synthetic seismograms for linear velocity converted waves.
//! Port of su/main/synthetics_waveforms_testpatterns/susynlvcw.c (Revision 1.18).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Based on Dave Hale's code susynlv, modified by Mohammed Alfaraj to handle
//! mode conversion. Date of modification: 01/07/92.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{max, nint};
use oxwave_f_convolution::convolve_cwp;
use oxwave_f_modeling::{
    addsinc, breakReflectors, decodeReflectors, makeref, makericker, raylv2, Reflector, Wavelet,
};
use oxwave_finite_difference::mkhdiff;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSYNLVCW - SYNthetic seismograms for Linear Velocity function\t",
    " \t\tfor mode Converted Waves\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " susynlvcw >outfile [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nt=101\t\tnumber of time samples\t\t\t\t",
    " dt=0.04\t\ttime sampling interval (sec)\t\t\t",
    " ft=0.0\t\tfirst time (sec)\t\t\t\t",
    " nxo=1\t\t\tnumber of source-receiver offsets\t\t",
    " dxo=0.05\t\toffset sampling interval (km)\t\t\t",
    " fxo=0.0\t\tfirst offset (km, see notes below)\t\t",
    " xo=fxo,fxo+dxo,...\tarray of offsets (use only for non-uniform offsets)",
    " nxm=101\t\tnumber of midpoints (see notes below)\t\t",
    " dxm=0.05\t\tmidpoint sampling interval (km)\t\t",
    " fxm=0.0\t\tfirst midpoint (km)\t\t\t\t",
    " nxs=101\t\tnumber of shotpoints (see notes below)\t\t",
    " dxs=0.05\t\tshotpoint sampling interval (km)\t\t",
    " fxs=0.0\t\tfirst shotpoint (km)\t\t\t\t",
    " x0=0.0\t\tdistance x at which v00 is specified\t\t",
    " z0=0.0\t\tdepth z at which v00 is specified\t\t",
    " v00=2.0\t\tvelocity at x0,z0 (km/sec)\t\t\t",
    " gamma=1.0\t\tvelocity ratio, upgoing/downgoing\t\t",
    " dvdx=0.0\t\tderivative of velocity with distance x (dv/dx)\t",
    " dvdz=0.0\t\tderivative of velocity with depth z (dv/dz)\t",
    " fpeak=0.2/dt\t\tpeak frequency of symmetric Ricker wavelet (Hz)\t",
    " ref=\"1:1,2;4,2\"\treflector(s):  \"amplitude:x1,z1;x2,z2;x3,z3;...\"",
    " smooth=0\t\t=1 for smooth (piecewise cubic spline) reflectors",
    " er=0\t\t\t=1 for exploding reflector amplitudes\t\t",
    " ls=0\t\t\t=1 for line source; default is point source\t",
    " ob=1\t\t\t=1 to include obliquity factors\t\t",
    " sp=1\t\t\t=1 to account for amplitude spreading\t\t",
    " \t\t\t=0 for constant amplitudes throught out\t\t",
    " tmin=10.0*dt\t\tminimum time of interest (sec)\t\t\t",
    " ndpfz=5\t\tnumber of diffractors per Fresnel zone\t\t",
    " verbose=0\t\t=1 to print some useful information\t\t",
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
    " Note that gamma<1 implies P-SV mode conversion, gamma>1 implies SV-P,\t",
    " and gamma=1 implies no mode conversion.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

// Parameters for the half-derivative filter.
const LHD: i32 = 20;
const NHD: usize = 1 + 2 * LHD as usize;

#[allow(unused_assignments)] // C reads midpoint defaults again, even in shot mode.
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut nt = 101;
    getparint("nt", &mut nt);
    CHECK_NT("nt", nt);
    let mut dt = 0.04f32;
    let mut ft = 0.0f32;
    getparfloat("dt", &mut dt);
    getparfloat("ft", &mut ft);
    let mut nxo = countparval("xo");
    let mut xo;
    if nxo != 0 {
        xo = vec![0.0f32; nxo as usize];
        getparfloat("xo", &mut xo);
    } else {
        let mut dxo = 0.05;
        let mut fxo = 0.0;
        if getparint("nxo", &mut nxo) == 0 {
            nxo = 1;
        }
        getparfloat("dxo", &mut dxo);
        getparfloat("fxo", &mut fxo);
        xo = Vec::<f32>::new();
        // ealloc1float receives a size_t: a negative count becomes an
        // impossible allocation. Preserve its ENOMEM error instead of
        // panicking on Vec's capacity overflow (which does not set errno).
        if xo.try_reserve_exact(nxo as usize).is_err() {
            err!(
                "%s: malloc failed (%s)",
                "ealloc.c",
                oxwave_core::err::strerror(12)
            );
        }
        xo.resize(nxo as usize, 0.0);
        for ixo in 0..nxo {
            xo[ixo as usize] = fxo + ixo as f32 * dxo;
        }
    }
    let mut nxs = 0;
    let mut dxs = 0.0;
    let mut fxs = 0.0;
    let mut nxm = 0;
    let mut dxm = 0.0;
    let mut fxm = 0.0;
    let shots = getparint("nxs", &mut nxs) != 0
        || getparfloat("dxs", &mut dxs) != 0
        || getparfloat("fxs", &mut fxs) != 0;
    let midpoints = getparint("nxm", &mut nxm) != 0
        || getparfloat("dxm", &mut dxm) != 0
        || getparfloat("fxm", &mut fxm) != 0;
    if shots && midpoints {
        err!("cannot specify both shot and midpoint sampling!\n");
    }

    let nxsm;
    let dxsm;
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
    let mut x0 = 0.0f32;
    if getparfloat("x0", &mut x0) == 0 {
        x0 = 0.0;
    }
    let mut z0 = 0.0f32;
    if getparfloat("z0", &mut z0) == 0 {
        z0 = 0.0;
    }
    let mut v00 = 0.0f32;
    if getparfloat("v00", &mut v00) == 0 {
        v00 = 2.0;
    }
    let mut gamma = 0.0f32;
    if getparfloat("gamma", &mut gamma) == 0 {
        gamma = 1.0;
    }
    let mut dvdx = 0.0f32;
    if getparfloat("dvdx", &mut dvdx) == 0 {
        dvdx = 0.0;
    }
    let mut dvdz = 0.0f32;
    if getparfloat("dvdz", &mut dvdz) == 0 {
        dvdz = 0.0;
    }
    let mut fpeak = 0.0f32;
    if getparfloat("fpeak", &mut fpeak) == 0 {
        fpeak = (0.2f64 / dt as f64) as f32;
    }
    let mut ls = 0;
    if getparint("ls", &mut ls) == 0 {
        ls = 0;
    }
    let mut er = 0;
    if getparint("er", &mut er) == 0 {
        er = 0;
    }
    let mut ob = 0;
    if getparint("ob", &mut ob) == 0 {
        ob = 1;
    }
    let mut sp = 0;
    if getparint("sp", &mut sp) == 0 {
        sp = 1;
    }
    let mut tmin = 0.0f32;
    if getparfloat("tmin", &mut tmin) == 0 {
        tmin = (10.0f64 * dt as f64) as f32;
    }
    let mut ndpfz = 0;
    if getparint("ndpfz", &mut ndpfz) == 0 {
        ndpfz = 5;
    }
    let mut smooth = 0;
    if getparint("smooth", &mut smooth) == 0 {
        smooth = 0;
    }
    let mut verbose = 0;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    let mut input = decodeReflectors();
    if smooth == 0 {
        input = breakReflectors(&input);
    }

    // Convert velocity v(x0,z0) to v(0,0).
    v00 -= dvdx * x0 + dvdz * z0;

    // Determine minimum velocity and minimum reflection time.
    let mut vmin = f32::MAX;
    let mut tminr = f32::MAX;
    for ir in 0..input.len() {
        for ixz in 0..input[ir].x.len() {
            let x = input[ir].x[ixz];
            let z = input[ir].z[ixz];
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

    // Determine maximum reflector segment length, then make reflectors.
    tmin = max(tmin, max(ft, dt));
    let dsmax =
        ((vmin / 2_i32.wrapping_mul(ndpfz) as f32) as f64 * ((tmin / fpeak) as f64).sqrt()) as f32;
    let r = makeref(dsmax, &input);

    // Count reflector segments.
    let mut ns = 0_i32;
    for ir in 0..r.len() {
        ns = ns.wrapping_add(r[ir].ns);
    }

    // Make wavelet.
    let w = makericker(fpeak, dt);

    if verbose != 0 {
        warn!("\nSYNLVCW:");
        warn!(
            "Minimum possible reflection time (assuming sources\nand receivers are at the surface Z=0) is %g s.\nYou may want to adjust the \"minimum time of \ninterest\" parameter.",
            tminr
        );
        warn!("Total number of small reflecting\nsegments is %d.", ns);
    }

    // Set constant segy trace header parameters.
    let mut tr = Trace::new();
    tr.trid = 1;
    tr.counit = 1;
    tr.ns = nt as u16;
    // C converts through a machine integer before storing the short fields.
    tr.dt = (1.0e6 * dt as f64) as u32 as u16;
    tr.delrt = (1.0e3 * ft as f64) as i32 as i16;
    let mut hd = None;

    // Loop over shots or midpoints, then offsets.
    let mut tracl = 0_i32;
    for ixsm in 0..nxsm {
        for ixo in 0..nxo {
            let xo = xo[ixo as usize];
            let xs = if shots {
                fxs + ixsm as f32 * dxs
            } else {
                ((fxm + ixsm as f32 * dxm) as f64 - 0.5 * xo as f64) as f32
            };
            let zs = 0.0;
            let xg = xs + xo;
            let zg = 0.0;

            // Set segy trace header parameters.
            tracl = tracl.wrapping_add(1);
            tr.tracl = tracl;
            tr.tracr = tracl;
            if shots {
                tr.fldr = 1_i32.wrapping_add(ixsm);
                tr.tracf = 1_i32.wrapping_add(ixo);
            } else {
                tr.cdp = 1_i32.wrapping_add(ixsm);
                tr.cdpt = 1_i32.wrapping_add(ixo);
            }
            tr.offset = nint(1000.0 * if dxsm > 0.0 { xo as f64 } else { -xo as f64 });
            tr.sx = nint(1000.0 * xs as f64);
            tr.gx = nint(1000.0 * xg as f64);

            makeone(
                v00,
                dvdx,
                dvdz,
                ls,
                er,
                ob,
                sp,
                &w,
                gamma,
                xs,
                zs,
                xg,
                zg,
                &r,
                nt,
                dt,
                ft,
                &mut tr.data,
                &mut hd,
            );
            puttr(&mut tr);
        }
    }
    exit(CWP_Exit());
}

/// Make one synthetic seismogram for v(x,z) = v00+dvdx*x+dvdz*z.
/// `hd` carries the C function's static half-derivative filter between traces.
#[allow(clippy::too_many_arguments)]
fn makeone(
    v00: f32,
    dvdx: f32,
    dvdz: f32,
    ls: i32,
    er: i32,
    ob: i32,
    sp: i32,
    w: &Wavelet,
    gamma: f32,
    xs: f32,
    zs: f32,
    xg: f32,
    zg: f32,
    r: &[Reflector],
    nt: i32,
    dt: f32,
    ft: f32,
    trace: &mut [f32],
    hd: &mut Option<[f32; NHD]>,
) {
    let tos = (2.0 / (1.0 + 1.0 / gamma as f64)) as f32;

    // If half-derivative filter not yet made, make it.
    let hd = hd.get_or_insert_with(|| {
        let mut hd = [0.0; NHD];
        mkhdiff(dt, LHD, &mut hd);
        hd
    });

    // Negative nt makes the C loops empty, then fails workspace allocation.
    let trace = &mut trace[..nt.max(0) as usize];
    for it in 0..nt {
        trace[it as usize] = 0.0;
    }

    // Velocities at source and receiver.
    let vs = v00 + dvdx * xs + dvdz * zs;
    let vg = (v00 + dvdx * xg + dvdz * zg) * gamma;

    // Loop over reflectors and their diffracting segments.
    for ir in 0..r.len() {
        let ar = r[ir].a;
        let ns = r[ir].ns;
        let ds = r[ir].ds;
        let rs = &r[ir].rs;
        for is in 0..ns as usize {
            let xd = rs[is].x;
            let zd = rs[is].z;
            let cd = rs[is].c;
            let sd = rs[is].s;

            // Velocity at diffractor.
            let vd = (v00 + dvdx * xd + dvdz * zd) * tos;

            // Rays from shot and receiver to diffractor.
            let (cs, ss, ts, qs) = raylv2(v00, dvdx, dvdz, xs, zs, xd, zd);
            let (cg, sg, tg, qg) = raylv2(gamma * v00, gamma * dvdx, gamma * dvdz, xg, zg, xd, zd);

            // Cosines of incidence and reflection angles.
            let (ci, cr) = if ob != 0 {
                (cd * cs + sd * ss, cd * cg + sd * sg)
            } else {
                (1.0, 1.0)
            };
            if ci < 0.0 || cr < 0.0 {
                continue;
            }

            // Two-way time and amplitude.
            let time = ts + tg;
            let amp = if sp != 0 {
                let mut amp = if er != 0 {
                    ((vg * vd / qg) as f64).sqrt() as f32
                } else if ls != 0 {
                    ((vs * vd * vd * vg / (qs * qg)) as f64).sqrt() as f32
                } else {
                    ((vs * vd * vd * vg / (qs * qg * (qs + qg))) as f64).sqrt() as f32
                };
                amp *= (ci + cr) * ar * ds;
                amp
            } else {
                (1.0 / (time as f64).sqrt()) as f32
            };

            // modeling.c:520-540 truncates negative tn toward zero, so an
            // overlapping wavelet can read a negative sinc row. Match the
            // bounds-checked C's clean error before calling the shared routine.
            let tn = (time - ft) / dt;
            let jt = tn as i32;
            let jsinc = ((tn - jt as f32) * 100.0) as i32;
            if jt.wrapping_add(4) >= 0
                && jt.wrapping_sub(3) < nt
                && !(0..101).contains(&jsinc)
            {
                err!("addsinc: sinc index out of bounds");
            }
            addsinc(time, amp, dt, ft, trace);
        }
    }

    // Allocate workspace, preserving ealloc1float's negative-count failure.
    let mut temp = Vec::new();
    if temp.try_reserve_exact(nt as usize).is_err() {
        err!(
            "%s: malloc failed (%s)",
            "ealloc.c",
            oxwave_core::err::strerror(12)
        );
    }
    temp.resize(nt as usize, 0.0);

    // Apply the half-derivative filter, then convolve the wavelet with trace.
    convolve_cwp(-LHD, hd, 0, trace, 0, &mut temp);
    convolve_cwp(w.iw, &w.wv, 0, &temp, 0, trace);
}
