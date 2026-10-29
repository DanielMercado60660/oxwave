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

//! SUSYNLV - synthetic seismograms for a linear velocity function.
//! Port of su/main/synthetics_waveforms_testpatterns/susynlv.c (Revision 1.22).
//! Credits: CWP Dave Hale, 09/17/91, Colorado School of Mines.
//! UTulsa Chris Liner, 5/22/03, added kilounits flag.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{max, nint};
use oxwave_f_convolution::*;
use oxwave_f_modeling::*;
use oxwave_finite_difference::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSYNLV - SYNthetic seismograms for Linear Velocity function\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " susynlv >outfile [optional parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nt=101                 number of time samples\t\t\t\t",
    " dt=0.04                time sampling interval (sec)\t\t\t",
    " ft=0.0                 first time (sec)\t\t\t\t",
    " kilounits=1            input length units are km or kilo-feet\t\t",
    "\t\t\t =0 for m or ft\t\t\t\t\t",
    "                        Note: Output (sx,gx,offset) are always m or ft ",
    " nxo=1                  number of source-receiver offsets\t\t",
    " dxo=0.05               offset sampling interval (kilounits)\t\t",
    " fxo=0.0                first offset (kilounits, see notes below)\t",
    " xo=fxo,fxo+dxo,...     array of offsets (use only for non-uniform offsets)",
    " nxm=101                number of midpoints (see notes below)\t\t",
    " dxm=0.05               midpoint sampling interval (kilounits)\t\t",
    " fxm=0.0                first midpoint (kilounits)\t\t\t",
    " nxs=101                number of shotpoints (see notes below)\t\t",
    " dxs=0.05               shotpoint sampling interval (kilounits)\t",
    " fxs=0.0                first shotpoint (kilounits)\t\t\t",
    " x0=0.0                 distance x at which v00 is specified\t\t",
    " z0=0.0                 depth z at which v00 is specified\t\t",
    " v00=2.0                velocity at x0,z0 (kilounits/sec)\t\t",
    " dvdx=0.0               derivative of velocity with distance x (dv/dx)\t",
    " dvdz=0.0               derivative of velocity with depth z (dv/dz)\t",
    " fpeak=0.2/dt           peak frequency of symmetric Ricker wavelet (Hz)",
    " ref=\"1:1,2;4,2\"        reflector(s):  \"amplitude:x1,z1;x2,z2;x3,z3;...\"",
    " smooth=0               =1 for smooth (piecewise cubic spline) reflectors",
    " er=0                   =1 for exploding reflector amplitudes\t\t",
    " ls=0                   =1 for line source; default is point source\t",
    " ob=1                   =1 to include obliquity factors\t\t",
    " tmin=10.0*dt           minimum time of interest (sec)\t\t\t",
    " ndpfz=5                number of diffractors per Fresnel zone\t\t",
    " verbose=0              =1 to print some useful information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Notes:\t\t\t\t\t\t\t\t",
    "Offsets are signed - may be positive or negative.  Receiver locations\t",
    "are computed by adding the signed offset to the source location.\t",
    "\t\t\t\t\t\t\t\t\t",
    "Specify either midpoint sampling or shotpoint sampling, but not both.\t",
    "If neither is specified, the default is the midpoint sampling above.\t",
    "\t\t\t\t\t\t\t\t\t",
    "More than one ref (reflector) may be specified. Do this by putting\t",
    "additional ref= entries on the commandline. When obliquity factors\t",
    "are included, then only the left side of each reflector (as the x,z\t",
    "reflector coordinates are traversed) is reflecting.  For example, if x\t",
    "coordinates increase, then the top side of a reflector is reflecting.\t",
    "Note that reflectors are encoded as quoted strings, with an optional\t",
    "reflector amplitude: preceding the x,z coordinates of each reflector.\t",
    "Default amplitude is 1.0 if amplitude: part of the string is omitted.\t",
];

// Parameters for the half-derivative filter.
const LHD: i32 = 20;
const NHD: usize = 1 + 2 * LHD as usize;

#[allow(unused_assignments)] // Retain C's repeated nxm read and midpoints assignment.
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    // Get parameters in the same order as susynlv.c, including the
    // short-circuit geometry detection and subsequent repeated reads.
    let mut nt = 0;
    if getparint("nt", &mut nt) == 0 {
        nt = 101;
    }
    CHECK_NT("nt", nt);
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = 0.04;
    }
    let mut ft = 0.0f32;
    if getparfloat("ft", &mut ft) == 0 {
        ft = 0.0;
    }
    let mut kilounits = 0;
    if getparint("kilounits", &mut kilounits) == 0 {
        kilounits = 1;
    }
    let mut nxo = countparval("xo");
    // C leaves dxo/fxo uninitialized for explicit xo in shot mode; its
    // resulting d2/f2 headers are unstable. Use zero for those unset fields.
    let mut dxo = 0.0f32;
    let mut fxo = 0.0f32;
    let mut xo = Vec::<f32>::new();
    if nxo != 0 {
        xo.resize(nxo as usize, 0.0);
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
        for ixo in 0..nxo {
            xo.push(fxo + ixo as f32 * dxo);
        }
    }
    let mut nxs = 0;
    let mut dxs = 0.0f32;
    let mut fxs = 0.0f32;
    let mut nxm = 0;
    let mut dxm = 0.0f32;
    let mut fxm = 0.0f32;
    let shots = getparint("nxs", &mut nxs) != 0
        || getparfloat("dxs", &mut dxs) != 0
        || getparfloat("fxs", &mut fxs) != 0;
    let mut midpoints = getparint("nxm", &mut nxm) != 0
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
        midpoints = true;
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
        fpeak = (0.2 / dt as f64) as f32;
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
    let mut tmin = 0.0f32;
    if getparfloat("tmin", &mut tmin) == 0 {
        tmin = (10.0 * dt as f64) as f32;
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

    let mut reflectors = decodeReflectors();
    if smooth == 0 {
        reflectors = breakReflectors(&reflectors);
    }

    // Convert velocity v(x0,z0) to v(0,0).
    v00 -= dvdx * x0 + dvdz * z0;

    // Determine minimum velocity and minimum reflection time.
    let mut vmin = f32::MAX;
    let mut tminr = f32::MAX;
    for ir in 0..reflectors.len() {
        for ixz in 0..reflectors[ir].x.len() {
            let x = reflectors[ir].x[ixz];
            let z = reflectors[ir].z[ixz];
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

    // Determine maximum reflector segment length.
    tmin = max(tmin, max(ft, dt));
    let dsmax =
        ((vmin / 2i32.wrapping_mul(ndpfz) as f32) as f64 * ((tmin / fpeak) as f64).sqrt()) as f32;

    // Make reflectors and count reflector segments.
    let r = makeref(dsmax, &reflectors);
    let mut ns = 0i32;
    for ir in 0..r.len() {
        ns = ns.wrapping_add(r[ir].ns);
    }

    // Make wavelet.
    let w = makericker(fpeak, dt);

    if verbose != 0 {
        warn!("\nSYNLV:");
        warn!(
            concat!(
                "Minimum possible reflection time (assuming sources\n",
                "and receivers are at the surface Z=0) is %g s.\n",
                "You may want to adjust the \"minimum time of \n",
                "interest\" parameter."
            ),
            tminr
        );
        warn!("Total number of small reflecting\nsegments is %d.\n", ns);
    }

    // susynlv.c does not call checkpars().
    // Set constant segy trace header parameters.
    let mut tr = Trace::new();
    tr.trid = 1;
    tr.counit = 1;
    tr.ns = nt as u16;
    // C converts to an integer before narrowing to the 16-bit header field.
    tr.dt = (1.0e6 * dt as f64) as i32 as u16;
    tr.delrt = (1.0e3 * ft as f64) as i32 as i16;

    // Own the makeone static filter here, initializing it on the first trace.
    let mut hd = None;
    let mut tracl = 0i32;
    for ixsm in 0..nxsm {
        for ixo in 0..nxo {
            // Compute source and receiver coordinates.
            let xs = if shots {
                fxs + ixsm as f32 * dxs
            } else {
                ((fxm + ixsm as f32 * dxm) as f64 - 0.5 * xo[ixo as usize] as f64) as f32
            };
            let zs = 0.0;
            let xg = xs + xo[ixo as usize];
            let zg = 0.0;

            // Set segy trace header parameters.
            tracl = tracl.wrapping_add(1);
            tr.tracl = tracl;
            tr.tracr = tracl;
            if shots {
                tr.fldr = 1i32.wrapping_add(ixsm);
                tr.tracf = 1i32.wrapping_add(ixo);
                tr.d2 = dxo;
                tr.f2 = fxo;
            } else {
                tr.cdp = 1i32.wrapping_add(ixsm);
                tr.cdpt = 1i32.wrapping_add(ixo);
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

            // Make and write one trace.
            makeone(
                v00,
                dvdx,
                dvdz,
                ls,
                er,
                ob,
                &w,
                xs,
                zs,
                xg,
                zg,
                &r,
                nt,
                dt,
                ft,
                &mut tr.data[..nt as usize],
                &mut hd,
            );
            puttr(&mut tr);
        }
    }
    exit(CWP_Exit());
}

/// Make one synthetic seismogram for v(x,z) = v00+dvdx*x+dvdz*z.
#[allow(clippy::too_many_arguments)]
fn makeone(
    v00: f32,
    dvdx: f32,
    dvdz: f32,
    ls: i32,
    er: i32,
    ob: i32,
    w: &Wavelet,
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
    // If the half-derivative filter is not yet made, make it.
    let hd = hd.get_or_insert_with(|| {
        let mut hd = [0.0; NHD];
        mkhdiff(dt, LHD, &mut hd);
        hd
    });

    for it in 0..nt as usize {
        trace[it] = 0.0;
    }

    // Velocities at source and receiver.
    let vs = v00 + dvdx * xs + dvdz * zs;
    let vg = v00 + dvdx * xg + dvdz * zg;

    for ir in 0..r.len() {
        // Amplitude, number of segments, segment length.
        let ar = r[ir].a;
        let ns = r[ir].ns;
        let ds = r[ir].ds;
        let rs = &r[ir].rs;

        for is in 0..ns as usize {
            // Diffractor midpoint, unit-normal, and length.
            let xd = rs[is].x;
            let zd = rs[is].z;
            let cd = rs[is].c;
            let sd = rs[is].s;
            let vd = v00 + dvdx * xd + dvdz * zd;

            // Rays from shot and receiver to diffractor.
            let (cs, ss, ts, qs) = raylv2(v00, dvdx, dvdz, xs, zs, xd, zd);
            let (cg, sg, tg, qg) = raylv2(v00, dvdx, dvdz, xg, zg, xd, zd);

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
            let mut amp = if er != 0 {
                ((vg * vd / qg) as f64).sqrt() as f32
            } else if ls != 0 {
                ((vs * vd * vd * vg / (qs * qg)) as f64).sqrt() as f32
            } else {
                ((vs * vd * vd * vg / (qs * qg * (qs + qg))) as f64).sqrt() as f32
            };
            amp *= (ci + cr) * ar * ds;

            // Add sinc wavelet to trace.
            addsinc(time, amp, dt, ft, trace);
        }
    }

    // Apply the half-derivative filter, then convolve the wavelet with trace.
    let mut temp = vec![0.0; nt as usize];
    convolve_cwp(-LHD, hd, 0, trace, 0, &mut temp);
    convolve_cwp(w.iw, &w.wv[..w.lw as usize], 0, &temp, 0, trace);
}
