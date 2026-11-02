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

//! SUDMOVZ - DMO for V(Z) media for common-offset gathers.
//! Port of su/main/dip_moveout/sudmovz.c (Revision 1.6).
//! Credits: CWP: Craig Artley. Ray tracing: Dave Hale; modified by Craig Artley.
//! Reference: Artley, C., 1992, Dip-moveout processing for depth-variable
//! velocity, M.S. Thesis, Colorado School of Mines; CWP-115.

#![allow(clippy::too_many_arguments)]
use oxwave_core::prelude::*;
use oxwave_cwp::{
    complex::*,
    cwp::{abs, nint, PI},
};
use oxwave_inverse_interp::*;
use oxwave_linear_interp::*;
use oxwave_lu::*;
use oxwave_pfa_fft::*;
use oxwave_sinc_interp::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUDMOVZ - DMO for V(Z) media for common-offset gathers\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sudmovz <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [...]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " cdpmin         minimum cdp (integer number) for which to apply DMO\t",
    " cdpmax         maximum cdp (integer number) for which to apply DMO\t",
    " dxcdp          distance between adjacent cdp bins (m)\t\t\t",
    " noffmix        number of offsets to mix (see notes)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " vfile=         binary (non-ascii) file containing interval velocities (m/s)",
    " tdmo=0.0       times corresponding to interval velocities in vdmo (s)\t",
    " vdmo=1500.0    interval velocities corresponding to times in tdmo (m/s)",
    " fmax=0.5/dt    maximum frequency in input traces (Hz)\t\t\t",
    " smute=1.5      stretch mute used for NMO correction\t\t\t",
    " speed=1.0      twist this knob for speed (and aliasing)\t\t",
    " verbose=0      =1 for diagnostic print\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input traces should be sorted into common-offset gathers.  One common-",
    " offset gather ends and another begins when the offset field of the trace",
    " headers changes.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
    " the cdp location expressed in units of meters or feet.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The number of offsets to mix (noffmix) should typically equal the ratio of",
    " the shotpoint spacing to the cdp spacing.  This choice ensures that every",
    " cdp will be represented in each offset mix.  Traces in each mix will\t",
    " contribute through DMO to other traces in adjacent cdps within that mix.",
    "\t\t\t\t\t\t\t\t\t",
    " vfile should contain the regularly sampled interval velocities as a\t",
    " function of time.  If vfile is not supplied, the interval velocity\t",
    " function is defined by linear interpolation of the values in the tdmo\t",
    " and vdmo arrays.  The times in tdmo must be monotonically increasing.\t",
    "\t\t\t\t\t\t\t\t\t",
    " For each offset, the minimum time to process is determined using the\t",
    " smute parameter.  The DMO correction is not computed for samples that\t",
    " have experienced greater stretch during NMO.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed:  nt, dt, delrt, offset, cdp.\t\t",
];

fn dmoapply(
    nx: usize,
    dx: f32,
    nt: usize,
    dt: f32,
    np: usize,
    dp2: f32,
    h: f32,
    tntable: &[f32],
    it0min: &[usize],
    nt0: &[usize],
    tmin: f32,
    fmax: f32,
    qx: &mut [f32],
) {
    let nxfft = npfar(nx as i32) as usize;
    let nk = nxfft / 2 + 1;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let fk = 0.0;
    let ft = 0.0;
    let mut qk = vec![Complex::default(); nk * nt];
    pfa2rc(-1, 2, nt as i32, nxfft as i32, qx, &mut qk);
    let mut k = fk;
    for ik in 0..nk {
        jakub1k(
            k,
            h,
            np,
            dp2,
            nt,
            dt,
            ft,
            tntable,
            it0min,
            nt0,
            fmax,
            &mut qk[ik * nt..(ik + 1) * nt],
        );
        k += dk;
    }
    pfa2cr(1, 2, nt as i32, nxfft as i32, &qk, qx);
    let scale = (1.0 / nxfft as f64) as f32;
    let itmin = (nint((tmin / dt) as f64) as usize).min(nt);
    for ix in 0..nx {
        for it in 0..itmin {
            qx[ix * nt + it] = 0.0;
        }
        for it in itmin..nt {
            qx[ix * nt + it] *= scale;
        }
    }
}

fn jakub1k(
    k: f32,
    h: f32,
    np: usize,
    dp2: f32,
    nt: usize,
    dt: f32,
    ft: f32,
    tntable: &[f32],
    it0min: &[usize],
    nt0: &[usize],
    fmax: f32,
    qq: &mut [Complex],
) {
    let ntfft = npfa(nt as i32) as usize;
    let nw = ntfft;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = (-PI / dt as f64 + 0.0 * h as f64) as f32;
    let czero = cmplx(0.0, 0.0);
    let mut qn = vec![czero; ntfft];
    let mut q0 = vec![czero; ntfft];
    let mut p = 0.0f32;
    let mut pl = 0.0f32;
    let mut ph2 = dp2;
    for ip in 0..np {
        let ph = (ph2 as f64).sqrt() as f32;
        let wl = k / ph;
        let wh = if pl == 0.0 {
            (2.0 * PI * fmax as f64) as f32
        } else {
            let a = (k / pl) as f64;
            let b = 2.0 * PI * fmax as f64;
            (if a < b { a } else { b }) as f32
        };
        let iwlm = ((-wl - fw) / dw) as i32;
        let mut iwhm = ((-wh - fw) / dw) as i32;
        if iwhm < 0 {
            iwhm = 0;
        }
        iwhm += 1;
        let mut iwlp = ((wl - fw) / dw) as i32;
        let mut iwhp = ((wh - fw) / dw) as i32;
        if iwhp >= nw as i32 {
            iwhp = nw as i32 - 1;
        }
        iwlp += 1;
        if iwhm <= iwlm || iwlp <= iwhp {
            for it in 0..it0min[ip] {
                qn[it] = czero;
            }
            // Repeated extrapolation of one-sample mappings can exceed
            // intt8c's integer coordinate range. C then indexes outside its
            // sinc table; reject that undefined input before calling it.
            let xouts = (1.0 / dt as f64) as f32;
            let xoutb = (8.0 - (ft * xouts) as f64) as f32;
            for &tn in &tntable[ip * nt..ip * nt + nt0[ip]] {
                let xoutn = xoutb + tn * xouts;
                if !xoutn.is_finite() || xoutn < i32::MIN as f32 || xoutn >= -(i32::MIN as f32) {
                    err!("DMO mapping exceeds sinc interpolation sample range");
                }
            }
            ints8c(
                dt,
                ft,
                qq,
                czero,
                czero,
                &tntable[ip * nt..ip * nt + nt0[ip]],
                &mut qn[it0min[ip]..it0min[ip] + nt0[ip]],
            );
            for it in it0min[ip] + nt0[ip]..ntfft {
                qn[it] = czero;
            }
            for it in (1..ntfft).step_by(2) {
                qn[it].r = -qn[it].r;
                qn[it].i = -qn[it].i;
            }
            pfacc(1, &mut qn);
            for iw in iwhm..=iwlm {
                q0[iw as usize].r += qn[iw as usize].r;
                q0[iw as usize].i += qn[iw as usize].i;
            }
            for iw in iwlp..=iwhp {
                q0[iw as usize].r += qn[iw as usize].r;
                q0[iw as usize].i += qn[iw as usize].i;
            }
        }
        pl = p;
        p = ph;
        ph2 += dp2;
    }
    pfacc(-1, &mut q0);
    let scale = (0.5 / ntfft as f64) as f32;
    for it in (0..nt).step_by(2) {
        qq[it].r = q0[it].r * scale;
        qq[it].i = q0[it].i * scale;
    }
    for it in (1..nt).step_by(2) {
        qq[it].r = -q0[it].r * scale;
        qq[it].i = -q0[it].i * scale;
    }
}

fn tmute(vrms: &[f32], nt: usize, dt: f32, ft: f32, h: f32, smute: f32) -> f32 {
    let mut t = vec![0.0f32; nt];
    let mut tn = ft;
    for it in 0..nt {
        t[it] = ((tn * tn) as f64 + 4.0 * h as f64 * h as f64 / (vrms[it] * vrms[it]) as f64).sqrt()
            as f32;
        tn += dt;
    }
    let mut it = nt - 1;
    while it > 0 {
        if dt / (t[it] - t[it - 1]) > smute {
            break;
        }
        it -= 1;
    }
    ft + it as f32 * dt
}

fn get_x(tsg: f32, h: f32, p0: f32, v: f32, vars: &mut [f32; 5]) {
    let h2 = h * h;
    let a2 = ((v * v * tsg * tsg) as f64 / 4.0) as f32;
    let d = a2 - h2;
    let b2 = if d > 0.0 { d } else { 0.0 };
    let sine = (0.5 * p0 as f64 * v as f64) as f32;
    let d = 1.0 - (sine * sine) as f64;
    let tangent = (sine as f64 / (if d > 0.0 { d } else { 0.0 }).sqrt()) as f32;
    let x = (tangent as f64 * ((a2 / (b2 + tangent * tangent)) as f64).sqrt()) as f32;
    let d = 1.0 - (x * x / a2) as f64;
    let z = (b2 as f64 * if d > 0.0 { d } else { 0.0 }).sqrt() as f32;
    let x0 = h2 / a2 * x;
    let tg = ((((h - x) * (h - x) + z * z) as f64).sqrt() / v as f64) as f32;
    let t0 = ((((x0 - x) * (x0 - x) + z * z) as f64).sqrt() * 2.0 / v as f64) as f32;
    let ps = (x + h) / (v * v * (tsg - tg));
    let pg = (x - h) / (v * v * tg);
    vars[0] = x0;
    vars[1] = ps;
    vars[2] = pg;
    vars[3] = tg;
    vars[4] = t0;
}

fn get_fj(
    x: &[Vec<f32>],
    a: &[Vec<f32>],
    tau: &[Vec<f32>],
    nt: usize,
    np: usize,
    dt: f32,
    dp: f32,
    tsg: f32,
    h: f32,
    p0: f32,
    vars: &[f32; 5],
    eqns: &mut [f32; 5],
    jacob: &mut [[f32; 5]; 5],
) {
    let x0 = vars[0];
    let ps = vars[1];
    let pg = vars[2];
    let tg = vars[3];
    let t0 = vars[4];
    let ts = if tsg - tg > 0.0 { tsg - tg } else { 0.0 };
    let ft = 0.0;
    let fp = 0.0;
    for i in 0..5 {
        eqns[i] = 0.0;
    }
    for j in 0..5 {
        for i in 0..5 {
            jacob[j][i] = 0.0;
        }
    }
    eqns[0] = (eqns[0] as f64 + 2.0 * h as f64) as f32;
    eqns[1] += h - x0;
    jacob[0][1] = (jacob[0][1] as f64 - 1.0) as f32;

    let (ftp, dfdt, dfdp) = intder(x, nt, dt, ft, np, dp, fp, (2.0 * ts as f64) as f32, abs(ps));
    let sign = if ps < 0.0 { -1.0f32 } else { 1.0f32 };
    eqns[0] -= sign * ftp;
    jacob[1][0] -= dfdp;
    jacob[3][0] = (jacob[3][0] as f64 + 2.0 * sign as f64 * dfdt as f64) as f32;

    let (ftp, dfdt, dfdp) = intder(x, nt, dt, ft, np, dp, fp, (2.0 * tg as f64) as f32, abs(pg));
    let sign = if pg < 0.0 { -1.0f32 } else { 1.0f32 };
    eqns[0] += sign * ftp;
    eqns[1] += sign * ftp;
    jacob[2][0] += dfdp;
    jacob[2][1] += dfdp;
    jacob[3][0] = (jacob[3][0] as f64 + 2.0 * sign as f64 * dfdt as f64) as f32;
    jacob[3][1] = (jacob[3][1] as f64 + 2.0 * sign as f64 * dfdt as f64) as f32;

    let (ftp, dfdt, _) = intder(x, nt, dt, ft, np, dp, fp, t0, abs(p0));
    let sign = if p0 < 0.0 { -1.0f32 } else { 1.0f32 };
    eqns[1] -= sign * ftp;
    jacob[4][1] -= sign * dfdt;

    let (ftp, dfdt, dfdp) = intder(
        tau,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        (2.0 * ts as f64) as f32,
        abs(ps),
    );
    let sign = if ps < 0.0 { -1.0f32 } else { 1.0f32 };
    eqns[2] -= ftp;
    jacob[1][2] -= sign * dfdp;
    jacob[3][2] = (jacob[3][2] as f64 + 2.0 * dfdt as f64) as f32;

    let (ftp, dfdt, dfdp) = intder(
        tau,
        nt,
        dt,
        ft,
        np,
        dp,
        fp,
        (2.0 * tg as f64) as f32,
        abs(pg),
    );
    let sign = if pg < 0.0 { -1.0f32 } else { 1.0f32 };
    eqns[2] += ftp;
    eqns[3] += ftp;
    jacob[2][2] += sign * dfdp;
    jacob[2][3] += sign * dfdp;
    jacob[3][2] = (jacob[3][2] as f64 + 2.0 * dfdt as f64) as f32;
    jacob[3][3] = (jacob[3][3] as f64 + 2.0 * dfdt as f64) as f32;

    let (ftp, dfdt, _) = intder(tau, nt, dt, ft, np, dp, fp, t0, abs(p0));
    eqns[3] -= ftp;
    jacob[4][3] -= dfdt;

    let (ftp, dfdt, dfdp) = intder(a, nt, dt, ft, np, dp, fp, (2.0 * ts as f64) as f32, abs(ps));
    let sign = if ps < 0.0 { -1.0f32 } else { 1.0f32 };
    eqns[4] += sign * ftp;
    jacob[1][4] += dfdp;
    jacob[3][4] = (jacob[3][4] as f64 - 2.0 * sign as f64 * dfdt as f64) as f32;

    let (ftp, dfdt, dfdp) = intder(a, nt, dt, ft, np, dp, fp, (2.0 * tg as f64) as f32, abs(pg));
    let sign = if pg < 0.0 { -1.0f32 } else { 1.0f32 };
    eqns[4] += sign * ftp;
    jacob[2][4] += dfdp;
    jacob[3][4] = (jacob[3][4] as f64 + 2.0 * sign as f64 * dfdt as f64) as f32;

    let (ftp, dfdt, _) = intder(a, nt, dt, ft, np, dp, fp, t0, abs(p0));
    let sign = if p0 < 0.0 { -1.0f32 } else { 1.0f32 };
    eqns[4] = (eqns[4] as f64 - 2.0 * sign as f64 * ftp as f64) as f32;
    jacob[4][4] = (jacob[4][4] as f64 - 2.0 * sign as f64 * dfdt as f64) as f32;
}

fn intder(
    f: &[Vec<f32>],
    nx: usize,
    dx: f32,
    fx: f32,
    ny: usize,
    dy: f32,
    fy: f32,
    x: f32,
    y: f32,
) -> (f32, f32, f32) {
    let lx = fx + (nx - 1) as f32 * dx;
    let ly = fy + (ny - 1) as f32 * dy;
    let (ix0, ix1, fracx) = if x <= fx {
        (0, 0, 0.0)
    } else if x >= lx {
        (nx - 1, nx - 1, 1.0)
    } else {
        let x0i = (x - fx) / dx;
        let ix0 = x0i as usize;
        (ix0, ix0 + 1, x0i - ix0 as f32)
    };
    let (iy0, iy1, fracy) = if y <= fy {
        (0, 0, 0.0)
    } else if y >= ly {
        (ny - 1, ny - 1, 1.0)
    } else {
        let y0i = (y - fy) / dy;
        let iy0 = y0i as usize;
        (iy0, iy0 + 1, y0i - iy0 as f32)
    };
    let f00 = f[iy0][ix0];
    let f01 = f[iy1][ix0];
    let f10 = f[iy0][ix1];
    let f11 = f[iy1][ix1];
    let fxy = (1.0f32 - fracx) * (1.0f32 - fracy) * f00
        + (1.0f32 - fracx) * fracy * f01
        + fracx * (1.0f32 - fracy) * f10
        + fracx * fracy * f11;
    let dfdx = ((1.0f32 - fracy) * (f10 - f00) + fracy * (f11 - f01)) / dx;
    let dfdy = ((1.0f32 - fracx) * (f01 - f00) + fracx * (f11 - f10)) / dy;
    (fxy, dfdx, dfdy)
}

fn rayvt(p: f32, nt: usize, dt: f32, v: &[f32], tau: &mut [f32], x: &mut [f32], a: &mut [f32]) {
    let mut vt = v[0];
    let mut dvt = v[1] - v[0];
    let pvt = p * vt;
    let mut at = (if 1.0 < pvt as f64 { 1.0 } else { pvt as f64 }).asin() as f32;
    tau[0] = 0.0;
    x[0] = 0.0;
    a[0] = at;
    let mut it = 1;
    while it < nt {
        let cosat = (at as f64).cos() as f32;
        tau[it] = tau[it - 1] + dt * cosat;
        if tau[it] <= 0.0 {
            break;
        }
        x[it] = x[it - 1] + dt * p * vt * vt;
        at += p * dvt;
        a[it] = at;
        if it < nt - 1 {
            let taui = tau[it] / dt;
            let itau = taui as usize;
            let frac = taui - itau as f32;
            vt = ((1.0 - frac as f64) * v[itau] as f64 + (frac * v[itau + 1]) as f64) as f32;
            dvt = v[itau + 1] - v[itau];
        }
        it += 1;
    }
    for jt in it..nt {
        tau[jt] = tau[jt - 1];
        x[jt] = x[jt - 1];
        a[jt] = a[jt - 1];
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let mut ft = (tr.delrt as f64 / 1000.0) as f32;
    if ft != 0.0 {
        // The reference executable prints zero for this argument for every
        // nonzero delrt (including negative and maximum i16 delays).
        err!("ft=%f --- nonzero ft not yet implemented!", 0.0);
    }
    let mut offset = tr.offset as f32;
    let (mut cdpmin, mut cdpmax, mut noffmix) = (0, 0, 0);
    let (mut dx, mut fmax, mut speed, mut smute) = (0.0f32, 0.0f32, 1.0f32, 1.5f32);
    let mut verbose = 0;
    if getparint("cdpmin", &mut cdpmin) == 0 {
        err!("must specify cdpmin");
    }
    if getparint("cdpmax", &mut cdpmax) == 0 {
        err!("must specify cdpmax");
    }
    if cdpmin > cdpmax {
        err!("cdpmin must not be greater than cdpmax");
    }
    if getparfloat("dxcdp", &mut dx) == 0 {
        err!("must specify dxcdp");
    }
    if getparint("noffmix", &mut noffmix) == 0 {
        err!("must specify noffmix");
    }
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (0.5 / dt as f64) as f32;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    if getparfloat("speed", &mut speed) == 0 {
        speed = 1.0;
    }
    if getparfloat("smute", &mut smute) == 0 {
        smute = 1.5;
    }
    let nttab = 2 * nt - 1;
    // rayvt reads v[1] even when nt=1; reserve a zero guard for that read.
    let mut vt = vec![0.0f32; nttab.max(2)];
    let mut vrms = vec![0.0f32; nttab];
    let mut vfile = String::new();
    if getparstring("vfile", &mut vfile) == 0 {
        let ntdmo = countparval("tdmo").max(1) as usize;
        let mut tdmo = vec![0.0f32; ntdmo];
        if getparfloat("tdmo", &mut tdmo) == 0 {
            tdmo[0] = 0.0;
        }
        let nvdmo = countparval("vdmo").max(1) as usize;
        if nvdmo != ntdmo {
            err!("number of tdmo and vdmo must be equal");
        }
        let mut vdmo = vec![0.0f32; nvdmo];
        if getparfloat("vdmo", &mut vdmo) == 0 {
            vdmo[0] = 1500.0;
        }
        for itdmo in 1..ntdmo {
            if tdmo[itdmo] <= tdmo[itdmo - 1] {
                err!("tdmo must increase monotonically");
            }
        }
        let mut t = 0.0f32;
        for it in 0..nt {
            intlin(
                &tdmo,
                &vdmo,
                vdmo[0],
                vdmo[ntdmo - 1],
                &[t],
                &mut vt[it..it + 1],
            );
            t += dt;
        }
    } else {
        let fp = fopen(&vfile, "r");
        if fp.as_ref().map_or(0, |fp| fread(&mut vt[..nt], fp)) != nt {
            err!("cannot read %d velocities from file %s", nt, vfile);
        }
    }
    checkpars();
    for it in nt..nttab {
        vt[it] = vt[nt - 1];
    }
    vrms[0] = vt[0];
    let mut sum = 0.0f32;
    for it in 1..nttab {
        sum += vt[it - 1] * vt[it - 1];
        vrms[it] = ((sum / it as f32) as f64).sqrt() as f32;
    }
    for it in 0..nttab {
        vt[it] = (vt[it] as f64 * 0.5) as f32;
    }
    ft = 0.0;
    let fp = 0.0;
    let lp = (1.0 / vt[0] as f64) as f32;
    let nptab = 401usize;
    let dptab = lp / (nptab - 1) as f32;
    let mut tau = vec![vec![0.0f32; nttab]; nptab];
    let mut x = tau.clone();
    let mut a = tau.clone();
    if verbose != 0 {
        fprintf!(stderr(), "Computing ray tables...\n");
    }
    let mut p = fp;
    for ip in 0..nptab {
        rayvt(p, nttab, dt, &vt, &mut tau[ip], &mut x[ip], &mut a[ip]);
        p += dptab;
    }
    if verbose != 0 {
        fprintf!(stderr(), "Ray tables ready...\n");
    }
    drop(vt);
    let nx = (cdpmax - cdpmin + 1) as usize;
    let nxpad = (0.5 * abs(offset / dx) as f64) as i32;
    let nxfft = npfar(nx as i32 + nxpad) as usize;
    let mut q = vec![0.0f32; nt * (nxfft + 2)];
    let mut qmix = vec![0.0f32; nt * nx];
    let hfp = etmpfile();
    let mut ntrace = 0;
    let mut tro = Trace::new();
    let mut oldoffset = offset;
    let mut gottrace = true;
    let mut noff = 0;
    loop {
        if gottrace {
            offset = tr.offset as f32;
        }
        if (gottrace && offset != oldoffset) || !gottrace {
            dmooff(
                &x, &a, &tau, nttab, dt, ft, nptab, dptab, fp, lp, nx, dx, oldoffset, nt, &vrms,
                fmax, smute, speed, verbose, &mut q,
            );
            for ix in 0..nx {
                for it in 0..nt {
                    qmix[ix * nt + it] += q[ix * nt + it];
                }
            }
            noff += 1;
            if offset != oldoffset {
                let nxpad = (0.5 * abs(offset / dx) as f64) as i32;
                let nxfft = npfar(nx as i32 + nxpad) as usize;
                q = vec![0.0f32; nt * (nxfft + 2)];
            }
        }
        if noff == noffmix || !gottrace {
            efseeko(&hfp, 0, SEEK_SET);
            for _itrace in 0..ntrace {
                efread(std::slice::from_mut(&mut tro.hdr), &hfp);
                let ix = (tro.cdp - cdpmin) as usize;
                tro.data[..nt].copy_from_slice(&qmix[ix * nt..(ix + 1) * nt]);
                puttr(&mut tro);
            }
            if verbose != 0 {
                fprintf!(
                    stderr(),
                    "\tCompleted mix of %d offsets with %d traces\n",
                    noff,
                    ntrace
                );
            }
            if !gottrace {
                break;
            }
            efseeko(&hfp, 0, SEEK_SET);
            noff = 0;
            ntrace = 0;
            qmix.fill(0.0);
        }
        if tr.cdp >= cdpmin && tr.cdp <= cdpmax {
            efwrite(std::slice::from_ref(&tr.hdr), &hfp);
            ntrace += 1;
            oldoffset = offset;
            let ix = (tr.cdp - cdpmin) as usize;
            q[ix * nt..(ix + 1) * nt].copy_from_slice(&tr.data[..nt]);
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    exit(CWP_Exit());
}

fn dmooff(
    x: &[Vec<f32>],
    a: &[Vec<f32>],
    tau: &[Vec<f32>],
    nttab: usize,
    dt: f32,
    ft: f32,
    nptab: usize,
    dptab: f32,
    fp: f32,
    lp: f32,
    nx: usize,
    dx: f32,
    offset: f32,
    nt: usize,
    vrms: &[f32],
    fmax: f32,
    smute: f32,
    speed: f32,
    verbose: i32,
    q: &mut [f32],
) {
    let tmin = tmute(vrms, nt, dt, ft, (0.5 * offset as f64) as f32, smute);
    let lt = (nt - 1) as f32 * dt;
    if tmin >= lt {
        err!("tmin=%f >= lt=%f:  no data survive the mute!\n", tmin, lt);
    }
    let dtn = (10.0 * dt as f64) as f32;
    let ftn = tmin;
    let ntn = (nint(((lt - ftn) / dtn) as f64) + 1) as usize;
    let mut dp2 = (tmin as f64 / (fmax as f64 * 0.25 * offset as f64 * offset as f64)) as f32;
    dp2 *= speed;
    let np = nint((lp * lp / dp2) as f64).wrapping_add(1);
    // A negative dimension becomes a huge size_t in C's ealloc1int.
    if np < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let np = np as usize;
    let mut it0min = vec![0usize; np];
    let mut nt0 = vec![0usize; np];
    // C reads one past a row when no mapping was found. The tables are
    // contiguous and zeroed, so supply that next row's zero explicitly.
    let mut t0table = vec![vec![0.0f32; ntn + 1]; np];
    // ealloc2float has contiguous rows. A nonmonotone mapping can make
    // the initial nt0 negative, so extrapolation starts in a previous row.
    let mut tntable = vec![0.0f32; nt * np];
    let h = (offset as f64 / 2.0).abs() as f32;
    let mut tn = ftn;
    for itn in 0..ntn {
        dmomap(
            x,
            a,
            tau,
            nttab,
            dt,
            ft,
            nptab,
            dptab,
            np,
            dp2,
            fp,
            lp,
            itn,
            tn,
            h,
            vrms,
            verbose,
            &mut t0table,
        );
        tn += dtn;
    }
    let mut ip = 0;
    while ip < np {
        let mut itn = 0;
        while itn < ntn && !(t0table[ip][itn] > 0.0) {
            itn += 1;
        }
        let itnmin = itn;
        it0min[ip] = (nint((t0table[ip][itnmin] / dt) as f64) as usize).min(nt - 1);
        let mut itnmax = ntn as i32 - 1;
        while itnmax >= 0 && !(t0table[ip][itnmax as usize] > 0.0) {
            itnmax -= 1;
        }
        if itnmax <= itnmin as i32 {
            break;
        }
        let itnmax = itnmax as usize;
        let it0max = (nint((t0table[ip][itnmax] / dt) as f64) as usize).min(nt - 1);
        let n = it0max as i32 - it0min[ip] as i32 + 1;
        for itn in itnmin + 1..itnmax {
            if t0table[ip][itn] <= t0table[ip][itn - 1] && verbose != 0 {
                fprintf!(
                    stderr(),
                    "t0table not monotonically increasing---extrapolating\n"
                );
            }
        }
        yxtoxy(
            dtn,
            ftn + itnmin as f32 * dtn,
            &t0table[ip][itnmin..itnmax + 1],
            dt,
            it0min[ip] as f32 * dt,
            ftn + itnmin as f32 * dtn,
            ftn + itnmax as f32 * dtn,
            // With ny <= 0, C's yxtoxy does not write any samples.
            &mut tntable[ip * nt..ip * nt + n.max(0) as usize],
        );
        for it in n..nt as i32 {
            let index = (ip * nt) as i64 + it as i64;
            // Reads within preceding rows are valid in the contiguous
            // allocation; reads before its beginning are undefined in C.
            if index < 2 {
                err!("DMO mapping has fewer than two samples for extrapolation");
            }
            let index = index as usize;
            tntable[index] = (2.0 * tntable[index - 1] as f64 - tntable[index - 2] as f64) as f32;
        }
        nt0[ip] = nt - it0min[ip];
        ip += 1;
    }
    let npdmo = ip.min(np);
    dmoapply(
        nx, dx, nt, dt, npdmo, dp2, h, &tntable, &it0min, &nt0, tmin, fmax, q,
    );
}

fn dmomap(
    x: &[Vec<f32>],
    a: &[Vec<f32>],
    tau: &[Vec<f32>],
    nt: usize,
    dt: f32,
    ft: f32,
    nptab: usize,
    dptab: f32,
    np: usize,
    dp2: f32,
    fp: f32,
    lp: f32,
    itn: usize,
    tn: f32,
    h: f32,
    vrms: &[f32],
    verbose: i32,
    t0table: &mut [Vec<f32>],
) {
    let lt = ft + (nt - 1) as f32 * dt + (0 * verbose) as f32;
    let mut vars = [0.0f32; 5];
    let mut eqns = [0.0f32; 5];
    let mut jacob = [[0.0f32; 5]; 5];
    let mut work = [0.0f32; 5];
    let mut ipvt = [0i32; 5];
    let vtn = vrms[(tn / dt) as usize];
    let tsg = ((tn * tn + 4.0f32 * h * h / (vtn * vtn)) as f64).sqrt() as f32;
    get_x(tsg, h, fp, vtn, &mut vars);
    let mut p2 = fp;
    // C leaves rcond uninitialized if its first trial already converges.
    // Any nonzero value has the same effect as the subsequent LU estimates.
    let mut rcond = 1.0f32;
    for ip in 0..np {
        let p0 = (p2 as f64).sqrt() as f32;
        for _iter in 0..10 {
            get_fj(
                x, a, tau, nt, nptab, dt, dptab, tsg, h, p0, &vars, &mut eqns, &mut jacob,
            );
            let mut error = 0.0f32;
            for i in 0..5 {
                error += abs(eqns[i]);
            }
            if error <= 1.0e-4f32 {
                break;
            }
            {
                let mut cols: Vec<&mut [f32]> = jacob.iter_mut().map(|r| &mut r[..]).collect();
                sgeco(&mut cols, 5, &mut ipvt, &mut rcond, &mut work);
            }
            let cols: Vec<&[f32]> = jacob.iter().map(|r| &r[..]).collect();
            sgesl(&cols, 5, &ipvt, &mut eqns, 0);
            if rcond != 0.0 {
                for i in 0..5 {
                    vars[i] -= eqns[i];
                }
            } else {
                break;
            }
            if abs(vars[1]) > lp {
                vars[1] = if vars[1] < 0.0 { -lp } else { lp };
            }
            if abs(vars[2]) > lp {
                vars[2] = if vars[2] < 0.0 { -lp } else { lp };
            }
            if vars[3] < 0.0 {
                vars[3] = 0.0;
            }
            if vars[3] > tsg {
                vars[3] = tsg;
            }
            if vars[4] < 0.0 {
                vars[4] = 0.0;
            }
            if vars[4] > lt {
                vars[4] = lt;
            }
        }
        let (tautp, _, _) = intder(tau, nt, dt, ft, nptab, dptab, fp, vars[4], abs(p0));
        if tautp < dt {
            break;
        }
        if rcond == 0.0 {
            break;
        }
        let x0 = vars[0];
        let t0 = vars[4];
        t0table[ip][itn] = t0 + p0 * x0;
        p2 += dp2;
    }
}
