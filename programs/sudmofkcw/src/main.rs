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

//! SUDMOFKCW - converted-wave DMO via F-K domain (log-stretch) method.
//! Port of su/main/dip_moveout/sudmofkcw.c (Revision 1.8).
//! Credits: CWP: Mohamed Alfaraj, Dave Hale.
//! Technical references: Mohammed Alfaraj, Ph.D. thesis (1993), Colorado
//! School of Mines; Dave Hale, Dip-Moveout Processing, SEG Course Notes (1988).

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::complex::cmplx;
use oxwave_cwp::cwp::{abs, max, min, nint, PI};
use oxwave_f_cubicspline::cmonot;
use oxwave_inverse_interp::yxtoxy;
use oxwave_pfa_fft::{npfa, npfar, pfa2cr, pfa2rc, pfacc};
use oxwave_sinc_interp::ints8c;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUDMOFKCW - converted-wave DMO via F-K domain (log-stretch) method for",
    " \t\tcommon-offset gathers\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sudmofkcw <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [...]\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " cdpmin\t\t  minimum cdp (integer number) for which to apply DMO",
    " cdpmax\t\t  maximum cdp (integer number) for which to apply DMO",
    " dxcdp\t\t   distance between adjacent cdp bins (m)\t\t",
    " noffmix\t\t number of offsets to mix (see notes)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " tdmo=0.0\t\ttimes corresponding to rms velocities in vdmo (s)",
    " vdmo=1500.0\t\trms velocities corresponding to times in tdmo (m/s)",
    " gamma=0.5\t\t velocity ratio, upgoing/downgoing\t\t",
    " ntable=1000\t\t number of tabulated z/h and b/h (see notes)\t",
    " sdmo=1.0\t\tDMO stretch factor; try 0.6 for typical v(z)\t",
    " flip=0\t\t =1 for negative shifts and exchanging s1 and s2",
    " \t\t\t (see notes below)\t\t\t\t",
    " fmax=0.5/dt\t\tmaximum frequency in input traces (Hz)\t\t",
    " verbose=0\t\t=1 for diagnostic print\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input traces should be sorted into common-offset gathers.  One common-",
    " offset gather ends and another begins when the offset field of the trace",
    " headers changes.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
    " the cdp location expressed in units of meters or feet.\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The number of offsets to mix (noffmix) should typically equal the ratio of",
    " the shotpoint spacing to the cdp spacing.  This choice ensures that every",
    " cdp will be represented in each offset mix.  Traces in each mix will\t",
    " contribute through DMO to other traces in adjacent cdps within that mix.",
    " \t\t\t\t\t\t\t\t\t",
    " The tdmo and vdmo arrays specify a velocity function of time that is\t",
    " used to implement a first-order correction for depth-variable velocity.",
    " The times in tdmo must be monotonically increasing. The velocity function",
    " is assumed to have been gotten by traditional NMO. \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " For each offset, the minimum time at which a non-zero sample exists is",
    " used to determine a mute time.  Output samples for times earlier than this",
    " mute time will be zeroed.  Computation time may be significantly reduced",
    " if the input traces are zeroed (muted) for early times at large offsets.",
    " \t\t\t\t\t\t\t\t\t",
    " z/h is horizontal-reflector depth normalized to half source-reciver offset",
    " h.  Normalized shift of conversion point is b/h.  The code now does not",
    " support signed offsets, therefore it is recommended that only end-on data,",
    " not split-spread, be used as input (of course after being sorted into\t",
    " common-offset gathers).  z/h vs b/h depends on gamma (see Alfaraj's Ph.D.",
    " thesis, 1993).\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Flip factor = 1 implies positive shift of traces (in the increasing CDP",
    " bin number direction).  When processing split-spread data, for example,",
    " if one side of the spread is processed with flip=0, then the other side",
    " of the spread should be processed with flip=1.  The flip factor also\t",
    " determines the actions of the factors s1 and s2, i.e., stretching or\t",
    " squeezing.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed:  nt, dt, delrt, offset, cdp.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get information from the first header */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut offset = tr.offset as f32;

    /* get parameters */
    let mut cdpmin = 0i32;
    let mut cdpmax = 0i32;
    let mut dx = 0.0f32;
    let mut noffmix = 0i32;
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
    let mut ntdmo = countparval("tdmo");
    if ntdmo == 0 {
        ntdmo = 1;
    }
    let mut tdmo = vec![0.0f32; ntdmo as usize];
    if getparfloat("tdmo", &mut tdmo) == 0 {
        tdmo[0] = 0.0;
    }
    let mut nvdmo = countparval("vdmo");
    if nvdmo == 0 {
        nvdmo = 1;
    }
    if nvdmo != ntdmo {
        err!("number of tdmo and vdmo must be equal");
    }
    let mut vdmo = vec![0.0f32; nvdmo as usize];
    if getparfloat("vdmo", &mut vdmo) == 0 {
        vdmo[0] = 1500.0;
    }
    for itdmo in 1..ntdmo as usize {
        if tdmo[itdmo] <= tdmo[itdmo - 1] {
            err!("tdmo must increase monotonically");
        }
    }
    let mut gamma = 0.0f32;
    let mut ntable = 0i32;
    let mut sdmo = 0.0f32;
    let mut flip = 0i32;
    let mut fmax = 0.0f32;
    let mut verbose = 0i32;
    if getparfloat("gamma", &mut gamma) == 0 {
        gamma = 0.5;
    }
    if getparint("ntable", &mut ntable) == 0 {
        ntable = 1000;
    }
    if getparfloat("sdmo", &mut sdmo) == 0 {
        sdmo = 1.0;
    }
    if getparint("flip", &mut flip) == 0 {
        flip = 0;
    }
    let sign = if flip != 0 { -1.0 } else { 1.0 };
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (0.5 / dt as f64) as f32;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    checkpars();

    /* allocate and generate tables of b/h and z/h if gamma not equal 1 */
    let mut zoh = Vec::new();
    let mut boh = Vec::new();
    if gamma != 1.0 {
        zoh = vec![0.0f32; ntable as usize];
        boh = vec![0.0f32; ntable as usize];
        table(ntable as usize, gamma, &mut zoh, &mut boh);
    }

    /* make uniformly sampled rms velocity function of time */
    let mut vrms = vec![0.0f32; nt];
    mkvrms(&tdmo, &vdmo, nt, dt, ft, &mut vrms);
    let nx = cdpmax.wrapping_sub(cdpmin).wrapping_add(1) as usize;
    let nxpad = (0.5 * abs(offset / dx) as f64) as i32;
    let nxfft = npfar((nx as i32).wrapping_add(nxpad)) as usize;
    // Flat storage preserves C's contiguous p[0] layout for the spatial FFT.
    let mut p = vec![0.0f32; nt * (nxfft + 2)];
    let mut q = vec![0.0f32; nt * nx];
    let hfp = etmpfile();

    let mut oldoffset = offset;
    let mut gottrace = true;
    let mut ntrace = 0;
    let mut noff = 0;
    let (mut s1, mut s2) = stretchfactor(sdmo, gamma);
    if flip != 0 {
        std::mem::swap(&mut s1, &mut s2);
    }
    if verbose != 0 {
        fprintf!(stderr(), "stretching factors: s1=%f s2=%f\n", s1, s2);
    }

    loop {
        if gottrace {
            offset = tr.offset as f32;
        }
        if (gottrace && offset != oldoffset) || !gottrace {
            dmooff(
                oldoffset, fmax, nx, dx, nt, dt, ft, &vrms, &mut p, gamma, &boh, &zoh, s1, s2, sign,
            );
            for ix in 0..nx {
                for it in 0..nt {
                    q[ix * nt + it] += p[ix * nt + it];
                }
            }
            noff += 1;
            p.clear();
            if offset != oldoffset {
                let nxpad = (0.5 * abs(offset / dx) as f64) as i32;
                let nxfft = npfar((nx as i32).wrapping_add(nxpad)) as usize;
                p = vec![0.0; nt * (nxfft + 2)];
            }
        }
        if noff == noffmix || !gottrace {
            efseeko(&hfp, 0, SEEK_SET);
            for _itrace in 0..ntrace {
                efread(tro.hdr.as_bytes_mut(), &hfp);
                let ix = (tro.cdp - cdpmin) as usize;
                tro.data[..nt].copy_from_slice(&q[ix * nt..(ix + 1) * nt]);
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
            for ix in 0..nx {
                for it in 0..nt {
                    q[ix * nt + it] = 0.0;
                }
            }
        }
        if tr.cdp >= cdpmin && tr.cdp <= cdpmax {
            efwrite(tr.hdr.as_bytes(), &hfp);
            ntrace += 1;
            oldoffset = offset;
            let ix = (tr.cdp - cdpmin) as usize;
            p[ix * nt..(ix + 1) * nt].copy_from_slice(&tr.data[..nt]);
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    exit(CWP_Exit());
}

/* make uniformly sampled vrms(t) for DMO; Dave Hale, 10/03/91. */
fn mkvrms(tdmo: &[f32], vdmo: &[f32], nt: usize, dt: f32, ft: f32, vrms: &mut [f32]) {
    let mut vdmod = vec![[0.0f32; 4]; tdmo.len()];
    cmonot(tdmo, vdmo, &mut vdmod);
    let mut t = ft;
    for it in 0..nt {
        intcub(0, tdmo, &vdmod, &[t], &mut vrms[it..it + 1]);
        t += dt;
    }
}

/* perform dmo in f-k domain for one offset; Dave Hale, 08/08/91. */
#[allow(clippy::too_many_arguments)]
fn dmooff(
    offset: f32,
    fmax: f32,
    nx: usize,
    dx: f32,
    nt: usize,
    dt: f32,
    ft: f32,
    vrms: &[f32],
    ptx: &mut [f32],
    gamma: f32,
    boh: &[f32],
    zoh: &[f32],
    s1: f32,
    s2: f32,
    sign: f32,
) {
    let nxfft = npfar((nx as i32).wrapping_add((0.5 * abs(offset / dx) as f64) as i32)) as usize;
    let mut itmin = nt;
    for ix in 0..nx {
        let mut it = 0;
        while it < itmin && ptx[ix * nt + it] == 0.0 {
            it += 1;
        }
        itmin = it;
    }
    if itmin >= nt {
        return;
    }
    // C's nt=1 path reads uoft[1] beyond its logical allocation. Its sole
    // temporal FFT coefficient is DC, which is explicitly zeroed below;
    // preserve the observed zero result without the out-of-bounds read.
    if nt == 1 {
        ptx[..nxfft].fill(0.0);
        return;
    }
    let (uoft, nu, du, fu, tofu, tcon) = maketu(offset, itmin, fmax, nt, dt, ft, vrms);

    let g1 = (2.0 * (gamma as f64).sqrt() / (1.0 + gamma as f64)) as f32;
    let h = (offset as f64 / 2.0) as f32;
    let lnt = nt - 1;
    let lt = (nt - 1) as f32 * dt;
    let vmin = (0.5
        * min(
            (1.0 + 1.0 / gamma as f64) * vrms[0] as f64,
            (1.0 + gamma as f64) * vrms[0] as f64,
        )) as f32;
    let (bboh, _itn) = if gamma != 1.0 {
        getbbohitn(offset, itmin, nt, dt, vrms, boh, zoh, gamma)
    } else {
        (Vec::new(), Vec::new())
    };
    let os1 = (1.0 / s1 as f64) as f32;
    let os2 = (1.0 / s2 as f64) as f32;
    let nupad = (1.5 * ((s1 + s2) as f64 / 2.0) * tcon as f64 / du as f64) as i32;
    let nufft = npfa((nu as i32).wrapping_add(nupad)) as usize;
    let nw = nufft;
    let dw = (2.0 * PI / (nufft as f32 * du) as f64) as f32;
    let mut pu = vec![cmplx(0.0, 0.0); nufft];
    let nk = nxfft / 2 + 1;
    let dk = (2.0 * PI / abs(nxfft as f32 * dx) as f64) as f32;
    let kmax = (PI / abs(dx) as f64) as f32;
    let ikmax = nint((kmax / dk) as f64) as usize;
    let mut ptk = vec![cmplx(0.0, 0.0); nk * nt];
    let fftscl = nk as f32 / (ikmax + 1) as f32 / (nufft as i32).wrapping_mul(nxfft as i32) as f32;
    pfa2rc(-1, 2, nt as i32, nxfft as i32, ptx, &mut ptk);

    let czero = cmplx(0.0, 0.0);
    let mut k = 0.0f32;
    for ik in 0..=ikmax {
        let hk = sign * h * k;
        let mut t = lt;
        for it in (itmin..=lnt).rev() {
            let phase = if gamma != 1.0 { bboh[it] * hk } else { 0.0 };
            let cphase = (phase as f64).cos() as f32;
            let sphase = (phase as f64).sin() as f32;
            let fr = ptk[ik * nt + it].r;
            let fi = ptk[ik * nt + it].i;
            ptk[ik * nt + it].r = fr * cphase + fi * sphase;
            ptk[ik * nt + it].i = -fr * sphase + fi * cphase;
            t -= dt;
        }
        let _ = t;
        ints8c(
            dt,
            ft,
            &ptk[ik * nt..(ik + 1) * nt],
            czero,
            czero,
            &tofu,
            &mut pu,
        );
        for iu in nu..nufft {
            pu[iu] = czero;
        }
        pfacc(1, &mut pu);
        let wmin = abs(0.5 * vmin as f64 * k as f64) as f32;
        let wmax = abs(PI / du as f64) as f32;
        let iwmin = max(1, min((nw as i32 - 1) / 2, nint((wmin / dw) as f64))) as usize;
        let iwmax = max(0, min((nw as i32 - 1) / 2, nint((wmax / dw) as f64))) as usize;
        let wwscl = (os1 as f64 * ((g1 * hk / tcon) as f64).powf(2.0)) as f32;
        let wwscl2 = wwscl * os2 / os1;
        pu[0] = czero;
        for iw in 1..iwmin {
            pu[iw] = czero;
            pu[nw - iw] = czero;
        }
        let mut w = iwmin as f32 * dw;
        for iw in iwmin..=iwmax {
            let iwn = nw - iw;
            let scales = (1.0 + (wwscl / (w * w)) as f64) as f32;
            let scale = (scales as f64).sqrt() as f32;
            let phase = ((s1 * w * tcon) as f64 * (scale as f64 - 1.0)) as f32;
            let amp = (fftscl as f64 * (1.0 - s1 as f64 + (s1 / scale) as f64)) as f32;
            let fr = (amp as f64 * (phase as f64).cos()) as f32;
            let fi = (amp as f64 * (phase as f64).sin()) as f32;
            let pwr = pu[iw].r;
            let pwi = pu[iw].i;
            pu[iw].r = pwr * fr - pwi * fi;
            pu[iw].i = pwr * fi + pwi * fr;

            let scales = (1.0 + (wwscl2 / (w * w)) as f64) as f32;
            let scale = (scales as f64).sqrt() as f32;
            let phase = ((s2 * w * tcon) as f64 * (scale as f64 - 1.0)) as f32;
            let amp = (fftscl as f64 * (1.0 - s2 as f64 + (s2 / scale) as f64)) as f32;
            let fr = (amp as f64 * (phase as f64).cos()) as f32;
            let fi = (amp as f64 * (phase as f64).sin()) as f32;
            let pwr = pu[iwn].r;
            let pwi = pu[iwn].i;
            pu[iwn].r = pwr * fr + pwi * fi;
            pu[iwn].i = pwi * fr - pwr * fi;
            w += dw;
        }
        for iw in iwmax + 1..=nw / 2 {
            pu[iw] = czero;
            pu[nw - iw] = czero;
        }
        pfacc(-1, &mut pu);
        ints8c(
            du,
            fu,
            &pu[..nu],
            czero,
            czero,
            &uoft,
            &mut ptk[ik * nt..(ik + 1) * nt],
        );
        k += dk;
    }
    for ik in ikmax + 1..nk {
        for it in 0..nt {
            ptk[ik * nt + it] = czero;
        }
    }
    pfa2cr(1, 2, nt as i32, nxfft as i32, &ptk, ptx);
}

/* make stretch and compress functions t(u) and u(t); Dave Hale, 08/08/91. */
fn maketu(
    offset: f32,
    itmin: usize,
    fmax: f32,
    nt: usize,
    dt: f32,
    ft: f32,
    vrms: &[f32],
) -> (Vec<f32>, usize, f32, f32, Vec<f32>, f32) {
    let numax = (500.0 + (nt as f32 as f64).ln() * (nt - 1) as f32 as f64) as i32;
    let mut uoft = vec![0.0f32; nt];
    let tmin = ft + itmin as f32 * dt;
    let et = ft + (nt - 1) as f32 * dt;
    let mut t1 = min(et, max(ft + dt, tmin));
    let mut t2 = if offset != 0.0 {
        max(
            t1 as f64,
            1.0 / (1.0 / et as f64
                + 0.2 * dt as f64 * vrms[0] as f64 * vrms[0] as f64 / (offset * offset) as f64),
        ) as f32
    } else {
        t1
    };
    t1 = ft + nint((t1 / dt) as f64) as f32 * dt;
    t2 = ft + nint((t2 / dt) as f64) as f32 * dt;

    let mut v2 = vrms[0];
    let mut v22 = v2 * v2;
    let mut v44 = v22 * v22;
    let mut gamma = 1.0f32;
    let mut t = ft;
    for it in 0..nt {
        let v2m = v2;
        let v22m = v22;
        let v44m = v44;
        let gammam = gamma;
        if t > 0.0 {
            v2 = vrms[it];
            v22 = v2 * v2;
            let vi2 = (t * v22 - (t - dt) * v22m) / dt;
            let vi4 = vi2 * vi2;
            v44 = (dt * vi4 + (t - dt) * v44m) / t;
        } else {
            v2 = v2m;
            v22 = v22m;
            v44 = v44m;
        }
        let dv2 = (v2 - v2m) / dt;
        let v24 = v22 * v22;
        gamma = (1.5 * v44 as f64 / v24 as f64 - (t * dv2 / v2) as f64 - 0.5) as f32;
        if t <= t1 {
            uoft[it] = t - t1;
        } else if t > t1 && t <= t2 {
            let du = (t1 as f64
                * (gamma as f64 * (t as f64 / (t as f64 - 0.5 * dt as f64)).ln()
                    - gammam as f64 * ((t - dt) as f64 / (t as f64 - 0.5 * dt as f64)).ln()))
                as f32;
            let dumin = (0.1 * dt as f64 * t1 as f64 / t as f64) as f32;
            // Match the bounds-checked C: negative delrt can start this
            // recurrence before any previous u(t) sample exists.
            if it < 1 {
                err!("maketu: u(t) requires a previous sample");
            }
            uoft[it] = uoft[it - 1] + max(dumin, du);
        } else if t > t2 {
            if it < 2 {
                err!("maketu: u(t) requires two previous samples");
            }
            uoft[it] = (2.0 * uoft[it - 1] as f64 - uoft[it - 2] as f64) as f32;
        }
        t += dt;
    }
    let mut dumin = uoft[1] - uoft[0];
    for it in 1..nt {
        dumin = min(dumin, uoft[it] - uoft[it - 1]);
    }
    let fu = 0.0f32;
    let eu = uoft[nt - 1];
    let mut du = (dumin as f64 / min(1.0, 2.0 * fmax as f64 * dt as f64)) as f32;
    let mut nu = 1 + nint(((eu - fu) / du) as f64);
    if nu > numax {
        nu = numax;
        du = (eu - fu) / (nu - 1) as f32;
    }
    // C passes nu to ealloc1float(size_t). A negative stretched length
    // (e.g. a positive delrt longer than a short trace) fails with ENOMEM;
    // use a fallible allocation so this path reports SU's error, not a panic.
    let mut tofu = Vec::<f32>::new();
    if tofu.try_reserve_exact(nu as usize).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    tofu.resize(nu as usize, 0.0);
    yxtoxy(dt, ft, &uoft, du, fu, ft, et, &mut tofu);
    (uoft, nu as usize, du, fu, tofu, t1)
}

/* table of z/h vs b/h; Mohammed Alfaraj, 01/08/92. */
fn table(ntable: usize, gamma: f32, zoh: &mut [f32], boh: &mut [f32]) {
    let f = ((1.0 - gamma as f64) / (1.0 + gamma as f64)) as f32;
    let alpha = (1.0 + 1.0 / (gamma * gamma) as f64) as f32;
    let beta = (1.0 - 1.0 / (gamma * gamma) as f64) as f32;
    let (bstart, db) = if gamma < 1.0 {
        let bstart = 1.0f32;
        (bstart, -(bstart - f) / (ntable as i32 + 2) as f32)
    } else {
        let bstart = -1.0f32;
        (bstart, (f - bstart) / (ntable as i32 + 2) as f32)
    };
    let mut b = bstart;
    for i in 0..ntable {
        let z = (((1.0f32 + b) * (1.0f32 - b * b)) as f64
            * (2.0 - alpha as f64 - (beta * b) as f64)
            / (alpha as f64 - 2.0 + b as f64 * ((beta - alpha) as f64 - 2.0)
                - (beta * b * b) as f64)) as f32;
        if z < 0.0 {
            err!("negative sqrt when generating table");
        }
        zoh[i] = (z as f64).sqrt() as f32;
        boh[i] = b;
        b += db;
    }
}

/* convert b/h to a function of NMO time; Mohammed Alfaraj, 01/08/92. */
#[allow(clippy::too_many_arguments)]
fn getbbohitn(
    x: f32,
    itmin: usize,
    nt: usize,
    dt: f32,
    v: &[f32],
    boh: &[f32],
    zoh: &[f32],
    gamma: f32,
) -> (Vec<f32>, Vec<i32>) {
    let mut bboh = vec![0.0f32; nt];
    let mut itn = vec![0i32; nt];
    let alpha = (1.0 + 1.0 / (gamma * gamma) as f64) as f32;
    let beta = (1.0 - 1.0 / (gamma * gamma) as f64) as f32;
    let tossq = (2.0 / ((1.0 + 1.0 / gamma as f64) * (1.0 + 1.0 / gamma as f64))) as f32;
    let tog = (2.0 / gamma as f64) as f32;
    let xsq = x * x;
    let mut j = 0;
    let mut t = itmin as f32 * dt;
    for it in itmin..nt {
        let xov = x / v[it];
        let nz = t / xov;
        let mut gotit = false;
        for i in j..boh.len() {
            if zoh[i] >= nz {
                bboh[it] = boh[i];
                let temp = (t * t / (1.0f32 - boh[i] * boh[i])
                    - xov * xov * (tog / (alpha + beta * boh[i]) - 1.0f32))
                    * (alpha + beta * boh[i])
                    * tossq;
                itn[it] = if temp < 0.0 {
                    0
                } else {
                    nint((temp as f64).sqrt() / dt as f64)
                };
                j = i;
                gotit = true;
                break;
            }
        }
        if gotit {
            t += dt;
            continue;
        }
        for j in it..nt {
            bboh[j] = boh[boh.len() - 1];
            let temp = (t * t / (1.0f32 - bboh[j] * bboh[j])
                - (xsq / (v[j] * v[j])) * (tog / (alpha + beta * bboh[j]) - 1.0f32))
                * (alpha + beta * bboh[j])
                * tossq;
            itn[j] = nint((temp as f64).sqrt() / dt as f64);
            t += dt;
        }
        break;
    }
    (bboh, itn)
}

/* calculate stretch/squeeze factors; Mohammed Alfaraj, 01/08/92. */
fn stretchfactor(sdmo: f32, gamma: f32) -> (f32, f32) {
    let (mut s1, mut s2) = if gamma <= 1.0 {
        (
            (0.62 + 0.633 * ((gamma as f64).powf(-3.0) - 1.0)) as f32,
            (0.62 + 0.615 * ((2.0 - gamma as f64).powf(-7.0) - 1.0)) as f32,
        )
    } else {
        (
            (0.62 + 0.615 * ((2.0 - 1.0 / gamma as f64).powf(-7.0) - 1.0)) as f32,
            (0.62 + 0.633 * ((1.0 / gamma as f64).powf(-3.0) - 1.0)) as f32,
        )
    };
    s1 *= sdmo;
    s2 *= sdmo;
    (s1, s2)
}
