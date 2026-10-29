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

//! SUDMOFK - DMO via F-K domain (log-stretch) for common-offset gathers.
//! Port of su/main/dip_moveout/sudmofk.c (Revision 1.20).
//! Credits (from the C source): CWP: Dave.
//! Technical reference: Dip-Moveout Processing - SEG Course Notes,
//! Dave Hale, 1988.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::cwp::{abs, max, min, nint, PI};
use oxwave_f_cubicspline::cmonot;
use oxwave_inverse_interp::yxtoxy;
use oxwave_pfa_fft::{npfa, npfar, pfa2cr, pfa2rc, pfacc, Complex};
use oxwave_sinc_interp::ints8c;
use oxwave_tempfile::temporary_filename;
use std::ffi::CString;
use std::sync::atomic::{AtomicPtr, Ordering};
use std::sync::OnceLock;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUDMOFK - DMO via F-K domain (log-stretch) method for common-offset gathers",
    "\t\t\t\t\t\t\t\t\t",
    " sudmofk <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [...]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " cdpmin\tminimum cdp (integer number) for which to apply DMO\t",
    " cdpmax\tmaximum cdp (integer number) for which to apply DMO\t",
    " dxcdp\t\tdistance between adjacent cdp bins (m)\t\t\t",
    " noffmix\tnumber of offsets to mix (see notes)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " tdmo=0.0\ttimes corresponding to rms velocities in vdmo (s)\t",
    " vdmo=1500.0\trms velocities corresponding to times in tdmo (m/s)\t",
    " sdmo=1.0\tDMO stretch factor; try 0.6 for typical v(z)\t\t",
    " fmax=0.5/dt\tmaximum frequency in input traces (Hz)\t\t\t",
    " verbose=0\t=1 for diagnostic print\t\t\t\t\t",
    " tmpdir=\tif non-empty, use the value as a directory path\tprefix\t",
    "\t\tfor storing temporary files; else if the CWP_TMPDIR\t",
    "\t\tenvironment variable is set use\tits value for the path;\t",
    "\t\telse use tmpfile()\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input traces should be sorted into common-offset gathers.  One common- ",
    " offset gather ends and another begins when the offset field of the trace",
    " headers changes.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
    " the cdp location expressed in units of meters or feet.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The number of offsets to mix (noffmix) should typically no smaller than",
    " the ratio of the shotpoint spacing to the cdp spacing.  This choice\t",
    " ensures that every cdp will be represented in each offset mix.  Traces ",
    " in each mix will contribute through DMO to other traces in adjacent cdps",
    " within that mix. (Values of noffmix 2 or 3 times the ratio of shotpoint",
    " spacing to the cdp spacing often yield better results.)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The tdmo and vdmo arrays specify a velocity function of time that is\t",
    " used to implement a first-order correction for depth-variable velocity.",
    " The times in tdmo must be monotonically increasing.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For each offset, the minimum time at which a non-zero sample exists is ",
    " used to determine a mute time.  Output samples for times earlier than this",
    " mute time will be zeroed.  Computation time may be significantly reduced",
    " if the input traces are zeroed (muted) for early times at large offsets.",
];

// The C globals exist so closefiles can remove the header file on a signal.
static HEADER_FILE: OnceLock<String> = OnceLock::new();
static HEADER_FP: AtomicPtr<CFile> = AtomicPtr::new(std::ptr::null_mut());

unsafe extern "C" {
    fn access(path: *const std::ffi::c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

extern "C" fn closefiles(_: i32) {
    let fp = HEADER_FP.load(Ordering::Acquire);
    if !fp.is_null() {
        // The single main thread installs a leaked handle before using it.
        // As in C, the handler closes it and immediately exits.
        unsafe {
            efclose(&*fp);
        }
    }
    if let Some(name) = HEADER_FILE.get() {
        eremove(name);
    }
    exit(EXIT_FAILURE);
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get information from the first header.
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut offset = tr.offset as f32;

    let mut cdpmin = 0i32;
    if getparint("cdpmin", &mut cdpmin) == 0 {
        err!("must specify cdpmin");
    }
    let mut cdpmax = 0i32;
    if getparint("cdpmax", &mut cdpmax) == 0 {
        err!("must specify cdpmax");
    }
    if cdpmin > cdpmax {
        err!("cdpmin must not be greater than cdpmax");
    }
    let mut dx = 0.0f32;
    if getparfloat("dxcdp", &mut dx) == 0 {
        err!("must specify dxcdp");
    }
    let mut noffmix = 0i32;
    if getparint("noffmix", &mut noffmix) == 0 {
        err!("must specify noffmix");
    }
    let ntdmo = max(countparval("tdmo"), 1) as usize;
    let mut tdmo = vec![0.0f32; ntdmo];
    if getparfloat("tdmo", &mut tdmo) == 0 {
        tdmo[0] = 0.0;
    }
    let nvdmo = max(countparval("vdmo"), 1) as usize;
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
    let mut sdmo = 0.0f32;
    if getparfloat("sdmo", &mut sdmo) == 0 {
        sdmo = 1.0;
    }
    let mut fmax = 0.0f32;
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (0.5 / dt as f64) as f32;
    }
    let mut verbose = 0i32;
    getparint("verbose", &mut verbose);

    // Look for user-supplied tmpdir.
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let directory = CString::new(tmpdir.split('\0').next().unwrap()).unwrap();
        if unsafe { access(directory.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();

    let mut vrms = vec![0.0f32; nt];
    mkvrms(&tdmo, &vdmo, nt, dt, ft, &mut vrms);

    let nx = cdpmax.wrapping_sub(cdpmin).wrapping_add(1) as usize;
    let nxpad = (0.5 * abs(offset / dx) as f64) as i32;
    let nxfft = npfar((nx as i32).wrapping_add(nxpad)) as usize;
    let mut p = vec![0.0f32; nt * (nxfft + 2)];
    let mut q = vec![0.0f32; nt * nx];

    // Open temporary file for headers.
    let headerfp;
    let mut headerfile = String::new();
    if tmpdir.is_empty() {
        headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        headerfile = temporary_filename(&tmpdir);
        let _ = HEADER_FILE.set(headerfile.clone());
        // These four signal numbers are shared by Darwin and Linux.
        for sig in [2, 3, 1, 15] {
            unsafe {
                signal(sig, closefiles);
            }
        }
        headerfp = efopen(&headerfile, "w+");
        HEADER_FP.store(Box::into_raw(Box::new(headerfp.clone())), Ordering::Release);
        if verbose != 0 {
            warn!("putting temporary header file in %s", tmpdir);
        }
    }

    let mut oldoffset = offset;
    let mut gottrace = true;
    let mut ntrace = 0i32;
    let mut noff = 0i32;

    loop {
        if gottrace {
            offset = tr.offset as f32;
        }

        // If an offset is complete, correct it and add it to the mix.
        if (gottrace && offset != oldoffset) || !gottrace {
            dmooff(oldoffset, fmax, sdmo, nx, dx, nt, dt, ft, &vrms, &mut p);
            for ix in 0..nx {
                for it in 0..nt {
                    q[ix * nt + it] += p[ix * nt + it];
                }
            }
            noff += 1;
            p = Vec::new();
            if offset != oldoffset {
                let nxpad = (0.5 * abs(offset / dx) as f64) as i32;
                let nxfft = npfar((nx as i32).wrapping_add(nxpad)) as usize;
                p = vec![0.0f32; nt * (nxfft + 2)];
            }
        }

        // If a mix of offsets is complete, output its saved headers and data.
        if noff == noffmix || !gottrace {
            erewind(&headerfp);
            for _itrace in 0..ntrace {
                efread(tro.hdr.as_bytes_mut(), &headerfp);
                let ix = tro.cdp.wrapping_sub(cdpmin) as usize;
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
            erewind(&headerfp);
            noff = 0;
            ntrace = 0;
            for ix in 0..nx {
                for it in 0..nt {
                    q[ix * nt + it] = 0.0;
                }
            }
        }

        if tr.cdp >= cdpmin && tr.cdp <= cdpmax {
            efwrite(tr.hdr.as_bytes(), &headerfp);
            ntrace += 1;
            oldoffset = offset;
            let ix = tr.cdp.wrapping_sub(cdpmin) as usize;
            p[ix * nt..(ix + 1) * nt].copy_from_slice(&tr.data[..nt]);
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }

    efclose(&headerfp);
    if !headerfile.is_empty() {
        eremove(&headerfile);
    }
    exit(CWP_Exit());
}

// Make uniformly sampled vrms(t) for DMO.
fn mkvrms(tdmo: &[f32], vdmo: &[f32], nt: usize, dt: f32, ft: f32, vrms: &mut [f32]) {
    let mut vdmod = vec![[0.0f32; 4]; tdmo.len()];
    cmonot(tdmo, vdmo, &mut vdmod);
    let mut t = ft;
    for it in 0..nt {
        intcub(0, tdmo, &vdmod, &[t], &mut vrms[it..it + 1]);
        t += dt;
    }
}

// Perform DMO in the F-K domain for one offset. ptx is a contiguous,
// zero-padded (nxfft + 2)-by-nt array, as in the C allocation.
#[allow(clippy::too_many_arguments)]
fn dmooff(
    offset: f32,
    fmax: f32,
    mut sdmo: f32,
    nx: usize,
    dx: f32,
    nt: usize,
    dt: f32,
    ft: f32,
    vrms: &[f32],
    ptx: &mut [f32],
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

    let (uoft, nu, du, fu, tofu, tcon) = maketu(offset, itmin, fmax, nt, dt, ft, vrms);
    // Adjust for the nominal log-stretch error, then take its inverse.
    sdmo = (sdmo as f64 * 0.62) as f32;
    let osdmo = (1.0 / sdmo as f64) as f32;
    let nupad = (1.5 * sdmo as f64 * tcon as f64 / du as f64) as i32;
    let nufft = npfa((nu as i32).wrapping_add(nupad)) as usize;
    let nw = nufft;
    let dw = (2.0 * PI / (nufft as f32 * du) as f64) as f32;
    let mut pu = vec![Complex::default(); nufft];

    let nk = nxfft / 2 + 1;
    let dk = (2.0 * PI / abs(nxfft as f32 * dx) as f64) as f32;
    let kmax = (PI / abs(dx) as f64) as f32;
    let ikmax = nint((kmax / dk) as f64);
    // The safe FFT API uses separate buffers for C's aliased ptx/ptk storage.
    let mut ptk = vec![Complex::default(); nk * nt];
    let fftscl = nk as f32 / (ikmax + 1) as f32 / (nufft * nxfft) as f32;
    pfa2rc(-1, 2, nt as i32, nxfft as i32, ptx, &mut ptk);

    let czero = Complex::default();
    let mut k = 0.0f32;
    let mut ik = 0i32;
    while ik <= ikmax {
        let ki = ik as usize * nt;
        ints8c(dt, ft, &ptk[ki..ki + nt], czero, czero, &tofu, &mut pu);
        for iu in nu..nufft {
            pu[iu].r = 0.0;
            pu[iu].i = 0.0;
        }
        pfacc(1, &mut pu);

        let wmin = abs(0.5 * vrms[0] as f64 * k as f64) as f32;
        let wmax = abs(PI / du as f64) as f32;
        let iwmin = max(1, min(((nw - 1) / 2) as i32, nint((wmin / dw) as f64)));
        let iwmax = max(0, min(((nw - 1) / 2) as i32, nint((wmax / dw) as f64)));
        let wwscl =
            (osdmo as f64 * (k as f64 * 0.5 * offset as f64 / tcon as f64).powf(2.0)) as f32;
        pu[0].r = 0.0;
        pu[0].i = 0.0;
        for iw in 1..iwmin as usize {
            let iwn = nw - iw;
            pu[iw].r = 0.0;
            pu[iw].i = 0.0;
            pu[iwn].r = 0.0;
            pu[iwn].i = 0.0;
        }

        let mut w = iwmin as f32 * dw;
        for iw in iwmin..=iwmax {
            let iw = iw as usize;
            let iwn = nw - iw;
            let scales = (1.0 + (wwscl / (w * w)) as f64) as f32;
            let scale = (scales as f64).sqrt() as f32;
            let phase = ((sdmo * w * tcon) as f64 * (scale as f64 - 1.0)) as f32;
            let amp = (fftscl as f64 * (1.0 - sdmo as f64 + (sdmo / scale) as f64)) as f32;
            let fr = (amp as f64 * (phase as f64).cos()) as f32;
            let fi = (amp as f64 * (phase as f64).sin()) as f32;
            let pwr = pu[iw].r;
            let pwi = pu[iw].i;
            pu[iw].r = pwr * fr - pwi * fi;
            pu[iw].i = pwr * fi + pwi * fr;
            let pwr = pu[iwn].r;
            let pwi = pu[iwn].i;
            pu[iwn].r = pwr * fr + pwi * fi;
            pu[iwn].i = pwi * fr - pwr * fi;
            w += dw;
        }
        for iw in (iwmax + 1) as usize..=nw / 2 {
            let iwn = nw - iw;
            pu[iw].r = 0.0;
            pu[iw].i = 0.0;
            pu[iwn].r = 0.0;
            pu[iwn].i = 0.0;
        }
        pfacc(-1, &mut pu);
        ints8c(
            du,
            fu,
            &pu[..nu],
            czero,
            czero,
            &uoft[..nt],
            &mut ptk[ki..ki + nt],
        );
        ik += 1;
        k += dk;
    }
    for ik in ik as usize..nk {
        for it in 0..nt {
            ptk[ik * nt + it].r = 0.0;
            ptk[ik * nt + it].i = 0.0;
        }
    }
    pfa2cr(1, 2, nt as i32, nxfft as i32, &ptk, ptx);
}

// Make stretch and compress functions t(u) and u(t).
#[allow(clippy::too_many_arguments)]
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
    // C reads uoft[1] even for nt=1. Supply the zero-filled allocation tail
    // explicitly so this degenerate case does not cause a Rust panic.
    let mut uoft = vec![0.0f32; max(nt, 2)];
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
            // A negative ft can round t1 below ft, making C read uoft[-1]
            // on the first sample. Use zero for that out-of-allocation value
            // to match the reference output explicitly instead of panicking.
            let previous = if it == 0 { 0.0 } else { uoft[it - 1] };
            uoft[it] = previous + max(dumin, du);
        } else if t > t2 {
            // The same negative-time boundary can reach the extrapolation
            // branch before two samples exist (for example at zero offset).
            let previous = if it == 0 { 0.0 } else { uoft[it - 1] };
            let previous2 = if it < 2 { 0.0 } else { uoft[it - 2] };
            uoft[it] = (2.0 * previous as f64 - previous2 as f64) as f32;
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
    let mut nu = 1i32.wrapping_add(nint(((eu - fu) / du) as f64));
    if nu > numax {
        nu = numax;
        du = (eu - fu) / (nu - 1) as f32;
    }
    let mut tofu = vec![0.0f32; nu as usize];
    yxtoxy(dt, ft, &uoft[..nt], du, fu, ft, et, &mut tofu);
    (uoft, nu as usize, du, fu, tofu, t1)
}
