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

//! SUSTOLT - Stolt migration for stacked data or common-offset gathers.
//! Port of su/main/migration_inversion/sustolt.c (Revision 1.22).
//! Credits: CWP: Dave Hale c. 1990.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::complex::Complex;
use oxwave_cwp::cwp::{max, min, nint, PI};
use oxwave_f_cubicspline::cmonot;
use oxwave_inverse_interp::yxtoxy;
use oxwave_pfa_fft::{npfao, npfar, pfa2cr, pfa2rc, pfacc};
use oxwave_sinc_interp::{ints8c, ints8r};
use oxwave_tempfile::temporary_filename;
use std::sync::atomic::{AtomicPtr, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSTOLT - Stolt migration for stacked data or common-offset gathers\t",
    "\t\t\t\t\t\t\t\t\t",
    " sustolt <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [...]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " cdpmin=\t\t  minimum cdp (integer number) in dataset\t",
    " cdpmax=\t\t  maximum cdp (integer number) in dataset\t",
    " dxcdp=\t\t  distance between adjacent cdp bins (m)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " noffmix=1\t\tnumber of offsets to mix (for unstacked data only)",
    " tmig=0.0\t\ttimes corresponding to rms velocities in vmig (s)",
    " vmig=1500.0\t\trms velocities corresponding to times in tmig (m/s)",
    " smig=1.0\t\tstretch factor (0.6 typical if vrms increasing)",
    " vscale=1.0\t\tscale factor to apply to velocities\t\t",
    " fmax=Nyquist\t\tmaximum frequency in input data (Hz)\t\t",
    " lstaper=0\t\tlength of side tapers (# of traces)\t\t",
    " lbtaper=0\t\tlength of bottom taper (# of samples)\t\t",
    " verbose=0\t\t=1 for diagnostic print\t\t\t\t",
    " tmpdir=\t\tif non-empty, use the value as a directory path\t",
    "\t\t\tprefix for storing temporary files; else if the\t",
    "\t\t\tthe CWP_TMPDIR environment variable is set use\t",
    "\t\t\tits value for the path; else use tmpfile()\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " If unstacked traces are input, they should be NMO-corrected and sorted",
    " into common-offset  gathers.  One common-offset gather ends and another",
    " begins when the offset field of the trace headers changes. If both\t",
    " NMO and DMO are applied, then this is equivalent to prestack time \t",
    " migration (though the velocity profile is assumed v(t), only).\t",
    "\t\t\t\t\t\t\t\t\t",
    " The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
    " the cdp location expressed in units of meters or feet.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The number of offsets to mix (noffmix) should be specified for\t",
    " unstacked data only.\tnoffmix should typically equal the ratio of the\t",
    " shotpoint spacing to the cdp spacing.\t This choice ensures that every\t",
    " cdp will be represented in each offset mix.  Traces in each mix will\t",
    " contribute through migration to other traces in adjacent cdps within\t",
    " that mix.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The tmig and vmig arrays specify a velocity function of time that is\t",
    " used to implement Stolt's stretch for depth-variable velocity.  The\t",
    " stretch factor smig is often referred to as the \"W\" factor.\t\t",
    " The times in tmig must be monotonically increasing.\t\t\t",
];

// The C globals are only used by the named-file signal cleanup handler.
static HEADERFP: AtomicPtr<CFile> = AtomicPtr::new(std::ptr::null_mut());
static HEADERFILE: AtomicPtr<String> = AtomicPtr::new(std::ptr::null_mut());

// OS FFI preserves the original program's named-file interrupt cleanup.
unsafe extern "C" {
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

// Storage for the C px row pointers and the offset-mix accumulator q.
struct Gathers {
    px: Vec<f32>,
    q: Vec<f32>,
    nxfft: usize,
    stride: usize,
    nt: usize,
}

impl Gathers {
    fn trace(&mut self, ix: usize) -> &mut [f32] {
        if ix < self.nxfft {
            &mut self.px[ix * self.stride..(ix + 1) * self.stride]
        } else {
            // Negative velocity/spacing can make nxfft < nx. The reference
            // overreads the px pointer table into the immediately following
            // q pointer table: px[nxfft + j] aliases q[j]. Preserve the
            // observed alias (including sequential accumulation) safely.
            let j = ix - self.nxfft;
            &mut self.q[j * self.nt..(j + 1) * self.nt]
        }
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get information from the first header */
    if gettr(&mut tr) == 0 {
        err!("cannot get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    if ft != 0.0 {
        err!("cannot handle non-zero time of first sample");
    }
    let mut offset = tr.offset as f32;

    /* get parameters */
    let mut cdpmin = 0i32;
    let mut cdpmax = 0i32;
    let mut dx = 0.0f32;
    let mut noffmix = 1i32;
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
        noffmix = 1;
    }

    /* get times and velocities */
    let ntmig = (countparval("tmig") + 1) as usize;
    let mut tmig = vec![0.0f32; ntmig];
    if getparfloat("tmig", &mut tmig[..]) == 0 {
        tmig[0] = 0.0;
    }
    tmig[ntmig - 1] = nt as f32 * dt + ft;
    let nvmig = (countparval("vmig") + 1) as usize;
    if nvmig != ntmig {
        err!("number of tmig and vmig must be equal");
    }
    let mut vmig = vec![0.0f32; nvmig];
    if getparfloat("vmig", &mut vmig[..]) == 0 {
        vmig[0] = 1500.0;
    }
    // With no vmig= the C reads vmig[-1]. In the reference build this is
    // zero: the default golden is byte-identical to tmig=0 vmig=0.
    // Use the observed value rather than an out-of-bounds Rust read.
    vmig[nvmig - 1] = if nvmig >= 2 { vmig[nvmig - 2] } else { 0.0 };
    for itmig in 1..ntmig {
        if tmig[itmig] <= tmig[itmig - 1] {
            err!("tmig must increase monotonically");
        }
    }
    let mut smig = 1.0f32;
    let mut fmax = 0.0f32;
    let mut vscale = 1.0f32;
    let mut lstaper = 0i32;
    let mut lbtaper = 0i32;
    let mut verbose = 0i32;
    if getparfloat("smig", &mut smig) == 0 {
        smig = 1.0;
    }
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (0.5 / dt as f64) as f32;
    }
    fmax = min(fmax as f64, 0.5 / dt as f64) as f32;
    if getparfloat("vscale", &mut vscale) == 0 {
        vscale = 1.0;
    }
    if getparint("lstaper", &mut lstaper) == 0 {
        lstaper = 0;
    }
    if getparint("lbtaper", &mut lbtaper) == 0 {
        lbtaper = 0;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* Look for user-supplied tmpdir */
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    unsafe extern "C" {
        fn access(path: *const std::ffi::c_char, mode: i32) -> i32;
    }
    if !tmpdir.is_empty() {
        let path = std::ffi::CString::new(tmpdir.as_bytes()).unwrap();
        // WRITE_OK == W_OK == 2. This is the OS permission check, not a
        // replacement for a CWP library routine.
        if unsafe { access(path.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();

    /* make uniformly sampled rms velocity function of time */
    let (v, vmin, vmax) = makev(&tmig, &vmig, vscale, nt, dt, ft);
    let vstolt = vmin;
    let (ut, nu, du, tu) = makeut(vstolt, fmax, &v, nt, dt);
    drop(v);
    let nx = cdpmax.wrapping_sub(cdpmin).wrapping_add(1) as usize;
    let nxpad = (0.5 * vmax as f64 * nt as f64 * dt as f64 / dx as f64) as i32;
    let nxfft = npfar((nx as i32).wrapping_add(nxpad)) as usize;
    let nk = nxfft / 2 + 1;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;

    // C aliases px and pk in one allocation. Separate buffers keep the same
    // flattened layout and transform inputs, with no simultaneous aliasing.
    let stride = max(nu, nt);
    let mut pk = vec![Complex::default(); stride * nk];
    let mut gathers = Gathers {
        px: vec![0.0f32; stride * nxfft],
        q: vec![0.0f32; nt * nx],
        nxfft,
        stride,
        nt,
    };

    /* open temporary file for headers */
    let headerfile;
    let headerfp;
    if tmpdir.is_empty() {
        headerfile = String::new();
        headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        headerfile = temporary_filename(&tmpdir);
        headerfp = efopen(&headerfile, "w+");
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
    }
    if !tmpdir.is_empty() {
        // These two locals remain in place until exit; only the signal
        // handler reads their pointers. Like C, the handler uses stdio.
        HEADERFP.store(&headerfp as *const CFile as *mut CFile, Ordering::Relaxed);
        HEADERFILE.store(
            &headerfile as *const String as *mut String,
            Ordering::Relaxed,
        );
        unsafe {
            signal(2, closefiles); // SIGINT
            signal(3, closefiles); // SIGQUIT
            signal(1, closefiles); // SIGHUP
            signal(15, closefiles); // SIGTERM
        }
    }
    let mut oldoffset = offset;
    let mut gottrace = true;
    let mut ntrace = 0usize;
    let mut noff = 0i32;

    /* loop over traces */
    loop {
        if gottrace {
            offset = tr.offset as f32;
        }
        if (gottrace && offset != oldoffset) || !gottrace {
            let scale = (1.0 / nxfft as f64) as f32;
            for ix in 0..nx {
                taper(
                    lstaper,
                    lbtaper,
                    nx as i32,
                    ix as i32,
                    nt,
                    gathers.trace(ix),
                );
            }
            if nu != nt {
                if let Some(tu) = &tu {
                    let mut temp = vec![0.0f32; nu];
                    for ix in 0..nx {
                        let trace = gathers.trace(ix);
                        ints8r(dt, 0.0, &trace[..nt], 0.0, 0.0, tu, &mut temp);
                        for iu in 0..nu {
                            trace[iu] = temp[iu];
                        }
                    }
                }
            }
            /* Fourier transform p(u,x) to p(u,k) */
            pfa2rc(-1, 2, stride as i32, nxfft as i32, &gathers.px, &mut pk);
            let mut k = dk;
            for ik in 1..nk {
                stolt1k(
                    k,
                    vstolt,
                    smig,
                    fmax,
                    stride,
                    du,
                    &mut pk[ik * stride..(ik + 1) * stride],
                );
                k += dk;
            }
            pfa2cr(1, 2, stride as i32, nxfft as i32, &pk, &mut gathers.px);
            for jx in 0..nx {
                let trace = gathers.trace(jx);
                for iu in 0..nu {
                    trace[iu] *= scale;
                }
            }
            if nu != nt {
                if let Some(ut) = &ut {
                    let mut temp = vec![0.0f32; nt];
                    for ix in 0..nx {
                        let trace = gathers.trace(ix);
                        ints8r(du, 0.0, &trace[..nu], 0.0, 0.0, ut, &mut temp);
                        for it in 0..nt {
                            trace[it] = temp[it];
                        }
                    }
                }
            }
            for ix in 0..nx {
                for it in 0..nt {
                    let sample = gathers.trace(ix)[it];
                    gathers.q[ix * nt + it] += sample;
                }
            }
            for ix in 0..nxfft {
                for iu in 0..nu {
                    gathers.px[ix * stride + iu] = 0.0;
                }
            }
            noff += 1;
        }
        if noff == noffmix || !gottrace {
            efseeko(&headerfp, 0, SEEK_SET);
            for _itrace in 0..ntrace {
                efread(&mut tro.hdr.as_bytes_mut()[..], &headerfp);
                let ix = tro.cdp.wrapping_sub(cdpmin) as usize;
                for it in 0..nt {
                    tro.data[it] = gathers.q[ix * nt + it];
                }
                puttr(&mut tro);
            }
            if verbose != 0 {
                fprintf!(
                    stderr(),
                    "\tCompleted mix of %d offsets with %d traces\n",
                    noff,
                    ntrace as i32
                );
            }
            if !gottrace {
                break;
            }
            efseeko(&headerfp, 0, SEEK_SET);
            noff = 0;
            ntrace = 0;
            gathers.q.fill(0.0);
        }
        if verbose != 0 && oldoffset != offset {
            warn!("migrating offset = %f", offset);
        }
        if tr.cdp >= cdpmin && tr.cdp <= cdpmax {
            efwrite(&tr.hdr.as_bytes()[..], &headerfp);
            ntrace += 1;
            oldoffset = offset;
            let ix = tr.cdp.wrapping_sub(cdpmin) as usize;
            let trace = gathers.trace(ix);
            for it in 0..nt {
                trace[it] = tr.data[it];
            }
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

/* make uniformly sampled rms velocity function v(t) for migration */
fn makev(
    tmig: &[f32],
    vmig: &[f32],
    vscale: f32,
    nt: usize,
    dt: f32,
    ft: f32,
) -> (Vec<f32>, f32, f32) {
    let mut vmigd = vec![[0.0f32; 4]; tmig.len()];
    cmonot(tmig, vmig, &mut vmigd);
    let mut vel = vec![0.0f32; nt];
    let mut t = ft;
    for it in 0..nt {
        intcub(0, tmig, &vmigd, &[t], &mut vel[it..it + 1]);
        t += dt;
    }
    for it in 0..nt {
        vel[it] *= vscale;
    }
    let mut velmin = vel[0];
    let mut velmax = vel[0];
    for it in 1..nt {
        velmin = min(velmin, vel[it]);
        velmax = max(velmax, vel[it]);
    }
    (vel, velmin, velmax)
}

/* Compute u(t) and t(u) for Stolt stretch. */
fn makeut(
    vstolt: f32,
    fmax: f32,
    vrms: &[f32],
    nt: usize,
    dt: f32,
) -> (Option<Vec<f32>>, usize, f32, Option<Vec<f32>>) {
    let mut it = 1usize;
    while it < nt {
        if vrms[it] != vrms[0] {
            break;
        }
        it += 1;
    }
    if it == nt {
        return (None, nt, dt, None);
    }
    let mut u = vec![0.0f32; nt];
    makeu(vstolt, vrms, nt, dt, &mut u);
    let mut duu = f32::MAX;
    for it in 1..nt {
        let delu = u[it] - u[it - 1];
        if delu < duu {
            duu = delu;
        }
    }
    let umax = u[nt - 1];
    duu = (duu as f64 / (2.0 * fmax as f64 * dt as f64)) as f32;
    let nuu = 1i32.wrapping_add(nint((umax / duu) as f64)) as usize;
    let mut t = vec![0.0f32; nuu];
    yxtoxy(
        dt,
        0.0,
        &u,
        duu,
        0.0,
        0.0,
        (nt as i32 - 1) as f32 * dt,
        &mut t,
    );
    (Some(u), nuu, duu, Some(t))
}

/* Compute u(t) via the trapezoidal rule. */
fn makeu(vstolt: f32, v: &[f32], nt: usize, dt: f32, u: &mut [f32]) {
    let scale = (2.0 / (vstolt * vstolt) as f64) as f32;
    let mut sum = 0.0f32;
    u[0] = 0.0;
    let mut t = dt;
    for it in 1..nt {
        sum = (sum as f64
            + 0.5 * dt as f64 * (t * v[it] * v[it] + (t - dt) * v[it - 1] * v[it - 1]) as f64)
            as f32;
        u[it] = ((scale * sum) as f64).sqrt() as f32;
        t += dt;
    }
}

/* Stolt's migration for one wavenumber k (p and q equivalenced). */
fn stolt1k(k: f32, v: f32, mut s: f32, fmax: f32, nt: usize, dt: f32, p: &mut [Complex]) {
    if s != 1.0 {
        s = (2.0 - s as f64) as f32;
    }
    let vko2s = (0.25 * v as f64 * v as f64 * k as f64 * k as f64) as f32;
    let wmax = (2.0 * PI * min(fmax as f64, 0.5 / dt as f64)) as f32;
    let nw = ((nt as i32 * 2) as f64 / 0.6) as i32;
    let nw = npfao(nw, nw.wrapping_mul(2)) as usize;
    let dw = (2.0 * PI / (nw as f32 * dt) as f64) as f32;
    let fw = (-PI / dt as f64) as f32;
    let ntau = nt;
    let dtau = dt;
    let nwtau = npfao(ntau as i32, (ntau as i32).wrapping_mul(2)) as usize;
    let dwtau = (2.0 * PI / (nwtau as f32 * dtau) as f64) as f32;
    let mut fwtau = (-PI / dtau as f64) as f32;
    fwtau = (fwtau as f64 + 0.001 * dwtau as f64) as f32;
    let wtauh = max(0.0f64, (wmax * wmax - s * vko2s) as f64).sqrt() as f32;
    let iwtauh = max(
        0,
        min(nwtau as i32 - 1, nint(((wtauh - fwtau) / dwtau) as f64)),
    ) as usize;
    let iwtaul = max(
        0,
        min(nwtau as i32 - 1, nint(((-wtauh - fwtau) / dwtau) as f64)),
    ) as usize;
    let _wtauh = fwtau + iwtauh as f32 * dwtau;
    let wtaul = fwtau + iwtaul as f32 * dwtau;
    let czero = Complex::default();
    let mut pp = vec![czero; nw];
    let mut qq = vec![czero; nwtau];
    let mut wwtau = vec![0.0f32; nwtau];
    for it in (0..nt).step_by(2) {
        pp[it] = p[it];
    }
    for it in (1..nt).step_by(2) {
        pp[it].r = -p[it].r;
        pp[it].i = -p[it].i;
    }
    pfacc(1, &mut pp);
    pp[0] = czero;
    let mut iwtau = iwtaul;
    let mut wtau = wtaul;
    if s == 1.0 {
        while wtau < 0.0 {
            wwtau[iwtau] = -((wtau * wtau + vko2s) as f64).sqrt() as f32;
            iwtau += 1;
            wtau += dwtau;
        }
        while iwtau <= iwtauh {
            wwtau[iwtau] = ((wtau * wtau + vko2s) as f64).sqrt() as f32;
            iwtau += 1;
            wtau += dwtau;
        }
    } else {
        let a = (1.0 / s as f64) as f32;
        let b = (1.0 - a as f64) as f32;
        while wtau < 0.0 {
            wwtau[iwtau] =
                ((b * wtau) as f64 - a as f64 * ((wtau * wtau + s * vko2s) as f64).sqrt()) as f32;
            iwtau += 1;
            wtau += dwtau;
        }
        while iwtau <= iwtauh {
            wwtau[iwtau] =
                ((b * wtau) as f64 + a as f64 * ((wtau * wtau + s * vko2s) as f64).sqrt()) as f32;
            iwtau += 1;
            wtau += dwtau;
        }
    }
    ints8c(
        dw,
        fw,
        &pp,
        czero,
        czero,
        &wwtau[iwtaul..iwtauh + 1],
        &mut qq[iwtaul..iwtauh + 1],
    );
    let fftscl = (1.0 / nwtau as f64) as f32;
    let mut wtau = wtaul;
    if s == 1.0 {
        for iwtau in iwtaul..=iwtauh {
            let scale = fftscl * wtau / wwtau[iwtau];
            qq[iwtau].r *= scale;
            qq[iwtau].i *= scale;
            wtau += dwtau;
        }
    } else {
        let a = (1.0 / (s * s) as f64) as f32;
        let b = (1.0 - 1.0 / s as f64) as f32;
        for iwtau in iwtaul..=iwtauh {
            let scale = fftscl * (b + a * wtau / (wwtau[iwtau] - b * wtau));
            qq[iwtau].r *= scale;
            qq[iwtau].i *= scale;
            wtau += dwtau;
        }
    }
    for iwtau in 0..iwtaul {
        qq[iwtau] = czero;
    }
    for iwtau in iwtauh + 1..nwtau {
        qq[iwtau] = czero;
    }
    pfacc(-1, &mut qq);
    for itau in (0..ntau).step_by(2) {
        p[itau] = qq[itau];
    }
    for itau in (1..ntau).step_by(2) {
        p[itau].r = -qq[itau].r;
        p[itau].i = -qq[itau].i;
    }
}

/* Taper traces near left and right sides of trace window. */
fn taper(lxtaper: i32, lbtaper: i32, nx: i32, ix: i32, nt: usize, trace: &mut [f32]) {
    let xtaper = if ix < lxtaper {
        (0.54 + 0.46 * (PI * (lxtaper - ix) as f64 / lxtaper as f64).cos()) as f32
    } else if ix >= nx - lxtaper {
        (0.54 + 0.46 * (PI * (lxtaper + ix + 1 - nx) as f64 / lxtaper as f64).cos()) as f32
    } else {
        1.0
    };
    if xtaper != 1.0 {
        for it in 0..nt {
            trace[it] *= xtaper;
        }
    }
    for it in max(0, nt as i32 - lbtaper)..nt as i32 {
        trace[it as usize] = (trace[it as usize] as f64
            * (0.54 + 0.46 * (PI * (lbtaper + it + 1 - nt as i32) as f64 / lbtaper as f64).cos()))
            as f32;
    }
}

/* for graceful interrupt termination */
extern "C" fn closefiles(_sig: i32) {
    let headerfp = HEADERFP.load(Ordering::Relaxed);
    let headerfile = HEADERFILE.load(Ordering::Relaxed);
    // Registered only after both locals have been initialized in main.
    unsafe {
        efclose(&*headerfp);
        eremove(&*headerfile);
    }
    exit(1);
}
