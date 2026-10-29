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

//! SUGAZMIGQ - Gazdag phase-shift migration of zero-offset data with attenuation Q.
//! Port of su/main/migration_inversion/sugazmig.c (Revision 1.21).
//! Credits: Constant Q attenuation correction by Chuck Oden, 5 May 2004;
//! CWP John Stockwell, 12 Oct 1992, based on a constant v version by Dave Hale.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cadd, cmplx, cmul, cneg, crmul, Complex};
use oxwave_cwp::cwp::PI;
use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::{npfa, npfar, pfa2cr, pfa2rc, pfacc};
use oxwave_tempfile::temporary_filename;
use std::ffi::{c_char, CString};
use std::sync::atomic::{AtomicBool, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUGAZMIGQ - SU version of Jeno GAZDAG's phase-shift migration \t",
    "\t     for zero-offset data, with attenuation Q.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sugazmig <infile >outfile vfile= [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt=from header(dt) or\t.004\ttime sampling interval\t\t\t",
    " dx=from header(d2) or 1.0\tmidpoint sampling interval\t\t",
    " ft=0.0\t\t\tfirst time sample\t\t\t",
    " ntau=nt(from data)\tnumber of migrated time samples\t\t\t",
    " dtau=dt(from header)\tmigrated time sampling interval\t\t\t",
    " ftau=ft\t\tfirst migrated time sample\t\t\t",
    " tmig=0.0\t\ttimes corresponding to interval velocities in vmig",
    " vmig=1500.0\tinterval velocities corresponding to times in tmig\t",
    " vfile=\t\tname of file containing velocities\t\t",
    " Q=1e6\t\t\tquality factor\t\t\t\t\t",
    " ceil=1e6\t\tgain ceiling beyond which migration ceases\t",
    "\t\t\t\t\t\t\t\t\t",
    " verbose=0\tverbose = 1 echoes information\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: ray bending effects not accounted for in this version.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The tmig and vmig arrays specify an interval velocity function of time.",
    " Linear interpolation and constant extrapolation is used to determine\t",
    " interval velocities at times not specified.  Values specified in tmig\t",
    " must increase monotonically.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Alternatively, interval velocities may be stored in a binary file\t",
    " containing one velocity for every time sample in the data that is to be",
    " migrated.  If vfile is specified, then the tmig and vmig arrays are ignored.",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat: Adding Q is a first attempt to address GPR issues.\t\t",
];

// Signals only set this lock-free flag. Stream operations and cleanup must
// run in normal execution, after any CFile borrows have been released.
static INTERRUPTED: AtomicBool = AtomicBool::new(false);

struct TempFiles {
    tracefp: CFile,
    headerfp: CFile,
    tracefile: String,
    headerfile: String,
}
unsafe extern "C" {
    fn access(path: *const c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn raise(sig: i32) -> i32;
}

extern "C" fn record_signal(_: i32) {
    INTERRUPTED.store(true, Ordering::Relaxed);
}

impl TempFiles {
    fn check_signal(&self) {
        if !INTERRUPTED.load(Ordering::Relaxed) {
            return;
        }
        efclose(&self.headerfp);
        efclose(&self.tracefp);
        if !self.headerfile.is_empty() {
            eremove(&self.headerfile);
        }
        if !self.tracefile.is_empty() {
            eremove(&self.tracefile);
        }
        exit(EXIT_FAILURE);
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get info from first trace.
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        if tr.dt != 0 {
            dt = (tr.dt as f64 / 1000000.0) as f32;
        } else {
            dt = 0.004;
            warn!("tr.dt not set, assuming dt=0.004");
        }
    }
    let mut dx = 0.0f32;
    if getparfloat("dx", &mut dx) == 0 {
        if tr.d2 != 0.0 {
            dx = tr.d2;
        } else {
            dx = 1.0;
            warn!("tr.d2 not set, assuming d2=1.0");
        }
    }

    let mut ft = 0.0f32;
    if getparfloat("ft", &mut ft) == 0 {
        ft = 0.0;
    }
    let mut ntau = 0i32;
    if getparint("ntau", &mut ntau) == 0 {
        ntau = nt;
    }
    CHECK_NT("ntau", ntau);
    let mut dtau = 0.0f32;
    if getparfloat("dtau", &mut dtau) == 0 {
        dtau = dt;
    }
    let mut ftau = 0.0f32;
    if getparfloat("ftau", &mut ftau) == 0 {
        ftau = ft;
    }
    let mut qual = 0.0f32;
    if getparfloat("Q", &mut qual) == 0 {
        qual = 1.0e6;
    }
    let mut ceil = 0.0f32;
    if getparfloat("ceil", &mut ceil) == 0 {
        ceil = 1.0e6;
    }
    let mut verbose = 0;
    if verbose != 0 {
        warn!("Q=%f ceil=%f", qual, ceil);
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty()
        && unsafe { access(CString::new(tmpdir.as_str()).unwrap().as_ptr(), 2) } != 0
    {
        err!("you can't write in %s (or it doesn't exist)", tmpdir);
    }

    // Store traces and headers in tempfiles while getting a count.
    let (tracefp, headerfp);
    let istmpdir = !tmpdir.is_empty();
    let (mut tracefile, mut headerfile) = (String::new(), String::new());
    if !istmpdir {
        tracefp = etmpfile();
        headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        tracefile = temporary_filename(&tmpdir);
        headerfile = temporary_filename(&tmpdir);
        unsafe {
            signal(2, record_signal); // SIGINT
            signal(3, record_signal); // SIGQUIT
            signal(1, record_signal); // SIGHUP
            signal(15, record_signal); // SIGTERM
        }
        tracefp = efopen(&tracefile, "w+");
        if tracefile == headerfile {
            // The tempfile crate returns the unchanged template when its
            // fixed /tmp allocation is denied. Keep the two streams separate
            // in that environment using an anonymous header file. Opening
            // tracefile above still preserves errors in the requested path.
            headerfp = etmpfile();
            headerfile.clear();
        } else {
            headerfp = efopen(&headerfile, "w+");
        }
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
    }
    let tempfiles = TempFiles {
        tracefp: tracefp.clone(),
        headerfp: headerfp.clone(),
        tracefile: tracefile.clone(),
        headerfile: headerfile.clone(),
    };
    let mut nx = 0i32;
    loop {
        tempfiles.check_signal();
        nx += 1;
        efwrite(std::slice::from_ref(&tr.hdr), &headerfp);
        efwrite(&tr.data[..nt as usize], &tracefp);
        tempfiles.check_signal();
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    tempfiles.check_signal();
    erewind(&tracefp);
    erewind(&headerfp);

    // Determine wavenumber sampling (for real to complex FFT).
    let nxfft = npfar(nx);
    let nk = nxfft / 2 + 1;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let fk = 0.0f32;
    let mut p = vec![0.0f32; nt as usize * nxfft as usize];
    // Negative alloc2 sizes fail in C; the later ealloc1float reports it.
    let mut q = vec![0.0f32; ntau.max(0) as usize * nxfft as usize];
    let mut cp = vec![Complex::default(); nt as usize * nk as usize];
    let mut cq = vec![Complex::default(); ntau.max(0) as usize * nk as usize];

    efread(&mut p[..(nt * nx) as usize], &tracefp);
    efclose(&tracefp);

    // Determine velocity function v(t).
    if ntau < 0 {
        err!("%s: malloc failed (Cannot allocate memory)", "ealloc.c");
    }
    let mut vt = vec![0.0f32; ntau as usize];
    let mut vfile = String::new();
    if getparstring("vfile", &mut vfile) == 0 {
        let mut ntmig = countparval("tmig");
        if ntmig == 0 {
            ntmig = 1;
        }
        let mut tmig = vec![0.0f32; ntmig as usize];
        if getparfloat("tmig", &mut tmig) == 0 {
            tmig[0] = 0.0;
        }
        let mut nvmig = countparval("vmig");
        if nvmig == 0 {
            nvmig = 1;
        }
        if nvmig != ntmig {
            err!("number of tmig and vmig must be equal");
        }
        let mut vmig = vec![0.0f32; nvmig as usize];
        if getparfloat("vmig", &mut vmig) == 0 {
            vmig[0] = 1500.0;
        }
        for itmig in 1..ntmig as usize {
            if tmig[itmig] <= tmig[itmig - 1] {
                err!("tmig must increase monotonically");
            }
        }
        let mut t = 0.0f32;
        for it in 0..ntau as usize {
            intlin(
                &tmig,
                &vmig,
                vmig[0],
                vmig[ntmig as usize - 1],
                &[t],
                &mut vt[it..it + 1],
            );
            t += dt;
        }
    } else {
        // C allocates ntau entries but reads nt: reserve the complete read so
        // a short output axis does not corrupt Rust's allocator metadata.
        vt.resize((ntau.max(nt)) as usize, 0.0);
        if let Some(fp) = fopen(&vfile, "r") {
            if fread(&mut vt[..nt as usize], &fp) != nt as usize {
                err!("cannot read %d velocities from file %s", nt, vfile);
            }
        } else {
            // C passes fopen's NULL result straight to fread.
            unsafe {
                raise(11);
            }
            exit(139);
        }
    }
    checkpars();
    tempfiles.check_signal();

    for ix in nx..nxfft {
        for it in 0..nt {
            p[(ix * nt + it) as usize] = 0.0;
        }
    }
    pfa2rc(-1, 2, nt, nxfft, &p, &mut cp);
    let mut k = fk;
    for ik in 0..nk as usize {
        gazdagvt(
            k,
            nt,
            dt,
            ft,
            ntau,
            dtau,
            ftau,
            &vt,
            &cp[ik * nt as usize..(ik + 1) * nt as usize],
            &mut cq[ik * ntau as usize..(ik + 1) * ntau as usize],
            qual,
            ceil,
            &tempfiles,
        );
        k += dk;
    }
    // There are no transforms when the fast dimension is zero. The Rust
    // FFT library asserts on this dimension; C's call does no work.
    if ntau != 0 {
        pfa2cr(1, 2, ntau, nxfft, &cq, &mut q);
    }
    tempfiles.check_signal();
    for ix in 0..nx {
        for itau in 0..ntau {
            q[(ix * ntau + itau) as usize] /= nxfft as f32;
        }
    }
    for ix in 0..nx as usize {
        tempfiles.check_signal();
        efread(std::slice::from_mut(&mut tr.hdr), &headerfp);
        tr.ns = ntau as u16;
        tr.dt = (dtau as f64 * 1000000.0) as u32 as u16;
        tr.delrt = (ftau as f64 * 1000.0) as i32 as i16;
        tr.data[..ntau as usize].copy_from_slice(&q[ix * ntau as usize..(ix + 1) * ntau as usize]);
        puttr(&mut tr);
    }
    tempfiles.check_signal();
    efclose(&headerfp);
    if !headerfile.is_empty() {
        eremove(&headerfile);
    }
    if istmpdir {
        eremove(&tracefile);
    }
    exit(if INTERRUPTED.load(Ordering::Relaxed) {
        EXIT_FAILURE
    } else {
        CWP_Exit()
    });
}

/// Gazdag's phase-shift zero-offset migration for one wavenumber,
/// adapted to a v(tau) velocity profile.
#[allow(clippy::too_many_arguments)]
fn gazdagvt(
    k: f32,
    nt: i32,
    dt: f32,
    ft: f32,
    ntau: i32,
    dtau: f32,
    ftau: f32,
    vt: &[f32],
    p: &[Complex],
    q: &mut [Complex],
    qual: f32,
    gainceil: f32,
    tempfiles: &TempFiles,
) {
    tempfiles.check_signal();
    let ntfft = npfa(nt);
    let nw = ntfft;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = (-PI / dt as f64) as f32;
    let tmax = ft + (nt - 1) as f32 * dt;
    let mut pp = vec![Complex::default(); nw as usize];
    let mut cumgain = vec![1.0f32; nw as usize];
    for it in 0..nt as usize {
        pp[it] = if it % 2 != 0 { cneg(p[it]) } else { p[it] };
    }
    for it in nt..ntfft {
        pp[it as usize] = cmplx(0.0, 0.0);
    }
    pfacc(1, &mut pp);

    // Preserve the C loop's comparison of an integer index with ftau.
    let mut itau = 0usize;
    while (itau as f32) < ftau {
        tempfiles.check_signal();
        let mut w = fw;
        for iw in 0..nw as usize {
            if w == 0.0 {
                w = (1e-10 / dt as f64) as f32;
            }
            let coss = (1.0 - (0.5 * vt[itau] as f64 * k as f64 / w as f64).powf(2.0)) as f32;
            if coss as f64 >= ((ftau / tmax) as f64).powf(2.0) {
                let phase = (w as f64 * (ft as f64 - ftau as f64 * (coss as f64).sqrt())) as f32;
                let cshift = cmplx((phase as f64).cos() as f32, (phase as f64).sin() as f32);
                pp[iw] = cmul(pp[iw], cshift);
            } else {
                pp[iw] = cmplx(0.0, 0.0);
            }
            w += dw;
        }
        itau += 1;
    }

    let mut tau = ftau;
    for itau in 0..ntau as usize {
        tempfiles.check_signal();
        q[itau] = cmplx(0.0, 0.0);
        let mut w = fw;
        for iw in 0..nw as usize {
            q[itau] = cadd(q[itau], pp[iw]);
            if w == 0.0 {
                w = (1e-10 / dt as f64) as f32;
            }
            let coss = (1.0 - (0.5 * vt[itau] as f64 * k as f64 / w as f64).powf(2.0)) as f32;
            if coss as f64 >= ((tau / tmax) as f64).powf(2.0) {
                let phase = ((-w * dtau) as f64 * (coss as f64).sqrt()) as f32;
                let cshift = cmplx((phase as f64).cos() as f32, (phase as f64).sin() as f32);
                if cumgain[iw] < gainceil {
                    let alpha = (w as f64 / (2.0 * vt[itau] as f64 * qual as f64)) as f32;
                    let gain = (0.5 * vt[itau] as f64 * dtau as f64 * alpha as f64)
                        .abs()
                        .exp() as f32;
                    pp[iw] = cmul(pp[iw], crmul(cshift, gain));
                    cumgain[iw] *= gain;
                } else {
                    pp[iw] = cmplx(0.0, 0.0);
                }
            } else {
                pp[iw] = cmplx(0.0, 0.0);
            }
            w += dw;
        }
        q[itau] = crmul(q[itau], (1.0 / nw as f64) as f32);
        tau += dtau;
    }
}
