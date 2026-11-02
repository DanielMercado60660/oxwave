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

//! SUMIGTK - migration via the T-K domain method for CMP stacked data.
//! Port of su/main/migration_inversion/sumigtk.c (Revision 1.17).
//! Credits (from the C source): CWP: Dave Hale.
//! mig1k: Dave Hale, Colorado School of Mines, 11/05/90.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cmplx, Complex};
use oxwave_cwp::cwp::{abs, max, min, PI};

use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::{npfaro, pfa2cr, pfa2rc};
use oxwave_sinc_interp::{ress8c, ress8r};
use oxwave_tempfile::{temporary_filename, temporary_stream};
use std::ffi::{c_char, CString};
use std::sync::atomic::{AtomicPtr, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGTK - MIGration via T-K domain method for common-midpoint stacked data",
    "\t\t\t\t\t\t\t\t\t",
    " sumigtk <stdin >stdout dxcdp= [optional parms]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " dxcdp                   distance between successive cdps\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " fmax=Nyquist            maximum frequency\t\t\t\t",
    " tmig=0.0                times corresponding to interval velocities in vmig",
    " vmig=1500.0             interval velocities corresponding to times in tmig",
    " vfile=                  binary (non-ascii) file containing velocities v(t)",
    " nxpad=0                 number of cdps to pad with zeros before FFT\t",
    " ltaper=0                length of linear taper for left and right edges",
    " verbose=0               =1 for diagnostic print\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input traces must be sorted by either increasing or decreasing cdp.\t",
    " \t\t\t\t\t\t\t\t\t",
    " The tmig and vmig arrays specify an interval velocity function of time.",
    " Linear interpolation and constant extrapolation is used to determine\t",
    " interval velocities at times not specified.  Values specified in tmig\t",
    " must increase monotonically.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Alternatively, interval velocities may be stored in a binary file\t",
    " containing one velocity for every time sample.  If vfile is specified,",
    " then the tmig and vmig arrays are ignored.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The time of first sample is assumed to be zero, regardless of the value",
    " of the trace header field delrt.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The migration is a reverse time migration in the (t,k) domain. In the\t",
    " first step, the data g(t,x) are Fourier transformed x->k into the\t",
    " the time-wavenumber domain g(t,k).\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Then looping over wavenumbers, the data are then reverse-time\t\t",
    " finite-difference migrated, wavenumber by wavenumber.  The resulting\t",
    " migrated data m(tau,k), now in the tau (migrated time) and k domain,\t",
    " are inverse fourier transformed back into m(tau,xout) and written out.",
    "\t\t\t\t\t\t\t\t\t",
];

/* Globals (so can trap signal) defining temporary disk files */
struct TemporaryFiles {
    tracefile: String,
    headerfile: String,
    tracefp: CFile,
    headerfp: CFile,
    istmpdir: bool,
}

static TEMPORARY_FILES: AtomicPtr<TemporaryFiles> = AtomicPtr::new(std::ptr::null_mut());

unsafe extern "C" {
    fn access(path: *const c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get information from the first header */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;

    /* get parameters */
    let mut dx = 0.0f32;
    if getparfloat("dxcdp", &mut dx) == 0 {
        err!("must specify dxcdp");
    }
    let mut fmax = 0.0f32;
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (0.5 / dt as f64) as f32;
    }
    let mut nxpad = 0i32;
    if getparint("nxpad", &mut nxpad) == 0 {
        nxpad = 0;
    }
    let mut ltaper = 0i32;
    if getparint("ltaper", &mut ltaper) == 0 {
        ltaper = 0;
    }
    let mut speed = 0.0f32;
    if getparfloat("speed", &mut speed) == 0 {
        speed = 1.0;
    }
    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* look for user-supplied tmpdir */
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    let path = CString::new(tmpdir.split('\0').next().unwrap()).unwrap();
    // WRITE_OK is W_OK (2): test this user's access, including path traversal.
    if !tmpdir.is_empty() && unsafe { access(path.as_ptr(), 2) } != 0 {
        err!("you can't write in %s (or it doesn't exist)", tmpdir);
    }

    /* determine velocity function v(t) */
    let mut vt = vec![0.0f32; nt as usize];
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
        for it in 0..nt as usize {
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
        // C passes an unchecked fopen to efread. A failed open terminates
        // with SIGSEGV; preserve its shell exit status without a C call.
        let fp = match fopen(&vfile, "r") {
            Some(fp) => fp,
            None => exit(139),
        };
        if efread(&mut vt, &fp) != nt as usize {
            err!("cannot read %d velocities from file %s", nt, vfile);
        }
    }
    checkpars();

    /* copy traces and headers to temporary files */
    let (tracefp, headerfp, tracefile, headerfile, istmpdir);
    let temporary_files;
    if tmpdir.is_empty() {
        tracefp = etmpfile();
        headerfp = etmpfile();
        tracefile = String::new();
        headerfile = String::new();
        istmpdir = false;
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        let directory = &tmpdir;
        tracefile = temporary_filename(directory);
        headerfile = temporary_filename(directory);
        if tracefile == headerfile {
            // The shared naming routine returns the same fallback name when
            // /tmp is unavailable (e.g. in a filesystem sandbox). Use its
            // directory-local streams to keep traces and headers separate.
            tracefp = CFile::from_file(temporary_stream(Some(directory)).unwrap_or_else(|| {
                err!("%s: efopen: fopen failed", "subcalls.c");
            }));
            headerfp = CFile::from_file(temporary_stream(Some(directory)).unwrap_or_else(|| {
                err!("%s: efopen: fopen failed", "subcalls.c");
            }));
            istmpdir = false; // temporary_stream unlinks its files immediately.
        } else {
            tracefp = efopen(&tracefile, "w+");
            headerfp = efopen(&headerfile, "w+");
            istmpdir = true;
        }
        temporary_files = Box::new(TemporaryFiles {
            tracefile: tracefile.clone(),
            headerfile: headerfile.clone(),
            tracefp: tracefp.clone(),
            headerfp: headerfp.clone(),
            istmpdir,
        });
        TEMPORARY_FILES.store(
            std::ptr::from_ref(temporary_files.as_ref()).cast_mut(),
            Ordering::Relaxed,
        );
        // Trap signals after opening both streams. C's handler during a
        // partially completed fopen sequence could fclose a NULL stream.
        unsafe {
            signal(2, closefiles); // SIGINT
            signal(3, closefiles); // SIGQUIT
            signal(1, closefiles); // SIGHUP
            signal(15, closefiles); // SIGTERM
        }
        if verbose != 0 {
            warn!("putting temporary files in %s", directory);
        }
    }

    /* count the traces */
    let mut nx = 0i32;
    loop {
        nx += 1;
        efwrite(&[tr.hdr], &headerfp);
        efwrite(&tr.data[..nt as usize], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    erewind(&headerfp);
    erewind(&tracefp);
    if verbose != 0 {
        fprintf!(stderr(), "\t%d traces input\n", nx);
    }

    /* determine wavenumber sampling */
    let nxfft = npfaro(
        nx.wrapping_add(nxpad),
        2i32.wrapping_mul(nx.wrapping_add(nxpad)),
    );
    let nk = nxfft / 2 + 1;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;

    /* allocate space for Fourier transform */
    let mut gtk = vec![Complex::default(); (nt * nk) as usize];
    // C aliases gtk as real storage. Separate buffers fit the library's safe
    // API while preserving the same x-major/time-minor array layout.
    let mut gtx = vec![0.0f32; (nt * nxfft) as usize];

    /* read and apply fft scaling to traces and pad with zeros */
    let fftscl = (1.0 / nxfft as f64) as f32;
    for ix in 0..nx {
        let row = &mut gtx[(ix * nt) as usize..((ix + 1) * nt) as usize];
        efread(row, &tracefp);
        for it in 0..nt as usize {
            row[it] *= fftscl;
        }
        if ix < ltaper {
            let taper = (ix + 1) as f32 / (ltaper + 1) as f32;
            for it in 0..nt as usize {
                row[it] *= taper;
            }
        } else if ix >= nx - ltaper {
            let taper = (nx - ix) as f32 / (ltaper + 1) as f32;
            for it in 0..nt as usize {
                row[it] *= taper;
            }
        }
    }
    for ix in nx..nxfft {
        for it in 0..nt {
            gtx[(ix * nt + it) as usize] = 0.0;
        }
    }

    /* Fourier transform g(t,x) to g(t,k) */
    pfa2rc(-1, 2, nt, nxfft, &gtx, &mut gtk);
    if verbose != 0 {
        fprintf!(stderr(), "\tFourier transform done\n");
    }

    /* loop over wavenumbers */
    let mut k = 0.0f32;
    for ik in 0..nk {
        if verbose != 0 && ik % (if nk / 10 > 0 { nk / 10 } else { 1 }) == 0 {
            fprintf!(stderr(), "\t%d of %d wavenumbers done\n", ik, nk);
        }
        mig1k(
            k,
            fmax,
            speed,
            nt,
            dt,
            &vt,
            &mut gtk[(ik * nt) as usize..((ik + 1) * nt) as usize],
        );
        k += dk;
    }

    /* Fourier transform g(t,k) to g(t,x) */
    pfa2cr(1, 2, nt, nxfft, &gtk, &mut gtx);
    if verbose != 0 {
        fprintf!(stderr(), "\tinverse Fourier transform done\n");
    }

    /* output migrated traces with headers */
    for ix in 0..nx {
        efread(std::slice::from_mut(&mut tr.hdr), &headerfp);
        tr.data[..nt as usize].copy_from_slice(&gtx[(ix * nt) as usize..((ix + 1) * nt) as usize]);
        puttr(&mut tr);
    }

    /* clean up */
    efclose(&headerfp);
    if istmpdir {
        eremove(&headerfile);
    }
    efclose(&tracefp);
    if istmpdir {
        eremove(&tracefile);
    }
    exit(CWP_Exit());
}

/// Migration in t-k for one wavenumber. C calls with p == q; pq retains that
/// in-place contract while gs holds the resampled input.
fn mig1k(k: f32, fmax: f32, speed: f32, nt: i32, dt: f32, v: &[f32], pq: &mut [Complex]) {
    let nfac = 5i32;
    let fac = 0.75f32;
    let czero = cmplx(0.0, 0.0);

    /* determine time sampling to avoid excessive grid dispersion */
    let mut vmin = v[0];
    let mut vmax = v[0];
    for it in 1..nt as usize {
        if v[it] < vmin {
            vmin = v[it];
        }
        if v[it] > vmax {
            vmax = v[it];
        }
    }
    let mut ds;
    if k != 0.0 {
        ds = (fac as f64
            * max(
                0.4 * PI / abs(vmax * k) as f64,
                0.1 * vmin as f64 / (vmax * fmax) as f64,
            )) as f32;
        ds *= speed;
        ds = min(ds, dt);
    } else {
        ds = dt;
    }
    let mut ns = (1.0f32 + (nt - 1) as f32 * dt / ds) as i32;
    ns = min(ns, nfac * nt);
    ds = (nt - 1) as f32 * dt / (ns - 1) as f32;
    fprintf!(stderr(), "ns=%d\n", ns);

    /* allocate workspace */
    // Negative ns converts to a huge size_t in C and malloc fails. Avoid
    // Rust's capacity-overflow panic and preserve ealloc.c's diagnostic.
    if ns < 0 {
        err!(
            "%s: malloc failed (%s)",
            "ealloc.c",
            oxwave_core::err::strerror(12)
        );
    }
    let mut vs = vec![0.0f32; ns as usize];
    let mut cs = vec![0.0f32; (ns * 2) as usize];
    let mut gs = vec![Complex::default(); ns as usize];

    /* resample v(t) and p(t) */
    ress8r(dt, 0.0, v, v[0], v[nt as usize - 1], ds, 0.0, &mut vs);
    ress8c(dt, 0.0, pq, czero, czero, ds, 0.0, &mut gs);

    /* compute finite-difference coefficients */
    for is in 0..ns as usize {
        let mut temp = (0.125 * vs[is] as f64 * k as f64 * ds as f64) as f32;
        temp = temp * temp;
        temp = ((1.0 - temp as f64) / (1.0 + temp as f64)) as f32;
        cs[2 * is] = temp;
        cs[2 * is + 1] = temp;
    }

    /* loop over t2 = (tau-t)/sqrt(2) */
    for i2 in 2 - ns..ns {
        let i1stop = if i2 <= 0 { 1 - i2 } else { i2 };
        let mut g00r = 0.0f32;
        let mut g00i = 0.0f32;
        let mut g01r = 0.0f32;
        let mut g01i = 0.0f32;
        let mut gi = ns - 1;
        let mut cp = i2 + ns - 2;
        let mut ct = cs[cp as usize];
        cp -= 1;

        /* loop over t1 = (tau+t)/sqrt(2) */
        for _i1 in (i1stop..=ns - 1).rev() {
            /* update real part */
            let g10r = g00r;
            let g11r = g01r;
            g00r = gs[gi as usize].r;
            g01r = ct * (g11r + g00r) - g10r;
            gs[gi as usize].r = g01r;

            /* update imaginary part */
            let g10i = g00i;
            let g11i = g01i;
            g00i = gs[gi as usize].i;
            g01i = ct * (g11i + g00i) - g10i;
            gs[gi as usize].i = g01i;

            /* update pointers and finite-difference coefficient */
            gi -= 1;
            // C's last, unused coefficient fetch can be cs[-1]. Load only
            // when another star update will actually consume the value.
            if _i1 > i1stop {
                ct = cs[cp as usize];
            }
            cp -= 1;
        }
    }

    /* resample q(t) */
    ress8c(ds, 0.0, &gs, czero, czero, dt, 0.0, pq);
}

/* for graceful interrupt termination */
extern "C" fn closefiles(_sig: i32) {
    // Installed only once both streams exist; main retains the Box until exit.
    // As in C, the handler closes stdio streams before removing their names.
    let files = unsafe { &*TEMPORARY_FILES.load(Ordering::Relaxed) };
    efclose(&files.headerfp);
    efclose(&files.tracefp);
    if files.istmpdir {
        eremove(&files.headerfile);
        eremove(&files.tracefile);
    }
    exit(EXIT_FAILURE);
}
