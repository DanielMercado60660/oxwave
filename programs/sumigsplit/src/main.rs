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

//! SUMIGSPLIT - Split-step depth migration for zero-offset data.
//! Port of su/main/migration_inversion/sumigsplit.c (Revision 1.9).
//! Credits (from the C source): CWP Baoniu Han, July 21th, 1997.
//! Reference: Stoffa, Fokkema, Freire and Kessinger (1990), Geophysics 55, 410-421.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::PI;

use oxwave_pfa_fft::{npfa, npfar, pfa2cc, pfa2rc};
use oxwave_tempfile::temporary_filename;
use std::sync::atomic::{AtomicPtr, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUMIGSPLIT - Split-step depth migration for zero-offset data.         ",
    "                                                                       ",
    " sumigsplit <infile >outfile vfile= [optional parameters]              ",
    "                                                                       ",
    " Required Parameters:                                                  ",
    " nz=                   number of depth samples                         ",
    " dz=                   depth sampling interval                         ",
    " vfile=                name of file containing velocities              ",
    "                                                                       ",
    " Optional Parameters:                                                  ",
    " dt=from header(dt) or .004    time sampling interval                  ",
    " dx=from header(d2) or 1.0     midpoint sampling interval              ",
    " ft=0.0                        first time sample                       ",
    " fz=                           first depth sample                      ",
    "                                                                       ",
    " tmpdir=        if non-empty, use the value as a directory path        ",
    "                prefix for storing temporary files; else if the        ",
    "                the CWP_TMPDIR environment variable is set use         ",
    "                its value for the path; else use tmpfile()             ",
    "                                                                       ",
    "                                                                       ",
    " Notes:                                                                ",
    " The input velocity file 'vfile' consists of C-style binary floats.  ",
    " The structure of this file is vfile[iz][ix]. Note that this means that",
    " the x-direction is the fastest direction instead of z-direction! Such a",
    " structure is more convenient for the downward continuation type       ",
    " migration algorithm than using z as fastest dimension as in other SU  ",
    " programs.                                                     \t",
    "                                                                       ",
    " Because most of the tools in the SU package (such as  unif2, unisam2, ",
    " and makevel) produce output with the structure vfile[ix][iz], you will",
    " need to transpose the velocity files created by these programs. You may",
    " use the SU program 'transp' in SU to transpose such files into the  ",
    " required vfile[iz][ix] structure.                                     ",
    " (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  \t\t",
    " denotes a v(z,x) array, the opposite of what Matlab and Fortran\t",
    " programmers may expect.)\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// The C globals used by closefiles. Keep the streams and names in one stable
// allocation owned by main until exit, including after tracefp is closed.
struct TemporaryFiles {
    tracefile: String,
    headerfile: String,
    tracefp: CFile,
    headerfp: CFile,
}

static TEMPORARY_FILES: AtomicPtr<TemporaryFiles> = AtomicPtr::new(std::ptr::null_mut());

unsafe extern "C" {
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get info from first trace */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;

    /* let user give dt and/or dx from command line */
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

    /* get optional parameters */
    let mut ft = 0.0f32;
    if getparfloat("ft", &mut ft) == 0 {
        ft = 0.0;
    }
    let mut nz = 0i32;
    if getparint("nz", &mut nz) == 0 {
        err!("nz must be specified");
    }
    let mut dz = 0.0f32;
    if getparfloat("dz", &mut dz) == 0 {
        err!("dz must be specified");
    }
    let mut fz = 0.0f32;
    if getparfloat("fz", &mut fz) == 0 {
        fz = 0.0;
    }
    let mut vfile = String::new();
    if getparstring("vfile", &mut vfile) == 0 {
        err!("vfile must be specified");
    }
    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    // ft is parsed but never used by the original.
    let _ = ft;

    /* Look for user-supplied tmpdir */
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        // Call the same POSIX access routine as C: mode bits alone do not
        // account for the real user/group IDs, ACLs, or read-only mounts.
        unsafe extern "C" {
            fn access(path: *const std::ffi::c_char, mode: std::ffi::c_int) -> std::ffi::c_int;
        }
        const WRITE_OK: std::ffi::c_int = 2;
        let path = std::ffi::CString::new(tmpdir.split('\0').next().unwrap()).unwrap();
        // SAFETY: path is NUL-terminated and remains alive for this call.
        if unsafe { access(path.as_ptr(), WRITE_OK) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }

    /* store traces and headers in tempfiles while getting a count */
    let tracefp;
    let headerfp;
    let temporary_files;
    if tmpdir.is_empty() {
        temporary_files = None;
        tracefp = etmpfile();
        headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        let directory = tmpdir.clone();
        let tracefile = temporary_filename(&directory);
        let headerfile = temporary_filename(&directory);
        tracefp = efopen(&tracefile, "w+");
        headerfp = efopen(&headerfile, "w+");
        let files = Box::new(TemporaryFiles {
            tracefile,
            headerfile,
            tracefp: tracefp.clone(),
            headerfp: headerfp.clone(),
        });
        TEMPORARY_FILES.store(
            std::ptr::from_ref(files.as_ref()).cast_mut(),
            Ordering::Relaxed,
        );
        temporary_files = Some(files);
        // Trap signals after the streams are initialized. The C handler's
        // fclose(NULL) during fopen would be undefined. These signal numbers
        // are the same on macOS and Linux.
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

    checkpars();

    let mut nx = 0i32;
    loop {
        nx += 1;
        efwrite(std::slice::from_ref(&tr.hdr), &headerfp);
        efwrite(&tr.data[..nt as usize], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    erewind(&tracefp);
    erewind(&headerfp);

    /* determine frequency sampling (for real to complex FFT) */
    let ntfft = npfar(nt);
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = 0.0f32;

    /* determine wavenumber sampling (for complex to complex FFT) */
    let nxfft = npfa(nx);
    let nk = nxfft;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let fk = (-PI / dx as f64) as f32;

    /* allocate space: flat arrays preserve the C contiguous layout */
    let mut p = vec![0.0f32; ntfft as usize * nx as usize];
    let mut cp = vec![Complex::default(); nw as usize * nx as usize];
    let mut cq = vec![Complex::default(); nw as usize * nk as usize];
    let mut cq1 = vec![Complex::default(); nw as usize * nk as usize];
    let mut cresult = vec![0.0f32; nz as usize * nx as usize];
    let mut v = vec![0.0f32; nx as usize * nz as usize];

    /* load traces into the zero-offset array and close tmpfile */
    for ix in 0..nx as usize {
        efread(
            &mut p[ix * ntfft as usize..ix * ntfft as usize + nt as usize],
            &tracefp,
        );
    }
    efclose(&tracefp);

    /* load velocity file */
    let vfp = efopen(&vfile, "r");
    efread(&mut v, &vfp);
    efclose(&vfp);

    /* pad with zeros and Fourier transform x to k */
    for ix in 0..nx as usize {
        for it in nt as usize..ntfft as usize {
            p[ix * ntfft as usize + it] = 0.0;
        }
    }
    pfa2rc(1, 1, ntfft, nx, &p, &mut cp);

    /* loops over depth */
    let mut tz = fz;
    for iz in 0..nz as usize {
        for ix in 0..nx as usize {
            cresult[ix * nz as usize + iz] = 0.0;
            for iw in 0..nw as usize {
                cresult[ix * nz as usize + iz] += cp[ix * nw as usize + iw].r / ntfft as f32;
            }
            for iw in 0..nw as usize {
                cp[ix * nw as usize + iw].r -= cresult[ix * nz as usize + iz] / ntfft as f32;
            }
        }

        let mut vmin = 0.0f32;
        for ix in 0..nx as usize {
            vmin = (vmin as f64 + 1.0 / v[iz * nx as usize + ix] as f64 / nx as f64) as f32;
        }
        vmin = (1.0 / vmin as f64 * 0.5) as f32;

        for ik in 0..nx as usize {
            for iw in 0..nw as usize {
                cq[ik * nw as usize + iw] = if ik % 2 != 0 {
                    cneg(cp[ik * nw as usize + iw])
                } else {
                    cp[ik * nw as usize + iw]
                };
            }
        }
        for ik in nx as usize..nk as usize {
            for iw in 0..nw as usize {
                cq[ik * nw as usize + iw] = cmplx(0.0, 0.0);
            }
        }

        /* FFT to W-K domain */
        pfa2cc(-1, 2, nw, nk, &mut cq);
        let v1 = vmin;
        let mut k = fk;
        for ik in 0..nk as usize {
            let mut w = fw;
            for iw in 0..nw as usize {
                if w == 0.0 {
                    w = (1.0e-10 / dt as f64) as f32;
                }
                let kz1 = 1.0 - ((v1 * k / w) as f64).powf(2.0);
                if kz1 > 0.0 {
                    let phase1 = -w as f64 * kz1.sqrt() * dz as f64 / v1 as f64;
                    let cshift1 = cmplx(phase1.cos() as f32, phase1.sin() as f32);
                    cq1[ik * nw as usize + iw] = cmul(cq[ik * nw as usize + iw], cshift1);
                } else {
                    let cshift1 = cmplx(0.0, 0.0);
                    cq1[ik * nw as usize + iw] = cmul(cq[ik * nw as usize + iw], cshift1);
                }
                w += dw;
            }
            k += dk;
        }

        pfa2cc(1, 2, nw, nk, &mut cq1);
        for ix in 0..nx as usize {
            let mut w = fw;
            for iw in 0..nw as usize {
                let index = ix * nw as usize + iw;
                cq1[index] = crmul(cq1[index], (1.0 / nxfft as f64) as f32);
                cq1[index] = if ix % 2 != 0 {
                    cneg(cq1[index])
                } else {
                    cq1[index]
                };
                let kz2 = (1.0 / v1 as f64 - 2.0 / v[iz * nx as usize + ix] as f64)
                    * w as f64
                    * dz as f64;
                let cshift2 = cmplx(kz2.cos() as f32, kz2.sin() as f32);
                cp[index] = cmul(cq1[index], cshift2);
                w += dw;
            }
        }
        tz += dz;
    }
    // tz is advanced but never used by the original.
    let _ = tz;

    /* restore header fields and write output */
    for ix in 0..nx as usize {
        efread(std::slice::from_mut(&mut tr.hdr), &headerfp);
        tr.ns = nz as u16;
        tr.d1 = dz;
        tr.data[..nz as usize].copy_from_slice(&cresult[ix * nz as usize..(ix + 1) * nz as usize]);
        puttr(&mut tr);
    }

    /* Clean up */
    efclose(&headerfp);
    if let Some(files) = &temporary_files {
        eremove(&files.headerfile);
        eremove(&files.tracefile);
    }
    exit(CWP_Exit());
}

/* for graceful interrupt termination */
extern "C" fn closefiles(_sig: i32) {
    // The handler is installed only after the Box is initialized, and main
    // keeps it alive until exit. Like the C helper, this uses stdio in the
    // handler; signals interrupting a temp-stream operation are not reentrant.
    let files = unsafe { &*TEMPORARY_FILES.load(Ordering::Relaxed) };
    efclose(&files.headerfp);
    efclose(&files.tracefp);
    eremove(&files.headerfile);
    eremove(&files.tracefile);
    exit(EXIT_FAILURE);
}
