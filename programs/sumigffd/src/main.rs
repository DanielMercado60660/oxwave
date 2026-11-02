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

//! SUMIGFFD - Fourier finite difference depth migration for zero-offset data.
//! Port of su/main/migration_inversion/sumigffd.c (Revision 1.9).
//! Credits: CWP Baoniu Han, July 21th, 1997.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::PI;
use oxwave_pfa_fft::{npfa, npfar, pfa2cc, pfa2rc};
use oxwave_tempfile::temporary_filename;
use std::ffi::{c_char, CString};
use std::sync::atomic::{AtomicPtr, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGFFD - Fourier finite difference depth migration for\t\t",
    "\t    zero-offset data. This method is a hybrid migration which\t",
    "\t    combines the advantages of phase shift and finite difference",
    "\t    migrations.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sumigffd <infile >outfile vfile= [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t  \t",
    " nz=\t\t   number of depth samples\t\t\t \t",
    " dz=\t\t   depth sampling interval\t\t\t \t",
    " vfile=\t\tname of file containing velocities\t      \t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t  \t",
    " dt=from header(dt) or .004    time sampling interval\t\t  \t",
    " dx=from header(d2) or 1.0     midpoint sampling interval\t  \t",
    " ft=0.0\t\t\tfirst time sample\t\t\t",
    " fz=0.0\t\t\tfirst depth sample\t\t      \t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir=\tif non-empty, use the value as a directory path\t\t",
    "\t\tprefix for storing temporary files; else if the\t\t",
    "\t\tthe CWP_TMPDIR environment variable is set use\t\t",
    "\t\tits value for the path; else use tmpfile()\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The input velocity file 'vfile' consists of C-style binary floats.  ",
    " The structure of this file is vfile[iz][ix]. Note that this means that",
    " the x-direction is the fastest direction instead of z-direction! Such a",
    " structure is more convenient for the downward continuation type\t",
    " migration algorithm than using z as fastest dimension as in other SU  ",
    " programs. (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  ",
    " denotes a v(z,x) array, the opposite of what Matlab and Fortran\t",
    " programmers may expect.)\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Because most of the tools in the SU package (such as  unif2, unisam2, ",
    " and makevel) produce output with the structure vfile[ix][iz], you will",
    " need to transpose the velocity files created by these programs. You may",
    " use the SU program 'transp' in SU to transpose such files into the  ",
    " required vfile[iz][ix] structure.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// Immutable C strings remain alive in main until exit. The callback only
// reads their pointers and calls async-signal-safe Unix functions; CFile
// cannot be borrowed while a signal has interrupted one of its operations.
static TRACEFILE: AtomicPtr<c_char> = AtomicPtr::new(std::ptr::null_mut());
static HEADERFILE: AtomicPtr<c_char> = AtomicPtr::new(std::ptr::null_mut());

unsafe extern "C" {
    fn access(path: *const c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn unlink(path: *const c_char) -> i32;
    fn _exit(status: i32) -> !;
}

extern "C" fn closefiles(_sig: i32) {
    // C closes both storage streams, removes their names, and exits 1.
    // Unlink the names immediately, even when gettr is blocked on a pipe.
    // _exit closes the descriptors without reentering interrupted CFile
    // buffers, whose temporary contents are discarded on this exit path.
    unsafe {
        unlink(HEADERFILE.load(Ordering::Acquire));
        unlink(TRACEFILE.load(Ordering::Acquire));
        _exit(EXIT_FAILURE);
    }
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
    let nt = tr.ns as usize;

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
    let mut verbose = 1i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    let path = CString::new(tmpdir.split('\0').next().unwrap()).unwrap();
    // WRITE_OK is W_OK (2) on Unix; access also tests directory traversal.
    if !tmpdir.is_empty() && unsafe { access(path.as_ptr(), 2) } != 0 {
        err!("you can't write in %s (or it doesn't exist)", tmpdir);
    }
    checkpars();

    let (tracefp, headerfp);
    let mut tracefile = String::new();
    let mut headerfile = String::new();
    let istmpdir = !tmpdir.is_empty();
    let (tracepath, headerpath);
    if !istmpdir {
        tracefp = etmpfile();
        headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        tracefile = temporary_filename(&tmpdir);
        headerfile = temporary_filename(&tmpdir);
        // Keep the two storage files distinct if the filename library
        // supplies identical fallback names after a creation failure.
        if tracefile == headerfile {
            headerfile.push_str(".headers");
        }
        tracepath = CString::new(tracefile.as_bytes()).unwrap();
        headerpath = CString::new(headerfile.as_bytes()).unwrap();
        TRACEFILE.store(tracepath.as_ptr().cast_mut(), Ordering::Release);
        HEADERFILE.store(headerpath.as_ptr().cast_mut(), Ordering::Release);
        unsafe {
            signal(2, closefiles); // SIGINT
            signal(3, closefiles); // SIGQUIT
            signal(1, closefiles); // SIGHUP
            signal(15, closefiles); // SIGTERM
        }
        tracefp = efopen(&tracefile, "w+");
        headerfp = efopen(&headerfile, "w+");
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
    }
    let mut nx = 0usize;
    loop {
        nx += 1;
        efwrite(std::slice::from_ref(&tr.hdr), &headerfp);
        efwrite(&tr.data[..nt], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    erewind(&tracefp);
    erewind(&headerfp);

    /* determine frequency sampling (for real to complex FFT) */
    let ntfft = npfar(nt as i32) as usize;
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = 0.0f32;
    let nxfft = npfa(nx as i32) as usize;
    let nk = nxfft;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let fk = (-PI / dx as f64) as f32;

    /* allocate space */
    let nz = nz as usize;
    let mut p = vec![0.0f32; ntfft * nx];
    let mut cp = vec![Complex::default(); nw * nx];
    let mut cq = vec![Complex::default(); nw * nk];
    let mut cresult = vec![vec![0.0f32; nz]; nx];
    let mut v = vec![0.0f32; nx * nz];
    for ix in 0..nx {
        efread(&mut p[ix * ntfft..ix * ntfft + nt], &tracefp);
    }
    efclose(&tracefp);
    let vfp = efopen(&vfile, "r");
    efread(&mut v, &vfp);
    efclose(&vfp);
    for ix in 0..nx {
        for it in nt..ntfft {
            p[ix * ntfft + it] = 0.0;
        }
    }
    pfa2rc(1, 1, ntfft as i32, nx as i32, &p, &mut cp);

    /* loop over depth */
    let mut tz = fz;
    for iz in 0..nz {
        for ix in 0..nx {
            cresult[ix][iz] = 0.0;
            for iw in 1..nw {
                cresult[ix][iz] += cp[ix * nw + iw].r / ntfft as f32;
            }
        }
        let mut vmin = v[iz * nx];
        let mut vmax = vmin;
        for ix in 0..nx {
            if v[iz * nx + ix] >= vmax {
                vmax = v[iz * nx + ix];
            }
            if v[iz * nx + ix] <= vmin {
                vmin = v[iz * nx + ix];
            }
        }
        for ik in 0..nx {
            for iw in 0..nw {
                cq[ik * nw + iw] = if ik % 2 != 0 {
                    cneg(cp[ik * nw + iw])
                } else {
                    cp[ik * nw + iw]
                };
            }
        }
        for ik in nx..nk {
            for iw in 0..nw {
                cq[ik * nw + iw] = cmplx(0.0, 0.0);
            }
        }
        pfa2cc(-1, 2, nw as i32, nk as i32, &mut cq);
        let v1 = (vmin as f64 * 0.5) as f32;
        let mut k = fk;
        for ik in 0..nk {
            let mut w = fw + dw;
            for iw in 1..nw {
                if w == 0.0 {
                    w = (1.0e-10 / dt as f64) as f32;
                }
                let kz1 = 1.0 - ((v1 * k / w) as f64).powf(2.0);
                let cshift1 = if kz1 > 0.0 {
                    let phase1 = -(w as f64) * kz1.sqrt() * dz as f64 / v1 as f64;
                    cmplx(phase1.cos() as f32, phase1.sin() as f32)
                } else {
                    cmplx(0.0, 0.0)
                };
                cq[ik * nw + iw] = cmul(cq[ik * nw + iw], cshift1);
                w += dw;
            }
            k += dk;
        }
        pfa2cc(1, 2, nw as i32, nk as i32, &mut cq);
        for ix in 0..nx {
            let mut w = fw + dw;
            for iw in 1..nw {
                cq[ix * nw + iw] = crmul(cq[ix * nw + iw], (1.0 / nxfft as f64) as f32);
                cp[ix * nw + iw] = if ix % 2 != 0 {
                    cneg(cq[ix * nw + iw])
                } else {
                    cq[ix * nw + iw]
                };
                w += dw;
            }
            let _ = w; // The C loop increments w but does not use it here.
        }
        for ix in 0..nx {
            let mut w = fw + dw;
            for iw in 1..nw {
                let kz2 = (1.0 / v1 as f64 - 2.0 / v[iz * nx + ix] as f64) * w as f64 * dz as f64;
                let cshift2 = cmplx(kz2.cos() as f32, kz2.sin() as f32);
                cp[ix * nw + iw] = cmul(cp[ix * nw + iw], cshift2);
                w += dw;
            }
        }
        fdmig(
            &mut cp,
            nx,
            nw,
            &v[iz * nx..(iz + 1) * nx],
            fw,
            dw,
            dz,
            dx,
            dt,
            v1,
        );
        tz += dz;
    }
    // ft is read but unused; tz (from fz) and vmax are also unused in C.
    let _ = (ft, tz);

    /* restore header fields and write output */
    for ix in 0..nx {
        efread(std::slice::from_mut(&mut tr.hdr), &headerfp);
        tr.ns = nz as u16;
        tr.d1 = dz;
        tr.data[..nz].copy_from_slice(&cresult[ix]);
        puttr(&mut tr);
    }
    efclose(&headerfp);
    if istmpdir {
        eremove(&headerfile);
        eremove(&tracefile);
    }
    exit(CWP_Exit());
}

#[allow(clippy::too_many_arguments)]
fn fdmig(
    cp: &mut [Complex],
    nx: usize,
    nw: usize,
    v: &[f32],
    fw: f32,
    dw: f32,
    dz: f32,
    dx: f32,
    dt: f32,
    vc: f32,
) {
    let a0 = 2.0f32;
    let trick = 0.14f32;
    let mut p = vec![0.0f32; nx];
    let mut s1 = vec![0.0f32; nx];
    let mut s2 = vec![0.0f32; nx];
    // C reads data[1], data[2] and data[nx-3] even with fewer than
    // three traces. Keep zeroed guard slots for those allocator reads.
    let mut data = vec![Complex::default(); nx.max(3)];
    let mut d = vec![Complex::default(); nx];
    let mut a = vec![Complex::default(); nx];
    let mut b = vec![Complex::default(); nx];
    let mut c = vec![Complex::default(); nx];
    for ix in 0..nx {
        p[ix] = (vc as f64 * 2.0 / v[ix] as f64) as f32;
        p[ix] = (0.5 * ((p[ix] * p[ix] + p[ix]) as f64 + 1.0)) as f32;
    }
    let v1 = v[0];
    let vn = v[nx - 1];
    let mut w = fw + dw;
    for iw in 1..nw {
        if w == 0.0 {
            w = (1.0e-10 / dt as f64) as f32;
        }
        for ix in 0..nx {
            s1[ix] = ((p[ix] * v[ix] * v[ix]) as f64
                / (4.0 * a0 as f64 * dx as f64 * dx as f64 * w as f64 * w as f64)
                + trick as f64) as f32;
            s2[ix] = ((1.0 - vc as f64 * 2.0 / v[ix] as f64) * v[ix] as f64 * dz as f64 * 0.5
                / (2.0 * w as f64 * a0 as f64 * dx as f64 * dx as f64)) as f32;
        }
        for ix in 0..nx {
            data[ix] = cp[ix * nw + iw];
        }
        let cp2 = data[1];
        let cp3 = data[2];
        let cpnm1 = if nx >= 2 {
            data[nx - 2]
        } else {
            Complex::default()
        };
        let cpnm2 = if nx >= 3 {
            data[nx - 3]
        } else {
            Complex::default()
        };
        let mut a1 = crmul(cmul(cp2, conjg(cp3)), 2.0);
        let b1 = cadd(cmul(cp2, conjg(cp2)), cmul(cp3, conjg(cp3)));
        if b1.r == 0.0 && b1.i == 0.0 {
            a1 = cwp_cexp(cmplx(0.0, (-((w * dx) as f64) * 0.5 / v1 as f64) as f32));
        } else {
            a1 = cdiv(a1, b1);
        }
        if a1.i > 0.0 {
            a1 = cwp_cexp(cmplx(0.0, (-((w * dx) as f64) * 0.5 / v1 as f64) as f32));
        }
        let mut a2 = crmul(cmul(cpnm1, conjg(cpnm2)), 2.0);
        let b2 = cadd(cmul(cpnm1, conjg(cpnm1)), cmul(cpnm2, conjg(cpnm2)));
        if b2.r == 0.0 && b2.i == 0.0 {
            a2 = cwp_cexp(cmplx(0.0, (-((w * dx) as f64) * 0.5 / vn as f64) as f32));
        } else {
            a2 = cdiv(a2, b2);
        }
        if a2.i > 0.0 {
            a2 = cwp_cexp(cmplx(0.0, (-((w * dx) as f64) * 0.5 / vn as f64) as f32));
        }
        for ix in 0..nx {
            a[ix] = cmplx(s1[ix], -s2[ix]);
            b[ix] = cmplx(
                (1.0 - 2.0 * s1[ix] as f64) as f32,
                (2.0 * s2[ix] as f64) as f32,
            );
        }
        d[0] = cadd(
            cmul(cadd(b[0], cmul(a[0], a1)), data[0]),
            // C indexes a[1] out of bounds for nx=1; keep the explicit
            // zero boundary without evaluating that invalid array access.
            cmul(data[1], if nx > 1 { a[1] } else { Complex::default() }),
        );
        d[nx - 1] = cadd(
            cmul(cadd(b[nx - 1], cmul(a[nx - 1], a2)), data[nx - 1]),
            if nx > 1 {
                cmul(data[nx - 2], a[nx - 2])
            } else {
                Complex::default()
            },
        );
        for ix in 1..nx - 1 {
            d[ix] = cadd(
                cadd(cmul(data[ix + 1], a[ix + 1]), cmul(data[ix - 1], a[ix - 1])),
                cmul(data[ix], b[ix]),
            );
        }
        for ix in 0..nx {
            data[ix] = cmplx(s1[ix], s2[ix]);
            b[ix] = cmplx(
                (1.0 - 2.0 * s1[ix] as f64) as f32,
                (-2.0 * s2[ix] as f64) as f32,
            );
        }
        let endl = cadd(b[0], cmul(data[0], a1));
        let endr = cadd(b[nx - 1], cmul(data[nx - 1], a2));
        for ix in 1..nx - 1 {
            a[ix] = data[ix + 1];
            c[ix] = data[ix - 1];
        }
        a[0] = data[1];
        c[nx - 1] = if nx > 1 {
            data[nx - 2]
        } else {
            Complex::default()
        };
        retris(&mut data, &a, &c, &b, endl, endr, nx, &d);
        for ix in 0..nx {
            cp[ix * nw + iw] = data[ix];
        }
        w += dw;
    }
}

#[allow(clippy::too_many_arguments)]
fn retris(
    data: &mut [Complex],
    a: &[Complex],
    c: &[Complex],
    b: &[Complex],
    endl: Complex,
    endr: Complex,
    nx: usize,
    d: &[Complex],
) {
    let mut e = vec![Complex::default(); nx];
    let mut f = vec![Complex::default(); nx];
    e[0] = cdiv(cneg(a[0]), endl);
    f[0] = cdiv(d[0], endl);
    for ix in 1..nx - 1 {
        let den = cadd(b[ix], cmul(c[ix], e[ix - 1]));
        e[ix] = cdiv(cneg(a[ix]), den);
        f[ix] = cdiv(csub(d[ix], cmul(f[ix - 1], c[ix])), den);
    }
    let fprev = if nx >= 2 {
        f[nx - 2]
    } else {
        Complex::default()
    };
    let eprev = if nx >= 2 {
        e[nx - 2]
    } else {
        Complex::default()
    };
    data[nx - 1] = cdiv(
        csub(d[nx - 1], cmul(fprev, c[nx - 1])),
        cadd(endr, cmul(c[nx - 1], eprev)),
    );
    for ix in (0..nx - 1).rev() {
        data[ix] = cadd(cmul(data[ix + 1], e[ix]), f[ix]);
    }
}
