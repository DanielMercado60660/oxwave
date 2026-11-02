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

//! SUMIGPSPI - Gazdag's phase-shift plus interpolation depth migration.
//! Port of su/main/migration_inversion/sumigpspi.c (Revision 1.12).
//! Credits: CWP, Baoniu Han, April 20th, 1998.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cmplx, cmul, cneg, crmul, cwp_cexp, Complex};
use oxwave_cwp::cwp::PI;
use oxwave_pfa_fft::{npfa, npfar, pfa2cc, pfa2rc};
use oxwave_tempfile::temporary_filename;
use std::ffi::{c_char, CString};
use std::sync::atomic::{AtomicPtr, Ordering};

unsafe extern "C" {
    fn access(path: *const c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUMIGPSPI - Gazdag's phase-shift plus interpolation depth migration   ",
    "            for zero-offset data, which can handle the lateral         ",
    "            velocity variation.                                        ",
    "                                                                       ",
    " sumigpspi <infile >outfile vfile= [optional parameters]               ",
    "                                                                       ",
    " Required Parameters:\t\t\t\t\t\t\t",
    " nz=\t\tnumber of depth samples\t\t\t\t\t",
    " dz=\t\tdepth sampling interval\t\t\t\t\t",
    " vfile=\tname of file containing velocities\t\t\t",
    "\t\t(Please see Notes below concerning the format of vfile)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:                                                  ",
    " dt=from header(dt) or .004    time sampling interval                  ",
    " dx=from header(d2) or 1.0     midpoint sampling interval              ",
    "                                                                       ",
    " tmpdir=        if non-empty, use the value as a directory path        ",
    "                prefix for storing temporary files; else if the        ",
    "                the CWP_TMPDIR environment variable is set use         ",
    "                its value for the path; else use tmpfile()             ",
    "                                                                       ",
    " Notes:\t\t\t\t\t\t\t\t",
    " The input velocity file 'vfile' consists of C-style binary floats.\t",
    " The structure of this file is vfile[iz][ix]. Note that this means that",
    " the x-direction is the fastest direction instead of z-direction! Such a",
    " structure is more convenient for the downward continuation type\t",
    " migration algorithm than using z as fastest dimension as in other SU\t",
    " programs. (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]\t",
    " denotes a v(z,x) array, the opposite of what Matlab and Fortran\t",
    " programmers may expect.)\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Because most of the tools in the SU package (such as  unif2, unisam2,\t",
    " and makevel) produce output with the structure vfile[ix][iz], you will",
    " need to transpose the velocity files created by these programs. You may",
    " use the SU program 'transp' in SU to transpose such files into the\t",
    " required vfile[iz][ix] structure.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

/* Globals (so can trap signal) defining temporary disk files */
#[derive(Default)]
struct TempFiles {
    tracefile: String,
    headerfile: String,
    tracefp: Option<CFile>,
    headerfp: Option<CFile>,
}

static FILES: AtomicPtr<TempFiles> = AtomicPtr::new(std::ptr::null_mut());

// Preserve C loop variables even where their values are unused.
#[allow(non_snake_case, unused_variables, unused_assignments)]
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut dt = 0.0f32;
    let mut dx = 0.0f32;
    let mut nz = 0i32;
    let mut dz = 0.0f32;
    let mut vfile = String::new();
    let mut verbose = 1i32;
    let mut tmpdir = String::new();

    /* hook up getpar to handle the parameters */
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get info from first trace */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;

    /* let user give dt and/or dx from command line */
    if getparfloat("dt", &mut dt) == 0 {
        if tr.dt != 0 {
            dt = (tr.dt as f64 / 1000000.0) as f32;
        } else {
            dt = 0.004;
            warn!("tr.dt not set, assuming dt=0.004");
        }
    }
    if getparfloat("dx", &mut dx) == 0 {
        if tr.d2 != 0.0 {
            dx = tr.d2;
        } else {
            dx = 1.0;
            warn!("tr.d2 not set, assuming d2=1.0");
        }
    }

    /* get optional parameters */
    if getparint("nz", &mut nz) == 0 {
        err!("nz must be specified");
    }
    if getparfloat("dz", &mut dz) == 0 {
        err!("dz must be specified");
    }
    if getparstring("vfile", &mut vfile) == 0 {
        err!("vfile must be specified");
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let path = CString::new(tmpdir.as_bytes()).unwrap();
        // Match access(tmpdir, WRITE_OK), which checks this user's access.
        if unsafe { access(path.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();

    /* store traces and headers in tempfiles while getting a count */
    let files = Box::leak(Box::new(TempFiles::default()));
    FILES.store(files as *mut TempFiles, Ordering::Relaxed);
    let mut istmpdir = false;
    if tmpdir.is_empty() {
        files.tracefp = Some(etmpfile());
        files.headerfp = Some(etmpfile());
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        files.tracefile = temporary_filename(&tmpdir);
        files.headerfile = temporary_filename(&tmpdir);
        /* Trap signals so can remove temp files */
        for sig in [2, 3, 1, 15] {
            unsafe { signal(sig, closefiles) };
        }
        files.tracefp = Some(efopen(&files.tracefile, "w+"));
        files.headerfp = Some(efopen(&files.headerfile, "w+"));
        istmpdir = true;
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
    }
    let tracefp = files.tracefp.as_ref().unwrap().clone();
    let headerfp = files.headerfp.as_ref().unwrap().clone();
    let mut nx = 0usize;
    loop {
        nx += 1;
        // C writes one HDRBYTES-sized item; preserve efwrite's error counts.
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

    /* determine wavenumber sampling (for complex to complex FFT) */
    let nxfft = npfa(nx as i32) as usize;
    let nk = nxfft;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let fk = (-PI / dx as f64) as f32;

    /* allocate space: contiguous arrays, as required by pfa2rc/pfa2cc */
    // Negative C allocation dimensions do not prevent the later vfile open.
    // Preserve that error ordering without an enormous Rust allocation.
    let nz = if nz < 0 { 0 } else { nz as usize };
    let mut p = vec![0.0f32; ntfft * nx];
    let mut cp = vec![Complex::default(); nw * nx];
    let mut cq = vec![Complex::default(); nw * nk];
    let mut cq1 = vec![Complex::default(); nw * nk];
    let mut cresult = vec![0.0f32; nz * nx];
    let mut v = vec![0.0f32; nx * nz];

    /* load traces into the zero-offset array and close tmpfile */
    for ix in 0..nx {
        efread(&mut p[ix * ntfft..ix * ntfft + nt], &tracefp);
    }
    efclose(&tracefp);

    /* load velocity file */
    let vfp = efopen(&vfile, "r");
    efread(&mut v, &vfp);
    efclose(&vfp);
    // C reads v[0][0] even for nz=0, but never uses it in that case.
    let mut vmax = v.first().copied().unwrap_or(0.0);
    let mut vmin = vmax;
    for iz in 0..nz {
        for ix in 0..nx {
            if v[iz * nx + ix] >= vmax {
                vmax = v[iz * nx + ix];
            }
            if v[iz * nx + ix] <= vmin {
                vmin = v[iz * nx + ix];
            }
        }
    }

    /* pad with zeros and Fourier transform t to w */
    for ix in 0..nx {
        for it in nt..ntfft {
            p[ix * ntfft + it] = 0.0;
        }
    }
    pfa2rc(1, 1, ntfft as i32, nx as i32, &p, &mut cp);

    let L = 10usize;
    let mut c = [0.0f32; 41];
    let mut P = [0.0f32; 40];
    let mut Y = [0.0f32; 41];
    let mut tz = 0.0f32;
    /* loops over depth */
    for iz in 0..nz {
        for ix in 0..nx {
            cresult[ix * nz + iz] = 0.0;
            for iw in 0..nw {
                cresult[ix * nz + iz] += cp[ix * nw + iw].r / ntfft as f32;
            }
        }
        let mut lvmax = v[iz * nx];
        let mut lvmin = v[iz * nx];
        for ix in 0..nx {
            if v[iz * nx + ix] >= lvmax {
                lvmax = v[iz * nx + ix];
            }
            if v[iz * nx + ix] <= lvmin {
                lvmin = v[iz * nx + ix];
            }
        }
        for ik in 0..nx {
            let mut w = fw;
            for iw in 0..nw {
                cp[ik * nw + iw] = cmul(
                    cp[ik * nw + iw],
                    cwp_cexp(cmplx(
                        0.0,
                        ((-w * dz) as f64 * 2.0 / v[iz * nx + ik] as f64) as f32,
                    )),
                );
                cq[ik * nw + iw] = if ik % 2 != 0 {
                    cneg(cp[ik * nw + iw])
                } else {
                    cp[ik * nw + iw]
                };
                w += dw;
            }
        }
        for ik in nx..nk {
            for iw in 0..nw {
                cq[ik * nw + iw] = cmplx(0.0, 0.0);
            }
        }

        /* FFT to W-K domain */
        pfa2cc(-1, 2, nw as i32, nk as i32, &mut cq);
        /* The second time phase shift */
        let v1 = (lvmin as f64 * 0.5) as f32;
        let v2 = (lvmax as f64 * 0.5) as f32;
        if (((v2 - v1) / v1) as f64) < 0.01 {
            let mut k = fk;
            for ik in 0..nk {
                let mut w = fw;
                for iw in 0..nw {
                    if w == 0.0 {
                        w = (1.0e-10 / dt as f64) as f32;
                    }
                    let kz1 = 1.0 - (((v1 * k) / w) as f64).powf(2.0);
                    let cshift1;
                    if kz1 > 0.0 {
                        let phase1 =
                            -(w as f64) * kz1.sqrt() * dz as f64 / v1 as f64 + (w * dz / v1) as f64;
                        cshift1 = cmplx(phase1.cos() as f32, phase1.sin() as f32);
                    } else {
                        let phase1 = -(w as f64) * (-kz1).sqrt() * dz as f64 / v1 as f64;
                        cshift1 = cwp_cexp(cmplx(phase1 as f32, w * dz / v1));
                    }
                    cq1[ik * nw + iw] = cmul(cq[ik * nw + iw], cshift1);
                    w += dw;
                }
                k += dk;
            }
            pfa2cc(1, 2, nw as i32, nk as i32, &mut cq1);
            for ix in 0..nx {
                for iw in 0..nw {
                    cq1[ix * nw + iw] = crmul(cq1[ix * nw + iw], (1.0 / nxfft as f64) as f32);
                    cp[ix * nw + iw] = if ix % 2 != 0 {
                        cneg(cq1[ix * nw + iw])
                    } else {
                        cq1[ix * nw + iw]
                    };
                }
            }
        } else {
            for ik in 0..=L {
                c[ik] = (vmin as f64 + ik as f64 * 1.0 * (vmax - vmin) as f64 / (L as f64 * 1.0))
                    as f32;
            }
            for ik in 0..L {
                P[ik] = 0.0;
            }
            for ix in 0..nx {
                for ik in 0..L {
                    if (v[iz * nx + ix] >= c[ik] && v[iz * nx + ix] < c[ik + 1])
                        || (ik == L - 1 && v[iz * nx + ix] == vmax)
                    {
                        P[ik] = (P[ik] as f64 + 1.0 / nx as f64) as f32;
                        break;
                    }
                }
            }
            let mut Sz = 0.0f32;
            for ik in 0..L {
                if P[ik] != 0.0 {
                    Sz = (Sz as f64 - P[ik] as f64 * (P[ik] as f64).ln()) as f32;
                }
            }
            let Bz = ((Sz as f64).exp() + 0.5) as i32 as usize;
            Y[0] = 0.0;
            Y[L] = 1.0;
            for ik in 1..L {
                Y[ik] = 0.0;
                for ix in 0..ik {
                    Y[ik] += P[ix];
                }
            }
            let mut V = vec![0.0f32; Bz + 1];
            V[0] = vmin;
            for ix in 1..=Bz {
                for ik in 0..L {
                    if ix as f64 * 1.0 / Bz as f64 > Y[ik] as f64
                        && ix as f64 * 1.0 / Bz as f64 <= Y[ik + 1] as f64
                    {
                        V[ix] = (c[ik] as f64
                            + (ix as f64 * 1.0 / Bz as f64 - Y[ik] as f64)
                                * (c[ik + 1] - c[ik]) as f64
                                / (Y[ik + 1] - Y[ik]) as f64)
                            as f32;
                        break;
                    }
                }
            }
            V[Bz] = (V[Bz] as f64 * 1.005) as f32;
            let mut cq2 = vec![Complex::default(); nw * nk * (Bz + 1)];
            for ix in 0..=Bz {
                let mut w = fw;
                for iw in 0..nw {
                    let mut k = fk;
                    for ik in 0..nk {
                        if w == 0.0 {
                            w = (1.0e-10 / dt as f64) as f32;
                        }
                        let kz1 = 1.0 - (V[ix] as f64 / 2.0 * k as f64 / w as f64).powf(2.0);
                        let cshift1;
                        if kz1 >= 0.0 {
                            let phase1 = -(w as f64) * kz1.sqrt() * dz as f64 * 2.0 / V[ix] as f64
                                + (w * dz) as f64 * 2.0 / V[ix] as f64;
                            cshift1 = cwp_cexp(cmplx(0.0, phase1 as f32));
                        } else {
                            let phase1 =
                                -(w as f64) * (-kz1).sqrt() * dz as f64 * 2.0 / V[ix] as f64;
                            cshift1 = cwp_cexp(cmplx(
                                phase1 as f32,
                                ((w * dz) as f64 * 2.0 / V[ix] as f64) as f32,
                            ));
                        }
                        cq2[(ix * nk + ik) * nw + iw] = cmul(cq[ik * nw + iw], cshift1);
                        k += dk;
                    }
                    w += dw;
                }
                pfa2cc(
                    1,
                    2,
                    nw as i32,
                    nk as i32,
                    &mut cq2[ix * nk * nw..(ix + 1) * nk * nw],
                );
                for ik in 0..nx {
                    let mut w = fw;
                    for iw in 0..nw {
                        let a = 0.015f32;
                        let mut g = 1.0f32;
                        let i = 20i32;
                        let ik_i = ik as i32;
                        let nx_i = nx as i32;
                        if ik_i <= i {
                            g = ((-a * (i - ik_i) as f32 * (i - ik_i) as f32) as f64).exp() as f32;
                        }
                        if ik_i >= nx_i - i {
                            g = ((-a * (-nx_i + i + ik_i) as f32 * (-nx_i + i + ik_i) as f32)
                                as f64)
                                .exp() as f32;
                        }
                        let index = (ix * nk + ik) * nw + iw;
                        cq2[index] = crmul(cq2[index], (g as f64 * 1.0 / nxfft as f64) as f32);
                        cq2[index] = if ik % 2 != 0 {
                            cneg(cq2[index])
                        } else {
                            cq2[index]
                        };
                        w += dw;
                    }
                }
            }
            for ix in 0..nx {
                for ik in 0..Bz {
                    if v[iz * nx + ix] >= V[ik] && v[iz * nx + ix] < V[ik + 1] {
                        let v1 = V[ik];
                        let v2 = V[ik + 1];
                        let mut w = fw;
                        for iw in 0..nw {
                            let a1 = cq2[(ik * nk + ix) * nw + iw].r as f64;
                            let a2 = cq2[((ik + 1) * nk + ix) * nw + iw].r as f64;
                            let theta1 = cq2[(ik * nk + ix) * nw + iw].i as f64;
                            let theta2 = cq2[((ik + 1) * nk + ix) * nw + iw].i as f64;
                            let a = a1 * (v2 - v[iz * nx + ix]) as f64 / (v2 - v1) as f64
                                + a2 * (v[iz * nx + ix] - v1) as f64 / (v2 - v1) as f64;
                            let theta = theta1 * (v2 - v[iz * nx + ix]) as f64 / (v2 - v1) as f64
                                + theta2 * (v[iz * nx + ix] - v1) as f64 / (v2 - v1) as f64;
                            cp[ix * nw + iw] = cmplx(a as f32, theta as f32);
                            w += dw;
                        }
                        break;
                    }
                }
            }
        }
        tz += dz;
    }

    /* restore header fields and write output */
    for ix in 0..nx {
        efread(std::slice::from_mut(&mut tr.hdr), &headerfp);
        tr.ns = nz as u16;
        tr.d1 = dz;
        tr.data[..nz].copy_from_slice(&cresult[ix * nz..(ix + 1) * nz]);
        puttr(&mut tr);
    }
    /* Clean up */
    efclose(&headerfp);
    if istmpdir {
        eremove(&files.headerfile);
    }
    if istmpdir {
        eremove(&files.tracefile);
    }
    let _ = tz;
    exit(CWP_Exit());
}

extern "C" fn closefiles(_sig: i32) {
    // main installs the handlers only after publishing the leaked file state.
    let files = unsafe { &*FILES.load(Ordering::Relaxed) };
    if let Some(fp) = &files.headerfp {
        efclose(fp);
    }
    if let Some(fp) = &files.tracefp {
        efclose(fp);
    }
    eremove(&files.headerfile);
    eremove(&files.tracefile);
    exit(EXIT_FAILURE);
}
