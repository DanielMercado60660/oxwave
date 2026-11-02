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

//! SUMIGFD - 45-90 degree finite difference depth migration for zero-offset data.
//! Port of su/main/migration_inversion/sumigfd.c (Revision 1.13).
//! Credits (from the C source): CWP Baoniu Han, April 20th, 1998.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::PI;
use oxwave_pfa_fft::{npfar, pfa2rc};
use oxwave_tempfile::{temporary_filename, temporary_stream};
use std::sync::atomic::{AtomicBool, AtomicPtr, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGFD - 45-90 degree Finite difference depth migration for\t\t",
    "           zero-offset data.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   sumigfd <infile >outfile vfile= [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " nz=\t\tnumber of depth samples\t\t\t\t\t",
    " dz=\t\tdepth sampling interval\t\t\t\t\t",
    " vfile=\tname of file containing velocities\t\t\t",
    " \t\t(see Notes below concerning format of this file)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt=from header(dt) or .004    time sampling interval\t\t\t",
    " dx=from header(d2) or 1.0\tmidpoint sampling interval\t\t",
    " dip=45,65,79,80,87,89,90  \tMaximum angle of dip reflector\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir=\tif non-empty, use the value as a directory path\t\t",
    "\t\tprefix for storing temporary files; else if the\t\t",
    "\t\tthe CWP_TMPDIR environment variable is set use\t\t",
    "\t\tits value for the path; else use tmpfile()\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The computation cost by dip angle is 45=65=79<80<87<89<90\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The input velocity file 'vfile' consists of C-style binary floats.\t",
    " The structure of this file is vfile[iz][ix]. Note that this means that",
    " the x-direction is the fastest direction instead of z-direction! Such a",
    " structure is more convenient for the downward continuation type\t",
    " migration algorithm than using z as fastest dimension as in other SU\t",
    " programs. (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  ",
    " denotes a v(z,x) array, the opposite of what Matlab and Fortran\t",
    " programmers may expect.)\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Because most of the tools in the SU package (such as  unif2, unisam2,\t",
    " and makevel) produce output with the structure vfile[ix][iz], you will",
    " need to transpose the velocity files created by these programs. You may",
    " use the SU program 'transp' in SU to transpose such files into the\t",
    " required vfile[iz][ix] structure.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

/* Temporary disk files are owned by main; paths are shared with the handler. */
struct TempFiles {
    tracefile: String,
    headerfile: String,
    tracefp: CFile,
    headerfp: CFile,
}
static INTERRUPTED: AtomicBool = AtomicBool::new(false);
static READING_TRACES: AtomicBool = AtomicBool::new(false);
static TRACEFILE: AtomicPtr<std::ffi::c_char> = AtomicPtr::new(std::ptr::null_mut());
static HEADERFILE: AtomicPtr<std::ffi::c_char> = AtomicPtr::new(std::ptr::null_mut());

// POSIX filesystem and signal APIs; all SU numerical routines use Rust crates.
unsafe extern "C" {
    fn access(path: *const std::ffi::c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
    fn unlink(path: *const std::ffi::c_char) -> i32;
    fn _exit(status: i32) -> !;
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

    let mut nz = 0i32;
    let mut dz = 0.0f32;
    let mut vfile = String::new();
    let mut verbose = 1i32;
    let mut dip = 0i32;
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
    if getparint("dip", &mut dip) == 0 {
        dip = 65;
    }

    /* Look for user-supplied tmpdir */
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let directory = std::ffi::CString::new(tmpdir.as_bytes()).unwrap();
        // SAFETY: NUL-terminated directory; POSIX WRITE_OK is 2.
        if unsafe { access(directory.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();

    // C copies nz samples into a fixed-size trace. Reject dimensions that
    // would overflow it (or turn negative allocation sizes into huge ones).
    if nz < 0 || nz as usize > tr.data.len() {
        err!("nz must be between 0 and %d", tr.data.len());
    }
    let nz = nz as usize;

    /* store traces and headers in tempfiles while getting a count */
    let istmpdir = !tmpdir.is_empty();
    let files = if !istmpdir {
        let tracefp = etmpfile();
        let headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
        TempFiles {
            tracefile: String::new(),
            headerfile: String::new(),
            tracefp,
            headerfp,
        }
    } else {
        let tracefile = temporary_filename(&tmpdir);
        let headerfile = temporary_filename(&tmpdir);
        // SAFETY: install the four handlers used in the original program.
        unsafe {
            signal(2, record_signal as *const () as usize);
            signal(3, record_signal as *const () as usize);
            signal(1, record_signal as *const () as usize);
            signal(15, record_signal as *const () as usize);
        }
        // When /tmp is unavailable, the library returns the same fallback
        // name twice. Use its directory-local streams in that case so the
        // header and sample files remain independent, as in the original.
        let (tracefp, headerfp, tracefile, headerfile) = if tracefile == headerfile {
            let tracefp = temporary_stream(Some(&tmpdir))
                .map(CFile::from_file)
                .unwrap_or_else(|| efopen(&tracefile, "w+"));
            let headerfp = temporary_stream(Some(&tmpdir))
                .map(CFile::from_file)
                .unwrap_or_else(|| efopen(&headerfile, "w+"));
            (tracefp, headerfp, String::new(), String::new())
        } else {
            let tracefp = efopen(&tracefile, "w+");
            let headerfp = efopen(&headerfile, "w+");
            (tracefp, headerfp, tracefile, headerfile)
        };
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
        TempFiles {
            tracefile,
            headerfile,
            tracefp,
            headerfp,
        }
    };
    let tracefp = &files.tracefp;
    let headerfp = &files.headerfp;
    // Keep these NUL-terminated paths alive until exit. Publishing the input
    // phase last lets the handler safely remove named files during gettr.
    let tracepath = std::ffi::CString::new(files.tracefile.as_bytes()).unwrap();
    let headerpath = std::ffi::CString::new(files.headerfile.as_bytes()).unwrap();
    TRACEFILE.store(tracepath.as_ptr().cast_mut(), Ordering::Relaxed);
    HEADERFILE.store(headerpath.as_ptr().cast_mut(), Ordering::Relaxed);
    READING_TRACES.store(true, Ordering::Release);
    let mut nx = 0usize;
    loop {
        check_signal(&files);
        nx += 1;
        efwrite(tr.hdr.as_bytes(), headerfp);
        efwrite(&tr.data[..nt], tracefp);
        check_signal(&files);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    READING_TRACES.store(false, Ordering::Release);
    check_signal(&files);
    erewind(tracefp);
    erewind(headerfp);

    /* determine frequency sampling (for real to complex FFT) */
    let ntfft = npfar(nt as i32) as usize;
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = 0.0f32;

    /* allocate and zero space; C arrays are contiguous */
    let mut p = vec![0.0f32; ntfft * nx];
    let mut cp = vec![Complex::default(); nw * nx];
    let mut cresult = vec![0.0f32; nz * nx];
    let mut v = vec![0.0f32; nz * nx];

    for ix in 0..nx {
        check_signal(&files);
        efread(&mut p[ix * ntfft..ix * ntfft + nt], tracefp);
    }
    efclose(tracefp);

    let vfp = efopen(&vfile, "r");
    efread(&mut v, &vfp);
    efclose(&vfp);
    check_signal(&files);

    pfa2rc(1, 1, ntfft as i32, nx as i32, &p, &mut cp);
    check_signal(&files);

    /* loop over depth */
    let mut tz = 0.0f32;
    for iz in 0..nz {
        check_signal(&files);
        for ix in 0..nx {
            cresult[ix * nz + iz] = 0.0;
            v[iz * nx + ix] = (v[iz * nx + ix] as f64 / 2.0) as f32;
            for iw in 1..nw {
                cresult[ix * nz + iz] += cp[ix * nw + iw].r / ntfft as f32;
            }
        }
        for ix in 0..nx {
            let mut w = fw + dw;
            for iw in 1..nw {
                let kz2 = -1.0 / v[iz * nx + ix] as f64 * w as f64 * dz as f64;
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
            dip,
            &files,
        );
        tz += dz;
    }
    let _ = tz;

    /* restore header fields and write output */
    for ix in 0..nx {
        check_signal(&files);
        efread(tr.hdr.as_bytes_mut(), headerfp);
        tr.ns = nz as u16;
        tr.d1 = dz;
        tr.data[..nz].copy_from_slice(&cresult[ix * nz..(ix + 1) * nz]);
        puttr(&mut tr);
    }
    efclose(headerfp);
    if istmpdir && !files.headerfile.is_empty() {
        eremove(&files.headerfile);
    }
    if istmpdir && !files.tracefile.is_empty() {
        eremove(&files.tracefile);
    }
    exit(if INTERRUPTED.load(Ordering::Relaxed) {
        1
    } else {
        CWP_Exit()
    });
}

extern "C" fn record_signal(_sig: i32) {
    INTERRUPTED.store(true, Ordering::Relaxed);
    if READING_TRACES.load(Ordering::Acquire) {
        // No output exists yet. Unlike deferred cleanup, these async-signal-
        // safe calls also work while gettr blocks with a CFile borrow held.
        // Process termination closes both temporary file descriptors.
        unsafe {
            unlink(HEADERFILE.load(Ordering::Relaxed));
            unlink(TRACEFILE.load(Ordering::Relaxed));
            _exit(1);
        }
    }
}

fn check_signal(files: &TempFiles) {
    if INTERRUPTED.load(Ordering::Relaxed) {
        // Called only after I/O routines have released their RefCell borrows.
        efclose(&files.headerfp);
        efclose(&files.tracefp);
        if !files.headerfile.is_empty() {
            eremove(&files.headerfile);
        }
        if !files.tracefile.is_empty() {
            eremove(&files.tracefile);
        }
        exit(1);
    }
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
    mut dip: i32,
    files: &TempFiles,
) {
    let mut s1 = vec![0.0f32; nx];
    let mut s2 = vec![0.0f32; nx];
    let mut data = vec![Complex::default(); nx];
    let mut d = vec![Complex::default(); nx];
    let mut a = vec![Complex::default(); nx];
    let mut b = vec![Complex::default(); nx];
    let mut c = vec![Complex::default(); nx];
    let trick = 0.1f32;
    let mut coefa = [0.0f32; 5];
    let mut coefb = [0.0f32; 5];
    let mut step = 1;

    if ![45, 65, 79, 80, 87, 89, 90].contains(&dip) {
        dip = 79;
    }
    if dip == 45 {
        coefa[0] = 0.5;
        coefb[0] = 0.25;
        step = 1;
    }
    if dip == 65 {
        coefa[0] = 0.478242060;
        coefb[0] = 0.376369527;
        step = 1;
    }
    if dip == 79 {
        coefa[0] = 0.4575;
        coefb[0] = 0.4575;
        step = 1;
    }
    if dip == 80 {
        coefa[1] = 0.040315157;
        coefb[1] = 0.873981642;
        coefa[0] = 0.457289566;
        coefb[0] = 0.222691983;
        step = 2;
    }
    if dip == 87 {
        coefa[2] = 0.00421042;
        coefb[2] = 0.972926132;
        coefa[1] = 0.081312882;
        coefb[1] = 0.744418059;
        coefa[0] = 0.414236605;
        coefb[0] = 0.150843924;
        step = 3;
    }
    if dip == 89 {
        coefa[3] = 0.000523275;
        coefb[3] = 0.994065088;
        coefa[2] = 0.014853510;
        coefb[2] = 0.919432661;
        coefa[1] = 0.117592008;
        coefb[1] = 0.614520676;
        coefa[0] = 0.367013245;
        coefb[0] = 0.105756624;
        step = 4;
    }
    if dip == 90 {
        coefa[4] = 0.000153427;
        coefb[4] = 0.997370236;
        coefa[3] = 0.004172967;
        coefb[3] = 0.964827992;
        coefa[2] = 0.033860918;
        coefb[2] = 0.824918565;
        coefa[1] = 0.143798076;
        coefb[1] = 0.483340757;
        coefa[0] = 0.318013812;
        coefb[0] = 0.073588213;
        step = 5;
    }
    let v1 = v[0];
    let vn = v[nx - 1];
    loop {
        step -= 1;
        let mut w = fw;
        for iw in 0..nw {
            check_signal(files);
            if (w as f64).abs() <= 1.0e-10 {
                w = (1.0e-10 / dt as f64) as f32;
            }
            for ix in 0..nx {
                s1[ix] = (v[ix] * v[ix]) * coefb[step] / (dx * dx * w * w) + trick;
                s2[ix] = ((-v[ix] * dz * coefa[step] / (w * dx * dx)) as f64 * 0.5) as f32;
            }
            for ix in 0..nx {
                data[ix] = cp[ix * nw + iw];
            }
            let cp2 = data[0];
            // For nx=1 the C routine reads both neighboring allocator slots.
            // They are zero in this build; provide those guard values explicitly.
            let cp3 = if nx == 1 { Complex::default() } else { data[1] };
            let cpnm1 = data[nx - 1];
            let cpnm2 = if nx == 1 {
                Complex::default()
            } else {
                data[nx - 2]
            };
            let mut a1 = cmul(cp2, conjg(cp3));
            let b1 = cmul(cp3, conjg(cp3));
            if b1.r == 0.0 && b1.i == 0.0 {
                a1 = cwp_cexp(cmplx(0.0, ((-w * dx) as f64 * 0.5 / v1 as f64) as f32));
            } else {
                a1 = cdiv(a1, b1);
            }
            if a1.i > 0.0 {
                a1 = cwp_cexp(cmplx(0.0, ((-w * dx) as f64 * 0.5 / v1 as f64) as f32));
            }
            let mut a2 = cmul(cpnm1, conjg(cpnm2));
            let b2 = cmul(cpnm2, conjg(cpnm2));
            if b2.r == 0.0 && b2.i == 0.0 {
                a2 = cwp_cexp(cmplx(0.0, ((-w * dx) as f64 * 0.5 / vn as f64) as f32));
            } else {
                a2 = cdiv(a2, b2);
            }
            if a2.i > 0.0 {
                a2 = cwp_cexp(cmplx(0.0, ((-w * dx) as f64 * 0.5 / vn as f64) as f32));
            }
            for ix in 0..nx {
                a[ix] = cmplx(s1[ix], s2[ix]);
                b[ix] = cmplx(
                    (1.0 - 2.0 * s1[ix] as f64) as f32,
                    (-2.0 * s2[ix] as f64) as f32,
                );
            }
            for ix in 1..nx - 1 {
                d[ix] = cadd(
                    cadd(cmul(data[ix + 1], a[ix + 1]), cmul(data[ix - 1], a[ix - 1])),
                    cmul(data[ix], b[ix]),
                );
            }
            if nx == 1 {
                // The right endpoint assignment overwrites the left one.
                d[0] = cadd(
                    cmul(cadd(b[0], cmul(a[0], a2)), data[0]),
                    cmul(Complex::default(), Complex::default()),
                );
            } else {
                d[0] = cadd(
                    cmul(cadd(b[0], cmul(a[0], a1)), data[0]),
                    cmul(data[1], a[1]),
                );
                d[nx - 1] = cadd(
                    cmul(cadd(b[nx - 1], cmul(a[nx - 1], a2)), data[nx - 1]),
                    cmul(data[nx - 2], a[nx - 2]),
                );
            }
            for ix in 0..nx {
                data[ix] = cmplx(s1[ix], -s2[ix]);
                b[ix] = cmplx(
                    (1.0 - 2.0 * s1[ix] as f64) as f32,
                    (2.0 * s2[ix] as f64) as f32,
                );
            }
            let endl = cadd(b[0], cmul(data[0], a1));
            let endr = cadd(b[nx - 1], cmul(data[nx - 1], a2));
            for ix in 1..nx - 1 {
                a[ix] = data[ix + 1];
                c[ix] = data[ix - 1];
            }
            a[0] = if nx == 1 { Complex::default() } else { data[1] };
            c[nx - 1] = if nx == 1 {
                Complex::default()
            } else {
                data[nx - 2]
            };
            retris(&mut data, &a, &c, &b, endl, endr, nx, &d);
            for ix in 0..nx {
                cp[ix * nw + iw] = data[ix];
            }
            w += dw;
        }
        if step == 0 {
            break;
        }
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
    if nx == 1 {
        // The C boundary solve reads zero guard slots at index -1.
        data[0] = cdiv(
            csub(d[0], cmul(Complex::default(), Complex::default())),
            cadd(endr, cmul(Complex::default(), Complex::default())),
        );
        return;
    }
    for ix in 1..nx - 1 {
        let den = cadd(b[ix], cmul(c[ix], e[ix - 1]));
        e[ix] = cdiv(cneg(a[ix]), den);
        f[ix] = cdiv(csub(d[ix], cmul(f[ix - 1], c[ix])), den);
    }
    // Preserve c[nx-2], exactly as in the original boundary solve.
    data[nx - 1] = cdiv(
        csub(d[nx - 1], cmul(f[nx - 2], c[nx - 2])),
        cadd(endr, cmul(c[nx - 2], e[nx - 2])),
    );
    for ix in (0..nx - 1).rev() {
        data[ix] = cadd(cmul(data[ix + 1], e[ix]), f[ix]);
    }
}
