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

//! SUMIGPS - Migration by phase shift with turning rays.
//! Port of su/main/migration_inversion/sumigps.c (Revision 1.15).
//! Credits (from the C source): CWP: Dave Hale (originally called supsmig.c).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
use oxwave_cwp::cwp::{abs, max, min, nint, PI};
use oxwave_index_search::xindex;
use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::{npfao, npfaro, pfa2cr, pfa2rc, pfacc};
use oxwave_tempfile::{temporary_filename, temporary_stream};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGPS - MIGration by Phase Shift with turning rays\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sumigps <stdin >stdout [optional parms]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " \tNone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt=from header(dt) or .004\ttime sampling interval\t\t\t",
    " dx=from header(d2) or 1.0\tdistance between sucessive cdp's\t",
    " ffil=0,0,0.5/dt,0.5/dt  trapezoidal window of frequencies to migrate\t",
    " tmig=0.0\t\ttimes corresponding to interval velocities in vmig",
    " vmig=1500.0\t\tinterval velocities corresponding to times in tmig",
    " vfile=\t\tbinary (non-ascii) file containing velocities v(t)",
    " nxpad=0\t\tnumber of cdps to pad with zeros before FFT\t",
    " ltaper=0\t\tlength of linear taper for left and right edges",
    " verbose=0\t\t=1 for diagnostic print\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input traces must be sorted by either increasing or decreasing cdp.\t",
    "\t\t\t\t\t\t\t\t\t",
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
];

struct TATable {
    np: usize,
    p: Vec<f32>,
    taumax: Vec<f32>,
    taumin: Vec<f32>,
    tturn: Vec<f32>,
    t: Vec<Vec<f32>>,
    a: Vec<Vec<f32>>,
}

/* Globals (so can trap signal) defining temporary disk files. */
struct TempFiles {
    tracefile: String,
    headerfile: String,
    tracefp: CFile,
    headerfp: CFile,
}
static TEMPFILES: std::sync::atomic::AtomicPtr<TempFiles> =
    std::sync::atomic::AtomicPtr::new(std::ptr::null_mut());

// POSIX OS APIs only. Numerical SU/CWP routines come from the OxWave Rust
// library crates; this executable does not link the original SU libraries or cref.
unsafe extern "C" {
    fn access(path: *const std::ffi::c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
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
    let nt = tr.ns as usize;

    /* let user give dt and/or dx from command line */
    let mut dt: f32 = 0.0;
    if getparfloat("dt", &mut dt) == 0 {
        if tr.dt != 0 {
            dt = (tr.dt as f64 / 1000000.0) as f32;
        } else {
            dt = 0.004;
            warn!("tr.dt not set, assuming dt=0.004");
        }
    }
    let mut dx: f32 = 0.0;
    if getparfloat("dx", &mut dx) == 0 {
        if tr.d2 != 0.0 {
            dx = tr.d2;
        } else {
            dx = 1.0;
            warn!("tr.d2 not set, assuming d2=1.0");
        }
    }

    /* get parameters */
    let mut ffil = Vec::<f32>::new();
    let nfil = getparfloat("ffil", &mut ffil);
    if nfil == 0 {
        ffil = vec![0.0, 0.0, (0.5 / dt as f64) as f32, (0.5 / dt as f64) as f32];
    } else if nfil != 4 {
        err!("if ffil specified, exactly 4 values must be provided");
    }
    let mut nxpad = 0i32;
    if getparint("nxpad", &mut nxpad) == 0 {
        nxpad = 0;
    }
    let mut ltaper = 0i32;
    if getparint("ltaper", &mut ltaper) == 0 {
        ltaper = 0;
    }
    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    let mut np = 200i32;
    if getparint("np", &mut np) == 0 {
        np = 200;
    }
    let mut ntflag = 1i32;
    if getparint("ntflag", &mut ntflag) == 0 {
        ntflag = 1;
    }

    /* look for user-supplied tmpdir */
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let directory = std::ffi::CString::new(tmpdir.as_bytes()).unwrap();
        // SAFETY: directory is NUL terminated; WRITE_OK is POSIX W_OK (2).
        if unsafe { access(directory.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }

    /* determine velocity function v(t) */
    let mut vt = vec![0.0f32; nt];
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
        for it in 0..nt {
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
        // The C uses unchecked fopen here; NULL subsequently faults in fread.
        let fp = match fopen(&vfile, "r") {
            Some(fp) => fp,
            None => {
                unsafe extern "C" {
                    fn raise(sig: i32) -> i32;
                }
                // SAFETY: reproduce the original NULL FILE* failure (SIGSEGV).
                unsafe {
                    signal(11, 0);
                    raise(11);
                }
                exit(139);
            }
        };
        if efread(&mut vt, &fp) != nt {
            err!("cannot read %d velocities from file %s", nt, vfile);
        }
    }
    checkpars();

    /* copy traces and headers to temporary files */
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
        // SAFETY: install the same four termination handlers as the C source.
        unsafe {
            signal(2, closefiles as *const () as usize);
            signal(3, closefiles as *const () as usize);
            signal(1, closefiles as *const () as usize);
            signal(15, closefiles as *const () as usize);
        }
        // If /tmp is inaccessible, the library returns the same fallback
        // name twice; C mkstemp still changes its template on failure.
        // Use the library's independent, already unlinked streams in the
        // requested directory in that case. Keep normal named-file handling.
        let (tracefp, headerfp, tracefile, headerfile) = if tracefile == headerfile {
            let tracefp = temporary_stream(Some(&tmpdir))
                .map(CFile::from_file)
                .unwrap_or_else(|| err!("subcalls.c: efopen: fopen failed"));
            let headerfp = temporary_stream(Some(&tmpdir))
                .map(CFile::from_file)
                .unwrap_or_else(|| err!("subcalls.c: efopen: fopen failed"));
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
    let files = Box::leak(Box::new(files));
    TEMPFILES.store(files, std::sync::atomic::Ordering::Release);
    let tracefp = &files.tracefp;
    let headerfp = &files.headerfp;

    /* count the traces */
    let mut nx = 0i32;
    loop {
        nx += 1;
        efwrite(tr.hdr.as_bytes(), headerfp);
        efwrite(&tr.data[..nt], tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    erewind(headerfp);
    erewind(tracefp);
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

    /* C overlays gtx on gtk; distinct buffers preserve the same layout. */
    let mut gtk = vec![Complex::default(); nt * nk as usize];
    let mut gtx = vec![0.0f32; nt * nxfft as usize];
    let fftscl = (1.0 / nxfft as f64) as f32;
    for ix in 0..nx as usize {
        efread(&mut gtx[ix * nt..(ix + 1) * nt], tracefp);
        for it in 0..nt {
            gtx[ix * nt + it] *= fftscl;
        }
        if (ix as i32) < ltaper {
            let taper = (ix as i32 + 1) as f32 / ltaper.wrapping_add(1) as f32;
            for it in 0..nt {
                gtx[ix * nt + it] *= taper;
            }
        } else if ix as i32 >= nx.wrapping_sub(ltaper) {
            let taper = (nx - ix as i32) as f32 / ltaper.wrapping_add(1) as f32;
            for it in 0..nt {
                gtx[ix * nt + it] *= taper;
            }
        }
    }
    for ix in nx as usize..nxfft as usize {
        for it in 0..nt {
            gtx[ix * nt + it] = 0.0;
        }
    }

    pfa2rc(-1, 2, nt as i32, nxfft, &gtx, &mut gtk);
    if verbose != 0 {
        fprintf!(stderr(), "\tFourier transform done\n");
    }

    let table = tableta(np as usize, nt, dt, &vt, nt, dt, verbose);
    if verbose != 0 {
        fprintf!(stderr(), "\tTime/amplitude table built\n");
    }

    let trig = TrigTable::new();
    let mut k = 0.0f32;
    for ik in 0..nk {
        if verbose != 0 && ik % (if nk / 10 > 0 { nk / 10 } else { 1 }) == 0 {
            fprintf!(stderr(), "\t%d of %d wavenumbers done\n", ik, nk);
        }
        mig1k(
            &table,
            &trig,
            k,
            &ffil,
            ntflag,
            nt,
            dt,
            &mut gtk[ik as usize * nt..(ik as usize + 1) * nt],
            nt,
            dt,
        );
        k += dk;
    }
    pfa2cr(1, 2, nt as i32, nxfft, &gtk, &mut gtx);
    if verbose != 0 {
        fprintf!(stderr(), "\tinverse Fourier transform done\n");
    }

    /* output migrated traces with headers */
    for ix in 0..nx as usize {
        efread(tr.hdr.as_bytes_mut(), headerfp);
        tr.data[..nt].copy_from_slice(&gtx[ix * nt..(ix + 1) * nt]);
        puttr(&mut tr);
    }
    efclose(headerfp);
    if !files.headerfile.is_empty() {
        eremove(&files.headerfile);
    }
    efclose(tracefp);
    if !files.tracefile.is_empty() {
        eremove(&files.tracefile);
    }
    exit(CWP_Exit());
}

/* The static sine/cosine table in mig1k is owned by main in Rust. */
struct TrigTable {
    ctab: Vec<f32>,
    stab: Vec<f32>,
    fntab: f32,
    opi2: f32,
}

impl TrigTable {
    fn new() -> Self {
        let ntab = 1025;
        let mut ctab = vec![0.0f32; ntab];
        let mut stab = vec![0.0f32; ntab];
        let dangle = (2.0 * PI / (ntab - 1) as f64) as f32;
        let mut angle = 0.0f32;
        for itab in 0..ntab {
            ctab[itab] = (angle as f64).cos() as f32;
            stab[itab] = (angle as f64).sin() as f32;
            angle += dangle;
        }
        Self {
            ctab,
            stab,
            fntab: ntab as f32,
            opi2: (1.0 / (PI * 2.0)) as f32,
        }
    }
}

/* Phase shift migration for one wavenumber k via tabulated times and amps. */
#[allow(clippy::too_many_arguments)]
fn mig1k(
    table: &TATable,
    trig: &TrigTable,
    k: f32,
    ffil: &[f32],
    ntflag: i32,
    nt: usize,
    dt: f32,
    cp: &mut [Complex],
    ntau: usize,
    dtau: f32,
) {
    let np = table.np;
    let p = &table.p;
    let taumax = &table.taumax;
    let taumin = &table.taumin;
    let tturn = &table.tturn;
    let t = &table.t;
    let a = &table.a;
    let ntfft = npfao(2 * nt as i32, 4 * nt as i32);
    let nw = ntfft;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let mut pw = vec![Complex::default(); ntfft as usize];
    let mut qn = vec![Complex::default(); ntau];
    let mut qt = vec![Complex::default(); ntau];
    let wnyq = (PI / dt as f64) as f32;
    let iwnyq = if ntfft % 2 != 0 {
        ntfft / 2
    } else {
        ntfft / 2 + 1
    };
    let iwfil0 = max(0, min(iwnyq, nint(2.0 * PI * ffil[0] as f64 / dw as f64)));
    let iwfil1 = max(0, min(iwnyq, nint(2.0 * PI * ffil[1] as f64 / dw as f64)));
    let iwfil2 = max(0, min(iwnyq, nint(2.0 * PI * ffil[2] as f64 / dw as f64)));
    let iwfil3 = max(0, min(iwnyq, nint(2.0 * PI * ffil[3] as f64 / dw as f64)));

    for it in 0..nt {
        pw[it] = cp[it];
    }
    for it in nt..ntfft as usize {
        pw[it] = Complex { r: 0.0, i: 0.0 };
    }
    pfacc(1, &mut pw);
    for iw in 0..nw as usize {
        pw[iw].r = (pw[iw].r as f64 * (1.0 / ntfft as f64)) as f32;
        pw[iw].i = (pw[iw].i as f64 * (1.0 / ntfft as f64)) as f32;
    }
    for itau in 0..ntau {
        qn[itau] = Complex { r: 0.0, i: 0.0 };
        qt[itau] = Complex { r: 0.0, i: 0.0 };
    }
    let wmin = abs(k) / p[np - 1];
    let iwmin = max(max(1, iwfil0), (wmin / dw) as i32);
    let wmax = if p[0] <= abs(k) / wnyq {
        wnyq
    } else {
        abs(k) / p[0]
    };
    let iwmax = min(min(iwnyq - 1, iwfil3), (wmax / dw) as i32);
    let mut w = iwmin as f32 * dw;
    for iw in iwmin..=iwmax {
        let kow = abs(k) / w;
        if kow > p[np - 1] {
            w += dw;
            continue;
        }
        let ampf = if iwfil0 <= iw && iw < iwfil1 {
            (iw - iwfil0) as f32 / (iwfil1 - iwfil0) as f32
        } else if iwfil2 < iw && iw <= iwfil3 {
            (iwfil3 - iw) as f32 / (iwfil3 - iwfil2) as f32
        } else {
            1.0f32
        };

        let mut ip1 = 0i32;
        xindex(p, kow, &mut ip1);
        let ip1 = min(ip1, np as i32 - 2) as usize;
        let ip2 = ip1 + 1;
        let p1 = p[ip1];
        let p2 = p[ip2];
        let a1 = (p2 * p2 - kow * kow) / (p2 * p2 - p1 * p1);
        let a2 = (kow * kow - p1 * p1) / (p2 * p2 - p1 * p1);
        let wa1 = w * a1;
        let wa2 = w * a2;
        let taumaxi = a1 * taumax[ip1] + a2 * taumax[ip2];
        let itaumax = min(ntau as i32 - 1, nint((taumaxi / dtau) as f64));
        let taumini = a1 * taumin[ip1] + a2 * taumin[ip2];
        let itaumin = max(0, nint((taumini / dtau) as f64));
        let phst = (2.0 * (wa1 * tturn[ip1] + wa2 * tturn[ip2]) as f64) as f32;
        let pwsr = ampf * (pw[iw as usize].r + pw[(nw - iw) as usize].r);
        let pwsi = ampf * (pw[iw as usize].i + pw[(nw - iw) as usize].i);
        let pwdr = ampf * (pw[iw as usize].r - pw[(nw - iw) as usize].r);
        let pwdi = ampf * (pw[iw as usize].i - pw[(nw - iw) as usize].i);
        if ntflag & 1 != 0 {
            for itau in 0..itaumax {
                let itau = itau as usize;
                let ampi = a1 * a[ip1][itau] + a2 * a[ip2][itau];
                let phsi = wa1 * t[ip1][itau] + wa2 * t[ip2][itau];
                let phsn = phsi * trig.opi2;
                let itab = (trig.fntab * (phsn - (phsn as i32) as f32)) as i32 as usize;
                qn[itau].r = qn[itau].r + ampi * (pwsr * trig.ctab[itab] + pwdi * trig.stab[itab]);
                qn[itau].i = qn[itau].i + ampi * (pwsi * trig.ctab[itab] - pwdr * trig.stab[itab]);
            }
        }
        if ntflag & 2 != 0 {
            for itau in itaumin + 1..itaumax {
                let itau = itau as usize;
                let ampi = a1 * a[ip1][itau] + a2 * a[ip2][itau];
                let phsi = phst - (wa1 * t[ip1][itau] + wa2 * t[ip2][itau]);
                let phsn = phsi * trig.opi2;
                let itab = (trig.fntab * (phsn - (phsn as i32) as f32)) as i32 as usize;
                qt[itau].r = qt[itau].r + ampi * (pwsr * trig.stab[itab] - pwdi * trig.ctab[itab]);
                qt[itau].i = qt[itau].i + ampi * (pwsi * trig.stab[itab] + pwdr * trig.ctab[itab]);
            }
        }
        w += dw;
    }
    for itau in 0..ntau {
        cp[itau].r = qn[itau].r + qt[itau].r;
        cp[itau].i = qn[itau].i + qt[itau].i;
    }
}

/* Tabulate time shifts and amplitudes for use in phase-shift migration. */
#[allow(unused_assignments)] // vel updates are retained verbatim from C.
fn tableta(
    np: usize,
    ntau: usize,
    dtau: f32,
    vtau: &[f32],
    nt: usize,
    dt: f32,
    verbose: i32,
) -> TATable {
    // C leaves these extrema uninitialized if the first ray step is negative
    // (possible at the last slope), or if nt=1 skips ray tracing altogether.
    // Zero initialization matches the original's observed fresh allocations;
    // arbitrary values from reused C heap storage cannot be reproduced.
    let mut table = TATable {
        np,
        p: vec![0.0; np],
        taumax: vec![0.0; np],
        taumin: vec![0.0; np],
        tturn: vec![0.0; np],
        t: vec![vec![0.0; ntau]; np],
        a: vec![vec![0.0; ntau]; np],
    };
    let TATable {
        p,
        taumax,
        taumin,
        tturn,
        t,
        a,
        ..
    } = &mut table;
    for jp in 0..np {
        let pj =
            (2.0 / vtau[0] as f64 * ((jp as f32 / (np as i32 - 1) as f32) as f64).sqrt()) as f32;
        p[jp] = pj;
        t[jp][0] = 0.0;
        a[jp][0] = 1.0;
        let mut ktau = 0usize;
        let mut taui = 0.0f32;
        let mut ti = 0.0f32;
        let mut ai = 1.0f32;
        let mut vel = (0.5 * vtau[0] as f64) as f32;
        // The original reads vtau[1] even for nt=1; no ray steps then use it.
        let mut dvel = if ntau > 1 {
            (0.5 * (vtau[1] - vtau[0]) as f64) as f32
        } else {
            0.0
        };
        let mut angle = min(1.0f64, (pj * vel) as f64).asin() as f32;
        for _it in 1..nt {
            let taul = taui;
            let tl = ti;
            let al = ai;
            let cosa = (angle as f64).cos() as f32;
            taui += dt * cosa;
            ti += dt * cosa * cosa;
            ai = 1.0;
            if taui < 0.0 {
                break;
            }
            if taui >= taul {
                tturn[jp] = ti;
                taumax[jp] = taui;
                taumin[jp] = taui;
            } else {
                taumin[jp] = taui;
            }
            let itau = (taui / dtau) as i32;
            let ltau = (taul / dtau) as i32;
            for jtau in ltau + 1..=min(itau, ntau as i32 - 1) {
                let tauj = jtau as f32 * dtau;
                let frac = (tauj - taul) / (taui - taul);
                ktau += 1;
                t[jp][ktau] = ((1.0 - frac as f64) * tl as f64 + (frac * ti) as f64) as f32;
                a[jp][ktau] = ((1.0 - frac as f64) * al as f64 + (frac * ai) as f64) as f32;
            }
            angle += pj * dvel;
            if itau < ntau as i32 - 1 {
                let frac = (taui - itau as f32 * dtau) / dtau;
                dvel = (0.5 * (vtau[itau as usize + 1] - vtau[itau as usize]) as f64) as f32;
                vel = (0.5
                    * ((1.0 - frac as f64) * vtau[itau as usize] as f64
                        + (frac * vtau[itau as usize + 1]) as f64)) as f32;
            } else {
                dvel = (0.5 * (vtau[ntau - 1] - vtau[ntau - 2]) as f64) as f32;
                vel = (0.5 * vtau[ntau - 1] as f64
                    + ((taui - (ntau - 1) as f32 * dtau) * dvel / dtau) as f64)
                    as f32;
            }
        }
        for jtau in ktau + 1..ntau {
            t[jp][jtau] = t[jp][ktau];
            a[jp][jtau] = a[jp][ktau];
        }
        if verbose != 0 {
            fprintf!(
                stderr(),
                "p=%g tturn=%g taumax=%g taumin=%g\n",
                p[jp],
                tturn[jp],
                taumax[jp],
                taumin[jp]
            );
        }
    }
    table
}

extern "C" fn closefiles(_sig: i32) {
    let files = TEMPFILES.load(std::sync::atomic::Ordering::Acquire);
    if !files.is_null() {
        // SAFETY: main leaks this allocation for the lifetime of the process.
        // Like the C program, this executable is single threaded.
        let files = unsafe { &*files };
        efclose(&files.headerfp);
        efclose(&files.tracefp);
        if !files.headerfile.is_empty() {
            eremove(&files.headerfile);
        }
        if !files.tracefile.is_empty() {
            eremove(&files.tracefile);
        }
    }
    exit(EXIT_FAILURE);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;
    use std::process::{Command, Stdio};

    // Compare the program-local ray tracer with the original diagnostics.
    // Only this test invokes the original executable, as a differential oracle.
    #[test]
    fn tableta_matches_original_diagnostics() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
        let bin = root.join("vendor/su-build/bin");
        let target = std::env::var_os("CARGO_TARGET_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| root.join(".targets/prog-sumigps"));
        std::fs::create_dir_all(&target).unwrap();
        let vfile = target.join(format!("sumigps-tableta-{}.bin", std::process::id()));

        for (nt, np, dt, first, step) in [
            (2, 5, 0.004f32, 1500.0f32, 0.0f32),
            (9, 11, 0.002, 2000.0, 0.0),
            (33, 17, 0.004, 1500.0, 16.0),
            (65, 31, 0.003, 800.0, 64.0),
            (47, 23, 0.004, 2500.0, -8.0),
            (129, 200, 0.001, 1500.0, 10.0),
        ] {
            let vt: Vec<f32> = (0..nt).map(|it| first + it as f32 * step).collect();
            let bytes: Vec<u8> = vt.iter().flat_map(|v| v.to_ne_bytes()).collect();
            std::fs::write(&vfile, bytes).unwrap();
            let input = Command::new(bin.join("suspike"))
                .args([format!("nt={nt}"), "ntr=1".to_owned()])
                .output()
                .unwrap();
            assert!(input.status.success());
            let mut original = Command::new(bin.join("sumigps"))
                .args([
                    "dx=20".to_owned(),
                    "verbose=1".to_owned(),
                    "tmpdir=".to_owned(),
                    format!("dt={dt}"),
                    format!("np={np}"),
                    format!("vfile={}", vfile.display()),
                ])
                .stdin(Stdio::piped())
                .stdout(Stdio::null())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            original
                .stdin
                .take()
                .unwrap()
                .write_all(&input.stdout)
                .unwrap();
            let output = original.wait_with_output().unwrap();
            assert!(output.status.success());
            let diagnostic = String::from_utf8(output.stderr).unwrap();
            let rows: Vec<&str> = diagnostic.lines().filter(|s| s.starts_with("p=")).collect();
            assert_eq!(rows.len(), np);
            let table = tableta(np, nt, dt, &vt, nt, dt, 0);
            // At the final slope, the first ray step can be negative after
            // rounding cos(pi/2); the C then prints uninitialized extrema.
            for jp in 0..np - 1 {
                let actual = sprintf!(
                    "p=%g tturn=%g taumax=%g taumin=%g",
                    table.p[jp],
                    table.tturn[jp],
                    table.taumax[jp],
                    table.taumin[jp]
                );
                assert_eq!(actual, rows[jp], "nt={nt} np={np} dt={dt} slope={jp}");
            }
        }
        std::fs::remove_file(vfile).unwrap();
    }
}
