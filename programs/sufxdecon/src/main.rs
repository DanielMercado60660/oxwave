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

//! SUFXDECON - random noise attenuation by FX-DECONvolution.
//! Port of su/main/decon_shaping/sufxdecon.c (Revision 1.10).
//! Credits: CWP: Carlos E. Theodoro (10/07/97).
//! References: Canales (1984), Gulunay (1986), Galbraith (1991).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cadd, cmplx, cmul, conjg, Complex};
use oxwave_cwp::cwp::nint;
use oxwave_f_axb::{backward_substitution, LU_decomposition};
use oxwave_pfa_fft::{npfaro, pfacr, pfarc};
use oxwave_tempfile::temporary_filename;
use std::ffi::{c_char, CString};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUFXDECON - random noise attenuation by FX-DECONvolution              ",
    "\t\t\t\t\t\t\t\t\t",
    " sufxdecon <stdin >stdout [...]\t                                ",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " taper=.1\tlength of taper                                         ",
    " fmin=6.       minimum frequency to process in Hz  (accord to twlen)   ",
    " fmax=.6/(2*dt)  maximum frequency to process in Hz                    ",
    " twlen=entire trace  time window length (minimum .3 for lower freqs)   ",
    " ntrw=10       number of traces in window                              ",
    " ntrf=4        number of traces for filter (smaller than ntrw)         ",
    " verbose=0\t=1 for diagnostic print\t\t\t\t\t",
    " tmpdir=\tif non-empty, use the value as a directory path\tprefix\t",
    "\t\tfor storing temporary files; else, if the CWP_TMPDIR\t",
    "\t\tenvironment variable is set, use its value for the path;",
    "\t\telse use tmpfile()\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes: Each trace is transformed to the frequency domain.             ",
    "        For each frequency, Wiener filtering, with unity prediction in ",
    "        space, is used to predict the next sample.                     ",
    "        At the end of the process, data is mapped back to t-x domain.  ",
    "\t\t\t\t\t\t\t\t\t",
    "                                                                       ",
];

const LOOKFAC: i32 = 2;
const PFA_MAX: i32 = 720720;

unsafe extern "C" {
    fn access(path: *const c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
}

// Globals in C, retained for graceful interrupt termination of named files.
struct TempFiles {
    tracefp: CFile,
    headerfp: CFile,
    tracefile: String,
    headerfile: String,
}
static mut TEMP_FILES: *const TempFiles = std::ptr::null();

extern "C" fn closefiles(_: i32) {
    // The box in main remains alive until program termination.
    unsafe {
        if let Some(files) = TEMP_FILES.as_ref() {
            efclose(&files.headerfp);
            efclose(&files.tracefp);
            eremove(&files.headerfile);
            eremove(&files.tracefile);
        }
    }
    exit(1);
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
    let nsp = tr.ns as i32;
    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    if verbose != 0 {
        warn!("nsp = %i", nsp);
    }
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    if dt == 0.0 {
        dt = 0.004;
        if verbose != 0 {
            warn!("dt not set, assumed to be .004");
        }
    }
    if verbose != 0 {
        warn!("dt = %f", dt);
    }
    let mut taper = 0.0f32;
    if getparfloat("taper", &mut taper) == 0 {
        taper = 0.1;
    }
    if taper == 0.0 {
        taper = 0.004;
    }
    let mut fmin = 0.0f32;
    if getparfloat("fmin", &mut fmin) == 0 {
        fmin = 6.0;
    }
    if fmin == 0.0 && verbose != 0 {
        warn!("using fmin=6 Hz");
    }
    let mut fmax = 0.0f32;
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (0.6 / (2.0f32 * dt) as f64) as f32;
    }
    if fmax == 0.0 && verbose != 0 {
        warn!("using fmax=.6/(2*dt) Hz");
    }
    let mut twlen = 0.0f32;
    if getparfloat("twlen", &mut twlen) == 0 {
        twlen = (nsp - 1) as f32 * dt;
    }
    if (twlen as f64) < 0.3 {
        twlen = 0.3;
        if verbose != 0 {
            warn!("twlen cannot be less than .3s, using .3s");
        }
    }

    // Setting taper and spatial and temporal windows.
    let mut nstaper = nint((taper / dt) as f64);
    nstaper = if nstaper % 2 != 0 {
        nstaper + 1
    } else {
        nstaper
    };
    let ntw = nint(((nsp - 1) as f32 * dt / twlen) as f64);
    if ntw == 1 {
        taper = 0.0;
    }
    let nspws = nint((twlen / dt) as f64) + 1;
    let nspwf = nspws + nstaper / 2;
    let nspwi = nspws + nstaper;
    let mut ntrw = 0i32;
    if getparint("ntrw", &mut ntrw) == 0 {
        ntrw = 10;
    }
    if ntrw == 0 {
        ntrw = 10;
        if verbose != 0 {
            warn!("ntrw cannot be zero, using 10 traces");
        }
    }
    if verbose != 0 {
        warn!("ntrw = %i", ntrw);
    }
    let mut ntrf = 0i32;
    if getparint("ntrf", &mut ntrf) == 0 {
        ntrf = 4;
    }
    if ntrf == 0 && verbose != 0 {
        warn!("using ntrf=4");
    }

    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let path = CString::new(tmpdir.as_bytes()).unwrap();
        if unsafe { access(path.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();

    let istmpdir = !tmpdir.is_empty();
    let files = if !istmpdir {
        let tracefp = etmpfile();
        let headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
        Box::new(TempFiles {
            tracefp,
            headerfp,
            tracefile: String::new(),
            headerfile: String::new(),
        })
    } else {
        let tracefile = temporary_filename(&tmpdir);
        let headerfile = temporary_filename(&tmpdir);
        let tracefp = efopen(&tracefile, "w+");
        let headerfp = efopen(&headerfile, "w+");
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
        Box::new(TempFiles {
            tracefp,
            headerfp,
            tracefile,
            headerfile,
        })
    };
    if istmpdir {
        unsafe {
            TEMP_FILES = &*files;
            for sig in [2, 3, 1, 15] {
                signal(sig, closefiles as *const () as usize);
            }
        }
    }

    let mut ntrc = 0i32;
    loop {
        ntrc += 1;
        efwrite(tr.hdr.as_bytes(), &files.headerfp);
        efwrite(&tr.data[..nsp as usize], &files.tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    erewind(&files.tracefp);
    erewind(&files.headerfp);
    if verbose != 0 {
        warn!("%i traces input", ntrc);
    }
    let nwx = ntrc / ntrw;
    let nfft = npfaro(2 * nspwf, LOOKFAC * 2 * nspwf);
    if nfft >= PFA_MAX {
        err!("Padded nsp=%d--too big", nfft);
    }
    let nf = nfft / 2 + 1;
    let d1 = (1.0 / (nfft as f32 * dt) as f64) as f32;

    // Flat trace arrays preserve C's contiguous ealloc2 storage.
    if ntrw < 0 || ntrf < 0 {
        // C converts negative allocation dimensions to size_t and ealloc fails.
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let zero = cmplx(0.0, 0.0);
    let mut ct = vec![zero; nfft as usize];
    let mut tidataw = vec![0.0f32; nfft as usize];
    let mut ttidataw = vec![0.0f32; nfft as usize];
    let mut fdata = vec![vec![zero; nf as usize]; ntrc as usize];
    let mut fdataw = vec![vec![zero; nf as usize]; (2 * ntrw + 2 * ntrf) as usize];
    // Rounded time-window counts can extend an intermediate window past nsp.
    // C then reads/writes the next contiguous trace, or overruns the allocation
    // on the last trace. One zeroed slack row reproduces its observed result
    // without out-of-bounds access (DESIGN rule 5); only real traces are read
    // from the temporary file and emitted below.
    let trace_samples = (nsp * ntrc) as usize;
    let mut tdatain = vec![0.0f32; trace_samples + nsp as usize];
    let mut tdataout = vec![0.0f32; trace_samples + nsp as usize];
    let mut sfreq = vec![zero; (2 * ntrw + 2 * ntrf) as usize];
    let mut sfreqcj = vec![zero; (2 * ntrw + 2 * ntrf) as usize];
    let mut autocorr = vec![zero; (ntrf + 1) as usize];
    let mut rautoc = vec![0.0f32; (ntrf + 1) as usize];
    let mut iautoc = vec![0.0f32; (ntrf + 1) as usize];
    let mut rmatrix = vec![vec![0.0f32; (2 * ntrf) as usize]; (2 * ntrf) as usize];
    let mut gvector = vec![0.0f32; (2 * ntrf) as usize];
    let mut fvector = vec![zero; (2 * ntrf + 1) as usize];
    // C leaves unprocessed bins uninitialized; fresh allocations are normally zero.
    let mut ffreq = vec![vec![zero; nf as usize]; (2 * ntrw) as usize];
    let mut sfreqout = vec![zero; (2 * ntrw) as usize];
    let mut freqv = vec![zero; nfft as usize];
    let mut todataw = vec![0.0f32; nfft as usize];
    // A final window shorter than nstaper starts the last overlap loop at a
    // negative it. C reads before ttodataw there. A zeroed prefix reproduces
    // the observed near-zero padding without forming an invalid index.
    let ttodataw_start = nstaper.max(0);
    let mut ttodataw = vec![0.0f32; (nfft + ttodataw_start) as usize];
    // C allocates 4*ntrw entries but LU writes 2*ntrf pivots. Provide slack
    // for oversized filters to preserve the solve without C's heap overrun.
    let solver_size = (4 * ntrw).max(2 * ntrf) as usize;
    let mut ipvt = vec![0i32; solver_size];
    let mut info = vec![0.0f32; solver_size];

    tdataout.fill(0.0);
    efread(&mut tdatain[..trace_samples], &files.tracefp);
    efclose(&files.tracefp);
    if dt != 0.0 && verbose != 0 {
        warn!("d1=%f", 1.0 / (nfft as f32 * dt) as f64);
    }
    if verbose != 0 {
        warn!("nf=%i, d1=%f, nfft=%i, nstaper=%i", nf, d1, nfft, nstaper);
    }

    // Loop over time windows.
    for iw in 0..ntw {
        let nspw = if iw > 0 && iw < ntw - 1 {
            nspwi
        } else if iw == 0 {
            if ntw > 1 {
                nspwf
            } else {
                nsp
            }
        } else {
            nsp - nspws * iw + nstaper / 2
        };
        if verbose != 0 {
            warn!("iw=%i, ntw=%i, nspw=%i, twlen=%f", iw, ntw, nspw, twlen);
        }
        for row in &mut fdata {
            row.fill(zero);
        }
        for jw in 0..ntrc {
            if iw > 0 {
                for it in 0..nspw {
                    tidataw[it as usize] =
                        tdatain[(jw * nsp + it + iw * nspws - nstaper / 2) as usize];
                }
            } else {
                for it in 0..nspw {
                    tidataw[it as usize] = tdatain[(jw * nsp + it) as usize];
                }
            }
            tidataw[nspw as usize..].fill(0.0);
            ct.fill(zero);
            for it in 0..nfft {
                ttidataw[it as usize] = if it % 2 != 0 {
                    -tidataw[it as usize]
                } else {
                    tidataw[it as usize]
                };
            }
            pfarc(1, &ttidataw, &mut ct);
            for ifq in 0..nf {
                fdata[jw as usize][ifq as usize] = ct[(nf - 1 - ifq) as usize];
            }
        }

        // Loop over space windows.
        for jx in 0..nwx {
            let ntrwu = if ntrc < jx * ntrw + 2 * ntrw {
                ntrc - jx * ntrw
            } else {
                ntrw
            };
            if verbose != 0 {
                warn!("jx=%i, ntrc=%i, ntrw=%i", jx, ntrc, ntrw);
                warn!("ntrwu=%i, nwx=%i, ntrf=%i", ntrwu, nwx, ntrf);
            }
            for jw in 0..ntrw + 2 * ntrf {
                fdataw[jw as usize].fill(zero);
            }
            for jw in 0..ntrwu + 2 * ntrf {
                for ifq in 0..nf {
                    let source = if jx > 0 && jx < nwx - 1 {
                        jw + jx * ntrw - ntrf
                    } else if jx == 0 {
                        if jw >= ntrf && jw < ntrw + ntrf {
                            jw - ntrf
                        } else if jw < ntrf {
                            0
                        } else if nwx > 1 {
                            jw - ntrf
                        } else {
                            ntrc - 1
                        }
                    } else if jw < ntrwu + ntrf {
                        jw + jx * ntrw - ntrf
                    } else {
                        ntrc - 1
                    };
                    fdataw[jw as usize][ifq as usize] = fdata[source as usize][ifq as usize];
                }
            }

            for ifq in 0..nf {
                if ifq as f32 * d1 >= fmin && ifq as f32 * d1 <= fmax {
                    sfreq[..(ntrwu + 2 * ntrf) as usize].fill(zero);
                    sfreqcj[..(ntrwu + 2 * ntrf) as usize].fill(zero);
                    for jw in 0..ntrwu + 2 * ntrf {
                        sfreq[jw as usize] = fdataw[jw as usize][ifq as usize];
                        sfreqcj[jw as usize] = conjg(fdataw[jw as usize][ifq as usize]);
                    }
                    autocorr.fill(zero);
                    cxcor(
                        ntrwu,
                        0,
                        &sfreq,
                        ntrwu,
                        0,
                        &sfreq,
                        ntrf + 1,
                        0,
                        &mut autocorr,
                    );
                    rautoc.fill(0.0);
                    iautoc.fill(0.0);
                    for jw in 0..ntrf + 1 {
                        rautoc[jw as usize] = autocorr[jw as usize].r;
                        iautoc[jw as usize] = autocorr[jw as usize].i;
                    }
                    gvector.fill(0.0);
                    fvector.fill(zero);
                    for row in &mut rmatrix {
                        row.fill(0.0);
                    }
                    for ir in 0..ntrf {
                        for jr in 0..ntrf {
                            rmatrix[ir as usize][jr as usize] = if ir >= jr {
                                autocorr[(ir - jr) as usize].r
                            } else {
                                autocorr[(jr - ir) as usize].r
                            };
                        }
                    }
                    for ir in ntrf..2 * ntrf {
                        for jr in 0..ntrf {
                            rmatrix[ir as usize][jr as usize] = if ir - ntrf < jr {
                                -autocorr[(jr - ir + ntrf) as usize].i
                            } else {
                                autocorr[(ir - jr - ntrf) as usize].i
                            };
                        }
                    }
                    for ir in ntrf..2 * ntrf {
                        for jr in ntrf..2 * ntrf {
                            rmatrix[ir as usize][jr as usize] =
                                rmatrix[(ir - ntrf) as usize][(jr - ntrf) as usize];
                        }
                    }
                    for ir in 0..ntrf {
                        for jr in ntrf..2 * ntrf {
                            rmatrix[ir as usize][jr as usize] =
                                -rmatrix[(ir + ntrf) as usize][(jr - ntrf) as usize];
                        }
                    }
                    for ig in 0..2 * ntrf {
                        gvector[ig as usize] = if ig < ntrf {
                            autocorr[(ig + 1) as usize].r
                        } else {
                            autocorr[(ig - ntrf + 1) as usize].i
                        };
                    }

                    // Zero-energy slices also reach the solver in C. The
                    // library emits singular-row diagnostics and continues,
                    // retaining the previous pivot when NaN comparisons fail.
                    // Keep this unconditional to preserve those diagnostics
                    // and the resulting nonfinite filter coefficients.
                    LU_decomposition(&mut rmatrix, &mut ipvt, &mut info[0]);
                    backward_substitution(&rmatrix, &ipvt, &mut gvector);

                    for ifv in 0..ntrf {
                        let ig = ntrf - 1 - ifv;
                        fvector[ifv as usize] = conjg(cmplx(
                            (gvector[ig as usize] as f64 / 2.0) as f32,
                            (gvector[(ig + ntrf) as usize] as f64 / 2.0) as f32,
                        ));
                    }
                    for ifv in ntrf + 1..2 * ntrf + 1 {
                        let ig = ifv - ntrf - 1;
                        fvector[ifv as usize] = cmplx(
                            (gvector[ig as usize] as f64 / 2.0) as f32,
                            (gvector[(ig + ntrf) as usize] as f64 / 2.0) as f32,
                        );
                    }
                    sfreqout[..ntrwu as usize].fill(zero);
                    cconv(
                        ntrwu + 2 * ntrf,
                        -ntrf,
                        &sfreq,
                        2 * ntrf + 1,
                        -ntrf,
                        &fvector,
                        ntrwu,
                        0,
                        &mut sfreqout,
                    );
                    for jw in 0..ntrwu {
                        ffreq[jw as usize][ifq as usize] = sfreqout[jw as usize];
                    }
                }
            }

            for jw in 0..ntrwu {
                for ifq in 0..nf {
                    freqv[ifq as usize] = ffreq[jw as usize][(nf - 1 - ifq) as usize];
                }
                freqv[nf as usize..].fill(zero);
                todataw.fill(0.0);
                pfacr(-1, &freqv, &mut todataw);
                for it in 0..nsp.min(nfft) {
                    todataw[it as usize] /= nfft as f32;
                }
                for it in 0..nsp.min(nfft) {
                    ttodataw[(it + ttodataw_start) as usize] = if it % 2 != 0 {
                        -todataw[it as usize]
                    } else {
                        todataw[it as usize]
                    };
                }
                let base = (jx * ntrw + jw) * nsp;
                if ntw > 1 {
                    if iw > 0 {
                        for it in 0..nstaper {
                            tdataout[(base + it + iw * nspws - nstaper / 2) as usize] +=
                                ttodataw[(it + ttodataw_start) as usize] * (it as f32 * dt / taper);
                        }
                    } else {
                        for it in 0..nstaper {
                            tdataout[(base + it) as usize] =
                                ttodataw[(it + ttodataw_start) as usize];
                        }
                    }
                    if iw > 0 {
                        for it in nstaper..nspw - nstaper {
                            tdataout[(base + it + iw * nspws - nstaper / 2) as usize] =
                                ttodataw[(it + ttodataw_start) as usize];
                        }
                    } else {
                        for it in nstaper..nspw - nstaper {
                            tdataout[(base + it) as usize] =
                                ttodataw[(it + ttodataw_start) as usize];
                        }
                    }
                    if iw > 0 && iw < ntw - 1 {
                        for it in nspw - nstaper..nspw {
                            let pos = (base + it + iw * nspws - nstaper / 2) as usize;
                            tdataout[pos] = (tdataout[pos] as f64
                                + ttodataw[(it + ttodataw_start) as usize] as f64
                                    * (1.0 - ((it - nspw + nstaper) as f32 * dt / taper) as f64))
                                as f32;
                        }
                    } else if iw == ntw - 1 {
                        for it in nspw - nstaper..nspw {
                            tdataout[(base + it + iw * nspws - nstaper / 2) as usize] =
                                ttodataw[(it + ttodataw_start) as usize];
                        }
                    } else {
                        for it in nspw - nstaper..nspw {
                            let pos = (base + it) as usize;
                            tdataout[pos] = (tdataout[pos] as f64
                                + ttodataw[(it + ttodataw_start) as usize] as f64
                                    * (1.0 - ((it - nspw + nstaper) as f32 * dt / taper) as f64))
                                as f32;
                        }
                    }
                } else {
                    for it in 0..nsp {
                        tdataout[(base + it) as usize] = ttodataw[(it + ttodataw_start) as usize];
                    }
                }
            }
        }
    }

    for jx in 0..ntrc {
        efread(tr.hdr.as_bytes_mut(), &files.headerfp);
        tr.data[..nsp as usize]
            .copy_from_slice(&tdataout[(jx * nsp) as usize..((jx + 1) * nsp) as usize]);
        puttr(&mut tr);
    }
    efclose(&files.headerfp);
    if istmpdir {
        eremove(&files.headerfile);
        eremove(&files.tracefile);
        unsafe {
            TEMP_FILES = std::ptr::null();
        }
    }
    exit(CWP_Exit());
}

// Complex correlation: local helper from sufxdecon.c.
#[allow(clippy::too_many_arguments)]
fn cxcor(
    lx: i32,
    ifx: i32,
    x: &[Complex],
    ly: i32,
    ify: i32,
    y: &[Complex],
    lz: i32,
    ifz: i32,
    z: &mut [Complex],
) {
    let mut xr = vec![Complex::default(); lx as usize];
    for i in 0..lx {
        xr[i as usize] = conjg(x[(lx - 1 - i) as usize]);
    }
    cconv(lx, 1 - ifx - lx, &xr, ly, ify, y, lz, ifz, z);
}

// Offsets implement C's pointer shifts without forming out-of-bounds pointers.
#[allow(clippy::too_many_arguments)]
fn cconv(
    lx: i32,
    ifx: i32,
    x: &[Complex],
    ly: i32,
    ify: i32,
    y: &[Complex],
    lz: i32,
    ifz: i32,
    z: &mut [Complex],
) {
    let ilx = ifx + lx - 1;
    let ily = ify + ly - 1;
    let ilz = ifz + lz - 1;
    for i in ifz..=ilz {
        let mut jlow = i - ily;
        if jlow < ifx {
            jlow = ifx;
        }
        let mut jhigh = i - ify;
        if jhigh > ilx {
            jhigh = ilx;
        }
        let mut sum = cmplx(0.0, 0.0);
        for j in jlow..=jhigh {
            sum = cadd(sum, cmul(x[(j - ifx) as usize], y[(i - j - ify) as usize]));
        }
        z[(i - ifz) as usize] = sum;
    }
}
