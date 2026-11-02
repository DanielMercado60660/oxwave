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

//! SUTAUP - forward and inverse T-X and F-K global slant stacks.
//! Port of su/main/transforms/sutaup.c (Revision 1.15).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: CWP: Gabriel Alvarez, 1995.
//! Reference: Levin, F., editor, 1991, Slant-Stack Processing,
//! Geophysics Reprint Series #14, SEG Press, Tulsa.

use oxwave_core::prelude::*;
use oxwave_f_taup::{fwd_FK_sstack, fwd_tx_sstack, inv_FK_sstack, inv_tx_sstack};
use oxwave_tempfile::temporary_filename;
use std::ffi::{c_char, CString};
use std::sync::atomic::{AtomicBool, AtomicPtr, Ordering};

unsafe extern "C" {
    fn access(path: *const c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn unlink(path: *const c_char) -> i32;
    fn _exit(status: i32) -> !;
}

// The C globals allow closefiles to remove named files on interruption.
static INTERRUPTED: AtomicBool = AtomicBool::new(false);
static READING_INPUT: AtomicBool = AtomicBool::new(false);
static TRACEFILE: AtomicPtr<c_char> = AtomicPtr::new(std::ptr::null_mut());
static HEADERFILE: AtomicPtr<c_char> = AtomicPtr::new(std::ptr::null_mut());

extern "C" fn note_signal(_: i32) {
    INTERRUPTED.store(true, Ordering::Relaxed);
    if READING_INPUT.load(Ordering::Relaxed) {
        // gettr may be blocked while borrowing a CFile. There is no output
        // yet, so discard temporary buffers and use only async-signal-safe
        // operations. The immutable CStrings live until main exits; process
        // exit closes the descriptors after unlinking both files.
        unsafe {
            unlink(HEADERFILE.load(Ordering::Relaxed));
            unlink(TRACEFILE.load(Ordering::Relaxed));
            _exit(EXIT_FAILURE);
        }
    }
}

fn closefiles(headerfp: &CFile, tracefp: Option<&CFile>, headerfile: &str, tracefile: &str) {
    if !INTERRUPTED.load(Ordering::Relaxed) {
        return;
    }
    // Outside input gathering, defer stream operations until the interrupted
    // CFile borrow has ended. Do not close an already closed trace stream.
    efclose(headerfp);
    if let Some(fp) = tracefp {
        efclose(fp);
    }
    eremove(headerfile);
    if !tracefile.is_empty() {
        eremove(tracefile);
    }
    exit(EXIT_FAILURE);
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUTAUP - forward and inverse T-X and F-K global slant stacks\t\t",
    "                                                                       ",
    "    sutaup <infile >outfile  [optional parameters]                 \t",
    "                                                                       ",
    " Optional Parameters:                                                  ",
    " option=1\t\t\t=1 for forward F-K domian computation\t",
    "\t\t\t\t=2 for forward T-X domain computation\t",
    "\t\t\t\t=3 for inverse F-K domain computation\t",
    "\t\t\t\t=4 for inverse T-X domain computation\t",
    " dt=tr.dt (from header) \ttime sampling interval (secs)           ",
    " nx=ntr   (counted from data)\tnumber of horizontal samples (traces)\t",
    " dx=1\t\t\t\thorizontal sampling interval (m)\t",
    " npoints=71\t\t\tnumber of points for rho filter\t\t",
    " pmin=0.0\t\t\tminimum slope for Tau-P transform (s/m)\t",
    " pmax=.006\t\t\tmaximum slope for Tau-P transform (s/m)\t",
    " np=nx\t\t\t\tnumber of slopes for Tau-P transform\t",
    " ntau=nt\t\t\tnumber of time samples in Tau-P domain  ",
    " fmin=3\t\t\tminimum frequency of interest \t        ",
    " xmin=0\t\t\toffset on first trace\t \t        ",
    "                                                                       ",
    " verbose=0\tverbose = 1 echoes information\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:                                                                ",
    " The cascade of a forward and inverse  tau-p transform preserves the\t",
    " relative amplitudes in a data panel, but not the absolute amplitudes  ",
    " meaning that a scale factor must be applied to data output by such a  ",
    " a cascade before the output may be compared to the original data.\t",
    " This is a characteristic of the algorithm employed in this program.\t",
    " (Suradon does not have this problem.)\t\t\t\t\t",
    "                                                                       ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();

    /* hook up getpar to handle the parameters */
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* Look for user-supplied tmpdir */
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let path = CString::new(tmpdir.as_str()).unwrap();
        if unsafe { access(path.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }

    /* get info from first trace */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let dtheader = (tr.dt as f64 / 1000000.0) as f32;

    /* Store traces in tmpfile while getting a count */
    let mut tracefile = String::new();
    let mut headerfile = String::new();
    let tracepath;
    let headerpath;
    let istmpdir = !tmpdir.is_empty();
    let (tracefp, headerfp) = if !istmpdir {
        let tracefp = etmpfile();
        let headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
        (tracefp, headerfp)
    } else {
        tracefile = temporary_filename(&tmpdir);
        headerfile = temporary_filename(&tmpdir);
        tracepath = CString::new(tracefile.as_str()).unwrap();
        headerpath = CString::new(headerfile.as_str()).unwrap();
        TRACEFILE.store(tracepath.as_ptr().cast_mut(), Ordering::Relaxed);
        HEADERFILE.store(headerpath.as_ptr().cast_mut(), Ordering::Relaxed);
        READING_INPUT.store(true, Ordering::Relaxed);
        // SIGINT, SIGQUIT, SIGHUP, SIGTERM, in the same order as C.
        for sig in [2, 3, 1, 15] {
            unsafe {
                signal(sig, note_signal);
            }
        }
        let tracefp = efopen(&tracefile, "w+");
        let headerfp = efopen(&headerfile, "w+");
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
        (tracefp, headerfp)
    };
    let mut ntr = 0i32;
    loop {
        ntr += 1;
        efwrite(tr.hdr.as_bytes(), &headerfp);
        efwrite(&tr.data[..nt as usize], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    READING_INPUT.store(false, Ordering::Relaxed);
    closefiles(&headerfp, Some(&tracefp), &headerfile, &tracefile);

    /* get general flags and parameters and set defaults */
    let mut np = 0i32;
    let mut nx = 0i32;
    let mut option = 0i32;
    let mut pmin = 0.0f32;
    let mut pmax = 0.0f32;
    let mut xmin = 0.0f32;
    let mut dx = 0.0f32;
    let mut dt = 0.0f32;
    let mut fmin = 0.0f32;
    let mut npoints = 0i32;
    if getparint("np", &mut np) == 0 {
        np = ntr;
    }
    if getparint("nx", &mut nx) == 0 {
        nx = ntr;
    }
    if getparint("option", &mut option) == 0 {
        option = 1;
    }
    if getparfloat("pmin", &mut pmin) == 0 {
        pmin = 0.0;
    }
    if getparfloat("pmax", &mut pmax) == 0 {
        pmax = 0.006;
    }
    if getparfloat("xmin", &mut xmin) == 0 {
        xmin = 0.0;
    }
    if getparfloat("dx", &mut dx) == 0 {
        dx = 1.0;
    }
    if getparfloat("dt", &mut dt) == 0 {
        dt = dtheader;
    }
    if getparfloat("fmin", &mut fmin) == 0 {
        fmin = 3.0;
    }
    if getparint("npoints", &mut npoints) == 0 {
        npoints = 71;
    }

    checkpars();
    if dt == 0.0 {
        err!("header field dt not set, must be getparred");
    }

    /* allocate space */
    let mut in_traces = vec![vec![0.0f32; nt as usize]; ntr as usize];
    let mut out_traces = vec![vec![0.0f32; nt as usize]; np as usize];
    let dp = (pmax - pmin) / np.wrapping_sub(1) as f32;

    /* load traces into an array and close temp file */
    erewind(&headerfp);
    erewind(&tracefp);
    for ix in 0..ntr as usize {
        fread(&mut in_traces[ix], &tracefp);
    }
    closefiles(&headerfp, Some(&tracefp), &headerfile, &tracefile);
    efclose(&tracefp);
    if istmpdir {
        eremove(&tracefile);
        tracefile.clear();
    }

    /* do requested operation */
    if option == 1 {
        /* compute F-K forward slant stack */
        fwd_FK_sstack(
            dt,
            nt as usize,
            xmin,
            dx,
            pmin,
            dp,
            fmin,
            // A negative C nx contributes no input traces. The F-K library
            // also needs to support an empty slice for this signed count.
            &in_traces[..if nx > 0 { nx as usize } else { 0 }],
            &mut out_traces,
        );
    } else if option == 2 {
        /* compute t-x forward slant stack */
        // C starts its slope loop at pmin without multiplying by dp. With
        // np=1, dp is infinite (or NaN), but is never used for a sample.
        // The Rust API evaluates pmin + ip*dp even for ip=0, so supply the
        // unused step as zero for this case; keep the original dp in tr.d2.
        let stack_dp = if np == 1 { 0.0 } else { dp };
        fwd_tx_sstack(
            dt,
            nt as usize,
            xmin,
            dx,
            pmin,
            stack_dp,
            // A negative C nx runs no iterations of the input-trace loop.
            &in_traces[..if nx > 0 { nx as usize } else { 0 }],
            &mut out_traces,
        );
    } else if option == 3 {
        /* compute F-K inverse slant stack */
        inv_FK_sstack(
            dt,
            nt as usize,
            xmin,
            dx,
            pmin,
            dp,
            fmin,
            &in_traces[..np as usize],
            &mut out_traces[..nx as usize],
        );
    } else if option == 4 {
        /* compute t-x inverse slant stack */
        // The Rust rho_filter API asserts its preconditions; retain the C
        // diagnostics here before entering the library (same validation order).
        let nph = npoints / 2;
        if 2 * nph == npoints {
            err!("npoints must be odd!\n");
        }
        if nph > nt {
            err!("filter length larger than number of time samples!");
        }
        if nt > 720720 {
            err!("number of time samples too large for npfa!");
        }
        inv_tx_sstack(
            dt,
            nt as usize,
            npoints as usize,
            xmin,
            dx,
            pmin,
            dp,
            &in_traces[..np as usize],
            &mut out_traces[..nx as usize],
        );
    } else {
        err!("option flag has to be between 1 and 4");
    }

    /* write output traces */
    closefiles(&headerfp, None, &headerfile, &tracefile);
    erewind(&headerfp);
    for itr in 0..np as usize {
        closefiles(&headerfp, None, &headerfile, &tracefile);
        efread(tr.hdr.as_bytes_mut(), &headerfp);
        tr.tracl = 1 + itr as i32;
        tr.tracr = 1 + itr as i32;
        tr.dt = (dt as f64 * 1000000.0) as i32 as u16;
        tr.ns = nt as u16;
        if option == 1 || option == 2 {
            tr.d2 = dp;
            tr.f2 = pmin;
        } else {
            tr.d2 = 0.0;
            tr.f2 = 0.0;
        }
        for it in 0..nt as usize {
            tr.data[it] = out_traces[itr][it];
        }
        puttr(&mut tr);
    }
    closefiles(&headerfp, None, &headerfile, &tracefile);
    efclose(&headerfp);
    if istmpdir {
        eremove(&headerfile);
    }
    exit(CWP_Exit());
}
