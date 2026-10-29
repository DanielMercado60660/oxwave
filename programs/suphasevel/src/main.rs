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

//! SUPHASEVEL - Multi-mode phase velocity dispersion map from shot records.
//! Port of su/main/transforms/suphasevel.c (Revision 1.5).
//! Original copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: UHouston: Chris Liner June2008 (cloned from suspecfk).
//! Reference: Park, Miller, and Xia (1998, SEG Abstracts).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cadd, cmplx, cmul, crmul, cwp_cexp, rcabs};
use oxwave_cwp::cwp::PI;
use oxwave_pfa_fft::*;
use oxwave_tempfile::*;
use std::ffi::CString;
use std::sync::atomic::{AtomicBool, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUPHASEVEL - Multi-mode PHASE VELocity dispersion map computed",
    "              from shot record(s)\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " suphasevel <infile >outfile [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " fv=330\tminimum phase velocity (m/s)\t\t\t",
    " nv=100\tnumber of phase velocities\t\t\t",
    " dv=25\t\tphase velocity step (m/s)\t\t\t",
    " fmax=50\tmaximum frequency to process (Hz)\t\t",
    "\t\t=0 process to nyquist\t\t\t\t",
    " norm=0\tdo not normalize by amplitude spectrum\t\t",
    "\t\t=1 normalize by amplitude spectrum\t\t",
    " verbose=0\tverbose = 1 echoes information\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Notes:  Offsets read from headers.\t\t\t \t",
    "  1. output is complex frequency data\t\t\t\t",
    "  2. offset header word must be set (signed offset is ok)\t",
    "  3. norm=1 tends to blow up aliasing and other artifacts\t",
    "  4. For correct suifft later, use fmax=0\t\t\t",
    "  5. For later processing outtrace.dt=domega\t\t\t",
    "  6. works for 2D or 3D shots in any offset order\t\t",
    "\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Using this program:\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " First: use \t\t\t\t\t\t\t",
    " \tsuspecfx < shotrecord.su | suximage\t\t\t",
    " to see what the maximum bandwidth is in your data. This will\t",
    " give you an idea about the possible value for fmax.\t\t",
    "\t\t\t\t\t\t\t\t",
    " Second: Plot your data or some subset of your data via:\t",
    "     suxwigb < shotrecord.su key=offset\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " You can then estimate the range of phase velocities by looking",
    " at the maximum and minimum slopes of arrivals in your data.\t",
    " This will allow you do set first velocity fv and the increment",
    " in velocity dv, that make sense for your data.\t\t",
    " \t\t\t\t\t\t\t\t",
    " You can pick values of offset and time by placing the cursor  ",
    " on the desired location on the plot and pressing the 's' key",
    " The picks will appear in your terminal window. \t\t",
    "\t\t\t\t\t\t\t\t",
    " When displaying, don't forget to use suamp to compute the\t",
    " modulus of the complex values that this program puts out.\t",
    "\t\t\t\t\t\t\t\t",
    "   suphasevel < shotrecord.su [parameters] | suamp | suximage\t",
    "\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

const PFA_MAX: i32 = 720720;
const MAX_OFFS: usize = 10000;
const FUNPACKNYQ: i16 = 111;

extern "C" {
    fn access(path: *const std::ffi::c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

struct TempFile {
    tracefile: String,
    tracefp: Option<CFile>,
}
static INTERRUPTED: AtomicBool = AtomicBool::new(false);

// AtomicBool stores are lock-free on the supported targets. The handler must
// neither borrow main's file state nor call stdio, allocation, or exit routines.
extern "C" fn note_signal(_sig: i32) {
    INTERRUPTED.store(true, Ordering::Relaxed);
}

impl TempFile {
    // Perform the C closefiles helper's cleanup only in normal control flow,
    // between I/O calls, and do not close a stream that main already closed.
    fn check_signal(&mut self) {
        if INTERRUPTED.load(Ordering::Relaxed) {
            self.close();
            exit(EXIT_FAILURE);
        }
    }

    fn close(&mut self) {
        if let Some(fp) = self.tracefp.take() {
            efclose(&fp);
        }
        if !self.tracefile.is_empty() {
            eremove(&self.tracefile);
            self.tracefile.clear();
        }
    }
}

#[allow(non_snake_case)]
fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get info from first trace. */
    if gettr(&mut intrace) == 0 {
        err!("can't get first trace");
    }
    let nt = intrace.ns as usize;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        if intrace.dt != 0 {
            dt = (intrace.dt as f64 / 1000000.0) as f32;
        } else {
            err!("tr.dt not set, stop.");
        }
    }
    warn!("dt=%f", dt);
    let (mut fv, mut dv) = (330.0f32, 25.0f32);
    let (mut nv, mut out, mut norm) = (100i32, 0i32, 0i32);
    let mut fmax = 50.0f32;
    let mut verbose = 0i32;
    if getparfloat("fv", &mut fv) == 0 {
        fv = 330.0;
    }
    if getparfloat("dv", &mut dv) == 0 {
        dv = 25.0;
    }
    if getparint("nv", &mut nv) == 0 {
        nv = 100;
    }
    if getparint("out", &mut out) == 0 {
        out = 0;
    }
    let _ = out; // Parsed but unused in the original program.
    if getparint("norm", &mut norm) == 0 {
        norm = 0;
    }
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = 50.0;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* Look for user-supplied tmpdir. */
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let directory = CString::new(tmpdir.as_str()).unwrap();
        // SAFETY: path is NUL terminated; WRITE_OK is POSIX W_OK (2).
        if unsafe { access(directory.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();

    /* Set up tmpfile. */
    let mut files = TempFile {
        tracefile: String::new(),
        tracefp: None,
    };
    let istmpdir = !tmpdir.is_empty();
    let tracefp;
    if !istmpdir {
        tracefp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        files.tracefile = temporary_filename(&tmpdir);
        for sig in [2, 3, 1, 15] {
            // SAFETY: these POSIX signals accept a handler with this ABI; it
            // only stores to the lock-free flag above.
            unsafe {
                signal(sig, note_signal);
            }
        }
        tracefp = efopen(&files.tracefile, "w+");
        files.tracefp = Some(tracefp.clone());
        files.check_signal();
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
    }

    /* We have to allocate offs(nx) before we know nx. */
    let mut offs = vec![0.0f32; MAX_OFFS];
    let mut ix = 0usize;
    let mut nx = 0i32;
    let mut xmax = 0.0f32;
    loop {
        files.check_signal();
        // The C fixed-size offset array overflows on the 10,001st trace.
        if ix == MAX_OFFS {
            files.close();
            err!("more than %d input traces, stop.", MAX_OFFS);
        }
        nx += 1;
        efwrite(&intrace.data[..nt], &tracefp);
        offs[ix] = intrace.offset as f32;
        if intrace.offset.wrapping_abs() as f32 > xmax {
            xmax = intrace.offset.wrapping_abs() as f32;
        }
        ix += 1;
        if gettr(&mut intrace) == 0 {
            break;
        }
    }
    files.check_signal();
    if xmax == 0.0 {
        err!("tr.offset not set, stop.");
    }

    /* Determine lengths for prime-factor FFTs. */
    let ntfft = npfar(nt as i32);
    if ntfft >= SU_NFLTS as i32 || ntfft >= PFA_MAX {
        err!("Padded nt=%d--too big", ntfft);
    }
    let mut nF = ntfft / 2 + 1;
    let fft_nF = nF as usize;
    let onfft = (1.0 / ntfft as f64) as f32;
    let twopi = (2.0 * PI) as f32;
    let domega = twopi * onfft / dt;

    /* Allocate space, load traces and close tmpfile. */
    let mut ct = vec![Complex::default(); fft_nF * nx as usize];
    let mut rt = vec![0.0f32; ntfft as usize * nx as usize];
    erewind(&tracefp);
    for ix in 0..nx as usize {
        files.check_signal();
        let start = ix * ntfft as usize;
        efread(&mut rt[start..start + nt], &tracefp);
        for it in nt..ntfft as usize {
            rt[start + it] = 0.0;
        }
    }
    efclose(&tracefp);
    files.tracefp = None;
    files.check_signal();
    pfa2rc(1, 1, ntfft, nx, &rt, &mut ct);
    files.check_signal();

    /* Set nF for processing. */
    if fmax == 0.0 {
        nF = ntfft / 2 + 1;
    } else {
        nF = (twopi * fmax / domega) as i32;
    }
    // Check the actual truncated bin count, preserving safe fmax values just
    // above Nyquist. The C reads past ct when processing more FFT bins, and
    // writes past outtrace.data if twice the count exceeds its capacity.
    if (norm == 1 || nv > 0) && nF > fft_nF as i32 {
        files.close();
        err!(
            "fmax selects %d frequency bins, only %d available, stop.",
            nF,
            fft_nF
        );
    }
    if nv > 0 && nF > (SU_NFLTS / 2) as i32 {
        files.close();
        err!(
            "output trace needs %d samples, maximum is %d, stop.",
            nF as i64 * 2,
            SU_NFLTS
        );
    }
    // A negative C allocation fails, but no iteration dereferences it.
    let mut cDisp = vec![Complex::default(); nF.max(0) as usize];
    if norm == 1 {
        for iF in 0..nF {
            files.check_signal();
            let _omega = iF as f32 * domega;
            for ix in 0..nx as usize {
                let index = ix * fft_nF + iF as usize;
                let amp = rcabs(ct[index]);
                let oamp = (1.0 / amp as f64) as f32;
                ct[index] = crmul(ct[index], oamp);
            }
        }
    }

    /* Set global output trace headers. */
    outtrace.ns = nF.wrapping_mul(2) as u16;
    // arm64 converts the double to unsigned int before narrowing to ushort.
    outtrace.dt = (dt as f64 * 1000000.0) as u32 as u16;
    outtrace.trid = FUNPACKNYQ;
    outtrace.d1 = (1.0 / (ntfft as f32 * dt) as f64) as f32;
    outtrace.f1 = 0.0;
    outtrace.d2 = dv;
    outtrace.f2 = fv;

    /* Loop over phase velocities. */
    for iv in 0..nv {
        files.check_signal();
        let v = fv + iv as f32 * dv;
        for iF in 0..nF {
            files.check_signal();
            let omega = iF as f32 * domega;
            let phi = omega / v;
            let iF = iF as usize;
            cDisp[iF] = cmplx(0.0, 0.0);
            /* Sum over abs offset (this is ok for 3D, too). */
            for ix in 0..nx as usize {
                let x = offs[ix].abs();
                let arg = -phi * x;
                let cExp = cwp_cexp(crmul(cmplx(0.0, 1.0), arg));
                cDisp[iF] = cadd(cDisp[iF], cmul(ct[ix * fft_nF + iF], cExp));
            }
        }
        outtrace.tracl = iv + 1;
        for iF in 0..nF {
            let iF = iF as usize;
            outtrace.data[2 * iF] = cDisp[iF].r;
            outtrace.data[2 * iF + 1] = cDisp[iF].i;
        }
        files.check_signal();
        puttr(&mut outtrace);
    }
    files.check_signal();
    files.close();
    exit(CWP_Exit());
}
