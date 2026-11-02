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
//! SUSPECFK - F-K Fourier spectrum. Port of su/main/transforms/suspecfk.c (Revision 1.21).
//! Credits: CWP: Dave (algorithm), Jack (reformatting for SU).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::rcabs;
use oxwave_pfa_fft::{npfa, npfar, pfa2rc, pfa2cc, Complex};
use oxwave_tempfile::temporary_filename;
use std::ffi::CString;
use std::sync::atomic::{AtomicPtr, Ordering};
unsafe extern "C" {
    fn access(path: *const i8, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn unlink(path: *const i8) -> i32;
    fn _exit(status: i32) -> !;
}

// The CString in main owns this path for the lifetime of the signal handlers.
static TRACEFILE: AtomicPtr<i8> = AtomicPtr::new(std::ptr::null_mut());

/* for graceful interrupt termination */
extern "C" fn closefiles(_: i32) {
    let tracefile = TRACEFILE.load(Ordering::Relaxed);
    unsafe {
        if !tracefile.is_null() { unlink(tracefile); }
        // Termination closes the removed temporary file. Avoid buffered Rust
        // I/O here: the signal may have interrupted a borrow of that stream.
        _exit(1);
    }
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUSPECFK - F-K Fourier SPECtrum of data set\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " suspecfk <infile >outfile [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " dt=from header\t\ttime sampling interval\t\t",
    " dx=from header(d2) or 1.0\tspatial sampling interval\t",
    "\t\t\t\t\t\t\t\t",
    " verbose=0\tverbose = 1 echoes information\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path",
    "\t\t prefix for storing temporary files; else if the",
    "\t         the CWP_TMPDIR environment variable is set use\t",
    "\t         its value for the path; else use tmpfile()\t",
    " \t\t\t\t\t\t\t\t",
    " Note: To facilitate further processing, the sampling intervals",
    "       in frequency and wavenumber as well as the first\t",
    "\tfrequency (0) and the first wavenumber are set in the\t",
    "\toutput header (as respectively d1, d2, f1, f2).\t\t",
    " \t\t\t\t\t\t\t\t",
    " Note: The relation: w = 2 pi F is well known, but there\t",
    "\tdoesn't\tseem to be a commonly used letter corresponding\t",
    "\tto F for the spatial conjugate transform variable.  We\t",
    "\tuse K for this.  More specifically we assume a phase:\t",
    "\t\ti(w t - k x) = 2 pi i(F t - K x).\t\t",
    "\tand F, K define our notion of frequency, wavenumber.\t",
    " \t\t\t\t\t\t\t\t",
];

fn workspace_len(n1: usize, n2: usize) -> usize {
    n1.checked_mul(n2)
        .unwrap_or_else(|| err!("FFT workspace size overflow"))
}

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut intrace) == 0 { err!("can't get first trace"); }
    let nt = intrace.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        if intrace.dt != 0 { dt = (intrace.dt as f64 / 1000000.0) as f32; }
        else { dt = 0.004; warn!("tr.dt not set, assuming dt=0.004"); }
    }
    let mut dx = 0.0f32;
    if getparfloat("dx", &mut dx) == 0 {
        if intrace.d2 != 0.0 { dx = intrace.d2; }
        else { dx = 1.0; warn!("tr.d2 not set, assuming d2=1.0"); }
    }
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
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
    checkpars();
    let mut tracefile = String::new();
    let signal_path;
    let tracefp = if tmpdir.is_empty() {
        let fp = etmpfile();
        if verbose != 0 { warn!("using tmpfile() call"); }
        fp
    } else {
        tracefile = temporary_filename(&tmpdir);
        signal_path = CString::new(tracefile.as_str()).unwrap();
        TRACEFILE.store(signal_path.as_ptr().cast_mut(), Ordering::Relaxed);
        /* Trap signals so can remove temp files */
        unsafe {
            for sig in [2, 3, 1, 15] { signal(sig, closefiles); }
        }
        let fp = efopen(&tracefile, "w+");
        if verbose != 0 { warn!("putting temporary files in %s", tmpdir); }
        fp
    };
    let mut nx = 0i32;
    loop {
        nx += 1;
        efwrite(&intrace.data[..nt as usize], &tracefp);
        if gettr(&mut intrace) == 0 { break; }
    }
    let ntfft = npfar(nt);
    let nxfft = npfa(nx);
    if ntfft >= SU_NFLTS as i32 || ntfft >= 720720 { err!("Padded nt=%d--too big", ntfft); }
    if nxfft >= SU_NFLTS as i32 || nxfft >= 720720 { err!("Padded nx=%d--too big", nxfft); }
    let d1 = (1.0 / ((ntfft as f32 * dt) as f64)) as f32;
    let d2 = (1.0 / ((nxfft as f32 * dx) as f64)) as f32;
    let f2 = (-1.0 / ((2.0f32 * dx) as f64)) as f32;
    let nf = ntfft / 2 + 1;
    let nk = nxfft;
    // C alloc2 uses size_t dimensions; the product can exceed i32::MAX.
    let ct_len = workspace_len(nf as usize, nk as usize);
    let rt_len = workspace_len(ntfft as usize, nxfft as usize);
    let mut ct = Vec::new();
    if ct.try_reserve_exact(ct_len).is_err() {
        err!("cannot allocate complex FFT workspace");
    }
    ct.resize(ct_len, Complex::default());
    let mut rt = Vec::new();
    if rt.try_reserve_exact(rt_len).is_err() {
        err!("cannot allocate real FFT workspace");
    }
    rt.resize(rt_len, 0.0f32);
    erewind(&tracefp);
    for ix in 0..nx as usize {
        let base = ix * ntfft as usize;
        efread(&mut rt[base..base + nt as usize], &tracefp);
        if ix % 2 != 0 {
            for it in 0..nt as usize { rt[base + it] = -rt[base + it]; }
        }
    }
    efclose(&tracefp);
    pfa2rc(1, 1, ntfft, nx, &rt, &mut ct);
    pfa2cc(-1, 2, nf, nxfft, &mut ct);
    for ik in 0..nk as usize {
        for i_f in 0..nf as usize { outtrace.data[i_f] = rcabs(ct[ik * nf as usize + i_f]); }
        outtrace.tracl = ik as i32 + 1;
        outtrace.ns = nf as u16;
        outtrace.dt = 0;
        outtrace.trid = 122;
        outtrace.d1 = d1;
        outtrace.f1 = 0.0;
        outtrace.d2 = d2;
        outtrace.f2 = f2;
        puttr(&mut outtrace);
    }
    if !tracefile.is_empty() { eremove(&tracefile); }
    exit(CWP_Exit());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn padded_workspace_exceeds_signed_int_range() {
        // Both are valid PFA dimensions below SU_NFLTS. Check the actual
        // length calculation without allocating a multi-gigabyte buffer.
        let ntfft = npfar(60060);
        let nxfft = npfa(60060);
        assert_eq!(ntfft, 60060);
        assert_eq!(nxfft, 60060);
        let len = workspace_len(ntfft as usize, nxfft as usize);
        assert_eq!(len, 3_607_203_600);
        assert!(len > i32::MAX as usize);
        assert_eq!(workspace_len((ntfft / 2 + 1) as usize, nxfft as usize), 1_803_661_860);
    }
}
