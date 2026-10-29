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
//! SUSPECK1K2 - 2D spatial Fourier amplitude spectrum.
//! Port of su/main/transforms/suspeck1k2.c (Revision 1.13).
//! Credits: CWP: John Stockwell, based on code by Dave Hale and Jack Cohen.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{Complex, rcabs};
use oxwave_pfa_fft::{npfa, npfar, pfa2cc, pfa2rc};
use oxwave_tempfile::temporary_filename;
use std::ffi::CString;

unsafe extern "C" {
    fn access(path: *const i8, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn unlink(path: *const i8) -> i32;
    fn _exit(status: i32) -> !;
}

static mut TRACEFILE: *const i8 = std::ptr::null();

extern "C" fn closefiles(_: i32) {
    unsafe {
        if !TRACEFILE.is_null() { unlink(TRACEFILE); }
        _exit(1);
    }
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUSPECK1K2 - 2D (K1,K2) Fourier SPECtrum of (x1,x2) data set\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suspeck1k2 <infile >outfile [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " d1=from header(d1) or 1.0\tspatial sampling interval in first (fast)",
    "\t\t\t\t   dimension\t\t\t\t",
    " d2=from header(d2) or 1.0\tspatial sampling interval in second\t",
    "\t\t\t\t (slow)  dimension\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " verbose=0\t\tverbose = 1 echoes information\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t \tif non-empty, use the value as a directory path",
    "\t\t \tprefix for storing temporary files; else if the",
    "\t         \tthe CWP_TMPDIR environment variable is set use\t",
    "\t         \tits value for the path; else use tmpfile()\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Because the data are assumed to be purely spatial (i.e. non-seismic), ",
    " the data are assumed to have trace id (30), corresponding to (z,x) data",
    "\t\t\t\t\t\t\t\t\t",
    " To facilitate further processing, the sampling intervals in wavenumber",
    " as well as the first frequency (0) and the first wavenumber are set in",
    " the output header (as respectively d1, d2, f1, f2).\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The relation: w = 2 pi F is well known for frequency, but there\t",
    " doesn't seem to be a commonly used letter corresponding to F for the\t",
    " spatial conjugate transform variables.  We use K1 and K2 for this.\t",
    " More specifically we assume a phase:\t\t\t\t\t",
    "\t\t-i(k1 x1 + k2 x2) = -2 pi i(K1 x1 + K2 x2).\t\t",
    " and K1, K2 define our respective wavenumbers.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    if gettr(&mut intrace) == 0 { err!("can't get first trace"); }
    if intrace.trid != 30 { warn!("tr.trid = %d", intrace.trid); }
    let nx1 = intrace.ns as i32;
    let mut dx1 = 0.0f32;
    let mut dx2 = 0.0f32;
    if getparfloat("d1", &mut dx1) == 0 {
        if intrace.d1 != 0.0 { dx1 = intrace.d1; }
        else { dx1 = 1.0; warn!("tr.d1 not set, assuming d1=1.0"); }
    }
    if getparfloat("d2", &mut dx2) == 0 {
        if intrace.d2 != 0.0 { dx2 = intrace.d2; }
        else { dx2 = 1.0; warn!("tr.d2 not set, assuming d2=1.0"); }
    }
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() && unsafe { access(CString::new(tmpdir.as_str()).unwrap().as_ptr(), 2) } != 0 {
        err!("you can't write in %s (or it doesn't exist)", tmpdir);
    }
    let mut tracefile = String::new();
    let signal_path;
    let tracefp = if tmpdir.is_empty() {
        let fp = etmpfile();
        if verbose != 0 { warn!("using tmpfile() call"); }
        fp
    } else {
        tracefile = temporary_filename(&tmpdir);
        signal_path = CString::new(tracefile.as_str()).unwrap();
        unsafe {
            TRACEFILE = signal_path.as_ptr();
            for sig in [2, 3, 1, 15] { signal(sig, closefiles); }
        }
        let fp = efopen(&tracefile, "w+");
        if verbose != 0 { warn!("putting temporary files in %s", tmpdir); }
        fp
    };
    checkpars();
    let mut nx2 = 0i32;
    loop {
        nx2 += 1;
        efwrite(&intrace.data[..nx1 as usize], &tracefp);
        if gettr(&mut intrace) == 0 { break; }
    }
    let nx1fft = npfar(nx1);
    let nx2fft = npfa(nx2);
    if nx1fft >= SU_NFLTS as i32 || nx1fft >= 720720 { err!("Padded nx1=%d--too big", nx1fft); }
    if nx2fft >= SU_NFLTS as i32 || nx2fft >= 720720 { err!("Padded nx2=%d--too big", nx2fft); }
    let d1 = (2.0f64 / ((nx1fft as f32 * dx1) as f64)) as f32;
    let d2 = (2.0f64 / ((nx2fft as f32 * dx2) as f64)) as f32;
    let f1 = (-1.0f64 / dx1 as f64 - d1 as f64 / 2.0) as f32;
    let f2 = (-1.0f64 / dx2 as f64 + d2 as f64) as f32;
    let nk1 = nx1fft / 2 + 1;
    let nk2 = nx2fft;
    let rt_len = (nx1fft as usize).checked_mul(nx2fft as usize)
        .unwrap_or_else(|| err!("FFT array size overflow"));
    let ct_len = (nk1 as usize).checked_mul(nk2 as usize)
        .unwrap_or_else(|| err!("FFT array size overflow"));
    let mut rt = Vec::new();
    rt.try_reserve_exact(rt_len).unwrap_or_else(|_| err!("cannot allocate FFT array"));
    rt.resize(rt_len, 0.0f32);
    let mut ct = Vec::new();
    ct.try_reserve_exact(ct_len).unwrap_or_else(|_| err!("cannot allocate FFT array"));
    ct.resize(ct_len, Complex::default());
    erewind(&tracefp);
    for ix2 in 0..nx2 as usize {
        let row = &mut rt[ix2 * nx1fft as usize..ix2 * nx1fft as usize + nx1 as usize];
        efread(row, &tracefp);
        if ix2 % 2 != 0 { for x in row { *x = -*x; } }
    }
    efclose(&tracefp);
    pfa2rc(-1, 1, nx1fft, nx2, &rt, &mut ct);
    pfa2cc(-1, 2, nk1, nk2, &mut ct);
    for ik2 in 0..nk2 as usize {
        for ik1 in 0..nk1 as usize {
            outtrace.data[ik1] = rcabs(ct[(nk2 as usize - 1 - ik2) * nk1 as usize + nk1 as usize - 1 - ik1]);
        }
        for ik1 in nk1 as usize..nx1fft as usize {
            let row = if ik2 < nk2 as usize - 1 { ik2 + 1 } else { 0 };
            outtrace.data[ik1] = rcabs(ct[row * nk1 as usize + ik1 - nk1 as usize]);
        }
        outtrace.tracl = ik2 as i32 + 1;
        outtrace.ns = nx1fft as u16;
        outtrace.dt = 0;
        outtrace.trid = 122;
        outtrace.d1 = d1;
        outtrace.f1 = f1;
        outtrace.d2 = d2;
        outtrace.f2 = f2;
        puttr(&mut outtrace);
    }
    if !tmpdir.is_empty() { eremove(&tracefile); }
    exit(CWP_Exit());
}
