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

//! SUPSIMAGE - PostScript image plot. Port of su/graphics/psplot/supsimage.c (Revision 1.25).
//! Credits: Dave Hale, Zhiming Li, Jack Cohen and John Stockwell.

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;
use std::ffi::CString;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUPSIMAGE - PostScript IMAGE plot of a segy data set\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " supsimage <stdin [optional parameters] | ...\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t \t",
    "\t\t\t\t\t\t\t\t\t",
    " n2=tr.ntr or number of traces in the data set (ntr is an alias for n2)",
    "\t\t\t\t\t\t\t\t\t",
    " d1=tr.d1 or tr.dt/10^6\tsampling interval in the fast dimension\t",
    "   =.004 for seismic \t\t(if not set)\t\t\t\t",
    "   =1.0 for nonseismic\t\t(if not set)\t\t\t\t",
    " \t\t\t\t\t\t\t        \t",
    " d2=tr.d2\t\t\tsampling interval in the slow dimension\t",
    "   =1.0 \t\t\t(if not set)\t\t\t\t",
    " \t\t\t\t\t\t\t        \t",
    " f1=tr.f1 or tr.delrt/10^3 or 0.0  first sample in the fast dimension\t",
    " \t\t\t\t\t\t\t        \t",
    " f2=tr.f2 or tr.tracr or tr.tracl  first sample in the slow dimension\t",
    "   =1.0 for seismic\t\t    (if not set)\t\t\t",
    "   =d2 for nonseismic\t\t    (if not set)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " verbose=0              =1 to print some useful information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir=\t \tif non-empty, use the value as a directory path\t",
    "\t\t \tprefix for storing temporary files; else if the\t",
    "\t         \tthe CWP_TMPDIR environment variable is set use\t",
    "\t         \tits value for the path; else use tmpfile()\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note that for seismic time domain data, the \"fast dimension\" is\t",
    " time and the \"slow dimension\" is usually trace number or range.\t",
    " Also note that \"foreign\" data tapes may have something unexpected\t",
    " in the d2,f2 fields, use segyclean to clear these if you can afford\t",
    " the processing time or use d2= f2= to override the header values if\t",
    " not.\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " See the psimage selfdoc for the remaining parameters.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " On NeXT:\tsupsimage < infile [optional parameters]  | open\t",
    "\t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns, ntr, tracr, tracl, delrt, trid,     ",
    "\tdt, d1, d2, f1, f2\t\t\t\t\t\t",
];

unsafe extern "C" {
    fn access(path: *const std::ffi::c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn unlink(path: *const std::ffi::c_char) -> i32;
    fn _exit(status: i32) -> !;
}

// The original removes its named temporary file on these interrupts.
static mut TRACEFILE: [u8; 8192] = [0; 8192];
extern "C" fn closefiles(_: i32) {
    unsafe {
        unlink(std::ptr::addr_of!(TRACEFILE).cast());
        _exit(1);
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let nt = tr.ns as usize;
    let mut ntr = tr.ntr;
    let mut have_ntr = ntr != 0;
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut d1 = 0.0f32;
    if getparfloat("d1", &mut d1) == 0 {
        if tr.d1 != 0.0 { d1 = tr.d1; }
        else if tr.dt != 0 { d1 = (tr.dt as f64 / 1000000.0) as f32; }
        else if seismic { d1 = 0.004; warn!("tr.dt not set, assuming dt=0.004"); }
        else { d1 = 1.0; warn!("tr.d1 not set, assuming d1=1.0"); }
    }
    let mut d2 = 0.0f32;
    if getparfloat("d2", &mut d2) == 0 { d2 = if tr.d2 != 0.0 { tr.d2 } else { 1.0 }; }
    let mut f1 = 0.0f32;
    if getparfloat("f1", &mut f1) == 0 {
        if tr.f1 != 0.0 { f1 = tr.f1; }
        else if tr.delrt != 0 { f1 = (tr.delrt as f64 / 1000.0) as f32; }
    }
    let mut f2 = 0.0f32;
    if getparfloat("f2", &mut f2) == 0 {
        if tr.f2 != 0.0 { f2 = tr.f2; }
        else if tr.tracr != 0 { f2 = tr.tracr as f32; }
        else if tr.tracl != 0 { f2 = tr.tracl as f32; }
        else if seismic { f2 = 1.0; }
    }
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 { tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default(); }
    if !tmpdir.is_empty() && unsafe { access(CString::new(tmpdir.as_str()).unwrap().as_ptr(), 2) } != 0 {
        err!("you can't write in %s (or it doesn't exist)", tmpdir);
    }
    let bindir = match std::env::var("CWPBIN") {
        Ok(s) => s,
        Err(_) => {
            let cwproot = std::env::var("CWPROOT").unwrap_or_default();
            if cwproot.is_empty() {
                warn!("CWPROOT environment variable is not set! ");
                err!("Set CWPROOT in shell environment as per instructions in CWP/SU Installation README files");
            }
            format!("{}/bin", cwproot)
        }
    };
    let mut trbuf = vec![0.0f32; nt];
    if getparint("n2", &mut ntr) != 0 || getparint("ntr", &mut ntr) != 0 { have_ntr = true; }
    let mut tracefile = String::new();
    let mut tracefp = None;
    if !have_ntr {
        if verbose != 0 { warn!("n2 not getparred and ntr header field not set"); warn!("....  counting traces"); }
        let fp = if tmpdir.is_empty() {
            let fp = etmpfile();
            if verbose != 0 { warn!("using tmpfile() call"); }
            fp
        } else {
            tracefile = temporary_filename(&tmpdir);
            let name = tracefile.as_bytes();
            if name.len() < 8192 {
                unsafe {
                    std::ptr::copy_nonoverlapping(name.as_ptr(), std::ptr::addr_of_mut!(TRACEFILE).cast(), name.len());
                    for sig in [2, 3, 1, 15] { signal(sig, closefiles); }
                }
            }
            let fp = efopen(&tracefile, "w+");
            if verbose != 0 { warn!("putting temporary files in %s", tmpdir); }
            fp
        };
        ntr = 0;
        loop {
            ntr += 1;
            efwrite(&tr.data[..nt], &fp);
            if gettr(&mut tr) == 0 { break; }
        }
        tracefp = Some(fp);
    }
    let mut plotcmd = sprintf!("%s/psimage n1=%d n2=%d d1=%f d2=%f f1=%f f2=%f", bindir, nt as i32, ntr, d1 as f64, d2 as f64, f1 as f64, f2 as f64);
    for arg in argv.iter().skip(1) {
        if !arg.starts_with("d1=") && !arg.starts_with("d2=") && !arg.starts_with("f1=") && !arg.starts_with("f2=") {
            plotcmd.push_str(" \""); plotcmd.push_str(arg); plotcmd.push('"');
        }
    }
    let plotfp = epopen(&plotcmd, "w");
    if let Some(fp) = tracefp {
        erewind(&fp);
        for _ in 0..ntr { efread(&mut trbuf, &fp); efwrite(&trbuf, &plotfp); }
        epclose(&plotfp);
        efclose(&fp);
        if !tmpdir.is_empty() { eremove(&tracefile); }
    } else {
        loop {
            efwrite(&tr.data[..nt], &plotfp);
            if gettr(&mut tr) == 0 { break; }
        }
        epclose(&plotfp);
    }
    exit(CWP_Exit());
}
