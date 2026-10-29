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
//! SUSTRIP - remove SEGY headers from traces.
//! Port of su/main/headers/sustrip.c (Revision 1.20).
//! Credits: SEP: Einar Kjartansson; CWP: Jack K. Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUSTRIP - remove the SEGY headers from the traces\t\t",
    " \t\t\t\t\t\t\t\t",
    " sustrip <stdin >stdout head=/dev/null outpar=/dev/tty ftn=0\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \thead=/dev/null\t\tfile to save headers in\t\t",
    " \t\t\t\t\t\t\t\t",
    " \toutpar=/dev/tty\t\toutput parameter file, contains:",
    " \t\t\t\tnumber of samples (n1=)\t\t",
    " \t\t\t\tnumber of traces (n2=)\t\t",
    " \t\t\t\tsample rate in seconds (d1=)\t",
    " \t\t\t\t\t\t\t\t",
    " \tftn=0\t\t\tFortran flag\t\t\t",
    " \t\t\t\t0 = write unformatted for C\t",
    " \t\t\t\t1 = ... for Fortran\t\t",
    " \t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " Invoking head=filename will write the trace headers into filename.",
    " You may paste the headers back onto the traces with supaste\t",
    " See:  sudoc  supaste \t for more information \t\t\t",
    " Related programs: supaste, suaddhead\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut ntr = 0i32;
    initargs(&argv);
    requestdoc(1, SDOC);
    match filestat(STDOUT) {
        FileType::BadFileType => { warn!("stdout is illegal filetype"); pagedoc(SDOC); }
        FileType::Tty => { warn!("stdout can't be tty"); pagedoc(SDOC); }
        _ => {}
    }
    let mut head = String::new();
    let mut outpar = String::new();
    let mut ftn = 0i32;
    if getparstring("head", &mut head) == 0 { head = "/dev/null".to_string(); }
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    if getparint("ftn", &mut ftn) == 0 { ftn = 0; }
    if ftn != 0 && ftn != 1 { err!("ftn=%d must be 0 or 1", ftn); }
    checkpars();
    let headfp = efopen(&head, "w");
    let outparfp = efopen(&outpar, "w");
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let ns = tr.ns as usize;
    let nsbytes = (ns * 4) as u32;
    loop {
        if ftn != 0 { efwrite(&[nsbytes], &stdout()); }
        match tr.trid {
            CHARPACK => { efwrite(&tr.data_bytes()[..ns], &stdout()); }
            SHORTPACK => { efwrite(&tr.data_bytes()[..ns * 2], &stdout()); }
            _ => { efwrite(&tr.data[..ns], &stdout()); }
        }
        if ftn != 0 { efwrite(&[nsbytes], &stdout()); }
        efwrite(tr.hdr.as_bytes(), &headfp);
        ntr = ntr.wrapping_add(1);
        if gettr(&mut tr) == 0 { break; }
    }
    fprintf!(outparfp, "n1=%d n2=%d d1=%f\nnt=%d ntr=%d dt=%f\nns=%d\n",
        tr.ns as i32, ntr, tr.dt as f64 / 1000000.0,
        tr.ns as i32, ntr, tr.dt as f64 / 1000000.0, tr.ns as i32);
    exit(CWP_Exit());
}
