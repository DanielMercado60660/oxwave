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
//! SUFLIP - flip a data set in various ways.
//! Port of su/main/operations/suflip.c (Revision 1.31).
//! Credits: CWP: Chris Liner, Jack K. Cohen, John Stockwell.

use oxwave_core::prelude::*;
use oxwave_tempfile::{temporary_filename, temporary_stream};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUFLIP - flip a data set in various ways\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " suflip <data1 >data2 flip=1 verbose=0\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tflip=1 \trotational sense of flip\t\t\t",
    " \t\t\t+1  = flip 90 deg clockwise\t\t",
    " \t\t\t-1  = flip 90 deg counter-clockwise\t",
    " \t\t\t 0  = transpose data\t\t\t",
    " \t\t\t 2  = flip right-to-left\t\t",
    " \t\t\t 3  = flip top-to-bottom\t\t",
    " \ttmpdir=\t if non-empty, use the value as a directory path",
    "\t\t prefix for storing temporary files; else if\t",
    "\t         the CWP_TMPDIR environment variable is set use\t",
    "\t         its value for the path; else use tmpfile()\t",
    " \t\t\t\t\t\t\t\t",
    " \tverbose=0\tverbose = 1 echoes flip info\t\t",
    " \t\t\t\t\t\t\t\t",
    " NOTE:  tr.dt header field is lost if flip=-1,+1.  It can be\t",
    "        reset using sushw.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " EXAMPLE PROCESSING SEQUENCES:\t\t\t\t\t",
    "   1.\tsuflip flip=-1 <data1 | sushw key=dt a=4000 >data2\t",
    " \t\t\t\t\t\t\t\t",
    "   2.\tsuflip flip=2 <data1 | suflip flip=2 >data1_again\t",
    " \t\t\t\t\t\t\t\t",
    "   3.\tsuflip tmpdir=/scratch <data1 | ...\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Caveat:  may fail on large files.\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    let mut flip = 1;
    getparint("verbose", &mut verbose);
    getparint("flip", &mut flip);
    if ![-1, 0, 1, 2, 3].contains(&flip) {
        err!("flip = %d, flag must be -1, 0, 1, 2, or 3", flip);
    }
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        if std::fs::metadata(&tmpdir).map_or(true, |m| m.permissions().readonly()) {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let trdt = tr.dt;
    let (tracefp, headerfp, tracefile, headerfile) = if tmpdir.is_empty() {
        let t = etmpfile();
        let h = etmpfile();
        if verbose != 0 { warn!("using tmpfile() call"); }
        (t, h, String::new(), String::new())
    } else {
        let tname = temporary_filename(&tmpdir);
        let hname = temporary_filename(&tmpdir);
        // The library's fixed /tmp template can fail in a sandbox and return
        // the same unexpanded name twice. Use its unlinked streams then.
        let (t, h, tname, hname) = if tname == hname {
            let t = temporary_stream(Some(&tmpdir)).unwrap_or_else(|| err!("cannot create temporary file"));
            let h = temporary_stream(Some(&tmpdir)).unwrap_or_else(|| err!("cannot create temporary file"));
            (CFile::from_file(t), CFile::from_file(h), String::new(), String::new())
        } else {
            (efopen(&tname, "w+"), efopen(&hname, "w+"), tname, hname)
        };
        if verbose != 0 { warn!("putting temporary files in %s", tmpdir); }
        (t, h, tname, hname)
    };
    let mut ntr = 0usize;
    loop {
        ntr += 1;
        let hdr = tr.hdr.as_bytes();
        efwrite(hdr, &headerfp);
        efwrite(&tr.data[..nt], &tracefp);
        if gettr(&mut tr) == 0 { break; }
    }
    erewind(&tracefp);
    erewind(&headerfp);
    let mut data = vec![0.0f32; nt * ntr];
    if flip == -1 || flip == 0 || flip == 1 { CHECK_NT("ntr", ntr as i32); }
    let mut flipdata = vec![0.0f32; nt * ntr];
    efread(&mut data, &tracefp);
    efclose(&tracefp);
    if !tracefile.is_empty() { eremove(&tracefile); }
    if verbose != 0 {
        fprintf!(stderr(), "\nSUFLIP:          flip = %d\n", flip);
        fprintf!(stderr(), "input:   samples/trace = %d", nt);
        fprintf!(stderr(), "    traces = %d \n", ntr);
        if flip == -1 || flip == 0 || flip == 1 {
            fprintf!(stderr(), "output:  samples/trace = %d", ntr);
            fprintf!(stderr(), "    traces = %d \n\n", nt);
        } else {
            fprintf!(stderr(), "output:  samples/trace = %d", nt);
            fprintf!(stderr(), "    traces = %d \n\n", ntr);
        }
    }
    for irow in 0..ntr {
        for icol in 0..nt {
            let dest = match flip {
                -1 => icol * ntr + ntr - 1 - irow,
                0 => icol * ntr + irow,
                1 => (nt - 1 - icol) * ntr + irow,
                2 => (ntr - 1 - irow) * nt + icol,
                3 => irow * nt + nt - 1 - icol,
                _ => unreachable!(),
            };
            flipdata[dest] = data[irow * nt + icol];
        }
    }
    let rotated = flip == -1 || flip == 0 || flip == 1;
    let (outtr, outnt) = if rotated { (nt, ntr) } else { (ntr, nt) };
    for i in 0..outtr {
        tr.data[..outnt].copy_from_slice(&flipdata[i * outnt..(i + 1) * outnt]);
        let hdr = tr.hdr.as_bytes_mut();
        efread(hdr, &headerfp);
        tr.ns = outnt as u16;
        if !rotated { tr.dt = trdt; }
        tr.tracl = (i + 1) as i32;
        puttr(&mut tr);
    }
    efclose(&headerfp);
    if !headerfile.is_empty() { eremove(&headerfile); }
    exit(CWP_Exit());
}
