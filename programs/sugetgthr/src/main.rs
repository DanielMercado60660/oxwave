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

//! SUGETGTHR - get SU files from a directory. Port of
//! su/main/windowing_sorting_muting/sugetgthr.c (Revision 1.5).
//! Credits: none listed in the C source.

use oxwave_core::prelude::*;
use std::io::Read;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUGETGTHR - Gets su files from a directory and put them               ",
    "             throught the unix pipe. This creates continous data flow.\t",
    " \t\t\t\t\t\t\t\t\t",
    "  sugetgthr  <stdin >sdout   \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dir=            Name of directory to fetch data from \t\t\t",
    " \t          Every file in the directory is treated as an su file\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " verbose=0\t\t=1 more chatty\t\t\t\t\t",
    " vt=0\t\t\t=1 allows gathers with variable length traces\t",
    " \t\t\tno header checking is done!\t\t\t",
    " ns=\t\t\tmust be specified if vt=1; number of samples to read",
    " \t\t\t\t\t\t\t\t\t",
];

fn raw_read(fp: &mut std::fs::File, tr: &mut [u8], ns: i32) -> usize {
    let nheader = fp.read(&mut tr[..240]).unwrap_or(0);
    let input_ns = u16::from_ne_bytes([tr[114], tr[115]]) as usize;
    let end = (ns.max(input_ns as i32) as usize).min(65535);
    tr[240 + input_ns * 4..240 + end * 4].fill(0);
    let ndata = fp.read(&mut tr[240..240 + input_ns * 4]).unwrap_or(0);
    nheader + ndata
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut dir = String::new();
    let mut verbose = 0;
    let mut vt = 0;
    let mut ns = 0;
    initargs(&argv);
    requestdoc(1, SDOC);
    mustgetparstring("dir", &mut dir);
    getparint("verbose", &mut verbose);
    getparint("vt", &mut vt);
    if vt != 0 { mustgetparint("ns", &mut ns); }
    checkpars();

    let dp = match std::fs::read_dir(&dir) {
        Ok(dp) => dp,
        Err(_) => err!(" %s directory not found\n", dir),
    };
    let mut raw = vec![0u8; 240 + 65535 * 4];
    for entry in dp.flatten() {
        let fname = entry.file_name();
        let ffname = format!("{}/{}", dir, fname.to_string_lossy());
        if verbose == 1 { warn!("%s", ffname); }
        let size = std::fs::metadata(&ffname).map(|s| s.len()).unwrap_or(0);
        if size > 0 {
            if vt != 0 {
                if let Ok(mut fp) = std::fs::File::open(&ffname) {
                    // The C write can include part of the header for negative ns.
                    // Outside this range it accesses beyond the trace object.
                    if !(-60..=65535).contains(&ns) {
                        err!("ns=%d is outside the supported trace buffer", ns);
                    }
                    raw_read(&mut fp, &mut raw, ns);
                    loop {
                        raw[114..116].copy_from_slice(&(ns as u16).to_ne_bytes());
                        efwrite(&raw[..(240 + ns * 4) as usize], &stdout());
                        if raw_read(&mut fp, &mut raw, ns) == 0 { break; }
                    }
                }
            } else {
                let fp = efopen(&ffname, "r");
                fgettr(&fp, &mut tr);
                loop {
                    puttr(&mut tr);
                    if fgettr(&fp, &mut tr) == 0 { break; }
                }
                efclose(&fp);
            }
        } else {
            warn!(" File %s has zero size, skipped.\n", ffname);
        }
    }
    exit(CWP_Exit());
}
