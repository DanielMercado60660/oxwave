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


//! SUDUMPTRACE - print selected header values and data.
//! Port of su/main/headers/sudumptrace.c (Revision 1.3).
//! Credits: MTU: David Forel, Jan 2005.

use oxwave_core::prelude::*;
use oxwave_tempfile::{temporary_filename, temporary_stream};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                   ",
    " SUDUMPTRACE - print selected header values and data.              ",
    "               Print first num traces.                             ",
    "               Use SUWIND to skip traces.                          ",
    "                                                                   ",
    " sudumptrace < stdin [> ascii_file]                                ",
    "                                                                   ",
    " Optional parameters:                                              ",
    "     num=4                    number of traces to dump             ",
    "     key=key1,key2,...        key(s) to print above trace values   ",
    "     hpf=0                    header print format is float         ",
    "                              =1 print format is exponential       ",
    "                                                                   ",
    " Examples:                                                         ",
    "   sudumptrace < inseis.su            PRINTS: 4 traces, no headers ",
    "   sudumptrace < inseis.su key=tracf,offset                        ",
    "   sudumptrace < inseis.su num=7 key=tracf,offset > info.txt       ",
    "   sudumptrace < inseis.su num=7 key=tracf,offset hpf=1 > info.txt ",
    "                                                                   ",
    " Related programs: suascii, sugethw                                ",
    "                                                                   ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let writable = match std::fs::metadata(&tmpdir) {
            Ok(metadata) if metadata.is_dir() => temporary_stream(Some(&tmpdir)).is_some(),
            Ok(metadata) => !metadata.permissions().readonly(),
            Err(_) => false,
        };
        if !writable {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let mut dt = (tr.dt as f64 / 1000000.0) as f32;
    if dt == 0.0 { getparfloat("dt", &mut dt); }
    if dt == 0.0 { mustgetparfloat("dt", &mut dt); }
    let delrt = (tr.delrt as f64 / 1000.0) as f32;
    let mut num = 0;
    let mut numtr = 4;
    if getparint("num", &mut num) != 0 { numtr = num; }
    let nkeys = countparval("key");
    let mut key = vec![String::new(); SU_NKEYS];
    if nkeys != 0 { getparstringarray("key", &mut key); }
    let mut hedr = vec![0.0f32; (nkeys as usize).saturating_mul(numtr.max(1) as usize)];
    let mut hpf = 0;
    if getparint("hpf", &mut hpf) == 0 { hpf = 0; }
    checkpars();
    let (tracefp, headerfp, tracefile, headerfile) = if tmpdir.is_empty() {
        (etmpfile(), etmpfile(), String::new(), String::new())
    } else {
        let tracefile = temporary_filename(&tmpdir);
        let mut headerfile = temporary_filename(&tmpdir);
        // The temporary-name library can return its fixed fallback twice if
        // /tmp is unavailable; keep the two storage files distinct.
        if headerfile == tracefile { headerfile.push_str(".header"); }
        let tracefp = efopen(&tracefile, "w+");
        let headerfp = efopen(&headerfile, "w+");
        (tracefp, headerfp, tracefile, headerfile)
    };
    let mut ntr = 0;
    loop {
        ntr += 1;
        efwrite(&[0u8; HDRBYTES], &headerfp);
        efwrite(&tr.data[..nt as usize], &tracefp);
        for ikey in 0..nkeys as usize {
            let mut val = Value::default();
            gethdval(&tr, &key[ikey], &mut val);
            hedr[(ntr - 1) as usize * nkeys as usize + ikey] = vtof(hdtype(&key[ikey]), val);
        }
        if ntr >= numtr || gettr(&mut tr) == 0 { break; }
    }
    erewind(&tracefp);
    erewind(&headerfp);
    let mut data = vec![0.0f32; (nt * ntr) as usize];
    efread(&mut data, &tracefp);
    efclose(&tracefp);
    if !tmpdir.is_empty() { eremove(&tracefile); }
    dump(&data, dt, &hedr, &key, delrt, nkeys, ntr, nt, hpf);
    efclose(&headerfp);
    if !tmpdir.is_empty() { eremove(&headerfile); }
    exit(CWP_Exit());
}

fn dump(data: &[f32], dt: f32, hedr: &[f32], key: &[String], delrt: f32,
        nkeys: i32, ntr: i32, ntime: i32, hpf: i32) {
    printf!("\nnum traces = %d    num samples = %d \n\n", ntr, ntime);
    for k in 0..nkeys as usize {
        printf!("%17s  ", key[k]);
        for m in 0..ntr as usize {
            if hpf == 0 { printf!("%11.4f\t", hedr[m * nkeys as usize + k]); }
            else { printf!("%11.4e\t", hedr[m * nkeys as usize + k]); }
        }
        putchar(b'\n' as i32);
    }
    putchar(b'\n' as i32);
    printf!("\nCounter    Time     Values\n");
    for i in 1..=ntime {
        printf!(" %6d ", i);
        printf!(" %8.3f   ", dt * i as f32 + delrt);
        for j in 1..=ntr {
            printf!("%11.4e\t", data[((j - 1) * ntime + i - 1) as usize]);
        }
        putchar(b'\n' as i32);
    }
    putchar(b'\n' as i32);
}

