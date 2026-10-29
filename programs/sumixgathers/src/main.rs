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

//! SUMIXGATHERS - mix two gathers.
//! Port of su/main/windowing_sorting_muting/sumixgathers.c (June 2000).
//! Credits: Daniel Trad, UBC.
//! Copyright (c) University of British Columbia, 1999. All rights reserved.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIXGATHERS - mix two gathers\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sumixgathers file1 file2 > stdout [optional parameters]\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " ntr=tr.ntr\tif ntr header field is not set, then ntr is mandatory\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters: none\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes: Both files have to be sorted by offset\t\t\t\t",
    " Mixes two gathers keeping only the traces of the first file\t\t",
    " if the offset is the same. The purpose is to substitute only\t\t",
    " traces non existing in file1 by traces interpolated store in file2. \t",
    "\t\t\t\t\t\t\t\t\t",
    " Example. If file1 is original data file and file 2 is obtained by\t",
    " resampling with Radon transform, then the output contains original  \t",
    " traces with gaps filled\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tr2 = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let fp1 = efopen(argv.get(1).map(String::as_str).unwrap_or(""), "r");
    let fp2 = efopen(argv.get(2).map(String::as_str).unwrap_or(""), "r");
    let tracefp = etmpfile();
    let headerfp = etmpfile();
    let mut scaling = 0;
    getparint("scaling", &mut scaling);
    if fgettr(&fp1, &mut tr) == 0 { err!("can't read first trace"); }
    let nt = tr.ns as usize;
    let mut ntr = 0;
    if getparint("ntr", &mut ntr) == 0 { ntr = tr.ntr; }
    checkpars();
    if ntr == 0 { err!("ntr neither set in header nor getparred!"); }
    let mut h = Vec::<f32>::new();
    let mut headers = Vec::new();
    loop {
        headers.push(tr.hdr.clone());
        efwrite(&tr.data[..nt], &tracefp);
        h.push(tr.offset as f32);
        if fgettr(&fp1, &mut tr) == 0 { break; }
    }
    erewind(&tracefp);
    erewind(&headerfp);
    let nh = h.len();
    warn!("nh=%d", nh as i32);
    while fgettr(&fp2, &mut tr2) != 0 {
        let mut flag = 1;
        for ih in 0..nh {
            if h[ih] >= 0.0 {
                if (tr2.offset as f64) > 0.999 * h[ih] as f64 && (tr2.offset as f64) < 1.001 * h[ih] as f64 { flag = 0; }
            }
            if h[ih] < 0.0 {
                if (tr2.offset as f64) < 0.99 * h[ih] as f64 && (tr2.offset as f64) > 1.01 * h[ih] as f64 { flag = 0; }
            }
        }
        if flag == 1 {
            if scaling != 0 && tr2.offset.wrapping_abs() > 0 {
                let scale = (1.0 + 0.03 * (tr2.offset.wrapping_abs() as f64 / 1000.0)) as f32;
                fprintf!(stderr(), "tr2.offset=%d,scale=%f\n", tr2.offset, scale);
                for it in 0..nt { tr2.data[it] *= scale; }
            }
            puttr(&mut tr2);
        }
    }
    for header in headers {
        tr.hdr = header;
        efread(&mut tr.data[..nt], &tracefp);
        puttr(&mut tr);
    }
    efclose(&fp1);
    efclose(&fp2);
    efclose(&tracefp);
    efclose(&headerfp);
    exit(CWP_Exit());
}
