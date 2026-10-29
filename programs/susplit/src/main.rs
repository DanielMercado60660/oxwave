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

//! SUSPLIT - Split traces into different output files by keyword value.
//! Port of su/main/windowing_sorting_muting/susplit.c.
//! Credits: Geocon: Garry Perratt hacked together from various other codes.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSPLIT - Split traces into different output files by keyword value\t",
    "\t\t\t\t\t\t\t\t\t",
    "     susplit <stdin >stdout [options]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    "\tkey=cdp\t\tKey header word to split on (see segy.h)\t",
    "\tstem=split_\tStem name for output files\t\t\t",
    "\tmiddle=key\tmiddle of name of output files\t\t\t",
    "\tsuffix=.su\tSuffix for output files\t\t\t\t",
    "\tnumlength=7\tLength of numeric part of filename\t\t",
    "\tverbose=0\t=1 to echo filenames, etc.\t\t\t",
    "\tclose=1\t\t=1 to close files before opening new ones\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The most efficient way to use this program is to presort the input data",
    " into common keyword gathers, prior to using susplit.\t\t\t\t\t\t\t\t\t\t\t\t",
    " Use \"suputgthr\" to put SU data into SU data directory format.\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key = String::new();
    let mut stem = String::new();
    let mut middle = String::new();
    let mut suffix = String::new();
    let mut numlength = 7;
    let mut verbose = 0;
    let mut close = 1;
    if getparstring("key", &mut key) == 0 { key = "cdp".into(); }
    if getparstring("stem", &mut stem) == 0 { stem = "split_".into(); }
    if getparstring("middle", &mut middle) == 0 { middle = key.clone(); }
    if getparstring("suffix", &mut suffix) == 0 { suffix = ".su".into(); }
    getparint("numlength", &mut numlength);
    getparint("verbose", &mut verbose);
    getparint("close", &mut close);
    checkpars();

    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let type_ = hdtype(&key);
    let index = getindex(&key);
    if verbose != 0 { warn!("key is %s", key); }
    let mut lastval_i = 0;
    let mut outfp = None;
    // Retain unclosed streams, as fopen does when close=0.
    let mut open_files = Vec::new();
    loop {
        let mut val = Value::default();
        gethval(&tr, index, &mut val);
        let val_i = vtoi(type_, val);
        if val_i != lastval_i || outfp.is_none() {
            if let Some(fp) = outfp.take() {
                if close != 0 { efclose(&fp); } else { open_files.push(fp); }
            }
            let format = sprintf!("%%s%%s%%0%dd%%s", numlength);
            let filename = sprintf!(&format, stem, middle, val_i, suffix);
            if verbose != 0 { warn!("output file is %s", filename); }
            outfp = Some(efopen(&filename, "ab"));
        }
        fputtr(outfp.as_ref().unwrap(), &mut tr);
        lastval_i = val_i;
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
