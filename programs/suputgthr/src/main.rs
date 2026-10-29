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

//! SUPUTGTHR - split traces into gathers.
//! Port of su/main/windowing_sorting_muting/suputgthr.c (Revision 1.5).
//! Credits: Balazs Nemeth, Potash Corporation, Saskatoon Saskatchewan; given to CWP in 2008.

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUPUTGTHR - split the stdout flow to gathers on the bases of given\t",
    " \t\tkey parameter. \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "\tsuputgthr <stdin   dir= [Optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " dir=\t\tName of directory where to put the gathers\t\t",
    " Optional parameters: \t\t\t\t\t\t\t",
    " key=ep\t\theader key word to watch   \t\t\t",
    " suffix=\".hsu\"\textension of the output files\t\t\t",
    " verbose=0\t\tverbose = 1 echos information\t\t\t",
    " numlength=7\t\tLength of numeric part of filename\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t    \t\t\t\t\t",
    " The name of the file is constructed from the key parameter. Traces\t",
    " are put into a temporary disk file, and renamed when key parameter\t",
    " changes in the input flow to \"key.suffix\". The result is that the\t",
    " directory \"dir\" contains separate files by \"key\" ensemble. \t",
    " \t\t\t\t\t\t\t\t\t",
    " Header field modified:  ntr  to be the number of traces in a given \t",
    " ensemble.\t\t\t\t\t\t\t\t",
    " Related programs: sugetgthr, susplit \t\t\t\t\t",
];

fn open_file(name: &str) -> CFile {
    match fopen(name, "w+") {
        Some(fp) => fp,
        None => err!("Cannot open file %s\n", name),
    }
}

fn output_name(directory: &str, key: &str, val_i: i32, suffix: &str, numlength: i32) -> String {
    let format = sprintf!("%%s%%s%%0%dd%%s", numlength);
    sprintf!(&format, format!("{}/", directory), format!("{}=", key), val_i, suffix)
}

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut tmptr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut dir = String::new();
    let mut key = String::new();
    let mut suffix = String::new();
    let mut verbose = 0i32;
    let mut numlength = 7i32;
    mustgetparstring("dir", &mut dir);
    getparint("verbose", &mut verbose);
    if getparstring("key", &mut key) == 0 { key = "ep".to_string(); }
    if getparstring("suffix", &mut suffix) == 0 { suffix = ".hsu".to_string(); }
    getparint("numlength", &mut numlength);
    checkpars();

    let type_ = hdtype(&key);
    let index = getindex(&key);
    let nsegy = gettr(&mut intrace);
    let mut tempfilename = temporary_filename(&dir);
    let mut tmpfp = open_file(&tempfilename);
    if verbose != 0 { warn!(" temporary filename = %s", tempfilename); }
    let mut val = Value::default();
    gethval(&intrace, index, &mut val);
    let mut ntr = 0i32;
    loop {
        ntr = ntr.wrapping_add(1);
        let mut valnew = Value::default();
        gethval(&intrace, index, &mut valnew);
        if valcmp(type_, val, valnew) != 0 || nsegy == 0 {
            let val_i = vtoi(type_, val);
            let outfile = output_name(&dir, &key, val_i, &suffix, numlength);
            if verbose != 0 { warn!(" outfile = %s", outfile); }
            erewind(&tmpfp);
            let outfp = open_file(&outfile);
            while fgettr(&tmpfp, &mut tmptr) != 0 {
                tmptr.ntr = ntr;
                fputtr(&outfp, &mut tmptr);
            }
            fclose(&tmpfp);
            fclose(&outfp);
            let _ = std::fs::remove_file(&tempfilename);
            if verbose != 0 { warn!("val= %i", val_i); }
            tempfilename = temporary_filename(&dir);
            tmpfp = open_file(&tempfilename);
            val = valnew;
            ntr = 0;
        }
        fputtr(&tmpfp, &mut intrace);
        if gettr(&mut intrace) == 0 { break; }
    }
    erewind(&tmpfp);
    let val_i = vtoi(type_, val);
    let outfile = output_name(&dir, &key, val_i, &suffix, numlength);
    let outfp = open_file(&outfile);
    while fgettr(&tmpfp, &mut tmptr) != 0 {
        tmptr.ntr = ntr;
        fputtr(&outfp, &mut tmptr);
    }
    fclose(&tmpfp);
    fclose(&outfp);
    let _ = std::fs::remove_file(&tempfilename);
    if verbose != 0 { warn!(" outfile = %s", outfile); }
    if verbose != 0 { warn!("val= %i", val_i); }
    exit(CWP_Exit());
}
