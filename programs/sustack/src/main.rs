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
//! SUSTACK - stack adjacent traces with the same header key.
//! Port of su/main/stacking/sustack.c (Revision 1.27).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: SEP: Einar Kjartansson; CWP: Jack K. Cohen, Dave Hale;
//! CENPET: Werner M. Heigl (repeat trace functionality).

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUSTACK - stack adjacent traces having the same key header word",
    " \t\t\t\t\t\t\t\t",
    "     sustack <stdin >stdout [Optional parameters]\t\t",
    " \t\t\t\t\t\t\t        ",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t        ",
    " Optional parameters: \t\t\t\t\t\t",
    " \tkey=cdp\t\theader key word to stack on\t\t",
    " \tnormpow=1.0\teach sample is divided by the\t\t",
    "\t\t\tnormpow'th number of non-zero values\t",
    "\t\t\tstacked (normpow=0 selects no division)\t",
    "\trepeat=0\t=1 repeats the stack trace nrepeat times",
    "\tnrepeat=10\trepeats stack trace nrepeat times in\t",
    "\t          \toutput file\t\t\t\t",
    " \tverbose=0\tverbose = 1 echos information\t\t",
    " \t\t\t\t\t\t\t        ",
    " Notes:\t\t\t\t\t\t\t",
    " ------\t\t\t\t\t\t\t",
    " The offset field is set to zero on the output traces, unless\t",
    " the user is stacking with key=offset. In that case, the value ",
    " of the offset field is left unchanged. \t\t        ",
    " \t\t\t\t\t\t\t        ",
    " Sushw can be used afterwards if this is not acceptable.\t",
    "\t\t\t\t\t\t\t\t",
    " For VSP users:\t\t\t\t\t\t",
    " The stack trace appears ten times in the output file when\t",
    " setting repeat=1 and nrepeat=10. Corridor stacking can be\t",
    " achieved by properly muting the upgoing data with SUMUTE\t",
    " before stacking.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    let mut normpow = 1.0f32;
    let mut key = String::new();
    let mut repeat = 0;
    let mut nrepeat = 10;
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    if getparfloat("normpow", &mut normpow) == 0 { normpow = 1.0; }
    if getparstring("key", &mut key) == 0 { key = "cdp".to_string(); }
    if getparint("repeat", &mut repeat) == 0 { repeat = 0; }
    if getparint("nrepeat", &mut nrepeat) == 0 { nrepeat = 10; }
    checkpars();
    let is_offset = key == "offset";
    let type_ = hdtype(&key);
    let indx = getindex(&key);
    let mut nsegy = gettr(&mut intrace);
    let nt = intrace.ns as usize;
    let mut outtrace = intrace.clone();
    let mut nnz = vec![0i32; nt];
    for i in 0..nt { nnz[i] = if intrace.data[i] != 0.0 { 1 } else { 0 }; }
    let mut fold = if intrace.nhs != 0 { intrace.nhs as i32 } else { 1 };
    let mut newtracl = 1i32;
    let mut val = Value::default();
    let mut valnew = Value::default();
    gethval(&intrace, indx, &mut val);
    while nsegy != 0 {
        nsegy = gettr(&mut intrace);
        gethval(&intrace, indx, &mut valnew);
        if valcmp(type_, val, valnew) != 0 || nsegy == 0 {
            if verbose != 0 {
                fprintf!(stderr(), "val=");
                fprintfval(&stderr(), type_, val);
                fprintf!(stderr(), "\tfold=%d\n", fold);
            }
            outtrace.nhs = fold as i16;
            outtrace.tracl = newtracl;
            newtracl = newtracl.wrapping_add(1);
            if !is_offset { outtrace.offset = 0; }
            if normpow != 0.0 && fold != 1 {
                for i in 0..nt {
                    let nnzi = nnz[i] as f32;
                    if nnzi != 0.0 {
                        outtrace.data[i] = ((outtrace.data[i] as f64) / (nnzi as f64).powf(normpow as f64)) as f32;
                    }
                }
            }
            if repeat != 0 {
                for _ in 0..nrepeat { puttr(&mut outtrace); }
            } else { puttr(&mut outtrace); }
            if nsegy != 0 { outtrace = intrace.clone(); }
            for i in 0..nt { nnz[i] = if intrace.data[i] != 0.0 { 1 } else { 0 }; }
            fold = if intrace.nhs != 0 { intrace.nhs as i32 } else { 1 };
            val = valnew;
        } else {
            for i in 0..nt {
                let datum = intrace.data[i];
                if datum != 0.0 { nnz[i] = nnz[i].wrapping_add(1); }
                outtrace.data[i] += datum;
            }
            fold = fold.wrapping_add(1);
        }
    }
    exit(CWP_Exit());
}
