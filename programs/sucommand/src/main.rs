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
// Copyright (c) Virginia Tech. All rights reserved.

//! SUCOMMAND - pipe traces having the same key header word to command.
//! Port of su/main/windowing_sorting_muting/sucommand.c (Revision 1.6).
//! Credits: VT: Matthias Imhof.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUCOMMAND - pipe traces having the same key header word to command\t",
    " \t\t\t\t\t\t\t\t\t",
    "     sucommand <stdin >stdout [Optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters: \t\t\t\t\t\t\t",
    " \tverbose=0\t\twordy output\t\t\t\t",
    " \tkey=cdp\t\t\theader key word to pipe on\t\t",
    " \tcommand=\"suop nop\"    command piped into\t\t\t",
    "\tdir=0\t\t0:  change of header key\t\t\t",
    "\t\t\t-1: break in monotonic decrease of header key\t",
    "\t\t\t+1: break in monotonic increase of header key\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "Notes:\t\t\t\t\t\t\t\t\t",
    " This program permits limited parallel processing capability by opening",
    " pipes for processes, signaled by a change in key header word value.\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut trace = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    let mut command = String::new();
    let mut key = String::new();
    let mut dir = 0i32;
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    if getparstring("command", &mut command) == 0 { command = "suop nop".into(); }
    if getparstring("key", &mut key) == 0 { key = "cdp".into(); }
    if getparint("dir", &mut dir) == 0 { dir = 0; }
    checkpars();
    dir = dir.wrapping_neg();
    let type_ = hdtype(&key);
    let indx = getindex(&key);
    let mut pipe = epopen(&command, "w");
    let mut fold = 1i32;
    let mut val = Value::default();
    let mut valnew = Value::default();
    gettr(&mut trace);
    gethval(&trace, indx, &mut val);
    fputtr(&pipe, &mut trace);
    loop {
        if gettr(&mut trace) == 0 { break; }
        gethval(&trace, indx, &mut valnew);
        let cmp = valcmp(type_, val, valnew);
        val = valnew;
        if (dir == 0 && cmp != 0) || (dir != 0 && dir != cmp) {
            if verbose != 0 {
                fprintf!(stderr(), "val=");
                fprintfval(&stderr(), type_, val);
                fprintf!(stderr(), "\tgather size=%d\n", fold);
            }
            epclose(&pipe);
            pipe = epopen(&command, "w");
            fold = 1;
            fputtr(&pipe, &mut trace);
        } else {
            fputtr(&pipe, &mut trace);
            fold = fold.wrapping_add(1);
        }
    }
    if verbose != 0 {
        fprintf!(stderr(), "val=");
        fprintfval(&stderr(), type_, val);
        fprintf!(stderr(), "\tgather size=%d\n", fold);
    }
    epclose(&pipe);
    exit(CWP_Exit());
}
