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
//! SUPASTE - paste existing SU headers on existing binary data.
//! Port of su/main/headers/supaste.c (Revision 1.17).
//! Credits: CWP: Jack K. Cohen, November 1990.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUPASTE - paste existing SU headers on existing binary data\t",
    " \t\t\t\t\t\t\t\t",
    " supaste <bare_data >segys  ns= head=headers ftn=0\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required parameter:\t\t\t\t\t\t",
    "\tns=the number of samples per trace\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \thead=headers\tfile with segy headers\t\t\t",
    "\tftn=0\t\tFortran flag\t\t\t\t",
    "\t\t\t0 = unformatted data from C\t\t",
    "\t\t\t1 = ... from Fortran\t\t\t",
    "\tverbose=0\t1= echo number of traces pasted\t\t",
    " Caution:\t\t\t\t\t\t\t",
    "\tAn incorrect ns field will munge subsequent processing.\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " This program is used when the option head=headers is used in\t",
    " sustrip. See:   sudoc sustrip    for more details. \t\t",
    "\t\t\t\t\t\t\t\t",
    " Related programs:  sustrip, suaddhead\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut ns = 0i32;
    let mut head = String::new();
    let mut ftn = 0i32;
    let mut verbose = 0i32;
    let mut junk = [0u8; 4];
    initargs(&argv);
    requestdoc(1, SDOC);
    if getparint("ns", &mut ns) == 0 { err!("must specify ns="); }
    if getparstring("head", &mut head) == 0 { head = "headers".to_string(); }
    if getparint("ftn", &mut ftn) == 0 { ftn = 0; }
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    if ftn != 0 && ftn != 1 { err!("ftn=%d must be 0 or 1", ftn); }
    checkpars();
    let headfp = efopen(&head, "r");
    let input = stdin();
    let mut ntr = 0i32;
    while !(feof(&headfp) || ferror(&headfp) || feof(&input) || ferror(&input)) {
        // C reads one HDRBYTES-sized item, not HDRBYTES byte items.
        if efread(std::slice::from_mut(&mut tr.hdr), &headfp) == 0 {
            if verbose != 0 { warn!("ntr=%d", ntr); }
            break;
        }
        tr.ns = ns as u16;
        if ftn != 0 { efread(std::slice::from_mut(&mut junk), &input); }
        let nread = if (ns as usize) <= tr.data.len() {
            efread(&mut tr.data[..ns as usize], &input)
        } else {
            // C converts ns to size_t for fread. Even an oversized request
            // is safe if EOF arrives within tr.data; keep those short reads
            // and the u16 conversion of tr.ns. Reject an actual overflow.
            let nread = fread(&mut tr.data, &input);
            let overflow = nread == tr.data.len() && fgetc(&input) != EOF;
            if ferror(&input) {
                // Retain the original request in efread's diagnostic.
                err!("%s: efread: fread only %d items of %d", "subcalls.c", nread, ns);
            }
            if overflow {
                err!("ns=%d exceeds trace buffer capacity %d", ns, SU_NFLTS);
            }
            nread
        };
        if nread == 0 {
            warn!("header without data for trace #%d", ntr + 1);
            exit(1);
        }
        puttr(&mut tr);
        ntr += 1;
        if ftn != 0 { efread(std::slice::from_mut(&mut junk), &input); }
    }
    exit(CWP_Exit());
}
