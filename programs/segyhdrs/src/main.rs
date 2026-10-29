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


//! SEGYHDRS - make SEG-Y ascii and binary headers for segywrite.
//! Port of su/main/data_conversion/segyhdrs.c (Revision 1.31).
//! Credits: CWP: Jack K. Cohen, John Stockwell; MOBIL: Stew Levin.

use oxwave_core::prelude::*;
use oxwave_core::segy::Bhed;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SEGYHDRS - make SEG-Y ascii and binary headers for segywrite\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " segyhdrs [ < sudata ] [optional parameters] [ > copy of sudata ]      ",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tns=  if no input trace header\t\t\t\t\t",
    "\tdt=  if no input trace header\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " \tns=tr.ns from header    number of samples on input traces\t",
    " \tdt=tr.dt from header\tsample rate (microseconds) from traces\t",
    " \tbfile=binary\t\tname of file containing binary block\t",
    " \thfile=header\t\tname of file containing ascii block\t",
    "   Some binary header fields are set:\t\t\t\t\t",
    " \tjobid=1\t\t\tjob id field\t\t\t\t",
    " \tlino=1\t\t\tline number (only one line per reel)\t",
    " \treno=1\t\t\treel number\t\t\t\t",
    " \tformat=1\t\tdata format\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " All other fields are set to 0, by default.\t\t\t\t",
    " To set any binary header field, use sukeyword to find out\t\t",
    " the appropriate keyword, then use the getpar form:\t\t\t",
    " \tkeyword=value\tto set keyword to value\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The header file is created as ascii and is translated to ebcdic\t",
    " by segywrite before being written to tape.  Its contents are\t\t",
    " formal but can be edited after creation as long as the forty\t\t",
    " line format is maintained.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat: This program has not been tested under XDR for machines       ",
    "\t not having a 2 byte unsigned short integral data type.\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    let redin = !matches!(filestat(STDIN), FileType::Tty | FileType::Directory | FileType::BadFileType);
    let redout = !matches!(filestat(STDOUT), FileType::Tty | FileType::Directory | FileType::BadFileType);
    if !redin { requestdoc(2, SDOC); }
    let mut hns = 0u16;
    let mut hdt = 0u16;
    if redin {
        if gettr(&mut tr) == 0 { err!("can't get first trace"); }
        if getparushort("ns", &mut hns) == 0 { hns = tr.ns; }
        if getparushort("dt", &mut hdt) == 0 { hdt = tr.dt; }
        if hns == 0 { warn!("hns not set in binary header, consider using ns="); }
        if hdt == 0 { warn!("hdt not set in binary header, consider using dt="); }
    } else {
        if getparushort("ns", &mut hns) == 0 { err!("need ns="); }
        if getparushort("dt", &mut hdt) == 0 { err!("need dt="); }
    }
    let mut hfile = String::new();
    let mut bfile = String::new();
    if getparstring("hfile", &mut hfile) == 0 { hfile = "header".into(); }
    if getparstring("bfile", &mut bfile) == 0 { bfile = "binary".into(); }
    let mut bh = Bhed::zeroed();
    bh.hns = hns;
    bh.hdt = hdt;
    bh.jobid = 1; bh.lino = 1; bh.reno = 1; bh.format = 1;
    getparint("jobid", &mut bh.jobid);
    getparint("lino", &mut bh.lino);
    getparint("reno", &mut bh.reno);
    getparshort("ntrpr", &mut bh.ntrpr);
    getparshort("nart", &mut bh.nart);
    getparushort("dto", &mut bh.dto);
    getparushort("nso", &mut bh.nso);
    getparshort("format", &mut bh.format);
    getparshort("fold", &mut bh.fold);
    getparshort("tsort", &mut bh.tsort);
    getparshort("vscode", &mut bh.vscode);
    getparshort("hsfs", &mut bh.hsfs);
    getparshort("hsfe", &mut bh.hsfe);
    getparshort("hslen", &mut bh.hslen);
    getparshort("hstyp", &mut bh.hstyp);
    getparshort("schn", &mut bh.schn);
    getparshort("hstas", &mut bh.hstas);
    getparshort("hstae", &mut bh.hstae);
    getparshort("htatyp", &mut bh.htatyp);
    getparshort("hcorr", &mut bh.hcorr);
    getparshort("bgrcv", &mut bh.bgrcv);
    getparshort("rcvm", &mut bh.rcvm);
    getparshort("mfeet", &mut bh.mfeet);
    getparshort("polyt", &mut bh.polyt);
    getparshort("vpol", &mut bh.vpol);
    checkpars();
    let headerfp = efopen(&hfile, "w");
    let binaryfp = efopen(&bfile, "w");
    for line in ["C      This tape was made at the", "C", "C      Center for Wave Phenomena", "C      Colorado School of Mines", "C      Golden, CO, 80401"] {
        let buf = sprintf!("%-79.79s\n", line);
        efwrite(buf.as_bytes(), &headerfp);
    }
    for _ in 0..35 {
        let buf = sprintf!("%-79.79s\n", "C");
        efwrite(buf.as_bytes(), &headerfp);
    }
    efwrite(&bh.to_xdr(), &binaryfp);
    efclose(&headerfp);
    efclose(&binaryfp);
    if redin {
        loop {
            if redout { puttr(&mut tr); }
            if gettr(&mut tr) == 0 { break; }
        }
    }
    exit(CWP_Exit());
}
