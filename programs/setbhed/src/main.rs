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


//! SETBHED - set fields in a SEG-Y binary tape header file.
//! Port of su/main/data_conversion/setbhed.c (Revision 1.8).
//! Credits (from the C source): CWP: John Stockwell, 11 Nov 1994.

use oxwave_core::prelude::*;
use oxwave_core::segy::Bhed;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SETBHED - SET the fields in a SEGY Binary tape HEaDer file, as would be",
    " \t    produced by segyread and segyhdrs\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " setbhed par= [optional parameters]\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameter:\t\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tbfile=binary\t\toutput binary tape header file\t\t",
    "\tpar=\t\t\t=parfile\t\t\t\t",
    " Set field by field, if desired:\t\t\t\t\t",
    " \tjobid=\t\t\tjob id field\t\t\t\t",
    " \tlino=\t\t\tline number (only one line per reel)\t",
    " \treno=\t\t\treel number\t\t\t\t",
    " \tformat=\t\t\tdata format\t\t\t\t",
    " ... etc....\t\t\t\t\t\t\t\t",
    " To set any binary header field, use sukeyword to find out\t\t",
    " the appropriate keyword, then use the getpar form:\t\t\t",
    " \tkeyword=value\tto set keyword to value\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " As with all other programs in the CWP/SU package that use getpars, \t",
    " (GET PARameters from the command line) a file filled with such\t",
    " statments may be included via option par=parfile. In particular, a\t",
    " parfile created by   \"bhedtopar\"  may be used as input for the program",
    " \"setbhed\".\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The binary header file that results from running segyread may have the",
    " wrong byte order. You will need to use \"swapbhed\" to change the byte, order before applying this program. \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    "   segyread tape=yourdata.segy bfile=yourdata.b > yourdata.su\t\t",
    " If  \t\t\t\t\t\t\t\t\t",
    "   bhedtopar < yourdata.b | more \t\t\t\t\t",
    " shows impossible values, then apply \t\t\t\t\t",
    "   swapbhed < yourdata.b > swapped.b\t\t\t\t\t",
    " then apply \t\t\t\t\t\t\t\t",
    "   bhedtopar < swapped.b | more \t\t\t\t\t",
    "   bhedtopar < swapped.b outpar=parfile\t\t\t\t",
    " hand edit parfile, and then apply \t\t\t\t\t",
    "  setbhed par=parfile bfile=swapped.b > new.b\t\t\t\t",
    " then apply \t\t\t\t\t\t\t\t",
    "   segywrite tape=fixeddata.segy bfile=new.b < yourdata.su\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat: This program breaks if a \"short\" isn't 2 bytes since\t",
    "         the SEG-Y standard demands a 2 byte integer for ns.\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut bfile = String::new();
    if getparstring("bfile", &mut bfile) == 0 { bfile = "binary".to_string(); }
    let mut jobid: i32 = 0;
    getparint("jobid", &mut jobid);
    let mut lino: i32 = 0;
    getparint("lino", &mut lino);
    let mut reno: i32 = 0;
    getparint("reno", &mut reno);
    let mut ntrpr: i32 = 0;
    getparint("ntrpr", &mut ntrpr);
    let mut hdt: i32 = 0;
    getparint("hdt", &mut hdt);
    let mut nart: i32 = 0;
    getparint("nart", &mut nart);
    let mut dto: i32 = 0;
    getparint("dto", &mut dto);
    let mut hns: i32 = 0;
    getparint("hns", &mut hns);
    let mut nso: i32 = 0;
    getparint("nso", &mut nso);
    let mut format: i32 = 0;
    getparint("format", &mut format);
    let mut fold: i32 = 0;
    getparint("fold", &mut fold);
    let mut tsort: i32 = 0;
    getparint("tsort", &mut tsort);
    let mut vscode: i32 = 0;
    getparint("vscode", &mut vscode);
    let mut hsfs: i32 = 0;
    getparint("hsfs", &mut hsfs);
    let mut hsfe: i32 = 0;
    getparint("hsfe", &mut hsfe);
    let mut hslen: i32 = 0;
    getparint("hslen", &mut hslen);
    let mut hstyp: i32 = 0;
    getparint("hstyp", &mut hstyp);
    let mut schn: i32 = 0;
    getparint("schn", &mut schn);
    let mut hstas: i32 = 0;
    getparint("hstas", &mut hstas);
    let mut hstae: i32 = 0;
    getparint("hstae", &mut hstae);
    let mut htatyp: i32 = 0;
    getparint("htatyp", &mut htatyp);
    let mut hcorr: i32 = 0;
    getparint("hcorr", &mut hcorr);
    let mut bgrcv: i32 = 0;
    getparint("bgrcv", &mut bgrcv);
    let mut rcvm: i32 = 0;
    getparint("rcvm", &mut rcvm);
    let mut mfeet: i32 = 0;
    getparint("mfeet", &mut mfeet);
    let mut polyt: i32 = 0;
    getparint("polyt", &mut polyt);
    let mut vpol: i32 = 0;
    getparint("vpol", &mut vpol);
    checkpars();
    let binaryfp = efopen(&bfile, "w");
    let mut bh = Bhed::zeroed();
    bh.jobid = jobid;
    bh.lino = lino;
    bh.reno = reno;
    bh.ntrpr = ntrpr as i16;
    bh.hdt = hdt as u16;
    bh.nart = nart as i16;
    bh.dto = dto as u16;
    bh.hns = hns as u16;
    bh.nso = nso as u16;
    bh.format = format as i16;
    bh.fold = fold as i16;
    bh.tsort = tsort as i16;
    bh.vscode = vscode as i16;
    bh.hsfs = hsfs as i16;
    bh.hsfe = hsfe as i16;
    bh.hslen = hslen as i16;
    bh.hstyp = hstyp as i16;
    bh.schn = schn as i16;
    bh.hstas = hstas as i16;
    bh.hstae = hstae as i16;
    bh.htatyp = htatyp as i16;
    bh.hcorr = hcorr as i16;
    bh.bgrcv = bgrcv as i16;
    bh.rcvm = rcvm as i16;
    bh.mfeet = mfeet as i16;
    bh.polyt = polyt as i16;
    bh.vpol = vpol as i16;
    efwrite(bh.as_bytes(), &binaryfp);
    efclose(&binaryfp);
    exit(CWP_Exit());
}
