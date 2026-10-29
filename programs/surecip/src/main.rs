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
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! SURECIP - sum opposing offsets in prepared data.
//! Port of su/main/stacking/surecip.c (Revision 1.14).
//! Credits: SEP: Shuki Ronen; CWP: Jack Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                \t\t\t\t",
    " SURECIP - sum opposing offsets in prepared data (see below)\t",
    "                                \t\t\t\t",
    " surecip <stdin >stdout\t \t\t               \t",
    "                                \t\t\t\t",
    " Sum traces with equal positive and negative offsets (i.e. assume",
    " reciprocity holds). \t\t\t\t\t\t",
    "                                \t\t\t\t",
    " Usage:\t\t\t\t\t\t\t",
    "\tsuabshw <data >absdata\t\t\t\t\t",
    "\tsusort cdp offset <absdata | surecip >sumdata\t\t",
    "                                \t\t\t\t",
    " Note that this processing stream can be simply evoked by:\t",
    "                                \t\t\t\t",
    "\trecip data sumdata\t\t\t\t\t",
    "                                \t\t\t\t",
];

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    let mut cdpval = Value::default();
    let mut offval = Value::default();
    let mut cdpvalnew = Value::default();
    let mut offvalnew = Value::default();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut norm = 1;
    getparint("norm", &mut norm);
    checkpars();
    let cdpindex = getindex("cdp");
    let offindex = getindex("offset");
    let mut nsegy = gettr(&mut intrace);
    let nt = intrace.ns as usize;
    gethval(&intrace, cdpindex, &mut cdpval);
    gethval(&intrace, offindex, &mut offval);
    if nsegy != 0 {
        outtrace.hdr = intrace.hdr;
        outtrace.data[..nt].copy_from_slice(&intrace.data[..nt]);
    }
    let mut newtracl: i32 = 1;
    let mut fold: i32 = 1;
    while nsegy != 0 {
        nsegy = gettr(&mut intrace);
        gethval(&intrace, cdpindex, &mut cdpvalnew);
        gethval(&intrace, offindex, &mut offvalnew);
        if valcmp("l", cdpval, cdpvalnew) != 0 || nsegy == 0 {
            outtrace.nhs = fold as i16;
            outtrace.tracl = newtracl;
            newtracl = newtracl.wrapping_add(1);
            if outtrace.sx > outtrace.gx {
                let itmp = outtrace.sx;
                outtrace.sx = outtrace.gx;
                outtrace.gx = itmp;
            }
            if norm != 0 {
                let ffold = fold as f32;
                if fold != 1 {
                    for i in 0..nt { outtrace.data[i] /= ffold; }
                }
            }
            puttr(&mut outtrace);
            if nsegy != 0 {
                outtrace.hdr = intrace.hdr;
                outtrace.data[..nt].copy_from_slice(&intrace.data[..nt]);
            }
            fold = 1;
            cdpval = cdpvalnew;
            offval = offvalnew;
        } else if valcmp("l", offval, offvalnew) != 0 {
            outtrace.nhs = fold as i16;
            outtrace.tracl = newtracl;
            newtracl = newtracl.wrapping_add(1);
            if outtrace.sx > outtrace.gx {
                let itmp = outtrace.sx;
                outtrace.sx = outtrace.gx;
                outtrace.gx = itmp;
            }
            if norm != 0 {
                let ffold = fold as f32;
                if fold != 1 {
                    for i in 0..nt { outtrace.data[i] /= ffold; }
                }
            }
            puttr(&mut outtrace);
            outtrace.hdr = intrace.hdr;
            outtrace.data[..nt].copy_from_slice(&intrace.data[..nt]);
            fold = 1;
            offval = offvalnew;
        } else {
            for i in 0..nt { outtrace.data[i] += intrace.data[i]; }
            fold = fold.wrapping_add(1);
        }
    }
    exit(CWP_Exit());
}
