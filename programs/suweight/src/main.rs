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
//! SUWEIGHT - weight traces by header parameter.
//! Port of su/main/amplitudes/suweight.c (Revision 1.8).
//! Credits: CWP: John Stockwell; inv option: Garry Perratt (Geocon).

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUWEIGHT - weight traces by header parameter, such as offset\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   suweight < stdin > stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t   \t\t",
    "   <none>\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t   \t\t",
    " key=offset\tkeyword of header field to weight traces by \t\t",
    " a=1.0\t\tconstant weighting parameter (see notes below)\t\t",
    " b=.0005\tvariable weighting parameter (see notes below)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "... or use values of a header field for the weighting ...\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " key2=\t\tkeyword of header field to draw weights from\t\t",
    " scale=.0001\tscale factor to apply to header field values\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " inv=0\t\tweight by header value\t\t\t \t\t",
    " \t\t=1 weight by inverse of header value\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t \t",
    " This code is initially written with offset weighting in mind, but may\t",
    " be used for other, user-specified schemes.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The rationale for this program is to correct for unwanted linear\t",
    " amplitude trends with offset prior to either CMP stacking or AVO work.",
    " The code has to be edited should other functions of a keyword be required.",
    "\t\t\t\t\t\t\t\t\t",
    " The default form of the weighting is to multiply the amplitudes of the",
    " traces by a factor of:    ( a + b*keyword).\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " If key2=  header field is  set then this program uses the weighting\t",
    " values read from that header field, instead. Note, that because most\t",
    " header fields are integers, the scale=.0001 permits 10001 in the header",
    " to represent 1.0001.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " To see the list of available keywords, type:    sukeyword  -o  <CR>\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut key = String::new();
    let mut a = 1.0f32;
    let mut b = 0.0005f32;
    if getparstring("key", &mut key) == 0 { key = "offset".to_string(); }
    if getparfloat("a", &mut a) == 0 { a = 1.0; }
    if getparfloat("b", &mut b) == 0 { b = 0.0005; }
    let type_ = hdtype(&key);
    let index = getindex(&key);

    let mut key2 = String::new();
    let is_key2 = getparstring("key2", &mut key2) != 0;
    let mut type2 = "";
    let mut index2 = 0;
    let mut scale = 0.0001f32;
    if is_key2 {
        type2 = hdtype(&key2);
        index2 = getindex(&key2);
        if getparfloat("scale", &mut scale) == 0 { scale = 0.0001; }
    }
    let mut inv = 0i32;
    if getparint("inv", &mut inv) == 0 { inv = 0; }
    checkpars();

    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let ns = tr.ns as usize;
    loop {
        let mut val = Value::default();
        if !is_key2 {
            gethval(&tr, index, &mut val);
            let fval = vtof(type_, val);
            if inv == 0 {
                for i in 0..ns { tr.data[i] *= a + fval * b; }
            } else {
                for i in 0..ns { tr.data[i] /= a + fval * b; }
            }
        } else {
            gethval(&tr, index2, &mut val);
            let mut fval = vtof(type2, val);
            fval *= scale;
            if inv == 0 {
                for i in 0..ns { tr.data[i] *= fval; }
            } else {
                for i in 0..ns { tr.data[i] /= fval; }
            }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
