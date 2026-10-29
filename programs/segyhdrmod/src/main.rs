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

//! SEGYHDRMOD - replace the text header on a SEGY file.
//! Port of su/main/data_conversion/segyhdrmod.c (Revision 1.5).
//! Credits (from the C source): Reginald H. Beardsley.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SEGYHDRMOD - replace the text header on a SEGY file\t\t",
    "\t\t\t\t\t\t\t\t",
    "   segyhdrmod text=file data=file\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "   Required parameters:\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "   text=      name of file containing new 3200 byte text header",
    "   data=      name of file containing SEGY data set\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " This program simply does a replacement of the content of the first",
    " 3200 bytes of the SEGY file with the contents of the file specified",
    " by the text= parameter. If the text header in the SEGY standard",
    " ebcdic format, the user will need to supply an ebcdic format file",
    " as the text=  as input file. A text file may be converted from",
    " ascii to ebcdic via:\t\t\t\t\t\t",
    "   dd if=ascii_filename of=ebcdic_filename conv=ebcdic ibs=3200 count=1",
    " or from ebcdic to ascii via:\t\t\t\t\t",
    "   dd if=ebcdic_filename of=ascii_filename ibs=3200 conv=ascii count=1",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut text = String::new();
    let mut data = String::new();
    let mut buf = [0u8; 3200];

    initargs(&argv);
    requestdoc(2, SDOC);
    if getparstring("text", &mut text) == 0 {
        err!("missing text header filename");
    }
    if getparstring("data", &mut data) == 0 {
        err!("missing data filename");
    }
    checkpars();

    let txtfp = match fopen(&text, "rb") {
        Some(fp) => fp,
        None => err!("unable to open %s", text),
    };
    let datfp = match fopen(&data, "rb+") {
        Some(fp) => fp,
        None => err!("unable to open %s", data),
    };
    if fread(&mut buf, &txtfp) != buf.len() {
        err!("unable to read new text header");
    }
    if fwrite(&buf, &datfp) != buf.len() {
        err!("write of new text header failed!!!!!");
    }
    exit(CWP_Exit());
}
