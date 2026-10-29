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

//! PSEPSI - add an EPSI preview to an EPS file.
//! Port of psplot/main/psepsi.c (Revision 1.10).
//! Author: Dave Hale, Colorado School of Mines, 03/01/92.

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " PSEPSI - add an EPSI formatted preview bitmap to an EPS file\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " psepsi <epsfile >epsifile\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note:\t\t\t\t\t\t\t\t\t",
    " This application requires\t\t\t\t\t\t",
    " (1) that gs (the Ghostscript interpreter) exist, and\t\t\t",
    " (2) that the input EPS file contain a BoundingBox and EndComments.\t",
    " Ghostscript is used to build the preview bitmap, which is then\t\t",
    " merged with the input EPS file to make the output EPSI file.\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

const PREVIEW: &[&str] = &[
    "% <width> <height> <filename> makepreviewbitmap",
    "%",
    "% Replaces the current device with a memory device with the",
    "% specified width and height.  When showpage is executed, this",
    "% memory device will write an EPSI preview bitmap to a file with",
    "% the specified filename.  The file will contain the following:",
    "%",
    "% %%BeginPreview: width height 1 nlines",
    "% <nlines of hex data representing the bitmap>",
    "% %%EndPreview",
    "%",
    "% where nlines is the number of lines of hexadecimal characters",
    "% representing the bitmap.  The bitmap may then be merged into",
    "% the prolog of an EPSI file.",
    "%",
    "% Author:  Dave Hale, Colorado School of Mines, 02/29/92",
    "",
    "% make dictionary",
    "/pdict 25 dict def",
    "pdict begin",
    "",
    "\t% save the showpage operator",
    "\t/realshowpage /showpage load def",
    "",
    "\t% define a monochrome palette with black and white reversed",
    "\t/monopalette [ 1 1 1 rgbcolor 0 0 0 rgbcolor ] def",
    "",
    "\t% the main procedure",
    "\t/makepreviewbitmap {",
    "",
    "\t\t% get file, height, and width",
    "\t\t(w) file /pfile exch def",
    "\t\t/pheight exch def",
    "\t\t/pwidth exch def",
    "",
    "\t\t% make memory device",
    "\t\t[1 0 0 -1 0 pheight] pwidth pheight monopalette",
    "\t\tmakeimagedevice /pdevice exch def",
    "",
    "\t\t% set the device (which does an erasepage)",
    "\t\tpdevice setdevice",
    "",
    "\t\t% replace definition of showpage",
    "\t\tuserdict /showpage {pdict begin writepreview end} bind put",
    "\t} def",
    "",
    "\t% procedure to write the preview bitmap to file",
    "\t/writepreview {",
    "",
    "\t\t% length of row in bytes an number of rows in device",
    "\t\t/lrow pwidth 7 add 8 idiv def",
    "\t\t/nrow pheight def",
    "",
    "\t\t% length of line and number of lines in preview bitmap",
    "\t\t/lline 38 def",
    "\t\t/nlineperrow lrow lline 1 sub add lline idiv def",
    "\t\t/nline nrow nlineperrow mul def",
    "",
    "\t\t% write BeginPreview comment",
    "\t\tpfile (%%BeginPreview: ) writestring",
    "\t\tpfile pwidth write==only pfile ( ) writestring",
    "\t\tpfile pheight write==only pfile ( ) writestring",
    "\t\tpfile (1 ) writestring",
    "\t\tpfile nline write==only pfile (\\n) writestring",
    "",
    "\t\t% string to contain one row from device",
    "\t\t/row lrow string def",
    "",
    "\t\t% zero last byte in row (which may contain unused bits)",
    "\t\trow lrow 1 sub 0 put",
    "",
    "\t\t% loop over rows in device",
    "\t\t0 1 nrow 1 sub {",
    "",
    "\t\t\t% get row from device",
    "\t\t\tpdevice exch row copyscanlines pop",
    "",
    "\t\t\t% loop over lines within row",
    "\t\t\t0 1 nlineperrow 1 sub {",
    "",
    "\t\t\t\t% each line is a comment",
    "\t\t\t\tpfile (%) writestring",
    "",
    "\t\t\t\t% index of first byte in line",
    "\t\t\t\tlline mul",
    "",
    "\t\t\t\t% number of bytes in line; may be less",
    "\t\t\t\t% than lline if this is last line in row",
    "\t\t\t\tdup lrow exch sub dup lline lt",
    "\t\t\t\texch lline ifelse",
    "",
    "\t\t\t\t% row (before index and count)",
    "\t\t\t\trow 3 1 roll",
    "",
    "\t\t\t\t% write line",
    "\t\t\t\tgetinterval pfile exch writehexstring",
    "\t\t\t\tpfile (\\n) writestring",
    "\t\t\t} for",
    "\t\t} for",
    "",
    "\t\t% write EndPreview comment",
    "\t\tpfile (%%EndPreview\\n) writestring",
    "",
    "\t\t% flush file and prepare for next page",
    "\t\tpfile flushfile",
    "\t\terasepage initgraphics",
    "\t} bind def",
    "end",
    "",
    "/makepreviewbitmap { pdict begin makepreviewbitmap end } bind def",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let gsfp = epopen("gs -dNODISPLAY - >/dev/null", "w");
    let infp = etmpfile();
    let (mut llx, mut lly, mut urx, mut ury) = (0i32, 0i32, 0i32, 0i32);
    let mut bboxfound = false;
    while let Some(line) = fgets(2048, &stdin()) {
        if !bboxfound && line.starts_with("%%BoundingBox") {
            if line.contains("atend") { continue; }
            // Skip the first %s field, then use the core's C decimal parser
            // for each of the four %d conversions.
            let s = line.split('\0').next().unwrap();
            let mut rest = &s[s.find(|c: char| c.is_ascii_whitespace()).unwrap_or(s.len())..];
            let mut values = [0i32; 4];
            for value in &mut values {
                let (number, consumed) = strtol(rest);
                if consumed == 0 { err!("bad BoundingBox in input EPS!"); }
                *value = number as i32;
                rest = &rest[consumed..];
            }
            (llx, lly, urx, ury) = (values[0], values[1], values[2], values[3]);
            bboxfound = true;
        }
        if fputs(&line, &infp) == EOF { err!("cannot copy input EPS to temporary file!"); }
    }
    if !bboxfound { err!("input EPS(?) does not have a BoundingBox!"); }
    // The C program passes an uninitialized buffer as the prefix.
    let pbname = "";
    let pbnamefp = temporary_filename(pbname);
    let pbfp = match fopen(&pbnamefp, "w+") {
        Some(fp) => fp,
        None => err!("cannot create temporary file for preview bitmap!"),
    };
    for line in PREVIEW { fprintf!(gsfp, "%s\n", *line); }
    fprintf!(gsfp, "%d %d (%s) makepreviewbitmap\n", urx.wrapping_sub(llx), ury.wrapping_sub(lly), pbname);
    fputs("gsave\n", &gsfp);
    fprintf!(gsfp, "%d %d translate\n", llx.wrapping_neg(), lly.wrapping_neg());
    efseeko(&infp, 0, SEEK_SET);
    while let Some(line) = fgets(2048, &infp) { fputs(&line, &gsfp); }
    fputs("grestore\n", &gsfp);
    fflush(&gsfp);
    epclose(&gsfp);
    efseeko(&infp, 0, SEEK_SET);
    efseeko(&pbfp, 0, SEEK_SET);
    let mut pbwritten = false;
    while let Some(line) = fgets(2048, &infp) {
        if !pbwritten && line.starts_with("%%EndComments") {
            while let Some(line2) = fgets(2048, &pbfp) {
                if fputs(&line2, &stdout()) == EOF { err!("cannot write EPSI to stdout!"); }
            }
            pbwritten = true;
        }
        if fputs(&line, &stdout()) == EOF { err!("cannot write EPSI to stdout!"); }
    }
    if std::fs::remove_file(pbname).is_err() { warn!("cannot remove temporary file %s!", pbname); }
    exit(CWP_Exit());
}
