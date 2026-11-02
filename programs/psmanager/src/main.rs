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

//! Original source: Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! PSMANAGER - printer manager. Port of psplot/main/psmanager.c (Revision 1.6).
//! Author: John Stockwell, June 1995, October 1997.
//! Reverse engineered from NeXTStep printer manager output; reference: Adobe Systems
//! PostScript Printer Description File Format Specification, version 4.2.

use oxwave_core::prelude::*;
use oxwave_core::cio::fputs_bytes;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " PSMANAGER - printer MANAGER for HP 4MV and HP 5Si Mx Laserjet ",
    "                PostScript printing\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    "   psmanager < stdin  [optional parameters] > stdout \t\t",
    " \t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t",
    "  none \t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " papersize=0\tpaper size  (US Letter default)\t\t\t",
    " \t\t=1       US Legal\t\t\t\t",
    " \t\t=2\t A4\t\t\t\t\t",
    " \t\t=3     \t 11x17\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " orient=0\tpaper orientation (Portrait default)\t\t",
    "  \t\t=1   \tLandscape\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " tray=3        printing tray (Bottom tray default)\t\t",
    "  \t\t=1\ttray 1 (multipurpose slot)\t\t",
    "  \t\t=2\ttray 2 \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " manual=0\tno manual feed \t\t\t\t\t",
    "  \t\t=1     (Manual Feed)\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " media=0\tregular paper\t\t\t\t\t",
    "  \t\t=1     Transparency\t\t\t\t",
    "  \t\t=2     Letterhead\t\t\t\t",
    "  \t\t=3     Card Stock\t\t\t\t",
    "  \t\t=4     Bond\t\t\t\t\t",
    "  \t\t=5     Labels\t\t\t\t\t",
    "  \t\t=6     Prepunched\t\t\t\t",
    "  \t\t=7     Recyled\t\t\t\t\t",
    "  \t\t=8     Preprinted\t\t\t\t",
    "  \t\t=9     Color (printing on colored paper)\t",
    " \t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t",
    " The option manual=1 implies tray=1. The media options apply\t",
    " only to the HP LaserJet 5Si MX model printer.\t\t\t",
    "  \t\t\t\t\t\t\t\t",
    " Examples: \t\t\t\t\t\t\t",
    "   overheads:\t\t\t\t\t\t\t",
    "    psmanager <  postscript_file manual=1 media=1 | lpr\t",
    "   labels:\t\t\t\t\t\t\t",
    "    psmanager <  postscript_file manual=1 media=5 | lpr\t",
    " \t\t\t\t\t\t\t\t",
];

const PROLOG: &str = "/__CWPbdef {\n    1 index where {\n        pop pop pop\n    } {\n        bind def\n    } ifelse\n} bind def /__CWPRectPath {\n    4 2 roll moveto 1 index 0 rlineto 0 exch rlineto neg 0 rlineto closepath\n} __CWPbdef /__CWPProcessRectArgs {\n    1 index type /arraytype eq {\n        exch 0 4 2 index length 1 sub {\n            dup 3 add 1 exch {\n                1 index exch get exch\n            } for 5 1 roll 5 index exec\n        } for pop pop\n    } {\n        exec\n    } ifelse\n} __CWPbdef /rectfill {\n    gsave newpath {\n        __CWPRectPath fill\n    } __CWPProcessRectArgs grestore\n} __CWPbdef /rectclip {\n    newpath {\n        __CWPRectPath\n    } __CWPProcessRectArgs clip newpath\n} __CWPbdef\n%%EndComments\n\ngsave\n-20 -28 translate\n /__CWPbasematrix matrix currentmatrix def\ngrestore\n%%EndProlog\n";

// The first field is discarded by the original sscanf formats.
fn scan_values(s: &str, floats: bool, values: &mut [f64]) -> usize {
    let s = s.trim_start_matches(|c: char| c.is_ascii_whitespace());
    let first = s.find(|c: char| c.is_ascii_whitespace()).unwrap_or(s.len());
    let mut rest = &s[first..];
    let mut count = 0;
    for value in values {
        let (v, used) = if floats {
            strtod(rest)
        } else {
            let (v, used) = strtol(rest);
            (v as i32 as f64, used)
        };
        if used == 0 { break; }
        *value = v;
        rest = &rest[used..];
        count += 1;
    }
    count
}

fn main() {
    let argv = args(); initargs(&argv); requestdoc(1, SDOC);
    let (mut papersize, mut orient, mut tray, mut manual, mut media) = (0,0,3,0,0);
    getparint("papersize", &mut papersize); getparint("orient", &mut orient);
    getparint("tray", &mut tray); getparint("manual", &mut manual); getparint("media", &mut media); checkpars();
    if manual != 0 { tray = 1; }
    let (slot, position) = match tray { 1 => ("MP",3), 2 => ("Cassette",0), 3 => ("LargeCapacity",1), _ => err!("Unsupported tray, %d", tray) };
    let direction = match orient { 0 => "Portrait", 1 => "Landscape", _ => err!("Unsupported orientation, %d", orient) };
    let (size, mut dim) = match papersize { 0 => ("Letter",[612,792]), 1 => ("Legal",[612,1008]), 2 => ("A4",[595,842]), 3 => ("Tabloid",[792,1224]), _ => err!("Unsupported papersize selection, %d", papersize) };
    let dims = dim; if orient != 0 { dim.swap(0,1); }
    let mediatype = match media { 0 => "None", 1 => "Transparency", 2 => "Letterhead", 3 => "Card Stock", 4 => "Bond", 5 => "Labels", 6 => "Prepunched", 7 => "Recycled", 8 => "Preprinted", 9 => "Color", _ => err!("Unsupported media selection, %d", media) };
    let traystr = format!("<</DeferredMediaSelection true /MediaPosition {}>> setpagedevice", position);
    let mediastr = if media == 0 { "<</DeferredMediaSelection true /MediaType null>> setpagedevice".to_string() } else { format!("<</DeferredMediaSelection true /MediaType ({})>> setpagedevice", mediatype) };
    let pagestr = format!("<</DeferredMediaSelection true /PageRegion [{} {}] /ImagingBBox null>> setpagedevice", dims[0], dims[1]);
    let mut lines: Vec<Vec<u8>> = Vec::new(); let mut line = Vec::new();
    let (mut llx, mut lly, mut urx, mut ury, mut pages) = (0i32,0i32,0i32,0i32,1i32);
    let (mut fllx, mut flly, mut furx, mut fury) = (0f32,0f32,0f32,0f32);
    loop {
        let c = fgetc(&stdin());
        if c != -1 { line.push(c as u8); }
        if c == -1 || c == 10 || line.len() == 1023 {
            if !line.is_empty() {
                let end = line.iter().position(|&b| b == 0).unwrap_or(line.len());
                let visible = &line[..end];
                if visible.first() == Some(&b'%') {
                    let s = String::from_utf8_lossy(visible);
                    if s.contains("%%BoundingBox:") {
                        if s.contains("atend") { line.clear(); if c == -1 { break; } else { continue; } }
                        let mut ints = [llx as f64, lly as f64, urx as f64, ury as f64];
                        let n = scan_values(&s, false, &mut ints);
                        llx = ints[0] as i32; lly = ints[1] as i32; urx = ints[2] as i32; ury = ints[3] as i32;
                        if n != 4 {
                            let mut floats = [fllx as f64, flly as f64, furx as f64, fury as f64];
                            let n = scan_values(&s, true, &mut floats);
                            fllx = floats[0] as f32; flly = floats[1] as f32; furx = floats[2] as f32; fury = floats[3] as f32;
                            if n == 4 {
                                llx = nint(fllx as f64); lly = nint(flly as f64); urx = nint(furx as f64); ury = nint(fury as f64);
                            } else { err!("Error reading BoundingBox!"); }
                        }
                    }
                    if s.contains("%%Pages:") {
                        let mut value = [pages as f64];
                        scan_values(&s, false, &mut value);
                        pages = value[0] as i32;
                    }
                }
                lines.push(visible.to_vec()); line.clear();
            }
            if c == -1 { break; }
        }
    }
    let pllx = ((dim[0] - urx + llx) as f64 / 2.0) as f32;
    let plly = ((dim[1] - ury + lly) as f64 / 2.0) as f32;
    let ipllx = pllx as i32; let iplly = plly as i32; let ipurx = dim[0]-ipllx; let ipury = dim[1]-iplly;
    printf!("%%!PS-Adobe-2.0\n%%%%Title:\n%%%%Creator: psmanager\n%%%%CreationDate\n%%%%For: psmanager user\n%%%%DocumentFonts: (atend)\n");
    if pages == 1 { printf!("%%%%Pages: (atend) 1\n"); } else { printf!("%%%%Pages: %d 1\n", pages); }
    printf!("%%%%BoundingBox: (atend)\n%%%%DocumentPaperSizes: %s\n%%%%Orientation: %s\n", size, direction);
    if pages == 1 { fputs(PROLOG, &stdout()); } else { printf!("%%%%EndComments\n"); }
    printf!("%%%%BeginSetup\n");
    if manual != 0 { printf!("%%%%BeginFeature: *ManualFeed True\n\n\n<</ManualFeed true>> setpagedevice\n%%%%EndFeature\n"); }
    printf!("%%%%BeginFeature: *InputSlot %s\n\n%s\n%%%%EndFeature\n", slot, traystr);
    printf!("%%%%BeginFeature: *MediaType %s\n\n%s\n%%%%EndFeature\n", mediatype, mediastr);
    printf!("%%%%BeginFeature: *PageRegion %s\n\n%s\n%%%%EndFeature\n", size, pagestr);
    printf!("%%%%Feature: *Resolution 600dpi\n%%%%EndFeature\n%%%%EndSetup\n\n\n");
    if pages == 1 { printf!("%%%%Page: 1 1\n%%%%PageBoundingBox: %d %d %d %d\n%%%%PageFonts: (atend)\n%%%%BeginPageSetup\n", ipllx, iplly, ipurx, ipury); }
    printf!("%%%%BeginFeature: *InputSlot %s\n\n%s\n%%%%EndFeature\n", slot, traystr);
    if manual != 0 { printf!("%%%%BeginFeature: *ManualFeed True\n\n\n<</ManualFeed true>> setpagedevice%%%%EndFeature\n"); }
    printf!("%%%%BeginFeature: *PageRegion %s\n\n\n", size);
    if pages == 1 {
        printf!("%s\n%%%%EndFeature\n/__CWPsheetsavetoken save def\n", pagestr);
        if orient != 0 { printf!("-90 rotate\n%d 0 translate\n", -dim[0]); }
        printf!("%3.1f %3.1f translate\ngsave\n-20 -28 translate\n /__CWPbasematrix matrix currentmatrix def\ngrestore\ngsave\n", pllx as f64, plly as f64);
        printf!("0 0 %d %d rectclip\n%d %d translate\n0 0 translate\n%%%%EndPageSetup\ngsave\n", urx-llx, ury-lly, -llx, -lly);
        printf!("%d %d %d %d rectclip\ngsave 1 setgray %d %d %d %d rectfill grestore\n", llx, lly, urx-llx, ury-lly, llx, lly, urx-llx, ury-lly);
        printf!("userdict /_CWPPreviewSaveObject0 systemdict /save get exec put\n/showpage {\n\n} def\n%%%%BeginFile:\n");
    }
    // fputs writes each input chunk to the tempfile; the subsequent fgets
    // reads it again, so a truncated chunk can join the following line.
    let contents: Vec<u8> = lines.concat();
    let mut offset = 0;
    let mut count = 0;
    while offset < contents.len() {
        let limit = (offset + 1023).min(contents.len());
        let end = contents[offset..limit].iter().position(|&c| c == b'\n')
            .map_or(limit, |i| offset + i + 1);
        let line = &contents[offset..end];
        if pages != 1 && line.windows(7).any(|w| w == b"%%Page:") {
            count += 1; printf!("clear userdict /_CWPPreviewSaveObject0 known {\n   _CWPPreviewSaveObject0 systemdict /restore get exec\n} if\n");
            fputs_bytes(line, &stdout()); printf!("userdict /_CWPPreviewSaveObject0 systemdict /save get exec put\n%%Page: %d %d\n", count, count);
        } else { fputs_bytes(line, &stdout()); }
        offset = end;
    }
    if pages == 1 {
        printf!("%%%%EndFile\nclear userdict /_CWPPreviewSaveObject0 known {\n    _CWPPreviewSaveObject0 systemdict /restore get exec\n} if\ngrestore\ngrestore\nshowpage\n__CWPsheetsavetoken restore\n%%%%PageTrailer\n%%%%Trailer\n%%%%Pages: 1 1\n");
        if orient != 0 { printf!("%%%%BoundingBox:%d %d %d %d\n", ipllx, iplly, ipury, ipurx); } else { printf!("%%%%BoundingBox:%d %d %d %d\n", ipllx, iplly, ipurx, ipury); }
    }
    exit(CWP_Exit());
}
