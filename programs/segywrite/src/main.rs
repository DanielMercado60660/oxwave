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
//! SEGYWRITE - port of segywrite.c (Revision 1.55).
//! Credits: SEP: Einar Kjartansson; CWP: Jack, Brian, Chris, John Stockwell.
//! IBM conversion: CWP: Brian Sumner.
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_core::segy::{Bhed, float_to_ibm, zebc};
use std::io::Write;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SEGYWRITE - write an SEG-Y tape\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " segywrite <stdin tape=\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\ttape=\t\ttape device to use (see sudoc segyread)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t\t",
    " verbose=0\tsilent operation\t\t\t\t",
    "\t\t=1 ; echo every 'vblock' traces\t\t\t",
    " vblock=50\techo every 'vblock' traces under verbose option ",
    " buff=1\t\tfor buffered device (9-track reel tape drive)\t",
    "\t\t=0 possibly useful for 8mm EXABYTE drive\t",
    " conv=1\t\t=0 don't convert to IBM format\t\t\t",
    " ebcdic=1\tconvert text header to ebcdic, =0 leave as ascii\t",
    " hfile=header\tebcdic card image header file\t\t\t",
    " bfile=binary\tbinary header file\t\t\t\t",
    " trmin=1 first trace to write\t\t\t\t\t",
    " trmax=INT_MAX  last trace to write\t\t\t       ",
    " endian=(autodetected)\t=1 for big-endian and =0 for little-endian byte order",
    " errmax=0\tallowable number of consecutive tape IO errors\t",
    " format=\t\toverride value of format in binary header file\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: The header files may be created with  'segyhdrs'.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: For buff=1 (default) tape is accessed with 'write', for buff=0\t",
    "\ttape is accessed with fwrite. Try the default setting of buff=1 ",
    "\tfor all tape types.\t\t\t\t\t\t",
    " Caveat: may be slow on an 8mm streaming (EXABYTE) tapedrive\t\t",
    " Warning: segyread or segywrite to 8mm tape is fragile. Allow time\t",
    "\t   between successive reads and writes.\t\t\t\t",
    " Precaution: make sure tapedrive is set to read/write variable blocksize",
    "\t   tapefiles.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For more information, type:\tsudoc <segywrite>\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args(); initargs(&argv); requestdoc(1, SDOC);
    let mut tape = String::new(); mustgetparstring("tape", &mut tape);
    let mut hfile = String::from("header"); let mut bfile = String::from("binary");
    getparstring("hfile", &mut hfile); getparstring("bfile", &mut bfile);
    let (mut trmin, mut trmax, mut verbose, mut vblock) = (1, i32::MAX, 0, 50);
    let (mut buff, mut conv, mut ebcdic, mut endian, mut errmax, mut format) = (1, 1, 1, if cfg!(target_endian="big") {1} else {0}, 0, 0);
    getparint("trmin", &mut trmin); getparint("trmax", &mut trmax);
    getparint("verbose", &mut verbose); getparint("vblock", &mut vblock);
    getparint("buff", &mut buff); getparint("conv", &mut conv);
    getparint("ebcdic", &mut ebcdic); getparint("endian", &mut endian);
    getparint("errmax", &mut errmax); let format_set = getparint("format", &mut format) != 0;
    checkpars();
    if trmin < 1 || trmax < 1 || trmax < trmin { err!("bad trmin/trmax values, trmin = %d, trmax = %d", trmin, trmax); }
    let mut tr = Trace::new(); gettr(&mut tr);
    let mut fd = if buff != 0 { Some(std::fs::OpenOptions::new().write(true).create(true).truncate(true).open(&tape).unwrap_or_else(|_| syserr!("%s: eopen: open failed", "syscalls.c"))) } else { None };
    let fp = if buff == 0 { Some(efopen(&tape, "w")) } else { None };
    if verbose != 0 { warn!("tape opened successfully "); }
    let headerfp = efopen(&hfile, "r");
    if verbose != 0 { warn!("header file opened successfully"); }
    let binaryfp = efopen(&bfile, "r");
    if verbose != 0 { warn!("binary file opened successfully"); }
    let mut ebcbuf = [0u8; 3201];
    if ebcdic == 1 {
        for i in (0..3200).step_by(80) {
            let mut j = 0;
            while j < 80 { let c = fgetc(&headerfp); if c == EOF { break; } ebcbuf[i+j] = c as u8; j += 1; if c == 10 { break; } }
            let length = ebcbuf[i..i+j].iter().position(|&c| c == 0).unwrap_or(j);
            ebcbuf[i+length] = b' ';
        }
        let source = ebcbuf; zebc(&source[..3200], &mut ebcbuf[..3200], 3200);
    } else { fread(&mut ebcbuf[..3200], &headerfp); }
    efclose(&headerfp);
    if verbose != 0 { warn!("header file closed successfully"); }
    let mut errcount = 0;
    write_tape(&mut fd, &fp, &ebcbuf[..3200], &mut errcount, errmax, verbose, "ebcdic header", 0);
    let mut bhbytes = [0u8; 400];
    let nread = fread(&mut bhbytes, &binaryfp);
    // xdrbhdrsub continues after a failed read, retaining its last utemp
    // for each subsequent pair of shorts.
    if nread >= 16 && nread < bhbytes.len() {
        let end = nread / 4 * 4;
        let last: [u8; 4] = bhbytes[end-4..end].try_into().unwrap();
        for i in (end..bhbytes.len()).step_by(4) {
            bhbytes[i..i+4].copy_from_slice(&last);
        }
    }
    let mut bh = Bhed::from_xdr(&bhbytes);
    if format_set { bh.format = format as i16; } bh.ntrpr = 1;
    let ns = bh.hns as usize;
    if ns == 0 { err!("bh.hns not set in binary header"); }
    write_tape(&mut fd, &fp, &bh.to_xdr(), &mut errcount, errmax, verbose, "binary header", 0);
    let mut itr = 0i32;
    loop {
        itr += 1; tr.tracr = itr;
        if tr.ns as usize != ns { warn!("conflict: tr.ns = %d, bh.ns = %d: trace %d", tr.ns, ns, itr); }
        if itr >= trmin {
            let mut trbuf = Vec::with_capacity(240 + 4*ns); trbuf.extend_from_slice(&tr.hdr.to_xdr());
            if conv != 0 { float_to_ibm(&mut tr.data[..ns], ns, endian); for x in &tr.data[..ns] { trbuf.extend_from_slice(&x.to_bits().to_ne_bytes()); } }
            else { for x in &tr.data[..ns] { trbuf.extend_from_slice(&x.to_bits().to_be_bytes()); } }
            write_tape(&mut fd, &fp, &trbuf, &mut errcount, errmax, verbose, "trace", itr);
            // C's remainder by zero yields zero on the target architecture.
            if verbose != 0 && itr.checked_rem(vblock).unwrap_or(0) == 0 { warn!(" %d traces written to tape", itr); }
        }
        if gettr(&mut tr) == 0 || itr >= trmax { break; }
    }
    if let Some(fp) = fp { efclose(&fp); } else { drop(fd); }
    if verbose != 0 { warn!("tape closed successfully"); }
    efclose(&binaryfp);
    if verbose != 0 { warn!("binary file closed successfully"); }
    exit(CWP_Exit());
}

fn write_tape(fd: &mut Option<std::fs::File>, fp: &Option<CFile>, data: &[u8], errcount: &mut i32, errmax: i32, verbose: i32, what: &str, itr: i32) {
    let failed = if let Some(fd) = fd { fd.write(data).map_or(true, |n| n != data.len()) } else { let fp = fp.as_ref().unwrap(); fwrite(data, fp); ferror(fp) };
    if failed {
        if verbose != 0 {
            if what == "trace" { warn!("tape write error on trace %d", itr); }
            else if what == "ebcdic header" { warn!("tape write error on ebcdic header"); }
            else { warn!("tape write error on binary header"); }
        }
        *errcount += 1; if *errcount > errmax { err!("exceeded maximum io errors"); }
        if let Some(fp) = fp { clearerr(fp); }
    } else { *errcount = 0; }
}
