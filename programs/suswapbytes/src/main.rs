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

//! SUSWAPBYTES - swap trace byte order.
//! Port of su/main/data_conversion/suswapbytes.c (Revision 1.26).
//! Original source: Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: CWP: John Stockwell; Institute fur Geophysik, Hamburg: Jens Hartmann (June 1993).

use oxwave_core::prelude::*;
use oxwave_core::segy::swap_float_4;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSWAPBYTES - SWAP the BYTES in SU data to convert data from big endian",
    "               to little endian byte order, and vice versa\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suswapbytes < stdin [optional parameter] > sdtout\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " \tformat=0\t\tforeign to native\t\t\t",
    " \t\t\t\t=1 native to foreign\t\t\t",
    "\tswaphdr=1\t\tswap the header byte order\t\t",
    " \t\t\t\t=0 do not change the header byte order\t",
    "\tswapdata=1\t\tswap the data byte order\t\t",
    " \t\t\t\t=0 do not change the data byte order\t",
    " \tns=from header\t\tif ns not set in header, must be set by hand",
    " Notes:\t\t\t\t\t\t\t\t",
    "  The 'native'\tendian is the endian (byte order) of the machine you are",
    "  running this program on. The 'foreign' endian is the opposite byte order.",
    "   \t\t\t\t\t\t\t\t\t",
    " Examples of big endian machines are: IBM RS6000, SUN, NeXT\t\t",
    " Examples of little endian machines are: PCs, DEC\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat: this code has not been tested on DEC\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    match filestat(STDOUT) {
        FileType::Tty => err!("stdout can't be tty"),
        FileType::Directory => err!("stdout must be a file, not a directory"),
        FileType::BadFileType => err!("stdout is illegal filetype"),
        _ => {}
    }
    let mut format = 0;
    let mut swaphdr = 1;
    let mut swapdata = 1;
    getparint("format", &mut format);
    getparint("swaphdr", &mut swaphdr);
    getparint("swapdata", &mut swapdata);
    if format == 0 {
        efread(tr.hdr.as_bytes_mut(), &stdin());
        if swaphdr != 0 {
            for i in 0..SU_NKEYS { swaphval(&mut tr.hdr, i as i32); }
        }
        let mut ns = 0;
        if getparint("ns", &mut ns) == 0 { ns = tr.ns as i32; }
        checkpars();
        if ns == 0 { err!("ns not set in header, set ns by hand!"); }
        // Retain the guard against a negative count becoming a huge allocation.
        // The first trace's separate data allocation is not limited by SU_NFLTS.
        if ns < 0 {
            err!("ns=%d must be between 1 and %d", ns, SU_NFLTS as i32);
        }
        let databytes = ns as usize * 4;
        let mut data = vec![0.0f32; databytes];
        efread(&mut oxwave_core::cio::as_bytes_mut(&mut data)[..databytes], &stdin());
        if swapdata != 0 {
            for i in 0..ns { swap_float_4(&mut data[i as usize]); }
        }
        efwrite(tr.hdr.as_bytes(), &stdout());
        efwrite(&oxwave_core::cio::as_bytes(&data)[..databytes], &stdout());
        let trbytes = HDRBYTES + databytes;
        // The C trace is static; retain bytes not overwritten by a short read.
        let mut raw = vec![0u8; trbytes];
        raw[..HDRBYTES].copy_from_slice(tr.hdr.as_bytes());
        while efread(&mut raw, &stdin()) != 0 {
            // Only later traces use C's fixed-size segy object. Oversized later
            // records overrun it in C; keep the safety guard on this path only.
            if ns as usize > SU_NFLTS {
                err!("ns=%d must be between 1 and %d", ns, SU_NFLTS as i32);
            }
            tr.hdr.as_bytes_mut().copy_from_slice(&raw[..HDRBYTES]);
            if swaphdr != 0 {
                for i in 0..SU_NKEYS { swaphval(&mut tr.hdr, i as i32); }
            }
            raw[..HDRBYTES].copy_from_slice(tr.hdr.as_bytes());
            if swapdata != 0 {
                tr.data_bytes_mut()[..databytes].copy_from_slice(&raw[HDRBYTES..]);
                for i in 0..ns { swap_float_4(&mut tr.data[i as usize]); }
                raw[HDRBYTES..].copy_from_slice(&tr.data_bytes()[..databytes]);
            }
            efwrite(&raw, &stdout());
        }
    } else {
        while gettr(&mut tr) != 0 {
            let ns = tr.ns as usize;
            if swapdata != 0 {
                for i in 0..ns { swap_float_4(&mut tr.data[i]); }
            }
            if swaphdr != 0 {
                for i in 0..SU_NKEYS { swaphval(&mut tr.hdr, i as i32); }
            }
            efwrite(tr.hdr.as_bytes(), &stdout());
            efwrite(&tr.data_bytes()[..4 * ns], &stdout());
        }
    }
    exit(CWP_Exit());
}
