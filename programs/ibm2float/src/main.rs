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

//! IBM2FLOAT - convert IBM tape floats to native binary floats.
//! Port of par/main/data_conversion/ibm2float.c (Revision 1.4).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: CWP: John Stockwell, based on code by Jack K. Cohen.
//! Conversion: CWP: Brian Sumner; modifications by Toralf Foerster and M. J. Rutty.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " IBM2FLOAT - convert IBM tape FLOATS to native binary FLOATS\t",
    " \t\t\t\t\t\t\t\t",
    " ibm2float <stdin >stdout \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " endian=\tbyte order of your system (autodetected)\t",
    " outpar=/dev/tty output parameter file, contains the\t\t",
    "\t\t\tnumber of values (n=)\t\t\t",
    "\t\t       other choices for outpar are: /dev/tty, ",
    "\t\t       /dev/stderr, or a name of a disk file   ",
    " \t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " endian=1 (big endian) endian=0 (little endian) byte order \t",
    " You probably will not have to set this, as the byte order of  ",
    " your system is autodetected by the program. \t\t\t",
    " This program is usable for reading SEG Y files with the headers",
    " stripped off.\t\t\t\t\t\t\t",
];

// This program's local conversion differs from segyread's routine.
fn ibm_to_float(from: i32, endian: i32) -> i32 {
    let mut fconv = from;
    if endian == 0 { fconv = fconv.swap_bytes(); }
    if fconv != 0 {
        let mut fmant = 0x00ff_ffff & fconv;
        if fmant == 0 { warn!("mantissa is zero data may not be in IBM FLOAT Format !"); }
        let mut t = ((0x7f00_0000 & fconv) >> 22) - 130;
        while (fmant & 0x0080_0000) == 0 {
            t = t.wrapping_sub(1);
            fmant = fmant.wrapping_shl(1);
        }
        if t > 254 { fconv = (0x8000_0000u32 as i32 & fconv) | 0x7f7f_ffff; }
        else if t <= 0 { fconv = 0; }
        else { fconv = (0x8000_0000u32 as i32 & fconv) | (t << 23) | (0x007f_ffff & fmant); }
    }
    fconv
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    match filestat(STDOUT) {
        FileType::BadFileType => { warn!("stdout is illegal filetype"); pagedoc(SDOC); }
        FileType::Tty => { warn!("stdout can't be tty"); pagedoc(SDOC); }
        _ => {}
    }
    let mut endian = 0;
    if getparint("endian", &mut endian) == 0 {
        endian = if cfg!(target_endian = "big") { 1 } else { 0 };
    }
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let outparfp = efopen(&outpar, "w");
    checkpars();
    let mut x = [0i32; 1];
    let mut n1 = 0i32;
    while efread(&mut x, &stdin()) != 0 {
        let z = [ibm_to_float(x[0], endian)];
        fwrite(&z, &stdout());
        n1 = n1.wrapping_add(1);
    }
    fprintf!(outparfp, "n=%d\n", n1);
    exit(CWP_Exit());
}
