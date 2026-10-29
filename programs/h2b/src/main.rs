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

//! H2B - convert hexadecimal bytes to floats. Port of par/main/data_conversion/h2b.c (Revision 1.11).
//! Credits: CWP: Craig Artley, with small changes by John Stockwell, July 1993.

use oxwave_core::prelude::*;
use std::ffi::{c_char, c_int, c_uint};

unsafe extern "C" {
    fn sscanf(s: *const c_char, format: *const c_char, ...) -> c_int;
}
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " H2B - convert 8 bit hexidecimal floats to binary\t\t",
    " \t\t\t\t\t\t\t\t",
    " h2b <stdin >stdout outpar=/dev/tty \t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \toutpar=/dev/tty output parameter file, contains the\t",
    "\t\t\tnumber of lines (n=)\t\t\t",
    " \t\t\tother choices for outpar are: /dev/tty,\t",
    " \t\t\t/dev/stderr, or a name of a disk file\t",
    " \t\t\t\t\t\t\t\t",
    " Note: this code may be used to recover binary data from PostScript",
    " bitmaps. To do this, strip away all parts of the PSfile that\t",
    " are not the actual hexidecimal bitmap and run through h2b.\t",
    " \t\t\t\t\t\t\t\t",
    " Note: that the binary file may need to be transposed using\t",
    " \"transp\" to appear to be the same as input data.\t\t",
    " \t\t\t\t\t\t\t\t",
    " Note:\toutput will be floats with the values 0-255\t\t",
];

fn main() {
    // BUFSIZ is 1024 in the macOS C library used for the original binary.
    const BUFSIZ: usize = 1024;
    let mut fz = [0.0f32; BUFSIZ];
    let mut z = [0u32; BUFSIZ];
    let mut line = [0u8; BUFSIZ + 2];
    let mut nz: i32 = 0;
    let infp = stdin();
    let outfp = stdout();

    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    match filestat(STDOUT) {
        FileType::BadFileType => { warn!("stdout is illegal filetype"); pagedoc(SDOC); }
        FileType::Tty => { warn!("stdout can't be tty"); pagedoc(SDOC); }
        _ => {}
    }

    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let outparfp = efopen(&outpar, "w");
    checkpars();

    loop {
        let mut len = 0;
        while len < BUFSIZ - 1 {
            let c = fgetc(&infp);
            if c == EOF { break; }
            line[len] = c as u8;
            len += 1;
            if c == b'\n' as i32 { break; }
        }
        if len == 0 { break; }
        line[len] = 0;

        let mut iz = 0usize;
        let mut lp = 0usize;
        while lp < BUFSIZ && iz < BUFSIZ {
            // sscanf supplies the C %2x handling of whitespace, signs and prefixes.
            let matched = unsafe {
                sscanf(line.as_ptr().add(lp).cast(), c"%2x".as_ptr(), &mut z[iz] as *mut c_uint)
            };
            if matched != 1 { break; }
            iz += 1;
            nz = nz.wrapping_add(1);
            lp += 2;
        }
        for jz in 0..iz { fz[jz] = 255u32.wrapping_sub(z[jz]) as f32; }
        fwrite(&fz[..iz], &outfp);
    }
    fprintf!(&outparfp, "total number of values=%d\n", nz);
    exit(CWP_Exit());
}
