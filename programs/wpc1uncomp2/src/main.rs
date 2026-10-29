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


//! WPC1UNCOMP2 - uncompress a section trace by trace.
//! Port of comp/dwpt/1d/main/wpc1uncomp2.c (Revision 1.5).
//! Credits: Tong Chen.

use oxwave_core::prelude::*;
use oxwave_lib_wpc1::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " WPC1UNCOMP2 --- UNCOMPRESS a 2D seismic section, which has been\t",
    "  \t\t\tcompressed using Wavelet Packets\t\t",
    "                                                                       ",
    " wpc1uncomp2 < stdin > stdout               \t\t\t\t",
    "                                                                       ",
    " Required Parameters:                                                  ",
    " none                                                                  ",
    "                                                                       ",
    " Optional Parameters:                                                  ",
    " none                                                                  ",
    "                                                                       ",
    " Notes:                                                                ",
    "  No parameter is required for this program. All the information for\t",
    "  uncompression has been encoded in the header of the compressed data.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveats:\t\t\t\t\t\t\t\t",
    "  For the current implementation, the compressed data themselves are\t",
    "  NOT portable, i.e., the data compressed on one platform might not be\t",
    "  recognizable on another.\t\t\t\t\t\t",
    "                                                                       ",
];

// The C program dereferences NULL when its first record cannot be read.
unsafe extern "C" { fn raise(sig: i32) -> i32; fn signal(sig: i32, handler: usize) -> usize; }

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut input = std::io::stdin().lock();
    let mut wpc1 = match wpc1Read(&mut input) {
        Ok(w) => w,
        Err(_) => {
            unsafe { signal(11, 0); raise(11); }
            exit(1);
        }
    };
    let n1 = wpc1Numsmp(&wpc1);
    if n1 < 0 { fprintf!(stderr(), "Inconsistent data, n=%d\n", n1); }
    let mut f = vec![0.0f32; n1.max(0) as usize];
    let mut i = 0i32;
    while i >= 0 {
        if wpc1.numsmp > f.len() as i32 {
            // C assumes subsequent records have the first record's length.
            f.resize(wpc1.numsmp as usize, 0.0);
        }
        if wpc1Uncompress(&wpc1, &mut f) == 0 {
            fprintf!(stderr(), "Inconsistent data, trace # %d\n", i);
            break;
        }
        i = i.wrapping_add(1);
        fwrite(&f[..n1.max(0) as usize], &stdout());
        match wpc1Read(&mut input) {
            Ok(w) => wpc1 = w,
            Err(_) => break,
        }
    }
    exit(CWP_Exit());
}
