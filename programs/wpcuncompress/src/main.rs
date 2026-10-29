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

//! WPCUNCOMPRESS - uncompress a 2D section.
//! Port of comp/dwpt/2d/main/wpcuncompress.c (Revision 1.4).
//! Credits (from the C source): CWP: Tong Chen.

use oxwave_core::prelude::*;
use oxwave_lib_wpc2::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " WPCUNCOMPRESS --- UNCOMPRESS a 2D section \t\t\t\t",
    "                                                                       ",
    " wpc1uncompress < stdin > stdout            \t\t\t\t",
    "                                                                       ",
    " Required Parameters:                                                  ",
    "                                                                       ",
    " Optional Parameters:                                                  ",
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

unsafe extern "C" {
    fn signal(sig: i32, handler: usize) -> usize;
    fn raise(sig: i32) -> i32;
}

fn segfault() -> ! {
    fflush(&stderr());
    // Restore SIG_DFL: Rust installs a SIGSEGV handler for stack overflow.
    unsafe {
        signal(11, 0);
        raise(11);
    }
    exit(EXIT_FAILURE);
}

fn main() {
    let argv = args();
    let time_used = 0.0f32;

    initargs(&argv);
    requestdoc(1, SDOC);

    let wpc = match wpcRead(&mut std::io::stdin().lock()) {
        Ok(wpc) => wpc,
        Err(_) => {
            // The C reader leaves its header uninitialized on empty/short
            // input; this build prints this diagnostic and dies with SIGSEGV.
            fprintf!(stderr(), "Inconsistent data\n");
            segfault();
        }
    };

    let n1 = wpcN1(&wpc);
    let n2 = wpcN2(&wpc);

    if n1 < 0 || n2 < 0 {
        fprintf!(stderr(), "Inconsistent data\n");
    }

    // alloc2float supplies contiguous storage, used through f[0] in C.
    // With a negative row count (or negative width and nonzero rows), the
    // original allocation fails and its unchecked pointer access faults.
    if n2 < 0 || (n1 < 0 && n2 > 0) {
        segfault();
    }
    let mut f = vec![0.0f32; n1.wrapping_mul(n2) as usize];

    if wpcUncompress(&wpc, &mut f) == 0 {
        fprintf!(stderr(), "Inconsistent data\n");
    }

    fwrite(&f, &stdout());
    fprintf!(
        stderr(),
        "time used in uncompression is %f sec \n",
        time_used
    );

    exit(CWP_Exit());
}
