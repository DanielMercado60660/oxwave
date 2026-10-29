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

//! WPCCOMPRESS - compress a 2D section using wavelet packets.
//! Port of comp/dwpt/2d/main/wpccompress.c (Revision 1.5).
//! Credits: CWP: Tong Chen.

use oxwave_core::cio::as_bytes_mut;
use oxwave_core::prelude::*;
use oxwave_lib_wpc2::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " WPCCOMPRESS --- COMPRESS a 2D section using Wavelet Packets\t\t",
    "                                                                       ",
    " wpccompress < stdin n1= n2= [optional parameters ] > stdout          \t",
    "                                                                       ",
    " Required Parameters:                                                  ",
    " n1=                    number of samples in the 1st dimension\t\t",
    " n2=                    number of samples in the 2nd dimenstion\t",
    "                                                                       ",
    " Optional Parameters:                                                  ",
    " error=0.01              relative RMS allowed in compress\t\t",
    "                                                                       ",
    " Notes:                                                                ",
    "  This program is used to compress a 2D section. It compresses in both\t",
    "  directions, vertically and horizontally.\t\t\t\t",
    "                                                                       ",
    "  The parameter error is used control the allowable compression error,\t",
    "  and thus the compression ratio. The larger the error, the more \t",
    "  the more compression you can get. The amount of error depends on \t",
    "  the type of data and the application of the compression. From my \t",
    "  experience, error=0.01 is a safe choice even when the compressed data \t",
    "  are used for further processing. For some other applications, it \t",
    "  may be set higher to achieve larger compression.\t\t\t",
    "                                                                       ",
    " Caveats:\t\t\t\t\t\t\t\t",
    "  For the current implementation, the compressed data themselves are\t",
    "  NOT portable, i.e., the data compressed on one platform might not be\t",
    "  recognizable on another.\t\t\t\t\t\t",
    "                                                                       ",
];

fn main() {
    let argv = args();
    let mut error = 0.0f32;
    let mut n1 = 0;
    let mut n2 = 0;
    let time_used = 0.0f32;

    /* hook up getpars */
    initargs(&argv);
    requestdoc(0, SDOC);

    if getparint("n1", &mut n1) == 0 {
        err!("Must specify n1");
    }
    if getparint("n2", &mut n2) == 0 {
        err!("Must specify n2");
    }

    if getparfloat("error", &mut error) == 0 {
        error = 0.01;
    }

    checkpars();

    if n1 < 0 || n2 < 0 {
        // The original dereferences the failed allocation for negative dimensions.
        unsafe extern "C" {
            fn signal(sig: i32, handler: usize) -> usize;
            fn raise(sig: i32) -> i32;
        }
        unsafe {
            signal(11, 0);
            raise(11);
        }
    }
    if n1 == 0 || n2 == 0 {
        err!("%s: efread: fread given 0 item size", "subcalls.c");
    }
    // C alloc2float is uninitialized: a short read leaves unspecified samples.
    let mut f = vec![0.0f32; (n1 as usize).wrapping_mul(n2 as usize)];
    // C passes efread(f[0], n1*n2, FSIZE, stdin): four items of runtime size.
    // Read bytes through core so even a partial float counts toward a C item.
    // Convert the byte count to the original call's units for its diagnostic.
    let item_size = n1.wrapping_mul(n2) as usize;
    let bytes = as_bytes_mut(&mut f);
    let nread = fread(&mut bytes[..item_size.wrapping_mul(4)], &stdin()) / item_size;
    if nread != 4 && ferror(&stdin()) {
        err!(
            "%s: efread: fread only %d items of %d",
            "subcalls.c",
            nread,
            4
        );
    }

    // The supplied library retains the valid code buffer after realloc. C's
    // wpcCompress leaves wpc->code stale at both growth and final shrink.
    // Affected cases are excluded in goldens/wpccompress/EXCLUDE.txt.
    let wpc = wpcCompress(&f, n1 as usize, n2 as usize, error).unwrap();
    let nbytes = wpcSize(&wpc);

    let mut data = Vec::new();
    wpcWrite(&wpc, &mut data).unwrap();
    fwrite(&data, &stdout());

    fprintf!(stderr(), "# of bytes after compression = %d\n", nbytes);
    fprintf!(stderr(), "time used in compression is %f sec\n", time_used);

    exit(CWP_Exit());
}
