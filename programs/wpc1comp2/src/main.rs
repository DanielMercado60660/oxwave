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

//! WPC1COMP2 - compress a section trace by trace.
//! Port of comp/dwpt/1d/main/wpc1comp2.c (Revision 1.5).
//! Credits (from the C source): CWP: Tong Chen.

use oxwave_core::prelude::*;
use oxwave_lib_wpc1::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " WPC1COMP2 --- COMPress a 2D seismic section trace-by-trace using \t",
    "\t\tWavelet Packets\t\t\t\t\t\t",
    "                                                                       ",
    " wpc1comp2 < stdin n1= [optional parameters ] > stdout              \t",
    "                                                                       ",
    " Required Parameters:                                                  ",
    " n1=                    number of samples in each trace\t\t",
    "                                                                       ",
    " Optional Parameters:                                                  ",
    " error=0.01              relative RMS allowed in compress\t\t",
    "                                                                       ",
    " Notes:                                                                ",
    "  This program is used to compress a 2D section using 1D method.\t",
    "  It is not a true 2D compression algorithm, since the later diretion\t",
    "  is not compressed at all. It can be used in situations where random\t",
    "  accessing of each trace is desirable, or when there are very few\t",
    "  traces and normal 2D compression algorithm will not compress.\t",
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
    initargs(&argv);
    requestdoc(1, SDOC);

    // C leaves n1 uninitialized when omitted (including n1=). With empty
    // input this build exits quietly; use the same path as a failed
    // allocation for a negative sample count.
    let mut n1 = -1i32;
    let mut error = 0.0f32;
    getparint("n1", &mut n1);
    if getparfloat("error", &mut error) == 0 {
        error = 0.01;
    }
    checkpars();

    // A negative n1 produces a failed allocation in C; fread then returns
    // zero on empty input instead of entering the compression loop.
    if n1 < 0 {
        exit(CWP_Exit());
    }
    let mut f = vec![0.0f32; n1 as usize];
    let mut nbytes = 0;
    let mut output = std::io::stdout().lock();
    while fread(&mut f, &stdin()) == n1 as usize {
        let wpc1 = if n1 == 0 {
            // C compresses zero samples as zero quantization blocks. The
            // library's wpc1Compress rejects an empty slice, but its encoder
            // supports the empty block list and supplies the C metadata.
            let mut wpc1 = Compressed {
                numsmp: n1,
                overallstep: 0.0,
                overallave: 0.0,
                code: Vec::new(),
            };
            wpc1Encoder(&[], &mut wpc1);
            wpc1
        } else {
            wpc1Compress(&mut f, error, &mut nbytes)
        };
        let _ = wpc1Write(&wpc1, &mut output);
    }
    let _ = std::io::Write::flush(&mut output);
    exit(CWP_Exit());
}
