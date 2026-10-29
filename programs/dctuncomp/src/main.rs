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

//! DCTUNCOMP - Discrete Cosine Transform Uncompression.
//! Port of comp/dct/main/dctuncomp.c (Revision 1.5).
//! Credits (from the C source): Author: CWP: Tong Chen Dec 1995.

use oxwave_core::prelude::*;
use oxwave_lib_comp_dct::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " DCTUNCOMP - Discrete Cosine Transform Uncompression \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   dctuncomp < stdin >stdout \t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The input of this program is a file compressed by dctcomp.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    // The C reads are unchecked. Zero initialization also preserves the
    // observed empty/short-header path to "Inconsistent data".
    let mut nsize = 0i32;
    let mut n1 = 0i32;
    let mut n2 = 0i32;
    let mut blocksize1 = 0i32;
    let mut blocksize2 = 0i32;
    let mut ave = 0.0f32;
    let mut step = 0.0f32;

    initargs(&argv);
    requestdoc(1, SDOC);

    /* get the parameters */
    fread(std::slice::from_mut(&mut nsize), &stdin());
    fread(std::slice::from_mut(&mut n1), &stdin());
    fread(std::slice::from_mut(&mut n2), &stdin());
    fread(std::slice::from_mut(&mut blocksize1), &stdin());
    fread(std::slice::from_mut(&mut blocksize2), &stdin());
    fread(std::slice::from_mut(&mut ave), &stdin());
    fread(std::slice::from_mut(&mut step), &stdin());

    /* regular sizes */
    // ARM64 C integer division by zero produces zero.
    let nblock1 = n1
        .wrapping_sub(1)
        .checked_div(blocksize1)
        .unwrap_or(0)
        .wrapping_add(1);
    let nblock2 = n2
        .wrapping_sub(1)
        .checked_div(blocksize2)
        .unwrap_or(0)
        .wrapping_add(1);
    let npad1 = nblock1.wrapping_mul(blocksize1);
    let npad2 = nblock2.wrapping_mul(blocksize2);

    /* allocate spaces */
    // f[0] is passed to uniDequant as a contiguous allocation in C.
    let mut f = vec![0.0f32; npad1.wrapping_mul(npad2) as usize];
    let mut g = vec![vec![0.0f32; blocksize1 as usize]; blocksize2 as usize];
    let mut qx = vec![0i32; npad1.wrapping_mul(npad2) as usize];

    /* allocate buffers */
    let nmax = 2i32.wrapping_mul(npad1).wrapping_mul(npad2);
    if nsize < 0 {
        // dctcomp can write nsize=-1 when its bounded Huffman buffer fills.
        // C converts that to a huge allocation, then dereferences NULL in
        // the unchecked buffer read. Preserve its SIGSEGV instead of a
        // Rust capacity-overflow panic (which has different stderr/status).
        unsafe extern "C" {
            fn raise(sig: i32) -> i32;
        }
        unsafe { raise(11) };
        exit(139);
    }
    let mut ibuff = buffAlloc1(nsize as usize);
    let mut obuff = buffAlloc1(nmax as usize);

    /* read data */
    fread(&mut ibuff.code, &stdin());

    /* Huffman decoding */
    if huffDecompress(&mut ibuff, &mut obuff) == -1 {
        err!("Inconsistent data \n");
    }

    /* run-length decoding */
    let pos = obuff.pos;
    buffRealloc1(&mut obuff, pos);
    obuff.pos = 0; // buffRewind
    buffRealloc1(&mut ibuff, nmax as usize);
    ibuff.pos = 0; // buffRewind
    codeDesilence(&mut obuff, &mut ibuff);

    /* prefix decoding */
    let pos = ibuff.pos;
    buffRealloc1(&mut ibuff, pos);
    ibuff.pos = 0; // buffRewind
    if pDecode(&mut ibuff, &mut qx) == -1 {
        err!("Inconsistent data \n");
    }

    /* dequantization */
    uniDequant(&mut f, ave, step, &qx);

    /* generate the transform tables */
    let c1 = dctAlloc(blocksize1 as usize);
    let c2 = dctAlloc(blocksize2 as usize);

    /* inverse DCT for each block */
    let mut ibeg2 = 0i32;
    for _iblock2 in 0..nblock2 {
        let mut ibeg1 = 0i32;
        for _iblock1 in 0..nblock1 {
            for j2 in 0..blocksize2 {
                let i2 = ibeg2 + j2;
                for j1 in 0..blocksize1 {
                    let i1 = ibeg1 + j1;
                    g[j2 as usize][j1 as usize] = f[(i2 * npad1 + i1) as usize];
                }
            }

            dct_2(&mut g, &c1, &c2, 1);

            for j2 in 0..blocksize2 {
                let i2 = ibeg2 + j2;
                for j1 in 0..blocksize1 {
                    let i1 = ibeg1 + j1;
                    f[(i2 * npad1 + i1) as usize] = g[j2 as usize][j1 as usize];
                }
            }
            ibeg1 += blocksize1;
        }
        ibeg2 += blocksize2;
    }

    for i2 in 0..n2 {
        let start = (i2 * npad1) as usize;
        fwrite(&f[start..start + n1 as usize], &stdout());
    }

    exit(CWP_Exit());
}
