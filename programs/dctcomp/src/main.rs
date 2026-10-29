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

//! DCTCOMP - Compression by Discrete Cosine Transform.
//! Port of comp/dct/main/dctcomp.c (Revision 1.6).
//! Credits (from the C source): CWP: Tong Chen, Dec 1995.

use oxwave_core::prelude::*;
use oxwave_lib_comp_dct::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " DCTCOMP - Compression by Discrete Cosine Transform\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   dctcomp < stdin n1= n2=   [optional parameter] > sdtout\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1=\t\t\tnumber of samples in the fast (first) dimension\t",
    " n2=\t\t\tnumber of samples in the slow (second) dimension",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " blocksize1=16\t\tblocksize in direction 1\t\t\t",
    " blocksize2=16\t\tblocksize in direction 2\t\t\t",
    " error=0.01\t\tacceptable error\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get the parameters */
    let (mut n1, mut n2, mut blocksize1, mut blocksize2) = (0i32, 0i32, 0i32, 0i32);
    let mut error = 0.0f32;
    if getparint("n1", &mut n1) == 0 {
        err!("Must specify n1");
    }
    if getparint("n2", &mut n2) == 0 {
        err!("Must specify n2");
    }
    if getparint("blocksize1", &mut blocksize1) == 0 {
        blocksize1 = 16;
    }
    if getparint("blocksize2", &mut blocksize2) == 0 {
        blocksize2 = 16;
    }
    if getparfloat("error", &mut error) == 0 {
        error = 0.01;
    }
    checkpars();

    /* regular sizes; arm64 C integer division by zero yields zero. */
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

    // Negative padded dimensions/block sizes fault in the C allocator or padding
    // loops. A negative n1 can still have a positive padded dimension.
    if npad1 < 0 || n2 < 0 || blocksize1 < 0 || blocksize2 < 0 {
        unsafe extern "C" {
            fn signal(sig: i32, handler: usize) -> usize;
            fn raise(sig: i32) -> i32;
        }
        unsafe {
            signal(11, 0); // Remove Rust's stack-overflow signal handler.
            raise(11);
        } // SIGSEGV on the target platform.
        exit(139);
    }

    /* allocate spaces */
    // C starts the padding loop at n1, even when n1 is negative. Reserve a
    // prefix for those columns so the successful C path remains addressable.
    // The original may instead fault if its allocation starts at a page edge.
    let prefix = if n1 < 0 { -(n1 as i64) as usize } else { 0 };
    let mut f_storage = vec![0.0f32; prefix + npad1.wrapping_mul(npad2) as usize];
    let mut g = vec![vec![0.0f32; blocksize1 as usize]; blocksize2 as usize];
    let mut qx = vec![0i32; npad1.wrapping_mul(npad2) as usize];

    /* generate the transform tables */
    let c1 = dctAlloc(blocksize1 as usize);
    let c2 = dctAlloc(blocksize2 as usize);

    /* input data */
    for i2 in 0..n2 {
        let row = prefix + (i2 * npad1) as usize;
        if n1 < 0 {
            // Preserve C's signed-to-size_t conversion at this call. The safe
            // core wrapper cannot express a negative count; libc returns a
            // short read on empty stdin, just as it does for the original.
            unsafe extern "C" {
                static __stdinp: *mut std::ffi::c_void;
                #[link_name = "fread"]
                fn c_fread(
                    ptr: *mut std::ffi::c_void,
                    size: usize,
                    count: usize,
                    stream: *mut std::ffi::c_void,
                ) -> usize;
            }
            unsafe {
                c_fread(
                    f_storage.as_mut_ptr().add(row).cast(),
                    std::mem::size_of::<f32>(),
                    n1 as usize,
                    __stdinp,
                );
            }
        } else {
            fread(&mut f_storage[row..row + n1 as usize], &stdin());
        }
        for i1 in n1..npad1 {
            // Negative columns are inside f_storage's reserved prefix. Keep
            // C's loop bounds and writes without indexing before a Rust slice.
            f_storage[(row as i64 + i1 as i64) as usize] = 0.0;
        }
    }
    /* pad with zeroes */
    for i2 in n2..npad2 {
        for i1 in 0..npad1 {
            f_storage[(prefix as i64 + i2 as i64 * npad1 as i64 + i1 as i64) as usize] = 0.0;
        }
    }
    let f = &mut f_storage[prefix..];

    /* DCT for each block */
    let mut ibeg2 = 0;
    for _iblock2 in 0..nblock2 {
        let mut ibeg1 = 0;
        for _iblock1 in 0..nblock1 {
            for j2 in 0..blocksize2 {
                for j1 in 0..blocksize1 {
                    g[j2 as usize][j1 as usize] = f[((ibeg2 + j2) * npad1 + ibeg1 + j1) as usize];
                }
            }
            dct_2(&mut g, &c1, &c2, 0);
            for j2 in 0..blocksize2 {
                for j1 in 0..blocksize1 {
                    f[((ibeg2 + j2) * npad1 + ibeg1 + j1) as usize] = g[j2 as usize][j1 as usize];
                }
            }
            ibeg1 += blocksize1;
        }
        ibeg2 += blocksize2;
    }

    /* quantization */
    let mut step = -1.0f32;
    let mut ave = 0.0f32;
    uniQuant(f, error, &mut ave, &mut step, &mut qx);

    /* prefix encoding */
    let mut ibuff = pEncode(&qx);
    /* allocate out buffer */
    let mut obuff = buffAlloc1(ibuff.code.len());
    /* rewind the in buffer */
    ibuff.pos = 0;
    /* run-length coding */
    codeSilence(&mut ibuff, &mut obuff);
    /* Huffman coding */
    let pos = obuff.pos;
    buffRealloc1(&mut obuff, pos);
    obuff.pos = 0;
    ibuff.pos = 0;
    let nsize = huffCompress(&mut obuff, &mut ibuff);

    fprintf!(stderr(), "size after compression = %d bytes\n", nsize);
    fprintf!(
        stderr(),
        "compression ratio = %f\n",
        ((n1 as f32 * n2 as f32) * std::mem::size_of::<i32>() as f32) / nsize as f32
    );
    fwrite(&[nsize, n1, n2, blocksize1, blocksize2], &stdout());
    fwrite(&[ave, step], &stdout());
    fwrite(&ibuff.code[..ibuff.pos], &stdout());
    exit(CWP_Exit());
}
