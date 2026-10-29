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

//! WTUNCOMP - uncompress WT compressed data.
//! Port of comp/dct/main/wtuncomp.c (Revision 1.5).
//! Credits: CWP: Tong Chen, 1995.

use oxwave_core::prelude::*;
use oxwave_lib_comp_dct::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " WTUNCOMP - UNCOMPression of WT compressed data\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   wtuncomp < stdin > sdtout\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    "    none\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    "    none\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn segfault() -> ! {
    unsafe extern "C" {
        fn signal(sig: i32, handler: usize) -> usize;
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        signal(11, 0);
        raise(11);
    }
    exit(139);
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut header = [0i32; 6];
    for value in &mut header {
        fread(std::slice::from_mut(value), &stdin());
    }
    let [nsize, n1, n2, nfilter, nstage1, nstage2] = header;
    let mut ave = 0.0f32;
    let mut step = 0.0f32;
    fread(std::slice::from_mut(&mut ave), &stdin());
    fread(std::slice::from_mut(&mut step), &stdin());
    if n1 <= 0 || n2 <= 0 {
        segfault();
    }
    let npow1 = if n1 == 1 { 0 } else { (n1 - 1).ilog2() + 1 };
    let npow2 = if n2 == 1 { 0 } else { (n2 - 1).ilog2() + 1 };
    let npad1 = 1usize << npow1;
    let npad2 = 1usize << npow2;
    let mut f = vec![vec![0.0f32; npad1]; npad2];
    let mut g = vec![vec![0.0f32; npad1]; npad2];
    let mut qx = vec![0i32; npad1 * npad2];
    let nmax = 2 * npad1 * npad2;
    if nsize < 0 {
        // The compressor can write -1 when its Huffman buffer fills; C
        // subsequently dereferences the failed allocation.
        segfault();
    }
    let mut ibuff = buffAlloc1(nsize as usize);
    let mut obuff = buffAlloc1(nmax);
    let filter = waveGetfilter(nfilter);
    fread(&mut ibuff.code, &stdin());
    ibuff.pos = 0;
    if huffDecompress(&mut ibuff, &mut obuff) == -1 {
        err!("Inconsistent data \n");
    }
    let pos = obuff.pos;
    buffRealloc1(&mut obuff, pos);
    buffRealloc1(&mut ibuff, nmax);
    obuff.pos = 0;
    ibuff.pos = 0;
    codeDesilence(&mut obuff, &mut ibuff);
    let pos = ibuff.pos;
    buffRealloc1(&mut ibuff, pos);
    ibuff.pos = 0;
    if pDecode(&mut ibuff, &mut qx) == -1 {
        err!("Inconsistent data \n");
    }
    let mut values = vec![0.0f32; npad1 * npad2];
    uniDequant(&mut values, ave, step, &qx);
    for i2 in 0..npad2 {
        g[i2].copy_from_slice(&values[i2 * npad1..(i2 + 1) * npad1]);
    }
    let filter = filter.unwrap_or_else(|| segfault());
    waveTrans_2(&mut f, &mut g, &filter, npow1 as usize, npow2 as usize,
                nstage1 as usize, nstage2 as usize, 1);
    for i2 in 0..n2 as usize {
        fwrite(&f[i2][..n1 as usize], &stdout());
    }
    exit(CWP_Exit());
}
