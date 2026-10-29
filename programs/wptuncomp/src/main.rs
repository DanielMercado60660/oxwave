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

//! WPTUNCOMP - Uncompress WPT compressed data. Port of comp/dct/main/wptuncomp.c (Revision 1.5).
//! Author: CWP: Tong Chen, Dec 1995.

use oxwave_core::prelude::*;
use oxwave_lib_comp_dct::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " WPTUNCOMP - Uncompress  WPT compressed data\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   wptuncomp < stdin > sdtout\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    "  none\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    "  none\t\t\t\t\t\t\t\t\t",
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
    let (mut nsize, mut n1, mut n2) = (0i32, 0i32, 0i32);
    let (mut nfilter1, mut nfilter2, mut nstage1, mut nstage2) = (0i32, 0i32, 0i32, 0i32);
    let (mut ave, mut step) = (0f32, 0f32);
    let input = stdin();
    fread(std::slice::from_mut(&mut nsize), &input);
    fread(std::slice::from_mut(&mut n1), &input);
    fread(std::slice::from_mut(&mut n2), &input);
    fread(std::slice::from_mut(&mut nfilter1), &input);
    fread(std::slice::from_mut(&mut nfilter2), &input);
    fread(std::slice::from_mut(&mut nstage1), &input);
    fread(std::slice::from_mut(&mut nstage2), &input);
    fread(std::slice::from_mut(&mut ave), &input);
    fread(std::slice::from_mut(&mut step), &input);

    if n1 <= 0 || n2 <= 0 {
        segfault();
    }
    let npow1 = if n1 == 1 { 0 } else { (32 - (n1 - 1).leading_zeros()) as usize };
    let npow2 = if n2 == 1 { 0 } else { (32 - (n2 - 1).leading_zeros()) as usize };
    let npad1 = 1usize << npow1;
    let npad2 = 1usize << npow2;
    let mut f = vec![vec![0f32; npad1]; npad2];
    let mut g = vec![vec![0f32; npad1]; npad2];
    let mut qx = vec![0i32; npad1 * npad2];
    let filter1 = waveGetfilter(nfilter1);
    let filter2 = waveGetfilter(nfilter2);
    let nmax = 2 * npad1 * npad2;
    // The compressor can write -1 when its bounded output buffer fills.
    // The C decoder then dereferences a failed allocation.
    if nsize < 0 {
        segfault();
    }
    let mut ibuff = buffAlloc1(nsize as usize);
    let mut obuff = buffAlloc1(nmax);
    fread(&mut ibuff.code, &input);
    ibuff.pos = 0;
    if huffDecompress(&mut ibuff, &mut obuff) == -1 {
        err!("Inconsistent data \n");
    }
    let size = obuff.pos;
    buffRealloc1(&mut obuff, size);
    buffRealloc1(&mut ibuff, nmax);
    ibuff.pos = 0;
    obuff.pos = 0;
    codeDesilence(&mut obuff, &mut ibuff);
    let size = ibuff.pos;
    buffRealloc1(&mut ibuff, size);
    ibuff.pos = 0;
    if pDecode(&mut ibuff, &mut qx) == -1 {
        err!("Inconsistent data \n");
    }
    let mut values = vec![0f32; npad1 * npad2];
    uniDequant(&mut values, ave, step, &qx);
    for i2 in 0..npad2 {
        g[i2].copy_from_slice(&values[i2 * npad1..(i2 + 1) * npad1]);
    }
    let (Some(filter1), Some(filter2)) = (filter1, filter2) else {
        segfault();
    };
    wavePack_2(&mut f, &mut g, &filter1, &filter2, npow1, npow2, nstage1 as usize, nstage2 as usize, 1);
    for i2 in 0..n2 as usize {
        fwrite(&f[i2][..n1 as usize], &stdout());
    }
    exit(CWP_Exit());
}
