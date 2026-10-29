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

//! WTCOMP - Compression by Wavelet Transform.
//! Port of comp/dct/main/wtcomp.c (Revision 1.5).
//! Credits: CWP: Tong Chen, Dec 1995.

use oxwave_core::prelude::*;
use oxwave_lib_comp_dct::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " WTCOMP - Compression by Wavelet Transform\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   wtcomp < stdin n1= n2=   [optional parameter] > sdtout\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1=\t\t\tnumber of samples in the fast (first) dimension\t",
    " n2=\t\t\tnumber of samples in the slow (second) dimension",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nstage1=\t\tnumber of stages (set automatically based on n1)",
    " nstage2=\t\tnumber of stages (set automatically based on n2)",
    " nfilter=11\t\tnumber of filters\t\t\t\t",
    " error=0.01\t\tacceptable error\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn invalid_stage(nstage1: i32, nstage2: i32) -> Option<&'static str> {
    if nstage1 < 0 {
        Some("nstage1 must be nonnegative")
    } else if nstage2 < 0 {
        Some("nstage2 must be nonnegative")
    } else {
        None
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut n1, mut n2, mut nfilter) = (0, 0, 11);
    let mut error = 0.01f32;
    if getparint("n1", &mut n1) == 0 { err!("Must specify n1\n"); }
    if getparint("n2", &mut n2) == 0 { err!("Must specify n2\n"); }
    getparint("nfilter", &mut nfilter);
    getparfloat("error", &mut error);
    let mut npow1 = 0i32;
    if n1 != 1 { while (n1 - 1) >> npow1 != 0 { npow1 += 1; } }
    let npad1 = 1usize << npow1;
    let mut npow2 = 0i32;
    if n2 != 1 { while (n2 - 1) >> npow2 != 0 { npow2 += 1; } }
    let npad2 = 1usize << npow2;
    let (mut nstage1, mut nstage2) = (npow1, npow2);
    getparint("nstage1", &mut nstage1);
    getparint("nstage2", &mut nstage2);
    if let Some(message) = invalid_stage(nstage1, nstage2) {
        err!("%s", message);
    }
    checkpars();
    let mut f = vec![vec![0.0f32; npad1]; npad2];
    let mut g = vec![vec![0.0f32; npad1]; npad2];
    let mut qx = vec![0i32; npad1 * npad2];
    let filter = waveGetfilter(nfilter).unwrap();
    for row in f.iter_mut().take(n2 as usize) {
        fread(&mut row[..n1 as usize], &stdin());
    }
    waveTrans_2(&mut f, &mut g, &filter, npow1 as usize, npow2 as usize,
                nstage1 as usize, nstage2 as usize, 0);
    let flat: Vec<f32> = g.into_iter().flatten().collect();
    let (mut ave, mut step) = (0.0f32, -1.0f32);
    uniQuant(&flat, error, &mut ave, &mut step, &mut qx);
    let mut ibuff = pEncode(&qx);
    let mut obuff = buffAlloc1(ibuff.code.len());
    ibuff.pos = 0;
    codeSilence(&mut ibuff, &mut obuff);
    let pos = obuff.pos;
    buffRealloc1(&mut obuff, pos);
    obuff.pos = 0;
    ibuff.pos = 0;
    let nsize = huffCompress(&mut obuff, &mut ibuff);
    fprintf!(stderr(), "size after compression = %d bytes\n", nsize);
    fprintf!(stderr(), "compression ratio = %f\n", (n1 as f32 * n2 as f32 * 4.0f32 / nsize as f32));
    fwrite(&[nsize, n1, n2, nfilter, nstage1, nstage2], &stdout());
    fwrite(&[ave, step], &stdout());
    fwrite(&ibuff.code[..ibuff.pos], &stdout());
    exit(CWP_Exit());
}

#[cfg(test)]
mod tests {
    use super::invalid_stage;

    #[test]
    fn negative_stages_are_rejected() {
        assert_eq!(invalid_stage(-1, 0), Some("nstage1 must be nonnegative"));
        assert_eq!(invalid_stage(0, -1), Some("nstage2 must be nonnegative"));
        assert_eq!(invalid_stage(0, 0), None);
        assert_eq!(invalid_stage(3, 2), None);
    }
}
