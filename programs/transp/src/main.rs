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

//! TRANSP - transpose an n1 by n2 element matrix.
//! Port of par/main/resampling/transp.c (Revision 1.10).
//! Credits: Dave Hale, Colorado School of Mines, 07/07/89.

use oxwave_core::prelude::*;
use oxwave_f_bigmatrix::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " TRANSP - TRANSPose an n1 by n2 element matrix\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " transp <infile >outfile n1= [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1                     number of elements in 1st (fast) dimension of matrix",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " n2=all                 number of elements in 2nd (slow) dimension of matrix",
    " nbpe=sizeof(float)     number of bytes per matrix element\t\t",
    " verbose=0              =1 for diagnostic information\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut n1 = 0;
    let mut n2 = 0;
    let mut nbpe = 0;
    if getparint("n1", &mut n1) == 0 { err!("must specify n1!\n"); }
    if getparint("nbpe", &mut nbpe) == 0 { nbpe = 4; }
    let infp = stdin();
    let outfp = stdout();
    if getparint("n2", &mut n2) == 0 {
        if efseeko(&infp, 0, SEEK_END) == -1 { err!("input file size unknown; specify n2\n"); }
        // On arm64 C integer division by zero yields zero.
        n2 = eftello(&infp).checked_div(n1.wrapping_mul(nbpe) as i64).unwrap_or(0) as i32;
        efseeko(&infp, 0, SEEK_SET);
    }
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    if verbose != 0 { fprintf!(stderr(), "n1=%d  n2=%d  nbpe=%d\n", n1, n2, nbpe); }
    checkpars();
    let row_bytes = n1.wrapping_mul(nbpe);
    let col_bytes = n2.wrapping_mul(nbpe);
    let size = (if n1 > n2 { n1 } else { n2 }).wrapping_mul(nbpe);
    let mut v = vec![0u8; if size > 0 { size as usize } else { 0 }];
    let mut state = bmalloc(nbpe, n1, n2);
    if verbose != 0 { fprintf!(stderr(), "Reading input file\n"); }
    for i2 in 0..n2 {
        if row_bytes < 0 || fread(&mut v[..row_bytes as usize], &infp) != row_bytes as usize {
            err!("Error reading input file:  i2=%d\n", i2);
        }
        if row_bytes > 0 { bmwrite(state.as_mut().unwrap(), 1, 0, i2 as usize, &v[..row_bytes as usize]); }
    }
    if verbose != 0 { fprintf!(stderr(), "Writing output file\n"); }
    for i1 in 0..n1 {
        if col_bytes > 0 { bmread(state.as_ref().unwrap(), 2, i1 as usize, 0, &mut v[..col_bytes as usize]); }
        if col_bytes < 0 || fwrite(&v[..col_bytes as usize], &outfp) != col_bytes as usize {
            err!("Error writing output file:  i1=%d\n", i1);
        }
    }
    if let Some(state) = state { bmfree(state); }
    if verbose != 0 { fprintf!(stderr(), "Transpose done!\n"); }
    exit(CWP_Exit());
}
