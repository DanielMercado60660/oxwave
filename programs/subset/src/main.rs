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
//! SUBSET - select samples from a three-dimensional file.
//! Port of par/main/resampling/subset.c (Revision 1.12).
//! Credits: Dave Hale, Colorado School of Mines, 07/07/89.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUBSET - select a SUBSET of the samples from a 3-dimensional file\t",
    " \t\t\t\t\t\t\t\t\t",
    " subset <infile >outfile [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " n1=nfloats             number of samples in 1st dimension\t\t",
    " n2=nfloats/n1          number of samples in 2nd dimension\t\t",
    " n3=nfloats/(n1*n2)     number of samples in 3rd dimension\t\t",
    " id1s=1                 increment in samples selected in 1st dimension\t",
    " if1s=0                 index of first sample selected in 1st dimension",
    " n1s=1+(n1-if1s-1)/id1s number of samples selected in 1st dimension\t",
    " ix1s=if1s,if1s+id1s,...indices of samples selected in 1st dimension\t",
    " id2s=1                 increment in samples selected in 2nd dimension\t",
    " if2s=0                 index of first sample selected in 2nd dimension",
    " n2s=1+(n2-if2s-1)/id2s number of samples selected in 2nd dimension\t",
    " ix2s=if2s,if2s+id2s,...indices of samples selected in 2nd dimension\t",
    " id3s=1                 increment in samples selected in 3rd dimension\t",
    " if3s=0                 index of first sample selected in 3rd dimension",
    " n3s=1+(n3-if3s-1)/id3s number of samples selected in 3rd dimension\t",
    " ix3s=if3s,if3s+id3s,...indices of samples selected in 3rd dimension\t",
    " \t\t\t\t\t\t\t\t\t",
    " For the 1st dimension, output is selected from input as follows:\t",
    "       output[i1s] = input[ix1s[i1s]], for i1s = 0 to n1s-1\t\t",
    " Likewise for the 2nd and 3rd dimensions.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn selection(n: i32, axis: &str) -> Vec<i32> {
    let ixname = format!("ix{}s", axis);
    let mut ix = vec![0i32; countparval(&ixname) as usize];
    if getparint(&ixname, &mut ix) == 0 {
        let mut id = 1;
        let mut first = 0;
        let mut count = 0;
        getparint(&format!("id{}s", axis), &mut id);
        getparint(&format!("if{}s", axis), &mut first);
        if getparint(&format!("n{}s", axis), &mut count) == 0 {
            // ARM64 integer division by zero yields zero.
            count = 1i32.wrapping_add(n.wrapping_sub(first).wrapping_sub(1).checked_div(id).unwrap_or(0));
        }
        ix = vec![0; count.max(0) as usize];
        let mut i = first;
        for item in &mut ix {
            *item = i;
            i = i.wrapping_add(id);
        }
    }
    ix
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let infp = stdin();
    let outfp = stdout();
    let mut n1 = 0i32;
    let mut n2 = 0i32;
    let mut n3 = 0i32;
    if getparint("n1", &mut n1) == 0 {
        if efseeko(&infp, 0, SEEK_END) == -1 { err!("input file size unknown; specify n1\n"); }
        n1 = (eftello(&infp) / 4) as i32;
    }
    if getparint("n2", &mut n2) == 0 {
        if efseeko(&infp, 0, SEEK_END) == -1 { err!("input file size unknown; specify n2\n"); }
        // The sizeof operand makes this unsigned; ARM64 division by zero yields zero.
        n2 = (eftello(&infp) as u64).checked_div((n1 as u64).wrapping_mul(4)).unwrap_or(0) as i32;
    }
    if getparint("n3", &mut n3) == 0 {
        if efseeko(&infp, 0, SEEK_END) == -1 { err!("input file size unknown; specify n3\n"); }
        n3 = (eftello(&infp) as u64).checked_div((n2.wrapping_mul(n1) as u64).wrapping_mul(4)).unwrap_or(0) as i32;
    }
    let ix1s = selection(n1, "1");
    let ix2s = selection(n2, "2");
    let ix3s = selection(n3, "3");
    checkpars();
    for (i, &ix) in ix1s.iter().enumerate() {
        if ix < 0 || ix > n1.wrapping_sub(1) { err!("ix1s[%d]=%d is out of bounds!\n", i, ix); }
    }
    for (i, &ix) in ix2s.iter().enumerate() {
        if ix < 0 || ix > n2.wrapping_sub(1) { err!("ix2s[%d]=%d is out of bounds!\n", i, ix); }
    }
    for (i, &ix) in ix3s.iter().enumerate() {
        if ix < 0 || ix > n3.wrapping_sub(1) { err!("ix3s[%d]=%d is out of bounds!\n", i, ix); }
    }
    let mut p = vec![0.0f32; n1.max(0) as usize];
    let mut ps = vec![0.0f32; ix1s.len()];
    for (i3s, &i3) in ix3s.iter().enumerate() {
        for (i2s, &i2) in ix2s.iter().enumerate() {
            let offset = (i2 as i64).wrapping_add((i3 as i64).wrapping_mul(n2 as i64)).wrapping_mul(n1 as i64).wrapping_mul(4);
            efseeko(&infp, offset, SEEK_SET);
            if fread(&mut p, &infp) == n1 as usize {
                for (i1s, &i1) in ix1s.iter().enumerate() { ps[i1s] = p[i1 as usize]; }
            } else {
                err!("no input for ix2s[%d]=%d ix3s[%d]=%d!\n", i2s, i2, i3s, i3);
            }
            efwrite(&ps, &outfp);
        }
    }
    exit(CWP_Exit());
}
