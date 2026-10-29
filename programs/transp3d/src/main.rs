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
//! TRANSP3D - transpose a matrix. Port of par/main/resampling/transp3d.c (Revision 1.4).
//! Author: Stewart A. Levin, Landmark Graphics, 01/20/10.

use oxwave_core::prelude::*;
use oxwave_f_vnd::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " TRANSP3D - TRANSPose an n1 by n2 by n3 element matrix\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " transp3d <infile >outfile n1= n2= [optional parameters]\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1\t\tnumber of elements in 1st (fast) dimension of matrix\t",
    " n2\t\tnumber of elements in 2nd (middle) dimension of matrix\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " n3=all\tnumber of elements in 3rd (slow) dimension of matrix\t",
    " perm=231\t\tdesired output axis ordering\t\t   \t",
    " nbpe=sizeof(float)\tnumber of bytes per matrix element\t\t",
    " scratchdir=/tmp\tdirectory for scratch files\t\t\t",
    " scratchstem=foo\tstem prefix for scratch file names\t\t",
    " verbose=0\t\t=1 for diagnostic information\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut n1, mut n2, mut n3) = (0, 0, 0);
    if getparint("n1", &mut n1) == 0 { err!("must specify n1!"); }
    if getparint("n2", &mut n2) == 0 { err!("must specify n2!"); }
    let ntempdirs = countparval("scratchdir");
    let mut tempdirs = Vec::new();
    if ntempdirs > 0 { getparstringarray("scratchdir", &mut tempdirs); }
    else { tempdirs.push("/tmp".to_string()); }
    let mut tempstem = String::new();
    if getparstring("scratchstem", &mut tempstem) == 0 { tempstem = "foo".to_string(); }
    let (mut nbpe, mut perm) = (4, 231);
    getparint("nbpe", &mut nbpe);
    getparint("perm", &mut perm);
    let infp = stdin();
    let outfp = stdout();
    if getparint("n3", &mut n3) == 0 {
        if efseeko(&infp, 0, 2) == -1 { err!("input file size unknown; specify n3\n"); }
        // arm64 C integer division by zero yields zero.
        n3 = eftello(&infp).checked_div((n1 as i64).wrapping_mul(n2 as i64).wrapping_mul(nbpe as i64)).unwrap_or(0) as i32;
        efseeko(&infp, 0, 0);
    }
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    if verbose != 0 { warn!("n1=%d  n2=%d n3=%d nbpe=%d perm=%d\n", n1, n2, n3, nbpe, perm); }
    checkpars();
    let mut v = vec![0u8; (n1.wrapping_add(n2).wrapping_add(n3) as usize).wrapping_mul(nbpe as usize)];
    let n = [n1 as i64, n2 as i64, n3 as i64];
    let dirs: Vec<&str> = tempdirs.iter().map(String::as_str).collect();
    let mut handle = VNDop(2, 200000000, &n, 1, nbpe as i64, &tempstem, dirs.len() as i32, &dirs, 1);
    if verbose != 0 { warn!("Reading input file"); }
    for i3 in 0..n3 { for i2 in 0..n2 {
        let count = nbpe.wrapping_mul(n1) as usize;
        if fread(&mut v[..count], &infp) != count { err!("Error reading input file:  i2=%d\n", i2); }
        VNDrw(b'w', 0, handle.as_mut().unwrap(), 0, &[0, i2 as i64, i3 as i64], 0, &mut v, 0, 1, n1 as i64, 21, "");
    }}
    if verbose != 0 { warn!("Writing output file"); }
    match perm {
        123 | 132 => {
            let (outer, inner) = if perm == 123 { (n3, n2) } else { (n2, n3) };
            for a in 0..outer { for b in 0..inner {
                let (i2, i3) = if perm == 123 { (b, a) } else { (a, b) };
                VNDrw(b'r', 0, handle.as_mut().unwrap(), 0, &[0, i2 as i64, i3 as i64], 0, &mut v, 0, 1, n1 as i64, perm, "");
                let count = nbpe.wrapping_mul(n1) as usize;
                if fwrite(&v[..count], &outfp) != count { err!("Error writing output file:  i2=%d i3=%d", i2, i3); }
            }}
        }
        213 | 231 => {
            let (outer, inner) = if perm == 213 { (n3, n1) } else { (n1, n3) };
            for a in 0..outer { for b in 0..inner {
                let (i1, i3) = if perm == 213 { (b, a) } else { (a, b) };
                VNDrw(b'r', 0, handle.as_mut().unwrap(), 1, &[i1 as i64, 0, i3 as i64], 0, &mut v, 0, 1, n2 as i64, perm, "");
                let count = nbpe.wrapping_mul(n2) as usize;
                if fwrite(&v[..count], &outfp) != count { err!("Error writing output file:  i1=%d i3=%d", i1, i3); }
            }}
        }
        312 | 321 => {
            let (outer, inner) = if perm == 312 { (n2, n1) } else { (n1, n2) };
            for a in 0..outer { for b in 0..inner {
                let (i1, i2) = if perm == 312 { (b, a) } else { (a, b) };
                VNDrw(b'r', 0, handle.as_mut().unwrap(), 2, &[i1 as i64, i2 as i64, 0], 0, &mut v, 0, 1, n3 as i64, perm, "");
                let count = nbpe.wrapping_mul(n3) as usize;
                if fwrite(&v[..count], &outfp) != count { err!("Error writing output file:  i1=%d i2=%d", i1, i2); }
            }}
        }
        _ => err!("Unrecognized transpose permutation %d", perm),
    }
    if let Some(handle) = handle { VNDcl(handle, 1); }
    if verbose != 0 { warn!("Transpose done!"); }
    exit(CWP_Exit());
}
