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

//! Z2XYZ - convert binary floats to ASCII triples.
//! Port of par/main/resampling/z2xyz.c (Revision 1.11).
//! Credits: CWP: John Stockwell, based on "b2a" by Jack Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " Z2XYZ - convert binary floats representing Z-values to ascii\t",
    " \t   form in X Y Z ordered triples\t\t\t",
    " \t\t\t\t\t\t\t\t",
    "    z2xyz <stdin >stdout \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " n1=\t\tnumber of floats in 1st (fast) dimension\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " outpar=/dev/tty\t\t output parameter file\t\t",
    " \t\t\t\t\t\t\t\t",
    " Notes: This program is useful for converting panels of float\t",
    " data (representing evenly spaced z values) to the x y z\t",
    " ordered triples required for certain 3D plotting packages.\t",
    " \t\t\t\t\t\t\t\t",
    " Example of NXplot3d usage on a NeXT:\t\t\t\t",
    " suplane | sufilter | z2xyz n1=64 > junk.ascii\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Now open junk.ascii as a mesh data file with NXplot3d.\t",
    " (NXplot3d is a NeXTStep-only utility for viewing 3d data sets\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".to_string(); }
    let outparfp = efopen(&outpar, "w");
    let mut n1 = 0i32;
    mustgetparint("n1", &mut n1);
    checkpars();
    // A negative C int converts to a huge size_t allocation, which fails.
    if n1 < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut z = vec![0.0f32; n1 as usize];
    let mut n2 = 0usize;
    loop {
        let n1read = efread(&mut z, &stdin());
        if n1read == 0 { break; }
        if n1read != n1 as usize { err!("out of data in forming line #%d", n2 + 1); }
        for i1 in 0..n1 { printf!("%lu %d %11.4e \n", n2, i1, z[i1 as usize]); }
        n2 += 1;
    }
    fprintf!(outparfp, "n2=%lu\n", n2);
    exit(CWP_Exit());
}
