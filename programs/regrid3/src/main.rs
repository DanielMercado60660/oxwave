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
//! REGRID3 - rewrite a three-dimensional grid.
//! Port of par/main/resampling/regrid3.c (Revision 1.4).
//! Credits: CWP: Zhaobo Meng, 1996, Colorado School of Mines.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " REGRID3 - REwrite a [ni3][ni2][ni1] GRID to a [no3][no2][no1] 3-D grid",
    " \t\t\t\t\t\t\t\t\t",
    " regrid3 < oldgrid > newgrid [parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " ni1=1   fastest (3rd) dimension in input grid                         ",
    " ni2=1   second fastest (2nd) dimension in input grid\t\t\t",
    " ni3=1   slowest (1st) dimension in input grid                       \t",
    " no1=1   fastest (3rd) dimension in output grid                        ",
    " no2=1   second fastest (2nd) dimension in output grid                 ",
    " no3=1   slowest (1st) dimension in output grid                        ",
    " Optional Parameters:                                                  ",
    " verbose=0\t=1 print some useful information\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " REGRID3 can be used to span a 1-D grid to a 2-D grid,  or a 2-D grid  ",
    " to a 3-D grid; or to change grid parameters within the dimensions.\t",
    " Together with MUL and UNIF3, most 3-D velocity model can be \t\t",
    " constructed.\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut ni1, mut ni2, mut ni3) = (1i32, 1i32, 1i32);
    let (mut no1, mut no2, mut no3) = (1i32, 1i32, 1i32);
    let mut verbose = 0i32;
    getparint("ni1", &mut ni1);
    getparint("ni2", &mut ni2);
    getparint("ni3", &mut ni3);
    getparint("no1", &mut no1);
    getparint("no2", &mut no2);
    getparint("no3", &mut no3);
    getparint("verbose", &mut verbose);
    checkpars();
    if verbose != 0 {
        warn!("ni1,ni2,ni3=%d %d %d", ni1, ni2, ni3);
        warn!("no3,no2,no1=%d %d %d", no1, no2, no3);
    }
    if ni1 <= 0 || ni2 <= 0 || ni3 <= 0 || no1 < 0 || no2 <= 0 || no3 <= 0 {
        err!("ni1,ni2,ni3 and no1,no2,no3 must be all positive");
    }
    let mut gridi = vec![0.0f32; ni1 as usize * ni2 as usize * ni3 as usize];
    let mut grido = vec![0.0f32; no1 as usize * no2 as usize * no3 as usize];
    if fread(&mut gridi, &stdin()) != gridi.len() {
        err!("Can not read in gridi");
    }
    let index = |i3: i32, i2: i32, i1: i32| -> usize {
        (i3 as usize * ni2 as usize + i2 as usize) * ni1 as usize + i1 as usize
    };
    for io1 in 0..no1 {
        let xi1 = io1 as f32 * ni1 as f32 / no1 as f32;
        let ix1 = (ni1 - 1).min(xi1 as i32);
        let alpha1 = xi1 - ix1 as f32;
        let ix1plus = (ni1 - 1).min(ix1 + 1);
        for io2 in 0..no2 {
            let xi2 = io2 as f32 * ni2 as f32 / no2 as f32;
            let ix2 = (ni2 - 1).min(xi2 as i32);
            let alpha2 = xi2 - ix2 as f32;
            let ix2plus = (ni2 - 1).min(ix2 + 1);
            for io3 in 0..no3 {
                let xi3 = io3 as f32 * ni3 as f32 / no3 as f32;
                let ix3 = (ni3 - 1).min(xi3 as i32);
                let alpha3 = xi3 - ix3 as f32;
                let ix3plus = (ni3 - 1).min(ix3 + 1);
                grido[(io3 as usize * no2 as usize + io2 as usize) * no1 as usize + io1 as usize] =
                    alpha1 * (alpha2 *
                        (alpha3 * gridi[index(ix3plus, ix2plus, ix1plus)] +
                         (1.0 - alpha3) * gridi[index(ix3plus, ix2plus, ix1)]) +
                        (1.0 - alpha2) *
                        (alpha3 * gridi[index(ix3plus, ix2, ix1plus)] +
                         (1.0 - alpha3) * gridi[index(ix3plus, ix2, ix1)])) +
                    (1.0 - alpha1) * (alpha2 *
                        (alpha3 * gridi[index(ix3, ix2plus, ix1plus)] +
                         (1.0 - alpha3) * gridi[index(ix3, ix2plus, ix1)]) +
                        (1.0 - alpha2) *
                        (alpha3 * gridi[index(ix3, ix2, ix1plus)] +
                         (1.0 - alpha3) * gridi[index(ix3, ix2, ix1)]));
            }
        }
    }
    fwrite(&grido, &stdout());
    exit(CWP_Exit());
}
