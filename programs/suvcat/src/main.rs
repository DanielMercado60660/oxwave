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

//! SUVCAT - append corresponding traces from two data sets.
//! Port of su/main/operations/suvcat.c (Revision 1.18).
//! Credits: CWP: Jack K. Cohen, Michel Dietrich, Steven D. Sheaffer;
//! IfG Kiel: Thies Beilecke.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUVCAT -  append one data set to another, with or without an  ",
    "           overlapping\tregion.  Data in the overlap may be     ",
    "           determined by one of several methods.               ",
    " \t\t\t\t\t\t\t\t",
    " suvcat data1 data2 >stdout\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "        none\t\t\t\t\t\t\t",
    "                                                               ",
    " Optional parameters for overlapping sections:\t\t\t",
    "                                                               ",
    "  taplen=0    Length of overlap in integer number of           ",
    "                  samples.(Default is 0.)                      ",
    "                                                               ",
    "  taptype=0    Type of taper or combination method in the\t",
    "                  overlap region.  0 - average                 ",
    "                                   1 - maximum magnitude       ",
    "                                   2 - cosine scaled           ",
    "                                   3 - summation               ",
    "                                                               ",
    " Computational Notes:\t\t\t\t\t\t",
    "                                                               ",
    " This program vertically concatenates traces from data2 onto   ",
    " the end of the corresponding traces in data1, with a region   ",
    " of overlap, defined by taplen.  Data in the overlapping       ",
    " region is combined by the method specified by taptype. The    ",
    " currently available methods are:                              ",
    "                                                               ",
    "     taptype=0    output is assigned the unweighted average of ",
    "                  each point in the overlap                    ",
    "     taptype=1    output is assigned the value of the maximum  ",
    "                  absolute value of each point in the overlap  ",
    "     taptype=2    output is assigned the weighted average of   ",
    "                  each point in the overlap, where the output  ",
    "                  is the sum of cos(x) times the values on the ",
    "                  first section, and 1-cos(x) times the values ",
    "                  on the second section, where x is factor that",
    "                  goes from 0 to pi/2 across the overlap. This ",
    "                  favors the upper section in the upper part of",
    "                  the overlap, and favors the lower section in ",
    "                  the lower part of the overlap.               ",
    "     taptype=3    output is assigned the sum of the amplitudes ",
    "                  at each sample in the overlap                ",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut intrace1 = Trace::new();
    let mut intrace2 = Trace::new();
    let mut data2flag = 0;
    let mut itr = 0;
    let mut taplen = 0;
    let mut taptype = 0;
    initargs(&argv);
    requestdoc(2, SDOC);
    getparint("taplen", &mut taplen);
    getparint("taptype", &mut taptype);
    checkpars();
    if taptype > 3 { err!("taptype must be 0, 1, 2, or 3!\n"); }
    let fp1 = efopen(&argv[1], "r");
    let fp2 = efopen(&argv[2], "r");
    while fgettr(&fp1, &mut intrace1) != 0 && {
        data2flag = fgettr(&fp2, &mut intrace2);
        data2flag != 0
    } {
        let nt1 = intrace1.ns as i32;
        let nt2 = intrace2.ns as i32;
        // Preserve C integer arithmetic even for extreme invalid overlaps.
        let nt = (nt1 + nt2).wrapping_sub(taplen);
        if nt > SU_NFLTS as i32 {
            err!("nt=%d exceeds SU_NFLTS=%d", nt, SU_NFLTS as i32);
        }
        // Preserve the sample-limit error above before rejecting undefined C
        // accesses: taplen > nt1 starts before intrace1.data, and taplen > nt2
        // writes past outtrace. Negative overlaps also index before intrace2.
        if taplen < 0 || taplen > nt1 || taplen > nt2 {
            err!("taplen=%d must be between 0 and both input trace lengths (%d, %d)", taplen, nt1, nt2);
        }
        // A negative taper type leaves the C overlap buffer uninitialized.
        if taplen > 0 && taptype < 0 {
            err!("taptype must be 0, 1, 2, or 3!\n");
        }
        let mut outtrace = vec![0.0f32; nt as usize];
        for i in 0..nt1-taplen {
            outtrace[i as usize] = intrace1.data[i as usize];
        }
        if taptype == 0 {
            for i in 0..taplen {
                let hold = intrace1.data[(nt1-taplen+i) as usize] + intrace2.data[i as usize];
                outtrace[(nt1-taplen+i) as usize] = (hold as f64 * 0.5) as f32;
            }
        }
        if taptype == 1 {
            for i in 0..taplen {
                let a = intrace1.data[(nt1-taplen+i) as usize];
                let b = intrace2.data[i as usize];
                let hold = if a.abs() > b.abs() { a } else { b };
                outtrace[(nt1-taplen+i) as usize] = hold;
            }
        }
        if taptype == 2 {
            for i in 0..taplen {
                // C integer division on arm64 returns zero for a zero divisor.
                let ratio = i.checked_div(taplen-1).unwrap_or(0);
                let s1 = (std::f64::consts::PI / 2.0 * ratio as f64).cos() as f32;
                let s2 = 1.0f32 - s1;
                let hold = s1 * intrace1.data[(nt1-taplen+i) as usize] + s2 * intrace2.data[i as usize];
                outtrace[(nt1-taplen+i) as usize] = hold;
            }
        }
        if taptype == 3 {
            for i in 0..taplen {
                let hold = intrace1.data[(nt1-taplen+i) as usize] + intrace2.data[i as usize];
                outtrace[(nt1-taplen+i) as usize] = hold;
            }
        }
        for i in 0..nt2-taplen {
            outtrace[(nt1+i) as usize] = intrace2.data[(i+taplen) as usize];
        }
        intrace1.ns = nt as u16;
        for i in 0..nt {
            intrace1.data[i as usize] = outtrace[i as usize];
        }
        puttr(&mut intrace1);
        itr += 1;
    }
    if data2flag == 0 {
        warn!("%s still had traces when %s was exhausted", argv[1], argv[2]);
        warn!("processed %d pairs of traces before EOF", itr);
    } else if fgettr(&fp2, &mut intrace2) != 0 {
        warn!("%s still had traces when %s was exhausted", argv[2], argv[1]);
        warn!("processed %d pairs of traces before EOF", itr);
    }
    exit(CWP_Exit());
}
