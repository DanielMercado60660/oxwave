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

// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
// David Forel, May 2005.
//! TVNMOQC - Check tnmo-vnmo pairs; create t-v column files.
//! Port of par/main/parameter_file_utilities/tvnmoqc.c.
//! Credits: MTU: David Forel (adapted from SUNMO).

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                    ",
    " TVNMOQC - Check tnmo-vnmo pairs; create t-v column files           ",
    "                                                                    ",
    " tvnmoqc [parameters] cdp=... tnmo=... vnmo=...                     ",
    "                                                                    ",
    "   Example:                                                         ",
    " tvnmoqc mode=1 \\                                                  ",
    " cdp=15,35 \\                                                       ",
    " tnmo=0.0091,0.2501,0.5001,0.7501,0.9941 \\                         ",
    " vnmo=1497.0,2000.0,2500.0,3000.0,3500.0 \\                         ",
    " tnmo=0.0082,0.2402,0.4902,0.7402,0.9842 \\                         ",
    " vnmo=1495.0,1900.0,2400.0,2900.0,3400.0                            ",
    "                                                                    ",
    " Required Parameter:                                                ",
    "   prefix=        Prefix of output t-v file(s)                      ",
    "                  Required only for mode=2                          ",
    "                                                                    ",
    " Optional Parameter:                                                ",
    "   mode=1         1=qc: check that tnmo values increase             ",
    "                  2=qc and output t-v files                         ",
    "                                                                    ",
    " mode=1                                                             ",
    "   TVNMOQC checks that there is a tnmo and vnmo series for each CDP ",
    "     and it checks that each tnmo series increases in time.         ",
    "                                                                    ",
    " mode=2                                                             ",
    "   TVNMOQC does mode=1 checking, plus ...                           ",
    "                                                                    ",
    "   TVNMOQC converts par (MKPARFILE) values written as:              ",
    "                                                                    ",
    "          cdp=15,35,...,95 \\                                       ",
    "          tnmo=t151,t152,...,t15n \\                                ",
    "          vnmo=v151,v152,...,v15n \\                                ",
    "          tnmo=t351,t352,...,t35n \\                                ",
    "          vnmo=v351,v352,...,v35n \\                                ",
    "          tnmo=... \\                                               ",
    "          vnmo=... \\                                               ",
    "          tnmo=t951,t952,...,t95n \\                                ",
    "          vnmo=v951,v952,...,v95n \\                                ",
    "                                                                    ",
    "   to column format. The format of each output file is:             ",
    "                                                                    ",
    "          t1 v1                                                     ",
    "          t2 v2                                                     ",
    "           ...                                                      ",
    "          tn vn                                                     ",
    "                                                                    ",
    "   One file is output for each input pair of tnmo-vnmo series.      ",
    "                                                                    ",
    "   A CDP VALUE MUST BE SUPPLIED FOR EACH TNMO-VNMO ROW PAIR.        ",
    "                                                                    ",
    "   Prefix of each output file is the user-supplied value of         ",
    "     parameter PREFIX.                                              ",
    "   Suffix of each output file is the cdp value.                     ",
    "   For the example above, output files names are:                   ",
    "     PREFIX.15  PREFIX.35  ...  PREFIX.95                           ",
    "                                                                    ",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut mode = 1;
    if getparint("mode", &mut mode) == 0 { mode = 1; }
    let mut prefix = String::new();
    if mode == 2 && getparstring("prefix", &mut prefix) == 0 {
        err!("When mode=2, you must supply a prefix name.");
    }
    let ncdp = countparval("cdp");
    warn!("This file has %i CDPs.", ncdp);
    if ncdp > 0 {
        if countparname("vnmo") != ncdp { err!("A vnmo set must be specified for each cdp"); }
        if countparname("tnmo") != ncdp { err!("A tnmo set must be specified for each cdp"); }
    } else {
        err!("A cdp value must be supplied for each tnmo-vnmo set");
    }
    let mut cdp = vec![0i32; ncdp as usize];
    if getparint("cdp", &mut cdp) == 0 {
        err!("A cdp value must be supplied for each tnmo-vnmo set");
    }
    for icdp in 0..ncdp {
        let nvnmo = countnparval(icdp + 1, "vnmo");
        let ntnmo = countnparval(icdp + 1, "tnmo");
        if nvnmo != ntnmo { err!("number of vnmo and tnmo values must be equal"); }
        if nvnmo == 0 { err!("Each cdp must have at least one velocity"); }
        if ntnmo == 0 { err!("Each cdp must have at least one time"); }
        let mut vnmo = vec![0.0f32; nvnmo as usize];
        let mut tnmo = vec![0.0f32; nvnmo as usize];
        if getnparfloat(icdp + 1, "vnmo", &mut vnmo) == 0 { err!("Each cdp must have at least one velocity"); }
        if getnparfloat(icdp + 1, "tnmo", &mut tnmo) == 0 { err!("Each cdp must have at least one time"); }
        for k in 1..ntnmo as usize {
            if tnmo[k] <= tnmo[k-1] {
                warn!("tnmo values must increase for use in NMO");
                warn!("For cdp=%i, check times %g and %g.", cdp[icdp as usize], tnmo[k-1], tnmo[k]);
            }
        }
        if mode == 2 {
            let outfile = sprintf!("%s%s%i", prefix, ".", cdp[icdp as usize]);
            mktvfile(&outfile, ntnmo, &tnmo, &vnmo);
        }
    }
    warn!("End of cdp-tnmo-vnmo check.");
    exit(CWP_Exit());
}

fn mktvfile(outfile: &str, ntnmo: i32, tnmo: &[f32], vnmo: &[f32]) {
    let fptr = efopen(outfile, "w");
    for i in 0..ntnmo as usize {
        fprintf!(fptr, "%g   %g\n", tnmo[i], vnmo[i]);
    }
    efclose(&fptr);
}
