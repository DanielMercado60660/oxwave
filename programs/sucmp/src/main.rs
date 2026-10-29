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

//! SUCMP - compare two seismic data sets.
//! Port of su/main/attributes_parameter_estimation/sucmp.c (Revision 1.5).
//! Credits (from the C source): Reginald H. Beardsley, rhb@acm.org.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUCMP   - CoMPare two seismic data sets, returns 0 to the shell\t",
    "             if the same and 1 if different\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  sucmp file_A file_B\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "      none\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   Optional parameters:\t\t\t\t\t\t",
    "      limit=1.0e-4    normalized difference threshold value\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This program is the seismic equivalent of the Unix cmp(1)\t\t",
    " command.  However, unlike cmp(1), it understands seismic data\t\t",
    " and will consider files which have only small numerical\t\t",
    " differences to be the same.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Sucmp first checks that the number of traces and number of samples\t",
    " are the same. It then compares the trace headers bit for bit.\t\t",
    " Finally it checks that the fractional difference of A & B is\t\t",
    " less than limit.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " This program is intended as an aid in regression testing changes to\t",
    " seismic processing programs.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Expected usage is in shell scripts or Makefiles, e.g.\t\t\t",
    "   #!/bin/sh\t\t\t\t\t\t\t\t",
    "    #-------------------------------------------------------\t\t",
    "    # Run a test data set and verify the result is correct\t\t",
    "    # If the data doesn't match show the data on the screen.\t\t",
    "   #-------------------------------------------------------\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  ./fubar par=tst1.par\t\t\t\t\t\t\t",
    "   sucmp tst1.su ref/tst1.su\t\t\t\t\t\t",
    "   if [ $? ]\t\t\t\t\t\t\t\t",
    "      then\t\t\t\t\t\t\t\t",
    "      suxwigb <tst1.su &\t\t\t\t\t\t",
    "      suxwigb <ref/tst1.su &    fi\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut trace_a = Trace::new();
    let mut trace_b = Trace::new();

    /* Initialize */
    initargs(&argv);
    requestdoc(2, SDOC); /* two file args required */
    let mut limit: f32 = 0.0;
    if getparfloat("limit", &mut limit) == 0 {
        limit = 1.0e-4;
    }
    checkpars();

    // C promotes limit to double for the bounds, then stores each in float.
    let upper = (1.0f64 + limit as f64) as f32;
    let lower = (1.0f64 - limit as f64) as f32;

    /* Open input files */
    let fp_a = efopen(&argv[1], "rb");
    let fp_b = efopen(&argv[2], "rb");

    /* Compare files */
    let mut j: i32 = 0;
    while !feof(&fp_a) && !feof(&fp_b) {
        let n_a = fgettr(&fp_a, &mut trace_a);
        let n_b = fgettr(&fp_b, &mut trace_b);
        j = j.wrapping_add(1);
        if n_a != n_b {
            printf!("Files %s & %s differ at trace %d\n", argv[1], argv[2], j);
            exit(1);
        } else if trace_a.hdr.as_bytes() != trace_b.hdr.as_bytes() {
            // Like memcmp, compare all header bits, including unused fields.
            printf!(
                "Files %s & %s differ in headers at trace %d\n",
                argv[1],
                argv[2],
                j
            );
            exit(1);
        } else {
            let nt = trace_a.ns as usize;
            for i in 0..nt {
                let a = trace_a.data[i];
                let b = trace_b.data[i];
                // Products stay in float. Zero (or NaN) in A is ignored in C.
                if (a > 0.0 && (a < b * lower || a > b * upper))
                    || (a < 0.0 && (a > b * lower || a < b * upper))
                {
                    printf!("Files %s & %s differ at", argv[1], argv[2]);
                    printf!(" Trace: %d Sample: %d\n", j, i as i32);
                    printf!("   A: %15g   B: %15g\n", a, b);
                    exit(1);
                }
            }
        }
    }
    exit(CWP_Exit());
}
