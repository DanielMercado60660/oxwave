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

// Copyright (c) Colorado School of Mines, 2011.
// All rights reserved.

//! SUREFCON - convolution of forward and reverse refraction shots.
//! Port of su/main/convolution_correlation/surefcon.c (Revision 1.7).
//! Credits (from the C source): CWP: Jack K. Cohen, Michel Dietrich;
//! UNSW: D. Palmer, K.T. LEE. Based on suconv.
//! CAVEATS: no space-variable or time-variable capacity.
//! The more than one trace allowed in sufile is the beginning of a hook
//! to handle the spatially variant case.
//!
//! References (from the C source):
//! Palmer, D, 2001a, Imaging refractors with the convolution section,
//! Geophysics 66, 1582-1589.
//! Palmer, D, 2001b, Resolving refractor ambiguities with amplitudes,
//! Geophysics 66, 1590-1593.
//! Exploration Geophysics (2005) 36, 18-25; Butsuri-Tansa (Vol. 58, No.1);
//! Mulli-Tamsa (Vol. 8): A simple approach to refraction statics with the
//! Generalized Main Reciprocal Method and the Refraction Convolution
//! Section Heading, by Derecke Palmer and Leonie Jones.

use oxwave_core::prelude::*;
use oxwave_f_convolution::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUREFCON -  Convolution of user-supplied Forward and Reverse\t\t",
    "\t\trefraction shots using XY trace offset in reverse shot\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tsurefcon <forshot sufile=revshot  xy=trace offseted  >stdout\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t \t",
    " sufile=\tfile containing SU trace to use as reverse shot\t\t",
    " xy=\t\tNumber of traces offseted from the 1st trace in sufile\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t \t",
    " none\t\t\t\t\t\t\t\t \t",
    "\t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns\t\t\t\t\t",
    " Trace header fields modified: ns\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This code implements the Refraction Convolution Section (RCS)\tmethod\t",
    " of generalized reciprocal refraction traveltime analysis developed by ",
    " Derecke Palmer and Leoni Jones.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The time sampling interval on the output traces is half of that on the",
    " traces in the input files.\t\t  \t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t surefcon <DATA sufile=DATA xy=1 | ...\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Here, the su data file, \"DATA\", convolved the nth trace by\t\t",
    " (n+xy)th trace in the same file\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    let mut sutrace = Trace::new();

    /* Initialize */
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get info from first trace */
    if gettr(&mut intrace) == 0 {
        err!("can't get 1st reverse shot trace");
    }
    let nt = intrace.ns as i32;

    /* Default parameters; User-defined overrides */
    let mut xy = 0;
    if getparint("xy", &mut xy) == 0 {
        xy = 0;
    }

    /* Check xy values */
    if xy < 0 {
        err!("xy=%d should be positive", xy);
    }
    let mut sufile = String::new();
    if getparstring("sufile", &mut sufile) == 0 {
        err!("must specify sufile= desired forward shot");
    } else {
        checkpars();

        /* Get parameters and set up forshot array */
        let fp = efopen(&sufile, "r");
        for _itr in 0..=xy {
            if fgettr(&fp, &mut sutrace) == 0 {
                err!("can't get 1st requested forward trace");
            }
        }
        let nforshot = sutrace.ns as i32;
        let mut forshot = vec![0.0f32; nforshot as usize];

        /* Set output trace length */
        let ntout = nt + nforshot - 1;
        // The C passes ntout to a fixed SU_NFLTS output array. Reject an
        // oversized convolution before writing past it, like the bounds-fixed C.
        if ntout > SU_NFLTS as i32 {
            err!(
                "output trace length %d exceeds SU_NFLTS=%d",
                ntout,
                SU_NFLTS as i32
            );
        }

        /* Main loop over reverse shot traces */
        loop {
            warn!("rev==%d\t , for=%d", intrace.tracf, sutrace.tracf);
            forshot.copy_from_slice(&sutrace.data[..nforshot as usize]);

            /* Convolve for shot with reverse shot trace */
            convolve_cwp(
                0,
                &forshot,
                0,
                &intrace.data[..nt as usize],
                0,
                &mut outtrace.data[..ntout as usize],
            );

            /* Output convolved trace */
            outtrace.hdr = intrace.hdr;
            outtrace.ns = ntout as u16;
            outtrace.dt /= 2;
            /* outtrace.cdp = 2*intrace.tracf + xy; */
            fprintf!(
                stderr(),
                "out_cdp=%d\n",
                intrace.tracf.wrapping_mul(2).wrapping_add(xy)
            );
            puttr(&mut outtrace);

            if gettr(&mut intrace) == 0 || fgettr(&fp, &mut sutrace) == 0 {
                break;
            }
        }
    }

    exit(CWP_Exit());
}
