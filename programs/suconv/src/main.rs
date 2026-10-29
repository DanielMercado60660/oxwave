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

//! SUCONV - convolution with user-supplied filter.
//! Port of su/main/convolution_correlation/suconv.c (Revision 1.19).
//! Credits: CWP: Jack K. Cohen, Michel Dietrich.

use oxwave_core::prelude::*;
use oxwave_f_convolution::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUCONV - convolution with user-supplied filter\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suconv <stdin >stdout  filter= [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters: ONE of\t\t\t\t\t\t",
    " sufile=\t\tfile containing SU trace to use as filter\t",
    " filter=\t\tuser-supplied convolution filter (ascii)\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " panel=0\t\tuse only the first trace of sufile\t\t",
    " \t\t\t=1 convolve corresponding trace in sufile with\t",
    " \t\t\ttrace in input data\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns\t\t\t\t\t",
    " Trace header fields modified: ns\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: It is quietly assumed that the time sampling interval on the\t",
    " single trace and the output traces is the same as that on the traces\t",
    " in the input file.  The sufile may actually have more than one trace,\t",
    " but only the first trace is used in panel=0. In panel=1 the corresponding",
    " trace from the sufile are convolved with its counterpart in the data.\t",
    " Caveat, in panel=1 there have to be at least as many traces in sufile\t",
    " as in the input data. If not, a warning is returned, and later traces\t",
    " in the dataset are returned unchanged.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "\tsuplane | suwind min=12 max=12 >TRACE\t\t\t\t",
    "\tsuconv<DATA sufile=TRACE | ...\t\t\t\t\t",
    " Here, the su data file, \"DATA\", is convolved trace by trace with the",
    " the single su trace, \"TRACE\".\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "\tsuconv<DATA filter=1,2,1 | ...\t\t\t\t\t",
    " Here, the su data file, \"DATA\", is convolved trace by trace with the",
    " the filter shown.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    let mut sutrace = Trace::new();
    let mut filter = Vec::<f32>::new(); /* filter coefficients */
    let mut nfilter: i32 = 0; /* length of input wavelet in samples */
    let mut sufile = String::new();
    let mut fp = None;
    let mut delrtf: i32 = 0; /* delrt from traces in sufile */
    let mut issufile = false;

    /* Initialize */
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut panel = 0;
    if getparint("panel", &mut panel) == 0 {
        panel = 0;
    }

    /* Get info from first trace */
    if gettr(&mut intrace) == 0 {
        err!("can't get first trace");
    }
    let nt = intrace.ns as i32;

    /* Get parameters and set up filter array */
    if getparstring("sufile", &mut sufile) == 0 {
        nfilter = countparval("filter");
        if nfilter == 0 {
            err!("must specify filter= desired filter");
        }
        filter = vec![0.0f32; nfilter as usize];
        getparfloat("filter", &mut filter);
    } else if panel == 0 {
        /* if panel=0 only use the first trace */
        fp = Some(efopen(&sufile, "r"));
        fgettr(fp.as_ref().unwrap(), &mut sutrace);

        /* HD: get delrt from filter sufile */
        delrtf = sutrace.delrt as i32;

        nfilter = sutrace.ns as i32;
        filter = vec![0.0f32; nfilter as usize];
        filter.copy_from_slice(&sutrace.data[..nfilter as usize]);
        issufile = true;
    } else {
        /* if panel=1 use each trace */
        fp = Some(efopen(&sufile, "r"));
        issufile = true;
    }
    checkpars();

    /* Set output trace length */
    let mut ntout = nt + nfilter - 1;

    /* Main loop over traces */
    loop {
        /* if panel=1 and sufile is defined */
        if panel == 1 && issufile {
            // C ignores EOF, reusing the last filter trace on a short file.
            fgettr(fp.as_ref().unwrap(), &mut sutrace);

            /* HD: get delrt from filter sufile */
            delrtf = sutrace.delrt as i32;

            nfilter = sutrace.ns as i32;
            filter = vec![0.0f32; nfilter as usize];
            filter.copy_from_slice(&sutrace.data[..nfilter as usize]);
            ntout = nt + nfilter - 1;
        }

        /* Convolve filter with trace */
        convolve_cwp(
            0,
            &filter[..nfilter as usize],
            0,
            &intrace.data[..nt as usize],
            0,
            &mut outtrace.data[..ntout as usize],
        );

        /* Output filtered trace */
        outtrace.hdr = intrace.hdr;
        outtrace.ns = ntout as u16;

        /* HD: update delrt */
        outtrace.delrt = (intrace.delrt as i32 + delrtf) as i16;
        puttr(&mut outtrace);

        if gettr(&mut intrace) == 0 {
            break;
        }
    }

    exit(CWP_Exit());
}
