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

//! SUGETHW - sugethw writes the values of the selected key words.
//! Port of su/main/headers/sugethw.c (SUGETHW: Revision 1.18).
//!
//! Credits (from the C source):
//!	SEP: Shuki Ronen
//!	CWP: Jack K. Cohen
//!      CWP: John Stockwell, added geom stuff, and getparstringarray

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUGETHW - sugethw writes the values of the selected key words\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "   sugethw key=key1,... [output=] <infile [>outfile]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " key=key1,...\t\tAt least one key word.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " output=ascii\t\toutput written as ascii for display\t\t",
    " \t\t\t=binary for output as binary floats\t\t",
    " \t\t\t=geom   ascii output for geometry setting\t",
    " verbose=0 \t\tquiet\t\t\t\t\t\t",
    " \t\t\t=1 chatty\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Output is written in the order of the keys on the command\t\t",
    " line for each trace in the data set.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    " \tsugethw <stdin key=sx,gx\t\t\t\t\t",
    " writes sx, gx values as ascii trace by trace to the terminal.\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Comment: \t\t\t\t\t\t\t\t",
    " Users wishing to edit one or more header field (as in geometry setting)",
    " may do this via the following sequence:\t\t\t\t",
    "     sugethw < sudata output=geom key=key1,key2,... > hdrfile \t\t",
    " Now edit the ASCII file hdrfile with any editor, setting the fields\t",
    " appropriately. Convert hdrfile to a binary format via:\t\t",
    "     a2b < hdrfile n1=nfields > binary_file\t\t\t\t",
    " Then set the header fields via:\t\t\t\t\t",
    "     sushw < sudata infile=binary_file key=key1,key2,... > sudata.edited",
    " \t\t\t\t\t\t\t\t\t",
];

const ASCII: i32 = 0;
const BINARY: i32 = 1;
const GEOM: i32 = 2;

fn main() {
    let argv = args();
    let argc = argv.len();
    let mut key: Vec<String> = vec![String::new(); SU_NKEYS]; /* array of keywords */
    let mut nkeys: i32; /* number of keywords to be gotten */
    let mut countkey: usize = 0; /* counter of keywords in argc loop */
    let mut output = String::new(); /* string representing output format */
    let mut ioutput = ASCII; /* integer representing output format */
    let mut verbose: i32 = 0; /* verbose? */
    let mut tr = Trace::new();

    /* Initialize */
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get key values */
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    nkeys = countparval("key");
    if nkeys != 0 {
        getparstringarray("key", &mut key);
    } else {
        /* support old fashioned method for inputting key fields */
        /* as single arguments:  sugethw key1 key2 ... 		*/
        if argc == 1 {
            err!("must set at least one key value!");
        }

        for iarg in 1..argc {
            let keyword: &String = &argv[iarg]; /* keyword */

            if verbose != 0 {
                warn!("argv=%s", argv[iarg]);
            }
            /* get array of types and indexes to be set */
            if !keyword.as_bytes().starts_with(b"output=") {
                if countkey < key.len() {
                    key[countkey] = keyword.clone();
                } else {
                    key.push(keyword.clone());
                }
                countkey += 1;
            }

            if countkey == 0 {
                err!("must set at least one key value!");
            }
        }
        nkeys = countkey as i32;
    }

    /* Set and check output format; recall ioutput initialized to ASCII */
    if getparstring("output", &mut output) == 0 {
        output = "ascii".to_string();
    }

    checkpars();

    if output == "binary" {
        ioutput = BINARY;
    } else if output == "geom" {
        ioutput = GEOM;
    } else if output != "ascii" {
        err!("unknown format output=\"%s\", see self-doc", output);
    }

    /* Loop over traces writing selected header field values */
    while gettr(&mut tr) != 0 {
        for ikey in 0..nkeys as usize {
            let mut val = Value::default();
            let fval: f32;

            gethdval(&tr, &key[ikey], &mut val);

            match ioutput {
                ASCII => {
                    printf!("%6s=", key[ikey]);
                    printfval(hdtype(&key[ikey]), val);
                    putchar(b'\t' as i32);
                }
                BINARY => {
                    fval = vtof(hdtype(&key[ikey]), val);
                    efwrite(&[fval], &stdout());
                }
                GEOM => {
                    printfval(hdtype(&key[ikey]), val);
                    putchar(b' ' as i32);
                }
                _ => {}
            }
        }

        match ioutput {
            GEOM => printf!("\n"),
            ASCII => printf!("\n\n"),
            _ => {}
        }
    }

    exit(CWP_Exit());
}
