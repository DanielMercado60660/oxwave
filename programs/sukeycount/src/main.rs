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

//! SUKEYCOUNT - count consecutive traces with the same header value.
//! Port of su/main/headers/sukeycount.c (Revision 1.3).
//! Credits (from the C source): MTU: David Forel, Jan 2005.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                             ",
    " SUKEYCOUNT - sukeycount writes a count of a selected key    ",
    "                                                             ",
    "   sukeycount key=keyword < infile [> outfile]                  ",
    "                                                             ",
    " Required parameters:                                        ",
    " key=keyword      One key word.                                 ",
    "                                                             ",
    " Optional parameters:                                        ",
    " verbose=0  quiet                                            ",
    "        =1  chatty                                           ",
    "                                                             ",
    " Writes the key and the count to the terminal or a text      ",
    "   file when a change of key occurs. This does not provide   ",
    "   a unique key count (see SUCOUNTKEY for that).             ",
    " Note that for key values  1 2 3 4 2 5                       ",
    "   value 2 is counted once per occurrence since this program ",
    "   only recognizes a change of key, not total occurrence.    ",
    "                                                             ",
    " Examples:                                                   ",
    "    sukeycount < stdin key=fldr                              ",
    "    sukeycount < stdin key=fldr > out.txt                    ",
    "                                                             ",
];

fn main() {
    let argv = args();
    let argc = argv.len();
    let mut tr = Trace::new();
    let mut key = vec![String::new(); SU_NKEYS];
    let nkeys: i32;
    let mut countkey = 0usize;
    let mut verbose = 0i32;
    let mut sortold = -99999.0f32;
    let mut isortold = -99999i32;
    let mut gatherkount = 0i32;
    let mut itotal = 0i32;
    let mut type1 = "";

    initargs(&argv);
    requestdoc(1, SDOC);
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    let n = countparval("key");
    if n != 0 {
        nkeys = n;
        getparstringarray("key", &mut key);
    } else {
        if argc == 1 { err!("must set one key value!"); }
        for iarg in 1..argc {
            let keyword = &argv[iarg];
            if verbose != 0 { warn!("argv=%s", argv[iarg]); }
            if !keyword.as_bytes().starts_with(b"output=") {
                key[countkey] = keyword.clone();
                countkey += 1;
            }
            if countkey == 0 { err!("must set one key value!"); }
            if countkey > 1 { err!("must set only one key value!"); }
        }
        nkeys = countkey as i32;
    }
    checkpars();
    if nkeys > 1 { err!("must set only one key value!"); }
    printf!("\n");
    while gettr(&mut tr) != 0 {
        let mut vsort = Value::default();
        gethdval(&tr, &key[0], &mut vsort);
        type1 = hdtype(&key[0]);
        if type1 == "f" {
            let sort = vtof(type1, vsort);
            if itotal == 0 { sortold = sort; }
            if sort != sortold {
                printf!(" %8s = %f", key[0], sortold);
                printf!("      has  %d  trace(s)", gatherkount);
                printf!("\n");
                sortold = sort;
                gatherkount = 0;
            }
        } else {
            let isort = vtoi(type1, vsort);
            if itotal == 0 { isortold = isort; }
            if isort != isortold {
                printf!(" %8s = %d", key[0], isortold);
                printf!("      has  %d  trace(s)", gatherkount);
                printf!("\n");
                isortold = isort;
                gatherkount = 0;
            }
        }
        gatherkount = gatherkount.wrapping_add(1);
        itotal = itotal.wrapping_add(1);
    }
    if type1 == "f" {
        printf!(" %8s = %f", key[0], sortold);
    } else {
        printf!(" %8s = %d", key[0], isortold);
    }
    printf!("      has  %d  trace(s)", gatherkount);
    printf!("\n");
    printf!("\n        %d trace(s) read\n\n", itotal);
    exit(CWP_Exit());
}
