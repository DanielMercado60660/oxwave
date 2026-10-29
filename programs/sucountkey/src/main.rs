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

//! SUCOUNTKEY - count unique header values.
//! Port of su/main/headers/sucountkey.c (Revision 1.6).
//! Credits: Baoniu Han, November 2000.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUCOUNTKEY - COUNT the number of unique values for a given KEYword.\t",
    "\t\t\t\t\t\t\t\t\t",
    " sucountkey < input.su key=[sx,gx,cdp,...]\t\t\t\t",
    " Required parameter:\t\t\t\t\t\t\t",
    " key=\t\t\tarray of SU header keywords being counted\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " verbose=1\t\tchatty, =0 just print keyword number\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    "\t  suplane | sucountkey key=tracl,tracr,offset\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut key: Vec<String> = vec![String::new(); SU_NKEYS];
    let mut verbose: i32 = 0;

    initargs(&argv);
    requestdoc(1, SDOC);

    let nkeys = countparval("key");
    if nkeys != 0 {
        getparstringarray("key", &mut key);
    } else {
        key[0] = "cdp".to_string();
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 1;
    }
    checkpars();

    let buf_size = 100000usize;
    let mut types = Vec::new();
    let mut index = Vec::new();
    let mut count_table: Vec<Vec<Value>> = Vec::new();
    for ikey in 0..nkeys as usize {
        types.push(hdtype(&key[ikey]));
        index.push(getindex(&key[ikey]));
        count_table.push(Vec::with_capacity(buf_size));
    }

    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    for ikey in 0..nkeys as usize {
        let mut val = Value::default();
        gethval(&tr, index[ikey], &mut val);
        count_table[ikey].push(val);
    }

    while gettr(&mut tr) != 0 {
        for ikey in 0..nkeys as usize {
            let mut val0 = Value::default();
            gethval(&tr, index[ikey], &mut val0);
            let mut i = 0;
            while i < count_table[ikey].len() {
                if valcmp(types[ikey], val0, count_table[ikey][i]) == 0 {
                    break;
                }
                i += 1;
            }
            if i == count_table[ikey].len() {
                let mut val = Value::default();
                gethval(&tr, index[ikey], &mut val);
                count_table[ikey].push(val);
            }
        }
    }

    if verbose == 1 {
        warn!("The counting of key word is finished, the total number for specified header is as follow:");
    }
    for ikey in 0..nkeys as usize {
        if verbose == 1 {
            warn!(" key %s has %d unique entries in this data", key[ikey], count_table[ikey].len() as i32);
        } else {
            fprintf!(stdout(), "%s %d ", key[ikey], count_table[ikey].len() as i32);
        }
    }
    exit(CWP_Exit());
}
