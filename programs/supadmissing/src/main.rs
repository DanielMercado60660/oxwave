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

//! SUPADMISSING - Add dead traces for missing key values.
//! Port of su/main/windowing_sorting_muting/supadmissing.c (Revision 1.0).
//! Credits: Andre Latour, February 2024.

use oxwave_core::prelude::*;
use oxwave_header_case::{GetCase, fromhead, tohead};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUPADMISSING - Add Dead Traces For Missing Key Values.                ",
    "                                                                       ",
    " supadmissing <stdin >stdout                                           ",
    "                                                                       ",
    " keyloc=cdp  Input traces must be in increasing order of this key.     ",
    " \t\t\t\t\t\t\t\t\t",
    " incval=1.0  Increment. Must be greater than 0.0                       ",
    "             Input traces should have keyloc values at whole number    ",
    "             multiples of incval added to minval (minval + incval*N).  ",
    "             So if minval=2 and incval=5 then 2,7,12,17...are expected.",
    "             Any missing values will cause a padded trace to be output.",
    "             The padded trace has keyloc set to the missing value.     ",
    "             Padded traces have all amplitudes equal to 0.0            ",
    "             Padded traces have nhs=0 set but all other key values     ",
    "             are copied from the next higher input trace (or last).    ",
    "        ***  This means the padded header values are NOT interpolated. ",
    "                                                                       ",
    " minval=     Minimum keyloc value to output. Input traces with keyloc  ",
    "             less than this are not output. If first input trace is    ",
    "             greater than this value, padding occurs up to this value. ",
    "             Default is keyloc value from first input trace.           ",
    "        ***  This default is often NOT what you want to do.            ",
    " \t\t\t\t\t\t\t\t\t",
    " maxval=     Maximum keyloc value to output. Input traces with keyloc  ",
    "             greater than this are not output. If last input trace is  ",
    "             less than this value, padding occurs up to this value.    ",
    "             Default is keyloc value from last input trace.            ",
    "        ***  This default is often NOT what you want to do.            ",
    " \t\t\t\t\t\t\t\t\t",
    " NOTE: Input trace range can be completely before or completely after  ",
    "       the minval,maxval range in which case no input trace is output. ",
    "       But padded traces are still output from minval to maxval using  ",
    "       the header values from the first or last input trace.           ",
    " \t\t\t\t\t\t\t\t\t",
    " NOTE: This program does not honour scalco and scalel. So if you use   ",
    "       any key associated with them as keyloc, then specify incval,    ",
    "       minval, maxval based on the actual values in the trace headers. ",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut trz = Trace::new();
    let mut keyloc = String::new();
    let mut incval = 1.0f64;
    let mut minval = 0.0f64;
    let mut maxval = 0.0f64;
    initargs(&argv);
    requestdoc(1, SDOC);
    if countparval("keyloc") > 1 { err!("**** Error: keyloc cannot be multiple keys."); }
    if countparval("keyloc") > 0 { getparstring("keyloc", &mut keyloc); }
    else { keyloc = "cdp".to_string(); }
    let kaseloc = GetCase(&keyloc);
    if kaseloc < 1 { err!("**** Error: Specified keyloc name %s is not recognized.", keyloc); }
    if getpardouble("incval", &mut incval) == 0 { incval = 1.0; }
    if incval <= 0.0 { err!("**** Error: incval must be greater than 0.0 "); }
    let minhave = countparval("minval") > 0;
    if minhave { getpardouble("minval", &mut minval); }
    let maxhave = countparval("maxval") > 0;
    if maxhave { getpardouble("maxval", &mut maxval); }
    checkpars();
    if gettr(&mut tr) == 0 { err!("Error: cannot get first trace"); }
    trz.hdr = tr.hdr;
    for i in 0..trz.ns as usize { trz.data[i] = 0.0; }
    trz.nhs = 0;
    let mut dval = fromhead(&tr, kaseloc);
    if !minhave { minval = dval; }
    let mut dpre = dval;
    if !maxhave { maxval = 1.0e30; }
    if minval > maxval { err!("Error: minval cannot be greater than maxval"); }
    if minval > dval {
        while minval > dval {
            if gettr(&mut tr) == 0 {
                tr.nhs = 0;
                tohead(&mut tr, kaseloc, minval);
                for i in 0..tr.ns as usize { tr.data[i] = 0.0; }
                dval = minval;
                break;
            }
            dval = fromhead(&tr, kaseloc);
            if dval <= dpre { err!("Error: keyloc=%s is not increasing. Next trace=%f  Previous=%f", keyloc, dpre, dval); }
            dpre = dval;
        }
        trz.hdr = tr.hdr;
        trz.nhs = 0;
    } else if minval < dval {
        while minval < dval {
            tohead(&mut trz, kaseloc, minval);
            puttr(&mut trz);
            minval += incval;
            if minval > maxval { break; }
        }
    }
    while dval <= maxval {
        if minval < dval {
            trz.hdr = tr.hdr;
            trz.nhs = 0;
            while minval < dval {
                tohead(&mut trz, kaseloc, minval);
                puttr(&mut trz);
                minval += incval;
            }
        }
        puttr(&mut tr);
        if gettr(&mut tr) != 0 {
            dval = fromhead(&tr, kaseloc);
            minval += incval;
            if dval <= dpre { err!("Error: keyloc=%s is not increasing. Next trace=%f  Previous=%f", keyloc, dpre, dval); }
            dpre = dval;
        } else {
            if maxhave {
                trz.hdr = tr.hdr;
                trz.nhs = 0;
                while dval < maxval {
                    dval += incval;
                    tohead(&mut trz, kaseloc, dval);
                    puttr(&mut trz);
                }
            }
            break;
        }
    }
    exit(CWP_Exit());
}
