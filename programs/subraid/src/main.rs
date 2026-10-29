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

//! SUBRAID - Braid traces from files together and output.
//! Port of su/main/windowing_sorting_muting/subraid.c (Revision 1.1).
//! Credits: Andre Latour; looked at suop2.c code to get started (efopen,fgettr).
//! Original source copyright (c) Colorado School of Mines, 2022. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_header_case::{GetCase, fromhead};
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUBRAID - Braid traces from files together and output.                ",
    "           (combine multiple ordered inputs into 1 ordered output)     ",
    "                                                                       ",
    " subraid in1=filea.su in2=fileb.su in3=...   key=cdp >stdout           ",
    "                                                                       ",
    "   inMM= An input file name. Where MM is 1 or 2 digits (0 to 99).      ",
    "         Therefore 100 files can be input in one setup.                ",
    "         Examples: in0= in1= in33= in99=                               ",
    "         There must be at least 2 input files.                         ",
    "                                                                       ",
    "   key=  Keys that indicate the order of the traces in every input.    ",
    "         Traces must be in this order. Order is checked on-the-fly     ",
    "         and an error-halt will occur if input trace order is wrong.   ",
    "         This parameter can be a list. For instance:  key=cdp,offset   ",
    "         means the traces are in offset order within cdps.             ",
    "         A negative sign in front of the name means decreasing order   ",
    "         for the values of that key.                                   ",
    "         Note: If the key= values are equal, the traces from the       ",
    "               input file with lower MM number are output first.       ",
    " \t\t\t\t\t\t\t\t\t",
    "  trac=1 Reset tracl and tracr keys to the output trace count          ",
    "         (really only done because the SUSORT program resets them).    ",
    "      =0 Do not reset. Copy as-is to output.                           ",
    " \t\t\t\t\t\t\t\t\t",
    " This program can combine multiple cdp ordered input files into one    ",
    " cdp ordered output file. It can also combine multiple input files     ",
    " of any other order into one output file in same order. This is done   ",
    " by a process conceptually similar to braiding hair or braiding a rope.",
    "                                                                       ",
    " Basically, you tell this program the names of the su keys which       ",
    " indicate the input trace order. This program opens all input files    ",
    " and outputs whichever trace is lowest. It then reads another trace    ",
    " from that input file and again finds lowest trace for all input files.",
    " It repeats the process until no traces remain on any input file.      ",
    " (For traces with the same key values, the trace from lowest inMM= is  ",
    " is output first. All traces are always output - except on error-halt).",
    "                                                                       ",
    " So, for cdp order, if file1 has 3 traces in cdp 17 and file2 has      ",
    " 5 traces in cdp 17, then the output will have 8 traces in cdp 17.     ",
    " It does not matter if a cdp only exists in one file, and it does      ",
    " not matter if there are gaps in cdp numbers in the files.             ",
    "                                                                       ",
    " Notes:                                                                ",
    "   1. See tests (examples) in src/demos/Utilities/Subraid              ",
    "   2. This program only stores 1 trace at a time from each input file. ",
    "   3. Trace order is checked on-the-fly to make sure it matches        ",
    "      your key= list. If it does not match, error-halt occurs after    ",
    "      output of the last trace which IS in the correct order.          ",
    "   4. Only trace order is checked. No other checks at all. You can     ",
    "      easily cause bad results by inputting files where the traces have",
    "      been procssed differently. Or by processing data in pieces when  ",
    "      it should be kept together (surface consistant analysis, etc.).  ",
    "   5. For clarify, remember that being in order does not necessarily   ",
    "      mean that traces have been through a sort program. This is fine, ",
    "      subraid.c does not care how the su key values got into order.    ",
    "                                                                       ",
    " EXAMPLE USES:                                                         ",
    "   1. Imagine that your survey is so large you have trouble sorting it ",
    "      to cdp order. Just sort howevermany shot gathers you want into   ",
    "      cdp order. Then sort another set of shot gathers into cdp order. ",
    "      Keep doing that until you have sorted all your traces.           ",
    "      Then use this program to combine the multiple cdp ordered files  ",
    "      into one output file. Note for both 2d and 3d, these cdp sorted  ",
    "      input files have over-lapping cdps. Some cdps only exist in      ",
    "      one of the input files and some cdps in multiple input files.    ",
    "   2. You may want to deliberately partition your traces for different ",
    "      processing. This can be done by sucleave (and other programs).   ",
    "      You can cleave your traces into offset ranges, while retaining   ",
    "      cdp order in those offset-ranged files. Then use subraid to      ",
    "      bring them back together.                                        ",
    "   3. In streamer-cable-marine you may want to partition your shot     ",
    "      gathers by cables (or port-starboard shot, cable combinations).  ",
    "      This is sometimes called nominal subsurface line processing.     ",
    "      This program can help put these kinds of data partitions back    ",
    "      together while reducing resource utilization.                    ",
    "                                                                       ",
];

fn compdouble(p1: &[f64], p2: &[f64]) -> i32 {
    for n in 0..p1.len() {
        if p1[n] < p2[n] { return -1; }
        if p1[n] > p2[n] { return 1; }
    }
    0
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let nkeys = countparval("key");
    if nkeys < 1 { err!("key= list must have at least 1 name"); }
    let mut key = Vec::<String>::new();
    getparstringarray("key", &mut key);
    let mut kase = Vec::new();
    let mut ddir = Vec::new();
    let mut dprev = Vec::new();
    for k in &key {
        let (name, dir) = if let Some(name) = k.strip_prefix('+') { (name, 1.0) }
            else if let Some(name) = k.strip_prefix('-') { (name, -1.0) }
            else { (k.as_str(), 1.0) };
        let case = GetCase(name);
        if case < 1 { err!("key= list contains unknown name (%s)", k); }
        kase.push(case);
        ddir.push(dir);
        dprev.push(-1.0e99);
    }
    let mut itrac = 1;
    if getparint("trac", &mut itrac) == 0 { itrac = 1; }
    if itrac > 1 || itrac < 0 { err!("error: trac= is out of range."); }
    let mut fps = Vec::new();
    let mut atrace = Vec::new();
    let mut dvals: Vec<Vec<f64>> = Vec::new();
    let mut idf = Vec::new();
    for i in 0..100 {
        let pname = format!("in{}", i);
        if countparname(&pname) > 0 {
            let mut fname = String::new();
            getparstring(&pname, &mut fname);
            let fp = efopen(&fname, "r");
            let mut tr = Trace::new();
            if fgettr(&fp, &mut tr) == 0 {
                err!("could not read first trace in file input by %s=", "in0");
            }
            let mut vals = Vec::new();
            for j in 0..nkeys as usize {
                vals.push(ddir[j] * fromhead(&tr, kase[j]));
            }
            fps.push(fp);
            atrace.push(tr);
            dvals.push(vals);
            idf.push(i);
        }
    }
    checkpars();
    if fps.len() < 2 { err!("At least 2 input files are needed."); }
    let mut ntrac: i32 = 0;
    while !fps.is_empty() {
        let mut ismall = 0;
        for i in 1..fps.len() {
            if compdouble(&dvals[i], &dvals[ismall]) < 0 { ismall = i; }
        }
        if compdouble(&dprev, &dvals[ismall]) > 0 {
            err!("traces out-of-order in file input by %s=", format!("in{}", idf[ismall]));
        }
        dprev.copy_from_slice(&dvals[ismall]);
        if itrac > 0 {
            ntrac = ntrac.wrapping_add(1);
            atrace[ismall].tracl = ntrac;
            atrace[ismall].tracr = ntrac;
        }
        puttr(&mut atrace[ismall]);
        if fgettr(&fps[ismall], &mut atrace[ismall]) != 0 {
            for j in 0..nkeys as usize {
                dvals[ismall][j] = ddir[j] * fromhead(&atrace[ismall], kase[j]);
            }
        } else {
            efclose(&fps[ismall]);
            fps.remove(ismall);
            atrace.remove(ismall);
            dvals.remove(ismall);
            idf.remove(ismall);
        }
    }
    exit(CWP_Exit());
}
