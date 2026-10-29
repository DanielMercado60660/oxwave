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

//! SUSTACKUP - stack traces by header combinations.
//! Port of su/main/stacking/sustackup.c (Revision 1.0).
//! Credits: Andre Latour, March 2023.

use oxwave_core::prelude::*;
use oxwave_header_case::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUSTACKUP - Input Any Trace Order and Stack To Any Key Combination.   ",
    "                                                                       ",
    " sustackup <stdin >stdout keyloc=cdp                                   ",
    "                                                                       ",
    " By default this program inputs any trace order and makes a cdp stack  ",
    " using all input traces. But this program can also make other stacks.  ",
    " Stacks can be made for any user-specified combination of key values.  ",
    " All input traces with the same set of key values are stacked together.",
    " In particular, 2D and 3D shot or receiver stacks can be made.         ",
    " Additionally, two range limits can be specified to exclude traces     ",
    " from stacking. For instance, the offset range can be limited. Range   ",
    " limiting can also be done by suwind, of course. But an important      ",
    " difference here is that, by default, this program outputs a zero      ",
    " amplitude trace with full header values even when the range limits    ",
    " exclude all traces for that shot or receiver.                         ",
    "                                                                       ",
    " keyloc=cdp List of keys that define output stack trace locations.     ",
    "          A different stacked trace is output for each unique          ",
    "          combination of these key values. The output traces have the  ",
    "          header values of the first trace with lowest absolute offset ",
    "          encountered for each keyloc unique value combination.        ",
    "          Except the nhs key is reset to the number of traces within   ",
    "          the range limits (i.e. the stack fold).                      ",
    "    Note: Output header is copied from lowest absolute offset trace    ",
    "          even if it is not within the range limits (specified below). ",
    "      *** Offset is always absolute valued on output (non-negative).   ",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " keyabs=  Key name for absolute value range limit. If specified then   ",
    "          traces are included in stack if within the range limit(s).   ",
    "          Often, the offset key name is used here.                     ",
    " \t\t\t\t\t\t\t\t\t",
    " keyabs2= If specified, this is subtracted from keyabs and the result  ",
    "          absolute valued before applying minabs,maxabs range.         ",
    " \t\t\t\t\t\t\t\t\t",
    " minabs=  Default is no minimum limit. See range specification note.   ",
    " maxabs=  Default is no maximum limit. See range specification note.   ",
    " \t\t\t\t\t\t\t\t\t",
    "  Note RANGE SPECIFICATION:                                            ",
    "          If minabs is less than maxabs then the range to be included  ",
    "          in stack is GREATER OR EQUAL to minabs and LESS than maxabs. ",
    "     ***  But if minabs is greater or equal to maxabs then the range   ",
    "          is GREATER OR EQUAL to minabs -or- LESS than maxabs.         ",
    "          This backwards specification is intended for some situations ",
    "          where the key contains an azimuth or an angle. For instance, ",
    "          if minabs=340 and maxabs=20 then traces >=340 or <20 degrees ",
    "          are in stack range.                                          ",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " keysign= Key name for signed value range limit. If specified then     ",
    "          traces are included in stack if in the minsign,maxsign range.",
    "    Note: If you specify both keyabs and keysign, traces must be in    ",
    "          both ranges to be included in stack.                         ",
    " \t\t\t\t\t\t\t\t\t",
    " keysign2= If specified, this is subtracted from keysign before        ",
    "           applying the minsign,maxsign range.                         ",
    " \t\t\t\t\t\t\t\t\t",
    " minsign= Default is no minimum limit. See range specification note.   ",
    " maxsign= Default is no maximum limit. See range specification note.   ",
    " \t\t\t\t\t\t\t\t\t",
    "  Note RANGE SPECIFICATION:                                            ",
    "          If minsign is less than maxsign then range to be included in ",
    "          stack is GREATER OR EQUAL to minsign and LESS than maxsign.  ",
    "     ***  But if minsign is greater or equal to maxsign then the range ",
    "          is GREATER OR EQUAL to minsign -or- LESS than maxsign.       ",
    "          This backwards specification is intended for some situations ",
    "          where the key contains an azimuth or angle. For instance, if ",
    "          minsign=340 and maxsign=20 then traces >=340 or <20 degrees  ",
    "          are in stack range.                                          ",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " keep=1   Make stacked output traces for all unique keyloc combinations",
    "          even if the range limits eliminates all input traces. These  ",
    "          stacked output traces have the same header values as output  ",
    "          for other traces (but key nhs=0 and all amplitudes are 0.)   ",
    "     =0   Only output stacked traces if at least one trace is in range.",
    "    Note: keep=1 makes it easier to compare the result of different    ",
    "          ranges since the same amount of stacked traces are output.   ",
    "                                                                       ",
    " maxstack=100000 Maximum Amount of Stacked Output Traces. This is also ",
    "                 the maximum unique combinations of keyloc values.     ",
    "           Note: Initial memory obtained is about 24*maxstack bytes.   ",
    "                 Then main memory for stacked traces is only obtained  ",
    "                 as each unique keyloc combination is encountered      ",
    "                 (about 12*trace_samples bytes per unique combination).",
    "                                                                       ",
    "                                                                       ",
    "  Note On Output Headers:                                              ",
    "          Key nhs is set to the number of traces within the ranges of  ",
    "          the unique combination of keyloc values that identify the    ",
    "          stacked output trace (i.e. the stack fold).                  ",
    "          All other keys are set to the values of the first trace      ",
    "          encountered with the smallest absolute offset for the        ",
    "          unique combination of keyloc values that identify the        ",
    "          traces to stack together.                                    ",
    "          The reason to output the header from smallest offset trace   ",
    "          is it means shot-stacks and receiver-stacks have cdp numbers ",
    "          that correspond approximately to their cdp-stack locations.  ",
    "          It also means cdp-stacks from this program have shot and     ",
    "          receiver coordinates and elevations (etc) of near locations. ",
    "                                                                       ",
];

struct Stack { loc: Vec<f64>, trace: Trace, amp: Vec<f64> }

fn comparedd(a: &[f64], b: &[f64]) -> i32 {
    for i in 0..a.len() {
        if a[i] < b[i] { return -1; }
        if a[i] > b[i] { return 1; }
    }
    0
}

fn bhighdd(all: &[Stack], guy: &[f64]) -> usize {
    let (mut low, mut high) = (0, all.len());
    while low < high {
        let mid = low + (high - low) / 2;
        if comparedd(guy, &all[mid].loc) >= 0 { low = mid + 1; } else { high = mid; }
    }
    low
}

fn key(name: &str, kind: &str) -> i32 {
    let kase = GetCase(name);
    if kase < 1 { err!("**** Error: Specified %s name %s is not recognized.", kind, name); }
    kase
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut keyloc = Vec::<String>::new();
    if countparval("keyloc") > 0 { getparstringarray("keyloc", &mut keyloc); }
    else { keyloc.push("cdp".to_string()); }
    let kaseloc: Vec<i32> = keyloc.iter().map(|s| key(s, "keyloc")).collect();
    let (mut kaseabs, mut kaseabs2, mut kasesign, mut kasesign2) = (0, 0, 0, 0);
    let (mut minabs, mut maxabs, mut minsign, mut maxsign) = (-1.0e31f64, 1.0e31f64, -1.0e31f64, 1.0e31f64);
    let mut name = String::new();
    if countparval("keyabs") > 0 {
        getparstring("keyabs", &mut name); kaseabs = key(&name, "keyabs");
    } else {
        if countparval("keyabs2") > 0 { err!("**** Error: specifying keyabs2 requires specifying keyabs."); }
        if countparval("minabs") > 0 { err!("**** Error: specifying minabs requires specifying keyabs."); }
        if countparval("maxabs") > 0 { err!("**** Error: specifying maxabs requires specifying keyabs."); }
    }
    if countparval("keyabs2") > 0 { getparstring("keyabs2", &mut name); kaseabs2 = key(&name, "keyabs2"); }
    if countparval("minabs") > 0 { getpardouble("minabs", &mut minabs); }
    if countparval("maxabs") > 0 { getpardouble("maxabs", &mut maxabs); }
    if countparval("keysign") > 0 {
        getparstring("keysign", &mut name); kasesign = key(&name, "keysign");
    } else {
        if countparval("keysign2") > 0 { err!("**** Error: specifying keysign2 requires specifying keyabs."); }
        if countparval("minsign") > 0 { err!("**** Error: specifying minsign requires specifying keyabs."); }
        if countparval("maxsign") > 0 { err!("**** Error: specifying maxsign requires specifying keyabs."); }
    }
    if countparval("keysign2") > 0 { getparstring("keysign2", &mut name); kasesign2 = key(&name, "keysign2"); }
    if countparval("minsign") > 0 { getpardouble("minsign", &mut minsign); }
    if countparval("maxsign") > 0 { getpardouble("maxsign", &mut maxsign); }
    let (mut keep, mut maxstack) = (1, 100000);
    if countparval("keep") > 0 { getparint("keep", &mut keep); }
    if keep < 0 || keep > 1 { err!("**** Error: keep parameter is out-of-range."); }
    if countparval("maxstack") > 0 { getparint("maxstack", &mut maxstack); }
    if maxstack < 1 { err!("**** Error: maxstack must be greater than 0"); }
    checkpars();
    let mut tr = Trace::new();
    let mut all: Vec<Stack> = Vec::new();
    let (mut nproct, mut lensam) = (0i32, 0usize);
    while gettr(&mut tr) != 0 {
        if nproct < 1 { lensam = tr.ns as usize; }
        nproct += 1;
        let mut ilive = 1i32;
        if kaseabs > 0 {
            let mut dval = fromhead(&tr.hdr, kaseabs);
            if kaseabs2 > 0 { dval -= fromhead(&tr.hdr, kaseabs2); }
            dval = dval.abs();
            if minabs < maxabs { if dval < minabs || dval >= maxabs { ilive = 0; } }
            else if dval < minabs && dval >= maxabs { ilive = 0; }
        }
        if kasesign > 0 {
            let mut dval = fromhead(&tr.hdr, kasesign);
            if kasesign2 > 0 { dval -= fromhead(&tr.hdr, kasesign2); }
            if minsign < maxsign { if dval < minsign || dval >= maxsign { ilive = 0; } }
            else if dval < minsign && dval >= maxsign { ilive = 0; }
        }
        if ilive < 1 && keep < 1 { continue; }
        let guy: Vec<f64> = kaseloc.iter().map(|&k| fromhead(&tr.hdr, k)).collect();
        let ihere = bhighdd(&all, &guy);
        tr.offset = tr.offset.wrapping_abs();
        if ihere > 0 && comparedd(&guy, &all[ihere - 1].loc) == 0 {
            let stack = &mut all[ihere - 1];
            if tr.offset < stack.trace.offset {
                let fold = stack.trace.nhs; stack.trace.hdr = tr.hdr; stack.trace.nhs = fold;
            }
            if ilive > 0 {
                stack.trace.nhs = stack.trace.nhs.wrapping_add(1);
                for i in 0..lensam {
                    if tr.data[i] != 0.0 { stack.amp[i] += tr.data[i] as f64; stack.trace.data[i] += 1.0; }
                }
            }
        } else {
            if all.len() >= maxstack as usize {
                err!("error: At input trace=%d  number of unique keyloc combinations is greater than maxstack=%d", nproct, maxstack);
            }
            let mut trace = tr.clone();
            let mut amp = vec![0.0; lensam];
            for i in 0..lensam {
                if ilive < 1 || tr.data[i] == 0.0 { trace.data[i] = 0.0; }
                else { amp[i] = tr.data[i] as f64; trace.data[i] = 1.0; }
            }
            trace.nhs = ilive as i16;
            all.insert(ihere, Stack { loc: guy, trace, amp });
        }
    }
    let lastloc = all.len();
    for stack in &mut all {
        for i in 0..lensam {
            if stack.trace.data[i] > 0.0 { stack.trace.data[i] = (stack.amp[i] / stack.trace.data[i] as f64) as f32; }
        }
        puttr(&mut stack.trace);
    }
    warn!("Number of traces input=%d  Number stacked output=%d ", nproct, lastloc as i32);
    exit(CWP_Exit());
}
