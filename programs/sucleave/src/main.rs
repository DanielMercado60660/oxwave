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

//! SUCLEAVE - cleave traces into multiple output files.
//! Port of su/main/windowing_sorting_muting/sucleave.c (Revision 1.2).
//! Credits: Andre Latour; started by looking at suop2.c (efopen, fgettr).
//! Original source copyright (c) Colorado School of Mines, 2022. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_header_case::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUCLEAVE - Cleave traces from 1 input into multiple output files.     ",
    "           (split to multiple output files, retaining input order).    ",
    "                                                                       ",
    " sucleave <stdin key=offset low=50 size=100 abs=1                      ",
    "                                                                       ",
    " key=offset  name of key used to seperate traces to output files.      ",
    " \t\t\t\t\t\t\t\t\t",
    "    abs=1    Use absolute value of key.                                ",
    "       =0    Do not use absolute value of key.                         ",
    " \t\t\t\t\t\t\t\t\t",
    "  size=100   size of key ranges for output files. Must be > 0.         ",
    " \t\t\t\t\t\t\t\t\t",
    "  low=size/2 Base value in lowest magnitude range for output files.    ",
    "             The default is size/2. Must be <= size/2 and > 0.         ",
    "             For each output file the centre of the range of traces    ",
    "             output to it is low+size*N (where N is an integer).       ",
    "              ( The lower boundary is included in the N range,         ",
    "               but the higher boundary belongs to the N+1 range ).     ",
    " \t\t\t\t\t\t\t\t\t",
    " high=100000 base value in highest magnitude range for output files.   ",
    "             (This value is adjusted to a multiple of low+size*N).     ",
    "             If abs=1 the range limits are from  low  to high.         ",
    "             If abs=0 the range limits are from -high to high.         ",
    " \t\t\t\t\t\t\t\t\t",
    "  outbase=   This value is used as part of the output file names.      ",
    "             By default outbase is cleave_key_ where key is key= name. ",
    "             The nearest integer of the range base value (low+size*N)  ",
    "             is always appended after outbase. For example:            ",
    "               cleave_offset_50.su   cleave_offset_150.su              ",
    " \t\t\t\t\t\t\t\t\t",
    "  print=0    Do not print.                                             ",
    "       =1    Print filename, trace count, average key value for        ",
    "             each created output file. Note that ranges with 0 traces  ",
    "             will make no print.                                       ",
    " \t\t\t\t\t\t\t\t\t",
    "  ATTENTION: No traces are deleted. All traces less than the minimum   ",
    "             key range are written to one extra file, and all traces   ",
    "             greater than the maximum key range are written to another ",
    "             extra file.     WARNINGS are printed if this occurs.      ",
    "             These extra files can be recognized because the number in ",
    "             their names is outside the low,high ranges you specified. ",
    " \t\t\t\t\t\t\t\t\t",
    "  NOTES                                                                ",
    "   1. See simple tests (examples) in src/demos/Utilities/Sucleave      ",
    "   2. Output files are not created until a trace is written to them.   ",
    "      So, specifying excessive low= and high= values has little        ",
    "      consequences except for some small memory allocation herein.     ",
    "   3. Bad choice of size can create an output file for every trace.    ",
    " \t\t\t\t\t\t\t\t\t",
    " EXAMPLE USES:                                                         ",
    "   1. You may want to deliberately partition your traces for different ",
    "      processing. You can cleave your traces into offset ranges, while ",
    "      retaining cdp order in those offset-ranged files.                ",
    "      You can then use subraid to bring them back together.            ",
    "   2. In streamer-cable-marine you may want to partition your shot     ",
    "      gathers by cables (or port-starboard shot, cable combinations).  ",
    "      This is sometimes called nominal subsurface line processing.     ",
    "      You can then use subraid to bring them back together.            ",
    "                                                                       ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key = String::new();
    if getparstring("key", &mut key) == 0 { key = "offset".into(); }
    let kase = GetCase(&key);
    if kase < 1 { err!("error: key= is not a recognized header key name"); }
    let mut iabs = 1;
    getparint("abs", &mut iabs);
    if iabs > 1 || iabs < 0 { err!("error: abs= is out of range."); }
    let mut dsize = 100.0;
    getpardouble("size", &mut dsize);
    if dsize <= 0.0 { err!("error: size= must be greater than zero."); }
    let mut dlow = dsize / 2.0;
    getpardouble("low", &mut dlow);
    if dlow < 0.0 { err!("error: low= must be greater than zero."); }
    if dlow > dsize / 2.0 { err!("error: low= must be less than or equal to size/2."); }
    let mut dhigh = 100000.0;
    getpardouble("high", &mut dhigh);
    if dlow > dhigh { err!("error: high= must be greater than low=."); }
    let mut outbase = String::new();
    if getparstring("outbase", &mut outbase) == 0 { outbase = format!("cleave_{}_", key); }
    let mut iprint = 0;
    getparint("print", &mut iprint);
    if iprint > 1 || iprint < 0 { err!("error: print= is out of range."); }
    checkpars();
    let file_count = 3.50000000001 + (dhigh - dlow) / dsize;
    if !file_count.is_finite() || file_count < 1.0 || file_count >= i32::MAX as f64 + 1.0 {
        err!("error: too many output file ranges.");
    }
    let mut nfile = file_count as i32;
    if iabs == 1 { dlow -= dsize; } else {
        nfile -= 1;
        dlow -= dsize * nfile as f64;
        nfile = match nfile.checked_mul(2) {
            Some(n) if n > 0 => n,
            // The wrapped negative C count reaches ealloc1int as an
            // impossibly large allocation on the reference build.
            _ => err!("ealloc.c: malloc failed (Cannot allocate memory)"),
        };
    }
    let mut fps: Vec<Option<CFile>> = Vec::new();
    let mut icount = Vec::new();
    let mut dsum = Vec::new();
    let count = nfile as usize;
    if fps.try_reserve_exact(count).is_err() || icount.try_reserve_exact(count).is_err()
        || dsum.try_reserve_exact(count).is_err() {
        err!("error: cannot allocate output file ranges.");
    }
    fps.resize_with(count, || None);
    icount.resize(count, 0i32);
    dsum.resize(count, 0.0f64);
    let mut minf = 2147483645i32;
    let mut maxf = -2147483645i32;
    if iprint == 1 { warn!(" Output file         trace count      key average "); }
    while gettr(&mut tr) != 0 {
        // fromhead has no defined value for case numbers outside its header cases.
        if kase > 80 { err!("error: key= is not a recognized header key name"); }
        let mut dval = fromhead(&tr, kase);
        if dval < 0.0 && iabs > 0 { dval = 0.0 - dval; }
        let mut ifile = (0.50000000001 + (dval - dlow) / dsize) as i32;
        if ifile > nfile - 1 { ifile = nfile - 1; }
        else if ifile < 1 { ifile = 0; }
        let idx = ifile as usize;
        if icount[idx] < 1 {
            let n = (ifile as f64 * dsize + dlow + 0.00000000001).round_ties_even() as i64 as i32;
            let fname = format!("{}{}.su", outbase, n);
            fps[idx] = Some(efopen(&fname, "w"));
            if ifile < minf { minf = ifile; }
            if ifile > maxf { maxf = ifile; }
        }
        fputtr(fps[idx].as_ref().unwrap(), &mut tr);
        icount[idx] = icount[idx].wrapping_add(1);
        dsum[idx] += dval;
    }
    for i in minf..=maxf {
        let idx = i as usize;
        if icount[idx] > 0 {
            efclose(fps[idx].as_ref().unwrap());
            if iprint == 1 {
                dsum[idx] /= icount[idx] as f64;
                let n = (i as f64 * dsize + dlow + 0.00000000001).round_ties_even() as i64 as i32;
                let fname = format!("{}{}.su", outbase, n);
                warn!("%s       %d       %f", fname, icount[idx], dsum[idx]);
            }
        }
    }
    if minf == 0 {
        let n = dlow.round_ties_even() as i64 as i32;
        let fname = format!("{}{}.su", outbase, n);
        warn!("%d traces have key < low range. Output in extra file %s", icount[0], fname);
    }
    if maxf == nfile - 1 {
        let n = (maxf as f64 * dsize + dlow).round_ties_even() as i64 as i32;
        let fname = format!("{}{}.su", outbase, n);
        warn!("%d traces have key > high range. Output in extra file %s", icount[maxf as usize], fname);
    }
    exit(CWP_Exit());
}
