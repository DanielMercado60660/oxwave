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
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! SUXCOR - correlation with user-supplied filter.
//! Port of su/main/convolution_correlation/suxcor.c (Revision 1.19).
//! Credits: CWP: Jack K. Cohen, Michel Dietrich; CWP: Ttjan, panel correlation;
//! UTK: Rick Williams, vibroseis correlation option.

use oxwave_core::prelude::*;
use oxwave_correlation::xcor;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUXCOR - correlation with user-supplied filter\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suxcor <stdin >stdout  filter= [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters: ONE of\t\t\t\t\t\t",
    " sufile=\t\tfile containing SU traces to use as filter\t",
    " filter=\t\tuser-supplied correlation filter (ascii)\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " vibroseis=0\t\t=nsout for correlating vibroseis data\t\t",
    " first=1\t\tsupplied trace is default first element of\t",
    " \t\t\tcorrelation.  =0 for it to be second.\t\t",
    " panel=0\t\tuse only the first trace of sufile as filter \t",
    " \t\t\t=1 xcor trace by trace an entire gather\t\t",
    " ftwin=0\t\tfirst sample on the first trace of the window \t",
    " \t\t\t\t(only with panel=1)\t\t \t",
    " ltwin=0\t\tfirst sample on the last trace of the window \t",
    " \t\t\t\t(only with panel=1)\t\t \t",
    " ntwin=nt\t\tnumber of samples in the correlation window\t",
    " \t\t\t\t(only with panel=1)\t\t \t",
    " ntrc=48\t\tnumber of traces on a gather \t\t\t",
    " \t\t\t\t(only with panel=1)\t\t \t",
    " \t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns\t\t\t\t\t",
    " Trace header fields modified: ns\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: It is quietly assumed that the time sampling interval on the\t",
    " single trace and the output traces is the same as that on the traces\t",
    " in the input file.  The sufile may actually have more than one trace,\t",
    " but only the first trace is used when panel=0. When panel=1 the number",
    " of traces in the sufile MUST be the same as the number of traces in \t",
    " the input.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "\tsuplane | suwind min=12 max=12 >TRACE\t\t\t\t",
    "\tsuxcor<DATA sufile=TRACE |...\t\t\t\t\t",
    " Here, the su data file, \"DATA\", is correlated trace by trace with the",
    " the single su trace, \"TRACE\".\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "\tsuxcor<DATA filter=1,2,1 | ...\t\t\t\t\t",
    " Here, the su data file, \"DATA\", is correlated trace by trace with the",
    " the filter shown.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Correlating vibroseis data with a vibroseis sweep:\t\t\t",
    " suxcor < data sufile=sweep vibroseis=nsout  |...\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " is equivalent to, but more efficient than:\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suxcor < data sufile=sweep |\t\t\t\t\t\t",
    " suwind itmin=nsweep itmax=nsweep+nsout | sushw key=delrt a=0.0 |...   ",
    " \t\t\t\t\t\t\t\t\t",
    " sweep=vibroseis sweep in SU format, nsweep=number of samples on\t",
    " the vibroseis sweep, nsout=desired number of samples on output\t",
    " \t\t\t\t\t\t\t\t\t",
    " or\t\t\t\t\t\t\t\t\t",
    " suxcor < data sufile=sweep |\t\t\t\t\t\t",
    " suwind itmin=nsweep itmax=nsweep+nsout | sushw key=delrt a=0.0 |...   ",
    " \t\t\t\t\t\t\t\t\t",
    " tsweep=sweep length in seconds, tout=desired output trace length in seconds",
    " \t\t\t\t\t\t\t\t\t",
    " In the spatially variant case (panel=1), a window with linear slope \t",
    " can be defined:\t\t\t\t\t\t \t",
    " \tftwin is the first sample of the first trace in the gather,  \t",
    " \tltwin is the first sample of the last trace in the gather,\t",
    " \tntwin is the lengthe of the window, \t\t\t\t",
    " \tntrc is the the number of traces in a gather. \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \tIf the data consists of a number gathers which need to be \t",
    "\tcorrelated with the same number gathers in the sufile, ntrc\t",
    "\tassures that the correlating window re-starts for each gather.\t",
    " \t\t\t\t\t\t\t\t\t",
    "\tThe default window is non-sloping and takes the entire trace\t",
    "\tinto account (ftwin=ltwin=0, ntwin=nt).\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    let mut sutrace = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut intrace) == 0 { err!("can't get first trace"); }
    let mut nt = intrace.ns as i32;
    let mut vibroseis = 0;
    let mut first = 1;
    let mut panel = 0;
    let mut ntwin = nt;
    let mut ftwin = 0;
    let mut ltwin = 0;
    let mut ntrc = 48;
    getparint("vibroseis", &mut vibroseis);
    getparint("first", &mut first);
    getparint("panel", &mut panel);
    getparint("ntwin", &mut ntwin);
    getparint("ftwin", &mut ftwin);
    getparint("ltwin", &mut ltwin);
    getparint("ntrc", &mut ntrc);
    let mut sufile = String::new();
    let mut nfilter = 0;
    let mut filter = Vec::new();
    let mut fp = None;
    if getparstring("sufile", &mut sufile) == 0 {
        nfilter = countparval("filter");
        if nfilter == 0 { err!("must specify filter= desired filter"); }
        filter = vec![0.0f32; nfilter as usize];
        getparfloat("filter", &mut filter);
    } else {
        let file = efopen(&sufile, "r");
        if panel == 0 {
            fgettr(&file, &mut sutrace);
            nfilter = sutrace.ns as i32;
            filter = sutrace.data[..nfilter as usize].to_vec();
        }
        fp = Some(file);
    }
    checkpars();
    let mut dtwin = 0;
    if panel == 1 {
        if fp.is_none() { err!("sufile= is required when panel=1"); }
        if ntwin < 0 { err!("ntwin must not be negative"); }
        nfilter = ntwin;
        nt = ntwin;
        // ARM64 integer division by zero yields zero in the original build.
        dtwin = ltwin.wrapping_sub(ftwin).checked_div(ntrc).unwrap_or(0);
    } else { ntwin = nt; }
    let ntout = if vibroseis <= 0 { nt as i64 + nfilter as i64 - 1 } else { vibroseis as i64 };
    if ntout < 0 || ntout > SU_NFLTS as i64 {
        err!("output trace length is outside the supported range");
    }
    if panel == 1 {
        if ntwin as usize > SU_NFLTS { err!("ntwin exceeds the trace buffer"); }
        filter = vec![0.0; nfilter as usize];
    }
    let mut trace = vec![0.0f32; nt as usize];
    let mut tracl = 1i32;
    loop {
        if panel == 1 {
            let file = match fp.as_ref() {
                Some(file) => file,
                None => err!("sufile= is required when panel=1"),
            };
            if fgettr(file, &mut sutrace) == 0 {
                err!("number of traces in sufile are not the same as the number of traces on the input");
            }
            // ARM64 remainder with a zero divisor leaves the dividend.
            let mut tracr = tracl.checked_rem(ntrc).unwrap_or(tracl);
            if tracr == 0 { tracr = ntrc; }
            let fs = ftwin.wrapping_add(tracr.wrapping_mul(dtwin));
            if ntwin > 0 && (fs < 0 || fs as i64 + ntwin as i64 > SU_NFLTS as i64) {
                err!("panel window is outside the trace buffer");
            }
            for it in 0..ntwin {
                let it2 = (it + fs) as usize;
                trace[it as usize] = intrace.data[it2];
                filter[it as usize] = sutrace.data[it2];
            }
        } else {
            trace.copy_from_slice(&intrace.data[..nt as usize]);
        }
        let ifz = if vibroseis <= 0 { if first != 0 { -nfilter + 1 } else { -nt + 1 } } else { 0 };
        if first != 0 {
            xcor(0, &filter, 0, &trace, ifz, &mut outtrace.data[..ntout as usize]);
        } else {
            xcor(0, &trace, 0, &filter, ifz, &mut outtrace.data[..ntout as usize]);
        }
        outtrace.hdr = intrace.hdr.clone();
        if vibroseis > 0 { outtrace.delrt = 0; }
        outtrace.ns = ntout as u16;
        puttr(&mut outtrace);
        tracl += 1;
        if gettr(&mut intrace) == 0 { break; }
    }
    exit(CWP_Exit());
}
