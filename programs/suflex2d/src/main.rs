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

//! SUFLEX2D - Duplicate traces into CDPs from neighbouring CDPs.
//! Port of su/main/headers/suflex2d.c (Revision 1.1).
//! Credits: Andre Latour, Nov 2022.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                            ",
    " SUFLEX2D - Duplicate Some Traces Into CDPs From Neighbouring CDPs.         ",
    "                                                                            ",
    " suflex2d     <in.su >out.su  (other parameters)                            ",
    "                                                                            ",
    " This program assigns input traces of a cdp to offset bins (offset ranges). ",
    " If an output cdp does not have a user-specified number of traces in an     ",
    " offset bin, this program searches the neighbouring cdps and duplicates and ",
    " outputs traces from the same offset bin. The neighbouring traces that are  ",
    " closer to the output cdp centre have higher duplication priority. To know  ",
    " which traces are closer to cdp centre, this program uses igi key values.   ",
    " The igi key values should be the inline distance between trace midpoint    ",
    " coordinates and cdp profile coordinates, as computed by program sunearcsv. ",
    "                                                                            ",
    " Traces duplicated from neighbouring cdps have their cdp key value changed  ",
    " to the output cdp they are duplicated into. (Also their igi key is changed ",
    " as described in the igiadd parameter). Note that being duplicated and      ",
    " output does not mean that a trace is deleted from a neighbouring cdp.      ",
    " A trace may be deleted if option keepmore is 0, but that is unrelated to   ",
    " the duplication.                                                           ",
    "                                                                            ",
    " Input traces must be ordered by cdp and absolute value of offset. This is  ",
    " checked and an error-halt will occur at the first input trace where this   ",
    " is not the case. (See Note On Output for a slight exception to this).      ",
    "                                                                            ",
    " Parameter overview:                                                        ",
    "                                                                            ",
    " maxflex=1    Maximum Relative cdp to include in flexing. Must be equal to  ",
    "              zero or greater. This program uses the input traces of cdps   ",
    "              from minflex to maxflex (inclusive) around each output cdp.   ",
    "                                                                            ",
    " minflex= -maxflex  Minimum Relative cdp to include in flexing. Must be     ",
    "                    equal to zero or smaller. The default is the magnitude  ",
    "                    of maxflex with a negative sign.                        ",
    "                                                                            ",
    " maxfold=100  Maximum number of traces per input cdp gather.                ",
    "              Note that the output cdp gathers usually have more traces     ",
    "              than input cdps (and sometimes more traces than this value).  ",
    "                                                                            ",
    " binsize=25   Offset Bin Size.                                              ",
    "              Note: The offsets are absolute valued before use herein.      ",
    "              Note: Even though the offset key contains integers,           ",
    "                    this parameter allows floating point numbers.           ",
    "              Note: Offset bins contain their lower boundary but their      ",
    "                    higher boundary belongs to the next higher bin.         ",
    "                                                                            ",
    " binalign=binsize/2   Offset Bin Alignment. Default is half of binsize.     ",
    "                      Specify the offset value at the centre of any bin     ",
    "                      that you find convenient. For instance, usually, the  ",
    "                      nominal minimum shot-to-receiver distance produces    ",
    "                      a good alignment of the bin centres.                  ",
    "              Note: Even though the offset key contains integers,           ",
    "                    this parameter allows floating point numbers.           ",
    "              Note: This parameter does not restrict the offset range, the  ",
    "                    offset bins extend from 0 to infinity and all traces    ",
    "                    are put into some offset bin. Empty bins are ignored.   ",
    "              Note: Having an offset distance closer to the bin centre does ",
    "                    not give higher priority to a trace. Traces in the same ",
    "                    bin have priority only based on cdp and igi key values. ",
    "                                                                            ",
    " binfold=2    Desired Number of Traces Per Offset Bin. Typical split-spread ",
    "              surveys have two traces with the same offset in each cdp.     ",
    "                                                                            ",
    " binbest=0    Use IGI key values to determine which traces are better to    ",
    "              duplicate to fullfill binfold requirements. Typically, igi    ",
    "              values come from the sunearcsv program and contain the        ",
    "              signed inline distance (positive values means traces are      ",
    "              closer to next cdp, negative means closer to previous cdp).   ",
    "        =1    Use IGC key values but absolute value them. Typically, igc    ",
    "              values come from the sunearcsv program and contain the        ",
    "              signed crossline distance. But this option absolute values    ",
    "              them before use herein. So it does not matter which side of   ",
    "              profile traces are located, only how far away they are in     ",
    "              the crossline direction (smaller absolute values are better). ",
    "        =2    Use IGC key values but do not absolute value them.            ",
    "              Smaller values (including negatives) are better. Normally     ",
    "              you would not use this option for igc values directly from    ",
    "              sunuearcsv since it means prefering one crossline direction.  ",
    "        =-2   Same as option 2, but reverse the sign of IGC values first.   ",
    "              This makes larger values better.                              ",
    "        =3    Use OFFSET key values. Values closer to the centre of the     ",
    "              offset bin of the trace are better (see binalign parameter).  ",
    "              Typically you might use this option with the subinbigcsv      ",
    "              program (and set binfold=1) to select traces from a super-cdp ",
    "              to get a nice offset distribution for analysis/display/tests. ",
    "        Note: Remember that trace priority only applies when comparing      ",
    "              symmetric cdps around the output cdp. Input traces for each   ",
    "              cdp always have more priority than input traces of any        ",
    "              neighbouring cdps. And a close neighbouring cdp always has    ",
    "              more priority than a far neighbouring cdp. Only SYMMETRIC     ",
    "              neighbouring cdps have their igi or igc values compared.      ",
    "        Note: If igi or igc values are all zero (or any constant), traces   ",
    "              will still be duplicated from the closest neighbouring cdp to ",
    "              fullfill binfold requirements. So you can still use suflex2d  ",
    "              even if you have not used sunearcsv to set igi,igc.           ",
    "        Note: No option here uses igi or igc as true distances, only to     ",
    "              determine the relative priority for duplicating traces.       ",
    "              You can therefore set one of these keys yourself to some      ",
    "              priority that you have computed.                              ",
    "                                                                            ",
    " keepmore=0   Do not keep more than binfold traces per bin. Traces in target",
    "              cdp always have priority over traces of neighbouring cdps.    ",
    "              But sometimes the target cdp already has more input traces in ",
    "              an offset bin than binfold. Option 0 here means that excess   ",
    "              traces with the worst binbest priority values are not output. ",
    "         =1   Output those traces anyway.                                   ",
    "                                                                            ",
    " igiadd=100   Add this to the igi key values for output traces that are     ",
    "              duplicated into a cdp from neighbouring cdps. Multiply this   ",
    "              value by the difference in cdp numbers between current cdp    ",
    "              and whichever neighbouring cdp the trace originated from.     ",
    "              Usually you should use a value here that easily allows you to ",
    "              tell which output traces originated from neighbouring cdps.   ",
    "              But in some cases setting this to the distance between cdps   ",
    "              is a good choice (especially if you are attempting to cascade ",
    "              the flexing - see next note).                                 ",
    "              Note that this value is added to the igi key values even if   ",
    "              the binbest option is 1 or 2.                                 ",
    "                                                                            ",
    " Note On Output:                                                            ",
    "       Traces duplicated from neighbouring CDPs have their cdp key value    ",
    "       changed to the cdp they are duplicated into. (Also their igi key is  ",
    "       changed as described in the igiadd parameter). But duplicated traces ",
    "       are not output in true offset order, just in their bin order. This   ",
    "       means the output must be re-sorted by offset distance if you need    ",
    "       to have true offset order.                                           ",
    "       And, this program has a slight exception to requirement that input   ",
    "       be ordered by cdp and absolute value of offset. The traces only must ",
    "       be ordered by cdp and absolute value of offset bins. Within bins the ",
    "       offsets can be out-of-order. Thus output from this program can be    ",
    "       re-input without sorting (if binsize and binalign remain the same).  ",
    "           (Re-input allows cascading of the flexing, but requires          ",
    "            extensive understanding of the situation).                      ",
    " ------------------------------------------------------------------------   ",
    "                                                                            ",
];

struct Gather {
    cdp: i32,
    traces: Vec<Trace>,
    bins: Vec<i32>,
    best: Vec<i32>,
}

fn flexit(
    gs: &[Gather],
    target: i32,
    min: i32,
    max: i32,
    fold: i32,
    best: i32,
    keep: i32,
    add: i32,
) -> i32 {
    let mut lo = 2147483645;
    let mut hi = -2147483645;
    let mut used = Vec::new();
    let mut bins = Vec::new();
    let mut blo = vec![0usize; gs.len()];
    let mut bhi = blo.clone();
    for g in gs {
        if !g.traces.is_empty() {
            if lo > g.bins[0] {
                lo = g.bins[0];
            }
            if hi < *g.bins.last().unwrap() {
                hi = *g.bins.last().unwrap();
            }
        }
        used.push(vec![false; g.traces.len()]);
        let mut b = g.bins.clone();
        b.push(2147483644);
        b.push(2147483645);
        bins.push(b);
    }
    let mut out = 0;
    let mut bin = lo;
    while bin <= hi {
        let mut n = 0;
        for (i, g) in gs.iter().enumerate() {
            if target >= g.cdp.wrapping_add(min) && target <= g.cdp.wrapping_add(max) {
                while bins[i][blo[i]] < bin {
                    blo[i] += 1;
                }
                bhi[i] = blo[i] + 1;
                while bins[i][bhi[i]] == bins[i][blo[i]] {
                    bhi[i] += 1;
                }
            }
        }
        loop {
            let mut choice = None;
            let mut priority = 2147483645i32;
            for (i, g) in gs.iter().enumerate() {
                if target >= g.cdp.wrapping_add(min)
                    && target <= g.cdp.wrapping_add(max)
                    && bins[i][blo[i]] == bin
                {
                    for k in blo[i]..bhi[i] {
                        let d = g.cdp.wrapping_sub(target);
                        let p = if best == 0 {
                            g.best[k]
                                .wrapping_add(d.wrapping_mul(100000))
                                .wrapping_abs()
                        } else {
                            g.best[k].wrapping_add(d.wrapping_abs().wrapping_mul(100000))
                        };
                        if !used[i][k] && p < priority {
                            priority = p;
                            choice = Some((i, k));
                        }
                    }
                }
            }
            let Some((i, k)) = choice else {
                break;
            };
            if keep == 1 && n >= fold && priority > 50000 {
                break;
            }
            let g = &gs[i];
            let mut tr = g.traces[k].clone();
            if g.cdp != target {
                tr.cdp = target;
                tr.igi = (i32::from(tr.igi)
                    .wrapping_add(g.cdp.wrapping_sub(target).wrapping_mul(add)))
                    as i16;
            }
            puttr(&mut tr);
            used[i][k] = true;
            n += 1;
            out += 1;
            if keep == 0 && n == fold {
                break;
            }
        }
        bin = bin.wrapping_add(1);
    }
    out
}

fn bin_and_best(tr: &Trace, size: f64, align: f64, best: i32) -> (i32, i32) {
    let bin = (0.5000001 + (tr.offset.wrapping_abs() as f64 - align) / size) as i32;
    let p = match best {
        0 => tr.igi as i32,
        1 => (tr.igc as i32).abs(),
        2 => tr.igc as i32,
        -2 => -(tr.igc as i32),
        _ => ((bin as f64 * size + align) as i32)
            .wrapping_sub(tr.offset)
            .wrapping_abs(),
    };
    (bin, p)
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut maxflex = 1;
    getparint("maxflex", &mut maxflex);
    if maxflex < 0 {
        err!("error: maxflex must be equal to or greater than 0.");
    }
    let mut minflex = -maxflex;
    getparint("minflex", &mut minflex);
    if minflex > 0 {
        err!("error: minflex must be equal to or less than 0.");
    }
    let numgathers = (maxflex - minflex + 1) as usize;
    let mut maxfold = 100;
    getparint("maxfold", &mut maxfold);
    if maxfold < 1 {
        err!("error: maxfold must be greater than 0.");
    }
    let mut binfold = 2;
    getparint("binfold", &mut binfold);
    if binfold < 1 {
        err!("error: binfold must be greater than 0.");
    }
    let mut binbest = 0;
    getparint("binbest", &mut binbest);
    if binbest != -2 && binbest != 0 && binbest != 1 && binbest != 2 && binbest != 3 {
        err!("error: binbest parameter is not one of the options.");
    }
    let mut binsize = 25.0;
    getpardouble("binsize", &mut binsize);
    if binsize <= 0.0 {
        err!("error: binsize must be greater than 0.0");
    }
    let mut binalign = binsize / 2.0;
    getpardouble("binalign", &mut binalign);
    if binalign < 0.0 {
        err!("error: binalign must be equal to or greater than 0.0");
    }
    binalign %= binsize;
    let mut keepmore = 0;
    getparint("keepmore", &mut keepmore);
    let mut igiadd = 100;
    getparint("igiadd", &mut igiadd);
    let mut gs: Vec<Gather> = (0..numgathers)
        .map(|_| Gather {
            cdp: 0,
            traces: Vec::new(),
            bins: Vec::new(),
            best: Vec::new(),
        })
        .collect();
    checkpars();
    let mut tr = Trace::new();
    let mut nproct = 0;
    let mut mproct = 0;
    let mut iroll = 0usize;
    let mut target = 0i32;
    let mut maxstored = 0i32;
    loop {
        let has = gettr(&mut tr) != 0;
        let mut newcdp = !has;
        let iend = if has { 0 } else { maxflex * 2 };
        if has {
            if nproct == 0 {
                iroll = 0;
                target = tr.cdp + minflex - 1;
                for (i, g) in gs.iter_mut().enumerate() {
                    g.cdp = tr.cdp + minflex - 10 * (i as i32 + 1);
                }
                maxstored = gs[0].cdp;
            }
            nproct += 1;
            if gs[iroll].cdp == tr.cdp {
                let g = &mut gs[iroll];
                if g.traces.len() == maxfold as usize {
                    err!("error: At input trace count %d there are more than maxfold (%d) traces in cdp %d.", nproct, maxfold, tr.cdp);
                }
                let (bin, p) = bin_and_best(&tr, binsize, binalign, binbest);
                if bin < *g.bins.last().unwrap() {
                    err!("error: At input trace count %d (cdp %d) trace is not in offset (bin) order.", nproct, tr.cdp);
                }
                g.traces.push(tr.clone());
                g.bins.push(bin);
                g.best.push(p);
            } else {
                newcdp = true;
            }
        }
        if newcdp {
            while target < maxstored + 1 - maxflex + iend {
                mproct += flexit(
                    &gs, target, minflex, maxflex, binfold, binbest, keepmore, igiadd,
                );
                target += 1;
            }
            if gs[iroll].cdp > tr.cdp {
                err!(
                    "error: At input trace count %d (cdp %d) is not in cdp order",
                    nproct,
                    tr.cdp
                );
            }
            iroll = (iroll + 1) % numgathers;
            let g = &mut gs[iroll];
            g.cdp = tr.cdp;
            g.traces.clear();
            g.bins.clear();
            g.best.clear();
            let (bin, p) = bin_and_best(&tr, binsize, binalign, binbest);
            g.traces.push(tr.clone());
            g.bins.push(bin);
            g.best.push(p);
            maxstored = tr.cdp;
        }
        if !has {
            break;
        }
    }
    warn!(
        "Number of traces input=%d.  Number output=%d.",
        nproct, mproct
    );
    exit(CWP_Exit());
}
