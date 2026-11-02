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
//! SUCVS4FOWLER - constant velocity stacks for Fowler codes.
//! Port of su/main/stacking/sucvs4fowler.c (Revision 1.5).
//! Credits: John E. Anderson, Spring 1994.

use oxwave_core::prelude::*;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUCVS4FOWLER --compute constant velocity stacks for input to Fowler codes",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameter:\t\t\t\t\t\t\t",
    " ncdps=\t\tnumber of input cdp gathers\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " vminstack=1500.\tminimum velocity panel in m/s to output\t\t",
    " nvstack=180\t\tnumber of stacking velocity panels to compute\t",
    "\t\t\t( Let offmax be the maximum offset, fmax be\t",
    "\t\t\tthe maximum freq to preserve, and tmute be\t",
    "\t\t\tthe starting mute time in sec on offmax, then\t",
    "\t\t\tthe recommended value for nvstack would be\t",
    "\t\t\tnvstack = 4 +(offmax*offmax*fmax)/(0.6*vmin*vmin*tmute)",
    "\t\t\t---you may want to make do with less---)\t\t",
    " lmute=24\t\tlength of mute taper in ms\t\t\t",
    " nonhyp=1\t\t1 if do mute at 2*offset/vhyp to avoid\t\t",
    "\t\t\tnon-hyperbolic moveout, 0 otherwise\t\t",
    " vhyp=2500.\t\tvelocity to use for non-hyperbolic moveout mute\t",
    " lbtaper=0\t\tlength of bottom taper in ms\t\t\t",
    " lstaper=0\t\tlength of side taper in traces\t\t\t",
    " dtout=1.5*dt\t\toutput sample rate in s,\t\t\t",
    "\t\t\tnote: typically fmax=salias*0.5/dtout\t\t",
    " mxfold=120\t\tmaximum number of offsets/input cmp\t\t",
    " salias=0.8\t\tfraction of output frequencies to force within sloth",
    "\t\t\tantialias limit.  This controls muting by offset of",
    "\t\t\tthe input data prior to computing the cv stacks\t",
    "\t\t\tfor values of choose=1 or choose=2.\t\t",
    " Required trace header words on input are ns, dt, cdp, offset.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn taper(lxtaper: i32, lbtaper: i32, nx: i32, ix: i32, nt: usize, trace: &mut [f32]) {
    let x = if ix < lxtaper {
        (0.54 + 0.46 * (std::f64::consts::PI * (lxtaper - ix) as f64 / lxtaper as f64).cos()) as f32
    } else if ix >= nx - lxtaper {
        (0.54
            + 0.46 * (std::f64::consts::PI * (lxtaper + ix + 1 - nx) as f64 / lxtaper as f64).cos())
            as f32
    } else {
        1.0
    };
    if x != 1.0 {
        for t in &mut trace[..nt] {
            *t *= x;
        }
    }
    for it in 0.max(nt as i32 - lbtaper) as usize..nt {
        trace[it] = (trace[it] as f64
            * (0.54
                + 0.46
                    * (std::f64::consts::PI * (lbtaper + it as i32 + 1 - nt as i32) as f64
                        / lbtaper as f64)
                        .cos())) as f32;
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    let mut trout = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("Can't get first trace \n");
    }
    let nt = tr.ns as usize;
    let dt = (0.000001 * tr.dt as f64) as f32;
    let mut oldcmp = tr.cdp;
    let mut salias = 0.8f32;
    getparfloat("salias", &mut salias);
    if salias > 1.0 {
        salias = 1.0;
    }
    let mut dtout = (1.5 * dt as f64) as f32;
    getparfloat("dtout", &mut dtout);
    let ntout = (1.0 + (nt as f32 * dt / dtout) as f64) as i32;
    let fmax = (salias as f64 * 0.5 / dtout as f64) as f32;
    let mut vmin = 1500f32;
    getparfloat("vminstack", &mut vmin);
    let mut vhyp = 2500f32;
    getparfloat("vhyp", &mut vhyp);
    let mut nvstack = 180i32;
    getparint("nvstack", &mut nvstack);
    let mut lmute = 24i32;
    getparint("lmute", &mut lmute);
    lmute = (1.0 + 0.001 * lmute as f64 / dtout as f64) as i32;
    let mut ncmps = 0i32;
    if getparint("ncdps", &mut ncmps) == 0 {
        err!("must enter ncdps");
    }
    let mut lbtaper = 0i32;
    getparint("lbtaper", &mut lbtaper);
    let mut lstaper = 0i32;
    getparint("lstaper", &mut lstaper);
    let mut mxfold = 120i32;
    getparint("mxfold", &mut mxfold);
    let mut nonhyp = 1i32;
    getparint("nonhyp", &mut nonhyp);
    lbtaper = (lbtaper as f64 / (1000.0 * dt as f64)) as i32;
    checkpars();
    let nout = ntout.max(0) as usize;
    let mut p2stack = vec![0f32; nvstack.max(0) as usize];
    // The original writes one slot before checking for excess fold.
    let slots = mxfold.max(0) as usize + 1;
    let mut mute = vec![0f32; slots];
    let mut off = vec![0f32; slots];
    let mut rt = vec![0f32; nout];
    let mut ttn = vec![0f32; nout];
    let mut qtn = vec![0f32; nout];
    let mut fold = vec![0f32; nout];
    let mut gath = vec![0f32; nt * slots];
    let dp2 = (1.0 / (vmin * vmin * (nvstack - 5) as f32) as f64) as f32;
    for iv in 0..p2stack.len() {
        p2stack[iv] = iv as f32 * dp2;
    }
    let mut icmp = 0i32;
    let mut noff = 0usize;
    let mut ieod = false;
    loop {
        if tr.cdp != oldcmp || ieod {
            for iv in 0..p2stack.len() {
                rt.fill(0.0);
                fold.fill(0.0);
                for ioff in 0..noff {
                    let factor = off[ioff] * off[ioff] * p2stack[iv];
                    let itmute = (mute[ioff] / dtout) as i32;
                    let start = itmute.max(0).min(ntout.max(0)) as usize;
                    for it in start..nout {
                        let t = it as f32 * dtout;
                        ttn[it] = ((t * t + factor) as f64).sqrt() as f32;
                    }
                    ints8r(
                        dt,
                        0.0,
                        &gath[ioff * nt..(ioff + 1) * nt],
                        0.0,
                        0.0,
                        &ttn[start..],
                        &mut qtn[start..],
                    );
                    for it in start..(itmute + lmute).max(0).min(ntout.max(0)) as usize {
                        qtn[it] *= (it as i32 - itmute + 1) as f32 / lmute as f32;
                    }
                    for it in start..nout {
                        rt[it] += qtn[it];
                    }
                    for it in start..nout {
                        fold[it] += 1.0;
                    }
                }
                for it in 0..nout {
                    if fold[it] == 0.0 {
                        fold[it] = 1.0;
                    }
                }
                for it in 0..nout {
                    trout.data[it] = rt[it] / fold[it];
                }
                trout.ns = ntout as u16;
                trout.dt = (1000000f32 * dtout) as u16;
                trout.cdp = oldcmp;
                trout.tracf = iv as i32;
                puttr(&mut trout);
            }
            if icmp == 20 * (icmp / 20) {
                fprintf!(
                    stderr(),
                    "cvstacks: completed stacking velocity analysis for cmp %d\n",
                    icmp
                );
            }
            icmp += 1;
            if icmp >= ncmps {
                ieod = true;
            }
            oldcmp = tr.cdp;
            noff = 0;
        }
        if lbtaper > 0 || lstaper > 0 {
            taper(lstaper, lbtaper, ncmps, icmp, nt, &mut tr.data);
        }
        gath[noff * nt..(noff + 1) * nt].copy_from_slice(&tr.data[..nt]);
        off[noff] = tr.offset as f32;
        mute[noff] = fmax * off[noff] * off[noff] * dp2;
        if nonhyp != 0 {
            let v = 2f32 * off[noff] / vhyp;
            if mute[noff] <= v {
                mute[noff] = v;
            }
        }
        noff += 1;
        if noff as i32 > mxfold {
            err!("cvstacks: input cdp has more traces than mxfold");
        }
        if ieod {
            break;
        }
        if gettr(&mut tr) == 0 {
            ieod = true;
        }
    }
    fprintf!(
        stderr(),
        "total VND memory at end = %ld\n",
        // VNDemalloc rounds the unfreed gather allocation to sizeof(long).
        (nt as i64 * mxfold as i64 * 4 + 7) / 8 * 8 + 24
    );
    fprintf!(
        stderr(),
        "cvstacks: read and stacked %d cmp gathers\n",
        icmp
    );
    exit(CWP_Exit());
}
