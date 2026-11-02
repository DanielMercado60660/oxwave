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

//! SUADDSTATICS - add random statics. Port of su/main/statics/suaddstatics.c (Revision 1.6).
//! Credits: CWP Wences Gouveia, 11/07/94, Colorado School of Mines.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_f_franuni::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUADDSTATICS - ADD random STATICS on seismic data\t\t\t",
    "\t\t\t\t\t\t\t\t \t",
    " suaddstatics required parameters [optional parameters] > stdout\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " shift=\t\tthe static shift will be generated \t \t",
    "\t\t\trandomly in the interval [+shift,-shif] (ms)\t",
    " sources=\t\tnumber of source locations\t\t\t",
    " receivers=\t\tnumber of receiver locations\t\t\t",
    " cmps=\t\t\tnumber of common mid point locations\t\t",
    " maxfold=\t\tmaximum fold of input data\t\t\t",
    " datafile=\t\tname and COMPLETE path of the input file\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " dt=tr.dt\t\t\ttime sampling interval (ms)\t\t",
    " seed=getpid()\t\t seed for random number generator\t\t",
    " verbose=0\t\t\t=1 print useful information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input data should be sorted into cdp gathers.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " SUADDSTATICS applies static time shifts in a surface consistent way on",
    " seismic data sets. SUADDSTATICS writes the static time shifts in the  ",
    " header field TSTAT. To perform the actual shifts the user should use \t",
    " the program SUSTATIC after SUADDSTATICS. SUADDSTATICS outputs the\t",
    " corrupted data set to stdout.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Header field used by SUADDSTATICS: cdp, sx, gx, tstat, dt.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut nsources, mut nreceivers, mut ncmp, mut maxfold) = (0, 0, 0, 0);
    let mut datafile = String::new();
    let (mut seed, mut verbose) = (std::process::id() as i32, 0);
    let (mut dt, mut shift) = (0.0f32, 0.0f32);
    if getparint("sources", &mut nsources) == 0 {
        err!("Specify number of sources!");
    }
    if getparint("receivers", &mut nreceivers) == 0 {
        err!("Specify number of receivers!");
    }
    if getparint("cmps", &mut ncmp) == 0 {
        err!("Specify number of cmps!");
    }
    if getparint("maxfold", &mut maxfold) == 0 {
        err!("Specify maximum fold of the data!");
    }
    if getparstring("datafile", &mut datafile) == 0 {
        err!("Specify data file name!");
    }
    getparint("seed", &mut seed);
    getparint("verbose", &mut verbose);
    getparfloat("dt", &mut dt);
    if getparfloat("shift", &mut shift) == 0 {
        err!("Specify time shift (sec)!");
    }
    if verbose != 0 {
        warn!("Reading input data set");
    }
    checkpars();
    let _first_fp = fopen(&datafile, "r");
    let fp = match fopen(&datafile, "r") {
        Some(fp) => fp,
        None => err!("SUADDSTATICS: can't open input data file"),
    };
    if fgettr(&fp, &mut tr) == 0 {
        err!("can't get first trace");
    }
    let nsamples = tr.ns as usize;
    if dt == 0.0 {
        dt = (tr.dt as f64 / 1000.0) as f32;
    }
    erewind(&fp);
    let mut pointer = vec![vec![[0i32; 2]; maxfold as usize]; ncmp as usize];
    let mut vect = vec![0.0f32; (nsources + nreceivers) as usize];
    let mut source = vec![0i32; nsources as usize];
    let mut receiver = vec![0i32; nreceivers as usize];
    let mut foldcmp = vec![0i32; ncmp as usize];
    let mut data = vec![vec![vec![0.0f32; nsamples]; maxfold as usize]; ncmp as usize];

    // Save native headers separately, as in the C temporary header file.
    let hfp = etmpfile();
    let (mut isource, mut ireceiver) = (0usize, 0usize);
    let (mut itrace, mut ntraces) = (0usize, 0i32);
    let mut icmp = -1i32;
    let mut cmp_last = 0;
    while fgettr(&fp, &mut tr) != 0 {
        // The first gather must start even for CDP zero (C otherwise uses pointer[-1]).
        if icmp < 0 || tr.cdp != cmp_last {
            if icmp >= 0 {
                foldcmp[icmp as usize] = itrace as i32;
            }
            icmp += 1;
            itrace = 0;
            cmp_last = tr.cdp;
        }
        // Match the minimal C bounds fixes instead of panicking on undersized input geometry.
        if icmp >= ncmp {
            err!("number of CMP gathers exceeds cmps=%d", ncmp);
        }
        if itrace as i32 >= maxfold {
            err!("CMP fold exceeds maxfold=%d", maxfold);
        }
        let c = icmp as usize;
        pointer[c][itrace] = [tr.sx, tr.gx];
        if c == 0 && itrace == 0 {
            source[isource] = tr.sx;
        } else if source[isource] != tr.sx {
            let mut i = isource as i32 - 1;
            while i >= 0 && source[i as usize] != tr.sx {
                i -= 1;
            }
            if i < 0 {
                isource += 1;
                if isource >= nsources as usize {
                    break;
                }
                source[isource] = tr.sx;
            }
        }
        if c == 0 && itrace == 0 {
            receiver[ireceiver] = tr.gx;
        } else if receiver[ireceiver] != tr.gx {
            let mut i = ireceiver as i32 - 1;
            while i >= 0 && receiver[i as usize] != tr.gx {
                i -= 1;
            }
            if i < 0 {
                ireceiver += 1;
                if ireceiver >= nreceivers as usize {
                    break;
                }
                receiver[ireceiver] = tr.gx;
            }
        }
        efwrite(std::slice::from_ref(&tr.hdr), &hfp);
        data[c][itrace].copy_from_slice(&tr.data[..nsamples]);
        itrace += 1;
        ntraces += 1;
    }
    foldcmp[icmp as usize] = itrace as i32;
    erewind(&hfp);
    if verbose != 0 {
        fprintf!(
            stderr(),
            "%d TRACES were read in %d CMPs\n",
            ntraces,
            icmp + 1
        );
    }
    for c in 0..ncmp as usize {
        for t in 0..foldcmp[c] as usize {
            let mut s = 0;
            while source[s] != pointer[c][t][0] {
                s += 1;
            }
            pointer[c][t][0] = s as i32;
            let mut r = 0;
            while receiver[r] != pointer[c][t][1] {
                r += 1;
            }
            pointer[c][t][1] = r as i32;
        }
    }
    drop(source);
    drop(receiver);

    sranuni(seed);
    let min_stat = (-shift / dt) as f64;
    let max_stat = (shift / dt) as f64;
    for i in 0..(nsources + nreceivers) as usize {
        vect[i] = randdouble(min_stat, max_stat) as f32;
    }
    for c in 0..ncmp as usize {
        if verbose != 0 {
            fprintf!(stderr(), "Applying fake statics to CMP# %d\n", c as i32);
        }
        for t in 0..foldcmp[c] as usize {
            efread(std::slice::from_mut(&mut tr.hdr), &hfp);
            let source_statics = vect[pointer[c][t][0] as usize] * dt;
            let receiver_statics = vect[(pointer[c][t][1] + nsources) as usize] * dt;
            let trace_static = source_statics + receiver_statics;
            tr.tstat = nint(trace_static as f64) as i16;
            tr.data[..nsamples].copy_from_slice(&data[c][t]);
            puttr(&mut tr);
        }
    }
    exit(CWP_Exit());
}

// Randdouble is local to suaddstatics.c. Its arithmetic promotes the f32
// random draw to f64; main narrows its result back to f32 for vect.
fn randdouble(dlow: f64, dhigh: f64) -> f64 {
    let random = franuni();
    let drandom = random as f64 * (dhigh - dlow) + dlow;
    drandom
}
