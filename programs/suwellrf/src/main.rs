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

//! SUWELLRF - convert well logs to a time-sampled reflectivity function.
//! Port of su/main/well_logs/suwellrf.c (Revision 1.14).
//! Author: CWP: John Stockwell, Summer 2001, updated Summer 2002;
//! inspired by a project by GP grad student Leo Brown.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_sinc_interp::ints8r;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUWELLRF - convert WELL log depth, velocity, density data into a\t",
    "\tuniformly sampled normal incidence Reflectivity Function of time",
    " \t\t\t\t\t\t\t\t\t",
    " suwellrf [required parameters] [optional parameters] > [stdout]\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " dvrfile=\tfile containing depth, velocity, and density values\t",
    " ...or...\t\t\t\t\t\t\t\t",
    " dvfile=\tfile containing depth and velocity values\t\t",
    " drfile=\tfile containing depth and density values\t\t",
    " ...or...\t\t\t\t\t\t\t\t",
    " dfile=\tfile containing depth values\t\t\t\t",
    " vfile=\tfile containing velocity log values\t\t\t",
    " rhofile=\tfile containing density log values\t\t\t",
    " nval= \tnumber of triplets of d,v,r values if dvrfile is set,\t",
    " \t\tnumber of pairs of d,v and d,r values dvfile and drfile\t",
    "\t\tare set, or number of values if dfile, vfile, and rhofile",
    "\t\tare set.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dtout=.004\tdesired time sampling interval (sec) in output\t\t",
    " ntr=1         number of traces to output \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The format of the input file(s) is C-style binary float. These files\t",
    " may be constructed from ascii file via:   \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "       a2b n1=3 < dvrfile.ascii > dvrfile.bin\t\t\t\t",
    " ...or...\t\t\t\t\t\t\t\t",
    "       a2b n1=2 < dvfile.ascii > dvfile.bin\t\t\t\t",
    "       a2b n1=2 < drfile.ascii > drfile.bin\t\t\t\t",
    " ...or...\t\t\t\t\t\t\t\t",
    "       a2b n1=1 < dfile.ascii > dfile.bin\t\t\t\t",
    "       a2b n1=1 < vfile.ascii > dfile.bin\t\t\t\t",
    "       a2b n1=1 < rhofile.ascii > rhofile.bin\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " A raw normal-incidence impedence reflectivity as a function of time is",
    " is generated using the smallest two-way traveltime implied by the\t",
    " input velocities as the time sampling interval. This raw reflectivity\t",
    " trace is then resampled to the desired output time sampling interval\t",
    " via 8 point sinc interpolation. If the number of samples on the output",
    " exceeds SU_NFLTS the output trace will be truncated to that value.\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Caveat: \t\t\t\t\t\t\t\t",
    " This program is really only a first rough attempt at creating a well\t",
    " log utility. User input and modifications are welcome.\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " See also:  suresamp \t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn input_array(len: usize) -> Vec<f32> {
    let mut values = Vec::new();
    if values.try_reserve_exact(len).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    values.resize(len, 0.0);
    values
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut nval = 0;
    mustgetparint("nval", &mut nval);
    let mut dfile = String::new();
    let mut vfile = String::new();
    let mut rhofile = String::new();
    let mut dvfile = String::new();
    let mut drfile = String::new();
    let mut dvrfile = String::new();
    getparstring("dfile", &mut dfile);
    getparstring("vfile", &mut vfile);
    getparstring("rhofile", &mut rhofile);
    getparstring("dvfile", &mut dvfile);
    getparstring("drfile", &mut drfile);
    getparstring("dvrfile", &mut dvrfile);
    let mut ntr = 1;
    getparint("ntr", &mut ntr);
    if nval < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    // The original computes these input counts as signed integers.
    let ndvr = nval.checked_mul(3).unwrap_or_else(|| err!("nval is too large")) as usize;
    let n = nval as usize;
    let mut dvrinput = input_array(ndvr);
    if !dvrfile.is_empty() {
        let fp = efopen(&dvrfile, "r");
        efread(&mut dvrinput, &fp);
    } else if !dvfile.is_empty() {
        let npairs = n.checked_mul(2).unwrap_or_else(|| err!("nval is too large"));
        let mut temp1 = input_array(npairs);
        let mut temp2 = input_array(npairs);
        let dvfp = efopen(&dvfile, "r");
        let drfp = efopen(&drfile, "r");
        efread(&mut temp1, &dvfp);
        efread(&mut temp2, &drfp);
        for i in 0..n {
            dvrinput[3*i] = temp1[2*i];
            dvrinput[3*i+1] = temp1[2*i+1];
            dvrinput[3*i+2] = temp2[2*i+1];
        }
    } else if !dfile.is_empty() && !vfile.is_empty() && !rhofile.is_empty() {
        let mut temp1 = input_array(n);
        let mut temp2 = input_array(n);
        let mut temp3 = input_array(n);
        let dfp = efopen(&dfile, "r");
        let vfp = efopen(&vfile, "r");
        let rfp = efopen(&rhofile, "r");
        efread(&mut temp1, &dfp);
        efread(&mut temp2, &vfp);
        efread(&mut temp3, &rfp);
        for i in 0..n {
            dvrinput[3*i] = temp1[i];
            dvrinput[3*i+1] = temp2[i];
            dvrinput[3*i+2] = temp3[i];
        }
    }
    let mut tempd = input_array(n);
    let mut tempv = input_array(n);
    let mut tempr = input_array(n);
    let mut dtmin = f32::MAX;
    let mut tmax = 0.0f32;
    for i in 0..n {
        tempd[i] = dvrinput[3*i];
        tempv[i] = dvrinput[3*i+1];
        tempr[i] = dvrinput[3*i+2];
        if i > 0 {
            if tempd[i] - tempd[i-1] < 0.0 { err!("depth not monotonically increasing!"); }
            let tempdt = (2.0 * (tempd[i] - tempd[i-1]) as f64 / tempv[i-1] as f64) as f32;
            dtmin = if dtmin < tempdt { dtmin } else { tempdt };
            tmax += tempdt;
        }
    }
    let mut dtout = 0.0f32;
    if getparfloat("dtout", &mut dtout) == 0 { dtout = (dtmin as f64 / 2.0) as f32; }
    let mut tr = Trace::new();
    tr.dt = nint(dtout as f64 * 1000000.0) as u16;
    checkpars();
    let ntmax = nint((tmax / dtmin) as f64);
    // The original writes its last spike at index ntmax, outside the input.
    let mut maxtrace = vec![0.0f32; (ntmax.max(0) as usize) + 1];
    let mut ntout = nint((tmax / dtout) as f64);
    if ntout > SU_NFLTS as i32 {
        warn!("ntout=%d too big! Setting ntout to SU_NFLTS=%d", ntout, SU_NFLTS as i32);
        ntout = SU_NFLTS as i32;
    }
    tr.ns = ntout as u16;
    let mut ttotal = 0.0f32;
    for i in 1..n {
        ttotal = (ttotal as f64 + 2.0 * (tempd[i]-tempd[i-1]) as f64 / tempv[i-1] as f64) as f32;
        let itempt = nint(((ntmax as f32 * ttotal) / tmax) as f64) as usize;
        if itempt < maxtrace.len() {
            maxtrace[itempt] = (tempv[i]*tempr[i] - tempv[i-1]*tempr[i-1]) / (tempv[i]*tempr[i] + tempv[i-1]*tempr[i-1]);
        }
    }
    // ealloc1float receives a size_t; a negative count becomes an
    // impossible allocation in the original.
    if ntout < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let mut t = input_array(ntout as usize);
    let mut tvalue = 0.0f32;
    for v in &mut t { *v = tvalue; tvalue += dtout; }
    ints8r(dtmin, 0.0, &maxtrace[..ntmax.max(0) as usize], 0.0, 0.0, &t, &mut tr.data);
    for _ in 0..ntr { puttr(&mut tr); }
    exit(CWP_Exit());
}
