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
//! SUVELAN_UCCS - velocity analysis using unnormalized crosscorrelation.
//! Port of su/main/velocity_analysis/suvelan_uccs.c (Revision 1.7).
//! Credits: CWP: Valmore Celis, Sept 2002; based on suvelan.c by
//! Colorado School of Mines: Dave Hale c. 1989. Reference: Neidell and
//! Taner, Geophysics, 36, 498-509 (1971).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUVELAN_UCCS - compute stacking VELocity panel for cdp gathers\t     ",
    "\t\tusing UnNormalized CrossCorrelation Sum \t             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " suvelan_uccs <stdin >stdout [optional parameters]\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:\t\t\t\t\t\t\t     ",
    " nx=tr.cdpt              number of traces in cdp \t\t\t     ",
    " nv=50                   number of velocities\t\t\t\t     ",
    " dv=50.0                velocity sampling interval\t\t\t     ",
    " fv=1500.0               first velocity\t\t\t\t     ",
    " smute=1.5               samples with NMO stretch exceeding smute are zeroed",
    " dtratio=5               ratio of output to input time sampling intervals   ",
    " nsmooth=dtratio*2+1     length of smoothing window                         ",
    " verbose=0               =1 for diagnostic print on stderr\t\t     ",
    " pwr=1.0                 semblance value to the power      \t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "Notes:\t\t\t\t\t\t\t\t\t     ",
    " Unnormalized crosscorrelation sum: sum all possible crosscorrelation trace ",
    " pairs in a CMP gather for each trial velocity and zero-offset two-way      ",
    " travel time inside a time window. This unnormalized coherency measure      ",
    " produces large spectral amplitudes for strong reflections and small        ",
    " spectral amplitudes for weaker ones. If M is the number of traces in the   ",
    " CMP gather M(M-1)/2 is the total number of crosscorrelations for each trial",
    " velocity and zero-offset two-way traveltime.\t\t\t \t     ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut cdp = tr.cdp;
    let mut nx = tr.cdpt;
    let mut nv = 50;
    let mut dv = 100.0f32;
    let mut fv = 1500.0f32;
    let mut smute = 1.5f32;
    let mut dtratio = 5;
    getparint("nx", &mut nx);
    getparint("nv", &mut nv);
    getparfloat("dv", &mut dv);
    getparfloat("fv", &mut fv);
    getparfloat("smute", &mut smute);
    if smute <= 1.0 {
        err!("smute must be greater than 1.0");
    }
    getparint("dtratio", &mut dtratio);
    let mut nsmooth = dtratio.wrapping_mul(2).wrapping_add(1);
    getparint("nsmooth", &mut nsmooth);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut pwr = 1.0f32;
    getparfloat("pwr", &mut pwr);
    if pwr < 0.0 {
        err!("we are not looking for noise: pwr < 0");
    }
    if pwr == 0.0 {
        err!("we are creating an all-white semblance: pwr = 0");
    }
    checkpars();
    // These sizes lead to invalid allocations or indexing in the C program.
    if nx <= 0 {
        err!("nx must be greater than zero");
    }
    if dtratio <= 0 {
        err!("dtratio must be greater than zero");
    }
    if nsmooth < 0 {
        err!("nsmooth must not be negative");
    }
    // C's arm64 integer division by zero returns zero.
    let ntout = 1 + (nt - 1).checked_div(dtratio).unwrap_or(0);
    CHECK_NT("ntout", ntout);
    let dtout = dt * dtratio as f32;
    if verbose != 0 {
        fprintf!(stderr(), "\tnumber of output time samples is %d\n", ntout);
        fprintf!(stderr(), "\toutput time sampling interval is %g\n", dtout);
        fprintf!(stderr(), "\toutput time of first sample is %g\n", ft);
    }
    let mut aoff = vec![0.0f32; nx as usize];
    let mut datm = vec![vec![0.0f32; nt as usize]; nx as usize];
    let mut datn = datm.clone();
    let mut sem = vec![0.0f32; ntout as usize];
    let mut gottrace = true;
    let mut cdpprev = tr.cdp;
    let mut ix = 0usize;
    loop {
        if gottrace {
            cdp = tr.cdp;
        }
        if cdp != cdpprev || !gottrace {
            // Output uses the current trace as its header and overwrites its
            // samples. Keep the first trace of the next gather intact.
            let incoming = if gottrace { Some(tr.clone()) } else { None };
            let mut v = fv;
            for _iv in 0..nv {
                for ix in 0..nx as usize {
                    let offovs = (aoff[ix] * aoff[ix]) / (v * v);
                    let tnmute = ((offovs as f64 / ((smute * smute) as f64 - 1.0)).sqrt()) as f32;
                    let itmute = if tnmute > ft {
                        ((tnmute - ft) / dt) as i32
                    } else {
                        0
                    };
                    if itmute < 0 {
                        err!("invalid NMO mute index");
                    }
                    let mut tn = ft + itmute as f32 * dt;
                    for it in itmute..nt {
                        let ti =
                            ((((tn * tn + offovs) as f64).sqrt() - ft as f64) / dt as f64) as f32;
                        let iti = ti as i32;
                        if iti < nt - 1 {
                            if iti < 0 {
                                err!("invalid NMO interpolation index");
                            }
                            let frac = ti - iti as f32;
                            datn[ix][it as usize] = ((1.0 - frac as f64)
                                * datm[ix][iti as usize] as f64
                                + (frac * datm[ix][(iti + 1) as usize]) as f64)
                                as f32;
                        }
                        tn += dt;
                    }
                }
                for itout in 0..ntout {
                    let it = itout * dtratio;
                    let mut ismin = it - nsmooth / 2;
                    let mut ismax = it + nsmooth / 2;
                    if ismin < 0 {
                        ismin = 0;
                    }
                    if ismax > nt - 1 {
                        ismax = nt - 1;
                    }
                    let mut ff = 0.0f32;
                    for is in ismin..ismax {
                        let mut dsum = 0.0f32;
                        let mut fac = 0.0f32;
                        for ix in 0..nx as usize {
                            let nsum = datn[ix][is as usize];
                            dsum += nsum;
                            fac += nsum * dsum;
                        }
                        ff += fac;
                    }
                    sem[itout as usize] = ff;
                }
                tr.offset = 0;
                tr.cdp = cdpprev;
                tr.ns = ntout as u16;
                tr.dt = (dtout as f64 * 1000000.0) as u16;
                tr.data[..ntout as usize].copy_from_slice(&sem);
                puttr(&mut tr);
                v += dv;
            }
            if verbose != 0 {
                warn!("semblance output for cdp=%d", cdpprev);
            }
            if !gottrace {
                break;
            }
            cdpprev = cdp;
            tr = incoming.unwrap();
            ix = 0;
            aoff.fill(0.0);
            for row in &mut datm { row.fill(0.0); }
            for row in &mut datn { row.fill(0.0); }
        }
        if gottrace {
            if ix >= nx as usize {
                err!("number of traces in cdp exceeds nx");
            }
            aoff[ix] = tr.offset as f32;
            datm[ix].copy_from_slice(&tr.data[..nt as usize]);
            ix += 1;
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    exit(CWP_Exit());
}
