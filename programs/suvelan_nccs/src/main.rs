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
//! SUVELAN_NCCS - normalized cross-correlation velocity analysis.
//! Port of su/main/velocity_analysis/suvelan_nccs.c, revision 1.6.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: Valmore Celis; based on suvelan.c by Dave Hale, CWP, c. 1989.
//! Reference: Neidell and Taner, Geophysics 36, 498-509 (1971).

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUVELAN_NCCS - compute stacking VELocity panel for cdp gathers\t     ",
    "\t\tusing Normalized CrossCorrelation Sum \t                     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " suvelan_uccs <stdin >stdout [optional parameters]\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:\t\t\t\t\t\t\t     ",
    " nx=tr.cdpt              number of traces in cdp\t\t\t     ",
    " nv=50                   number of velocities\t\t\t\t     ",
    " dv=50.0                 velocity sampling interval\t\t\t     ",
    " fv=1500.0               first velocity\t\t\t\t     ",
    " smute=1.5               samples with NMO stretch exceeding smute are zeroed",
    " dtratio=5               ratio of output to input time sampling intervals   ",
    " nsmooth=dtratio*2+1     length of smoothing window                         ",
    " verbose=0               =1 for diagnostic print on stderr\t\t     ",
    " pwr=1.0                 semblance value to the power      \t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Notes:\t\t\t\t\t\t\t\t     ",
    " Normalized CrossCorrelation sum: sum all possible crosscorrelation\t     ",
    " trace pairs in a CMP gather for each trial velocity and zero-offset        ",
    " two-way travel time inside a time window. This coherence measure is        ",
    " normalized by dividing each crosscorrelation trace pair by the geometric   ",
    " mean of the energy, inside the chosen time window, of each trace pair      ",
    " involved in each crosscorrelation. Then, to achieve a maximum amplitude    ",
    " of unity, the result is multiplied by  2/(M(M-1)), which is the inverse    ",
    " of the total number of crosscorrelation. The normalization allows to\t     ",
    " bring out weak reflection as long as these reflections have moveouts close ",
    " to a hyperbola.\t\t\t\t\t\t\t     ",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut cdp = tr.cdp;
    let mut nx = tr.cdpt;
    getparint("nx", &mut nx);
    let inc = nx.wrapping_mul(nx.wrapping_sub(1)) / 2;
    let mut nv = 50; getparint("nv", &mut nv);
    let mut dv = 100.0f32; getparfloat("dv", &mut dv);
    let mut fv = 1500.0f32; getparfloat("fv", &mut fv);
    let mut smute = 1.5f32; getparfloat("smute", &mut smute);
    if smute <= 1.0 { err!("smute must be greater than 1.0"); }
    let mut dtratio = 5; getparint("dtratio", &mut dtratio);
    let mut nsmooth = dtratio * 2 + 1; getparint("nsmooth", &mut nsmooth);
    let mut verbose = 0; getparint("verbose", &mut verbose);
    let mut pwr = 1.0f32; getparfloat("pwr", &mut pwr);
    if pwr < 0.0 { err!("we are not looking for noise: pwr < 0"); }
    if pwr == 0.0 { err!("we are creating an all-white semblance: pwr = 0"); }
    checkpars();
    // arm64 C division by zero yields zero.
    let ntout = 1 + (nt - 1).checked_div(dtratio).unwrap_or(0);
    CHECK_NT("ntout", ntout);
    let dtout = dt * dtratio as f32;
    if verbose != 0 {
        fprintf!(stderr(), "\tnumber of output time samples is %d\n", ntout);
        fprintf!(stderr(), "\toutput time sampling interval is %g\n", dtout);
        fprintf!(stderr(), "\toutput time of first sample is %g\n", ft);
    }
    let n = nx.max(0) as usize;
    let mut norm = vec![0.0f32; n];
    // An extra slot accommodates the original's write before the CDP check.
    let mut aoff = vec![0.0f32; n + 1];
    let mut datm = vec![vec![0.0f32; nt as usize]; n + 1];
    let mut datn = vec![vec![0.0f32; nt as usize]; n];
    let mut sem = vec![0.0f32; ntout.max(0) as usize];
    let mut gottrace = true;
    let mut cdpprev = tr.cdp;
    let mut ix = 0usize;
    let mut temp1 = 0.0f32;
    loop {
        if gottrace {
            if ix >= datm.len() { datm.push(vec![0.0; nt as usize]); aoff.push(0.0); }
            aoff[ix] = tr.offset as f32;
            cdp = tr.cdp;
            datm[ix].copy_from_slice(&tr.data[..nt as usize]);
            ix += 1;
        }
        if cdp != cdpprev || !gottrace {
            let mut v = fv;
            for _ in 0..nv {
                for j in 0..n {
                    let offovs = (aoff[j] * aoff[j]) / (v * v);
                    let tnmute = ((offovs as f64 / ((smute * smute) as f64 - 1.0)).sqrt()) as f32;
                    let itmute = if tnmute > ft { ((tnmute - ft) / dt) as i32 } else { 0 };
                    let mut tn = ft + itmute as f32 * dt;
                    for it in itmute..nt {
                        let ti = ((((tn * tn + offovs) as f64).sqrt() - ft as f64) / dt as f64) as f32;
                        let iti = ti as i32;
                        if iti < nt - 1 && iti >= 0 && it >= 0 {
                            let frac = ti - iti as f32;
                            datn[j][it as usize] = ((1.0 - frac as f64) * datm[j][iti as usize] as f64 + (frac * datm[j][iti as usize + 1]) as f64) as f32;
                        }
                        tn += dt;
                    }
                }
                for itout in 0..ntout {
                    let it = itout * dtratio;
                    let ismin = (it - nsmooth / 2).max(0);
                    let ismax = (it + nsmooth / 2).min(nt - 1);
                    let mut ff = 0.0f32;
                    for j in 0..n {
                        for is in ismin..ismax { temp1 += datn[j][is as usize] * datn[j][is as usize]; }
                        norm[j] = ((temp1 as f64).sqrt() + 1.0) as f32;
                        temp1 = 0.0;
                    }
                    for is in ismin..ismax {
                        let mut dsum = 0.0f32;
                        let mut fac = 0.0f32;
                        for j in 0..n {
                            let nsum = datn[j][is as usize] / norm[j];
                            dsum += nsum;
                            fac += nsum * dsum;
                        }
                        ff += fac;
                    }
                    sem[itout as usize] = ff / inc as f32;
                }
                tr.offset = 0;
                tr.cdp = cdpprev;
                tr.ns = ntout as u16;
                tr.dt = (dtout as f64 * 1000000.0) as u16;
                tr.data[..ntout as usize].copy_from_slice(&sem);
                puttr(&mut tr);
                v += dv;
            }
            if verbose != 0 { warn!("semblance output for cdp=%d", cdpprev); }
            if !gottrace { break; }
            cdpprev = cdp;
            ix = n;
        }
        if gettr(&mut tr) == 0 { gottrace = false; }
    }
    exit(CWP_Exit());
}
