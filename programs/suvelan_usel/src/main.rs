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
//! SUVELAN_USEL, revision 1.7, velocity analysis.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Valmore Celis (2002); based on suvelan_.c by Dave Hale (1989).
//! References: Neidell and Taner, Geophysics 36, 498-509 (1971);
//! Celis, Selective-correlation velocity analysis, CSM thesis (2002).
use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUVELAN_USEL - compute stacking velocity panel for cdp gathers\t     ",
    "\t\tusing the UnNormalized Selective CrossCorrelation Sum\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " suvelan_usel <stdin >stdout [optional parameters]\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:\t\t\t\t\t\t\t     ",
    " nx=tr.cdpt              number of traces in cdp\t\t\t     ",
    " dx=tr.d2                offset increment\t\t\t\t     ",
    " nv=50                   number of velocities\t\t\t\t     ",
    " dv=50.0                 velocity sampling interval\t\t\t     ",
    " fv=1500.0               first velocity\t\t\t\t     ",
    " tau=0.5                 threshold for significance values                  ",
    " smute=1.5               samples with NMO stretch exceeding smute are zeroed",
    " dtratio=5               ratio of output to input time sampling intervals   ",
    " nsmooth=dtratio*2+1     length of smoothing window                         ",
    " verbose=0               =1 for diagnostic print on stderr\t\t     ",
    " pwr=1.0                 semblance value to the power      \t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Notes:\t\t\t\t\t\t\t\t     ",
    " UnNormalized Selective CrossCorrelation sum: is based on the coherence     ",
    " measure known as crosscorrelation sum. The difference is that the selective",
    " approach sum only crosscorrelation pairs with relatively large differential",
    " moveout, thus increasing the resolving power in the velocity spectra       ",
    " compared to that achieved by conventional methods.  The selection is made  ",
    " using a parabolic approximation of the differential moveout and imposing a ",
    " threshold for those differential moveouts.\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " That threshold is the parameter tau in this program, which varies between  ",
    " 0 to 1.  A value of tau=0, means conventional crosscorrelation sum is      ",
    " applied implying that all crosscorrelations are included in the sum. In    ",
    " contrast, a value of tau=1 (not recomended) means that only the            ",
    " crosscorrelation formed by the trace pair involving the shortest and longest",
    " offset is included in the sum. Intermediate values will produce percentages",
    " of the crosscorrelations included in the sum that will be shown in the     ",
    " screen before computing the velocity spectra. Typical values for tau are   ",
    " between 0.2 and 0.6, producing approximated percentages of crosscorrelations",
    " summed between 60% and 20%. The higher the value of tau the lower the     ",
    " percentage and higher the increase in the resolving power of velocity\t      ",
    " spectra.\t\t\t\t\t\t\t\t      ",
    "\t\t\t\t\t\t\t\t\t      ",
    " Keeping the percentage of crosscorrelations included in the sum between 20%",
    " and 60% will increase resolution and avoid the precense of artifacts in the",
    " results.  In data contaminated by random noise or statics distortions is    ",
    " recomended to mantaing the percentage of crosscorrelations included in the  ",
    " sum above 25%.  After computing the velocity spectra one might want to     ",
    " adjust the level and number of contours before velocity picking.  \t      ",
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
    let mut nx = tr.cdpt;
    let mut dx = tr.d2;
    let mut nv = 50;
    let mut dv = 50.0f32;
    let mut fv = 1500.0f32;
    let mut tau = 0.5f32;
    let mut smute = 5.0f32;
    let mut dtratio = 5;
    let mut verbose = 0;
    let mut pwr = 1.0f32;
    getparint("nx", &mut nx);
    getparfloat("dx", &mut dx);
    getparint("nv", &mut nv);
    getparfloat("dv", &mut dv);
    getparfloat("fv", &mut fv);
    getparfloat("tau", &mut tau);
    getparfloat("smute", &mut smute);
    if smute <= 1.0 {
        err!("smute must be greater than 1.0");
    }
    getparint("dtratio", &mut dtratio);
    let mut nsmooth = dtratio * 2 + 1;
    getparint("nsmooth", &mut nsmooth);
    getparint("verbose", &mut verbose);
    getparfloat("pwr", &mut pwr);
    if pwr < 0.0 {
        err!("we are not looking for noise: pwr < 0");
    }
    if pwr == 0.0 {
        err!("we are creating an all-white semblance: pwr = 0");
    }
    checkpars();
    // arm64 division by zero yields zero.
    let ntout = 1 + (nt - 1).checked_div(dtratio).unwrap_or(0);
    CHECK_NT("ntout", ntout);
    let dtout = dt * dtratio as f32;
    if verbose != 0 {
        fprintf!(stderr(), "\tnumber of output time samples is %d\n", ntout);
        fprintf!(stderr(), "\toutput time sampling interval is %g\n", dtout);
        fprintf!(stderr(), "\toutput time of first sample is %g\n", ft);
    }
    if nx <= 0 {
        err!("nx must be positive");
    }
    let mut aoff = vec![0.0f32; nx as usize];
    let mut datm = vec![vec![0.0f32; nt as usize]; nx as usize];
    let mut datn = datm.clone();
    let mut sem = vec![0.0f32; ntout as usize];
    let mut gottrace = true;
    let mut cdpprev = tr.cdp;
    let mut cdp = tr.cdp;
    let mut ix = 0usize;
    let mut m = 0i32;
    loop {
        if gottrace && tr.cdp == cdpprev {
            if ix >= nx as usize {
                err!("gather has more than nx traces");
            }
            aoff[ix] = tr.offset as f32;
            datm[ix].copy_from_slice(&tr.data[..nt as usize]);
            ix += 1;
        }
        if gottrace {
            cdp = tr.cdp;
        }
        if cdp != cdpprev || !gottrace {
            // The next CDP's trace must survive output header and sample edits.
            let pending = if gottrace { Some(tr.clone()) } else { None };
            let mut ol = (aoff[0] as f64).abs().floor() as i32;
            for x in 1..nx as usize {
                let oaux = (aoff[x] as f64).abs().floor() as i32;
                m = if oaux < ol { oaux } else { ol };
            }
            ol = (m as f32 / dx) as i32;
            let delta = (nx + ol)
                .wrapping_mul(nx + ol)
                .wrapping_sub(ol.wrapping_mul(ol));
            let j0t = (((tau * delta as f32 + ol.wrapping_mul(ol) as f32) as f64).sqrt()
                - ol as f64)
                .floor() as i32
                + 1;
            let mut beta = 0i32;
            for x in j0t..=nx {
                let arg = (x + ol).wrapping_mul(x + ol) as f32 - tau * delta as f32;
                beta = beta.wrapping_add(((arg as f64).sqrt() - ol as f64).floor() as i32 - 1);
            }
            let tncc = (nx - 1).wrapping_mul(nx) / 2;
            let pcc = beta as f32 * 100.0 / tncc as f32;
            warn!("\t---- the percentage of crosscorrelations included in the sum is = %6.2f ----\n\n",pcc);
            let mut v = fv;
            for _iv in 0..nv {
                for x in 0..nx as usize {
                    let offovs = aoff[x] * aoff[x] / (v * v);
                    let tnmute = (offovs as f64 / ((smute * smute) as f64 - 1.0)).sqrt() as f32;
                    let itmute = if tnmute > ft {
                        ((tnmute - ft) / dt) as i32
                    } else {
                        0
                    };
                    let mut tn = ft + itmute as f32 * dt;
                    for it in itmute..nt {
                        let ti = (((tn * tn + offovs) as f64).sqrt() - ft as f64) / dt as f64;
                        let ti = ti as f32;
                        let iti = ti as i32;
                        if iti < nt - 1 {
                            if it < 0 || iti < 0 {
                                err!("invalid NMO sample index");
                            }
                            let frac = ti - iti as f32;
                            datn[x][it as usize] = ((1.0 - frac as f64)
                                * datm[x][iti as usize] as f64
                                + (frac * datm[x][(iti + 1) as usize]) as f64)
                                as f32;
                        }
                        tn += dt;
                    }
                }
                for itout in 0..ntout {
                    let it = itout * dtratio;
                    let ismin = (it - nsmooth / 2).max(0);
                    let ismax = (it + nsmooth / 2).min(nt - 1);
                    let mut ff = 0.0f32;
                    for is in ismin..ismax {
                        if j0t < 0 || j0t >= nx || is < 0 || is >= nt {
                            err!("invalid correlation selection index");
                        }
                        let arg = (j0t + ol).wrapping_mul(j0t + ol) as f32 - tau * delta as f32;
                        let ktj = ((arg as f64).sqrt() - ol as f64).floor() as i32 - 1;
                        if ktj > nx {
                            err!("invalid correlation selection index");
                        }
                        let mut dsum = 0.0f32;
                        for ik in 0..ktj {
                            dsum += datn[ik as usize][is as usize];
                        }
                        let mut fac = datn[j0t as usize][is as usize] * dsum;
                        for x in j0t + 1..nx {
                            let arg1 =
                                (x - 1 + ol).wrapping_mul(x - 1 + ol) as f32 - tau * delta as f32;
                            let arg2 = (x + ol).wrapping_mul(x + ol) as f32 - tau * delta as f32;
                            let ktj1 = ((arg1 as f64).sqrt() - ol as f64).floor() as i32;
                            let ktj2 = ((arg2 as f64).sqrt() - ol as f64).floor() as i32 - 1;
                            if ktj1 <= ktj2 && (ktj1 < 0 || ktj2 >= nx) {
                                err!("invalid correlation selection index");
                            }
                            let mut temp = 0.0f32;
                            for ik in ktj1..=ktj2 {
                                temp += datn[ik as usize][is as usize];
                            }
                            dsum += temp;
                            let nsum = datn[x as usize][is as usize];
                            fac += nsum * dsum;
                        }
                        ff += fac;
                    }
                    if ff < 0.0 {
                        ff = 0.0;
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
            tr = pending.unwrap();
            ix = 0;
            aoff.fill(0.0);
            for row in &mut datm { row.fill(0.0); }
            for row in &mut datn { row.fill(0.0); }
            // Store the trace that triggered the boundary on the next iteration.
            continue;
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    exit(CWP_Exit());
}
