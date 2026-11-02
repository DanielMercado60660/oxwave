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
//! SUVELAN_NSEL - selective normalized cross-correlation velocity analysis.
//! Port of su/main/velocity_analysis/suvelan_nsel.c, revision 1.7.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: Valmore Celis (2002); based on suvelan.c by Dave Hale (1989).
//! References: Neidell and Taner, Geophysics 36 (1971); Celis, CSM thesis (2002).

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUVELAN_NSEL - compute stacking VELocity panel for cdp gathers\t     ",
    "\t\tusing the Normalized Selective CrossCorrelation sum\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " suvelan_usel <stdin >stdout [optional parameters]\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:\t\t\t\t\t\t\t     ",
    " nx=tr.cdpt              number of traces in cdp\t\t\t     ",
    " dx=tr.d2 \t          offset increment\t\t\t\t     ",
    " nv=50                   number of velocities\t\t\t\t     ",
    " dv=100.0                velocity sampling interval\t\t\t     ",
    " fv=1500.0               first velocity\t\t\t\t     ",
    " tau=0.5                 threshold for significance values                  ",
    " smute=1.5               samples with NMO stretch exceeding smute are zeroed",
    " dtratio=5               ratio of output to input time sampling intervals   ",
    " nsmooth=dtratio*2+1     length of smoothing window                         ",
    " verbose=0               =1 for diagnostic print on stderr\t\t     ",
    " pwr=1.0                 semblance value to the power      \t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Notes:\t\t\t\t\t\t\t\t     ",
    " Normalized Selective CrossCorrelation Sum: is based on the coherence       ",
    " measure known as crosscorrelation sum. The difference is that the selective",
    " approach sum only crosscorrelation pairs with relatively large differential",
    " moveout, thus increasing the resolving power in the velocity spectra       ",
    " compared to that achieved by conventional methods. The normalization is    ",
    " achieved in much the same way of normalizing the conventional              ",
    " crosscorrelation sum.\t\t\t\t\t\t             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Each crosscorrelation is divided by the geometric mean\t\t     ",
    " of the energy of the traces involved, and the multiplying by a constant to ",
    " achieve maximum amplitude of unity. The constant is just the inverse of the",
    " total number of crosscorrelations included in the sum.  The selection is   ",
    " made using a parabolic approximation of the differential moveout and       ",
    " imposing a threshold for those differential moveouts.\t\t   \t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " That threshold is the parameter tau in this program, which varies between 0",
    " to 1.\t A value of tau=0, means conventional crosscorrelation sum is applied",
    " implying that all crosscorrelations are included in the sum. In contrast,  ",
    " a value of tau=1 (not recomended) means that only the crosscorrelation     ",
    " formed by the trace pair involving the shortest and longest offset is      ",
    " included in the sum. Intermediate values will produce percentages of the   ",
    " crosscorrelations included in the sum that will be shown in the screen     ",
    " before computing the velocity spectra. Typical values for tau are between  ",
    " 0.2 and 0.6, producing approximated percentages of crosscorrelations summed",
    " between 60% and 20%. The higher the value of tau the lower the percentage",
    " and higher the increase in the resolving power of velocity spectra.        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Keeping the percentage of crosscorrelations included in the sum between 20%",
    " and 60% will increase resolution and avoid the precense of artifacts in   ",
    " the results.  In data contaminated by random noise or statics distortions   ",
    " is recomended to mantaing the percentage of crosscorrelations included in   ",
    " the sum above 25%. After computing the velocity spectra one might want to  ",
    " adjust the level  and number of contours before velocity picking.\t      ",
];

use oxwave_core::prelude::*;
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut nx = tr.cdpt;
    getparint("nx", &mut nx);
    let mut dx = tr.d2;
    getparfloat("dx", &mut dx);
    let mut nv = 50;
    getparint("nv", &mut nv);
    let mut dv = 50.0;
    getparfloat("dv", &mut dv);
    let mut fv = 1500.0;
    getparfloat("fv", &mut fv);
    let mut tau = 0.5;
    getparfloat("tau", &mut tau);
    let mut smute = 15.0;
    getparfloat("smute", &mut smute);
    if smute <= 1.0 {
        err!("smute must be greater than 1.0");
    }
    let mut dtratio = 5;
    getparint("dtratio", &mut dtratio);
    let mut nsmooth = dtratio + 1;
    getparint("nsmooth", &mut nsmooth);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut pwr = 1.0;
    getparfloat("pwr", &mut pwr);
    if pwr < 0.0 {
        err!("we are not looking for noise: pwr < 0");
    }
    if pwr == 0.0 {
        err!("we are creating an all-white semblance: pwr = 0");
    }
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
    if nx <= 0 {
        err!("nx must be positive");
    }
    let mut aoff = vec![0.0f32; nx as usize];
    let mut norm = vec![0.0f32; nx as usize];
    let mut datm = vec![vec![0.0f32; nt as usize]; nx as usize];
    let mut datn = datm.clone();
    let mut sem = vec![0.0f32; ntout as usize];
    let mut gottrace = true;
    let mut cdpprev = tr.cdp;
    let mut ix = 0usize;
    let mut m = 0i32;
    loop {
        let incoming_cdp = tr.cdp;
        let boundary = gottrace && tr.cdp != cdpprev;
        if gottrace {
            if ix >= nx as usize {
                err!("more than nx traces in gather");
            }
            aoff[ix] = tr.offset as f32;
            datm[ix].copy_from_slice(&tr.data[..nt as usize]);
            ix += 1;
        }
        if boundary || !gottrace {
            let mut ol = (aoff[0] as f64).abs().floor() as i32;
            for i in 1..nx as usize {
                let oaux = (aoff[i] as f64).abs().floor() as i32;
                m = if oaux < ol { oaux } else { ol };
            }
            ol = (m as f32 / dx) as i32;
            let diff = (nx + ol - 1) * (nx + ol - 1) - ol * ol;
            let j0t = (((tau * diff as f32 + (ol * ol) as f32) as f64).sqrt() - ol as f64).floor()
                as i32
                + 2;
            let mut beta = 0i32;
            for i in j0t..=nx {
                beta += (((((i + ol - 1) * (i + ol - 1)) as f32 - tau * diff as f32) as f64).sqrt()
                    - ol as f64)
                    .floor() as i32;
            }
            let tncc = (nx - 1) * nx / 2;
            let pcc = beta as f32 * 100.0 / tncc as f32;
            warn!("\t---- the percentage of crosscorrelations included in the sum is = %6.2f ----\n\n",pcc);
            let mut v = fv;
            for _ in 0..nv {
                for i in 0..nx as usize {
                    let offovs = aoff[i] * aoff[i] / (v * v);
                    let tnmute = (offovs as f64 / ((smute * smute) as f64 - 1.0)).sqrt() as f32;
                    let itmute = if tnmute > ft {
                        ((tnmute - ft) / dt) as i32
                    } else {
                        0
                    };
                    if itmute < 0 {
                        err!("NMO mute index is outside the trace");
                    }
                    let mut tn = ft + itmute as f32 * dt;
                    for it in itmute..nt {
                        let ti =
                            ((((tn * tn + offovs) as f64).sqrt() - ft as f64) / dt as f64) as f32;
                        let iti = ti as i32;
                        if iti < nt - 1 {
                            if iti < 0 {
                                err!("NMO interpolation index is outside the trace");
                            }
                            let frac = ti - iti as f32;
                            datn[i][it as usize] = ((1.0 - frac as f64)
                                * datm[i][iti as usize] as f64
                                + (frac * datm[i][(iti + 1) as usize]) as f64)
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
                    for i in 0..nx as usize {
                        let mut temp1 = 0.0f32;
                        for is in ismin..ismax {
                            temp1 += datn[i][is as usize] * datn[i][is as usize];
                        }
                        norm[i] = ((temp1 as f64).sqrt() + 1.0) as f32;
                    }
                    if ismin < ismax && (j0t < 0 || j0t >= nx) {
                        err!("correlation selection is outside the gather");
                    }
                    for is in ismin..ismax {
                        let ktj = (((((j0t + ol - 1) * (j0t + ol - 1)) as f32 - tau * diff as f32)
                            as f64)
                            .sqrt()
                            - ol as f64)
                            .floor() as i32;
                        if ktj < 0 || ktj > nx {
                            err!("correlation selection is outside the gather");
                        }
                        let mut dsum = 0.0f32;
                        for ik in 0..ktj {
                            dsum += datn[ik as usize][is as usize] / norm[ik as usize];
                        }
                        let mut fac = datn[j0t as usize][is as usize] / norm[j0t as usize] * dsum;
                        for i in j0t + 1..nx {
                            let ktj1 = ((((((i - 1) + ol - 1) * ((i - 1) + ol - 1)) as f32
                                - tau * diff as f32)
                                as f64)
                                .sqrt()
                                - ol as f64)
                                .floor() as i32
                                + 1;
                            let ktj2 = (((((i + ol - 1) * (i + ol - 1)) as f32 - tau * diff as f32)
                                as f64)
                                .sqrt()
                                - ol as f64)
                                .floor() as i32;
                            if ktj1 <= ktj2 && (ktj1 < 0 || ktj2 >= nx) {
                                err!("correlation selection is outside the gather");
                            }
                            let mut temp = 0.0f32;
                            for ik in ktj1..=ktj2 {
                                temp += datn[ik as usize][is as usize] / norm[ik as usize];
                            }
                            dsum += temp;
                            let nsum = datn[i as usize][is as usize] / norm[i as usize];
                            fac += nsum * dsum;
                        }
                        ff += fac;
                    }
                    if ff < 0.0 {
                        ff = 0.0;
                    }
                    sem[itout as usize] = ff / beta as f32;
                }
                tr.offset = 0;
                tr.cdp = cdpprev;
                tr.ns = ntout as u16;
                // Convert through u32 before narrowing, as in the defined C fix.
                tr.dt = (dtout as f64 * 1000000.0) as u32 as u16;
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
            cdpprev = incoming_cdp;
            // The C correlation loops leave ix at the end of the allocation.
            ix = if nv > 0 { nx as usize } else { nx as usize + 1 };
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    exit(CWP_Exit());
}
