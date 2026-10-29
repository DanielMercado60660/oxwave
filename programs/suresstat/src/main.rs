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

//! SURESSTAT - Surface consistent source and receiver statics calculation.
//! Port of su/main/statics/suresstat.c (Revision 1.20).
//! Credits: CWP: Timo Tjan, 4 October 1994.
//! Rewritten by Thomas Pratt, USGS, Feb. 2000.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SURESSTAT - Surface consistent source and receiver statics calculation",
    " \t\t\t\t\t\t\t\t\t",
    "   suresstat <stdin [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters: \t\t\t\t\t\t\t",
    " ssol=\t\toutput file source statics\t\t\t\t",
    " rsol=\t\toutput file receiver statics\t\t\t\t",
    " ntraces=\tnumber of traces in input data set (must be correct!)\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " ntpick=50 \tmaximum static shift (samples)         \t\t\t",
    " niter=5 \tnumber of iterations\t\t\t\t\t",
    " nshot=240 \tlargest shot number (fldr=1 to nshot)\t\t\t",
    " nr=335 \tlargest receiver number (tracf=1 to nr)\t\t\t",
    " nc=574 \tmaximum number of cmp's (for array allocation)\t\t",
    " sfold=96 \tmaximum shot gather fold\t\t\t\t",
    " rfold=96 \tmaximum receiver gather fold\t\t\t\t",
    " cfold=48 \tmaximum cmp gather fold\t\t\t\t\t",
    " sub=0 \tsubtract super trace 1 from super trace 2 (=1)\t\t",
    " \t\tsub=0 strongly biases static to a value of 0\t\t",
    " mode=0 \tuse global maximum in cross-correllation window\t\t",
    "\t\t=1 choose the peak perc=percent smaller than the global max.",
    " perc=10. \tpercent of global max (used only for mode=1)\t\t",
    " verbose=0 \tprint diagnostic output (verbose=1)                     ",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Estimates surface-consistent source and receiver statics, meaning that",
    " there is one static correction value estimated for each shot and receiver",
    " position.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The method employed here is based on the method of Ronen and Claerbout:",
    " Geophysics 50, 2759-2767 (1985).\t\t\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " The output files are binary files containing the source and receiver\t",
    " statics, as a function of shot number (trace header fldr) and      \t",
    " receiver station number (trace header tracf). \t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " The code builds a supertrace1 and supertrace2, which are subsequently\t",
    " cross-correllated. The program then picks the time lag associated with",
    " the largest peak in the cross-correllation according to two possible\t",
    " criteria set by the parameter \"mode\". If mode=0, the maximum of the\t",
    " cross-correllation window is chosen. If mode=1, the program will pick ",
    " a peak which is up to perc=percent smaller than the global maximum, but",
    " closer to zero lag than the global maximum.\t(Choosing mode=0 is\t",
    " recommended.)\t\t\t\t\t\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " The geometry can be irregular: the program simply computes a static \t",
    " correction for each shot record (fldr=1 to fldr=nshot), with any missing ",
    " shots being assigned a static of 0.  A static correction for each    \t",
    " receiver station (tracf=1 to tracf=nr) is calculated, with missing    ",
    " receivers again assigned a static of 0.                               ",
    "\t\t\t\t\t\t\t\t\t",
    " The ntracesces parameter must be equal to the number of prestack traces.",
    " The ntpick parameter sets the maximum allowable shift desired (in\t",
    "   samples NOT time).\t\t\t\t\t\t\t",
    " The niter parameter sets the number of iterations desired.\t\t",
    " The nshot parameter must be equal to the maximum fldr number in\t",
    "     the data. Note that this number might be different from the actual",
    "     number of shot records in the data (i.e., the maximum ep number).\t",
    "     For getting the correct maximum fldr number, you may use the surange",
    "     command.\t\t\t\t\t\t\t\t",
    " The nr parameter must be equal to the largest number of receivers\t",
    "     per shot in the whole data.\t\t\t\t\t",
    " The nc parameter must be equal to the number of prestack traces in\t",
    "     the data.\t\t\t\t\t\t\t\t",
    " The sfold parameter must be equal to the nr parameter.\t\t",
    " The rfold parameter must be equal to the maximum ep number.\t\t",
    " The cfold parameter must be equal to the maximum CDP fold,\t\t",
    "     which is equal to the maximum number under the cdpt entry in the\t",
    "     output of the surange command.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " To apply the static corrections, use sustatic with hdrs=3\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get info from first trace and store first header */
    if gettra(&mut tr, 0) == 0 {
        err!("can't get first trace");
    }
    let mut dt = (tr.dt as f64 / 1000.0) as f32;
    let mdt = (tr.dt as f64 / 1000.0) as f32;
    if dt == 0.0 {
        getparfloat("dt", &mut dt);
    }
    if dt == 0.0 {
        mustgetparfloat("dt", &mut dt);
    }
    let nt = tr.ns as i32;

    /* Get optional parameters */
    let mut ntcc = 250i32;
    getparint("ntcc", &mut ntcc);
    let mut icmpshift = 9i32;
    getparint("icmpshift", &mut icmpshift);
    if ntcc % 2 == 0 {
        ntcc = ntcc.wrapping_add(1);
    }
    let mut ntpick = 50i32;
    getparint("ntpick", &mut ntpick);
    let mut ntraces = 50i32;
    getparint("ntraces", &mut ntraces);
    let mut resamp = 4i32;
    getparint("resamp", &mut resamp);
    let mut n_o = 7i32;
    getparint("n_o", &mut n_o);
    let mut niter = 5i32;
    getparint("niter", &mut niter);
    let mut nshot = 240i32;
    getparint("nshot", &mut nshot);
    let mut nr = 335i32;
    getparint("nr", &mut nr);
    let mut nc = 574i32;
    getparint("nc", &mut nc);
    let mut sfold = 96i32;
    getparint("sfold", &mut sfold);
    let mut rfold = 96i32;
    getparint("rfold", &mut rfold);
    let mut cfold = 48i32;
    getparint("cfold", &mut cfold);
    let mut mode = 0i32;
    getparint("mode", &mut mode);
    let mut sub = 0i32;
    getparint("sub", &mut sub);
    let mut perc = 10.0f32;
    getparfloat("perc", &mut perc);
    let mut verbose = 0i32;
    getparint("verbose", &mut verbose);
    let mut ssol = String::new();
    let mut rsol = String::new();
    if getparstring("ssol", &mut ssol) == 0 {
        err!("must specify a source statics output file");
    }
    if getparstring("rsol", &mut rsol) == 0 {
        err!("must specify a receiver statics output file");
    }
    checkpars();

    /* Compute time windowing parameters */
    let nt_r = ntcc.wrapping_mul(resamp);
    if verbose == 1 {
        warn!("nt_r=%i", nt_r);
    }
    let ntout = 2i32.wrapping_mul(ntpick).wrapping_add(1);
    let nt_super = nt.wrapping_add(2i32.wrapping_mul(ntpick));
    let i = if nr > nshot { nr } else { nshot };
    let ntotal_super = nt_super.wrapping_mul(i);
    if verbose == 1 {
        warn!(
            "nt=%i, nt_super=%i, ntotal_super=%i",
            nt, nt_super, ntotal_super
        );
    }

    /* Allocate space and zero arrays. The C indexes samples 1..=nt even
     * though alloc2float has stride nt; likewise header has stride 3 but
     * uses columns 1..=3. Keep these strides in flat storage. One extra
     * element makes the final one-past-end access explicit and safe. Gather
     * locations also retain their C strides, including row aliasing when
     * a fold is understated. Uninitialized C storage starts as zero here. */
    let mut cmpntr = vec![0i32; (nc + 1) as usize];
    let mut recntr = vec![0i32; (nr + 1) as usize];
    let mut shotntr = vec![0i32; (nshot + 1) as usize];
    let mut g_trace = vec![0.0f32; ntotal_super as usize + 1];
    if verbose == 1 {
        warn!("nr=%i rfold=%i", nr, rfold);
    }
    let mut cmp_loc = vec![0i32; (cfold + 1) as usize * (nc + 1) as usize];
    let mut rec_loc = vec![0i32; (rfold + 1) as usize * (nr + 1) as usize];
    let mut shot_loc = vec![0i32; (sfold + 1) as usize * (nshot + 1) as usize];
    let mut header = vec![0i32; 3 * (ntraces + 1) as usize + 1];
    let mut model = vec![0.0f32; nt as usize * (nc + 1) as usize + 1];
    let mut model_trace = vec![0.0f32; ntotal_super as usize + 1];
    let mut corr_trace = vec![0.0f32; ntout as usize + 1];
    let mut data = vec![0.0f32; nt as usize * (ntraces + 1) as usize + 1];
    let mut tsstat = vec![0.0f32; (nshot + 1) as usize];
    let mut trstat = vec![0.0f32; (nr + 1) as usize];
    let mut sstat = vec![0i32; (nshot + 1) as usize];
    let mut rstat = vec![0i32; (nr + 1) as usize];

    /* Read rest of data */
    for itrace in 0..ntraces {
        if itrace > 0 {
            gettra(&mut tr, itrace);
        }
        let icmp = tr.cdp as usize;
        let irec = tr.tracf as usize;
        let ishot = tr.fldr as usize;
        cmpntr[icmp] += 1;
        recntr[irec] += 1;
        shotntr[ishot] += 1;
        let j = cmpntr[icmp] as usize;
        let k = recntr[irec] as usize;
        let l = shotntr[ishot] as usize;
        cmp_loc[icmp * (cfold + 1) as usize + j] = itrace + 1;
        rec_loc[irec * (rfold + 1) as usize + k] = itrace + 1;
        shot_loc[ishot * (sfold + 1) as usize + l] = itrace + 1;
        header[3 * (itrace + 1) as usize + 1] = icmp as i32;
        header[3 * (itrace + 1) as usize + 2] = irec as i32;
        header[3 * (itrace + 1) as usize + 3] = ishot as i32;
        for it in 1..=nt {
            /* At ns=SU_NFLTS the C reads beyond tr.data. Its value is
             * undefined; report the invalid read instead of panicking. */
            if it as usize >= tr.data.len() {
                err!(
                    "sample index %d exceeds trace buffer (%d samples)",
                    it,
                    tr.data.len()
                );
            }
            data[((itrace + 1) * nt + it) as usize] = tr.data[it as usize];
        }
    }

    /* Start iterations */
    for iter in 1..=niter {
        if verbose == 1 {
            fprintf!(stderr(), "iteration #= %i\n", iter);
        }
        /* Construct CMP stack */
        if verbose == 1 {
            fprintf!(stderr(), "constructing CMP stack\n");
        }
        for i in 1..=nc {
            for it in 1..=nt {
                model[(i * nt + it) as usize] = 0.0;
            }
            for j in 1..=cmpntr[i as usize] {
                for it in 1..=nt {
                    model[(i * nt + it) as usize] +=
                        data[(cmp_loc[(i * (cfold + 1) + j) as usize] * nt + it) as usize];
                }
            }
        }

        /* Loop over shots */
        for ishot in 1..=nshot {
            if verbose == 1 {
                fprintf!(stderr(), "\nstarting shot loop, shot=%i\n", ishot);
            }
            if shotntr[ishot as usize] > 1 {
                /* Construct shot and cmp super traces */
                if verbose == 1 {
                    fprintf!(stderr(), "creating supertraces\n");
                }
                for it in 1..=ntotal_super {
                    g_trace[it as usize] = 0.0;
                }
                for it in 1..=ntotal_super {
                    model_trace[it as usize] = 0.0;
                }
                for ichan in 1..=shotntr[ishot as usize] {
                    let j = (ichan - 1) * nt_super + ntpick;
                    let l =
                        header[3 * shot_loc[(ishot * (sfold + 1) + ichan) as usize] as usize + 1];
                    for it in 1..=nt {
                        g_trace[(j + it) as usize] = data
                            [(shot_loc[(ishot * (sfold + 1) + ichan) as usize] * nt + it) as usize];
                        model_trace[(j + it) as usize] = model[(l * nt + it) as usize];
                    }
                }
                /* Subtract shot supertrace from cmp supertrace if desired */
                if sub == 1 && verbose == 1 {
                    fprintf!(stderr(), "subtracting supertraces\n");
                }
                if sub == 1 {
                    for it in 1..=ntotal_super {
                        model_trace[it as usize] = model_trace[it as usize] - g_trace[it as usize];
                    }
                }
                /* Cross-correlate super trace 1 with super trace 2 */
                if verbose == 1 {
                    fprintf!(stderr(), "crosscorrelating supertraces\n");
                }
                for i in 1..=ntout {
                    corr_trace[i as usize] = 0.0;
                    let tlag = i - ntpick - 1;
                    for j in ntpick..=(ntotal_super - ntpick) {
                        corr_trace[i as usize] = corr_trace[i as usize]
                            + g_trace[(j + tlag) as usize] * model_trace[j as usize];
                    }
                }
                /* Pick cross-correlation peak */
                let tlag = max(&corr_trace, mode, perc, ntout);
                sstat[ishot as usize] = tlag;
                tsstat[ishot as usize] += tlag as f32;
                if verbose == 1 {
                    fprintf!(
                        stderr(),
                        "finished shot %i, sstat=%i, tsstat=%f (samples)\n",
                        ishot,
                        sstat[ishot as usize],
                        tsstat[ishot as usize]
                    );
                }
            }
        }
        /* Correct traces for shot statics (use g_trace as temp trace) */
        for i in 1..=ntraces {
            let ishot = header[3 * i as usize + 3] as usize;
            for it in 1..=nt {
                g_trace[it as usize] = 0.0;
            }
            let s = sstat[ishot];
            for it in (1 - if s > 0 { 0 } else { s })..=(nt - if s > 0 { s } else { 0 }) {
                g_trace[it as usize] = data[(i * nt + it + s) as usize];
            }
            for it in 1..=nt {
                data[(i * nt + it) as usize] = g_trace[it as usize];
            }
        }
        /* Construct CMP stack of corrected traces */
        if verbose == 1 {
            fprintf!(
                stderr(),
                "constructing new CMP stack using shot-corrected traces\n"
            );
        }
        for i in 1..=nc {
            for it in 1..=nt {
                model[(i * nt + it) as usize] = 0.0;
            }
            for j in 1..=cmpntr[i as usize] {
                for it in 1..=nt {
                    model[(i * nt + it) as usize] +=
                        data[(cmp_loc[(i * (cfold + 1) + j) as usize] * nt + it) as usize];
                }
            }
        }

        /* Loop over receivers */
        if verbose == 1 {
            fprintf!(stderr(), "\n\nstarting reciever loop\n");
        }
        for irec in 1..=nr {
            if recntr[irec as usize] > 1 {
                /* Construct receiver and cmp super traces */
                if verbose == 1 {
                    fprintf!(stderr(), "creating supertraces\n");
                }
                for it in 1..=ntotal_super {
                    g_trace[it as usize] = 0.0;
                }
                for it in 1..=ntotal_super {
                    model_trace[it as usize] = 0.0;
                }
                for ichan in 1..=recntr[irec as usize] {
                    let j = (ichan - 1) * nt_super + ntpick;
                    let l = header[3 * rec_loc[(irec * (rfold + 1) + ichan) as usize] as usize + 1];
                    for it in 1..=nt {
                        g_trace[(j + it) as usize] = data
                            [(rec_loc[(irec * (rfold + 1) + ichan) as usize] * nt + it) as usize];
                        model_trace[(j + it) as usize] = model[(l * nt + it) as usize];
                    }
                }
                /* Subtract receiver supertrace from cmp supertrace if desired */
                if sub == 1 && verbose == 1 {
                    fprintf!(stderr(), "subtracting supertraces\n");
                }
                if sub == 1 {
                    for it in 1..=ntotal_super {
                        model_trace[it as usize] = model_trace[it as usize] - g_trace[it as usize];
                    }
                }
                /* Cross-correlate super trace 1 with super trace 2 */
                if verbose == 1 {
                    fprintf!(stderr(), "crosscorrelating supertraces\n");
                }
                for i in 1..=ntout {
                    corr_trace[i as usize] = 0.0;
                    let tlag = i - ntpick - 1;
                    for j in ntpick..=(ntotal_super - ntpick) {
                        corr_trace[i as usize] = corr_trace[i as usize]
                            + g_trace[(j + tlag) as usize] * model_trace[j as usize];
                    }
                }
                let tlag = max(&corr_trace, mode, perc, ntout);
                rstat[irec as usize] = tlag;
                trstat[irec as usize] += tlag as f32;
                if verbose == 1 {
                    fprintf!(
                        stderr(),
                        "finished receiver %i, rstat=%i, trstat=%f (samples)\n",
                        irec,
                        rstat[irec as usize],
                        trstat[irec as usize]
                    );
                }
                if verbose == 1 {
                    fprintf!(stderr(), "finished one receiver\n");
                }
            }
        }
        /* Correct traces for receiver statics (use g_trace as temp trace) */
        for i in 1..=ntraces {
            let irec = header[3 * i as usize + 2] as usize;
            for it in 1..=nt {
                g_trace[it as usize] = 0.0;
            }
            let s = rstat[irec];
            for it in (1 - if s > 0 { 0 } else { s })..=(nt - if s > 0 { s } else { 0 }) {
                g_trace[it as usize] = data[(i * nt + it + s) as usize];
            }
            for it in 1..=nt {
                data[(i * nt + it) as usize] = g_trace[it as usize];
            }
        }
        if iter < niter {
            if verbose == 1 {
                fprintf!(
                    stderr(),
                    "constructing new CMP stack using shot-corrected traces\n"
                );
            }
            for i in 1..=nc {
                for it in 1..=nt {
                    model[(i * nt + it) as usize] = 0.0;
                }
                for j in 1..=cmpntr[i as usize] {
                    for it in 1..=nt {
                        model[(i * nt + it) as usize] +=
                            data[(cmp_loc[(i * (cfold + 1) + j) as usize] * nt + it) as usize];
                    }
                }
            }
        }
    }
    /* Output final statics */
    for it in 1..=nshot {
        tsstat[it as usize] = tsstat[it as usize] * mdt;
    }
    let fps = efopen(&ssol, "wb");
    efwrite(&tsstat, &fps);
    efclose(&fps);
    for it in 1..=nr {
        trstat[it as usize] = trstat[it as usize] * mdt;
    }
    let fpr = efopen(&rsol, "wb");
    efwrite(&trstat, &fpr);
    efclose(&fpr);
    exit(CWP_Exit());
}

#[allow(dead_code)]
fn window(trace: &[f32], nt: i32, mut nnt: i32, ntrace: &mut [f32]) {
    if nnt % 2 == 0 {
        nnt += 1;
    }
    let ft = (nt - 1) / 2 - (nnt - 1) / 2;
    for j in 0..nnt {
        ntrace[j as usize] = trace[(j + ft) as usize];
    }
}

fn max(trace: &[f32], mode: i32, mut perc: f32, mut nt: i32) -> i32 {
    let mut maxamp = 0.0f32;
    let mut locmax = 0i32;
    let mut lag = 0i32;
    if nt % 2 == 0 {
        nt += 1;
    }
    let zero = (nt + 1) / 2;
    let mut globmax = (nt + 1) / 2;
    perc = (1.0f64 - (perc / 100.0f32) as f64) as f32;
    for j in 1..=nt {
        if trace[j as usize] > maxamp {
            maxamp = trace[j as usize];
            globmax = j;
        }
    }
    if mode == 1 {
        for j in 1..=nt {
            if trace[j as usize] > perc * maxamp && (j - zero).abs() < (locmax - zero).abs() {
                locmax = j;
            }
        }
    }
    match mode {
        0 => lag = globmax - zero,
        1 => lag = locmax - zero,
        _ => {}
    }
    lag
}
