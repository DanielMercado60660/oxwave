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

//! SUPEF - Wiener (least squares) predictive error filtering.
//! Port of su/main/decon_shaping/supef.c (Revision 1.47).
//! Credits: CWP: Shuki Ronen, Jack K. Cohen, Ken Larner.
//! CWP: John Stockwell, added mixing feature (April 1998).
//! CSM: Tanya Slota (September 2005), added cdp feature.
//! Technical reference: A. Ziolkowski, "Deconvolution", page 91.

use oxwave_core::prelude::*;
use oxwave_correlation::xcor;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::cwp::nint;
use oxwave_f_cubicspline::{cakima, cmonot, csplin};
use oxwave_linear_interp::intlin;
use oxwave_toeplitz::stoepf;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUPEF - Wiener (least squares) predictive error filtering\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " supef <stdin >stdout  [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " dt is mandatory if not set in header\t\t\t \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " cdp= \t\t\tCDPs for which minlag, maxlag, pnoise, mincorr, ",
    "\t\t\tmaxcorr are set\t(see Notes)\t\t\t",
    " minlag=dt\t\tfirst lag of prediction filter (sec)\t\t",
    " maxlag=last\t\tlag default is (tmax-tmin)/20\t\t\t",
    " pnoise=0.001\t\trelative additive noise level\t\t\t",
    " mincorr=tmin\t\tstart of autocorrelation window (sec)\t\t",
    " maxcorr=tmax\t\tend of autocorrelation window (sec)\t\t",
    " wienerout=0\t\t=1 to show Wiener filter on each trace\t\t",
    " mix=1,...\t \tarray of weights (floats) for moving\t\t",
    "\t\t\t\taverage of the autocorrelations\t\t",
    " outpar=/dev/null\toutput parameter file, contains the Wiener filter",
    " \t\t\tif wienerout=1 is set\t\t\t\t",
    " method=linear\t for linear interpolation of cdp values\t\t\t",
    "\t\t       =mono for monotonic cubic interpolation of cdps\t",
    "\t\t       =akima for Akima's cubic interpolation of cdps\t",
    "\t\t       =spline for cubic spline interpolation of cdps\t",
    " \t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns, dt\t\t\t\t\t",
    " Trace header fields modified: none\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " 1) To apply spiking decon (Wiener filtering with no gap):\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Run the following command\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "    suacor < data.su | suximage perc=95\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " You will see horizontal stripe running across the center of your plot.",
    " This is the autocorrelation wavelet for each trace. The idea of spiking",
    " decon is to apply a Wiener filter with no gap to the data to collapse\t",
    " the waveform into a spike. The idea is to pick the width of the\t",
    " autocorrelation waveform _from beginning to end_ (not trough to trough)",
    " and use this time for MAXLAG_SPIKING:\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  supef < data.su maxlag=MAXLAG_SPIKING  > dataspiked.su\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " 2) Prediction Error Filter (i.e. gapped Wiener filtering)\t\t",
    " The purpose of gapped decon is to suppress repetitions in the data\t",
    " such as those caused by water bottom multiples.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " To look for the period of the repetitions\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "    suacor ntout=1000 < dataspiked.su | suximage perc=95\t\t",
    " or \t\t\t\t\t\t\t\t\t",
    "    suacor ntout=1000 < dataspiked.su | sustack key=dt |suxwigb\t",
    " \t\t\t\t\t\t\t\t\t",
    " The value of ntout must be larger than the default 100. The idea is\t",
    " to look for repetitions in the autocorrelation. These repetitions will",
    " appear as a family of parallel stripes above and below the main\t",
    " autocorrelation waveform. Or, if you stack the data, these will be\t",
    " repetitive spikes.  This repetition time is the GAP. We set \t\t",
    " MINLAG_PEF to the value of this repetition time.\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " We set the minlag to MINLAG_PEF = GAP\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " We set the maxlag to MAXLAG_PEF = GAP + MAXLAG_SPIKING\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  supef < dataspiked.su minlag=MINLAG_PEF maxlag=MAXLAG_PEF > datapef.su",
    " \t\t\t\t\t\t\t\t\t",
    " Some experimentation may be required to get a satisfactory result.\t",
    " In particular you may find that you need to reduce the value of the   ",
    " minlag \t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " 3) It may be effective to sort your data into cdp gathers with susort,",
    " and perform sunmo correction to the water speed with sunmo, prior to \t",
    " attempts to suppress water bottom multiples. After applying supef, the",
    " user should apply inverse nmo to undo the nmo to water speed prior to\t",
    " further processing. Or, do the predictive decon on fully nmo-corrected",
    " gathers.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " If you flatten your data with sunmo, then make sure that you turn off ",
    " the stretch mute by using smute=20\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  | sunmo vnmo=v1,v2,... tnmo=t1,t2,... smute=20 | supef ...\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For a filter expressed as a function of cdp, specify the array\t",
    "     cdp=cdp1,cdp2,...\t\t\t\t\t\t\t",
    " and for each cdp specified, specify the minlag and maxlag arrays as\t",
    "      minlag=min1,min2,...     maxlag=max1,max2,...   \t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " It is required that the number of minlag and maxlag values be equal to",
    " the number of cdp's specified.  If the number of\t\t\t",
    " values in these arrays does not equal the number of cdp's, only the first",
    " value will be used.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t\t",
    " The wienerout=1 option writes out the wiener filter to outpar, and   ",
    " the prediction error filter to stdout, which is\t\t\t",
    "     1,0,0,...,-wiener[0],...,-wiener[imaxlag-1] \t\t\t",
    " where the sample value of -wiener[0], is  iminlag in the pe-filter.\t",
    " The pe-filter is output as a SU format datafile, one pe-filter for each",
    " trace input.\t\t\t\t\t\t\t\t",
    "\t...| supef ... wienerout | suxwigb\t\t\t\t",
    " shows the prediction error filters\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get info from first trace */
    if gettr(&mut intrace) == 0 {
        err!("can't get first trace");
    }
    let nt = intrace.ns as i32;
    let mut dt = (intrace.dt as f64 / 1000000.0) as f32;
    if dt == 0.0 {
        mustgetparfloat("dt", &mut dt);
    }

    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 {
        outpar = "/dev/null".to_owned();
    }
    let outparfp = efopen(&outpar, "w");

    /* Count parameters and allocate space */
    let mut ncdp = countparval("cdp") as usize;
    let mut icdp = 0usize;
    let mut jcdp = 0usize;
    let (mut cdp, mut minlag, mut maxlag, mut iminlag, mut imaxlag);
    if ncdp > 0 {
        let nminlag = countparval("minlag") as usize;
        if nminlag != 0 {
            if nminlag != ncdp {
                err!("number of cdp values must equal number of minlag values");
            }
        } else {
            err!("cdp set, but minlag not getparred!");
        }
        let nmaxlag = countparval("maxlag") as usize;
        if nmaxlag != 0 {
            if nmaxlag != ncdp {
                err!("number of cdp values must equal number of maxlag values");
            }
        } else {
            err!("cdp set, but maxlag not getparred!");
        }

        cdp = vec![0.0f32; ncdp];
        minlag = vec![0.0f32; ncdp];
        maxlag = vec![0.0f32; ncdp];
        getparfloat("minlag", &mut minlag);
        getparfloat("maxlag", &mut maxlag);
        getparfloat("cdp", &mut cdp);
        for icdp in 0..ncdp - 1 {
            if cdp[icdp] > cdp[icdp + 1] {
                err!("cdp values must increase monotonically!");
            }
        }
        iminlag = vec![0i32; ncdp];
        imaxlag = vec![0i32; ncdp];
        while icdp < ncdp {
            iminlag[icdp] = nint((minlag[icdp] / dt) as f64);
            imaxlag[icdp] = nint((maxlag[icdp] / dt) as f64);
            if iminlag[icdp] < 1 {
                err!("minlag[%d]=%g too small", icdp, minlag[icdp]);
            }
            if imaxlag[icdp] >= nt {
                err!("maxlag[%d]=%g too large", icdp, maxlag[icdp]);
            }
            // C compares the allocation pointers here, not their elements.
            // Its successive lag allocations are in increasing address order;
            // no numerical minlag >= maxlag validation takes place.
            icdp += 1;
        }
    } else {
        ncdp = 1;
        cdp = vec![intrace.cdp as f32];
        minlag = vec![0.0f32; ncdp];
        maxlag = vec![0.0f32; ncdp];
        iminlag = vec![0i32; ncdp];
        imaxlag = vec![0i32; ncdp];
        if getparfloat("minlag", &mut minlag) != 0 {
            iminlag[0] = nint((minlag[0] / dt) as f64);
        } else {
            iminlag[0] = 1;
        }
        if getparfloat("maxlag", &mut maxlag) != 0 {
            imaxlag[0] = nint((maxlag[0] / dt) as f64);
        } else {
            imaxlag[0] = nint(0.05 * nt as f64);
        }
    }

    let mut wienerout = 0i32;
    getparint("wienerout", &mut wienerout);
    let mut pnoise = 0.001f32;
    getparfloat("pnoise", &mut pnoise);
    let mut mincorr = 0.0f32;
    let imincorr = if getparfloat("mincorr", &mut mincorr) != 0 {
        nint((mincorr / dt) as f64)
    } else {
        0
    };
    if imincorr < 0 {
        err!("mincorr=%g too small", mincorr);
    }
    let mut maxcorr = 0.0f32;
    let imaxcorr = if getparfloat("maxcorr", &mut maxcorr) != 0 {
        nint((maxcorr / dt) as f64)
    } else {
        nt
    };
    if imaxcorr > nt {
        err!("maxcorr=%g too large", maxcorr);
    }
    if imincorr >= imaxcorr {
        err!("mincorr=%g, maxcorr=%g", mincorr, maxcorr);
    }

    let mut nmix = countparval("mix") as usize;
    let mut mix;
    if nmix != 0 {
        mix = vec![0.0f32; nmix];
        getparfloat("mix", &mut mix);
    } else {
        nmix = 1;
        mix = vec![1.0f32; nmix];
    }
    for imix in 0..nmix {
        mix[imix] /= nmix as f32;
    }
    let mut method = "linear".to_owned();
    getparstring("method", &mut method);
    checkpars();

    let mut nlag = imaxlag[0].wrapping_sub(iminlag[0]).wrapping_add(1);
    let mut ncorr = imaxcorr.wrapping_sub(imincorr).wrapping_add(1);
    let mut lcorr = imaxlag[0].wrapping_add(1);
    let lagcount = nlag;
    let mixcount = lcorr as i64 * nmix as i64;
    let mut wiener = floats(nlag);
    let mut spiker = floats(nlag);
    let mut autocorr = floats(lcorr);
    let mut temp = floats(lcorr);
    let mut mixacorr = floats64(mixcount);
    let mut crosslag = iminlag[0];
    let mut oldcdp = intrace.cdp;
    let mut itr = 0i32;

    loop {
        itr += 1;
        /* If necessary, compute new filter parameters. The argument order
         * of the original interpolation calls is intentionally retained. */
        if intrace.cdp != oldcdp && ncdp > 1 && jcdp < ncdp - 1 {
            let cdptotal = (cdp[ncdp - 1] - intrace.cdp as f32) as i32;
            jcdp += 1;
            let mut cdpinterp = floats(cdptotal);
            let mut ominlag = floats(cdptotal);
            let mut omaxlag = floats(cdptotal);
            for k in 0..cdptotal as usize {
                cdpinterp[k] = intrace.cdp.wrapping_add(k as i32) as f32;
            }
            for k in 0..cdptotal as usize {
                ominlag[k] = sample(&minlag, k as i32) + k as f32;
            }
            for k in 0..cdptotal as usize {
                omaxlag[k] = sample(&maxlag, k as i32) + k as f32;
            }
            let first = method.as_bytes().first().copied();
            if first == Some(b'l') || ncdp == 1 {
                // C uses cdptotal input points and ncdp output points,
                // even when these exceed the allocations. Fresh malloc
                // padding is zero in the reference; retain that padding.
                let mut xin = floats(cdptotal);
                let mut yin = floats(cdptotal);
                for k in 0..cdptotal as usize {
                    xin[k] = sample(&minlag, k as i32);
                    yin[k] = sample(&cdp, k as i32);
                }
                cdpinterp.resize((cdptotal as usize).max(ncdp), 0.0);
                ominlag.resize((cdptotal as usize).max(ncdp), 0.0);
                omaxlag.resize((cdptotal as usize).max(ncdp), 0.0);
                if cdptotal > 0 {
                    intlin(
                        &xin,
                        &yin,
                        cdp[0],
                        cdp[ncdp - 1],
                        &cdpinterp[..ncdp],
                        &mut ominlag[..ncdp],
                    );
                    for k in 0..cdptotal as usize {
                        xin[k] = sample(&maxlag, k as i32);
                    }
                    intlin(
                        &xin,
                        &yin,
                        cdp[0],
                        cdp[ncdp - 1],
                        &cdpinterp[..ncdp],
                        &mut omaxlag[..ncdp],
                    );
                }
                icdp = 0;
                while icdp < ncdp {
                    iminlag[icdp] = nint((ominlag[icdp] / dt) as f64);
                    imaxlag[icdp] = nint((omaxlag[icdp] / dt) as f64);
                    icdp += 1;
                }
            } else if matches!(first, Some(b'm' | b'a' | b's')) {
                let mut zind = vec![[0.0f32; 4]; ncdp];
                // C allocates cdptotal entries but reads/writes ncdp entries
                // in intcub. Retain zero padding as in the linear branch.
                cdpinterp.resize((cdptotal as usize).max(ncdp), 0.0);
                ominlag.resize((cdptotal as usize).max(ncdp), 0.0);
                omaxlag.resize((cdptotal as usize).max(ncdp), 0.0);
                let coefficients = match first {
                    Some(b'm') => cmonot,
                    Some(b'a') => cakima,
                    _ => csplin,
                };
                coefficients(&minlag, &cdp, &mut zind);
                intcub(0, &minlag, &zind, &cdpinterp[..ncdp], &mut ominlag[..ncdp]);
                coefficients(&maxlag, &cdp, &mut zind);
                intcub(0, &maxlag, &zind, &cdpinterp[..ncdp], &mut omaxlag[..ncdp]);
                icdp = 0;
                while icdp < ncdp {
                    iminlag[icdp] = nint((ominlag[icdp] / dt) as f64);
                    imaxlag[icdp] = nint((omaxlag[icdp] / dt) as f64);
                    icdp += 1;
                }
            } else {
                err!("%s is an unknown interpolation method!\n", method);
            }

            nlag = imaxlag[jcdp].wrapping_sub(iminlag[jcdp]).wrapping_add(1);
            ncorr = imaxcorr.wrapping_sub(imincorr).wrapping_add(1);
            lcorr = imaxlag[jcdp].wrapping_add(1);
            wiener = floats(nlag);
            spiker = floats(nlag);
            autocorr = floats(lcorr);
            temp = floats(lcorr);
            mixacorr = floats64(lcorr as i64 * nmix as i64);
            crosslag = iminlag[jcdp];
            // The C memset still uses the initial mixbytes after resizing.
            mixacorr.resize(mixacorr.len().max(mixcount as usize), 0.0);
        }

        // C keeps the initial lagbytes even after reallocating the filters.
        wiener.resize(wiener.len().max(lagcount as usize), 0.0);
        spiker.resize(spiker.len().max(lagcount as usize), 0.0);
        wiener[..lagcount as usize].fill(0.0);
        spiker[..lagcount as usize].fill(0.0);
        autocorr.fill(0.0);
        temp.fill(0.0);

        /* Form autocorrelation vector. The data pointer is NOT advanced by
         * imincorr, and the default inclusive window includes data[nt]. */
        let input_len = ncorr as usize;
        if input_len <= intrace.data.len() {
            xcor(
                imincorr,
                &intrace.data[..input_len],
                imincorr,
                &intrace.data[..input_len],
                0,
                &mut autocorr,
            );
        } else {
            let mut input = intrace.data.clone();
            input.resize(input_len, 0.0);
            xcor(imincorr, &input, imincorr, &input, 0, &mut autocorr);
        }

        if sample(&autocorr, 0) == 0.0 {
            puttr(&mut intrace);
            if wienerout != 0 {
                warn!("NO Wiener filter, trace: %d", itr);
            }
            if gettr(&mut intrace) == 0 {
                break;
            }
            continue;
        }
        autocorr[0] = (autocorr[0] as f64 * (1.0 + pnoise as f64)) as f32;
        mixacorr[..lcorr as usize].copy_from_slice(&autocorr);
        for ilag in 0..lcorr as usize {
            for imix in 0..nmix {
                temp[ilag] += mixacorr[imix * lcorr as usize + ilag] * mix[imix];
            }
            autocorr[ilag] = temp[ilag];
        }
        for imix in (1..nmix).rev() {
            for ilag in 0..lcorr as usize {
                mixacorr[imix * lcorr as usize + ilag] =
                    mixacorr[(imix - 1) * lcorr as usize + ilag];
            }
        }

        /* Get inverse filter by Wiener-Levinson. Zero length in C writes
         * malloc padding but contributes no terms to the convolution. */
        if nlag > 0 {
            let mut r = floats(nlag);
            let mut g = floats(nlag);
            for k in 0..nlag as usize {
                r[k] = sample(&autocorr, k as i32);
                g[k] = sample(&autocorr, crosslag.wrapping_add(k as i32));
            }
            stoepf(&r, &g, &mut wiener, &mut spiker);
        }

        /* Convolve pefilter with trace - don't do zero multiplies */
        for i in 0..nt {
            let n = if i < imaxlag[jcdp] { i } else { imaxlag[jcdp] };
            let mut sum = intrace.data[i as usize];
            for j in iminlag[jcdp]..=n {
                sum -=
                    wiener[(j - iminlag[jcdp]) as usize] * sample(&intrace.data, i.wrapping_sub(j));
            }
            outtrace.data[i as usize] = sum;
        }

        if wienerout != 0 && autocorr[0] != 0.0 {
            warn!("Wiener filter, trace: %d", itr);
            // icdp is deliberately the loop counter left by setup, as in C.
            let output_maxlag = imaxlag.get(icdp).copied().unwrap_or(0);
            let output_minlag = iminlag.get(icdp).copied().unwrap_or(0);
            for i in 0..output_maxlag {
                fprintf!(
                    &outparfp,
                    "%10g%c ",
                    filter_sample(&wiener, &spiker, nlag, i),
                    ' ' as i32
                );
            }
            fprintf!(&outparfp, "\n");
            outtrace.hdr = intrace.hdr;
            outtrace.data[..nt as usize].fill(0.0);
            outtrace.data[0] = 1.0;
            for i in 0..output_maxlag {
                let k = output_minlag.wrapping_add(i);
                if k >= 0 && (k as usize) < outtrace.data.len() {
                    outtrace.data[k as usize] = -filter_sample(&wiener, &spiker, nlag, i);
                }
            }
            puttr(&mut outtrace);
        } else {
            outtrace.hdr = intrace.hdr;
            puttr(&mut outtrace);
        }
        oldcdp = intrace.cdp;
        if gettr(&mut intrace) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}

// Original allocation failures go through ealloc1float. Negative signed
// sizes convert to an enormous size_t, and malloc fails in the reference.
fn floats(n: i32) -> Vec<f32> {
    floats64(n as i64)
}

fn floats64(n: i64) -> Vec<f32> {
    if n < 0 {
        err!("%s: malloc failed (Cannot allocate memory)", "ealloc.c");
    }
    let mut v = Vec::new();
    if v.try_reserve_exact(n as usize).is_err() {
        err!("%s: malloc failed (Cannot allocate memory)", "ealloc.c");
    }
    v.resize(n as usize, 0.0);
    v
}

// Only used for C's reads outside the logical array. Fresh allocation
// padding and the static samples beyond ns are zero in the reference.
fn sample(data: &[f32], i: i32) -> f32 {
    if i < 0 {
        0.0
    } else {
        data.get(i as usize).copied().unwrap_or(0.0)
    }
}

// The C output loop reads imaxlag coefficients, although only nlag were
// allocated. In the reference's small malloc size classes the next allocation
// (spiker) starts at the following 16-byte boundary. Reproduce the observed
// overread explicitly: padding is zero, then the auxiliary Toeplitz solution.
fn filter_sample(wiener: &[f32], spiker: &[f32], nlag: i32, i: i32) -> f32 {
    if i < nlag {
        return sample(wiener, i);
    }
    let stride = ((nlag.max(1) + 3) / 4) * 4;
    if i >= stride {
        sample(spiker, i - stride)
    } else {
        0.0
    }
}
