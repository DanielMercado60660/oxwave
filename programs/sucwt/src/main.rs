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

//! SUCWT - continuous wavelet amplitude and regularity analysis.
//! Port of su/main/transforms/sucwt.c (Revision 1.12).
//! Credits (from the C source): CWP: John Stockwell, Nov 2004.
//! Inspired in part by "bhpcwt" in the BHP_SU package, code written by
//! BHP: Michael Glinsky, c. 2002, based loosely on a Matlab CWT function.
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_f_convolution::*;
use oxwave_hilbert::*;
use oxwave_linear_regression::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUCWT - generates Continous Wavelet Transform amplitude, regularity\t",
    "         analysis in the wavelet basis\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "     sucwt < stdin [Optional parameters ] > stdout\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " base=10\tBase value for wavelet transform scales\t\t\t",
    " first=-1\tFirst exponent value for wavelet transform scales\t",
    " expinc=0.01\tExponent increment for wavelet transform scales\t\t",
    " last=1.5\tLast exponent value for wavelet transform scales\t",
    "\t\t\t\t\t\t\t\t\t",
    " Wavelet Parameters:\t\t\t\t\t\t\t",
    " wtype=0\t\t2nd derivative of Gaussian (Mexican hat)\t",
    "\t\t\t=1 4th derivative of Gaussian (witch's hat)\t",
    "\t\t\t=2 6th derivative of Gaussian (wizard's hat)\t",
    " nwavelet=1024\t\tnumber of samples in the wavelet\t\t",
    " xmin=-20\t\tminimum x value wavelet is computed\t\t",
    " xcenter=0\t\tcenter x value  wavelet is computed \t\t",
    " xmax=20\t\tmaximum x value wavelet is computed\t\t",
    " sigma=1\t\tsharpness parameter ( sigma > 1 sharper)\t",
    "\t\t\t\t\t\t\t\t\t",
    " verbose=0\t\tsilent, =1 chatty\t\t\t\t",
    " holder=0\t\t=1 compute Holder regularity estimate\t\t",
    " divisor=1.0\t\ta floating point number >= 1.0 (see notes)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " This is the CWT version of the time frequency analysis notion that is ",
    " applied in sugabor.\t\t\t\t\t\t\t",
    " The parameter base is the base of the power that is applied to scale\t",
    " the wavelet. Some mathematical literature assume base 2. Base 10 works",
    " well here.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Default option yields an output similar to that of sugabor. With the  ",
    " parameter holder=1 an estimate of the instantaneous Holder regularity ",
    " (the Holder exponent) is output for each input data value. The result ",
    " is a Holder exponent trace for each corresponding input data trace.\t",
    " \t\t\t\t\t\t\t\t\t",
    " The strict definition of the Holder exponent is the maximum slope of  ",
    " the rise of the spectrum in the log(amplitude) versus log(scale) domain:",
    " \t\t\t\t\t\t\t\t\t",
    " divisor=1.0 means the exponent is computed simply by fitting a line   ",
    " through all of the values in the transform. A value of divisor>1.0    ",
    " indicates that the Holder exponent is determined as the max of slopes ",
    " found in (total scales)/divisor length segments.\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Some experimentation with the parameters nwavelet, first, last, and   ",
    " expinc may be necessary before a desirable output is obtained. The\t",
    " most effective way to proceed is to perform a number of tests with    ",
    " holder=0 to determine the range of first, last, and expinc that best  ",
    " represents the data in the wavelet domain. Then experimentation with  ",
    " holder=1 and values of divisor>=1.0 may proceed.\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut base = 10.0f32;
    let mut first = -1.0f32;
    let mut expinc = 0.01f32;
    let mut last = 1.5f32;
    let mut wtype = 0;
    let mut nwavelet = 1024;
    let mut xmin = -20.0f32;
    let mut xcenter = 0.0f32;
    let mut xmax = 20.0f32;
    let mut sigma = 1.0f32;
    let mut holder = 0;
    let mut divisor = 1.0f32;
    let mut verbose = 0;
    getparfloat("base", &mut base);
    getparfloat("first", &mut first);
    getparfloat("expinc", &mut expinc);
    getparfloat("last", &mut last);
    getparint("wtype", &mut wtype);
    getparint("nwavelet", &mut nwavelet);
    getparfloat("xmin", &mut xmin);
    if getparfloat("xcenter", &mut xcenter) == 0 {
        // Preserve the C assignment to xmin, rather than xcenter.
        xmin = 0.0;
    }
    getparfloat("xmax", &mut xmax);
    getparfloat("sigma", &mut sigma);
    getparint("holder", &mut holder);
    getparfloat("divisor", &mut divisor);
    getparint("verbose", &mut verbose);
    checkpars();
    if verbose != 0 {
        warn!(
            "base=%f, first=%f, expinc=%f, last=%f",
            base, first, expinc, last
        );
    }
    // C allocates before checking wtype or nwavelet. A negative count is
    // converted to size_t and fails allocation, rather than reaching the
    // MexicanHatFunction check. Rust capacity errors do not set errno, so
    // supply the corresponding ENOMEM text explicitly.
    let mut xvalues = Vec::<f32>::new();
    if xvalues.try_reserve_exact(nwavelet as usize).is_err() {
        err!(
            "%s: malloc failed (%s)",
            "ealloc.c",
            oxwave_core::err::strerror(12)
        );
    }
    let mut wavelet = Vec::<f32>::new();
    if wavelet.try_reserve_exact(nwavelet as usize).is_err() {
        err!(
            "%s: malloc failed (%s)",
            "ealloc.c",
            oxwave_core::err::strerror(12)
        );
    }
    xvalues.resize(nwavelet as usize, 0.0);
    wavelet.resize(nwavelet as usize, 0.0);
    if wtype == 0 {
        mexican_hat_function(nwavelet, xmin, xcenter, xmax, sigma, &mut wavelet);
    } else {
        err!("%d  type of wavelet not yet implemented", wtype);
    }
    let waveletinc = (xmax - xmin) / (nwavelet - 1) as f32;
    if verbose != 0 {
        warn!(
            "xmin=%f, xmax=%f, nwavelet=%d, waveletinc=%f",
            xmin, xmax, nwavelet, waveletinc
        );
    }
    let mut x = xmin;
    for v in &mut xvalues {
        *v = x;
        x += waveletinc;
    }
    xvalues[nwavelet as usize - 1] = xmax;
    let mut exponent = first;
    x = 0.0;
    let minscale = (base as f64).powf(first as f64) as f32;
    let maxscale = (base as f64).powf(last as f64) as f32;
    let mut scales = Vec::new();
    while x <= maxscale {
        x = (base as f64).powf(exponent as f64) as f32;
        exponent += expinc;
        scales.push(x);
        if scales.len() == 32767 {
            err!("Too many scales, change params and re-run\n");
        }
    }
    scales.pop();
    let nscales = scales.len();
    let mut nconv = vec![0i32; nscales];
    let mut index = vec![0i32; nwavelet as usize];
    let mut waveletsum = vec![0.0f32; nwavelet as usize];
    // ealloc2float(nwavelet, nscales) is one contiguous allocation in C.
    let mut filt = vec![0.0f32; nwavelet as usize * nscales];
    let mut f = vec![0.0f32; nwavelet as usize];
    for i in (0..nwavelet as usize).rev() {
        xvalues[i] -= xvalues[0];
    }
    let dx = xvalues[1];
    xmax = xvalues[nwavelet as usize - 1];
    if verbose != 0 {
        warn!("first xvalues=%f, last xvalues=%f", xvalues[0], xmax);
        warn!("dx=%f, xmax=%f", dx, xmax);
    }
    let mut fmin = 0.0f32;
    for i in 0..nwavelet as usize {
        fmin += wavelet[i];
        waveletsum[i] = fmin * dx;
    }
    for i in 0..nscales {
        nconv[i] = 1 + (scales[i] * xmax) as i32;
        for j in 0..nconv[i] as usize {
            // C adds one in float before converting the result to int.
            index[j] = (1.0f32 + j as f32 / (scales[i] * dx)) as i32;
        }
        for j in 0..nconv[i] as usize {
            f[j] = waveletsum[(index[j] - 1) as usize];
        }
        let mut k = nconv[i] - 1;
        for j in 0..nconv[i] as usize {
            filt[i * nwavelet as usize + j] = f[k as usize];
            k -= 1;
        }
    }
    if verbose != 0 {
        warn!("Convolution Lengths");
        for i in 0..nscales {
            warn!("%d ", nconv[i]);
        }
        warn!("%d scales will be used for transforms", scales.len() as i32);
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("Cannot get first trace\n");
    }
    let ns = tr.ns as usize;
    let mut outtr = Trace::new();
    let mut rt = vec![0.0f32; ns];
    let mut qt = vec![0.0f32; ns];
    let mut tmpdata = vec![vec![0.0f32; nscales]; ns];
    // For a one-sample filter C reads sucwt_buff[-1] for the first sample.
    // That read is zero in the reference build. Keep an explicit leading
    // zero and offset the convolution output by one to reproduce it safely.
    let mut sucwt_buff = vec![0.0f32; ns + nconv[nscales - 1] as usize + 2];

    loop {
        outtr.d2 = waveletinc;
        outtr.f2 = minscale;
        // The original immediately overwrites d2 and f2 with the input header.
        outtr.hdr = tr.hdr;

        for i in 0..nscales {
            sucwt_buff.fill(0.0);
            convolve_cwp(
                0,
                &tr.data[..ns],
                0,
                &filt[i * nwavelet as usize..i * nwavelet as usize + nconv[i] as usize],
                0,
                &mut sucwt_buff[1..ns + 1],
            );
            for j in 0..ns {
                rt[j] = sucwt_buff[(j as i32 + nconv[i] / 2) as usize];
            }
            for j in (1..ns).rev() {
                rt[j] = rt[j] - rt[j - 1];
            }
            for j in 0..ns {
                rt[j] = (-(scales[i] as f64).sqrt() * rt[j] as f64) as f32;
            }
            hilbert(&rt, &mut qt);
            if holder == 0 {
                for j in 0..ns {
                    outtr.data[j] = ((rt[j] * rt[j] + qt[j] * qt[j]) as f64).sqrt() as f32;
                }
                outtr.cdpt = i as i32 + 1;
                puttr(&mut outtr);
            } else {
                for j in 0..ns {
                    tmpdata[j][i] = ((rt[j] * rt[j] + qt[j] * qt[j]) as f64).sqrt() as f32;
                }
            }
        }

        if holder != 0 {
            let mut x = vec![0.0f32; nscales];
            let mut y = vec![0.0f32; nscales];
            let mut lrcoeff = [0.0f32; 4];
            for j in 0..ns {
                let mut icount = 0i32;
                x[0] = 0.0;
                for i in 1..nscales {
                    if i > 1 && (tmpdata[j][i - 1] - tmpdata[j][1]) as f64 > 0.0 {
                        let diff = tmpdata[j][i] - tmpdata[j][1];
                        let abs = if diff < 0.0 { -diff } else { diff };
                        y[icount as usize] = (abs as f64).ln() as f32;
                        x[icount as usize] = ((scales[i] - scales[1]) as f64).ln() as f32;
                        icount += 1;
                    }
                }
                icount -= 1;
                if icount > 10 && divisor == 1.0 {
                    linear_regression(&y[..icount as usize], &x[..icount as usize], &mut lrcoeff);
                    outtr.data[j] = lrcoeff[0];
                } else if icount > 10 && divisor > 1.0 {
                    let mut maxalpha = 0.0f32;
                    let interval = icount as f32 / divisor;
                    let mut k = interval as i32;
                    while k < icount {
                        linear_regression(&y[..k as usize], &x[..k as usize], &mut lrcoeff);
                        maxalpha = if lrcoeff[0] > maxalpha {
                            lrcoeff[0]
                        } else {
                            maxalpha
                        };
                        k = (k as f32 + interval) as i32;
                    }
                    outtr.data[j] = maxalpha;
                } else if icount < 10 && divisor >= 1.0 {
                    outtr.data[j] = 0.0;
                } else if divisor < 1.0 {
                    err!("divisor = %f < 1.0!", divisor);
                }
                // C leaves the previous output sample unchanged when icount == 10.
            }
            puttr(&mut outtr);
        }
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}

fn mexican_hat_function(
    nwavelet: i32,
    xmin: f32,
    xcenter: f32,
    xmax: f32,
    sigma: f32,
    wavelet: &mut [f32],
) {
    if nwavelet <= 1 {
        err!("nwavelet must be greater than 1!");
    }
    let dxwavelet = ((xmax - xmin) / (nwavelet - 1) as f32) as f64;
    let twopi = 2.0 * oxwave_cwp::cwp::PI;
    let mult = 1.0 / ((sigma * sigma * sigma) as f64 * twopi.sqrt());
    for i in 0..nwavelet as usize {
        let x = (i as f64 * dxwavelet - xcenter as f64) as f32;
        wavelet[i] = (mult
            * ((x * x / (sigma * sigma)) as f64 - 1.0)
            * (((-x * x) as f64) / (2.0 * sigma as f64 * sigma as f64)).exp())
            as f32;
    }
}
