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

//! SUCCEPSTRUM - compute the complex cepstrum of a seismic trace.
//! Port of su/main/transforms/succepstrum.c (Revision 1.5).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: John Stockwell, Dec 2010, based on sucepstrum.c.
//! Credits: Balazs Nemeth of Potash Corporation of Saskatchewan Inc.,
//! given to CWP in 2008.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{rcabs, Complex};
use oxwave_cwp::cwp::closeto;
use oxwave_pfa_fft::{npfao, pfacc, pfarc};
use oxwave_phase_unwrap::{oppenheim_unwrap_phase, simple_unwrap_phase};
use oxwave_sblas::sscal;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUCCEPSTRUM - Compute the complex CEPSTRUM of a seismic trace \t\t\t\t\t\t\t\t\t\t",
    "  sucepstrum < stdin > stdout\t\t\t\t\t   \t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t  \t",
    "\tnone\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t  \t",
    " sign1=1\t\tsign of real to complex transform\t\t",
    " sign2=-1\t\tsign of complex to complex (inverse) transform\t",
    "\t\t\t\t\t\t\t\t\t",
    " ...phase unwrapping .....\t\t\t\t\t\t",
    " mode=ouphase\t\tOppenheim's algorithm for phase unwrapping\t",
    "\t\t\t=suphase  simple unwrap phase\t\t\t",
    " unwrap=1\t |dphase| > pi/unwrap constitutes a phase wrapping\t",
    "\t\t\t(operative only for mode=suphase)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " trend=1\t\tderamp the phase, =0 do not deramp the phase\t",
    " zeromean=0\t\tassume phase starts at 0,  =1 phase is zero mean",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The cepstrum is defined as the fourier transform of the the decibel   ",
    " spectrum, as though it were a time domain signal.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " CC(t) = FT[ln[T(omega)] ] = FT[ ln|T(omega)| + i phi(omega) ]\t\t",
    "\tT(omega) = |T(omega)| exp(i phi(omega))\t\t\t\t",
    "       phi(omega) = unwrapped phase of T(omega)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Phase unwrapping:\t\t\t\t\t\t\t",
    " The mode=ouphase uses the phase unwrapping method of Oppenheim and\t",
    " Schaffer, 1975, which operates integrating the derivative of the phase",
    "\t\t\t\t\t\t\t\t\t",
    " The mode=suphase generates unwrapped phase assuming that jumps\t",
    " in phase larger than dphase=pi/unwrap constitute a phase wrapping. In this case",
    " the jump in phase is replaced with the average of the jumps in phase  ",
    " on either side of the location where the suspected phase wrapping occurs.",
    " \t\t\t\t\t\t\t\t\t",
    " In either mode, the user has the option of de-ramping the phase, by   ",
    " removing its linear trend via trend=1 and of deciding whether the \t",
    " phase starts at phase=0 or is of  zero mean via zeromean=1.\t\t",
];

const LOOKFAC: i32 = 4;
const PFA_MAX: i32 = 720720;
const SUPHASE: i32 = 5;
const OUPHASE: i32 = 6;

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut sign1 = 1;
    let mut sign2 = -1;
    let mut mode = String::from("ouphase");
    let mut unwrap = 1.0f32;
    let mut trend = 1;
    let mut zeromean = 0;
    /* Get parameters. */
    if getparint("sign1", &mut sign1) == 0 {
        sign1 = 1;
    }
    if getparint("sign2", &mut sign2) == 0 {
        sign2 = -1;
    }
    if getparstring("mode", &mut mode) == 0 {
        mode = "ouphase".into();
    }
    if getparfloat("unwrap", &mut unwrap) == 0 {
        unwrap = 1.0;
    }
    if getparint("trend", &mut trend) == 0 {
        trend = 1;
    }
    if getparint("zeromean", &mut zeromean) == 0 {
        zeromean = 0;
    }
    checkpars();
    let imode = if mode == "suphase" {
        SUPHASE
    } else if mode == "ouphase" {
        OUPHASE
    } else {
        err!("unknown mode=\"%s\", see self-doc", mode);
    };
    /* Get information from the first header. */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    /* Loop over traces. The C leaves ns unchanged. */
    loop {
        rcceps(
            sign1,
            sign2,
            imode,
            unwrap,
            trend,
            zeromean,
            nt,
            dt,
            &mut tr.data,
        );
        tr.trid = 111;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}

fn rcceps(
    sign1: i32,
    sign2: i32,
    imode: i32,
    unwrap: f32,
    trend: i32,
    zeromean: i32,
    nt: i32,
    dt: f32,
    data: &mut [f32],
) {
    /* Set up pfa fft. Preserve npfao, as used by the C program. */
    let nfftc = npfao(nt, LOOKFAC * nt);
    if nfftc >= SU_NFLTS as i32 || nfftc >= PFA_MAX {
        err!("Padded nt=%d--too big", nfftc);
    }
    let nf = (nfftc / 2 + 1) as usize;
    let snfftc = (1.0 / nfftc as f64) as f32;
    let df = snfftc / dt;
    // C's pfarc writes float slot nfftc+1, one past its allocation when
    // nfftc is odd. Reserve that slot, but retain nf for subsequent calls.
    let mut w = vec![Complex::default(); nf + (nfftc % 2) as usize];
    let mut xr = vec![0.0f32; nf];
    let mut xi = vec![0.0f32; nf];
    let mut a = vec![0.0f32; nfftc as usize];
    let mut p = vec![0.0f32; nfftc as usize];
    a[..nt as usize].copy_from_slice(&data[..nt as usize]);
    sscal(nt, snfftc, &mut a, 1);
    // pfarc retains C's packing and integer nfftc/2 even for odd lengths.
    pfarc(sign1, &a, &mut w);
    for iw in 0..nf {
        xr[iw] = w[iw].r;
        xi[iw] = w[iw].i;
        a[iw] = rcabs(w[iw]);
        p[iw] = (w[iw].i as f64).atan2(w[iw].r as f64) as f32;
    }
    /* Unwrap phase. */
    if imode == OUPHASE {
        p[0] = 0.0;
        oppenheim_unwrap_phase(trend, zeromean, df, &xr, &xi, &mut p[..nf]);
    } else {
        simple_unwrap_phase(trend, zeromean, unwrap, &mut p[..nf]);
    }
    /* Take the log of the amplitude. */
    for iw in 0..nf {
        if !closeto(a[iw], 0.0) {
            w[iw].r = (a[iw] as f64).ln() as f32;
            w[iw].i = p[iw];
        } else {
            w[iw].r = 0.0;
            w[iw].i = 0.0;
        }
    }
    // Keep the C transform length, including inadmissible nf=25 and 64.
    // The extra slot for odd pfarc output is not part of this transform.
    // C overruns w for repeated factors and incorporates neighboring buffers;
    // zero-padding or selecting a valid FFT length changes its output.
    // The frozen FFT kernels reject these lengths when the radix has no
    // modular inverse: (ifac, m) is (5, 5) for 25 and (16, 4) for 64.
    // pfamcc with a nonzero sample stride rejects them even with spare storage.
    // Matching C requires shared support for its out-of-bounds accesses and
    // neighboring allocation contents.
    pfacc(sign2, &mut w[..nf]);
    for iw in 0..nf {
        data[2 * iw] = w[iw].r;
        data[2 * iw + 1] = w[iw].i;
    }
}
