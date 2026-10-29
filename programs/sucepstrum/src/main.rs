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

//! SUCEPSTRUM - transform to the cepstral domain.
//! Port of su/main/transforms/sucepstrum.c (Revision 1.10).
//! Credits (from the C source):
//! CWP: John Stockwell, June 2013 based on sufft by:
//! CWP: Shuki Ronen, Chris Liner, Jack K. Cohen
//! CENPET: Werner M. Heigl - added well log support
//! U Montana: Bob Lankston - added m_unwrap_phase feature
//!
//! Note: leave dt set for later inversion.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
use oxwave_f_smooth::*;
use oxwave_pfa_fft::*;
use oxwave_phase_unwrap::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUCEPSTRUM - transform to the CEPSTRal domain\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  sucepstrum <stdin >sdout sign1=1 \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " sign1=1\t\t\tsign in exponent of fft\t\t\t",
    " sign2=-1\t\t\tsign in exponent of ifft\t\t",
    " dt=from header\t\tsampling interval\t\t\t",
    " verbose=1\t\t\t=0 to stop advisory messages\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " .... phase unwrapping options .....\t\t\t\t   \t",
    " mode=ouphase\tOppenheim's phase unwrapping\t\t\t\t",
    "\t\t=suphase simple jump detecting phase unwrapping\t\t",
    " unwrap=1       |dphase| > pi/unwrap constitutes a phase wrapping\t",
    " \t \t=0 no phase unwrapping\t(in mode=suphase  only)\t\t",
    " trend=1\tremove linear trend from the unwrapped phase\t   \t",
    " zeromean=0     assume phase(0)=0.0, else assume phase is zero mean\t",
    " smooth=0      apply damped least squares smoothing to unwrapped phase ",
    " r=10     ... damping coefficient, only active when smooth=1           ",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The complex log fft of a function F(t) is given by:\t\t\t",
    " clogfft(F(t)) = log(FFT(F(t)) = log|F(omega)| + iphi(omega)\t\t",
    " where phi(omega) is the unwrapped phase. Note that 0< unwrap<= 1.0 \t",
    " allows phase unwrapping to be tuned, if desired. \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The ceptrum is the inverse Fourier transform of the log fft of F(t) \t",
    " F(t_c) =cepstrum(F(t)) = INVFFT[log(FFT(F(t))]\t\t\t",
    "                        =INVFFT[log|F(omega)| + iphi(omega)]\t\t",
    " Here t_c is the cepstral time domain. \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " To facilitate further processing, the sampling interval\t\t",
    " in quefrency and first quefrency (0) are set in the\t\t\t",
    " output header.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveats: \t\t\t\t\t\t\t\t",
    " No check is made that the data ARE real time traces!\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Use suminphase to make minimum phase representations of signals \t",
];

const LOOKFAC: i32 = 2; // Look ahead factor for npfaro.
const PFA_MAX: i32 = 720720; // Largest allowed nfft.
const SIMPLE: i32 = 0;
const OPPENHEIM: i32 = 1;

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 1;
    getparint("verbose", &mut verbose);
    let mut trend = 1;
    getparint("trend", &mut trend);
    let mut zeromean = 0;
    getparint("zeromean", &mut zeromean);
    let mut mode = String::from("ouphase");
    getparstring("mode", &mut mode);
    let imode = if mode == "suphase" {
        SIMPLE
    } else if mode == "ouphase" {
        OPPENHEIM
    } else {
        err!("unknown mode=\"%s\", see self-doc", mode);
    };
    if verbose != 0 {
        warn!("mode = %s", mode);
    }
    let mut smooth = 0;
    getparint("smooth", &mut smooth);
    let mut r = 10.0f32;
    getparfloat("r", &mut r);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    // ISSEISMIC(tr.trid); the documented dt parameter is never read in C.
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let mut dt = if seismic {
        if verbose != 0 {
            warn!("input is seismic data, trid=%d", tr.trid);
        }
        (tr.dt as f64 / 1000000.0) as f32
    } else {
        if verbose != 0 {
            warn!("input is not seismic data, trid=%d", tr.trid);
        }
        tr.d1
    };
    if dt == 0.0 {
        dt = 0.004;
        if verbose != 0 {
            warn!("dt or d1 not set, assumed to be .004");
        }
    }
    // Set up pfa fft.
    let nfft = npfaro(nt as i32, LOOKFAC * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= PFA_MAX {
        err!("Padded nt=%d--too big", nfft);
    }
    let nfft = nfft as usize;
    let nf = nfft / 2 + 1;
    let d1 = (1.0 / (nfft as f32 * dt) as f64) as f32;
    let mut sign1 = 1;
    getparint("sign1", &mut sign1);
    let mut sign2 = -1;
    getparint("sign2", &mut sign2);
    let mut unwrap = 1.0f32;
    getparfloat("unwrap", &mut unwrap);
    if sign1 != 1 && sign1 != -1 {
        err!("sign1 = %d must be 1 or -1", sign1);
    }
    if sign2 != -1 && sign2 != 1 {
        err!("sign2 = %d must be -1 or 1", sign2);
    }
    checkpars();
    let mut rt = vec![0.0f32; nfft];
    let mut ct = vec![Complex::default(); nf];
    let mut phase = vec![0.0f32; nf];
    let mut real = vec![0.0f32; nf];
    let mut imag = vec![0.0f32; nf];
    let mut clogfft = vec![Complex::default(); nf];
    if dt != 0.0 && verbose != 0 {
        warn!("d1=%f", 1.0 / (nfft as f32 * dt) as f64);
    }
    loop {
        phase.fill(0.0);
        real.fill(0.0);
        imag.fill(0.0);
        clogfft.fill(Complex::default());
        rt[..nt].copy_from_slice(&tr.data[..nt]);
        rt[nt..].fill(0.0);
        pfarc(sign1, &rt, &mut ct);
        for i in 0..nf {
            let re = ct[i].r;
            let im = ct[i].i;
            let ampsqrd = re * re + im * im;
            if !(ampsqrd.abs() <= 0.0) {
                clogfft[i].r = (0.5 * (ampsqrd as f64).ln()) as f32;
                phase[i] = (im as f64).atan2(re as f64) as f32;
            } else {
                phase[i] = 0.0;
                tr.data[2 * i] = 0.0;
            }
            real[i] = ct[i].r;
            imag[i] = ct[i].i;
        }
        // Unwrap phase. This guard also bypasses smoothing in both modes.
        if unwrap != 0.0 {
            // Smoothing array of x values.
            let mut xout = vec![0.0f32; nf];
            if smooth != 0 {
                for i in 0..nf {
                    xout[i] = i as f32;
                }
            }
            if imode == SIMPLE {
                simple_unwrap_phase(trend, zeromean, unwrap, &mut phase);
            } else {
                phase.fill(0.0);
                oppenheim_unwrap_phase(trend, zeromean, d1, &real, &imag, &mut phase);
            }
            // Smooth the unwrapped phase.
            if smooth != 0 {
                smooth_1(&xout, &mut phase, r);
            }
        }
        for i in 0..nf {
            clogfft[i].i = phase[i];
        }
        let nifft = 2 * nf - 2;
        let onifft = 1.0f32 / nifft as f32;
        let newd1 = if tr.d1 != 0.0 {
            onifft / tr.d1
        } else if tr.dt != 0 {
            (tr.dt as f64 / 1000000.0) as f32
        } else {
            0.0
        };
        let mut cepstrum = vec![0.0f32; nifft];
        pfacr(sign2, &clogfft, &mut cepstrum);
        for i in 0..nifft {
            tr.data[i] = cepstrum[i] * onifft;
        }
        tr.ns = nifft as u16;
        tr.trid = 1;
        tr.d1 = newd1;
        tr.f1 = 0.0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
