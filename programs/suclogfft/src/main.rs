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

//! SUCLOGFFT - complex log frequency transform; port of suclogfft.c, Revision 1.12.
//! Credits: CWP: John Stockwell, Shuki Ronen, Chris Liner, Jack K. Cohen;
//! CENPET: Werner M. Heigl; U Montana: Bob Lankston.

use oxwave_core::prelude::*;
use oxwave_pfa_fft::{npfaro, pfarc, Complex};
use oxwave_phase_unwrap::{oppenheim_unwrap_phase, simple_unwrap_phase};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUCLOGFFT - fft real time traces to complex log frequency domain traces",
    " \t\t\t\t\t\t\t\t\t",
    " suclogftt <stdin >sdout sign=1 \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " sign=1\t\t\tsign in exponent of fft\t\t\t",
    " dt=from header\t\tsampling interval\t\t\t",
    " verbose=1\t\t=0 to stop advisory messages\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " .... phase unwrapping options .....\t\t\t\t   \t",
    " mode=suphase\tsimple jump detecting phase unwrapping\t\t\t",
    " \t\t=ouphase  Oppenheim's phase unwrapping\t\t\t",
    " unwrap=1       |dphase| > pi/unwrap constitutes a phase wrapping\t",
    " \t \t=0 no phase unwrapping\t(in mode=suphase  only)\t\t",
    " trend=1\tremove linear trend from the unwrapped phase\t   \t",
    " zeromean=0     assume phase(0)=0.0, else assume phase is zero mean\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " clogfft(F(t)) = log(FFT(F(t)) = log|F(omega)| + iphi(omega)\t\t",
    " where phi(omega) is the unwrapped phase. Note that 0< unwrap<= 1.0 \t",
    " allows phase unwrapping to be tuned, if desired. \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " To facilitate further processing, the sampling interval\t\t",
    " in frequency and first frequency (0) are set in the\t\t\t",
    " output header.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suclogfft unwrap=0 | suiclogfft is not quite a no-op since the trace\t",
    " length will usually be longer due to fft padding.\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveats: \t\t\t\t\t\t\t\t",
    " No check is made that the data ARE real time traces!\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Output is type complex. To view amplitude, phase or real, imaginary\t",
    " parts, use    suamp \t\t\t\t\t\t\t",
    " PI/unwrap = minimum dphase is assumed to constitute a wrap in phase\t",
    " for suphase unwrapping only \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Examples: \t\t\t\t\t\t\t\t",
    " suclogfft < stdin | suamp mode=real | .... \t\t\t\t",
    " suclogfft < stdin | suamp mode=imag | .... \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The real and imaginary parts of the complex log spectrum are the\t",
    " respective amplitude and phase (unwrapped) phase spectra of the input\t",
    " signal. \t\t\t\t\t\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Example:  Homomorphic wavelet estimation \t\t\t\t",
    " suclogfft < shotgather | suamp mode=real | sustack key=dt > real.su\t",
    " suclogfft < shotgather | suamp mode=imag | sustack key=dt > imag.su\t",
    " suop2 real.su imag.su op=zipper | suiclogfft | suminphase > wavelet.su",
    "  \t\t\t\t\t\t\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 1;
    let mut trend = 1;
    let mut zeromean = 0;
    getparint("verbose", &mut verbose);
    getparint("trend", &mut trend);
    getparint("zeromean", &mut zeromean);
    let mut mode = String::from("ouphase");
    getparstring("mode", &mut mode);
    let imode = if mode == "suphase" { 0 } else if mode == "ouphase" { 1 } else {
        err!("unknown mode=\"%s\", see self-doc", mode);
    };
    if verbose != 0 { warn!("mode = %s", mode); }
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let mut dt = if seismic {
        if verbose != 0 { warn!("input is seismic data, trid=%d", tr.trid); }
        (tr.dt as f64 / 1000000.0) as f32
    } else {
        if verbose != 0 { warn!("input is not seismic data, trid=%d", tr.trid); }
        tr.d1
    };
    if dt == 0.0 {
        dt = 0.004;
        if verbose != 0 { warn!("dt or d1 not set, assumed to be .004"); }
    }
    let nfft = npfaro(nt as i32, 2 * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    let nf = nfft as usize / 2 + 1;
    let d1 = (1.0 / (nfft as f32 * dt) as f64) as f32;
    let mut sign = 1;
    let mut unwrap = 1.0f32;
    getparint("sign", &mut sign);
    getparfloat("unwrap", &mut unwrap);
    if sign != 1 && sign != -1 { err!("sign = %d must be 1 or -1", sign); }
    checkpars();
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ct = vec![Complex::default(); nf];
    let mut phase = vec![0.0f32; nf];
    let mut real = vec![0.0f32; nf];
    let mut imag = vec![0.0f32; nf];
    if dt != 0.0 && verbose != 0 { warn!("d1=%f", 1.0 / (nfft as f32 * dt) as f64); }
    loop {
        phase.fill(0.0);
        real.fill(0.0);
        imag.fill(0.0);
        rt[..nt].copy_from_slice(&tr.data[..nt]);
        rt[nt..].fill(0.0);
        pfarc(sign, &rt, &mut ct);
        for i in 0..nf {
            let re = ct[i].r;
            let im = ct[i].i;
            let ampsqrd = re * re + im * im;
            if !((ampsqrd as f64).abs() <= 0.0) {
                tr.data[2*i] = (0.5 * (ampsqrd as f64).ln()) as f32;
                phase[i] = (im as f64).atan2(re as f64) as f32;
            } else {
                phase[i] = 0.0;
                tr.data[2*i] = 0.0;
            }
            real[i] = re;
            imag[i] = im;
        }
        if unwrap != 0.0 {
            if imode == 0 { simple_unwrap_phase(trend, zeromean, unwrap, &mut phase); }
            else { phase.fill(0.0); oppenheim_unwrap_phase(trend, zeromean, d1, &real, &imag, &mut phase); }
        }
        for i in 0..nf { tr.data[2*i+1] = phase[i]; }
        tr.ns = (2*nf) as u16;
        tr.trid = 111;
        tr.d1 = d1;
        tr.f1 = 0.0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
