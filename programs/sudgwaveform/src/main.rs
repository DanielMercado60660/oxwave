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
//
// Portions of this file are a Rust translation of code that is
// Copyright (c) 2007 by the Society of Exploration Geophysicists.
// For more information, go to http://software.seg.org/2007/0004 .
// You must read and accept usage terms at:
// http://software.seg.org/disclaimer.txt before use.
//
// Original code: Werner M. Heigl, Apache Corporation, E&P Technology, 2006-2007,
// accompanying: Heigl, W. M., 2007, Computing Gaussian derivative waveforms of
// any order: Geophysics, 72(4), H39-H42, doi:10.1190/1.2716624.
// The SEG-derived code is marked below between `SEG-BEGIN` and `SEG-END` comments.
//
// Modification notice (SEG terms of use, item 4):
//   Nature:      translated from C to Rust; same algorithm, statement order and
//                floating-point precision as the C; no functional change intended.
//   Extent:      the SEG-BEGIN/SEG-END blocks in this file.
//   Purpose:     OxWave, a Rust port of CWP/SU (Seismic Un*x).
//   When:        October 2026.
//   Who:         the OxWave project, maintained by Daniel Mercado;
//                translation produced with AI coding agents and checked against the
//                original C by differential tests.
//   Affiliation: OxWave project (independent open-source project).
// Keep this notice in all redistributions.

//! SUDGWAVEFORM - Gaussian derivative waveform.
//! Port of sudgwaveform.c, revision 1.4.
//! Credits: Werner M. Heigl, February 2007.
//! Copyright (c) 2007 by the Society of Exploration Geophysicists.
//! http://software.seg.org/2007/0004 ; usage terms: http://software.seg.org/disclaimer.txt .

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_waveforms::deriv_n_gauss;

/// Self-documentation from the C source.
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t   \t\t\t\t\t\t\t\t",
    " SUDGWAVEFORM - make Gaussian derivative waveform in SU format\t\t",
    " \t   \t\t\t\t\t\t\t\t",
    "  sudgwaveform >stdout  [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " n=2    \torder of derivative (n>=1)\t\t\t\t",
    " fpeak=35\tpeak frequency\t\t\t\t\t\t",
    " nfpeak=n*n\tmax. frequency = nfpeak * fpeak\t\t\t\t",
    " nt=128\tlength of waveform\t\t\t\t\t",
    " shift=0\tadditional time shift in s (used for plotting)\t\t",
    " sign=1\tuse =-1 to change sign\t\t\t\t\t",
    " verbose=0\t=0 don't display diagnostic messages\t\t\t",
    "               =1 display diagnostic messages\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This code computes a waveform that is the n-th order derivative of a\t",
    " Gaussian. The variance of the Gaussian is specified through its peak\t",
    " frequency, i.e. the frequency at which the amplitude spectrum of the\t",
    " Gaussian has a maximum. nfpeak is used to compute maximum frequency,\t",
    " which in turn is used to compute the sampling interval. Increasing\t",
    " nfpeak gives smoother plots. In order to have a (pseudo-) causal\t",
    " pulse, the program computes a time shift equal to sqrt(n)/fpeak. An\t",
    " additional shift can be applied with the parameter shift. A positive\t",
    " value shifts the waveform to the right.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    " 2-loop Ricker: dgwaveform n=1\t>ricker2.su\t\t\t\t",
    " 3-loop Ricker: dgwaveform n=2 >ricker3.su\t\t\t\t",
    " Sonic transducer pulse: dgwaveform n=10 fpeak=300 >sonic.su\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " To display use suxgraph. For example:\t\t\t\t\t",
    " dgwaveform n=10 fpeak=300 | suxgraph style=normal &\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For other seismic waveforms, please use \"suwaveform\".\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// SEG-BEGIN sudgwaveform.c:main
/// Generate a Gaussian derivative waveform.
/// Copyright (c) 2007 by the Society of Exploration Geophysicists.
/// For more information, go to http://software.seg.org/2007/0004 .
/// You must read and accept usage terms at:
/// http://software.seg.org/disclaimer.txt before use.
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut n: i32 = 2;
    let mut fpeak = 35.0f32;
    let mut sign = 1;
    let mut nfpeak = 0;
    let mut shift = 0.0f64;
    let mut verbose = 0;
    let mut nt = 0;
    getparint("n", &mut n);
    getparfloat("fpeak", &mut fpeak);
    getparint("sign", &mut sign);
    if getparint("nfpeak", &mut nfpeak) == 0 { nfpeak = n.wrapping_mul(n); }
    getpardouble("shift", &mut shift);
    getparint("verbose", &mut verbose);
    let fmax = nfpeak as f32 * fpeak;
    let dt = 0.5 / fmax as f64;
    let t0 = shift + (n as f64).sqrt() / fpeak as f64;
    if getparint("nt", &mut nt) == 0 { nt = (2.0 * t0 / dt + 1.0) as i32; }
    checkpars();
    warn!("n=%d fpeak=%.0f fmax=%.0f t0=%f nt=%d dt=%.12f", n, fpeak, fmax, t0, nt, dt);
    if dt < 1e-6 { err!("single-precision exceeded: reduce nfpeak or fpeak"); }
    // A negative length becomes an impossible size_t allocation in C.
    if nt < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut w = vec![0.0f64; nt as usize];
    if verbose != 0 { warn!("memory for waveform allocated and initialized"); }
    if n >= 1 {
        // The C library reads w[0] even for nt=0; no output samples are produced.
        if nt != 0 { deriv_n_gauss(dt, t0, fpeak, n, &mut w, sign, verbose); }
        for i in 0..nt as usize { tr.data[i] = w[i] as f32; }
        tr.tracl = 1;
        tr.ns = nt as u16;
        tr.trid = 1;
        tr.dt = nint(dt * 1000000.0) as u16;
        tr.ntr = 1;
        puttr(&mut tr);
        warn!("waveform written to stdout");
    } else { err!("specified n not >=1 !!"); }
    if verbose != 0 { warn!("memory freed"); }
    exit(CWP_Exit());
}
// SEG-END
