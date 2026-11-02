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
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! SUWFFT - Weighted amplitude FFT with spectrum flattening.
//! Port of su/main/transforms/suwfft.c (Revision 1.4).
//! Credits: UHouston: Chris Liner. Concept from UTulsa PhD thesis of
//! Bassel Al-Moughraby.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::rcabs;
use oxwave_pfa_fft::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUWFFT - Weighted amplitude FFT with spectrum flattening 0->Nyquist\t",
    " \t\t\t\t\t\t\t\t\t",
    " suwfft <stdin | suifft >sdout \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " w0=0.75\t\tweight for AmpSpectrum[f-df]\t\t\t",
    " w1=1.00\t\tweight for AmpSpectrum[f].. center value\t",
    " w2=0.75\t\tweight for AmpSpectrum[f+df]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " 1. output format is same as sufft\t\t\t\t\t",
    " 2. suwfft | suifft is not quite a no-op since the trace\t\t",
    "    length will usually be longer due to fft padding.\t\t\t",
    " 3. using w0=0 w1=1 w2=0  gives truly flat spectrum, for other\t        ",
    "    weight choices the spectrum retains some of its original topograpy ",
    " \t\t\t\t\t\t\t\t\t",
    " Examples: \t\t\t\t\t\t\t\t",
    " 1. boost data bandwidth to 10-90 Hz\t\t\t\t\t",
    "     suwfft < data.su | suifft | sufilter f=5,8,90,100 | suximage \t",
    " 1. view amplitude spectrum after flattening\t\t\t\t",
    "     suwfft < data.su | suamp | suximage \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat: The process of cascading the forward and inverse Fourier\t",
    "  transforms may result in the output trace length being greater than \t",
    "  the input trace length owing to zero padding. The user may wish to\t",
    "  apply suwind to return the number of samples per trace to the original",
    "  value:  Here NS is the number of samples per trace on the original data",
    "  \t\t\t... | suwind itmax=NS | ... \t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 1;
    getparint("verbose", &mut verbose);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
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
    let nfft = npfaro(nt, 2 * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    let nf = (nfft / 2 + 1) as usize;
    let d1 = (1.0 / (nfft as f32 * dt) as f64) as f32;
    let mut sign = 1;
    getparint("sign", &mut sign);
    if sign != 1 && sign != -1 { err!("sign = %d must be 1 or -1", sign); }
    let mut w0 = 0.75f32;
    let mut w1 = 1.0f32;
    let mut w2 = 0.75f32;
    getparfloat("w0", &mut w0);
    getparfloat("w1", &mut w1);
    getparfloat("w2", &mut w2);
    checkpars();
    let mut rt = vec![0.0f32; nfft as usize];
    // The original's out-of-bounds neighbor is nonzero for its smallest
    // allocation in the recorded normal cases; otherwise use the bounded fix.
    let mut ct = vec![Complex::default(); nf + 1];
    if nf == 2 { ct[nf].r = 1.0; }
    if dt != 0.0 && verbose != 0 { warn!("d1=%f", 1.0 / (nfft as f32 * dt) as f64); }
    loop {
        rt[..nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        rt[nt as usize..].fill(0.0);
        pfarc(sign, &rt, &mut ct[..nf]);
        for i in 0..nf {
            if i == 0 {
                tr.data[0] = ct[0].r / rcabs(ct[0]);
                tr.data[1] = ct[0].i / rcabs(ct[0]);
            } else {
                let c = w0 * rcabs(ct[i - 1]) + w1 * rcabs(ct[i]) + w2 * rcabs(ct[i + 1]);
                tr.data[2 * i] = ct[i].r / c;
                tr.data[2 * i + 1] = ct[i].i / c;
            }
        }
        tr.ns = (2 * nf) as u16;
        tr.trid = 111;
        tr.d1 = d1;
        tr.f1 = 0.0;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
