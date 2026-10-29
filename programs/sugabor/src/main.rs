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
//! SUGABOR - Gaussian multifilter time-frequency analysis.
//! Port of su/main/transforms/sugabor.c (Revision 1.20).
//! Credits: CWP: John Stockwell, Oct 1994; holder option Oct 2004.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUGABOR -  Outputs a time-frequency representation of seismic data via",
    "\t        the Gabor transform-like multifilter analysis technique ",
    "\t\tpresented by Dziewonski, Bloch and  Landisman, 1969.\t",
    "\t\t\t\t\t\t\t\t\t",
    "    sugabor <stdin >stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t \t\t",
    "\tif dt is not set in header, then dt is mandatory\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tdt=(from header)\ttime sampling interval (sec)\t\t",
    "\tfmin=0\t\t\tminimum frequency of filter array (hz)\t",
    "\tfmax=NYQUIST \t\tmaximum frequency of filter array (hz)\t",
    "\tbeta=3.0\t\tln[filter peak amp/filter endpoint amp]\t",
    "\tband=.05*NYQUIST\tfilter bandwidth (hz) \t\t\t",
    "\talpha=beta/band^2\tfilter width parameter\t\t\t",
    "\tverbose=0\t\t=1 supply additional info\t\t",
    "\tholder=0\t\t=1 output Holder regularity estimate\t",
    "\t\t\t\t=2 output linear regularity estimate\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes: This program produces a muiltifilter (as opposed to moving window)",
    " representation of the instantaneous amplitude of seismic data in the\t",
    " time-frequency domain. (With Gaussian filters, moving window and multi-",
    " filter analysis can be shown to be equivalent.)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " An input trace is passed through a collection of Gaussian filters\t",
    " to produce a collection of traces, each representing a discrete frequency",
    " range in the input data. For each of these narrow bandwidth traces, a ",
    " quadrature trace is computed via the Hilbert transform. Treating the narrow",
    " bandwidth trace and its quadrature trace as the real and imaginary parts",
    " of a \"complex\" trace permits the \"instantaneous\" amplitude of each",
    " narrow bandwidth trace to be compute. The output is thus a representation",
    " of instantaneous amplitude as a function of time and frequency.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Some experimentation with the \"band\" parameter may necessary to produce",
    " the desired time-frequency resolution. A good rule of thumb is to run ",
    " sugabor with the default value for band and view the image. If band is",
    " too big, then the t-f plot will consist of stripes parallel to the frequency",
    " axis. Conversely, if band is too small, then the stripes will be parallel",
    " to the time axis. \t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t\t",
    " The Gabor transform is not a wavelet transform, but rather are sharp\t",
    " frame basis. However, it is nearly a Morlet continuous wavelet transform",
    " so the concept of Holder regularity may have some meaning. If you are\t",
    " computing Holder regularity of, say, a migrated seismic section, then",
    " set band to 1/3 of the frequency band of your data.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "    suvibro | sugabor | suximage\t\t\t\t\t",
    "    suvibro | sugabor | suxmovie n1= n2= n3= \t\t\t\t",
    "     (because suxmovie scales it's amplitudes off of the first panel,  ",
    "      may have to experiment with the wclip and bclip parameters\t",
    "    suvibro | sugabor | supsimage | ... ( your local PostScript utility)",
    "\t\t\t\t\t\t\t\t\t",
];

use oxwave_core::prelude::*;
use oxwave_cwp::complex::crmul;
use oxwave_cwp::cwp::nint;
use oxwave_hilbert::hilbert;
use oxwave_linear_regression::linear_regression;
use oxwave_pfa_fft::{npfaro, pfarc, pfacr, Complex};

fn gaussian_filter(fcent: f32, dt: f32, nfft: usize, alpha: f32, band: f32, filter: &mut [f32]) {
    let nf = nfft / 2 + 1;
    let onfft = (1.0 / nfft as f64) as f32;
    let df = onfft / dt;
    let mut iflower = nint(((fcent - band / 2.0) / df) as f64);
    if iflower < 0 { iflower = 0; }
    let mut ifupper = nint(((fcent + band / 2.0) / df) as f64);
    if ifupper > nf as i32 { ifupper = nf as i32; }
    filter.fill(0.0);
    for i in 0..nf {
        let f = i as f32 * df - fcent;
        if iflower <= i as i32 && i as i32 <= ifupper {
            filter[i] = ((-4.0 * alpha as f64 * f as f64 * f as f64).exp() * onfft as f64) as f32;
        }
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let ntr = tr.ntr;
    if !(0..=17).contains(&tr.trid) && tr.trid != 109 { warn!("input is not seismic data, trid=%d", tr.trid); }
    let nt = tr.ns as usize;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if dt == 0.0 { err!("dt field is zero and not getparred"); }
    let nyq = (0.5 / dt as f64) as f32;
    let nfft = npfaro(nt as i32, 2 * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d -- too big", nfft); }
    let nfft = nfft as usize;
    let nf = nfft / 2 + 1;
    let mut verbose = 0i32; let mut holder = 0i32;
    let mut fmin = 0.0f32; let mut fmax = 0.0f32; let mut band = 0.0f32; let mut beta = 0.0f32; let mut alpha = 0.0f32;
    getparint("verbose", &mut verbose); getparint("holder", &mut holder);
    if getparfloat("fmin", &mut fmin) == 0 { fmin = (0.0 * nyq as f64) as f32; }
    if getparfloat("fmax", &mut fmax) == 0 { fmax = nyq; }
    if getparfloat("band", &mut band) == 0 { band = (0.05 * nyq as f64) as f32; }
    if getparfloat("beta", &mut beta) == 0 { beta = 3.0; }
    if getparfloat("alpha", &mut alpha) == 0 { alpha = beta / (band * band); }
    checkpars();
    let nfilt = nint(((fmax - fmin) / band) as f64);
    if verbose != 0 {
        warn!("            alpha = %f", alpha); warn!("             beta = %f", beta);
        warn!("             band = %f", band); warn!("number of filters = %d", nfilt);
    }
    let mut rt = vec![0.0f32; nfft]; let mut qt = vec![0.0f32; nfft];
    let mut ct = vec![Complex::default(); nf]; let mut fct = vec![Complex::default(); nf];
    let mut filter = vec![0.0f32; nf];
    // C allocates a contiguous nt by nfft array. Holder mode can access
    // beyond it, including an unpopulated column; provide space for those accesses.
    let extra = if holder != 0 { nfilt.max(0) as usize + 1 } else { 0 };
    let mut tmpdata = vec![0.0f32; nt * nfft + extra];
    loop {
        rt[..nt].copy_from_slice(&tr.data[..nt]); rt[nt..].fill(0.0);
        pfarc(1, &rt, &mut ct);
        let mut tracr = 0i32; let mut fcent = fmin;
        for i in 0..nfilt.max(0) as usize {
            gaussian_filter(fcent, dt, nfft, alpha, band, &mut filter);
            for j in 0..nf { fct[j] = crmul(ct[j], filter[j]); }
            pfacr(-1, &fct, &mut rt);
            hilbert(&rt[..nt], &mut qt[..nt]);
            if holder == 0 {
                for j in 0..nt { tr.data[j] = ((rt[j] * rt[j] + qt[j] * qt[j]) as f64).sqrt() as f32; }
                tracr += 1; tr.tracr = tracr; tr.sx = tr.tracl; tr.gx = tr.tracr;
                tr.d1 = dt; tr.f2 = fmin; tr.d2 = band; tr.trid = 118; tr.ntr = ntr.wrapping_mul(nfilt);
                puttr(&mut tr);
            } else {
                for j in 0..nt { tmpdata[j * nt + i] = ((rt[j] * rt[j] + qt[j] * qt[j]) as f64).sqrt() as f32; }
            }
            fcent += band;
        }
        if holder != 0 {
            let count = nfilt.max(0) as usize;
            let mut x = vec![0.0f32; count.max(4)]; let mut y = vec![0.0f32; count.max(4)]; let mut coeff = [0.0f32; 4];
            for j in 0..nt {
                x[0] = 0.0;
                for i in 1..count {
                    let diff = tmpdata[j * nt + i + 1] - tmpdata[j * nt];
                    let a = if diff < 0.0 { -diff } else { diff };
                    if holder == 1 { y[i] = (a as f64).ln() as f32; x[i] = (i as f64).ln() as f32; }
                    else if holder == 2 { y[i] = a; x[i] = i as f32; }
                    else { err!("holder must = 0,1,2!"); }
                }
                y[0] = y[1] - (y[3] - y[2]);
                if count > 0 { linear_regression(&y[..count], &x[..count], &mut coeff); }
                tr.data[j] = (coeff[0] as f64 + 0.5) as f32;
            }
            tracr += 1; tr.tracr = tracr; tr.sx = tr.tracl; tr.gx = tr.tracr;
            tr.d1 = dt; tr.f2 = fmin; tr.d2 = band; tr.trid = 118; tr.ntr = ntr.wrapping_mul(nfilt);
            puttr(&mut tr);
        }
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
