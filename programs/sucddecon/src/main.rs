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

//! SUCDDECON - complex division deconvolution.
//! Port of su/main/decon_shaping/sucddecon.c (Revision 1.18).
//! Credits: CWP: Ivan Vasconcelos; changes by John Stockwell.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_pfa_fft::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUCDDECON - DECONvolution with user-supplied filter by straightforward",
    " \t      Complex Division in the frequency domain\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sucddecon <stdin >stdout [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " filter= \t\tascii filter values separated by commas\t\t",
    " \t\t...or...\t\t\t\t\t\t",
    " sufile=\t\tfile containing SU traces to use as filter\t",
    "                       (must have same number of traces as input data\t",
    " \t\t\t for panel=1)\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " panel=0\t\tuse only the first trace of sufile as filter\t",
    " \t\t\t=1 decon trace by trace an entire gather\t",
    " pnoise=0.001\t\twhite noise factor for stabilizing results\t",
    "\t \t\t\t\t(see below)\t\t \t",
    " sym=0\t\tnot centered, =1 center the output on each trace",
    " verbose=0\t\tsilent, =1 chatty\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " For given time-domain input data I(t) (stdin) and deconvolution\t",
    " filter F(t), the frequency-domain deconvolved trace can be written as:",
    "\t\t\t\t\t\t\t\t\t",
    "\t I(f)\t\tI(f) * complex_conjugate[F(f)]\t\t\t",
    " D(f) = ----- ===> D(f) = ------------------------ \t\t\t",
    "\t F(f)\t\t|F(f)|^2 + delta\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The real scalar delta is introduced to prevent the resulting deconvolved",
    " trace to be dominated by frequencies at which the filter power is close",
    " to zero. As described above, delta is set to some fraction (pnoise) of ",
    " the mean of the filter power spectra. Time sampling rate must be the \t",
    " same in the input data and filter traces. If panel=1 the two input files",
    " must have the same number of traces. Data and filter traces don't need to",
    " necessarily have the same number of samples, but the filter trace length",
    " length be always equal or shorter than the data traces. \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The time delay 'delrt' of the sufile subtracted from the delrt of each",
    " input trace. This is designed to make this program closer to an inverse",
    " of suconv, which adds the time delay of the sufile to each trace.\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveats: \t\t\t\t\t\t\t\t",
    " Owing to zero padding by the Fourier transform, the number of samples ",
    " on the output traces will likely be greater than on input traces.\t",
    " The user may wish to use\t\t\t\t\t\t",
    "      sucddecon .... |  suwind itmin=0 itmax=NSminus1\t\t\t",
    " where NSminus1 is the numerical value of the number of samples on the ",
    " input trace minus 1.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The user may need to apply frequency filtering to get acceptable output",
    "   sucddecon  ...| sufilter f=f1,f2,f3,f4 \t\t\t\t",
    " where f1,f2,f3,f4 are an acceptable frequency range.\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The user may need to mute artifacts that appear in previously muted   ",
    " parts of the input data, as well.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns, delrt\t\t\t\t",
    " Trace header fields modified: none, delrt\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    let mut sutrace = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let nfft = npfaro(nt as i32, 5 * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d -- too big", nfft); }
    let nfft = nfft as usize;
    let mut panel = 0;
    let mut verbose = 0;
    let mut sym = 0;
    let mut pnoise = 0.0f32;
    getparint("panel", &mut panel);
    getparint("verbose", &mut verbose);
    getparint("sym", &mut sym);
    if getparfloat("pnoise", &mut pnoise) == 0 { pnoise = 0.001; }
    let mut sufile = String::new();
    let is_filter = getparstring("sufile", &mut sufile) == 0;
    let mut filter = Vec::<f32>::new();
    let mut fp = None;
    let nfilter;
    let mut delrtf = 0i32;
    if is_filter {
        nfilter = countparval("filter") as usize;
        if nfilter == 0 { warn!("must specify filter= desired filter"); err!(" or sufile= "); }
        filter.resize(nfilter, 0.0);
        getparfloat("filter", &mut filter);
        panel = 0;
    } else {
        let file = efopen(&sufile, "r");
        fgettr(&file, &mut sutrace);
        nfilter = sutrace.ns as usize;
        delrtf = sutrace.delrt as i32;
        fp = Some(file);
    }
    let nfiltby2 = nfilter / 2;
    checkpars();
    let size = (nfft + nfiltby2).max(nfilter).max(nt + nfiltby2);
    let mut rf = vec![0f32; size];
    let mut rt = vec![0f32; size];
    let mut cf = vec![Complex::default(); size.max(2)];
    let mut ct = vec![Complex::default(); size.max(2)];
    if is_filter { rf[..nfilter].copy_from_slice(&filter[..nfilter]); }
    else { rf[..nfilter].copy_from_slice(&sutrace.data[..nfilter]); }
    if verbose != 0 { warn!(" nt=%d nfilter=%d nfft=%d ", nt as i32, nfilter as i32, nfft as i32); }
    if let Some(ref f) = fp { erewind(f); }
    let mut ntr = 0;
    let mut sum_powspec = 0.0f32;
    loop {
        ntr += 1;
        if panel != 0 {
            if fgettr(fp.as_ref().unwrap(), &mut sutrace) == 0 {
                delrtf = sutrace.delrt as i32;
                if verbose != 0 { warn!("Out of traces in sufile at trace= %d !", ntr); }
            }
            rt[..nfft].fill(0.0);
            rf[..nfft].fill(0.0);
            rt[nfiltby2..nfiltby2+nt].copy_from_slice(&tr.data[..nt]);
            rf[..nfilter].copy_from_slice(&sutrace.data[..nfilter]);
        } else {
            rt[..nfft].fill(0.0);
            rt[nfiltby2..nfiltby2+nt].copy_from_slice(&tr.data[..nt]);
        }
        ct.fill(Complex::default());
        cf.fill(Complex::default());
        pfarc(1, &rt[..nfft], &mut ct);
        pfarc(1, &rf[..nfft], &mut cf);
        for i in 0..nfft { sum_powspec += rcabs(cf[i]); }
        let delta = pnoise * sum_powspec / nfft as f32;
        let cdelta = cmplx(delta, 0.0);
        for i in 0..nfft {
            let numer = cmul(ct[i], conjg(cf[i]));
            let denom = cadd(cmul(cf[i], conjg(cf[i])), cdelta);
            ct[i] = cdiv(numer, denom);
            if sym != 0 && i % 2 != 0 { ct[i].r = -ct[i].r; ct[i].i = -ct[i].i; }
        }
        pfacr(-1, &ct, &mut rt[..nfft]);
        tr.data[..nt].copy_from_slice(&rt[nfiltby2..nfiltby2+nt]);
        tr.delrt = (tr.delrt as i32 - delrtf) as i16;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
