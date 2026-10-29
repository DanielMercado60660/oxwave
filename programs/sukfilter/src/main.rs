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

//! SUKFILTER - radial wavenumber filtering.
//! Port of su/main/filters/sukfilter.c (Revision 1.7).
//! Credits: CWP: John Stockwell, June 1997.

use oxwave_core::prelude::*;
use oxwave_cwp::{complex::crmul, cwp::nint};
use oxwave_pfa_fft::{npfa, npfaro, pfa2cc, pfa2cr, pfa2rc, Complex};
use oxwave_polar_grid::polartorect;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUKFILTER - radially symmetric K-domain, sin^2-tapered, polygonal\t",
    "\t\t  filter\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "     sukfilter <infile >outfile [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " k=val1,val2,...\tarray of K filter wavenumbers\t\t\t",
    " amps=a1,a2,...\tarray of K filter amplitudes\t\t\t",
    " d1=tr.d1 or 1.0\tsampling interval in first (fast) dimension\t",
    " d2=tr.d1 or 1.0\tsampling interval in second (slow) dimension\t",
    "\t\t\t\t\t\t\t\t\t",
    " Defaults:\t\t\t\t\t\t\t\t",
    " k=.10*(nyq),.15*(nyq),.45*(nyq),.50*(nyq)\t\t\t\t",
    " amps=0.,1.,...,1.,0.  trapezoid-like bandpass filter\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The nyquist wavenumbers, nyq=sqrt(nyq1^2 + nyq2^2) is  computed\t",
    " internally.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The filter is assumed to be symmetric, to yield real output.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Because the data are assumed to be purely spatial (i.e. non-seismic), ",
    " the data are assumed to have trace id (30), corresponding to (z,x) data",
    "\t\t\t\t\t\t\t\t\t",
    " The relation: w = 2 pi F is well known for frequency, but there\t",
    " doesn't seem to be a commonly used letter corresponding to F for the\t",
    " spatial conjugate transform variables.  We use K1 and K2 for this.\t",
    " More specifically we assume a phase:\t\t\t\t\t",
    "\t\t-i(k1 x1 + k2 x2) = -2 pi i(K1 x1 + K2 x2).\t\t",
    " and K1, K2 define our respective wavenumbers.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

// su/include/segy.h defines 130; the selfdoc retains the historical value 30.
const TRID_DEPTH: i16 = 130;

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    if tr.trid != TRID_DEPTH {
        warn!("tr.trid = %d", tr.trid);
    }
    let nx1 = tr.ns as i32;
    let mut dx1 = 0.0f32;
    let mut dx2 = 0.0f32;
    if getparfloat("d1", &mut dx1) == 0 {
        if tr.d1 != 0.0 {
            dx1 = tr.d1;
        } else {
            dx1 = 1.0;
            warn!("tr.d1 not set, assuming d1=1.0");
        }
    }
    if getparfloat("d2", &mut dx2) == 0 {
        if tr.d2 != 0.0 {
            dx2 = tr.d2;
        } else {
            dx2 = 1.0;
            warn!("tr.d2 not set, assuming d2=1.0");
        }
    }
    let dx = nint(((dx1 * dx1 + dx2 * dx2) as f64).sqrt()) as f32;
    // C parses verbose but never uses it, and never calls checkpars().
    let mut verbose = 1;
    getparint("verbose", &mut verbose);
    let nyq = (0.5 / dx as f64) as f32;
    // Store traces in tmpfile while getting a count.
    let tracefp = etmpfile();
    let hfp = etmpfile();
    let mut nx2 = 0;
    loop {
        nx2 += 1;
        efwrite(std::slice::from_ref(&tr.hdr), &hfp);
        efwrite(&tr.data[..nx1 as usize], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    let nx1fft = npfaro(nx1, 2 * nx1);
    let nx2fft = npfa(nx2);
    if nx1fft >= SU_NFLTS as i32 || nx1fft >= 720720 {
        err!("Padded nx1=%d--too big", nx1fft);
    }
    if nx2fft >= 720720 {
        err!("Padded nx2=%d--too big", nx2fft);
    }
    // Determine number of wavenumbers in K1 and K2.
    let nk1 = nx1fft / 2 + 1;
    let nk2 = nx2fft / 2 + 1;
    let nk = nint((nk1.wrapping_mul(nk1).wrapping_add(nk2.wrapping_mul(nk2)) as f64).sqrt());
    let nxfft = 2 * (nk - 1); // Faked for polygonal filter in k.

    let mut npoly = countparval("k");
    let mut k = Vec::<f32>::new();
    if npoly != 0 {
        getparfloat("k", &mut k);
    } else {
        npoly = 4;
        k = vec![
            (0.10 * nyq as f64) as f32,
            (0.15 * nyq as f64) as f32,
            (0.45 * nyq as f64) as f32,
            (0.50 * nyq as f64) as f32,
        ];
    }
    if npoly < 2 {
        warn!("Only %d value defining filter", npoly);
    }
    for ik in 0..npoly as usize - 1 {
        if k[ik] < 0.0 || k[ik] > k[ik + 1] {
            err!("Bad filter parameters");
        }
    }
    let mut namps = countparval("amps");
    let mut amps = Vec::<f32>::new();
    if namps != 0 {
        getparfloat("amps", &mut amps);
    } else {
        namps = npoly;
        amps = vec![1.0; namps as usize];
        amps[0] = 0.0;
        amps[namps as usize - 1] = 0.0;
    }
    if namps != npoly {
        err!("number of k values must = number of amps values");
    }
    let mut icount = 0i32;
    for &amp in &amps {
        icount = (icount as f32 + amp) as i32;
        if amp < 0.0 {
            err!("amps values must be positive");
        }
    }
    if icount == 0 {
        err!("All amps values are zero");
    }
    icount = 0;
    for i in 0..namps as usize - 1 {
        if amps[i] != amps[i + 1] {
            icount += 1;
        }
    }
    if icount == 0 {
        warn!("All amps values are the same");
    }
    // Allocate and zero arrays; karray is filled by polartorect.
    let mut rt = vec![0.0f32; nx1fft as usize * nx2fft as usize];
    let mut ct = vec![Complex::default(); nk1 as usize * nx2fft as usize];
    let mut karray = vec![0.0f32; nx1fft as usize * nx2fft as usize];
    let mut kfilter = vec![0.0f32; nx1fft as usize * nx2fft as usize];
    let mut kfilt = vec![0.0f32; nk as usize];
    // ealloc2float(1,nK): the angular dimension is exactly one.
    let mut kphi = vec![0.0f32; nk as usize];

    // Build filter and k,phi filter.
    polygonal_filter(&k, &amps, nxfft, dx, &mut kfilt);
    for _iphi in 0..1 {
        for ik in 0..nk as usize {
            kphi[ik] = kfilt[ik] * kfilt[ik];
        }
    }

    // Convert polar to rectangular coordinates. The library infers dimensions
    // from row slices; keep the C allocation's single angle and zero da.
    {
        let kphi_rows: Vec<&[f32]> = kphi.chunks_exact(1).collect();
        let mut karray_rows: Vec<&mut [f32]> = karray.chunks_exact_mut(nx1fft as usize).collect();
        polartorect(
            0.0,
            0.0,
            1.0,
            0.0,
            &kphi_rows,
            1.0,
            -nk1 as f32,
            1.0,
            -nk2 as f32,
            &mut karray_rows,
        );
    }

    // Positive k1, positive k2.
    for ik2 in 0..nk2 as usize {
        for ik1 in 0..nk1 as usize {
            kfilter[ik2 * nx1fft as usize + ik1] =
                karray[(nk2 as usize - 1 - ik2) * nx1fft as usize + nk1 as usize - 1 - ik1];
        }
    }
    // Positive k1, negative k2.
    for ik2 in nk2 as usize..nx2fft as usize {
        for ik1 in 0..nk1 as usize {
            kfilter[ik2 * nx1fft as usize + ik1] =
                karray[(ik2 - nk2 as usize) * nx1fft as usize + nk1 as usize - 1 - ik1];
        }
    }

    // Load traces into FFT arrays and close tmpfile.
    erewind(&tracefp);
    for ix2 in 0..nx2 as usize {
        let start = ix2 * nx1fft as usize;
        efread(&mut rt[start..start + nx1 as usize], &tracefp);
    }
    efclose(&tracefp);

    // Fourier transform dimensions 1 and 2.
    pfa2rc(-1, 1, nx1fft, nx2, &rt, &mut ct);
    pfa2cc(-1, 2, nk1, nx2fft, &mut ct);

    // Apply filter.
    for ik2 in 0..nx2fft as usize {
        for ik1 in 0..nk1 as usize {
            let index = ik2 * nk1 as usize + ik1;
            ct[index] = crmul(ct[index], kfilter[ik2 * nx1fft as usize + ik1]);
        }
    }

    // Inverse Fourier transform dimensions 2 and 1.
    pfa2cc(1, 2, nk1, nx2fft, &mut ct);
    pfa2cr(1, 1, nx1fft, nx2, &ct, &mut rt);

    erewind(&hfp);
    // Output filtered traces.
    for ix2 in 0..nx2 as usize {
        efread(std::slice::from_mut(&mut tr.hdr), &hfp);
        for ix1 in 0..nx1 as usize {
            tr.data[ix1] = rt[ix2 * nx1fft as usize + ix1];
        }
        puttr(&mut tr);
    }
    efclose(&hfp);
    exit(CWP_Exit());
}

// This filter is defined in sukfilter.c itself, rather than in a CWP library.
// Author: CWP: John Stockwell, 1992.
// Its sine taper is squared when building kphi, as in the original program.
fn polygonal_filter(f: &[f32], amps: &[f32], nfft: i32, dt: f32, filter: &mut [f32]) {
    let nf = nfft / 2 + 1;
    let nfm1 = nf - 1;
    let onfft = (1.0f64 / nfft as f64) as f32;
    let df = onfft / dt;
    let mut intfr = vec![0i32; f.len()];
    for ifs in 0..f.len() {
        intfr[ifs] = nint((f[ifs] / df) as f64);
        if intfr[ifs] > nfm1 {
            intfr[ifs] = nfm1;
        }
    }
    for icount in 0..intfr[0] {
        filter[icount as usize] = amps[0] * onfft;
    }
    let mut taper = 0;
    for ifs in 0..f.len() - 1 {
        if amps[ifs] < amps[ifs + 1] {
            taper += 1;
            for icount in intfr[ifs]..=intfr[ifs + 1] {
                let c = (1.57079632679490f64 / (intfr[ifs + 1] - intfr[ifs] + 2) as f64) as f32;
                let s = ((c * (icount - intfr[ifs] + 1) as f32) as f64).sin() as f32;
                let adiff = amps[ifs + 1] - amps[ifs];
                filter[icount as usize] = (amps[ifs] + adiff * s) * onfft;
            }
        } else if amps[ifs] > amps[ifs + 1] {
            taper += 1;
            for icount in intfr[ifs]..=intfr[ifs + 1] {
                let c = (1.57079632679490f64 / (intfr[ifs + 1] - intfr[ifs] + 2) as f64) as f32;
                let s = ((c * (intfr[ifs + 1] - icount + 1) as f32) as f64).sin() as f32;
                let adiff = amps[ifs] - amps[ifs + 1];
                filter[icount as usize] = (amps[ifs + 1] + adiff * s) * onfft;
            }
        } else if taper == 0 {
            for icount in intfr[ifs]..=intfr[ifs + 1] {
                filter[icount as usize] = amps[ifs] * onfft;
            }
        } else {
            for icount in intfr[ifs] + 1..=intfr[ifs + 1] {
                filter[icount as usize] = amps[ifs] * onfft;
            }
        }
    }
    for icount in intfr[f.len() - 1] + 1..nf {
        filter[icount as usize] = amps[f.len() - 1] * onfft;
    }
}
