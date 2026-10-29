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

// Copyright (c) Colorado School of Mines, 2016. All rights reserved.

//! SUST - time-frequency representation via the Stockwell transform.
//! Port of su/main/transforms/sust.c (2016/05/16).
//! Credits (from the C source): BRGM: Adnand BITRI, 20/05/2016.
//! References: Stockwell, Mansinha, and Lowe (1996), Localization of the
//! complex spectrum: The S transform, IEEE Trans. Signal Process. 44,
//! No. 4, 2957-2962; Brown, Lauzon, and Frayen (2010), A general description
//! of linear time-frequency transforms and formulation of a fast, invertible
//! transform that samples the continuous S-transform spectrum nonredundantly,
//! IEEE Trans. Signal Process. 58, No. 1, 281-290.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;
use oxwave_pfa_fft::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUST -  Outputs a time-frequency representation of seismic data via",
    "\t   the Stockwell transform (S- transform)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "    sust <stdin >stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t \t\t",
    "\tif dt is not set in header, then dt is mandatory\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tdt=(from header)\ttime sampling interval (sec)\t\t",
    "\tfmin=0\t\t\tminimum frequency of filter array (hz)\t",
    "\tfmax=NYQUIST \t\tmaximum frequency of filter array (hz)\t",
    "\tverbose=0\t\t=1 supply additional info\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes: The S transform provide a time dependent frequency distribution ",
    " of the signal. It is similar to the Gabor transform which which utilizes a Gaussian window for for spectral location. In the S transform the    ",
    " Gaussian window is scalable whith the frequency which provides a better",
    "  time freqency resolution    \t\t\t                ",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    "    suvibro | sust | suximage\t\t\t\t\t",
    "    suvibro | sust | suxmovie n1= n2= n3= \t\t\t\t",
    "     (because suxmovie scales it's amplitudes off of the first panel,  ",
    "      may have to experiment with the wclip and bclip parameters\t",
    "    suvibro | sust | supsimage | ... ( your local PostScript utility)",
    "\t\t\t\t\t\t\t\t\t",
];

const FRAC0: f64 = 0.0;
const FRAC1: f64 = 1.0;
const LOOKFAC: i32 = 2;
const PFA_MAX: i32 = 720720;

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get info from first trace. */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    // ISSEISMIC from segy.h.
    let seismic = matches!(tr.trid, 0..=17 | 109);
    let _ntr = tr.ntr;
    if !seismic {
        warn!("input is not seismic data, trid=%d", tr.trid);
    }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    if dt == 0.0 {
        err!("dt field is zero and not getparred");
    }
    let nyq = (0.5 / dt as f64) as f32;

    /* Get parameters. */
    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    let _ = verbose; // Parsed by C but otherwise unused by this program.
    let mut fmin = 0.0f32;
    if getparfloat("fmin", &mut fmin) == 0 {
        fmin = (FRAC0 * nyq as f64) as f32;
    }
    let mut fmax = 0.0f32;
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (FRAC1 * nyq as f64) as f32;
    }

    /* Set up FFT parameters. */
    let nfft = npfaro(nt, LOOKFAC * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= PFA_MAX {
        err!("Padded nt=%d -- too big", nfft);
    }
    let d1 = (1.0 / (nfft as f32 * dt) as f64) as f32;
    let nflow = (fmin / d1 + 1.0f32) as i32;
    let nfupr = ((fmax / d1 + 1.0f32) as f64 + 0.5) as i32;
    let nf = nfupr.wrapping_sub(nflow).wrapping_add(1);
    checkpars();

    /* Allocate arrays. */
    let mut rt = vec![0.0f32; nfft as usize];
    if nf < 0 {
        // ealloc2float converts nf to size_t and its allocation fails.
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    // st subtracts nfft from n+i at most once. Every spectrum access is
    // in bounds exactly when 0 <= n <= nfft. Reject C's undefined reads,
    // but preserve empty ranges and the valid endpoint n == nfft.
    if nflow <= nfupr && (nflow < 0 || nfupr > nfft) {
        err!(
            "frequency-bin range [%d, %d] outside supported range [0, %d]",
            nflow,
            nfupr,
            nfft
        );
    }
    // C allocates rows with stride nt, but st writes nfft values per row.
    // Keep that stride and the overlapping writes, with extra backing storage
    // for the final row so padded transforms do not corrupt Rust's allocator.
    let rows = if nflow == 0 && nfupr == 0 {
        nfft / 2 + 1
    } else {
        nf
    };
    let vgr_len = rows as usize * nt as usize + (nfft - nt) as usize;
    let mut vgr = Vec::new();
    if vgr.try_reserve_exact(vgr_len).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    vgr.resize(vgr_len, 0.0f32);

    loop {
        /* Load trace into rt (zero-padded). */
        rt[..nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        rt[nt as usize..].fill(0.0);

        st(nfft, nflow, nfupr, &rt, &mut vgr, nt as usize);

        let mut tracr = 0;
        for i in 0..nf as usize {
            for j in 0..nt as usize {
                tr.data[j] = vgr[i * nt as usize + j];
            }
            tr.d1 = dt;
            tracr += 1;
            tr.tracr = tracr;
            tr.gx = tr.tracr;
            tr.f2 = (fmin as f64 * 1.0) as f32;
            tr.d2 = (d1 as f64 * 1.0) as f32;
            tr.trid = 118; // AMPLITUDE
            puttr(&mut tr);
        }
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}

#[allow(non_snake_case)]
fn st(len: i32, lo: i32, mut hi: i32, data: &[f32], result: &mut [f32], stride: usize) {
    /* Check for frequency defaults. */
    if lo == 0 && hi == 0 {
        hi = len / 2;
    }

    let mut p = vec![0.0f64; len as usize];
    let mut g = vec![0.0f64; len as usize];
    let mut h = vec![Complex::default(); len as usize];
    let mut G = vec![Complex::default(); len as usize];

    /* Convert the input to complex. Also compute the mean. */
    let mut s = 0.0f64;
    for i in 0..len as usize {
        h[i].r = data[i];
        s += data[i] as f64;
    }
    s /= len as f64;
    let _mean = s;
    // Preserve C's npfaro-selected length: nt=95 and nt=96 both select 96.
    // At its radix-16 stage, m=6 has no inverse modulo 16. C proceeds with
    // mm=96 and accesses h[96], beyond the logical transform. The frozen
    // pfacc, pfamcc, and pfa2cc APIs reject this for both signs, including
    // when the strided APIs receive extra storage. Both this forward call
    // and the inverse call below need library compatibility support;
    // choosing a different FFT length would also change d1 and the bins.
    pfacc(1, &mut h);

    /* Hilbert transform: double the upper half-circle and zero the lower. */
    let mut l2 = (len + 1) / 2;
    for i in 1..l2 as usize {
        h[i].r = (h[i].r as f64 * 2.0) as f32;
        h[i].i = (h[i].i as f64 * 2.0) as f32;
    }
    l2 = len / 2 + 1;
    for i in l2 as usize..len as usize {
        h[i].r = 0.0;
        h[i].i = 0.0;
    }

    let mut n = lo;
    /* Rows contain the inverse FFT of the spectrum multiplied with the FFT
     * of scaled gaussians. */
    while n <= hi {
        for i in 0..len as usize {
            p[i] = 0.0;
            g[i] = 0.0;
        }
        g[0] = gaussa(n, 0);
        l2 = len / 2 + 1;
        for i in 1..l2 as usize {
            let v = gaussa(n, i as i32);
            g[len as usize - i] = v;
            g[i] = v;
        }
        for i in 0..len as usize {
            s = g[i];
            let mut k = n.wrapping_add(i as i32);
            if k >= len {
                k -= len;
            }
            G[i].r = (h[k as usize].r as f64 * s) as f32;
            G[i].i = (h[k as usize].i as f64 * s) as f32;
        }

        /* Inverse FFT the result to get the next row. */
        pfacc(-1, &mut G);
        for i in 0..len as usize {
            let mut ab1 = (G[i].r / len as f32) as f64;
            ab1 = ab1 * ab1;
            let mut ab2 = (G[i].i / len as f32) as f64;
            ab2 = ab2 * ab2;
            p[i] += (ab1 + ab2).sqrt();
        }
        for i in 0..len as usize {
            result[(n - lo) as usize * stride + i] = p[i] as f32;
        }
        n += 1;
    }
}

/* Fourier transform of a Gaussian (static helper in sust.c). */
fn gaussa(n: i32, m: i32) -> f64 {
    (-2.0 * PI * PI * m as f64 * m as f64 / n.wrapping_mul(n) as f64).exp()
}
