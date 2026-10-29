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

//! Port of cwp/lib/hankel.c.
//! Credits: Dave Hale, Colorado School of Mines, 06/04/90.
//!
//! The C contract requires `nfft` to be selected by `npfar` or `npfaro`.
//! Each transform consumes and produces `nfft / 2 + 1` samples. In particular,
//! the smallest valid real FFT length is two, giving two input samples.

use oxwave_abel::{abel, abelalloc, abelfree, AbelTransformer};
use oxwave_cwp::complex::Complex;
use oxwave_pfa_fft::pfarc;

/// Workspace and precomputed Abel coefficients for a Hankel transform.
pub struct HankelTransformer {
    nfft: usize,
    gx: Vec<f32>,
    gk: Vec<Complex>,
    at: AbelTransformer,
}

/// Allocate a transformer for a valid even real FFT length of at least 2.
/// Select the length using `oxwave_pfa_fft::npfar` or `npfaro`.
/// The owned workspace can be reused for either transform and is released by
/// `hankelfree` or by dropping it. Rust initializes the scratch arrays to zero;
/// C leaves them uninitialized (see `hankel1`).
///
/// C prototype: `void *hankelalloc (int nfft);`
pub fn hankelalloc(nfft: i32) -> HankelTransformer {
    assert!(nfft >= 2 && nfft % 2 == 0);
    let nfft = nfft as usize;
    HankelTransformer {
        nfft,
        gx: vec![0.0; nfft],
        gk: vec![Complex::default(); nfft / 2 + 1],
        at: abelalloc((nfft / 2 + 1) as i32),
    }
}

/// Free a transformer by consuming it.
///
/// C prototype: `void hankelfree (void *ht);`
pub fn hankelfree(ht: HankelTransformer) {
    let HankelTransformer { gx, gk, at, .. } = ht;
    drop(gx);
    drop(gk);
    abelfree(at);
}

/// Compute the zeroth-order transform; slices need `nfft/2+1` elements.
///
/// C prototype: `void hankel0 (void *ht, float f[], float h[]);`
pub fn hankel0(ht: &mut HankelTransformer, f: &[f32], h: &mut [f32]) {
    let nfft = ht.nfft;
    let n = nfft / 2 + 1;
    assert!(f.len() >= n && h.len() >= n);
    abel(&ht.at, f, &mut ht.gx);
    for i in n..nfft {
        ht.gx[i] = ht.gx[nfft - i];
    }
    pfarc(-1, &ht.gx, &mut ht.gk);
    let scale = (1.0 / (2.0 * std::f64::consts::PI)) as f32;
    for i in 0..n {
        h[i] = scale * ht.gk[i].r;
    }
}

/// Compute the first-order transform; slices need `nfft/2+1` elements.
/// The input value `f[0]` is unused, as in C.
///
/// A fresh Rust transformer has a zero-initialized workspace. The C routine
/// reads an uninitialized `gx[0]` when called immediately after allocation,
/// then multiplies its Abel transform by zero. To obtain a defined C reference,
/// initialize that workspace with `hankel0` before calling `hankel1`.
///
/// C prototype: `void hankel1 (void *ht, float f[], float h[]);`
pub fn hankel1(ht: &mut HankelTransformer, f: &[f32], h: &mut [f32]) {
    let nfft = ht.nfft;
    let n = nfft / 2 + 1;
    assert!(f.len() >= n && h.len() >= n);
    for i in 1..n {
        ht.gx[i] = f[i] / i as f32;
    }
    // As in C, gx[0] retains its previous value; the Abel output there is
    // subsequently multiplied by zero. Rust initializes fresh storage to zero
    // rather than reading uninitialized memory on the first hankel1 call.
    // Copy because the safe Abel API uses separate input and output slices.
    let input = ht.gx[..n].to_vec();
    abel(&ht.at, &input, &mut ht.gx);
    for i in 0..n {
        ht.gx[i] *= i as f32;
    }
    for i in n..nfft {
        ht.gx[i] = -ht.gx[nfft - i];
    }
    ht.gx[nfft / 2] = 0.0;
    pfarc(-1, &ht.gx, &mut ht.gk);
    let scale = (-1.0 / (2.0 * std::f64::consts::PI)) as f32;
    for i in 0..n {
        h[i] = scale * ht.gk[i].i;
    }
}
