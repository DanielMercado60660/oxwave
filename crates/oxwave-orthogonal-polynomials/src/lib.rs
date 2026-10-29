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

//! Port of cwp/lib/orthpolynomials.c.
//! Author: Werner M. Heigl, Apache Corporation, E&P Technology, December 2006.

// SEG-BEGIN orthpolynomials.c:hermite_n_polynomial
/// Compute the generalized Hermite polynomial through the recurrence.
///
/// C prototype: `void hermite_n_polynomial(double *h, double *h0, double *h1, double *t, int nt, int n, double sigma)`.
/// Slices must have equal lengths; `h0` and `h1` are updated in place.
/// The recurrence executes at least once, including when `n` is less than two.
///
/// Note: Note that n in the function call is the order of the derivative and
/// j in the code below is the n in the recurrence relation.
///
/// Copyright (c) 2007 by the Society of Exploration Geophysicists.
/// For more information, go to http://software.seg.org/2007/0004 .
/// You must read and accept usage terms at:
/// http://software.seg.org/disclaimer.txt before use.
/// Author: Werner M. Heigl, Apache Corporation, E&P Technology, December 2006.
pub fn hermite_n_polynomial(h: &mut [f64], h0: &mut [f64], h1: &mut [f64], t: &[f64], n: i32, sigma: f64) {
    assert_eq!(h.len(), h0.len());
    assert_eq!(h.len(), h1.len());
    assert_eq!(h.len(), t.len());
    let mut j = 1i32;
    loop {
        for i in 0..h.len() {
            h[i] = (t[i] * h1[i] - (j as f64) * h0[i]) / sigma;
        }
        h0.copy_from_slice(h1);
        h1.copy_from_slice(h);
        j += 1;
        if j >= n { break; }
    }
}
// SEG-END
