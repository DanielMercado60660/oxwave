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

//! Port of `cwp/lib/differentiate.c`.
//! Author: John Stockwell, CWP, 1994.

/// Computes a two-point derivative with leading and lagging endpoint differences.
///
/// C prototype: `void differentiate(int n, float h, float *f, float *fprime);`
///
/// The input length supplies `n`; the output must be at least as long.
/// Output elements beyond `n` are left untouched. Zero sample intervals retain
/// C's IEEE floating-point division behavior (infinities and NaNs).
///
/// # Panics
///
/// Panics if fewer than two input samples are provided or the output is too
/// short. The C implementation reads outside its input for `n < 2`.
pub fn differentiate(h: f32, f: &[f32], fprime: &mut [f32]) {
    let n = f.len();
    assert!(n >= 2 && fprime.len() >= n);
    // C's `2*h` converts the integer literal to float, so this stays f32.
    let h2 = 2.0f32 * h;
    // The C source allocates an unused `temp` vector and leaks it; omit it.
    fprime[0] = (f[1] - f[0]) / h;
    for i in 1..n - 1 {
        fprime[i] = (f[i + 1] - f[i - 1]) / h2;
    }
    fprime[n - 1] = (f[n - 1] - f[n - 2]) / h;
}

/// Computes a double-precision two-point derivative with endpoint differences.
///
/// C prototype: `void ddifferentiate(int n, double h, double *f, double *fprime);`
///
/// The input length supplies `n`; the output must be at least as long.
/// Output elements beyond `n` are left untouched. Zero sample intervals retain
/// C's IEEE floating-point division behavior (infinities and NaNs).
///
/// # Panics
///
/// Panics if fewer than two input samples are provided or the output is too
/// short. The C implementation reads outside its input for `n < 2`.
pub fn ddifferentiate(h: f64, f: &[f64], fprime: &mut [f64]) {
    let n = f.len();
    assert!(n >= 2 && fprime.len() >= n);
    let h2 = 2.0f64 * h;
    // The C source allocates an unused `temp` vector and leaks it; omit it.
    fprime[0] = (f[1] - f[0]) / h;
    for i in 1..n - 1 {
        fprime[i] = (f[i + 1] - f[i - 1]) / h2;
    }
    fprime[n - 1] = (f[n - 1] - f[n - 2]) / h;
}
