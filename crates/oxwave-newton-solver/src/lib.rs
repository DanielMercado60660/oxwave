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

//! Newton solver ported from `cwp/lib/mnewt.c`.
//! Author: Dave Hale, Colorado School of Mines, 06/06/91.

use oxwave_lu::{sgefa, sgesl};

/// Solve a nonlinear system `f(x) = 0` using Newton's method.
///
/// C prototype:
/// ```c
/// int mnewt (int maxiter, float ftol, float dxtol, int n, float *x, void *aux,
///     void (*fdfdx)(int n, float *x, float *f, float **dfdx, void *aux));
/// ```
///
/// The equation count is `x.len()` and must be positive and fit in a C `int`.
/// `x` initially contains the guess and receives the solution. The callback
/// receives `(x, f, dfdx, aux)` and must fill every residual and Jacobian entry
/// on each call. The Jacobian uses C's column-major indexing:
/// `dfdx[j][i] = df_i / dx_j`. The solver overwrites this workspace with LU
/// factors between calls. `aux` replaces C's opaque auxiliary pointer; pass
/// `&mut ()` when no auxiliary data is needed.
///
/// Returns the zero-based iteration at convergence, or `-1` on a singular
/// Jacobian or exhaustion of `maxiter`. Residual and step convergence use
/// `sum(abs(value)) <= tolerance`. Like C, step convergence returns **before**
/// applying that step to `x`, and exhaustion retains the last updated guess.
/// A zero iteration limit returns `-1`; a negative limit returns `0`, as in C.
pub fn mnewt<A, F>(
    maxiter: i32,
    ftol: f32,
    dxtol: f32,
    x: &mut [f32],
    aux: &mut A,
    mut fdfdx: F,
) -> i32
where
    A: ?Sized,
    F: FnMut(&mut [f32], &mut [f32], &mut [&mut [f32]], &mut A),
{
    let n = x.len();
    assert!(n > 0, "mnewt requires at least one equation");
    let ni = i32::try_from(n).expect("mnewt equation count must fit in a C int");

    // Allocate workspace. C aliases dx and f to the same allocation.
    let mut ipvt = vec![0; n];
    let mut f = vec![0.0f32; n];
    let mut matrix = vec![0.0f32; n.checked_mul(n).expect("Jacobian size overflow")];
    let mut dfdx: Vec<&mut [f32]> = matrix.chunks_mut(n).collect();

    let mut niter = 0;
    while niter < maxiter {
        // Compute f(x) and f'(x).
        fdfdx(x, &mut f, &mut dfdx, aux);

        // Check for ftol convergence. C's ABS preserves negative zero/NaNs.
        let mut sum = 0.0f32;
        for i in 0..n {
            sum += if f[i] < 0.0 { -f[i] } else { f[i] };
        }
        if sum <= ftol {
            break;
        }

        // Compute LU decomposition; a singular Jacobian has no solution.
        let mut info = 0;
        sgefa(&mut dfdx, ni, &mut ipvt, &mut info);
        if info >= 0 {
            return -1;
        }

        // Solve for dx in place in f, just as with C's aliased dx/f pointers.
        let factors: Vec<&[f32]> = dfdx.iter().map(|column| &column[..]).collect();
        sgesl(&factors, ni, &ipvt, &mut f, 0);

        // Check for dxtol convergence before updating the solution.
        let mut sum = 0.0f32;
        for i in 0..n {
            sum += if f[i] < 0.0 { -f[i] } else { f[i] };
        }
        if sum <= dxtol {
            break;
        }

        for i in 0..n {
            x[i] -= f[i];
        }
        niter += 1;
    }

    if niter == maxiter {
        -1
    } else {
        niter
    }
}
