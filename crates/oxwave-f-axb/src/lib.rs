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

//! Clean-room port of cwp/lib/axb.c using Crout LU with scaled partial pivoting.
//! Credits: Adapted from discussions in Numerical Recipes, by Gabriel Alvarez (1995).
//! Implemented from the documented algorithm and checked against the linked C
//! library, without consulting the original C function bodies.
#![allow(non_snake_case)]

/// Decomposes a square matrix in place, recording pivots and determinant sign.
///
/// C: `void LU_decomposition(int nrows, float **matrix, int *index, float *d);`
///
/// The number of rows supplies `nrows`; each row must contain at least that
/// many entries. Columns beyond the square matrix are preserved. `index`
/// records each zero-based row exchange, and `d` is the permutation sign.
/// Pivot scores are scaled by each row's largest absolute entry; equal scores
/// select the last candidate row, as verified against the reference library.
/// Exact zero pivots are replaced with `1e-6`, as checked against C.
/// An entirely zero row prints `error, singular matrix in LU decomposition`
/// to stderr and retains infinite reciprocal scaling; factorization continues.
/// If every candidate score in a column is NaN, the previous column's pivot
/// row is reused (initially row zero), matching the reference library.
///
/// Panics for short rows or a short index vector.
pub fn LU_decomposition(matrix: &mut [Vec<f32>], index: &mut [i32], d: &mut f32) {
    let n = matrix.len();
    assert!(index.len() >= n && matrix.iter().all(|r| r.len() >= n));
    // Black-box zero-row cases establish infinite reciprocal scaling:
    // an all-zero row has a NaN pivot score, not a zero score.
    let mut scale = vec![f32::INFINITY; n];
    *d = 1.0;
    for i in 0..n {
        let mut largest = 0.0f32;
        for j in 0..n {
            let x = matrix[i][j].abs();
            if x > largest {
                largest = x;
            }
        }
        if largest == 0.0 {
            eprintln!("error, singular matrix in LU decomposition");
        } else {
            scale[i] = 1.0 / largest;
        }
    }
    // The C library remembers this row across columns, including columns
    // whose NaN scores never select a new candidate.
    let mut pivot = 0;
    for j in 0..n {
        for i in 0..j {
            let mut sum = matrix[i][j];
            for k in 0..i {
                sum -= matrix[i][k] * matrix[k][j];
            }
            matrix[i][j] = sum;
        }
        let mut largest = 0.0f32;
        for i in j..n {
            let mut sum = matrix[i][j];
            for k in 0..j {
                sum -= matrix[i][k] * matrix[k][j];
            }
            matrix[i][j] = sum;
            let candidate = scale[i] * sum.abs();
            if candidate >= largest {
                largest = candidate;
                pivot = i;
            }
        }
        if j != pivot {
            // Exchange only the active matrix entries, preserving each row's
            // allocation and any caller storage beyond the square matrix.
            for k in 0..n {
                let value = matrix[j][k];
                matrix[j][k] = matrix[pivot][k];
                matrix[pivot][k] = value;
            }
            *d = -*d;
            scale[pivot] = scale[j];
        }
        index[j] = pivot as i32;
        if matrix[j][j] == 0.0 {
            matrix[j][j] = 1.0e-6;
        }
        if j + 1 < n {
            let reciprocal = 1.0f32 / matrix[j][j];
            for i in j + 1..n {
                matrix[i][j] *= reciprocal;
            }
        }
    }
}

/// Solves a system using an LU decomposition; overwrites the right-hand side.
///
/// C: `void backward_substitution(int nrows, float **matrix, int *index, float *b);`
///
/// `matrix` and `index` must be the outputs of [`LU_decomposition`]. The row
/// count supplies `nrows`; only the first `nrows` entries of `b` are changed.
/// The factorization and permutation vector are preserved, including padding.
/// Panics for short rows or vectors, or a pivot outside the right-hand side.
pub fn backward_substitution(matrix: &[Vec<f32>], index: &[i32], b: &mut [f32]) {
    let n = matrix.len();
    assert!(index.len() >= n && b.len() >= n && matrix.iter().all(|r| r.len() >= n));
    let mut first = None;
    for i in 0..n {
        let pivot = index[i] as usize;
        let mut sum = b[pivot];
        b[pivot] = b[i];
        if let Some(start) = first {
            for j in start..i {
                sum -= matrix[i][j] * b[j];
            }
        } else if sum != 0.0 {
            first = Some(i);
        }
        b[i] = sum;
    }
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in i + 1..n {
            sum -= matrix[i][j] * b[j];
        }
        b[i] = sum / matrix[i][i];
    }
}

/// Replaces a square matrix by its inverse.
///
/// C: `void inverse_matrix(int nrows, float **matrix);`
///
/// A mathematical inverse requires a nonsingular matrix. Singular inputs
/// continue with the diagnostics and zero-pivot replacement described in
/// [`LU_decomposition`]; results may contain infinities or NaNs.
/// The number of rows supplies `nrows`, with row storage as described there.
/// Entries beyond the first `nrows` columns in each row are preserved.
pub fn inverse_matrix(matrix: &mut [Vec<f32>]) {
    let n = matrix.len();
    let mut lu = matrix.to_vec();
    let mut index = vec![0; n];
    let mut d = 1.0;
    LU_decomposition(&mut lu, &mut index, &mut d);
    for j in 0..n {
        let mut b = vec![0.0; n];
        b[j] = 1.0;
        backward_substitution(&lu, &index, &mut b);
        for i in 0..n {
            matrix[i][j] = b[i];
        }
    }
}

/// Computes the inverse of the first matrix times the second, preserving inputs.
///
/// C: `void inverse_matrix_multiply(int nrows1, float **matrix1, int ncols2, int nrows2, float **matrix2, float **out_matrix);`
///
/// `matrix1` must be square; a mathematical inverse requires it to be
/// nonsingular. Singular inputs continue as described in [`LU_decomposition`].
/// Both inputs and the output have
/// the same row count; `matrix2` supplies the column count. Solves one right-hand
/// side per column without explicitly forming an inverse. Panics for
/// incompatible dimensions.
/// Output entries beyond the column count of `matrix2` are preserved.
pub fn inverse_matrix_multiply(
    matrix1: &[Vec<f32>],
    matrix2: &[Vec<f32>],
    out_matrix: &mut [Vec<f32>],
) {
    let n = matrix1.len();
    assert_eq!(matrix2.len(), n);
    assert_eq!(out_matrix.len(), n);
    let cols = matrix2.first().map_or(0, Vec::len);
    assert!(matrix2.iter().all(|r| r.len() == cols) && out_matrix.iter().all(|r| r.len() >= cols));
    let mut lu = matrix1.to_vec();
    let mut index = vec![0; n];
    let mut d = 1.0;
    LU_decomposition(&mut lu, &mut index, &mut d);
    for j in 0..cols {
        let mut b: Vec<_> = matrix2.iter().map(|row| row[j]).collect();
        backward_substitution(&lu, &index, &mut b);
        for i in 0..n {
            out_matrix[i][j] = b[i];
        }
    }
}
