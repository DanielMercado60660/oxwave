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

//! CONJTEST - compute solution for x in Ax=b using conjugate gradients.
//! Port of cwp/main/conjtest.c.
//! Credits (from the C source): Copyright (c) Colorado School of Mines, 2011.

use oxwave_core::prelude::*;
use oxwave_f_conjugate_gradient::*;

// The C program has no sdoc array or requestdoc call. Running sdoc2rs.py
// on the source reports "no sdoc array found"; tty runs solve the system too.

fn main() {
    /* compute solution for x in Ax=b using conjugate gradients */
    let testvector: [f32; 5] = [3.00, 3.00, 5.00, 7.00, 9.00];
    let testmatrix: [[f32; 4]; 5] = [
        [1.00, 1.00, 1.00, 0.00],
        [1.00, 2.00, 0.00, 0.00],
        [1.00, 3.00, 1.00, 0.00],
        [1.00, 4.00, 0.00, 1.00],
        [1.00, 5.00, 1.00, 1.00],
    ];
    let m: i32 = 5; /* size of b vector */
    let n: i32 = 4; /* size of x vector */
    let niter: i32 = 10; /* do 4 iterations */

    /* allocate space */
    let mut x = vec![0.0f32; n as usize]; /* unknown vector being solved for in Ax=b */
    let mut b = vec![0.0f32; m as usize]; /* vector b in Ax=b */
    let mut a = vec![vec![0.0f32; n as usize]; m as usize]; /* matrix A in Ax=b */

    /* read data into arrays */
    for j in 0..m as usize {
        b[j] = testvector[j];
        for i in 0..n as usize {
            a[j][i] = testmatrix[j][i];
        }
    }

    /* call conjugate gradient solver to find x in Ax = b */
    let a_rows: Vec<&[f32]> = a.iter().map(|row| row.as_slice()).collect();
    simple_conj_gradient(&mut x, &b, &a_rows, niter);

    /* print out the results */
    printf!("\n");
    printf!("Solving Ax=b by the conjugate gradient method \n");
    printf!("\n");
    printf!("The 5x4  matrix A\n");
    for j in 0..m as usize {
        for i in 0..n as usize {
            printf!("%f ", a[j][i]);
        }
        printf!("\n");
    }
    printf!("\n");

    printf!("The vector b\n");
    for j in 0..m as usize {
        printf!("%f ", b[j]);
    }

    printf!("\n");
    printf!("\n");

    printf!("The output vector x \n");
    for i in 0..n as usize {
        printf!("%f ", x[i]);
    }
    printf!("\n");

    exit(CWP_Exit());
}
