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

//! Port of cwp/lib/yxtoxy.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

/// Compute the regularly sampled inverse of a monotonically increasing function.
///
/// C prototype: `void yxtoxy(int nx, float dx, float fx, float y[],
/// int ny, float dy, float fy, float xylo, float xyhi, float x[]);`
///
/// `y` must be nonempty and strictly increasing, and `dx` and `dy` must
/// be positive. The lengths of `y` and `x` replace C's `nx` and `ny`.
/// Values below and above the input range receive `xylo` and `xyhi`.
/// With one input sample, C returns `fx` at `y[0]` only when that output
/// sample is the last one; otherwise it returns `xyhi` there.
pub fn yxtoxy(dx: f32, fx: f32, y: &[f32], dy: f32, fy: f32, xylo: f32, xyhi: f32, x: &mut [f32]) {
    let nxi = y.len();
    let dxi = dx;
    let fxi = fx;
    let nyo = x.len();
    let dyo = dy;
    let fyo = fy;
    let fyi = y[0];

    // Loop over output y less than smallest input y.
    let mut jyo = 0;
    let mut yo = fyo;
    while jyo < nyo {
        if yo >= fyi {
            break;
        }
        x[jyo] = xylo;
        jyo += 1;
        yo += dyo;
    }
    // Loop over output y between smallest and largest input y.
    // jyo + 1 == nyo also preserves the C condition when nyo is zero.
    if jyo + 1 == nyo && yo == fyi {
        x[jyo] = fxi;
        jyo += 1;
        yo += dyo;
    }
    let mut jxi1 = 0;
    let mut jxi2 = 1;
    let mut xi1 = fxi;
    while jxi2 < nxi && jyo < nyo {
        let yi1 = y[jxi1];
        let yi2 = y[jxi2];
        if yi1 <= yo && yo <= yi2 {
            x[jyo] = xi1 + dxi * (yo - yi1) / (yi2 - yi1);
            jyo += 1;
            yo += dyo;
        } else {
            jxi1 += 1;
            jxi2 += 1;
            xi1 += dxi;
        }
    }
    // Loop over output y greater than largest input y.
    while jyo < nyo {
        x[jyo] = xyhi;
        jyo += 1;
    }
}
