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

//! Port of tetra/lib/tetra.c.
//! Credits: CWP: Zhaobo Meng, October 1997.

/// Find the volume of a tetrahedron defined by four points.
///
/// C prototype: `float tetra_volume(float x0[3], float x1[3], float x2[3], float x3[3]);`
/// Each point must contain at least three coordinates.
pub fn tetra_volume(x0: &[f32], x1: &[f32], x2: &[f32], x3: &[f32]) -> f32 {
    let mut vv = -(x1[0] * (x2[1] * x3[2] - x2[2] * x3[1])
        - x1[1] * (x2[0] * x3[2] - x2[2] * x3[0])
        + x1[2] * (x2[0] * x3[1] - x3[0] * x2[1]));
    vv += x0[0] * (x2[1] * x3[2] - x2[2] * x3[1])
        - x0[1] * (x2[0] * x3[2] - x2[2] * x3[0])
        + x0[2] * (x2[0] * x3[1] - x3[0] * x2[1]);
    vv -= x0[0] * (x1[1] * x3[2] - x1[2] * x3[1])
        - x0[1] * (x1[0] * x3[2] - x1[2] * x3[0])
        + x0[2] * (x1[0] * x3[1] - x3[0] * x1[1]);
    vv += x0[0] * (x1[1] * x2[2] - x1[2] * x2[1])
        - x0[1] * (x1[0] * x2[2] - x1[2] * x2[0])
        + x0[2] * (x1[0] * x2[1] - x2[0] * x1[1]);
    vv = ((if vv < 0.0 { -vv } else { vv }) as f64 / 6.0) as f32;
    vv
}

/// Find the area computed by the original routine for three 3D points.
///
/// C prototype: `float area3d(float x0[3], float x1[3], float x2[3]);`
/// Each point must contain at least three coordinates.
pub fn area3d(x0: &[f32], x1: &[f32], x2: &[f32]) -> f32 {
    let d1 = (((x0[0] - x1[0]) * (x0[0] - x1[0])
        + (x0[1] - x1[1]) * (x0[1] - x1[1])
        + (x0[2] - x1[2]) * (x0[2] - x1[2])) as f64).sqrt() as f32;
    let d2 = (((x0[0] - x2[0]) * (x0[0] - x2[0])
        + (x0[1] - x2[1]) * (x0[1] - x2[1])
        + (x0[2] - x2[2]) * (x0[2] - x2[2])) as f64).sqrt() as f32;
    let d3 = (((x1[0] - x2[0]) * (x1[0] - x2[0])
        + (x1[1] - x2[1]) * (x1[1] - x2[1])
        + (x1[2] - x2[2]) * (x1[2] - x2[2])) as f64).sqrt() as f32;
    ((((d1 + d2 + d3) * (d2 + d3) * (d1 + d3) * (d1 + d2)) as f64).sqrt() / 12.0) as f32
}
