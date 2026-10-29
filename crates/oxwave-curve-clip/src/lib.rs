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

//! Port of cwp/lib/yclip.c.
//! Author: Dave Hale, Colorado School of Mines, 07/03/89.

/// Clip a uniformly sampled, piecewise linear function to `[ymin, ymax]`.
///
/// C prototype: `int yclip(int nx, float dx, float fx, float y[], float ymin,
/// float ymax, float xc[], float yc[]);`
///
/// `y` must be nonempty, `ymin <= ymax`, and both output slices must have
/// room for twice the number of input samples. Returns the number written.
pub fn yclip(
    dx: f32,
    fx: f32,
    y: &[f32],
    ymin: f32,
    ymax: f32,
    xc: &mut [f32],
    yc: &mut [f32],
) -> i32 {
    assert!(!y.is_empty());
    assert!(ymin <= ymax);
    assert!(xc.len() >= 2 * y.len() && yc.len() >= 2 * y.len());
    let mut xix = fx;
    let mut yix = y[0];
    let mut nc = 0;
    xc[nc] = xix;
    if yix < ymin {
        yc[nc] = ymin;
    } else if yix > ymax {
        yc[nc] = ymax;
    } else {
        yc[nc] = yix;
    }
    nc += 1;
    for ix in 1..y.len() {
        let yixm1 = yix;
        xix += dx;
        yix = y[ix];
        if yixm1 < ymin {
            if yix >= ymin {
                xc[nc] = xix - dx * (yix - ymin) / (yix - yixm1);
                yc[nc] = ymin;
                nc += 1;
                if yix <= ymax {
                    xc[nc] = xix;
                    yc[nc] = yix;
                } else {
                    xc[nc] = xix - dx * (yix - ymax) / (yix - yixm1);
                    yc[nc] = ymax;
                }
                nc += 1;
            }
        } else if yixm1 > ymax {
            if yix <= ymax {
                xc[nc] = xix - dx * (yix - ymax) / (yix - yixm1);
                yc[nc] = ymax;
                nc += 1;
                if yix >= ymin {
                    xc[nc] = xix;
                    yc[nc] = yix;
                } else {
                    xc[nc] = xix - dx * (yix - ymin) / (yix - yixm1);
                    yc[nc] = ymin;
                }
                nc += 1;
            }
        } else if yix < ymin {
            xc[nc] = xix - dx * (yix - ymin) / (yix - yixm1);
            yc[nc] = ymin;
            nc += 1;
        } else if yix > ymax {
            xc[nc] = xix - dx * (yix - ymax) / (yix - yixm1);
            yc[nc] = ymax;
            nc += 1;
        } else {
            xc[nc] = xix;
            yc[nc] = yix;
            nc += 1;
        }
    }
    if yix < ymin {
        xc[nc] = xix;
        yc[nc] = ymin;
        nc += 1;
    } else if yix > ymax {
        xc[nc] = xix;
        yc[nc] = ymax;
        nc += 1;
    }
    nc as i32
}
