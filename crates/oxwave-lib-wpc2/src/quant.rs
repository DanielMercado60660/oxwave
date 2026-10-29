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

//! Port of comp/dwpt/2d/lib/wpcquant.c. Author: Tong Chen, 07/20/94.
//! Modifier: Tong Chen, 07/28/94, for API, without iteration.
use crate::*;

/// Owned quantization status corresponding to C `wpcQUANT`.
#[derive(Clone, Debug)]
pub struct WpcQuant {
    pub zeroflag: u8,
    pub flowflag: Vec<u8>,
    pub blockind: Vec<u8>,
    pub ave: f32,
    pub step: f32,
    pub qx: Vec<Vec<i32>>,
}
impl WpcQuant {
    pub(crate) fn new(c: &WpcConfig) -> Self {
        Self {
            zeroflag: 0,
            flowflag: vec![0; c.nblock as usize],
            blockind: vec![0; c.nblock as usize],
            ave: 0.0,
            step: 0.0,
            qx: vec![vec![0; c.lblock as usize]; c.nblock as usize],
        }
    }
}

/// C prototype: `void quantFixerror(float **x, float error, void *inconf, void *qstate)`.
pub fn quantFixerror(x: &mut [Vec<f32>], mut error: f32, c: &WpcConfig, q: &mut WpcQuant) {
    let (nblock, lblock) = (c.nblock as usize, c.lblock as usize);
    let nsz = c.tileszt.wrapping_mul(c.tileszx);
    let mut fmax = vec![0.0f32; nblock];
    let mut dev = vec![0.0f32; nblock];
    let mut indices = vec![0i32; nblock];
    let rn = (1.0 / (lblock as f32) as f64) as f32;
    let mut clip = 0.0f32;
    for b in 0..nblock - 1 {
        for i in 0..lblock {
            let v = x[b][i];
            let atmp = if v < 0.0 { -v } else { v };
            if fmax[b] < atmp {
                fmax[b] = atmp;
            }
            dev[b] += atmp * atmp;
        }
        clip += dev[b];
        dev[b] *= rn;
        dev[b] = (dev[b] as f64).sqrt() as f32;
    }
    let mut average = 0.0f32;
    for i in 0..lblock {
        average += x[nblock - 1][i];
    }
    average *= rn;
    q.ave = average;
    for i in 0..lblock {
        x[nblock - 1][i] -= average;
        let v = x[nblock - 1][i];
        let atmp = if v < 0.0 { -v } else { v };
        if fmax[nblock - 1] < atmp {
            fmax[nblock - 1] = atmp;
        }
        dev[nblock - 1] += atmp * atmp;
    }
    clip += dev[nblock - 1];
    dev[nblock - 1] *= rn;
    dev[nblock - 1] = (dev[nblock - 1] as f64).sqrt() as f32;
    clip /= nsz as f32;
    clip = (clip as f64).sqrt() as f32;
    if clip <= f32::MIN_POSITIVE {
        q.zeroflag = 1;
        return;
    }
    q.zeroflag = 0;
    error = (100.0 * error as f64) as f32;
    clip *= error;
    if clip < f32::MIN_POSITIVE {
        clip = f32::MIN_POSITIVE;
    }
    let ratio = 5.0f32;
    let bound = ratio * clip;
    let mut step_factor = (126.0 / bound as f64) as f32;
    if step_factor > f32::MAX {
        step_factor = f32::MAX;
    }
    q.step = step_factor;
    for b in 0..nblock {
        for i in 0..lblock {
            let atmp = x[b][i] * step_factor;
            q.qx[b][i] = if atmp > 0.0 {
                (atmp as f64 + 0.5) as i32
            } else {
                (atmp as f64 - 0.5) as i32
            };
        }
        let loading = fmax[b] * step_factor;
        q.flowflag[b] = if loading as f64 > 126.0 { 1 } else { 0 };
    }
    for b in 0..nblock {
        indices[b] = b as i32;
    }
    oxwave_sort::qkisort(&dev, &mut indices);
    for b in 0..nblock {
        q.blockind[b] = indices[b] as u8;
    }
}

/// C prototype: `void quantDe(float **x, void *inconf, void *qstate)`.
pub fn quantDe(x: &mut [Vec<f32>], c: &WpcConfig, q: &WpcQuant) {
    let (nblock, lblock) = (c.nblock as usize, c.lblock as usize);
    let rstep = (1.0 / q.step as f64) as f32;
    let ave = q.ave;
    if q.zeroflag == 1 {
        for b in 0..nblock {
            for i in 0..lblock {
                x[b][i] = 0.0;
            }
        }
        return;
    }
    for b in 0..nblock {
        for i in 0..lblock {
            x[b][i] = q.qx[b][i] as f32 * rstep;
        }
    }
    for i in 0..lblock {
        x[nblock - 1][i] += ave;
    }
}
