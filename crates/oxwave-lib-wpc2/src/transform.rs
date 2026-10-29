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

//! Port of comp/dwpt/2d/lib/wpcpack2.c. Author: Tong Chen, 05/25/94.
//! Modifier: Tong Chen, 07/28/94, for API.
use crate::*;

/// Wavelet filters corresponding to C `waveFILTER`.
pub struct WaveFilter<'a> {
    pub ishift: i32,
    pub filterh: &'a [f32],
    pub filterg: &'a [f32],
}
fn position(n: usize, shift: i32) -> usize {
    let mut pos = n as i32;
    while pos < shift.abs() {
        pos <<= 1;
    }
    (pos + shift) as usize
}

/// C prototype: `void wavePack1(float *x, float *y, waveFILTER *filter, int twopow, int stage, int type)`.
/// Filters must have equal, even lengths and even shift; `twopow >= 1`.
/// Like C, stage zero copies the input and then continues through the first stage.
pub fn wavePack1(
    x: &mut [f32],
    y: &mut [f32],
    filter: &WaveFilter<'_>,
    twopow: usize,
    stage: usize,
    kind: i32,
) {
    let n = 1 << twopow;
    if stage == 0 {
        for i in 0..n {
            if kind == 0 {
                y[i] = x[i];
            } else {
                x[i] = y[i];
            }
        }
    }
    let h = filter.filterh;
    let g = filter.filterg;
    let len = h.len();
    let nhalf = n >> 1;
    let mut tmp = vec![0.0f32; n];
    if kind == 0 {
        let mask = n - 1;
        let pos = position(n, filter.ishift);
        for j in 0..nhalf {
            tmp[j] = 0.0;
            y[j] = 0.0;
            for i in (1..len).step_by(2) {
                let k = (i + 2 * j - 1 + pos) & mask;
                tmp[j] += x[k] * h[i - 1] + x[k + 1] * h[i];
                y[j] += x[k] * g[i - 1] + x[k + 1] * g[i];
            }
        }
        for level in 1..stage {
            let nband = (1 << level) >> 1;
            let nmid = 1 << (twopow - level);
            let nmidh = nmid >> 1;
            let mask = nmid - 1;
            let pos = position(nmid, filter.ishift);
            for band in 0..nband {
                let first = (band * nmid) << 1;
                let mut out = first + nmid;
                for j in 0..nmidh {
                    let jj = out + j;
                    tmp[jj] = 0.0;
                    y[jj] = 0.0;
                    for i in (1..len).step_by(2) {
                        let k = ((i + 2 * j - 1 + pos) & mask) + first;
                        tmp[jj] += tmp[k] * h[i - 1] + tmp[k + 1] * h[i];
                        y[jj] += tmp[k] * g[i - 1] + tmp[k + 1] * g[i];
                    }
                }
                out += nmidh;
                for j in 0..nmidh {
                    let jj = out + j;
                    tmp[jj] = 0.0;
                    y[jj] = 0.0;
                    for i in (1..len).step_by(2) {
                        let k = ((i + 2 * j - 1 + pos) & mask) + first;
                        tmp[jj] += y[k] * g[i - 1] + y[k + 1] * g[i];
                        y[jj] += y[k] * h[i - 1] + y[k + 1] * h[i];
                    }
                }
                for j in 0..nmidh {
                    tmp[first + j] = tmp[out + j];
                    y[first + j] = y[out + j];
                }
            }
        }
        let nmid = 1 << (twopow - stage);
        let nmidh = nmid << 1;
        let nband = (1 << stage) >> 1;
        for band in 0..nband {
            let first = band * nmidh;
            let out = first + nmid;
            for i in 0..nmid {
                y[out + i] = tmp[first + i];
            }
        }
    } else {
        for i in 0..n {
            x[i] = y[i];
        }
        for level in (2..=stage).rev() {
            let nmid = 1 << (twopow - level);
            let nmidh = nmid << 1;
            let nband = (1 << level) >> 1;
            let mask = nmidh - 1;
            let pos = position(nmidh, filter.ishift);
            for band in 0..nband {
                let first = band * nmidh;
                let out = first + nmid;
                if band == (band >> 1) << 1 {
                    for i in 0..nmid {
                        tmp[first + i] = x[out + i];
                    }
                    for i in 0..nmid {
                        tmp[out + i] = x[first + i];
                    }
                } else {
                    for i in 0..nmidh {
                        tmp[first + i] = x[first + i];
                    }
                }
                for i in 0..nmidh {
                    x[first + i] = 0.0;
                }
                for j in 0..nmid {
                    let jj = first + j;
                    for i in 0..len {
                        let k = ((2 * j + i + pos) & mask) + first;
                        x[k] += tmp[jj] * g[i] + tmp[jj + nmid] * h[i];
                    }
                }
            }
        }
        for i in 0..n {
            tmp[i] = x[i];
            x[i] = 0.0;
        }
        let mask = n - 1;
        let pos = position(n, filter.ishift);
        for j in 0..nhalf {
            for i in 0..len {
                let k = (2 * j + i + pos) & mask;
                x[k] += tmp[j] * g[i] + tmp[j + nhalf] * h[i];
            }
        }
    }
}

// Private C wavePack_row, with the column loop in the original position.
fn wave_pack_row(
    x: &mut [Vec<f32>],
    y: &mut [Vec<f32>],
    filter: &WaveFilter<'_>,
    p1: usize,
    p2: usize,
    stage: usize,
    kind: i32,
) {
    let (n1, n2) = (1 << p1, 1 << p2);
    if stage == 0 {
        for i in 0..n2 {
            for j in 0..n1 {
                y[i][j] = x[i][j];
            }
        }
    }
    let h = filter.filterh;
    let g = filter.filterg;
    let len = h.len();
    let nhalf = n2 >> 1;
    let mut tmp = vec![vec![0.0f32; n1]; n2];
    if kind == 0 {
        let mask = n2 - 1;
        let pos = position(n2, filter.ishift);
        for j in 0..nhalf {
            for l in 0..n1 {
                tmp[j][l] = 0.0;
                y[j][l] = 0.0;
            }
            for i in (1..len).step_by(2) {
                let k = (i + 2 * j - 1 + pos) & mask;
                for l in 0..n1 {
                    tmp[j][l] += x[k][l] * h[i - 1] + x[k + 1][l] * h[i];
                    y[j][l] += x[k][l] * g[i - 1] + x[k + 1][l] * g[i];
                }
            }
        }
        for level in 1..stage {
            let nband = (1 << level) >> 1;
            let nmid = 1 << (p2 - level);
            let nmidh = nmid >> 1;
            let mask = nmid - 1;
            let pos = position(nmid, filter.ishift);
            for band in 0..nband {
                let first = (band * nmid) << 1;
                let mut out = first + nmid;
                for j in 0..nmidh {
                    let jj = out + j;
                    for l in 0..n1 {
                        tmp[jj][l] = 0.0;
                        y[jj][l] = 0.0;
                    }
                    for i in (1..len).step_by(2) {
                        let k = ((i + 2 * j - 1 + pos) & mask) + first;
                        for l in 0..n1 {
                            tmp[jj][l] += tmp[k][l] * h[i - 1] + tmp[k + 1][l] * h[i];
                            y[jj][l] += tmp[k][l] * g[i - 1] + tmp[k + 1][l] * g[i];
                        }
                    }
                }
                out += nmidh;
                for j in 0..nmidh {
                    let jj = out + j;
                    for l in 0..n1 {
                        tmp[jj][l] = 0.0;
                        y[jj][l] = 0.0;
                    }
                    for i in (1..len).step_by(2) {
                        let k = ((i + 2 * j - 1 + pos) & mask) + first;
                        for l in 0..n1 {
                            tmp[jj][l] += y[k][l] * g[i - 1] + y[k + 1][l] * g[i];
                            y[jj][l] += y[k][l] * h[i - 1] + y[k + 1][l] * h[i];
                        }
                    }
                }
                for j in 0..nmidh {
                    for l in 0..n1 {
                        tmp[first + j][l] = tmp[out + j][l];
                        y[first + j][l] = y[out + j][l];
                    }
                }
            }
        }
        let nmid = 1 << (p2 - stage);
        let nmidh = nmid << 1;
        let nband = (1 << stage) >> 1;
        for band in 0..nband {
            let first = band * nmidh;
            let out = first + nmid;
            for i in 0..nmid {
                for l in 0..n1 {
                    y[out + i][l] = tmp[first + i][l];
                }
            }
        }
    } else {
        for i in 0..n2 {
            for l in 0..n1 {
                x[i][l] = y[i][l];
            }
        }
        for level in (2..=stage).rev() {
            let nmid = 1 << (p2 - level);
            let nmidh = nmid << 1;
            let nband = (1 << level) >> 1;
            let mask = nmidh - 1;
            let pos = position(nmidh, filter.ishift);
            for band in 0..nband {
                let first = band * nmidh;
                let out = first + nmid;
                if band == (band >> 1) << 1 {
                    for i in 0..nmid {
                        for l in 0..n1 {
                            tmp[first + i][l] = x[out + i][l];
                        }
                    }
                    for i in 0..nmid {
                        for l in 0..n1 {
                            tmp[out + i][l] = x[first + i][l];
                        }
                    }
                } else {
                    for i in 0..nmidh {
                        for l in 0..n1 {
                            tmp[first + i][l] = x[first + i][l];
                        }
                    }
                }
                for i in 0..nmidh {
                    for l in 0..n1 {
                        x[first + i][l] = 0.0;
                    }
                }
                for j in 0..nmid {
                    let jj = first + j;
                    for i in 0..len {
                        let k = ((2 * j + i + pos) & mask) + first;
                        for l in 0..n1 {
                            x[k][l] += tmp[jj][l] * g[i] + tmp[jj + nmid][l] * h[i];
                        }
                    }
                }
            }
        }
        for i in 0..n2 {
            for l in 0..n1 {
                tmp[i][l] = x[i][l];
                x[i][l] = 0.0;
            }
        }
        let mask = n2 - 1;
        let pos = position(n2, filter.ishift);
        for j in 0..nhalf {
            for i in 0..len {
                let k = (2 * j + i + pos) & mask;
                for l in 0..n1 {
                    x[k][l] += tmp[j][l] * g[i] + tmp[j + nhalf][l] * h[i];
                }
            }
        }
    }
}

/// C prototype: `void wavePack2(float **x, float **y, void *inconf, int type)`.
pub fn wavePack2(x: &mut [Vec<f32>], y: &mut [Vec<f32>], config: &WpcConfig, kind: i32) {
    let (p1, p2) = (config.powt as usize, config.powx as usize);
    let (s1, s2) = (config.staget as usize, config.stagex as usize);
    let filter = WaveFilter {
        ishift: config.ishift,
        filterh: &config.filterh[..config.len as usize],
        filterg: &config.filterg[..config.len as usize],
    };
    let mut z = vec![vec![0.0; config.tileszt as usize]; config.tileszx as usize];
    if kind == 0 {
        for i in 0..config.tileszx as usize {
            wavePack1(&mut x[i], &mut z[i], &filter, p1, s1, kind);
        }
        wave_pack_row(&mut z, x, &filter, p1, p2, s2, kind);
        waveReorder(x, y, p1, p2, s1, s2, kind);
    } else {
        waveReorder(x, y, p1, p2, s1, s2, kind);
        for i in 0..config.tileszx as usize {
            wavePack1(&mut z[i], &mut x[i], &filter, p1, s1, kind);
        }
        wave_pack_row(x, &mut z, &filter, p1, p2, s2, kind);
    }
}
