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

//! Compression driver matching SU's buffer-full retry behavior.
//! Keep the transformed block when growing the output, as in wpc.c.
//! The shared crate restores its input on retry; this program needs the C
//! behavior. Reuse its transform, quantizer, encoder, and owned result types.

use oxwave_lib_wpc2::*;

pub(super) fn compress(f: &[f32], n1: usize, n2: usize, error: f32) -> Option<WpcPacked> {
    assert!(n1 > 0 && n2 > 0 && f.len() >= n1 * n2);
    let c = WpcConfig::default();
    let (nt, nx) = (c.tileszt as usize, c.tileszx as usize);
    let nblock = (n1 - 1) / nt + 1;
    let nstrip = (n2 - 1) / nx + 1;
    let n1p = nblock * nt;
    let n2l = n2 - (nstrip - 1) * nx;
    let mut strip = vec![vec![0.0; n1p]; nx];
    let mut buff = buffAlloc1(n1 * n2);
    let mut size = 0;
    let mut start = 0;
    for istrip in 0..nstrip {
        let rows = if istrip < nstrip - 1 { nx } else { n2l };
        for j in 0..rows {
            strip[j][..n1].copy_from_slice(&f[start..start + n1]);
            start += n1;
            for i in n1..n1p {
                strip[j][i] = 0.0;
            }
        }
        for j in rows..nx {
            for i in 0..n1p {
                strip[j][i] = 0.0;
            }
        }
        for it in (0..n1p).step_by(nt) {
            let mut block: Vec<Vec<f32>> =
                (0..nx).map(|i| strip[i][it..it + nt].to_vec()).collect();
            loop {
                let pos = buff.pos;
                let bsize = wpcEncoder(&mut block, error, &c, &mut buff);
                if bsize == -1 {
                    let bound = buff.code.len() * 2;
                    buffRealloc1(&mut buff, bound);
                    buff.pos = pos;
                    // C wpc.c:202-219,254-271 retries the transformed block.
                    continue;
                }
                size += bsize;
                break;
            }
            for i in 0..nx {
                strip[i][it..it + nt].copy_from_slice(&block[i]);
            }
        }
    }
    if size as usize != buff.pos {
        return None;
    }
    buff.code.truncate(size as usize);
    Some(WpcPacked {
        version: 1,
        n1: n1 as i32,
        n2: n2 as i32,
        config: c,
        code: buff.code,
    })
}
