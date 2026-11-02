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

//! Port of comp/dwpt/2d/lib/wpc.c and wpcendec.c.
//! Author: Tong Chen, 07/27/94 (WPC), 07/28/94 (block wrappers).
use crate::*;
use std::io::{Read, Write};

/// Configuration corresponding exactly to C `wpcCONFIG`.
#[repr(C)]
#[derive(Clone, Debug, PartialEq)]
pub struct WpcConfig {
    pub powt: i32,
    pub powx: i32,
    pub tileszt: i32,
    pub tileszx: i32,
    pub nblock: i32,
    pub lblock: i32,
    pub codenblock: i32,
    pub codelblock: i32,
    pub staget: i32,
    pub stagex: i32,
    pub ishift: i32,
    pub len: i32,
    pub filterh: [f32; 12],
    pub filterg: [f32; 12],
}
impl Default for WpcConfig {
    fn default() -> Self {
        let h = [
            0.016387336464,
            -0.041464936782,
            -0.067372554722,
            0.386110066823,
            0.812723635450,
            0.416998396426,
            -0.076488599079,
            -0.059434418647,
            0.023680171946,
            0.005611434819,
            -0.001823208871,
            -0.000720549445,
        ];
        let mut g = [0.0; 12];
        for i in 0..12 {
            g[i] = h[11 - i];
        }
        for i in (1..12).step_by(2) {
            g[i] = -g[i];
        }
        Self {
            powt: 8,
            powx: 8,
            tileszt: 256,
            tileszx: 256,
            nblock: 64,
            lblock: 1024,
            codenblock: 4,
            codelblock: 16,
            staget: 3,
            stagex: 3,
            ishift: -4,
            len: 12,
            filterh: h,
            filterg: g,
        }
    }
}

/// Owned C `wpcPTR` object. Size is the length of `code`.
#[derive(Clone, Debug)]
pub struct WpcPacked {
    pub version: i32,
    pub n1: i32,
    pub n2: i32,
    pub config: WpcConfig,
    pub code: Vec<u8>,
}

/// C prototype: `int wpcEncoder(float **block, float error, void *inconf, void *outbuff, void *interblock)`.
pub fn wpcEncoder(
    block: &mut [Vec<f32>],
    error: f32,
    c: &WpcConfig,
    output: &mut WpcBuffer,
) -> i32 {
    let mut x = vec![vec![0.0; c.lblock as usize]; c.nblock as usize];
    let mut q = WpcQuant::new(c);
    wavePack2(block, &mut x, c, 0);
    quantFixerror(&mut x, error, c, &mut q);
    codeEncoder(c, &q, output)
}

/// C prototype: `int wpcDecoder(float **block, void *inconf, void *outbuff, void *interblock)`.
/// Uses the flag and average restored by `codeDecoder`, including zero tiles.
pub fn wpcDecoder(block: &mut [Vec<f32>], c: &WpcConfig, input: &mut WpcBuffer) -> i32 {
    let mut x = vec![vec![0.0; c.lblock as usize]; c.nblock as usize];
    let mut q = WpcQuant::new(c);
    let retval = codeDecoder(c, &mut q, input);
    if retval != 0 {
        quantDe(&mut x, c, &q);
        wavePack2(block, &mut x, c, 1);
    }
    if retval != 0 {
        1
    } else {
        0
    }
}

/// C prototype: `void *wpcCompress(float *f, int n1, int n2, float error, int *nbytes)`.
/// Returns the owned result; `nbytes` becomes `wpcSize(&result)`.
/// C retries after a full buffer using the already transformed block; OxWave
/// retries from the untransformed block (DEVIATIONS.md). Rust keeps ownership when realloc moves.
pub fn wpcCompress(f: &[f32], n1: usize, n2: usize, error: f32) -> Option<WpcPacked> {
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
            // OxWave deviation from SU (DEVIATIONS.md): wpcEncoder transforms the block in
            // place, and SU retries a full buffer with the already transformed block, which
            // can push the round-trip error past 100%. Keep the input and retry from it.
            let original = block.clone();
            loop {
                let pos = buff.pos;
                let bsize = wpcEncoder(&mut block, error, &c, &mut buff);
                if bsize == -1 {
                    let bound = buff.code.len() * 2;
                    buffRealloc1(&mut buff, bound);
                    buff.pos = pos;
                    block.clone_from(&original);
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

/// C prototype: `int wpcUncompress(void *w, float *f)`.
pub fn wpcUncompress(w: &WpcPacked, f: &mut [f32]) -> i32 {
    if w.n1 <= 0 || w.n2 <= 0 || w.code.is_empty() {
        return 0;
    }
    let (n1, n2) = (w.n1 as usize, w.n2 as usize);
    let c = &w.config;
    let (nt, nx) = (c.tileszt as usize, c.tileszx as usize);
    let nblock = (n1 - 1) / nt + 1;
    let nstrip = (n2 - 1) / nx + 1;
    let n1p = nblock * nt;
    let n2l = n2 - (nstrip - 1) * nx;
    let mut strip = vec![vec![0.0; n1p]; nx];
    let mut buff = WpcBuffer {
        code: w.code.clone(),
        pos: 0,
    };
    let mut start = 0;
    for istrip in 0..nstrip {
        for it in (0..n1p).step_by(nt) {
            let mut block: Vec<Vec<f32>> =
                (0..nx).map(|i| strip[i][it..it + nt].to_vec()).collect();
            let retval = wpcDecoder(&mut block, c, &mut buff);
            for i in 0..nx {
                strip[i][it..it + nt].copy_from_slice(&block[i]);
            }
            if retval == 0 {
                return 0;
            }
        }
        let rows = if istrip < nstrip - 1 { nx } else { n2l };
        for j in 0..rows {
            f[start..start + n1].copy_from_slice(&strip[j][..n1]);
            start += n1;
        }
    }
    1
}

/// C prototype: `int wpcN1(void *wpc)`.
pub fn wpcN1(w: &WpcPacked) -> i32 {
    w.n1
}
/// C prototype: `int wpcN2(void *wpc)`.
pub fn wpcN2(w: &WpcPacked) -> i32 {
    w.n2
}
/// C prototype: `int wpcSize(void *wpc)`.
pub fn wpcSize(w: &WpcPacked) -> i32 {
    w.code.len() as i32
}
/// C prototype: `void wpcFree(void *wpc)`.
pub fn wpcFree(_w: WpcPacked) {}

/// C prototype: `void wpcWrite(void *w, FILE *fp)`.
/// Writes native-endian fields, matching C's machine-dependent format.
pub fn wpcWrite(w: &WpcPacked, fp: &mut impl Write) -> std::io::Result<()> {
    for v in [w.version, w.n1, w.n2, w.code.len() as i32] {
        fp.write_all(&v.to_ne_bytes())?;
    }
    let c = &w.config;
    for v in [
        c.powt,
        c.powx,
        c.tileszt,
        c.tileszx,
        c.nblock,
        c.lblock,
        c.codenblock,
        c.codelblock,
        c.staget,
        c.stagex,
        c.ishift,
        c.len,
    ] {
        fp.write_all(&v.to_ne_bytes())?;
    }
    for filter in [&c.filterh, &c.filterg] {
        for v in filter {
            fp.write_all(&v.to_ne_bytes())?;
        }
    }
    fp.write_all(&w.code)
}

/// C prototype: `void *wpcRead(FILE *fp)`.
/// C ignores short config/code reads. Missing bytes are zero-filled here;
/// defined C decoding must not consume those uninitialized bytes.
pub fn wpcRead(fp: &mut impl Read) -> std::io::Result<WpcPacked> {
    fn available(fp: &mut impl Read, mut b: &mut [u8]) -> std::io::Result<usize> {
        let len = b.len();
        while !b.is_empty() {
            match fp.read(b) {
                Ok(0) => break,
                Ok(n) => b = &mut b[n..],
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            }
        }
        Ok(len - b.len())
    }
    fn int(fp: &mut impl Read) -> std::io::Result<i32> {
        let mut b = [0; 4];
        fp.read_exact(&mut b)?;
        Ok(i32::from_ne_bytes(b))
    }
    let version = int(fp)?;
    let n1 = int(fp)?;
    let n2 = int(fp)?;
    let size = int(fp)?;
    if size < 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "negative WPC size",
        ));
    }
    let mut config_bytes = [0; 12 * 4 + 24 * 4];
    let config_len = available(fp, &mut config_bytes)?;
    if config_len < config_bytes.len() && n1 > 0 && n2 > 0 && size > 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "truncated WPC configuration",
        ));
    }
    let mut config_input = &config_bytes[..];
    let mut c = WpcConfig::default();
    c.powt = int(&mut config_input)?;
    c.powx = int(&mut config_input)?;
    c.tileszt = int(&mut config_input)?;
    c.tileszx = int(&mut config_input)?;
    c.nblock = int(&mut config_input)?;
    c.lblock = int(&mut config_input)?;
    c.codenblock = int(&mut config_input)?;
    c.codelblock = int(&mut config_input)?;
    c.staget = int(&mut config_input)?;
    c.stagex = int(&mut config_input)?;
    c.ishift = int(&mut config_input)?;
    c.len = int(&mut config_input)?;
    for filter in [&mut c.filterh, &mut c.filterg] {
        for v in filter {
            let mut b = [0; 4];
            config_input.read_exact(&mut b)?;
            *v = f32::from_ne_bytes(b);
        }
    }
    let mut code = vec![0; size as usize];
    available(fp, &mut code)?;
    Ok(WpcPacked {
        version,
        n1,
        n2,
        config: c,
        code,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_declared_code_buffer_can_contain_complete_tile() {
        let mut samples = vec![0.0; 256 * 256];
        samples[2 * 256 + 6] = 1.0;
        let packed = wpcCompress(&samples, 256, 256, 0.01).unwrap();
        let mut bytes = Vec::new();
        wpcWrite(&packed, &mut bytes).unwrap();
        let size = packed.code.len() + 17;
        bytes[12..16].copy_from_slice(&(size as i32).to_ne_bytes());
        let read = wpcRead(&mut bytes.as_slice()).unwrap();
        assert_eq!(read.code.len(), size);
        assert_eq!(&read.code[..packed.code.len()], packed.code);
        let mut expected = vec![0.0; samples.len()];
        let mut actual = expected.clone();
        assert_eq!(wpcUncompress(&packed, &mut expected), 1);
        assert_eq!(wpcUncompress(&read, &mut actual), 1);
        assert_eq!(actual, expected);
    }

    #[test]
    fn empty_data_does_not_require_unused_config() {
        for (n1, n2, size) in [(0i32, 256i32, 0i32), (256, 0, 0), (256, 256, 0)] {
            let bytes: Vec<u8> = [1, n1, n2, size]
                .iter()
                .flat_map(|v| v.to_ne_bytes())
                .collect();
            let packed = wpcRead(&mut bytes.as_slice()).unwrap();
            assert_eq!(packed.n1, n1);
            assert_eq!(packed.n2, n2);
            assert_eq!(wpcUncompress(&packed, &mut []), 0);
        }
    }
}
