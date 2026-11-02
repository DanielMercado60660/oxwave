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

//! Ports of comp/dct/lib/buffalloc.c, dct1.c, dct2.c, dctalloc.c,
//! huffman.c, lct1.c, pcoding.c, quant.c and rle.c.
//! Credits from the C sources:
//! Author: Tong Chen (buffalloc.c, rle.c).
//! Author: Tong Chen, 03/16/95 (dct1.c, dct2.c).
//! Modifier: Tong Chen, 06/01/95, use the pre-generated table.
//! Author: Tong Chen, 06/01/95 (dctalloc.c).
//! Author: Tong Chen, 08/01/94 (huffman.c).
//! The algorithm is based on that in Mark Nelson's "Data Compression" book.
//! Author: Tong Chen, 05/22/95 (lct1.c).
//! Author: Tong Chen, 05/31/95, modified from en/decoding (pcoding.c).
//! Author: Tong Chen, 07/15/94 (quant.c).
//! Modifier: Tong Chen, 05/30/95, for unified interface.
#![allow(non_snake_case)]

mod prediction;
mod wave;
pub use prediction::{lpred, LpFilter};
pub use wave::{
    waveFreefilter, waveGetfilter, wavePack_1, wavePack_2, waveTrans_1, waveTrans_2, WaveFilter,
};

/// C prototype: `float **dctivAlloc(int n);`
/// The transform length must be positive, as required by C's alloc2float.
pub fn dctivAlloc(n: usize) -> Vec<Vec<f32>> {
    assert!(n > 0);
    let c1 = (1.0 / (0.5 * n as f64).sqrt()) as f32;
    let c2 = (std::f64::consts::PI / n as f64) as f32;
    let mut c = vec![vec![0.0; n]; n];
    for j in 0..n {
        for i in 0..n {
            c[j][i] = (c1 as f64 * ((i as f64 + 0.5) * (j as f64 + 0.5) * c2 as f64).cos()) as f32;
        }
    }
    c
}
/// C prototype: `float **dctAlloc(int n);`
/// The transform length must be positive, as required by C's alloc2float.
pub fn dctAlloc(n: usize) -> Vec<Vec<f32>> {
    assert!(n > 0);
    let c1 = (1.0 / (n as f64).sqrt()) as f32;
    let c2 = (1.0 / (0.5 * n as f64).sqrt()) as f32;
    let c3 = (0.5 * std::f64::consts::PI / n as f64) as f32;
    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        c[0][i] = c1;
    }
    for j in 1..n {
        for i in 0..n {
            c[j][i] = (c2 as f64 * ((2.0 * i as f64 + 1.0) * j as f64 * c3 as f64).cos()) as f32;
        }
    }
    c
}
/// C prototype: `void dctiv(float *x, int n, float **c);`
pub fn dctiv(x: &mut [f32], c: &[Vec<f32>]) {
    dct(x, c, 0);
}
/// C prototype: `void dct(float *x, int n, float **c, int type);`
pub fn dct(x: &mut [f32], c: &[Vec<f32>], type_: i32) {
    let n = x.len();
    let mut tmp = vec![0.0f32; n];
    if type_ == 0 {
        for i in 0..n {
            for j in 0..n {
                tmp[i] += x[j] * c[i][j];
            }
        }
    } else {
        for j in 0..n {
            for i in 0..n {
                tmp[i] += x[j] * c[j][i];
            }
        }
    }
    x.copy_from_slice(&tmp);
}
/// C prototype: `void dctiv_2(float **x, int n1, int n2, float **c1, float **c2);`
pub fn dctiv_2(x: &mut [Vec<f32>], c1: &[Vec<f32>], c2: &[Vec<f32>]) {
    dct_2(x, c1, c2, 0);
}
/// C prototype: `void dct_2(float **x, int n1, int n2, float **c1, float **c2, int type);`
pub fn dct_2(x: &mut [Vec<f32>], c1: &[Vec<f32>], c2: &[Vec<f32>], type_: i32) {
    let n2 = x.len();
    if n2 == 0 {
        return;
    }
    let n1 = x[0].len();
    for row in x.iter_mut() {
        assert_eq!(row.len(), n1);
        dct(row, c1, type_);
    }
    let mut tmp = vec![vec![0.0f32; n1]; n2];
    if type_ == 0 {
        for i in 0..n2 {
            for j in 0..n2 {
                for k in 0..n1 {
                    tmp[i][k] += x[j][k] * c2[i][j];
                }
            }
        }
    } else {
        for j in 0..n2 {
            for i in 0..n2 {
                for k in 0..n1 {
                    tmp[i][k] += x[j][k] * c2[j][i];
                }
            }
        }
    }
    for i in 0..n2 {
        x[i].copy_from_slice(&tmp[i]);
    }
}
/// C prototype: `void uniQuant(float *x, int n, float error, float *ave, float *step, int *qx);`
pub fn uniQuant(x: &[f32], error: f32, ave: &mut f32, step: &mut f32, qx: &mut [i32]) {
    let n = x.len();
    assert!(qx.len() >= n);
    let rn = (1.0 / n as f64) as f32;
    let mut lave = 0.0f32;
    for &v in x {
        lave += v;
    }
    lave *= rn;
    let mut g = vec![0.0; n];
    let dev;
    if *step < 0.0 {
        let mut d = 0.0f32;
        for i in 0..n {
            g[i] = x[i] - lave;
            let a = if g[i] < 0.0 { -g[i] } else { g[i] };
            d += a * a;
        }
        d *= rn;
        dev = (d as f64).sqrt() as f32;
    } else {
        for i in 0..n {
            g[i] = x[i] - lave;
        }
        dev = *step;
    }
    let mut lstep = ((dev * error) as f64 * 3.4641) as f32;
    lstep = (1.0 / lstep as f64) as f32;
    oxwave_core::fprintf!(oxwave_core::prelude::stderr(), "lstep=%f\n", lstep);
    for i in 0..n {
        let a = g[i] * lstep;
        qx[i] = if a > 0.0 {
            (a as f64 + 0.5) as i32
        } else {
            (a as f64 - 0.5) as i32
        };
    }
    oxwave_core::fprintf!(oxwave_core::prelude::stderr(), "after quantization\n");
    *ave = lave;
    *step = lstep;
}
/// C prototype: `void uniDequant(float *x, int n, float ave, float step, int *qx);`
pub fn uniDequant(x: &mut [f32], ave: f32, step: f32, qx: &[i32]) {
    let rstep = (1.0 / step as f64) as f32;
    for i in 0..x.len() {
        x[i] = ave + qx[i] as f32 * rstep;
    }
}

// Ports of buffalloc.c, rle.c and huffman.c, author Tong Chen.
/// Memory buffer used by the coding functions.
pub use oxwave_lib_wpc2::WpcBuffer;
/// C prototype: `memBUFF *buffAlloc1(int n);`
pub fn buffAlloc1(n: usize) -> WpcBuffer {
    oxwave_lib_wpc2::buffAlloc1(n)
}
/// C prototype: `memBUFF **buffAlloc2(int n1, int n2);`
pub fn buffAlloc2(n1: usize, n2: usize) -> Vec<WpcBuffer> {
    oxwave_lib_wpc2::buffAlloc2(n1, n2)
}
/// C prototype: `void buffFree1(memBUFF *buff);`
pub fn buffFree1(buff: WpcBuffer) {
    oxwave_lib_wpc2::buffFree1(buff)
}
/// C prototype: `void buffFree2(memBUFF **buff, int n1);`
pub fn buffFree2(buff: Vec<WpcBuffer>) {
    oxwave_lib_wpc2::buffFree2(buff)
}
/// C prototype: `void buffRealloc1(memBUFF *buff, int n);`
pub fn buffRealloc1(buff: &mut WpcBuffer, n: usize) {
    oxwave_lib_wpc2::buffRealloc1(buff, n)
}
/// C prototype: `void codeSilence(void *inb, void *outb);`
pub fn codeSilence(input: &mut WpcBuffer, output: &mut WpcBuffer) {
    oxwave_lib_wpc2::codeSilence(input, output)
}
/// C prototype: `void codeDesilence(void *inb, void *outb);`
pub fn codeDesilence(input: &mut WpcBuffer, output: &mut WpcBuffer) {
    oxwave_lib_wpc2::codeDesilence(input, output)
}
/// C prototype: `int huffCompress(void *inb, void *outb);`
pub fn huffCompress(input: &mut WpcBuffer, output: &mut WpcBuffer) -> i32 {
    use oxwave_lib_wpc2::*;
    let mut bit = WpcBitBuffer {
        buff: output,
        mask: 0x80,
        rack: 0,
    };
    let mut counts = [0; 256];
    let mut nodes = [HuffNode::default(); 514];
    let mut codes = [HuffCode::default(); 257];
    huffCountbytes(input, &mut counts);
    huffScalecounts(&mut counts, &mut nodes);
    let mut nbytes = huffOutputcounts(&mut bit, &nodes);
    if nbytes == -1 {
        return -1;
    }
    let root = huffBuildtree(&mut nodes);
    huffConverttreetocode(&nodes, &mut codes, 0, 0, root);
    let retval = huffEncoder(input, &mut bit, &codes);
    // This diagnostic occurs in comp/dct/lib/huffman.c only, not wpchuff.c.
    oxwave_core::fprintf!(oxwave_core::prelude::stderr(), "retval=%d\n", retval);
    if retval == -1 {
        return -1;
    }
    nbytes += retval;
    if bit.mask != 0x80 {
        if bit.buff.pos == bit.buff.code.len() {
            return -1;
        }
        bit.buff.code[bit.buff.pos] = bit.rack;
        bit.buff.pos += 1;
        nbytes += 1;
    }
    nbytes
}
/// C prototype: `int huffDecompress(void *inb, void *outb);`
pub fn huffDecompress(input: &mut WpcBuffer, output: &mut WpcBuffer) -> i32 {
    oxwave_lib_wpc2::huffDecompress(input, output)
}

// Port of pcoding.c, author Tong Chen, 05/31/95.
/// C prototype: `memBUFF *pEncode(int *qx, int n);`
/// Preserves C's bounded 2*n allocation: prefix codes can be truncated when
/// more than half of the input values need three bytes.
pub fn pEncode(qx: &[i32]) -> WpcBuffer {
    let mut b = buffAlloc1(2 * qx.len());
    for &x in qx {
        if x > 126 || x < -126 {
            let s = x.unsigned_abs().min(0x7fff) * 2 + 1;
            for c in [if x > 0 { 127 } else { 129 }, (s >> 8) as u8, s as u8] {
                if b.pos < b.code.len() {
                    b.code[b.pos] = c;
                    b.pos += 1;
                }
            }
        } else if b.pos < b.code.len() {
            b.code[b.pos] = x as u8;
            b.pos += 1;
        }
    }
    let n = b.pos;
    buffRealloc1(&mut b, n);
    b
}
/// C prototype: `int pDecode(memBUFF *buff, int *qx, int n);`
/// The output must fit all decoded values, as in C. Its length also supplies
/// the expected count checked by C; a count mismatch returns -1.
pub fn pDecode(buff: &mut WpcBuffer, qx: &mut [i32]) -> i32 {
    let mut i = 0;
    let mut c1 = 0u8;
    let mut c2 = 0u8;
    while buff.pos < buff.code.len() {
        let sc = buff.code[buff.pos] as i8;
        buff.pos += 1;
        let x = if sc == 127 || sc == -127 {
            if buff.pos < buff.code.len() {
                c1 = buff.code[buff.pos];
                buff.pos += 1;
            }
            if buff.pos < buff.code.len() {
                c2 = buff.code[buff.pos];
                buff.pos += 1;
            }
            let s = (c1 as i32 * 256 + c2 as i32) >> 1;
            if sc == 127 {
                s
            } else {
                -s
            }
        } else {
            sc as i32
        };
        qx[i] = x;
        i += 1;
    }
    if i == qx.len() {
        0
    } else {
        -1
    }
}

// Port of lct1.c, author Tong Chen, 05/22/95.
/// Rising cutoff samples, indexed from `-rlast` to `rlast - 1`.
#[derive(Clone, Debug)]
pub struct RcFunct {
    pub rlast: usize,
    pub rc: Vec<f32>,
}
/// C prototype: `rcFunct *rcalloc(int hlen);`
pub fn rcalloc(hlen: usize) -> RcFunct {
    RcFunct {
        rlast: hlen,
        rc: vec![0.0; 2 * hlen],
    }
}
/// C prototype: `void rcfmidp(rcFunct *rcf);`
pub fn rcfmidp(rcf: &mut RcFunct) {
    let h = rcf.rlast;
    let mut x = (0.5 / h as f64) as f32;
    let dx = (1.0 / h as f64) as f32;
    fn sample(t: f32) -> f32 {
        if t <= -1.0 {
            0.0
        } else if t >= 1.0 {
            1.0
        } else {
            let tmp = (0.5 * std::f64::consts::PI * t as f64).sin() as f32;
            (0.25 * std::f64::consts::PI * (1.0 + tmp as f64)).sin() as f32
        }
    }
    for j in 0..h {
        rcf.rc[h + j] = sample(x);
        rcf.rc[h - j - 1] = sample(-x);
        x += dx;
    }
}

fn fold(sig: &mut [f32], at: usize, rcf: &RcFunct, inverse: bool) {
    let h = rcf.rlast;
    for k in 0..h {
        let a = sig[at + k];
        let b = sig[at - k - 1];
        let p = rcf.rc[h + k];
        let m = rcf.rc[h - k - 1];
        if inverse {
            sig[at + k] = p * a - m * b;
            sig[at - k - 1] = p * b + m * a;
        } else {
            sig[at + k] = p * a + m * b;
            sig[at - k - 1] = p * b - m * a;
        }
    }
}
/// C prototype: `void lct(float *sig, int n, rcFunct *rcf, float **c);`
/// `sig` includes the cutoff halo on each side; `start` selects the first sample.
pub fn lct(sig: &mut [f32], start: usize, n: usize, rcf: &RcFunct, c: &[Vec<f32>]) {
    fold(sig, start, rcf, false);
    fold(sig, start + n, rcf, false);
    dctiv(&mut sig[start..start + n], c);
}
/// C prototype: `void ilct(float *sig, int n, rcFunct *rcf, float **c);`
/// `sig` includes the cutoff halo on each side; `start` selects the first sample.
pub fn ilct(sig: &mut [f32], start: usize, n: usize, rcf: &RcFunct, c: &[Vec<f32>]) {
    dctiv(&mut sig[start..start + n], c);
    fold(sig, start, rcf, true);
    fold(sig, start + n, rcf, true);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_quantization_preserves_ieee_arithmetic() {
        let mut ave = 1.0;
        let mut step = -1.0;
        uniQuant(&[], 0.01, &mut ave, &mut step, &mut []);
        assert!(ave.is_nan());
        assert!(step.is_nan());

        // A supplied deviation is still used, even though the mean is NaN.
        step = 2.0;
        uniQuant(&[], 0.5, &mut ave, &mut step, &mut []);
        assert!(ave.is_nan());
        assert_eq!(step, (1.0 / 3.4641) as f32);
    }
}
