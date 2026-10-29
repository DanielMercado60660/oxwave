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

//! Port of comp/dwpt/2d/lib/{wpcbuffal,wpcsilence,wpchuff,wpcpack2}.c.
//! Credits: CWP: Tong Chen, 1994; Huffman algorithm based on Mark Nelson,
//! "Data Compression". Transform: Tong Chen, 05/25/94; API: 07/28/94.
#![allow(non_snake_case)]

/// Bounded byte buffer; `pos` is the current read or write position.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WpcBuffer {
    pub code: Vec<u8>,
    pub pos: usize,
}

impl WpcBuffer {
    fn get(&mut self) -> Option<u8> {
        if self.pos == self.code.len() {
            None
        } else {
            let c = self.code[self.pos];
            self.pos += 1;
            Some(c)
        }
    }
    fn put(&mut self, c: u8) {
        if self.pos != self.code.len() {
            self.code[self.pos] = c;
            self.pos += 1;
        }
    }
}

/// C prototype: `void *buffAlloc1(int n)`.
pub fn buffAlloc1(n: usize) -> WpcBuffer {
    WpcBuffer {
        code: vec![0; n],
        pos: 0,
    }
}

/// C prototype: `void buffRealloc1(void *inbuff, int n)`.
pub fn buffRealloc1(buff: &mut WpcBuffer, n: usize) {
    buff.code.resize(n, 0);
    if buff.pos > n {
        buff.pos = n;
    }
}

/// C prototype: `void buffFree1(void *inbuff)`.
pub fn buffFree1(_buff: WpcBuffer) {}

/// C prototype: `void **buffAlloc2(int n1, int n2)`.
pub fn buffAlloc2(n1: usize, n2: usize) -> Vec<WpcBuffer> {
    (0..n1).map(|_| buffAlloc1(n2)).collect()
}

/// C prototype: `void buffFree2(void **inbuff, int n1)`.
pub fn buffFree2(_buff: Vec<WpcBuffer>) {}

/// C prototype: `void codeSilence(void *inb, void *outb)`.
pub fn codeSilence(input: &mut WpcBuffer, output: &mut WpcBuffer) {
    let mut look = [-1i32; 8];
    let mut i = 0;
    while i < 8 {
        if let Some(c) = input.get() {
            look[i] = c as i32;
            i += 1;
        } else {
            break;
        }
    }
    if i < 8 {
        if output.pos + input.pos <= output.code.len() {
            output.code[output.pos..output.pos + input.pos]
                .copy_from_slice(&input.code[..input.pos]);
            output.pos += input.pos;
            input.pos = 0;
        }
        return;
    }
    let mut index = 0;
    loop {
        if look[index] == -1 {
            break;
        }
        if (0..5).all(|j| look[(index + j) & 7] == 0) {
            let mut run = 0;
            loop {
                look[index] = input.get().map_or(-1, |c| c as i32);
                index = (index + 1) & 7;
                run += 1;
                if run == 255 {
                    output.put(128);
                    output.put(255);
                    run = 0;
                }
                if look[index] != 0 {
                    break;
                }
            }
            if run > 0 {
                output.put(128);
                output.put(run);
            }
        }
        if look[index] == -1 {
            break;
        }
        if look[index] == 128 {
            look[index] += 1;
        }
        output.put(look[index] as u8);
        look[index] = input.get().map_or(-1, |c| c as i32);
        index = (index + 1) & 7;
    }
}

/// C prototype: `void codeDesilence(void *inb, void *outb)`.
pub fn codeDesilence(input: &mut WpcBuffer, output: &mut WpcBuffer) {
    let mut run_count = 0;
    while let Some(c) = input.get() {
        if c == 128 {
            if let Some(n) = input.get() {
                run_count = n as i32;
            }
            while run_count > 0 {
                output.put(0);
                run_count -= 1;
            }
            run_count -= 1;
        } else {
            output.put(c);
        }
    }
}

/// A node from the Huffman tree in `wpchuff.c` (Tong Chen).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HuffNode {
    pub count: u32,
    pub saved_count: u32,
    pub child_0: i32,
    pub child_1: i32,
}
/// A Huffman code from `wpchuff.c`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HuffCode {
    pub code: u32,
    pub code_bits: i32,
}

/// C prototype: `void huffCountbytes(wpcBUFF *buff, unsigned int *counts)`.
pub fn huffCountbytes(buff: &mut WpcBuffer, counts: &mut [u32]) {
    let pos = buff.pos;
    while let Some(c) = buff.get() {
        counts[c as usize] = counts[c as usize].wrapping_add(1);
    }
    buff.pos = pos;
}

/// C prototype: `void huffScalecounts(unsigned int *counts, huffNODE *nodes)`.
pub fn huffScalecounts(counts: &mut [u32], nodes: &mut [HuffNode]) {
    let mut max_count = 0;
    for &c in &counts[..256] {
        if c > max_count {
            max_count = c;
        }
    }
    if max_count == 0 {
        counts[0] = 1;
        max_count = 1;
    }
    max_count = max_count / 255 + 1;
    for i in 0..256 {
        nodes[i].count = counts[i] / max_count;
        if nodes[i].count == 0 && counts[i] != 0 {
            nodes[i].count = 1;
        }
    }
    nodes[256].count = 1;
}

/// C prototype: `int huffBuildtree(huffNODE *nodes)`.
pub fn huffBuildtree(nodes: &mut [HuffNode]) -> i32 {
    nodes[513].count = 0xffff;
    let mut next = 257;
    loop {
        let (mut a, mut b) = (513, 513);
        for i in 0..next {
            if nodes[i].count != 0 {
                if nodes[i].count < nodes[a].count {
                    b = a;
                    a = i;
                } else if nodes[i].count < nodes[b].count {
                    b = i;
                }
            }
        }
        if b == 513 {
            break;
        }
        nodes[next].count = nodes[a].count.wrapping_add(nodes[b].count);
        nodes[a].saved_count = nodes[a].count;
        nodes[a].count = 0;
        nodes[b].saved_count = nodes[b].count;
        nodes[b].count = 0;
        nodes[next].child_0 = a as i32;
        nodes[next].child_1 = b as i32;
        next += 1;
    }
    next -= 1;
    nodes[next].saved_count = nodes[next].count;
    next as i32
}

/// C prototype: `void huffConverttreetocode(huffNODE *nodes, huffCODE *codes, unsigned int code_so_far, int bits, int node)`.
pub fn huffConverttreetocode(
    nodes: &[HuffNode],
    codes: &mut [HuffCode],
    code: u32,
    bits: i32,
    node: i32,
) {
    if node <= 256 {
        codes[node as usize] = HuffCode {
            code,
            code_bits: bits,
        };
        return;
    }
    let code = code.wrapping_shl(1);
    huffConverttreetocode(nodes, codes, code, bits + 1, nodes[node as usize].child_0);
    huffConverttreetocode(
        nodes,
        codes,
        code | 1,
        bits + 1,
        nodes[node as usize].child_1,
    );
}

/// C prototype: `void waveReorder(float **x, float **y, int twopow1, int twopow2, int stage1, int stage2, int type)`.
pub fn waveReorder(
    x: &mut [Vec<f32>],
    y: &mut [Vec<f32>],
    twopow1: usize,
    twopow2: usize,
    stage1: usize,
    stage2: usize,
    kind: i32,
) {
    let nband1 = 1 << stage1;
    let nband2 = 1 << stage2;
    let lband1 = 1 << (twopow1 - stage1);
    let lband2 = 1 << (twopow2 - stage2);
    for j2 in 0..nband2 {
        for j1 in 0..nband1 {
            let block = j2 * nband1 + j1;
            let mut offset = 0;
            for i2 in j2 * lband2..(j2 + 1) * lband2 {
                for i1 in j1 * lband1..(j1 + 1) * lband1 {
                    if kind == 0 {
                        y[block][offset] = x[i2][i1];
                    } else {
                        x[i2][i1] = y[block][offset];
                    }
                    offset += 1;
                }
            }
        }
    }
}

mod codec;
mod coding;
mod huffman;
mod quant;
mod transform;
pub use codec::*;
pub use coding::*;
pub use huffman::*;
pub use quant::*;
pub use transform::*;

impl WpcBuffer {
    fn write(&mut self, data: &[u8]) -> bool {
        if self.pos + data.len() > self.code.len() {
            return false;
        }
        self.code[self.pos..self.pos + data.len()].copy_from_slice(data);
        self.pos += data.len();
        true
    }
    fn read<const N: usize>(&mut self) -> Option<[u8; N]> {
        if self.pos + N > self.code.len() {
            return None;
        }
        let mut bytes = [0; N];
        bytes.copy_from_slice(&self.code[self.pos..self.pos + N]);
        self.pos += N;
        Some(bytes)
    }
    fn put_checked(&mut self, c: u8) -> bool {
        if self.pos == self.code.len() {
            false
        } else {
            self.put(c);
            true
        }
    }
    fn merge(&mut self, other: &mut WpcBuffer) -> bool {
        if !self.write(&other.code[..other.pos]) {
            return false;
        }
        other.pos = 0;
        true
    }
}
