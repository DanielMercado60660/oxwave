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

//! Port of comp/dwpt/2d/lib/wpchuff.c. Author: Tong Chen, 08/01/94.
//! Algorithm based on Mark Nelson's "Data Compression" book.
use crate::*;

/// State of the C `wpcBITBUFF`; initialize mask to 0x80 and rack to zero.
pub struct WpcBitBuffer<'a> {
    pub buff: &'a mut WpcBuffer,
    pub mask: u8,
    pub rack: u8,
}
impl WpcBitBuffer<'_> {
    fn output(&mut self, code: u32, count: i32) -> bool {
        let mut lmask = 1u32.wrapping_shl((count - 1) as u32);
        while lmask != 0 {
            if lmask & code != 0 {
                self.rack |= self.mask;
            }
            self.mask >>= 1;
            if self.mask == 0 {
                if !self.buff.put_checked(self.rack) {
                    return false;
                }
                self.rack = 0;
                self.mask = 0x80;
            }
            lmask >>= 1;
        }
        true
    }
    fn input(&mut self) -> Option<bool> {
        if self.mask == 0x80 {
            self.rack = self.buff.get()?;
        }
        let value = self.rack & self.mask;
        self.mask >>= 1;
        if self.mask == 0 {
            self.mask = 0x80;
        }
        Some(value != 0)
    }
}

/// C prototype: `int huffOutputcounts(wpcBITBUFF *bit_buff, huffNODE *nodes)`.
pub fn huffOutputcounts(bit: &mut WpcBitBuffer<'_>, nodes: &[HuffNode]) -> i32 {
    let mut nbytes = 0;
    let mut first = 0;
    while first < 255 && nodes[first].count == 0 {
        first += 1;
    }
    while first < 256 {
        let mut last = first + 1;
        let next;
        loop {
            while last < 256 && nodes[last].count != 0 {
                last += 1;
            }
            last -= 1;
            let mut n = last + 1;
            while n < 256 && nodes[n].count == 0 {
                n += 1;
            }
            if n > 255 || n - last > 3 {
                next = n;
                break;
            }
            last = n;
        }
        if !bit.buff.put_checked(first as u8) {
            return -1;
        }
        nbytes += 1;
        if !bit.buff.put_checked(last as u8) {
            return -1;
        }
        nbytes += 1;
        for i in first..=last {
            if !bit.buff.put_checked(nodes[i].count as u8) {
                return -1;
            }
            nbytes += 1;
        }
        first = next;
    }
    if !bit.buff.put_checked(0) {
        return -1;
    }
    nbytes + 1
}

/// C prototype: `int huffInputcounts(wpcBITBUFF *bit_buff, huffNODE *nodes)`.
pub fn huffInputcounts(bit: &mut WpcBitBuffer<'_>, nodes: &mut [HuffNode]) -> i32 {
    for node in &mut nodes[..256] {
        node.count = 0;
    }
    let Some(mut first) = bit.buff.get() else {
        return -1;
    };
    let Some(mut last) = bit.buff.get() else {
        return -1;
    };
    loop {
        for i in first as usize..=last as usize {
            let Some(c) = bit.buff.get() else {
                return -1;
            };
            nodes[i].count = c as u32;
        }
        let Some(c) = bit.buff.get() else {
            return -1;
        };
        first = c;
        if first == 0 {
            break;
        }
        let Some(c) = bit.buff.get() else {
            return -1;
        };
        last = c;
    }
    nodes[256].count = 1;
    0
}

/// C prototype: `int huffEncoder(wpcBUFF *inbuff, wpcBITBUFF *bit_buff, huffCODE *codes)`.
pub fn huffEncoder(input: &mut WpcBuffer, bit: &mut WpcBitBuffer<'_>, codes: &[HuffCode]) -> i32 {
    let nbytes = bit.buff.pos;
    while let Some(c) = input.get() {
        if !bit.output(codes[c as usize].code, codes[c as usize].code_bits) {
            return -1;
        }
    }
    if !bit.output(codes[256].code, codes[256].code_bits) {
        return -1;
    }
    (bit.buff.pos - nbytes) as i32
}

/// C prototype: `int huffDecoder(wpcBITBUFF *bit_buff, wpcBUFF *outbuff, huffNODE *nodes, int root_node)`.
pub fn huffDecoder(
    bit: &mut WpcBitBuffer<'_>,
    output: &mut WpcBuffer,
    nodes: &[HuffNode],
    root: i32,
) -> i32 {
    loop {
        let mut node = root;
        loop {
            let Some(b) = bit.input() else {
                return -1;
            };
            node = if b {
                nodes[node as usize].child_1
            } else {
                nodes[node as usize].child_0
            };
            if node <= 256 {
                break;
            }
        }
        if node == 256 {
            break;
        }
        output.put(node as u8); // C deliberately ignores a full output buffer.
    }
    0
}

/// C prototype: `int huffCompress(void *inb, void *outb)`.
pub fn huffCompress(input: &mut WpcBuffer, output: &mut WpcBuffer) -> i32 {
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
    if retval == -1 {
        return -1;
    }
    nbytes += retval;
    if bit.mask != 0x80 {
        if !bit.buff.put_checked(bit.rack) {
            return -1;
        }
        nbytes += 1;
    }
    nbytes
}

/// C prototype: `int huffDecompress(void *inb, void *outb)`.
pub fn huffDecompress(input: &mut WpcBuffer, output: &mut WpcBuffer) -> i32 {
    let mut bit = WpcBitBuffer {
        buff: input,
        mask: 0x80,
        rack: 0,
    };
    let mut nodes = [HuffNode::default(); 514];
    if huffInputcounts(&mut bit, &mut nodes) == -1 {
        return -1;
    }
    let root = huffBuildtree(&mut nodes);
    huffDecoder(&mut bit, output, &nodes, root)
}
