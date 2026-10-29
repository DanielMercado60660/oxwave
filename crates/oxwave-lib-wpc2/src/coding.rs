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

//! Port of comp/dwpt/2d/lib/wpccoding.c. Author: Tong Chen, 07/29/94.
use crate::*;

/// C prototype: `int codeEncoder(void *inconf, void *qstate, void *obuff, void *interblock)`.
/// The C interblock argument is reserved and unused.
pub fn codeEncoder(c: &WpcConfig, q: &WpcQuant, output: &mut WpcBuffer) -> i32 {
    if !output.put_checked(if q.zeroflag == 1 { 1 } else { 0 }) {
        return -1;
    }
    if !output.write(&q.ave.to_ne_bytes()) {
        return -1;
    }
    if q.zeroflag == 1 {
        return 5;
    }
    if !output.write(&q.step.to_ne_bytes()) {
        return -1;
    }
    let mut outputsize = 9;
    let (nblock, lblock) = (c.nblock as usize, c.lblock as usize);
    let (codenblock, codelblock) = (c.codenblock as usize, c.codelblock as usize);
    let huffmbound = 4 * codelblock * lblock;
    let mut quantbuff = buffAlloc2(nblock, 4 * lblock);
    let mut rlebuff = buffAlloc2(codenblock, huffmbound);
    let mut huffbuff = buffAlloc2(codenblock, huffmbound);
    for b in 0..nblock {
        let p = &mut quantbuff[b];
        if q.flowflag[b] == 0 {
            p.put(0);
            for i in 0..lblock {
                p.put(q.qx[b][i] as u8);
            }
        } else {
            p.put(1);
            for i in 0..lblock {
                let v = q.qx[b][i];
                if v > 126 || v < -126 {
                    let mut s = if v > 0 {
                        v as u32
                    } else {
                        v.wrapping_neg() as u32
                    };
                    if s > 0x7fff {
                        s = 0x7fff;
                    }
                    s = (s << 1) + 1;
                    p.put(if v > 0 { 127 } else { (-127i8) as u8 });
                    for _ in 0..2 {
                        p.put((s & 0xff) as u8);
                        s >>= 8;
                    }
                } else {
                    p.put(v as u8);
                }
            }
        }
        buffRealloc1(p, p.pos);
    }
    for b in 0..nblock {
        huffbuff[0].put(q.blockind[b]);
    }
    for b in 0..nblock {
        huffbuff[0].write(&(quantbuff[b].code.len() as i32).to_ne_bytes());
    }
    for b in 0..codenblock {
        let p = &mut rlebuff[b];
        for i in 0..codelblock {
            let j = (codenblock - 1 - b) * codelblock + i;
            let bind = q.blockind[j] as usize;
            p.merge(&mut quantbuff[bind]);
        }
        buffRealloc1(p, p.pos);
        p.pos = 0;
    }
    for b in 0..codenblock {
        codeSilence(&mut rlebuff[b], &mut huffbuff[b]);
        let p = &mut huffbuff[b];
        buffRealloc1(p, p.pos);
    }
    for b in 0..codenblock {
        let p = &mut huffbuff[b];
        p.pos = 0;
        let nsize = huffCompress(p, output);
        if nsize == -1 {
            return -1;
        }
        outputsize += nsize;
    }
    outputsize
}

/// C prototype: `int codeDecoder(void *inconf, void *qstate, void *ibuff, void *interblock)`.
/// Stores the decoded zero flag and average, repairing C's missing assignments
/// so zero tiles do not dequantize using unset state and produce NaNs.
pub fn codeDecoder(c: &WpcConfig, q: &mut WpcQuant, input: &mut WpcBuffer) -> i32 {
    let Some(zeroflag) = input.get() else {
        return 0;
    };
    if zeroflag == 1 {
        let Some(ave) = input.read::<4>() else {
            return 0;
        };
        q.zeroflag = zeroflag;
        q.ave = f32::from_ne_bytes(ave);
        return 1;
    }
    if zeroflag != 0 {
        return 0;
    }
    let Some(ave) = input.read::<4>() else {
        return 0;
    };
    let Some(step) = input.read::<4>() else {
        return 0;
    };
    q.zeroflag = zeroflag;
    q.ave = f32::from_ne_bytes(ave);
    q.step = f32::from_ne_bytes(step);
    let (nblock, lblock) = (c.nblock as usize, c.lblock as usize);
    let (codenblock, codelblock) = (c.codenblock as usize, c.codelblock as usize);
    let huffmbound = 4 * codelblock * lblock;
    let mut rlebuff = buffAlloc2(codenblock, huffmbound);
    let mut huffbuff = buffAlloc2(codenblock, huffmbound);
    for b in 0..codenblock {
        let p = &mut huffbuff[b];
        if huffDecompress(input, p) == -1 {
            return 0;
        }
        buffRealloc1(p, p.pos);
        p.pos = 0;
    }
    for b in 0..nblock {
        let Some(v) = huffbuff[0].get() else {
            return 0;
        };
        q.blockind[b] = v;
        if v as usize > nblock {
            return 0;
        }
    }
    let mut bufflen = vec![0; nblock];
    for b in 0..nblock {
        let Some(v) = huffbuff[0].read::<4>() else {
            return 0;
        };
        let len = i32::from_ne_bytes(v);
        if len < 0 {
            return 0;
        }
        bufflen[b] = len as usize;
    }
    for b in 0..codenblock {
        codeDesilence(&mut huffbuff[b], &mut rlebuff[b]);
        let p = &mut rlebuff[b];
        buffRealloc1(p, p.pos);
    }
    let mut quantbuff: Vec<_> = bufflen.iter().map(|&n| buffAlloc1(n)).collect();
    for b in 0..codenblock {
        let p = &mut rlebuff[b];
        p.pos = 0;
        for i in 0..codelblock {
            let j = (codenblock - 1 - b) * codelblock + i;
            let bind = q.blockind[j] as usize;
            // C has no checks for corrupted indices/lengths before memcpy.
            // Reject those streams instead of accessing memory outside a slice.
            if bind >= nblock {
                return 0;
            }
            let len = bufflen[bind];
            if p.pos + len <= p.code.len() {
                quantbuff[bind]
                    .code
                    .copy_from_slice(&p.code[p.pos..p.pos + len]);
                quantbuff[bind].pos = 0;
                p.pos += len;
            }
        }
    }
    let mut last_c = 0u8;
    for b in 0..nblock {
        let p = &mut quantbuff[b];
        p.pos = 0;
        if let Some(v) = p.get() {
            q.flowflag[b] = v;
        }
        let mut i = 0;
        while let Some(v) = p.get() {
            let sc = v as i8;
            let value = if q.flowflag[b] != 0 && (sc == 127 || sc == -127) {
                if let Some(v) = p.get() {
                    last_c = v;
                }
                let mut s = last_c as u32;
                if let Some(v) = p.get() {
                    last_c = v;
                }
                s += (last_c as u32) << 8;
                if sc == 127 {
                    (s >> 1) as i32
                } else {
                    -((s >> 1) as i32)
                }
            } else {
                sc as i32
            };
            if i >= lblock {
                return 0;
            }
            q.qx[b][i] = value;
            i += 1;
        }
        if i != lblock {
            return 0;
        }
    }
    1
}
