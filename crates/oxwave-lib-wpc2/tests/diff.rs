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

//! Differential tests for comp/dwpt/2d/lib/wpc*.c. Author: CWP: Tong Chen, 1994.
//! Huffman algorithm based on Mark Nelson's "Data Compression" book.
#![allow(non_snake_case)]
use oxwave_cref::XorShift;
use oxwave_lib_wpc2 as rust;
use oxwave_lib_wpc2::*;
use std::ffi::c_void;
use std::sync::Mutex;
// Both C and oxwave-sort use persistent qkisort seeds. Pair the calls in
// the same order so that ties use the same pivot history.
static SORT_PAIR: Mutex<()> = Mutex::new(());

#[repr(C)]
struct CBitBuffer {
    buff: *mut CBuffer,
    mask: u8,
    rack: u8,
}
#[repr(C)]
struct CWaveFilter {
    len: i32,
    ishift: i32,
    filterh: *mut f32,
    filterg: *mut f32,
}
#[repr(C)]
struct CQuant {
    zeroflag: u8,
    flowflag: *mut u8,
    blockind: *mut u8,
    ave: f32,
    step: f32,
    qx: *mut *mut i32,
}
#[repr(C)]
struct CPacked {
    version: i32,
    n1: i32,
    n2: i32,
    size: i32,
    config: *mut WpcConfig,
    code: *mut u8,
}

#[repr(C)]
struct CBuffer {
    bound: i32,
    pos: i32,
    code: *mut u8,
}
#[link(name = "wpc", kind = "static")]
extern "C" {
    fn codeEncoder(c: *const WpcConfig, q: *mut CQuant, b: *mut CBuffer, inter: *mut c_void)
        -> i32;
    fn codeDecoder(c: *const WpcConfig, q: *mut CQuant, b: *mut CBuffer, inter: *mut c_void)
        -> i32;
    fn huffCompress(a: *mut CBuffer, b: *mut CBuffer) -> i32;
    fn huffDecompress(a: *mut CBuffer, b: *mut CBuffer) -> i32;
    fn huffEncoder(a: *mut CBuffer, b: *mut CBitBuffer, codes: *mut HuffCode) -> i32;
    fn huffDecoder(a: *mut CBitBuffer, b: *mut CBuffer, nodes: *mut HuffNode, root: i32) -> i32;
    fn huffInputcounts(b: *mut CBitBuffer, nodes: *mut HuffNode) -> i32;
    fn huffOutputcounts(b: *mut CBitBuffer, nodes: *mut HuffNode) -> i32;
    fn quantFixerror(x: *mut *mut f32, error: f32, c: *const WpcConfig, q: *mut CQuant);
    fn quantDe(x: *mut *mut f32, c: *const WpcConfig, q: *mut CQuant);
    fn wavePack1(x: *mut f32, y: *mut f32, f: *mut CWaveFilter, p: i32, s: i32, kind: i32);
    fn wavePack2(x: *mut *mut f32, y: *mut *mut f32, c: *const WpcConfig, kind: i32);
    fn wpcEncoder(
        x: *mut *mut f32,
        error: f32,
        c: *const WpcConfig,
        b: *mut CBuffer,
        inter: *mut c_void,
    ) -> i32;
    fn wpcDecoder(
        x: *mut *mut f32,
        c: *const WpcConfig,
        b: *mut CBuffer,
        inter: *mut c_void,
    ) -> i32;
    fn wpcCompress(f: *mut f32, n1: i32, n2: i32, error: f32, size: *mut i32) -> *mut CPacked;
    fn wpcUncompress(w: *mut CPacked, f: *mut f32) -> i32;
    fn wpcN1(w: *mut CPacked) -> i32;
    fn wpcN2(w: *mut CPacked) -> i32;
    fn wpcSize(w: *mut CPacked) -> i32;
    fn wpcFree(w: *mut CPacked);
    fn wpcRead(fp: *mut c_void) -> *mut CPacked;
    fn wpcWrite(w: *mut CPacked, fp: *mut c_void);
    fn fopen(path: *const std::ffi::c_char, mode: *const std::ffi::c_char) -> *mut c_void;
    fn rewind(fp: *mut c_void);
    fn fclose(fp: *mut c_void) -> i32;
    fn fread(p: *mut c_void, size: usize, n: usize, fp: *mut c_void) -> usize;
    fn fwrite(p: *const c_void, size: usize, n: usize, fp: *mut c_void) -> usize;
    fn free(p: *mut c_void);

    fn huffCountbytes(buff: *mut CBuffer, counts: *mut u32);
    fn huffScalecounts(counts: *mut u32, nodes: *mut HuffNode);
    fn huffBuildtree(nodes: *mut HuffNode) -> i32;
    fn huffConverttreetocode(
        nodes: *mut HuffNode,
        codes: *mut HuffCode,
        code: u32,
        bits: i32,
        node: i32,
    );
    fn waveReorder(
        x: *mut *mut f32,
        y: *mut *mut f32,
        p1: i32,
        p2: i32,
        s1: i32,
        s2: i32,
        kind: i32,
    );
    fn buffAlloc1(n: i32) -> *mut CBuffer;
    fn buffAlloc2(n: i32, m: i32) -> *mut *mut CBuffer;
    fn buffRealloc1(b: *mut CBuffer, n: i32);
    fn buffFree1(b: *mut CBuffer);
    fn buffFree2(b: *mut *mut CBuffer, n: i32);
    fn codeSilence(a: *mut c_void, b: *mut c_void);
    fn codeDesilence(a: *mut c_void, b: *mut c_void);
}

#[test]
fn buffAlloc1_matches() {
    for n in 0..3000 {
        unsafe {
            let c = buffAlloc1(n);
            let r = oxwave_lib_wpc2::buffAlloc1(n as usize);
            assert_eq!(((*c).bound, (*c).pos), (r.code.len() as i32, r.pos as i32));
            buffFree1(c);
        }
    }
}
#[test]
fn buffRealloc1_matches() {
    let mut rng = XorShift::new(79);
    for case in 0..3000 {
        unsafe {
            let old = [1, 2, 3, 7, 31, 127, 256, 257, 1024, 8191][case % 10];
            let n = [0, 1, 2, 7, 32, 127, 256, 257, 1024, 8192][(case / 10) % 10];
            let c = buffAlloc1(old as i32);
            let mut r = rust::buffAlloc1(old);
            for i in 0..old {
                let v = rng.next_u32() as u8;
                *(*c).code.add(i) = v;
                r.code[i] = v;
            }
            let pos = rng.next_u32() as usize % (old + 1);
            (*c).pos = pos as i32;
            r.pos = pos;
            buffRealloc1(c, n as i32);
            rust::buffRealloc1(&mut r, n);
            assert_eq!(((*c).bound, (*c).pos), (r.code.len() as i32, r.pos as i32));
            if old.min(n) > 0 {
                assert_eq!(
                    std::slice::from_raw_parts((*c).code, old.min(n)),
                    &r.code[..old.min(n)]
                );
            }
            buffFree1(c);
        }
    }
}
#[test]
fn buffFree1_matches() {
    for n in 0..3000 {
        unsafe {
            buffFree1(buffAlloc1(n));
        }
        oxwave_lib_wpc2::buffFree1(oxwave_lib_wpc2::buffAlloc1(n as usize));
    }
}
#[test]
fn buffAlloc2_matches() {
    for n in 0..3000 {
        unsafe {
            let count = [0, 1, 2, 3, 7, 16][n % 6];
            let c = buffAlloc2(count, n as i32);
            let r = rust::buffAlloc2(count as usize, n);
            for i in 0..count as usize {
                assert_eq!(
                    ((*(*c.add(i))).bound, (*(*c.add(i))).pos),
                    (r[i].code.len() as i32, r[i].pos as i32)
                );
            }
            buffFree2(c, count);
        }
    }
}
#[test]
fn buffFree2_matches() {
    for n in 0..3000 {
        let count = [0, 1, 2, 3, 7, 16][n % 6];
        unsafe {
            buffFree2(buffAlloc2(count, n as i32), count);
        }
        rust::buffFree2(rust::buffAlloc2(count as usize, n));
    }
}
fn compare(input: Vec<u8>, encode: bool) {
    unsafe {
        let n = input.len();
        let mut source = input.clone();
        let mut target = vec![0u8; n * if encode { 2 } else { 256 } + 512];
        let mut a = CBuffer {
            bound: n as i32,
            pos: 0,
            code: source.as_mut_ptr(),
        };
        let mut b = CBuffer {
            bound: target.len() as i32,
            pos: 0,
            code: target.as_mut_ptr(),
        };
        let mut ar = WpcBuffer {
            code: input,
            pos: 0,
        };
        let mut br = oxwave_lib_wpc2::buffAlloc1(target.len());
        if encode {
            codeSilence(&mut a as *mut _ as _, &mut b as *mut _ as _);
            oxwave_lib_wpc2::codeSilence(&mut ar, &mut br);
        } else {
            codeDesilence(&mut a as *mut _ as _, &mut b as *mut _ as _);
            oxwave_lib_wpc2::codeDesilence(&mut ar, &mut br);
        }
        assert_eq!(a.pos as usize, ar.pos);
        assert_eq!(b.pos as usize, br.pos);
        assert_eq!(&target[..b.pos as usize], &br.code[..br.pos]);
    }
}
#[test]
fn codeSilence_matches() {
    let mut rng = XorShift::new(73);
    for i in 0..3000 {
        let mut v = vec![0; [0, 1, 2, 5, 7, 8, 31, 255, 256, 257, 1024, 8191][i % 12]];
        for x in &mut v {
            let r = rng.next_u32();
            *x = if r % 5 == 0 { r as u8 } else { 0 };
        }
        if i % 3 == 0 {
            v.fill(0);
        }
        compare(v, true);
    }
}
#[test]
fn codeDesilence_matches() {
    let mut rng = XorShift::new(91);
    for i in 0..3000 {
        let v = (0..i % 300).map(|_| rng.next_u32() as u8).collect();
        compare(v, false);
    }
}

#[test]
fn huffCountbytes_matches() {
    let mut rng = XorShift::new(123);
    for n in 0..3000 {
        let mut data: Vec<u8> = (0..n % 500).map(|_| rng.next_u32() as u8).collect();
        let pos = n % (data.len() + 1);
        let mut c = CBuffer {
            bound: data.len() as i32,
            pos: pos as i32,
            code: data.as_mut_ptr(),
        };
        let mut r = WpcBuffer {
            code: data.clone(),
            pos,
        };
        let mut a = [0u32; 256];
        let mut b = a;
        unsafe {
            huffCountbytes(&mut c, a.as_mut_ptr());
        }
        oxwave_lib_wpc2::huffCountbytes(&mut r, &mut b);
        assert_eq!(a, b);
        assert_eq!(c.pos as usize, r.pos);
    }
}
fn huff_case(seed: u32) {
    let mut rng = XorShift::new(seed as u64 + 1);
    let mut counts = [0u32; 256];
    for x in &mut counts {
        let v = rng.next_u32();
        *x = if v % 4 == 0 { 0 } else { v % 100000 };
    }
    if seed == 0 {
        counts.fill(0);
    }
    let mut rc = counts;
    let mut cn = [HuffNode::default(); 514];
    let mut rn = cn;
    unsafe {
        huffScalecounts(counts.as_mut_ptr(), cn.as_mut_ptr());
    }
    oxwave_lib_wpc2::huffScalecounts(&mut rc, &mut rn);
    assert_eq!(counts, rc);
    assert_eq!(cn, rn);
    let root = unsafe { huffBuildtree(cn.as_mut_ptr()) };
    assert_eq!(root, oxwave_lib_wpc2::huffBuildtree(&mut rn));
    assert_eq!(cn, rn);
    let mut cc = [HuffCode::default(); 257];
    let mut rr = cc;
    unsafe {
        huffConverttreetocode(cn.as_mut_ptr(), cc.as_mut_ptr(), 0, 0, root);
    }
    oxwave_lib_wpc2::huffConverttreetocode(&rn, &mut rr, 0, 0, root);
    assert_eq!(cc, rr);
}
#[test]
fn huffScalecounts_matches() {
    for i in 0..3000 {
        huff_case(i);
    }
}
#[test]
fn huffBuildtree_matches() {
    for i in 3000..6000 {
        huff_case(i);
    }
}
#[test]
fn huffConverttreetocode_matches() {
    for i in 6000..9000 {
        huff_case(i);
    }
}

#[test]
fn waveReorder_matches() {
    let mut rng = XorShift::new(431);
    for case in 0..3000 {
        let p1 = 1 + case % 5;
        let p2 = 1 + (case / 5) % 5;
        let s1 = (case / 25) % (p1 + 1);
        let s2 = (case / 37) % (p2 + 1);
        let (n1, n2) = (1 << p1, 1 << p2);
        let (blocks, length) = (1 << (s1 + s2), 1 << (p1 + p2 - s1 - s2));
        let mut x: Vec<Vec<f32>> = (0..n2)
            .map(|_| (0..n1).map(|_| rng.next_u32() as i32 as f32).collect())
            .collect();
        let mut y: Vec<Vec<f32>> = (0..blocks)
            .map(|_| (0..length).map(|_| rng.next_u32() as i32 as f32).collect())
            .collect();
        let mut rx = x.clone();
        let mut ry = y.clone();
        let mut xp: Vec<_> = x.iter_mut().map(|r| r.as_mut_ptr()).collect();
        let mut yp: Vec<_> = y.iter_mut().map(|r| r.as_mut_ptr()).collect();
        let kind = (case % 2) as i32;
        unsafe {
            waveReorder(
                xp.as_mut_ptr(),
                yp.as_mut_ptr(),
                p1 as i32,
                p2 as i32,
                s1 as i32,
                s2 as i32,
                kind,
            );
        }
        oxwave_lib_wpc2::waveReorder(&mut rx, &mut ry, p1, p2, s1, s2, kind);
        assert_eq!(x, rx);
        assert_eq!(y, ry);
    }
}

fn cbuff(data: &mut [u8], pos: usize) -> CBuffer {
    CBuffer {
        bound: data.len() as i32,
        pos: pos as i32,
        code: data.as_mut_ptr(),
    }
}
fn byte_data(case: usize, rng: &mut XorShift) -> Vec<u8> {
    let sizes = [0, 1, 2, 3, 5, 7, 31, 127, 256, 257, 1024, 4096, 8191];
    (0..sizes[case % sizes.len()])
        .map(|_| {
            let v = rng.next_u32();
            match case % 5 {
                0 => 0,
                1 => 255,
                2 => (v % 3) as u8,
                3 => {
                    if v % 7 == 0 {
                        v as u8
                    } else {
                        0
                    }
                }
                _ => v as u8,
            }
        })
        .collect()
}
fn huff_stream_case(case: usize, rng: &mut XorShift) {
    let data = byte_data(case, rng);
    let mut source = data.clone();
    let mut a = cbuff(&mut source, 0);
    let mut out = vec![0xa5; data.len() * 2 + 1024];
    // Include full buffers and nonzero starting positions.
    if case % 7 == 0 {
        out.truncate(case % 263);
    }
    let pos = if out.is_empty() {
        0
    } else {
        case % out.len().min(7)
    };
    let mut b = cbuff(&mut out, pos);
    let mut ra = WpcBuffer {
        code: data.clone(),
        pos: 0,
    };
    let mut rb = WpcBuffer {
        code: out.clone(),
        pos,
    };
    let cc = unsafe { huffCompress(&mut a, &mut b) };
    assert_eq!(cc, rust::huffCompress(&mut ra, &mut rb));
    assert_eq!(a.pos as usize, ra.pos);
    assert_eq!(b.pos as usize, rb.pos);
    assert_eq!(out, rb.code);
    if cc != -1 {
        let mut packed = out[pos..b.pos as usize].to_vec();
        if case % 11 == 0 {
            packed.truncate(packed.len().saturating_sub(1));
        }
        let mut cinput = cbuff(&mut packed, 0);
        let mut rinput = WpcBuffer {
            code: packed.clone(),
            pos: 0,
        };
        let bound = if case % 3 == 0 {
            data.len() / 2
        } else {
            data.len() + 4
        };
        let mut decoded = vec![0xa5; bound];
        let mut rdecoded = WpcBuffer {
            code: decoded.clone(),
            pos: 0,
        };
        let mut coutput = cbuff(&mut decoded, 0);
        let cstatus = unsafe { huffDecompress(&mut cinput, &mut coutput) };
        assert_eq!(cstatus, rust::huffDecompress(&mut rinput, &mut rdecoded));
        assert_eq!(cinput.pos as usize, rinput.pos);
        assert_eq!(coutput.pos as usize, rdecoded.pos);
        assert_eq!(decoded, rdecoded.code);
    }
}
#[test]
fn huffCompress_matches() {
    let mut rng = XorShift::new(801);
    for i in 0..2100 {
        huff_stream_case(i, &mut rng);
    }
}
#[test]
fn huffDecompress_matches() {
    let mut rng = XorShift::new(802);
    for i in 2100..4200 {
        huff_stream_case(i, &mut rng);
    }
}

fn counts_case(case: usize, rng: &mut XorShift) {
    let mut nodes = [HuffNode::default(); 514];
    for node in &mut nodes[..256] {
        let r = rng.next_u32();
        node.count = if r % 5 == 0 { r % 255 + 1 } else { 0 };
    }
    if case % 13 == 0 {
        for n in &mut nodes[..256] {
            n.count = 0;
        }
        nodes[case % 256].count = 1;
    }
    let mut out = vec![0xa5; 1024];
    if case % 3 == 0 {
        out.truncate(case % 263);
    }
    let mut rb = WpcBuffer {
        code: out.clone(),
        pos: 0,
    };
    let mut cb = cbuff(&mut out, 0);
    let mask = 0x80 >> (case % 8);
    let mut bit = CBitBuffer {
        buff: &mut cb,
        mask,
        rack: 0x5a,
    };
    let mut rbit = WpcBitBuffer {
        buff: &mut rb,
        mask,
        rack: 0x5a,
    };
    let retval = unsafe { huffOutputcounts(&mut bit, nodes.as_mut_ptr()) };
    assert_eq!(retval, rust::huffOutputcounts(&mut rbit, &nodes));
    assert_eq!((bit.mask, bit.rack), (rbit.mask, rbit.rack));
    assert_eq!(cb.pos as usize, rb.pos);
    assert_eq!(out, rb.code);
    out.truncate(cb.pos as usize);
    let mut inp = cbuff(&mut out, 0);
    let mut rinp = WpcBuffer {
        code: out.clone(),
        pos: 0,
    };
    let mut bit = CBitBuffer {
        buff: &mut inp,
        mask,
        rack: 0x7b,
    };
    let mut rbit = WpcBitBuffer {
        buff: &mut rinp,
        mask,
        rack: 0x7b,
    };
    let mut cn = [HuffNode {
        count: 33,
        saved_count: 9,
        child_0: 3,
        child_1: 7,
    }; 514];
    let mut rn = cn;
    let ret = unsafe { huffInputcounts(&mut bit, cn.as_mut_ptr()) };
    assert_eq!(ret, rust::huffInputcounts(&mut rbit, &mut rn));
    assert_eq!(cn, rn);
    assert_eq!(inp.pos as usize, rbit.buff.pos);
    assert_eq!((bit.mask, bit.rack), (rbit.mask, rbit.rack));
}
#[test]
fn huffOutputcounts_matches() {
    let mut rng = XorShift::new(811);
    for i in 0..2100 {
        counts_case(i, &mut rng);
    }
}
#[test]
fn huffInputcounts_matches() {
    let mut rng = XorShift::new(812);
    for i in 2100..4200 {
        counts_case(i, &mut rng);
    }
}

fn huff_bits_case(case: usize, rng: &mut XorShift) {
    let data = byte_data(case, rng);
    let mut counts = [0u32; 256];
    for &v in &data {
        counts[v as usize] += 1;
    }
    let mut nodes = [HuffNode::default(); 514];
    rust::huffScalecounts(&mut counts, &mut nodes);
    let root = rust::huffBuildtree(&mut nodes);
    let mut codes = [HuffCode::default(); 257];
    rust::huffConverttreetocode(&nodes, &mut codes, 0, 0, root);
    let mut source = data.clone();
    let mut a = cbuff(&mut source, 0);
    let mut ra = WpcBuffer { code: data, pos: 0 };
    let mut bytes = vec![0xa5; source.len() * 2 + 32];
    let mut b = cbuff(&mut bytes, 0);
    let mut rb = WpcBuffer {
        code: bytes.clone(),
        pos: 0,
    };
    let mut bit = CBitBuffer {
        buff: &mut b,
        mask: 0x80,
        rack: 0,
    };
    let mut rbit = WpcBitBuffer {
        buff: &mut rb,
        mask: 0x80,
        rack: 0,
    };
    let ret = unsafe { huffEncoder(&mut a, &mut bit, codes.as_mut_ptr()) };
    assert_eq!(ret, rust::huffEncoder(&mut ra, &mut rbit, &codes));
    assert_eq!((bit.mask, bit.rack), (rbit.mask, rbit.rack));
    assert_eq!(a.pos as usize, ra.pos);
    assert_eq!(b.pos as usize, rbit.buff.pos);
    assert_eq!(bytes, rbit.buff.code);
    let mut len = b.pos as usize;
    if bit.mask != 0x80 {
        bytes[len] = bit.rack;
        len += 1;
    }
    bytes.truncate(len);
    if case % 9 == 0 {
        bytes.truncate(len.saturating_sub(1));
    }
    let mut input = cbuff(&mut bytes, 0);
    let mut ri = WpcBuffer {
        code: bytes.clone(),
        pos: 0,
    };
    let mut out = vec![0xa5; source.len() + 3];
    let mut output = cbuff(&mut out, 0);
    let mut ro = WpcBuffer {
        code: out.clone(),
        pos: 0,
    };
    let mut bit = CBitBuffer {
        buff: &mut input,
        mask: 0x80,
        rack: 0,
    };
    let mut rbit = WpcBitBuffer {
        buff: &mut ri,
        mask: 0x80,
        rack: 0,
    };
    let ret = unsafe { huffDecoder(&mut bit, &mut output, nodes.as_mut_ptr(), root) };
    assert_eq!(ret, rust::huffDecoder(&mut rbit, &mut ro, &nodes, root));
    assert_eq!((bit.mask, bit.rack), (rbit.mask, rbit.rack));
    assert_eq!(input.pos as usize, rbit.buff.pos);
    assert_eq!(output.pos as usize, ro.pos);
    assert_eq!(out, ro.code);
}
#[test]
fn huffEncoder_matches() {
    let mut rng = XorShift::new(821);
    for i in 0..2100 {
        huff_bits_case(i, &mut rng);
    }
}
#[test]
fn huffDecoder_matches() {
    let mut rng = XorShift::new(822);
    for i in 2100..4200 {
        huff_bits_case(i, &mut rng);
    }
}

fn config(case: usize) -> WpcConfig {
    let mut c = WpcConfig::default();
    c.powt = 1 + (case % 5) as i32;
    c.powx = 1 + ((case / 5) % 4) as i32;
    c.tileszt = 1 << c.powt;
    c.tileszx = 1 << c.powx;
    c.staget = 1 + ((case / 20) % c.powt as usize) as i32;
    c.stagex = 1 + ((case / 29) % c.powx as usize) as i32;
    c.nblock = 1 << (c.staget + c.stagex);
    c.lblock = c.tileszt * c.tileszx / c.nblock;
    // Coding requires room for 5*nblock bytes of header. Keep at least
    // two coefficients per block, and merge all blocks into one code block.
    if c.lblock < 2 {
        c.staget = 0;
        c.nblock = 1 << c.stagex;
        c.lblock = c.tileszt * c.tileszx / c.nblock;
    }
    c.codenblock = 1;
    c.codelblock = c.nblock;
    c
}
fn ptrs<T>(x: &mut [Vec<T>]) -> Vec<*mut T> {
    x.iter_mut().map(|v| v.as_mut_ptr()).collect()
}
fn floats(rng: &mut XorShift, n1: usize, n2: usize, case: usize) -> Vec<Vec<f32>> {
    (0..n2)
        .map(|_| {
            (0..n1)
                .map(|_| match case % 17 {
                    0 => 0.0,
                    1 => -0.0,
                    2 => -3.25,
                    3 => rng.uniform_f32(-1e-18, 1e-18),
                    4 => rng.uniform_f32(-1e12, 1e12),
                    _ => rng.uniform_f32(-100.0, 100.0),
                })
                .collect()
        })
        .collect()
}
fn same_floats(a: &[Vec<f32>], b: &[Vec<f32>]) {
    assert_eq!(a.len(), b.len());
    let mut max = 0.0f32;
    for row in a {
        for &v in row {
            if v.is_finite() && v.abs() > max {
                max = v.abs();
            }
        }
    }
    for (row, other) in a.iter().zip(b) {
        assert_eq!(row.len(), other.len());
        for (&x, &y) in row.iter().zip(other) {
            if x.is_nan() || y.is_nan() {
                assert!(x.is_nan() && y.is_nan(), "{x} != {y}");
            } else if x.is_infinite() || y.is_infinite() {
                assert_eq!(x.to_bits(), y.to_bits());
            } else {
                assert!(
                    x.to_bits() == y.to_bits() || (x - y).abs() <= max * 1e-5,
                    "{x} != {y}, max={max}"
                );
            }
        }
    }
}

fn quant_state(c: &WpcConfig) -> WpcQuant {
    WpcQuant {
        zeroflag: 9,
        flowflag: vec![7; c.nblock as usize],
        blockind: vec![11; c.nblock as usize],
        ave: 1.25,
        step: 7.0,
        qx: vec![vec![-123; c.lblock as usize]; c.nblock as usize],
    }
}
fn cq(q: &mut WpcQuant, p: &mut Vec<*mut i32>) -> CQuant {
    *p = ptrs(&mut q.qx);
    CQuant {
        zeroflag: q.zeroflag,
        flowflag: q.flowflag.as_mut_ptr(),
        blockind: q.blockind.as_mut_ptr(),
        ave: q.ave,
        step: q.step,
        qx: p.as_mut_ptr(),
    }
}
fn update_q(q: &mut WpcQuant, c: &CQuant) {
    q.zeroflag = c.zeroflag;
    q.ave = c.ave;
    q.step = c.step;
}
fn same_quant(a: &WpcQuant, b: &WpcQuant) {
    assert_eq!(a.zeroflag, b.zeroflag);
    assert_eq!(a.ave.to_bits(), b.ave.to_bits());
    assert_eq!(a.step.to_bits(), b.step.to_bits());
    assert_eq!(a.flowflag, b.flowflag);
    assert_eq!(a.blockind, b.blockind);
    assert_eq!(a.qx, b.qx);
}
#[test]
fn quantFixerror_matches() {
    let _lock = SORT_PAIR.lock().unwrap();
    let mut rng = XorShift::new(831);
    for i in 0..3000 {
        let mut c = config(i); // Include non-power-of-two coefficient block lengths.
        c.nblock = [1, 2, 3, 5, 7, 16, 64][i % 7];
        c.lblock = [1, 2, 7, 31, 64][(i / 7) % 5];
        c.tileszt = c.nblock * c.lblock;
        c.tileszx = 1;
        let mut x = floats(&mut rng, c.lblock as usize, c.nblock as usize, i);
        let mut rx = x.clone();
        let mut q = quant_state(&c);
        let mut rq = q.clone();
        let mut qp = vec![];
        let mut state = cq(&mut q, &mut qp);
        let mut xp = ptrs(&mut x);
        let error = [0.0, -0.01, 1e-6, 0.001, 0.01, 0.1, 1.0][i % 7];
        unsafe {
            quantFixerror(xp.as_mut_ptr(), error, &c, &mut state);
        }
        update_q(&mut q, &state);
        rust::quantFixerror(&mut rx, error, &c, &mut rq);
        same_quant(&q, &rq);
        same_floats(&x, &rx);
    }
}
#[test]
fn quantDe_matches() {
    let mut rng = XorShift::new(832);
    for i in 0..3000 {
        let mut c = config(i);
        c.nblock = [1, 2, 3, 7, 64][i % 5];
        c.lblock = [1, 2, 7, 31, 64][(i / 5) % 5];
        let mut q = quant_state(&c);
        q.zeroflag = (i % 3) as u8;
        q.ave = rng.uniform_f32(-1000.0, 1000.0);
        q.step = [0.0, -1.0, 1e-20, 1e20, 3.14159][i % 5];
        for row in &mut q.qx {
            for v in row {
                *v = rng.next_u32() as i32;
            }
        }
        let mut x = vec![vec![12.0; c.lblock as usize]; c.nblock as usize];
        let mut rx = x.clone();
        let mut xp = ptrs(&mut x);
        let mut qp = vec![];
        let mut state = cq(&mut q, &mut qp);
        unsafe {
            quantDe(xp.as_mut_ptr(), &c, &mut state);
        }
        rust::quantDe(&mut rx, &c, &q);
        same_floats(&x, &rx);
    }
}
#[test]
fn wavePack1_matches() {
    let mut rng = XorShift::new(841);
    for i in 0..3000 {
        let p = 1 + i % 11;
        let stage = (i / 11) % (p + 1);
        let kind = (i % 2) as i32;
        let mut c = WpcConfig::default();
        c.ishift = [-16, -4, 0, 2, 12][i % 5];
        let mut x = floats(&mut rng, 1 << p, 1, i).remove(0);
        let mut y = floats(&mut rng, 1 << p, 1, i + 3).remove(0);
        let mut rx = x.clone();
        let mut ry = y.clone();
        let mut f = CWaveFilter {
            len: 12,
            ishift: c.ishift,
            filterh: c.filterh.as_mut_ptr(),
            filterg: c.filterg.as_mut_ptr(),
        };
        let rf = WaveFilter {
            ishift: c.ishift,
            filterh: &c.filterh,
            filterg: &c.filterg,
        };
        unsafe {
            wavePack1(
                x.as_mut_ptr(),
                y.as_mut_ptr(),
                &mut f,
                p as i32,
                stage as i32,
                kind,
            );
        }
        rust::wavePack1(&mut rx, &mut ry, &rf, p, stage, kind);
        same_floats(&[x, y], &[rx, ry]);
    }
}
#[test]
fn wavePack2_matches() {
    let mut rng = XorShift::new(842);
    for i in 0..3000 {
        let c = config(i);
        let mut x = floats(&mut rng, c.tileszt as usize, c.tileszx as usize, i);
        let mut y = floats(&mut rng, c.lblock as usize, c.nblock as usize, i + 3);
        let mut rx = x.clone();
        let mut ry = y.clone();
        let mut xp = ptrs(&mut x);
        let mut yp = ptrs(&mut y);
        let kind = (i % 2) as i32;
        unsafe {
            wavePack2(xp.as_mut_ptr(), yp.as_mut_ptr(), &c, kind);
        }
        rust::wavePack2(&mut rx, &mut ry, &c, kind);
        same_floats(&x, &rx);
        same_floats(&y, &ry);
    }
}

fn coding_case(case: usize, rng: &mut XorShift) {
    let mut c = config(case);
    // The header needs 5*nblock bytes in addition to coefficient bytes.
    // lblock >= 8 satisfies the C fixed-buffer capacity assumptions.
    c.lblock = c.lblock.max(8);
    c.tileszt = c.nblock * c.lblock;
    c.tileszx = 1;
    let mut q = quant_state(&c);
    q.zeroflag = if case % 11 == 0 { 1 } else { 0 };
    q.ave = rng.uniform_f32(-100.0, 100.0);
    q.step = rng.uniform_f32(0.01, 100.0);
    for b in 0..c.nblock as usize {
        q.blockind[b] = b as u8;
        q.flowflag[b] = if b % 2 == 0 { 0 } else { 1 };
        for v in &mut q.qx[b] {
            *v = if b % 2 == 0 {
                rng.range_i32(-126, 126)
            } else {
                rng.range_i32(-100000, 100000)
            };
            if case % 7 == 0 {
                *v = 0;
            }
        }
    }
    for b in (1..c.nblock as usize).rev() {
        let j = (rng.next_u32() as usize) % (b + 1);
        q.blockind.swap(b, j);
    }
    let mut qp = vec![];
    let mut state = cq(&mut q, &mut qp);
    let mut out = vec![0xa5; c.tileszt as usize * c.tileszx as usize * 10 + 2048];
    if case % 5 == 0 {
        out.truncate(case % 32);
    }
    let pos = out.len().min(case % 4);
    let mut cb = cbuff(&mut out, pos);
    let mut rb = WpcBuffer {
        code: out.clone(),
        pos,
    };
    let ret = unsafe { codeEncoder(&c, &mut state, &mut cb, std::ptr::null_mut()) };
    assert_eq!(ret, rust::codeEncoder(&c, &q, &mut rb));
    assert_eq!(cb.pos as usize, rb.pos);
    assert_eq!(out, rb.code);
    if ret != -1 {
        let mut input = out[pos..cb.pos as usize].to_vec();
        let mut ci = cbuff(&mut input, 0);
        let mut ri = WpcBuffer {
            code: input.clone(),
            pos: 0,
        };
        let mut dq = quant_state(&c);
        let mut rdq = dq.clone();
        // C omits the decoded flag and the zero-path average assignments.
        // Seed only its state with the encoded values; Rust must recover them
        // from the stream, including when reusing the opposite flag state.
        dq.zeroflag = q.zeroflag;
        dq.ave = q.ave;
        rdq.zeroflag = 1 - q.zeroflag;
        let mut dp = vec![];
        let mut ds = cq(&mut dq, &mut dp);
        let ret = unsafe { codeDecoder(&c, &mut ds, &mut ci, std::ptr::null_mut()) };
        update_q(&mut dq, &ds);
        assert_eq!(ret, rust::codeDecoder(&c, &mut rdq, &mut ri));
        assert_eq!(ci.pos as usize, ri.pos);
        same_quant(&dq, &rdq);
        // Compare decoded coefficients to C rather than the original qx:
        // codeSilence remaps every byte 128 to 129, even inside short values.
    }
}
#[test]
fn codeEncoder_matches() {
    let mut rng = XorShift::new(851);
    for i in 0..2100 {
        coding_case(i, &mut rng);
    }
}
#[test]
fn codeDecoder_matches() {
    let mut rng = XorShift::new(852);
    for i in 2100..4200 {
        coding_case(i, &mut rng);
    }
}

#[test]
fn wpcEncoder_matches() {
    let _lock = SORT_PAIR.lock().unwrap();
    let mut rng = XorShift::new(861);
    for i in 0..2100 {
        let c = config(i);
        let mut x = floats(&mut rng, c.tileszt as usize, c.tileszx as usize, i);
        let mut rx = x.clone();
        let mut xp = ptrs(&mut x);
        let mut out = vec![0xa5; c.tileszt as usize * c.tileszx as usize * 8 + 2048];
        if i % 7 == 0 {
            out.truncate(i % 53);
        }
        let pos = out.len().min(i % 4);
        let mut cb = cbuff(&mut out, pos);
        let mut rb = WpcBuffer {
            code: out.clone(),
            pos,
        };
        let error = [0.001, 0.01, 0.1, 1.0][i % 4];
        let ret = unsafe { wpcEncoder(xp.as_mut_ptr(), error, &c, &mut cb, std::ptr::null_mut()) };
        assert_eq!(ret, rust::wpcEncoder(&mut rx, error, &c, &mut rb));
        assert_eq!(cb.pos as usize, rb.pos);
        assert_eq!(out, rb.code);
        same_floats(&x, &rx);
    }
}

fn non_dc_packet(c: &WpcConfig, rng: &mut XorShift) -> Vec<u8> {
    let mut q = quant_state(c);
    q.zeroflag = 0;
    q.ave = rng.uniform_f32(-5.0, 5.0);
    q.step = rng.uniform_f32(0.1, 2.0);
    for b in 0..c.nblock as usize {
        q.flowflag[b] = 0;
        q.blockind[b] = b as u8;
        for v in &mut q.qx[b] {
            *v = rng.range_i32(-126, 126);
        }
    }
    let mut out = rust::buffAlloc1(c.tileszt as usize * c.tileszx as usize * 8 + 2048);
    assert!(rust::codeEncoder(c, &q, &mut out) > 0);
    out.code.truncate(out.pos);
    out.code
}
#[test]
fn wpcDecoder_matches() {
    let mut rng = XorShift::new(862);
    for i in 0..2100 {
        let c = config(i);
        let mut packed = non_dc_packet(&c, &mut rng);
        if i % 7 == 0 {
            packed = vec![2 + (i % 254) as u8];
        }
        let mut ci = cbuff(&mut packed, 0);
        let mut ri = WpcBuffer {
            code: packed.clone(),
            pos: 0,
        };
        let mut x = floats(&mut rng, c.tileszt as usize, c.tileszx as usize, i);
        let mut rx = x.clone();
        let mut xp = ptrs(&mut x);
        let ret = unsafe { wpcDecoder(xp.as_mut_ptr(), &c, &mut ci, std::ptr::null_mut()) };
        assert_eq!(ret, rust::wpcDecoder(&mut rx, &c, &mut ri));
        assert_eq!(ci.pos as usize, ri.pos);
        same_floats(&x, &rx);
    }
}

// C FILE objects are used only as an independent oracle. Rust I/O uses Vec
// and Cursor. All temporary FILEs are closed, including after read/free tests.
unsafe fn temp_file() -> *mut c_void {
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let id = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../.targets/lib-lib_wpc2");
    let path = root.join(format!("wpc-diff-{}-{id}.tmp", std::process::id()));
    let name = std::ffi::CString::new(path.to_str().unwrap()).unwrap();
    let fp = fopen(name.as_ptr(), b"w+b\0".as_ptr().cast());
    assert!(!fp.is_null());
    std::fs::remove_file(&path).unwrap();
    fp
}
unsafe fn c_read_bytes(bytes: &[u8]) -> *mut CPacked {
    let fp = temp_file();
    assert!(!fp.is_null());
    assert_eq!(
        fwrite(bytes.as_ptr().cast(), 1, bytes.len(), fp),
        bytes.len()
    );
    rewind(fp);
    let w = wpcRead(fp);
    assert!(!w.is_null());
    fclose(fp);
    w
}
unsafe fn c_write_bytes(w: *mut CPacked) -> Vec<u8> {
    let fp = temp_file();
    assert!(!fp.is_null());
    wpcWrite(w, fp);
    rewind(fp);
    let mut bytes = vec![0; 16 + std::mem::size_of::<WpcConfig>() + (*w).size as usize];
    let n = bytes.len();
    assert_eq!(fread(bytes.as_mut_ptr().cast(), 1, n, fp), n);
    fclose(fp);
    bytes
}
fn packed_fixture(case: usize, rng: &mut XorShift) -> WpcPacked {
    WpcPacked {
        version: rng.next_u32() as i32,
        n1: rng.next_u32() as i32,
        n2: rng.next_u32() as i32,
        config: config(case),
        code: byte_data(case, rng),
    }
}
unsafe fn borrowed_c(w: &mut WpcPacked) -> CPacked {
    CPacked {
        version: w.version,
        n1: w.n1,
        n2: w.n2,
        size: w.code.len() as i32,
        config: &mut w.config,
        code: w.code.as_mut_ptr(),
    }
}
#[test]
fn wpcN1_matches() {
    let mut rng = XorShift::new(871);
    for i in 0..3000 {
        let mut w = packed_fixture(i, &mut rng);
        let mut cw = unsafe { borrowed_c(&mut w) };
        assert_eq!(unsafe { wpcN1(&mut cw) }, rust::wpcN1(&w));
    }
}
#[test]
fn wpcN2_matches() {
    let mut rng = XorShift::new(872);
    for i in 0..3000 {
        let mut w = packed_fixture(i, &mut rng);
        let mut cw = unsafe { borrowed_c(&mut w) };
        assert_eq!(unsafe { wpcN2(&mut cw) }, rust::wpcN2(&w));
    }
}
#[test]
fn wpcSize_matches() {
    let mut rng = XorShift::new(873);
    for i in 0..3000 {
        let mut w = packed_fixture(i, &mut rng);
        let mut cw = unsafe { borrowed_c(&mut w) };
        assert_eq!(unsafe { wpcSize(&mut cw) }, rust::wpcSize(&w));
    }
}
#[test]
fn wpcWrite_matches() {
    let mut rng = XorShift::new(874);
    for i in 0..2100 {
        let mut w = packed_fixture(i, &mut rng);
        let mut cw = unsafe { borrowed_c(&mut w) };
        let expected = unsafe { c_write_bytes(&mut cw) };
        let mut actual = vec![];
        rust::wpcWrite(&w, &mut actual).unwrap();
        assert_eq!(expected, actual);
    }
}
#[test]
fn wpcRead_matches() {
    let mut rng = XorShift::new(875);
    for i in 0..2100 {
        let mut w = packed_fixture(i, &mut rng);
        let mut cw = unsafe { borrowed_c(&mut w) };
        let bytes = unsafe { c_write_bytes(&mut cw) };
        let cr = unsafe { c_read_bytes(&bytes) };
        let rr = rust::wpcRead(&mut std::io::Cursor::new(&bytes)).unwrap();
        unsafe {
            assert_eq!(
                ((*cr).version, (*cr).n1, (*cr).n2, (*cr).size),
                (rr.version, rr.n1, rr.n2, rr.code.len() as i32)
            );
            assert_eq!(*(*cr).config, rr.config);
            assert_eq!(
                std::slice::from_raw_parts((*cr).code, (*cr).size as usize),
                rr.code
            );
            wpcFree(cr);
        }
    }
}
#[test]
fn wpcFree_matches() {
    let mut rng = XorShift::new(876);
    for i in 0..2100 {
        let w = packed_fixture(i, &mut rng);
        let mut bytes = vec![];
        rust::wpcWrite(&w, &mut bytes).unwrap();
        let cr = unsafe { c_read_bytes(&bytes) };
        unsafe {
            assert_eq!((*cr).size, rust::wpcSize(&w));
            wpcFree(cr);
        }
        rust::wpcFree(w);
    }
}
#[test]
fn wpcUncompress_matches() {
    let mut rng = XorShift::new(877);
    for i in 0..2100 {
        let c = config(i);
        let n1 = [1, 2, 3, 7, 17, 32, 33][i % 7];
        let n2 = [1, 2, 5, 16, 17][(i / 7) % 5];
        let mut code = vec![];
        for _ in 0..((n2 - 1) / c.tileszx as usize + 1) {
            for _ in 0..((n1 - 1) / c.tileszt as usize + 1) {
                code.extend(non_dc_packet(&c, &mut rng));
            }
        }
        let mut w = WpcPacked {
            version: 1,
            n1: n1 as i32,
            n2: n2 as i32,
            config: c,
            code,
        };
        if i % 11 == 0 {
            w.n1 = -(i as i32);
        }
        if i % 13 == 0 {
            w.code = vec![255];
        }
        let mut cw = unsafe { borrowed_c(&mut w) };
        let mut expected = vec![13.0; n1 * n2];
        let mut actual = expected.clone();
        let ret = unsafe { wpcUncompress(&mut cw, expected.as_mut_ptr()) };
        assert_eq!(ret, rust::wpcUncompress(&w, &mut actual));
        same_floats(&[expected], &[actual]);
    }
}

#[test]
fn codeDecoder_restores_zero_tile_header() {
    let mut rng = XorShift::new(879);
    for case in 0..3000 {
        let c = config(case);
        let ave = match case % 7 {
            0 => 0.0,
            1 => -0.0,
            2 => f32::MIN_POSITIVE,
            3 => -f32::MIN_POSITIVE,
            _ => rng.uniform_f32(-1000.0, 1000.0),
        };
        let pos = case % 4;
        let mut bytes = vec![0xa5; pos];
        bytes.push(1);
        bytes.extend_from_slice(&ave.to_ne_bytes());
        let mut ci = cbuff(&mut bytes, pos);
        let mut ri = WpcBuffer {
            code: bytes.clone(),
            pos,
        };
        let mut rq = quant_state(&c);
        // Neither the step nor coefficients are transmitted for a zero tile.
        // Poison the step so accidentally taking the nonzero path is visible.
        rq.step = f32::NAN;
        let mut expected = rq.clone();
        // The original C decoder omits these two assignments. Preinitialize
        // only the C state, as in coding_case, to compare all remaining state.
        expected.zeroflag = 1;
        expected.ave = ave;
        let mut qp = vec![];
        let mut cq = cq(&mut expected, &mut qp);
        assert_eq!(
            unsafe { codeDecoder(&c, &mut cq, &mut ci, std::ptr::null_mut()) },
            1
        );
        update_q(&mut expected, &cq);
        assert_eq!(rust::codeDecoder(&c, &mut rq, &mut ri), 1);
        assert_eq!(ri.pos, bytes.len());
        assert_eq!(ci.pos as usize, ri.pos);
        same_quant(&expected, &rq);

        // Keep quantDe's original zero-flag behavior: it clears every
        // coefficient, independently of the unused step and integer values.
        let mut x = vec![vec![f32::NAN; c.lblock as usize]; c.nblock as usize];
        rust::quantDe(&mut x, &c, &rq);
        assert!(x.iter().flatten().all(|&v| v.to_bits() == 0));
    }
}

#[test]
fn codeDecoder_rejects_incomplete_zero_tile_headers() {
    let c = WpcConfig::default();
    let mut packet = vec![1];
    packet.extend_from_slice(&(-3.25f32).to_ne_bytes());
    for len in 0..packet.len() {
        let mut bytes = packet[..len].to_vec();
        let mut ci = cbuff(&mut bytes, 0);
        let mut ri = WpcBuffer {
            code: bytes.clone(),
            pos: 0,
        };
        let mut q = quant_state(&c);
        let mut rq = q.clone();
        let mut qp = vec![];
        let mut cq = cq(&mut q, &mut qp);
        assert_eq!(
            unsafe { codeDecoder(&c, &mut cq, &mut ci, std::ptr::null_mut()) },
            0
        );
        assert_eq!(rust::codeDecoder(&c, &mut rq, &mut ri), 0);
        assert_eq!(ci.pos as usize, ri.pos);
        same_quant(&q, &rq);
    }
}

// Compare full compressed bytes through C's block encoder, avoiding the stale
// pointer returned by C wpcCompress. These fixtures fit the initial buffer, so
// both encoders make the same calls to the persistent qkisort generators.
// The caller holds SORT_PAIR for the entire comparison.
fn compress_zero_tile_fixture(f: &[f32], n1: usize, n2: usize) -> WpcPacked {
    let c = WpcConfig::default();
    let (nt, nx) = (c.tileszt as usize, c.tileszx as usize);
    let mut bytes = vec![0; n1 * n2];
    let mut cb = cbuff(&mut bytes, 0);
    for ix in (0..n2).step_by(nx) {
        for it in (0..n1).step_by(nt) {
            let mut block = vec![vec![0.0; nt]; nx];
            for j in 0..nx.min(n2 - ix) {
                let len = nt.min(n1 - it);
                let start = (ix + j) * n1 + it;
                block[j][..len].copy_from_slice(&f[start..start + len]);
            }
            let mut bp = ptrs(&mut block);
            assert!(
                unsafe { wpcEncoder(bp.as_mut_ptr(), 0.1, &c, &mut cb, std::ptr::null_mut()) } > 0
            );
        }
    }
    let packed = rust::wpcCompress(f, n1, n2, 0.1).unwrap();
    assert_eq!(packed.code, bytes[..cb.pos as usize]);
    packed
}

#[test]
fn wpcUncompress_preserves_muted_tile() {
    let _lock = SORT_PAIR.lock().unwrap();
    let mut rng = XorShift::new(880);
    let (n1, n2) = (512, 256);
    let mut input = vec![0.0; n1 * n2];
    for trace in input.chunks_exact_mut(n1) {
        for v in &mut trace[256..] {
            *v = rng.uniform_f32(-100.0, 100.0);
        }
    }
    let packed = compress_zero_tile_fixture(&input, n1, n2);
    // The first tile has only its zero flag and average; the next is nonzero.
    assert_eq!(&packed.code[..5], &[1, 0, 0, 0, 0]);
    assert_eq!(packed.code[5], 0);
    let mut decoded = vec![f32::NAN; input.len()];
    assert_eq!(rust::wpcUncompress(&packed, &mut decoded), 1);
    assert!(decoded.iter().all(|v| v.is_finite()));
    for trace in decoded.chunks_exact(n1) {
        assert!(trace[..256].iter().all(|&v| v == 0.0));
    }
    assert!(decoded
        .chunks_exact(n1)
        .any(|trace| trace[256..].iter().any(|&v| v != 0.0)));

    // Compare the live tile against the original C decoder independently;
    // its zero-tile output is undefined and is deliberately not an oracle.
    let mut bytes = packed.code[5..].to_vec();
    let mut cb = cbuff(&mut bytes, 0);
    let mut expected = vec![vec![0.0; 256]; 256];
    let mut ep = ptrs(&mut expected);
    assert_eq!(
        unsafe {
            wpcDecoder(
                ep.as_mut_ptr(),
                &packed.config,
                &mut cb,
                std::ptr::null_mut(),
            )
        },
        1
    );
    let actual: Vec<_> = decoded
        .chunks_exact(n1)
        .map(|trace| trace[256..].to_vec())
        .collect();
    assert!(expected.iter().flatten().all(|v| v.is_finite()));
    same_floats(&expected, &actual);
}

#[test]
fn wpcUncompress_preserves_zero_edge_tiles() {
    let _lock = SORT_PAIR.lock().unwrap();
    let mut rng = XorShift::new(881);
    for (n1, n2) in [(257, 259), (256, 512), (512, 512)] {
        let mut input = vec![0.0; n1 * n2];
        for trace in input.chunks_exact_mut(n1).take(256) {
            for v in &mut trace[..256] {
                *v = rng.uniform_f32(-100.0, 100.0);
            }
        }
        let packed = compress_zero_tile_fixture(&input, n1, n2);
        let mut decoded = vec![f32::NAN; input.len()];
        assert_eq!(rust::wpcUncompress(&packed, &mut decoded), 1);
        assert!(decoded.iter().all(|v| v.is_finite()));
        for (j, trace) in decoded.chunks_exact(n1).enumerate() {
            let zeros = if j < 256 { &trace[256..] } else { trace };
            assert!(zeros.iter().all(|&v| v == 0.0));
        }
        assert!(decoded.iter().any(|&v| v != 0.0));
    }
}

#[test]
fn wpcUncompress_all_zero_sections() {
    let _lock = SORT_PAIR.lock().unwrap();
    for (n1, n2) in [(1, 1), (2, 2), (7, 31), (256, 256), (512, 256), (257, 259)] {
        let input = vec![0.0; n1 * n2];
        // Zero-tile quantization returns before sorting, so these calls do
        // not advance the Rust qkisort seed or require a paired C compressor.
        let packed = rust::wpcCompress(&input, n1, n2, 0.1).unwrap();
        let mut decoded = vec![f32::NAN; input.len()];
        assert_eq!(rust::wpcUncompress(&packed, &mut decoded), 1);
        assert!(decoded.iter().all(|&v| v == 0.0));
    }
}

/// Smooth data compresses well, so the encoder never overflows its buffer and
/// never retries: here the Rust must match SU exactly. Inputs that overflow the
/// buffer are tested Rust-only in tests/roundtrip.rs: there OxWave deliberately
/// differs from SU (DEVIATIONS.md), and calling only one side would desynchronise
/// the paired qkisort seeds that the other tests in this file rely on.
#[test]
fn wpcCompress_matches() {
    let _lock = SORT_PAIR.lock().unwrap_or_else(|e| e.into_inner());
    let mut n = 0;
    for (k, &(n1, n2)) in [(256usize, 17usize), (257, 64), (512, 33), (300, 257)].iter().enumerate() {
        for &error in &[0.01f32, 0.1] {
            let mut f: Vec<f32> = (0..n1 * n2)
                .map(|idx| {
                    let (i, j) = ((idx % n1) as f32, (idx / n1) as f32);
                    (0.05 * i + 0.02 * j + k as f32).sin() + 0.5 * (0.013 * i * (1.0 + 0.01 * j)).cos()
                })
                .collect();
            let mut size = 0;
            let cw = unsafe { wpcCompress(f.as_mut_ptr(), n1 as i32, n2 as i32, error, &mut size) };
            let rw = rust::wpcCompress(&f, n1, n2, error).expect("Rust compressor consistency");
            assert!(!cw.is_null());
            assert!(rust::wpcSize(&rw) as usize <= n1 * n2, "smooth case {n1}x{n2} error={error} hit the retry path");
            unsafe {
                assert_eq!(size, rust::wpcSize(&rw), "size for {n1}x{n2} error={error}");
                assert_eq!(*(*cw).config, rw.config);
                // (*cw).code is stale after wpc.c's final realloc (SU bug), so it is neither read
                // nor freed; exact stream equality is covered by wpcEncoder_matches.
                free((*cw).config.cast());
                free(cw.cast());
            }
            let mut g = vec![0.0f32; n1 * n2];
            rust::wpcUncompress(&rw, &mut g);
            let num: f64 = f.iter().zip(&g).map(|(a, b)| ((a - b) as f64).powi(2)).sum();
            let den: f64 = f.iter().map(|a| (*a as f64).powi(2)).sum();
            assert!((num / den).sqrt() <= 3.0 * error as f64, "round trip {n1}x{n2} error={error}");
            rust::wpcFree(rw);
            n += 1;
        }
    }
    assert_eq!(n, 8);
}
