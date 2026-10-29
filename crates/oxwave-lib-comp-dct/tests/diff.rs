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

//! Differential tests against comp/dct/lib/{buffalloc,dct1,dct2,dctalloc,
//! getfilter,huffman,lct1,lpred,pcoding,quant,rle,wavepack1,wavepack2,
//! wavetrans1,wavetrans2}.c. Credits: Tong Chen (1994–1995), Huffman
//! algorithm based on Mark Nelson's "Data Compression" book.
#![allow(non_snake_case)]

use oxwave_cref::XorShift;
use oxwave_lib_comp_dct as rust;

#[repr(C)]
struct CBuffer {
    bound: i32,
    pos: i32,
    code: *mut u8,
}
#[repr(C)]
struct CRc {
    rlen: i32,
    rlast: i32,
    rc: *mut f32,
}
#[repr(C)]
struct CWave {
    id: i32,
    len: i32,
    ishift: i32,
    h: *mut f32,
    g: *mut f32,
}
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
struct CLp {
    hlfilter: i32,
    lfilter: i32,
    nslope: i32,
    nblock: i32,
    filter: *mut *mut *mut f32,
}

#[link(name = "comp", kind = "static")]
unsafe extern "C" {
    fn buffAlloc1(n: i32) -> *mut CBuffer;
    fn buffAlloc2(n1: i32, n2: i32) -> *mut *mut CBuffer;
    fn buffFree1(b: *mut CBuffer);
    fn buffFree2(b: *mut *mut CBuffer, n1: i32);
    fn buffRealloc1(b: *mut CBuffer, n: i32);
    fn codeSilence(input: *mut std::ffi::c_void, output: *mut std::ffi::c_void);
    fn codeDesilence(input: *mut std::ffi::c_void, output: *mut std::ffi::c_void);
    fn dct(x: *mut f32, n: i32, c: *mut *mut f32, typ: i32);
    fn dctAlloc(n: i32) -> *mut *mut f32;
    fn dct_2(x: *mut *mut f32, n1: i32, n2: i32, c1: *mut *mut f32, c2: *mut *mut f32, typ: i32);
    fn dctiv(x: *mut f32, n: i32, c: *mut *mut f32);
    fn dctivAlloc(n: i32) -> *mut *mut f32;
    fn dctiv_2(x: *mut *mut f32, n1: i32, n2: i32, c1: *mut *mut f32, c2: *mut *mut f32);
    fn huffCompress(input: *mut std::ffi::c_void, output: *mut std::ffi::c_void) -> i32;
    fn huffDecompress(input: *mut std::ffi::c_void, output: *mut std::ffi::c_void) -> i32;
    fn ilct(sig: *mut f32, n: i32, rc: *mut CRc, c: *mut *mut f32);
    fn lct(sig: *mut f32, n: i32, rc: *mut CRc, c: *mut *mut f32);
    fn lpred(
        f: *mut *mut f32,
        n1: i32,
        n2: i32,
        lblock: i32,
        hlfilter: i32,
        nslope: i32,
        lpf: *mut CLp,
        res: *mut *mut f32,
    );
    fn pDecode(b: *mut CBuffer, qx: *mut i32, n: i32) -> i32;
    fn pEncode(qx: *mut i32, n: i32) -> *mut CBuffer;
    fn rcalloc(n: i32) -> *mut CRc;
    fn rcfmidp(rc: *mut CRc);
    fn uniDequant(x: *mut f32, n: i32, ave: f32, step: f32, qx: *mut i32);
    fn uniQuant(x: *mut f32, n: i32, error: f32, ave: *mut f32, step: *mut f32, qx: *mut i32);
    fn waveFreefilter(filter: *mut CWave);
    fn waveGetfilter(id: i32) -> *mut CWave;
    fn wavePack_1(
        x: *mut f32,
        y: *mut f32,
        filter: *mut CWave,
        pow: i32,
        stage: i32,
        typ: i32,
    ) -> i32;
    fn wavePack_2(
        x: *mut *mut f32,
        y: *mut *mut f32,
        f1: *mut CWave,
        f2: *mut CWave,
        p1: i32,
        p2: i32,
        s1: i32,
        s2: i32,
        typ: i32,
    ) -> i32;
    fn waveTrans_1(
        x: *mut f32,
        y: *mut f32,
        filter: *mut CWave,
        pow: i32,
        stage: i32,
        typ: i32,
    ) -> i32;
    fn waveTrans_2(
        x: *mut *mut f32,
        y: *mut *mut f32,
        filter: *mut CWave,
        p1: i32,
        p2: i32,
        s1: i32,
        s2: i32,
        typ: i32,
    ) -> i32;
    fn free2float(p: *mut *mut f32);
    fn free(p: *mut std::ffi::c_void);
}

const CASES: usize = 2500;

fn close(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    let scale = b
        .iter()
        .filter(|x| x.is_finite())
        .fold(0.0f32, |m, x| m.max(x.abs()));
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        if x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()) {
            continue;
        }
        assert!(
            x.is_finite() && y.is_finite() && (x - y).abs() <= 1e-5 * scale + 1e-30,
            "sample {i}: {x} != {y}, scale {scale}"
        );
    }
}

fn matrix_close(a: &[Vec<f32>], b: &[Vec<f32>]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        close(x, y);
    }
}

fn ptrs(x: &mut [Vec<f32>]) -> Vec<*mut f32> {
    x.iter_mut().map(|r| r.as_mut_ptr()).collect()
}

fn size(rng: &mut XorShift, iter: usize) -> usize {
    const EDGES: [usize; 13] = [1, 2, 3, 7, 8, 16, 31, 64, 127, 257, 1024, 4093, 8192];
    if iter % 50 < EDGES.len() {
        EDGES[iter % 50]
    } else {
        rng.range_i32(1, 300) as usize
    }
}

fn signal(rng: &mut XorShift, n: usize, iter: usize) -> Vec<f32> {
    let scale = [1.0, 1e-12, 1e12, 300.0][iter % 4];
    let mut x = rng.vec_f32(n, -scale, scale);
    match iter % 19 {
        0 => x.fill(0.0),
        1 => x.fill(-scale),
        2 => {
            x.fill(0.0);
            x[n / 2] = -scale;
        }
        3 => {
            for (j, v) in x.iter_mut().enumerate() {
                *v = if j % 2 == 0 { scale } else { -scale };
            }
        }
        _ => (),
    }
    x
}

unsafe fn bytes<'a>(p: *mut CBuffer) -> &'a [u8] {
    // malloc(0) and realloc(p,0) may return null; do not construct a null slice.
    if (*p).pos == 0 {
        &[]
    } else {
        std::slice::from_raw_parts((*p).code, (*p).pos as usize)
    }
}

unsafe fn input_buffer(data: &[u8], pos: usize) -> *mut CBuffer {
    let p = buffAlloc1(data.len() as i32);
    if !data.is_empty() {
        std::ptr::copy_nonoverlapping(data.as_ptr(), (*p).code, data.len());
    }
    (*p).pos = pos as i32;
    p
}

// Allocate/free tests initialize all C storage before examining it: C malloc
// leaves bytes uninitialized, whereas the safe Rust allocation zeroes them.
fn buffer_one(reallocate: bool) {
    unsafe {
        let mut rng = XorShift::new(12);
        for iter in 0..CASES {
            let n = if iter % 29 == 0 {
                0
            } else {
                size(&mut rng, iter)
            };
            let p = buffAlloc1(n as i32);
            let mut r = rust::buffAlloc1(n);
            assert_eq!(((*p).bound, (*p).pos), (n as i32, 0));
            assert_eq!((r.code.len(), r.pos), (n, 0));
            for i in 0..n {
                let c = rng.next_u32() as u8;
                *(*p).code.add(i) = c;
                r.code[i] = c;
            }
            let pos = rng.range_i32(0, n as i32) as usize;
            (*p).pos = pos as i32;
            r.pos = pos;
            if reallocate {
                let nn = if iter % 3 == 0 {
                    0
                } else {
                    size(&mut rng, iter + 1)
                };
                buffRealloc1(p, nn as i32);
                rust::buffRealloc1(&mut r, nn);
                assert_eq!(
                    ((*p).bound as usize, (*p).pos as usize),
                    (r.code.len(), r.pos)
                );
                for i in 0..n.min(nn) {
                    assert_eq!(*(*p).code.add(i), r.code[i]);
                }
            }
            buffFree1(p);
            rust::buffFree1(r);
        }
    }
}

fn buffer_two() {
    unsafe {
        let mut rng = XorShift::new(15);
        for iter in 0..CASES {
            let n1 = rng.range_i32(0, 9) as usize;
            let n2 = if iter % 31 == 0 {
                0
            } else {
                size(&mut rng, iter)
            };
            let p = buffAlloc2(n1 as i32, n2 as i32);
            let mut r = rust::buffAlloc2(n1, n2);
            assert_eq!(r.len(), n1);
            for i in 0..n1 {
                let row = *p.add(i);
                assert_eq!(
                    ((*row).bound as usize, (*row).pos as usize),
                    (r[i].code.len(), r[i].pos)
                );
                for j in 0..n2 {
                    let c = rng.next_u32() as u8;
                    *(*row).code.add(j) = c;
                    r[i].code[j] = c;
                }
            }
            buffFree2(p, n1 as i32);
            rust::buffFree2(r);
        }
    }
}

#[test]
fn buffAlloc1_matches_c() {
    buffer_one(false);
}
#[test]
fn buffAlloc2_matches_c() {
    buffer_two();
}
#[test]
fn buffFree1_matches_c() {
    buffer_one(false);
}
#[test]
fn buffFree2_matches_c() {
    buffer_two();
}
#[test]
fn buffRealloc1_matches_c() {
    buffer_one(true);
}

fn dct_table(iv: bool) {
    unsafe {
        let mut rng = XorShift::new(16);
        for iter in 0..CASES {
            // Thousands of table allocations, including primes and a large table.
            let n = if iter == 0 {
                512
            } else {
                rng.range_i32(1, 65) as usize
            };
            let p = if iv {
                dctivAlloc(n as i32)
            } else {
                dctAlloc(n as i32)
            };
            let r = if iv {
                rust::dctivAlloc(n)
            } else {
                rust::dctAlloc(n)
            };
            for i in 0..n {
                close(&r[i], std::slice::from_raw_parts(*p.add(i), n));
            }
            free2float(p);
        }
    }
}
#[test]
fn dctAlloc_matches_c() {
    dct_table(false);
}
#[test]
fn dctivAlloc_matches_c() {
    dct_table(true);
}

fn dct_one(iv: bool) {
    unsafe {
        let mut rng = XorShift::new(17);
        for n in [1, 2, 3, 7, 16, 31, 64, 257, 1024] {
            let p = if iv {
                dctivAlloc(n as i32)
            } else {
                dctAlloc(n as i32)
            };
            let r = if iv {
                rust::dctivAlloc(n)
            } else {
                rust::dctAlloc(n)
            };
            let count = if n <= 64 { 360 } else { 4 };
            for iter in 0..count {
                let mut a = signal(&mut rng, n, iter);
                let mut b = a.clone();
                let typ = [0, 1, -1, 2][iter % 4];
                if iv {
                    dctiv(a.as_mut_ptr(), n as i32, p);
                    rust::dctiv(&mut b, &r);
                } else {
                    dct(a.as_mut_ptr(), n as i32, p, typ);
                    rust::dct(&mut b, &r, typ);
                }
                close(&b, &a);
            }
            free2float(p);
        }
    }
}
#[test]
fn dct_matches_c() {
    dct_one(false);
}
#[test]
fn dctiv_matches_c() {
    dct_one(true);
}

fn dct_two(iv: bool) {
    unsafe {
        let mut rng = XorShift::new(18);
        for (n1, n2) in [(1, 1), (1, 7), (2, 3), (7, 2), (16, 8), (31, 17), (127, 64)] {
            let p1 = if iv { dctivAlloc(n1) } else { dctAlloc(n1) };
            let p2 = if iv { dctivAlloc(n2) } else { dctAlloc(n2) };
            let r1 = if iv {
                rust::dctivAlloc(n1 as usize)
            } else {
                rust::dctAlloc(n1 as usize)
            };
            let r2 = if iv {
                rust::dctivAlloc(n2 as usize)
            } else {
                rust::dctAlloc(n2 as usize)
            };
            for iter in 0..if n1 > 64 { 4 } else { 420 } {
                let mut a: Vec<_> = (0..n2)
                    .map(|_| signal(&mut rng, n1 as usize, iter))
                    .collect();
                let mut b = a.clone();
                let mut pp = ptrs(&mut a);
                let typ = [0, 1, -1][iter % 3];
                if iv {
                    dctiv_2(pp.as_mut_ptr(), n1, n2, p1, p2);
                    rust::dctiv_2(&mut b, &r1, &r2);
                } else {
                    dct_2(pp.as_mut_ptr(), n1, n2, p1, p2, typ);
                    rust::dct_2(&mut b, &r1, &r2, typ);
                }
                matrix_close(&b, &a);
            }
            free2float(p1);
            free2float(p2);
        }
    }
}
#[test]
fn dct_2_matches_c() {
    dct_two(false);
}
#[test]
fn dctiv_2_matches_c() {
    dct_two(true);
}

#[test]
fn uniQuant_matches_c() {
    unsafe {
        let mut rng = XorShift::new(19);
        for iter in 0..CASES {
            let n = size(&mut rng, iter);
            let mut x = signal(&mut rng, n, iter);
            let error = [0.03, 0.001, 0.5, 2.0, -0.1, 0.0][iter % 6];
            let mut sa = if iter % 17 == 0 {
                0.0
            } else if iter % 2 == 0 {
                -1.0
            } else {
                rng.uniform_f32(0.1, 1000.0)
            };
            let mut sb = sa;
            let mut aa = -999.0;
            let mut ab = aa;
            let mut qa = vec![0; n];
            let mut qb = qa.clone();
            uniQuant(
                x.as_mut_ptr(),
                n as i32,
                error,
                &mut aa,
                &mut sa,
                qa.as_mut_ptr(),
            );
            rust::uniQuant(&x, error, &mut ab, &mut sb, &mut qb);
            assert_eq!(qa, qb, "case {iter}, n={n}");
            close(&[ab], &[aa]);
            close(&[sb], &[sa]);
        }
    }
}

#[test]
fn uniDequant_matches_c() {
    unsafe {
        let mut rng = XorShift::new(20);
        for iter in 0..CASES {
            let n = size(&mut rng, iter);
            let mut q: Vec<_> = (0..n).map(|_| rng.next_u32() as i32).collect();
            if iter % 11 == 0 {
                q.fill(0);
            }
            let ave = rng.uniform_f32(-10000.0, 10000.0);
            let step = [0.0, -0.0, 0.01, -0.03, 1e-20, 1e20, f32::INFINITY][iter % 7];
            let mut a = vec![0.0; n];
            let mut b = a.clone();
            uniDequant(a.as_mut_ptr(), n as i32, ave, step, q.as_mut_ptr());
            rust::uniDequant(&mut b, ave, step, &q);
            close(&b, &a);
        }
    }
}

fn cutoff(fill: bool) {
    unsafe {
        let mut rng = XorShift::new(21);
        for iter in 0..CASES {
            let n = if iter == 0 { 0 } else { size(&mut rng, iter) };
            let p = rcalloc(n as i32);
            let mut r = rust::rcalloc(n);
            assert_eq!(
                ((*p).rlen as usize, (*p).rlast as usize),
                (r.rc.len(), r.rlast)
            );
            if fill {
                rcfmidp(p);
                rust::rcfmidp(&mut r);
                if n != 0 {
                    close(&r.rc, std::slice::from_raw_parts((*p).rc, 2 * n));
                }
            } else {
                for i in 0..2 * n {
                    let v = rng.uniform_f32(-2.0, 2.0);
                    *(*p).rc.add(i) = v;
                    r.rc[i] = v;
                }
                if n != 0 {
                    close(&r.rc, std::slice::from_raw_parts((*p).rc, 2 * n));
                }
            }
            free((*p).rc.cast());
            free(p.cast());
        }
    }
}
#[test]
fn rcalloc_matches_c() {
    cutoff(false);
}
#[test]
fn rcfmidp_matches_c() {
    cutoff(true);
}

fn local_cosine(inverse: bool) {
    unsafe {
        let mut rng = XorShift::new(22);
        for n in [1usize, 2, 3, 7, 16, 31, 64, 257] {
            let p = dctivAlloc(n as i32);
            let r = rust::dctivAlloc(n);
            for iter in 0..if n > 64 { 4 } else { 360 } {
                // Include overlapping folding regions and zero-width halos.
                let h = rng.range_i32(0, n as i32) as usize;
                let rc = rcalloc(h as i32);
                rcfmidp(rc);
                let mut rr = rust::rcalloc(h);
                rust::rcfmidp(&mut rr);
                let mut a = signal(&mut rng, n + 2 * h, iter);
                let mut b = a.clone();
                if inverse {
                    ilct(a.as_mut_ptr().add(h), n as i32, rc, p);
                    rust::ilct(&mut b, h, n, &rr, &r);
                } else {
                    lct(a.as_mut_ptr().add(h), n as i32, rc, p);
                    rust::lct(&mut b, h, n, &rr, &r);
                }
                close(&b, &a);
                free((*rc).rc.cast());
                free(rc.cast());
            }
            free2float(p);
        }
    }
}
#[test]
fn lct_matches_c() {
    local_cosine(false);
}
#[test]
fn ilct_matches_c() {
    local_cosine(true);
}

fn prefix(decode: bool) {
    unsafe {
        let mut rng = XorShift::new(23);
        let edges = [
            i32::MIN,
            i32::MAX,
            -32768,
            -32767,
            -128,
            -127,
            -126,
            -1,
            0,
            1,
            126,
            127,
            128,
            32767,
            32768,
        ];
        for iter in 0..CASES {
            let n = if iter % 37 == 0 {
                0
            } else {
                size(&mut rng, iter)
            };
            let mut q: Vec<_> = (0..n)
                .map(|j| {
                    if iter % 4 == 0 {
                        edges[j % edges.len()]
                    } else if iter % 4 == 1 {
                        rng.range_i32(-126, 126)
                    } else {
                        rng.next_u32() as i32
                    }
                })
                .collect();
            let p = pEncode(q.as_mut_ptr(), n as i32);
            let mut r = rust::pEncode(&q);
            assert_eq!(
                ((*p).bound as usize, (*p).pos as usize),
                (r.code.len(), r.pos)
            );
            assert_eq!(bytes(p), r.code);
            if decode {
                (*p).pos = 0;
                r.pos = 0;
                let count = n + iter % 2; // also test the C consistency failure.
                let mut a = vec![i32::MIN; count + 1];
                let mut b = a.clone();
                let ca = pDecode(p, a.as_mut_ptr(), count as i32);
                let cb = rust::pDecode(&mut r, &mut b[..count]);
                assert_eq!(ca, cb);
                assert_eq!(a, b);
                assert_eq!((*p).pos as usize, r.pos);
            }
            buffFree1(p);
            if decode {
                // Arbitrary and truncated prefix streams exercise stale c1/c2
                // behavior; both implementations are given sufficient output.
                let count = rng.range_i32(1, 100) as usize;
                let mut data: Vec<_> = (0..count).map(|_| rng.next_u32() as u8).collect();
                data[count - 1] = if iter % 2 == 0 { 127 } else { 129 };
                let pos = iter % count;
                let p = input_buffer(&data, pos);
                let mut r = rust::WpcBuffer { code: data, pos };
                let mut a = vec![0x7654321; count + 1];
                let mut b = a.clone();
                assert_eq!(
                    pDecode(p, a.as_mut_ptr(), count as i32),
                    rust::pDecode(&mut r, &mut b[..count])
                );
                assert_eq!(a, b);
                assert_eq!((*p).pos as usize, r.pos);
                buffFree1(p);
            }
        }
    }
}
#[test]
fn pEncode_matches_c() {
    prefix(false);
}
#[test]
fn pDecode_matches_c() {
    prefix(true);
}

fn byte_input(rng: &mut XorShift, n: usize, iter: usize) -> Vec<u8> {
    (0..n)
        .map(|i| match iter % 6 {
            0 => 0,
            1 => 128,
            2 => (i % 256) as u8,
            3 => {
                if rng.next_u32() % 10 != 0 {
                    0
                } else {
                    rng.next_u32() as u8
                }
            }
            4 => 255,
            _ => rng.next_u32() as u8,
        })
        .collect()
}

fn silence(decode: bool) {
    unsafe {
        let mut rng = XorShift::new(24);
        for iter in 0..CASES {
            let n = if iter % 37 == 0 {
                0
            } else {
                size(&mut rng, iter)
            };
            let data = byte_input(&mut rng, n, iter);
            let pos = if n == 0 || iter % 2 == 0 {
                0
            } else {
                rng.range_i32(0, n as i32) as usize
            };
            let capacity = if iter % 3 == 0 {
                rng.range_i32(0, n as i32) as usize
            } else {
                n * 255 + 8
            };
            let outpos = rng.range_i32(0, capacity as i32) as usize;
            let p = input_buffer(&data, pos);
            let out = vec![0xa5; capacity];
            let cp = input_buffer(&out, outpos);
            let mut r = rust::WpcBuffer { code: data, pos };
            let mut ro = rust::WpcBuffer {
                code: out,
                pos: outpos,
            };
            if decode {
                codeDesilence(p.cast(), cp.cast());
                rust::codeDesilence(&mut r, &mut ro);
            } else {
                codeSilence(p.cast(), cp.cast());
                rust::codeSilence(&mut r, &mut ro);
            }
            assert_eq!((*p).pos as usize, r.pos);
            assert_eq!((*cp).pos as usize, ro.pos);
            assert_eq!(bytes(cp), &ro.code[..ro.pos]);
            if capacity > 0 {
                assert_eq!(std::slice::from_raw_parts((*cp).code, capacity), ro.code);
            }
            buffFree1(p);
            buffFree1(cp);
        }
    }
}
#[test]
fn codeSilence_matches_c() {
    silence(false);
}
#[test]
fn codeDesilence_matches_c() {
    silence(true);
}

fn huffman(decode: bool) {
    unsafe {
        let mut rng = XorShift::new(25);
        for iter in 0..CASES {
            let n = size(&mut rng, iter);
            let data = byte_input(&mut rng, n, iter);
            let pos = if iter % 3 == 0 {
                rng.range_i32(0, n as i32) as usize
            } else {
                0
            };
            let capacity = if !decode && iter % 3 == 0 {
                rng.range_i32(0, 500) as usize
            } else {
                n * 2 + 2048
            };
            let outpos = if iter % 5 == 0 { capacity / 2 } else { 0 };
            let p = input_buffer(&data, pos);
            let cp = input_buffer(&vec![0xa5; capacity], outpos);
            let mut r = rust::WpcBuffer { code: data, pos };
            let mut ro = rust::WpcBuffer {
                code: vec![0xa5; capacity],
                pos: outpos,
            };
            let ca = huffCompress(p.cast(), cp.cast());
            let cb = rust::huffCompress(&mut r, &mut ro);
            assert_eq!(ca, cb, "case {iter}, n={n}");
            assert_eq!((*p).pos as usize, r.pos);
            assert_eq!((*cp).pos as usize, ro.pos);
            assert_eq!(bytes(cp), &ro.code[..ro.pos]);
            if capacity > 0 {
                assert_eq!(std::slice::from_raw_parts((*cp).code, capacity), ro.code);
            }
            if decode {
                assert!(ca >= 0);
                let encoded = bytes(cp)[outpos..].to_vec();
                let ci = input_buffer(&encoded, 0);
                let mut ri = rust::WpcBuffer {
                    code: encoded,
                    pos: 0,
                };
                let outlen = if iter % 3 == 0 { n / 2 } else { n + 8 };
                let outstart = if iter % 7 == 0 { outlen } else { 0 };
                let co = input_buffer(&vec![0xa5; outlen], outstart);
                let mut rr = rust::WpcBuffer {
                    code: vec![0xa5; outlen],
                    pos: outstart,
                };
                assert_eq!(
                    huffDecompress(ci.cast(), co.cast()),
                    rust::huffDecompress(&mut ri, &mut rr)
                );
                assert_eq!((*ci).pos as usize, ri.pos);
                assert_eq!((*co).pos as usize, rr.pos);
                assert_eq!(bytes(co), &rr.code[..rr.pos]);
                buffFree1(ci);
                buffFree1(co);
                // A truncated prefix of a valid stream exercises EOF in both
                // the count table and the bit payload, without inventing an
                // invalid tree that could make the C routine access garbage.
                let cut = rng.range_i32(0, ri.code.len() as i32 - 1) as usize;
                ri.code.truncate(cut);
                ri.pos = 0;
                let ci = input_buffer(&ri.code, 0);
                let co = input_buffer(&vec![0xa5; outlen], 0);
                let mut rr = rust::WpcBuffer {
                    code: vec![0xa5; outlen],
                    pos: 0,
                };
                assert_eq!(
                    huffDecompress(ci.cast(), co.cast()),
                    rust::huffDecompress(&mut ri, &mut rr)
                );
                assert_eq!((*ci).pos as usize, ri.pos);
                assert_eq!((*co).pos as usize, rr.pos);
                assert_eq!(bytes(co), &rr.code[..rr.pos]);
                buffFree1(ci);
                buffFree1(co);
            }
            buffFree1(p);
            buffFree1(cp);
        }
    }
}
#[test]
fn huffCompress_matches_c() {
    huffman(false);
}
#[test]
fn huffDecompress_matches_c() {
    huffman(true);
}

fn filter_cases() {
    unsafe {
        let mut rng = XorShift::new(26);
        for iter in 0..(2 * CASES) {
            let id = if iter < 16 {
                iter as i32
            } else {
                rng.range_i32(0, 25)
            };
            let p = waveGetfilter(id);
            let r = rust::waveGetfilter(id);
            assert_eq!(p.is_null(), r.is_none());
            if let Some(r) = r {
                assert_eq!(
                    ((*p).id, (*p).len as usize, (*p).ishift),
                    (r.id, r.len, r.ishift)
                );
                // Coefficients and their signed reversals must be bit-exact.
                let h = std::slice::from_raw_parts((*p).h, r.len);
                let g = std::slice::from_raw_parts((*p).g, r.len);
                for i in 0..r.len {
                    assert_eq!(h[i].to_bits(), r.filterh[i].to_bits());
                    assert_eq!(g[i].to_bits(), r.filterg[i].to_bits());
                }
                waveFreefilter(p);
                rust::waveFreefilter(r);
            }
        }
    }
}
#[test]
fn waveGetfilter_matches_c() {
    filter_cases();
}
#[test]
fn waveFreefilter_matches_c() {
    filter_cases();
}

fn wave_one(packet: bool) {
    unsafe {
        let mut rng = XorShift::new(27);
        for iter in 0..CASES {
            let id = (iter % 15) as i32;
            let p = waveGetfilter(id);
            let r = rust::waveGetfilter(id).unwrap();
            let pow = if iter % 500 == 0 {
                12
            } else {
                rng.range_i32(0, 9) as usize
            };
            let n = 1 << pow;
            let stage = rng.range_i32(0, pow as i32 + 1) as usize;
            let mut a = signal(&mut rng, n, iter);
            let mut b = signal(&mut rng, n, iter + 7);
            let mut ra = a.clone();
            let mut rb = b.clone();
            let typ = [0, 1, -1, 2][iter % 4];
            let ca = if packet {
                wavePack_1(
                    a.as_mut_ptr(),
                    b.as_mut_ptr(),
                    p,
                    pow as i32,
                    stage as i32,
                    typ,
                )
            } else {
                waveTrans_1(
                    a.as_mut_ptr(),
                    b.as_mut_ptr(),
                    p,
                    pow as i32,
                    stage as i32,
                    typ,
                )
            };
            let cb = if packet {
                rust::wavePack_1(&mut ra, &mut rb, &r, pow, stage, typ)
            } else {
                rust::waveTrans_1(&mut ra, &mut rb, &r, pow, stage, typ)
            };
            assert_eq!(ca, cb, "packet={packet} case {iter}");
            close(&ra, &a);
            close(&rb, &b);
            // A forward transform followed by an inverse also uses C as oracle.
            if stage <= pow && typ == 0 {
                if packet {
                    assert_eq!(
                        wavePack_1(
                            a.as_mut_ptr(),
                            b.as_mut_ptr(),
                            p,
                            pow as i32,
                            stage as i32,
                            1
                        ),
                        rust::wavePack_1(&mut ra, &mut rb, &r, pow, stage, 1)
                    );
                } else {
                    assert_eq!(
                        waveTrans_1(
                            a.as_mut_ptr(),
                            b.as_mut_ptr(),
                            p,
                            pow as i32,
                            stage as i32,
                            1
                        ),
                        rust::waveTrans_1(&mut ra, &mut rb, &r, pow, stage, 1)
                    );
                }
                close(&ra, &a);
                close(&rb, &b);
            }
            waveFreefilter(p);
        }
    }
}
#[test]
fn wavePack_1_matches_c() {
    wave_one(true);
}
#[test]
fn waveTrans_1_matches_c() {
    wave_one(false);
}

fn wave_two(packet: bool) {
    unsafe {
        let mut rng = XorShift::new(28);
        for iter in 0..CASES {
            let id1 = (iter % 15) as i32;
            let id2 = rng.range_i32(0, 14);
            let p1 = waveGetfilter(id1);
            let r1 = rust::waveGetfilter(id1).unwrap();
            let p2 = waveGetfilter(id2);
            let r2 = rust::waveGetfilter(id2).unwrap();
            let pow1 = if iter % 500 == 0 {
                10
            } else {
                rng.range_i32(0, 5) as usize
            };
            let pow2 = rng.range_i32(0, 5) as usize;
            let s1 = rng.range_i32(0, pow1 as i32 + 1) as usize;
            let s2 = rng.range_i32(0, pow2 as i32 + 1) as usize;
            let n1 = 1 << pow1;
            let n2 = 1 << pow2;
            let mut a: Vec<_> = (0..n2).map(|_| signal(&mut rng, n1, iter)).collect();
            let mut b: Vec<_> = (0..n2).map(|_| signal(&mut rng, n1, iter + 7)).collect();
            let mut ra = a.clone();
            let mut rb = b.clone();
            let mut ap = ptrs(&mut a);
            let mut bp = ptrs(&mut b);
            let typ = [0, 1, -1][iter % 3];
            let ca = if packet {
                wavePack_2(
                    ap.as_mut_ptr(),
                    bp.as_mut_ptr(),
                    p1,
                    p2,
                    pow1 as i32,
                    pow2 as i32,
                    s1 as i32,
                    s2 as i32,
                    typ,
                )
            } else {
                waveTrans_2(
                    ap.as_mut_ptr(),
                    bp.as_mut_ptr(),
                    p1,
                    pow1 as i32,
                    pow2 as i32,
                    s1 as i32,
                    s2 as i32,
                    typ,
                )
            };
            let cb = if packet {
                rust::wavePack_2(&mut ra, &mut rb, &r1, &r2, pow1, pow2, s1, s2, typ)
            } else {
                rust::waveTrans_2(&mut ra, &mut rb, &r1, pow1, pow2, s1, s2, typ)
            };
            assert_eq!(ca, cb, "packet={packet} case {iter}");
            matrix_close(&ra, &a);
            matrix_close(&rb, &b);
            if s1 <= pow1 && s2 <= pow2 && typ == 0 {
                let ca = if packet {
                    wavePack_2(
                        ap.as_mut_ptr(),
                        bp.as_mut_ptr(),
                        p1,
                        p2,
                        pow1 as i32,
                        pow2 as i32,
                        s1 as i32,
                        s2 as i32,
                        1,
                    )
                } else {
                    waveTrans_2(
                        ap.as_mut_ptr(),
                        bp.as_mut_ptr(),
                        p1,
                        pow1 as i32,
                        pow2 as i32,
                        s1 as i32,
                        s2 as i32,
                        1,
                    )
                };
                let cb = if packet {
                    rust::wavePack_2(&mut ra, &mut rb, &r1, &r2, pow1, pow2, s1, s2, 1)
                } else {
                    rust::waveTrans_2(&mut ra, &mut rb, &r1, pow1, pow2, s1, s2, 1)
                };
                assert_eq!(ca, cb);
                matrix_close(&ra, &a);
                matrix_close(&rb, &b);
            }
            waveFreefilter(p1);
            waveFreefilter(p2);
        }
    }
}
#[test]
fn wavePack_2_matches_c() {
    wave_two(true);
}
#[test]
fn waveTrans_2_matches_c() {
    wave_two(false);
}

#[test]
fn lpred_matches_c() {
    unsafe {
        let mut rng = XorShift::new(29);
        for iter in 0..CASES {
            let n1 = if iter % 100 == 0 {
                257
            } else {
                rng.range_i32(1, 64) as usize
            };
            let n2 = rng.range_i32(2, 8) as usize;
            let lblock = rng.range_i32(1, 17) as usize;
            let hlfilter = rng.range_i32(0, 4) as usize;
            let nslope = rng.range_i32(1, n2 as i32 - 1) as usize;
            let padded = n1.div_ceil(lblock) * lblock;
            // C deliberately reads complete final blocks from f even past n1.
            // Allocate and initialize those samples instead of comparing garbage.
            let mut f: Vec<_> = (0..n2).map(|_| signal(&mut rng, padded, iter)).collect();
            let mut a = vec![vec![12345.0; padded + 1]; n2];
            let mut b = a.clone();
            let mut fp = ptrs(&mut f);
            let mut ap = ptrs(&mut a);
            let mut cp = CLp {
                hlfilter: 7,
                lfilter: 15,
                nslope: 3,
                nblock: 9,
                filter: std::ptr::null_mut(),
            };
            let original_cp = cp;
            let mut rp = rust::LpFilter {
                hlfilter: 7,
                lfilter: 15,
                nslope: 3,
                nblock: 9,
                filter: vec![vec![vec![42.0]]],
            };
            let original_rp = rp.clone();
            lpred(
                fp.as_mut_ptr(),
                n1 as i32,
                n2 as i32,
                lblock as i32,
                hlfilter as i32,
                nslope as i32,
                &mut cp,
                ap.as_mut_ptr(),
            );
            rust::lpred(&f, n1, lblock, hlfilter, nslope, &mut rp, &mut b);
            for i in 0..n2 {
                close(&b[i][..padded], &a[i][..padded]);
            }
            assert_eq!(cp, original_cp);
            assert_eq!(rp, original_rp);
            for row in a.iter().chain(&b) {
                assert_eq!(row[padded], 12345.0);
            }
            // Original C loses its locally allocated filter; no pointer is
            // returned to free. All buffers owned by this test are still released.
        }
    }
}
