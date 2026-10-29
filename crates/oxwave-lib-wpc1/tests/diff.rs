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

//! Differential tests for comp/dwpt/1d/lib/{wbuffalloc,wpc1,wpc1coding,wpc1quant,wpc1trans}.c.
//! Credits: CWP: Tong Chen, 1994–1995.
use oxwave_cref::XorShift;
use oxwave_lib_wpc1 as rust;
use std::ffi::c_void;

// wpc1Compress discards realloc's returned pointer (wpc1.c, after encoding).
// Keep shrinking allocations in place in this TEST BINARY so the unmodified C
// routine has a valid pointer, including dead traces whose encoded size is zero.
// C permits realloc to preserve the allocation when shrinking; no codec operation
// is replaced. Growing allocations follow malloc/copy/free. This policy applies
// to Rust allocations too, so each allocation still has one ordinary libc owner.
extern "C" {
    fn malloc(n: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    #[cfg(target_os = "macos")]
    fn malloc_size(p: *const c_void) -> usize;
    #[cfg(not(target_os = "macos"))]
    fn malloc_usable_size(p: *mut c_void) -> usize;
}

#[export_name = "realloc"]
unsafe extern "C" fn reference_realloc(p: *mut c_void, n: usize) -> *mut c_void {
    if p.is_null() {
        return malloc(n.max(1));
    }
    #[cfg(target_os = "macos")]
    let old = malloc_size(p);
    #[cfg(not(target_os = "macos"))]
    let old = malloc_usable_size(p);
    if n <= old {
        return p;
    }
    let new = malloc(n);
    if !new.is_null() {
        std::ptr::copy_nonoverlapping(p as *const u8, new as *mut u8, old);
        free(p);
    }
    new
}

#[repr(C)]
struct CConfig {
    powt: i32,
    tileszt: i32,
    staget: i32,
    nblock: i32,
    lblock: i32,
    ishift: i32,
    len: i32,
    h: [f32; 12],
    g: [f32; 12],
}
#[repr(C)]
struct CQuant {
    step: f32,
    ave: f32,
    numsmp: i32,
    nbit: i32,
    qx: *mut i32,
}
#[repr(C)]
struct CBuffer {
    mbound: i32,
    pos: i32,
    code: *mut u8,
}
#[repr(C)]
struct CFilter {
    len: i32,
    ishift: i32,
    h: *mut f32,
    g: *mut f32,
}
#[repr(C)]
struct CCompressed {
    numsmp: i32,
    nsize: i32,
    step: f32,
    ave: f32,
    code: *mut u8,
}

#[link(name = "wpc1", kind = "static")]
extern "C" {
    fn buffAlloc1(n: i32) -> *mut CBuffer;
    fn buffRealloc1(b: *mut CBuffer, n: i32);
    fn buffFree1(b: *mut CBuffer);
    fn buffAlloc2(n1: i32, n2: i32) -> *mut *mut CBuffer;
    fn buffFree2(b: *mut *mut CBuffer, n1: i32);
    fn wpc1Quant(x: *mut f32, n: i32, error: f32, c: *const c_void, q: *mut c_void);
    fn wpc1Dequant(x: *mut f32, n: i32, c: *const c_void, q: *mut c_void);
    fn wavePack1(x: *mut f32, y: *mut f32, f: *mut CFilter, pow: i32, stage: i32, typ: i32);
    fn wpc1Trans(x: *mut f32, c: *const c_void, typ: i32);
    fn bookCoder(buff: *mut c_void, book: *mut c_void, q: *mut c_void, n: i32);
    fn bookDecoder(buff: *mut c_void, book: *mut c_void, q: *mut c_void, n: i32) -> i32;
    fn entCoder(q: *mut c_void, n: i32, buff: *mut c_void);
    fn entDecoder(q: *mut c_void, n: i32, buff: *mut c_void) -> i32;
    fn wpc1Compress(f: *mut f32, n: i32, error: f32, bytes: *mut i32) -> *mut c_void;
    fn wpc1Uncompress(w: *mut c_void, f: *mut f32) -> i32;
    fn wpc1Numsmp(w: *mut c_void) -> i32;
    fn wpc1Size(w: *mut c_void) -> i32;
    fn wpc1Free(w: *mut c_void);
    fn wpc1Write(w: *mut c_void, fp: *mut c_void);
    fn wpc1Read(fp: *mut c_void) -> *mut c_void;
    fn wpc1Encoder(q: *mut c_void, n: i32, w: *mut c_void);
    fn wpc1Decoder(q: *mut c_void, n: i32, w: *mut c_void) -> i32;
}

extern "C" {
    fn tmpfile() -> *mut c_void;
    fn fclose(f: *mut c_void) -> i32;
    fn rewind(f: *mut c_void);
    fn fwrite(p: *const c_void, size: usize, n: usize, f: *mut c_void) -> usize;
    fn fread(p: *mut c_void, size: usize, n: usize, f: *mut c_void) -> usize;
}

// The original compression functions use static configuration pointers.
static C_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
const SIZES: [usize; 14] = [1, 2, 3, 17, 31, 64, 127, 128, 129, 251, 256, 257, 509, 4093];

fn bits_eq(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    for (i, (&a, &b)) in a.iter().zip(b).enumerate() {
        assert_eq!(a.to_bits(), b.to_bits(), "sample {i}: {a} versus {b}");
    }
}

fn config(pow: i32, stage: i32) -> (CConfig, rust::Config) {
    let h = [
        0.016387336464f32,
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
    (
        CConfig {
            powt: pow,
            tileszt: 1 << pow,
            staget: stage,
            nblock: 1 << stage,
            lblock: 1 << (pow - stage),
            ishift: -4,
            len: 12,
            h,
            g,
        },
        rust::Config {
            powt: pow,
            tileszt: 1 << pow,
            staget: stage,
            nblock: 1 << stage,
            lblock: 1 << (pow - stage),
            ishift: -4,
            len: 12,
            filterh: h,
            filterg: g,
        },
    )
}

fn samples(rng: &mut XorShift, n: usize, trial: usize) -> Vec<f32> {
    let scale = [1e-12, 0.01, 1.0, 1000.0, 1e12][trial % 5];
    let mut x: Vec<f32> = (0..n).map(|_| rng.uniform_f32(-scale, scale)).collect();
    match trial % 19 {
        0 => x.fill(0.0),
        1 => x.fill(-scale),
        2 => {
            x.fill(0.0);
            x[n / 2] = -scale;
        }
        3 => {
            for (i, v) in x.iter_mut().enumerate() {
                *v = if i % 2 == 0 { scale } else { -scale };
            }
        }
        _ => {}
    }
    x
}

#[test]
fn buff_alloc1_matches_c() {
    let mut rng = XorShift::new(11);
    unsafe {
        for trial in 0..3000 {
            let n = if trial % 11 == 0 {
                0
            } else {
                SIZES[trial % SIZES.len()] + (rng.next_u32() % 2) as usize
            };
            let c = buffAlloc1(n as i32);
            let r = rust::buffAlloc1(n);
            assert_eq!((*c).mbound as usize, r.code.len());
            assert_eq!((*c).pos as usize, r.pos);
            buffFree1(c);
            rust::buffFree1(r);
        }
    }
}

#[test]
fn buff_realloc1_matches_c() {
    let mut rng = XorShift::new(12);
    unsafe {
        for trial in 0..3000 {
            let n = SIZES[trial % SIZES.len()];
            let new = (rng.next_u32() as usize) % (n * 2 + 1);
            let c = buffAlloc1(n as i32);
            let mut r = rust::buffAlloc1(n);
            for i in 0..n {
                let byte = rng.next_u32() as u8;
                *(*c).code.add(i) = byte;
                r.code[i] = byte;
            }
            let pos = rng.next_u32() as usize % (n + 1);
            (*c).pos = pos as i32;
            r.pos = pos;
            buffRealloc1(c, new as i32);
            rust::buffRealloc1(&mut r, new);
            assert_eq!((*c).pos as usize, r.pos);
            assert_eq!((*c).mbound as usize, r.code.len());
            assert_eq!(
                std::slice::from_raw_parts((*c).code, new.min(n)),
                &r.code[..new.min(n)]
            );
            buffFree1(c);
            rust::buffFree1(r);
        }
    }
}

#[test]
fn buff_free1_matches_c() {
    let mut rng = XorShift::new(13);
    unsafe {
        for trial in 0..3000 {
            let n = SIZES[trial % SIZES.len()];
            let c = buffAlloc1(n as i32);
            let mut r = rust::buffAlloc1(n);
            for i in 0..n {
                let byte = rng.next_u32() as u8;
                *(*c).code.add(i) = byte;
                r.code[i] = byte;
            }
            assert_eq!(std::slice::from_raw_parts((*c).code, n), r.code);
            buffFree1(c);
            rust::buffFree1(r);
        }
    }
}

fn check_buffers2(seed: u64) {
    let mut rng = XorShift::new(seed);
    unsafe {
        for trial in 0..2500 {
            let n1 = rng.next_u32() as usize % 33;
            let n2 = if trial % 17 == 0 {
                0
            } else {
                SIZES[trial % SIZES.len()]
            };
            let c = buffAlloc2(n1 as i32, n2 as i32);
            let mut r = rust::buffAlloc2(n1, n2);
            for (i, r) in r.iter_mut().enumerate() {
                let b = *c.add(i);
                assert_eq!((*b).mbound as usize, r.code.len());
                assert_eq!((*b).pos as usize, r.pos);
                for j in 0..n2 {
                    let byte = rng.next_u32() as u8;
                    *(*b).code.add(j) = byte;
                    r.code[j] = byte;
                }
                assert_eq!(std::slice::from_raw_parts((*b).code, n2), r.code);
            }
            buffFree2(c, n1 as i32);
            rust::buffFree2(r);
        }
    }
}

#[test]
fn buff_alloc2_matches_c() {
    check_buffers2(14);
}
#[test]
fn buff_free2_matches_c() {
    check_buffers2(15);
}

#[test]
fn wpc1_quant_matches_c() {
    let mut rng = XorShift::new(21);
    for trial in 0..4000 {
        let (cc, rc) = config(8, 4);
        let transformed = trial % 3 == 0;
        let n = if transformed {
            256
        } else {
            SIZES[trial % SIZES.len()]
        };
        let cn = if transformed { 0 } else { n as i32 };
        let mut x = samples(&mut rng, n, trial);
        let mut y = x.clone();
        let mut cq = vec![0; n];
        let mut q = CQuant {
            step: 0.0,
            ave: 0.0,
            numsmp: 0,
            nbit: 0,
            qx: cq.as_mut_ptr(),
        };
        let mut rq = rust::Quant::default();
        let error = [0.0001, 0.01, 0.05, 0.1, 0.5, 1.0, 10.0][trial % 7];
        unsafe {
            wpc1Quant(
                x.as_mut_ptr(),
                cn,
                error,
                &cc as *const _ as _,
                &mut q as *mut _ as _,
            );
        }
        rust::wpc1Quant(&mut y, cn as usize, error, &rc, &mut rq);
        bits_eq(&x, &y);
        assert_eq!(q.nbit, rq.nbit);
        assert_eq!(q.numsmp as usize, rq.numsmp);
        assert_eq!(q.step.to_bits(), rq.step.to_bits());
        assert_eq!(q.ave.to_bits(), rq.ave.to_bits());
        assert_eq!(cq, rq.qx);
    }
}

#[test]
fn wpc1_dequant_matches_c() {
    let mut rng = XorShift::new(22);
    for trial in 0..4000 {
        let (cc, rc) = config(8, 4);
        let transformed = trial % 3 == 0;
        let n = if transformed {
            256
        } else {
            SIZES[trial % SIZES.len()]
        };
        let cn = if transformed { 0 } else { n as i32 };
        let mut qx: Vec<i32> = (0..n).map(|_| rng.range_i32(-32767, 32767)).collect();
        let step = match trial % 9 {
            0 => 0.0,
            1 => f32::MIN_POSITIVE,
            2 => -1.0,
            _ => rng.uniform_f32(0.001, 1000.0),
        };
        let ave = rng.uniform_f32(-1e5, 1e5);
        let mut q = CQuant {
            step,
            ave,
            numsmp: n as i32,
            nbit: 14,
            qx: qx.as_mut_ptr(),
        };
        let rq = rust::Quant {
            step,
            ave,
            numsmp: n,
            nbit: 14,
            qx,
        };
        let mut x = vec![123.0; n];
        let mut y = x.clone();
        unsafe {
            wpc1Dequant(
                x.as_mut_ptr(),
                cn,
                &cc as *const _ as _,
                &mut q as *mut _ as _,
            );
        }
        rust::wpc1Dequant(&mut y, cn as usize, &rc, &rq);
        bits_eq(&x, &y);
    }
}

fn check_transform(seed: u64, wrapper: bool) {
    let mut rng = XorShift::new(seed);
    for trial in 0..3000 {
        let pow = (trial % 11) as i32;
        let n = 1usize << pow;
        let stage = (rng.next_u32() % (pow as u32 + 1)) as i32;
        let typ = ((trial / 11) % 2) as i32;
        let (mut cc, mut rc) = config(pow, stage);
        if trial % 3 != 0 {
            cc.len = [2, 6, 12][trial % 3];
            rc.len = cc.len as usize;
            cc.ishift = 2 * rng.range_i32(-16, 16);
            rc.ishift = cc.ishift;
            for i in 0..rc.len {
                cc.h[i] = rng.uniform_f32(-0.5, 0.5);
                cc.g[i] = rng.uniform_f32(-0.5, 0.5);
            }
            rc.filterh = cc.h;
            rc.filterg = cc.g;
        }
        let mut cx = samples(&mut rng, n, trial);
        let mut cy = samples(&mut rng, n, trial + 1);
        let mut rx = cx.clone();
        let mut ry = cy.clone();
        if wrapper {
            unsafe {
                wpc1Trans(cx.as_mut_ptr(), &cc as *const _ as _, typ);
            }
            rust::wpc1Trans(&mut rx, &rc, typ);
        } else {
            let mut cf = CFilter {
                len: cc.len,
                ishift: cc.ishift,
                h: cc.h.as_mut_ptr(),
                g: cc.g.as_mut_ptr(),
            };
            let rf = rust::WaveFilter {
                ishift: rc.ishift,
                h: rc.filterh[..rc.len].to_vec(),
                g: rc.filterg[..rc.len].to_vec(),
            };
            unsafe {
                wavePack1(cx.as_mut_ptr(), cy.as_mut_ptr(), &mut cf, pow, stage, typ);
            }
            rust::wavePack1(&mut rx, &mut ry, &rf, pow, stage, typ);
            bits_eq(&cy, &ry);
        }
        bits_eq(&cx, &rx);
    }
}

#[test]
fn wave_pack1_matches_c() {
    check_transform(31, false);
}
#[test]
fn wpc1_trans_matches_c() {
    check_transform(32, true);
}

fn check_io(seed: u64) {
    let mut rng = XorShift::new(seed);
    for trial in 0..2000 {
        let n = if trial % 17 == 0 {
            0
        } else {
            SIZES[trial % SIZES.len()]
        };
        let mut code: Vec<u8> = (0..n).map(|_| rng.next_u32() as u8).collect();
        let w = rust::Compressed {
            numsmp: rng.next_u32() as i32,
            overallstep: rng.uniform_f32(-1e12, 1e12),
            overallave: rng.uniform_f32(-1e12, 1e12),
            code: code.clone(),
        };
        let mut c = CCompressed {
            numsmp: w.numsmp,
            nsize: n as i32,
            step: w.overallstep,
            ave: w.overallave,
            code: code.as_mut_ptr(),
        };
        unsafe {
            let f = tmpfile();
            assert!(!f.is_null());
            wpc1Write(&mut c as *mut _ as _, f);
            rewind(f);
            let mut bytes = vec![0; n + 16];
            assert_eq!(
                fread(bytes.as_mut_ptr() as _, 1, bytes.len(), f),
                bytes.len()
            );
            let mut actual = Vec::new();
            rust::wpc1Write(&w, &mut actual).unwrap();
            assert_eq!(bytes, actual);
            rewind(f);
            let cw = wpc1Read(f);
            assert!(!cw.is_null());
            assert_eq!(wpc1Numsmp(cw), rust::wpc1Numsmp(&w));
            assert_eq!(wpc1Size(cw), rust::wpc1Size(&w));
            let rw = rust::wpc1Read(&mut bytes.as_slice()).unwrap();
            assert_eq!(rw, w);
            let cd = &*(cw as *const CCompressed);
            assert_eq!(cd.step.to_bits(), rw.overallstep.to_bits());
            assert_eq!(cd.ave.to_bits(), rw.overallave.to_bits());
            assert_eq!(std::slice::from_raw_parts(cd.code, n), rw.code);
            wpc1Free(cw);
            rust::wpc1Free(rw);
            assert_eq!(fclose(f), 0);
        }
    }
}

#[test]
fn wpc1_write_matches_c() {
    check_io(41);
}
#[test]
fn wpc1_read_matches_c() {
    check_io(42);
    // C returns NULL on each of the four incomplete header fields.
    for len in 0..16 {
        let bytes = vec![0u8; len];
        unsafe {
            let f = tmpfile();
            assert!(!f.is_null());
            assert_eq!(fwrite(bytes.as_ptr() as _, 1, len, f), len);
            rewind(f);
            assert!(wpc1Read(f).is_null());
            assert!(rust::wpc1Read(&mut bytes.as_slice()).is_err());
            fclose(f);
        }
    }
}
#[test]
fn wpc1_numsmp_matches_c() {
    check_io(43);
}
#[test]
fn wpc1_size_matches_c() {
    check_io(44);
}
#[test]
fn wpc1_free_matches_c() {
    check_io(45);
}

#[repr(C)]
struct CEntCode {
    value: u32,
    bits: i32,
}
#[repr(C)]
struct CBook {
    maxrun: i32,
    bitmaxrun: i32,
    booksize: i32,
    code: *mut CEntCode,
    pref2bit: *mut i32,
}
#[repr(C)]
struct CBitBuffer {
    buff: *mut CBuffer,
    mask: u8,
    rack: u8,
}

// Own all C arrays, so every raw pointer remains valid until the call finishes.
struct CQuants {
    data: Vec<Vec<i32>>,
    records: Vec<CQuant>,
    pointers: Vec<*mut c_void>,
}

impl CQuants {
    fn new(q: &[rust::Quant]) -> Self {
        let mut data: Vec<Vec<i32>> = q.iter().map(|q| q.qx.clone()).collect();
        let mut records: Vec<CQuant> = q
            .iter()
            .zip(&mut data)
            .map(|(q, d)| CQuant {
                step: q.step,
                ave: q.ave,
                numsmp: q.numsmp as i32,
                nbit: q.nbit,
                qx: d.as_mut_ptr(),
            })
            .collect();
        let pointers = records
            .iter_mut()
            .map(|q| q as *mut _ as *mut c_void)
            .collect();
        Self {
            data,
            records,
            pointers,
        }
    }
    fn ptr(&mut self) -> *mut c_void {
        self.pointers.as_mut_ptr() as _
    }
    fn assert_eq(&self, q: &[rust::Quant]) {
        for ((c, d), r) in self.records.iter().zip(&self.data).zip(q) {
            assert_eq!(c.step.to_bits(), r.step.to_bits());
            assert_eq!(c.ave.to_bits(), r.ave.to_bits());
            assert_eq!(c.nbit, r.nbit);
            assert_eq!(c.numsmp as usize, r.numsmp);
            assert_eq!(d, &r.qx);
        }
    }
}

fn symbols(rng: &mut XorShift, n: usize, nbit: usize, trial: usize) -> Vec<i32> {
    // For books >= 3 the all-one amplitude is reserved as EOB. This is why
    // wpc1Quant reduces its scale at the upper boundary.
    let max = (1 << (nbit + 1)) - 1 - if nbit >= 3 { 1 } else { 0 };
    (0..n)
        .map(|i| match trial % 7 {
            0 => 0,
            1 => {
                if i % 2 == 0 {
                    max
                } else {
                    -max
                }
            }
            2 => {
                if i % 31 == 0 {
                    -1
                } else {
                    0
                }
            }
            _ => rng.range_i32(-max, max),
        })
        .collect()
}

fn quant_blocks(rng: &mut XorShift, trial: usize) -> Vec<rust::Quant> {
    let nblock = 1 + rng.next_u32() as usize % 8;
    (0..nblock)
        .map(|i| {
            let n = SIZES[(trial + i) % SIZES.len()];
            let nbit = (trial + i) % 15;
            let mut step = rng.uniform_f32(0.001, 1e4);
            let mut ave = rng.uniform_f32(-1e4, 1e4);
            if trial % 11 < 2 {
                step = f32::MIN_POSITIVE;
            }
            if trial % 11 == 0 || trial % 11 == 2 {
                ave = 0.0;
            }
            rust::Quant {
                step,
                ave,
                numsmp: n,
                nbit: nbit as i32,
                qx: symbols(rng, n, nbit, trial + i),
            }
        })
        .collect()
}

fn empty_quants(source: &[rust::Quant]) -> Vec<rust::Quant> {
    source
        .iter()
        .map(|q| rust::Quant {
            step: -123.0,
            ave: 0.0,
            numsmp: q.numsmp,
            nbit: -1,
            qx: vec![-999; q.numsmp],
        })
        .collect()
}

#[test]
fn book_coder_matches_c() {
    let books = test_books();
    let mut rng = XorShift::new(51);
    for trial in 0..3000 {
        let nbit = trial % 15;
        let book = &books[nbit];
        let n = SIZES[trial % SIZES.len()];
        let mut q = symbols(&mut rng, n, nbit, trial);
        let mut codes: Vec<CEntCode> = book
            .code
            .iter()
            .map(|c| CEntCode {
                value: c.value,
                bits: c.bits,
            })
            .collect();
        let mut pref = book.pref2bit.clone();
        let mut cb = CBook {
            maxrun: book.maxrun,
            bitmaxrun: book.bitmaxrun,
            booksize: codes.len() as i32,
            code: codes.as_mut_ptr(),
            pref2bit: pref.as_mut_ptr(),
        };
        let capacity = if trial % 13 == 0 {
            trial % 7
        } else {
            n * 4 + 32
        };
        let mut bytes = vec![0xa5; capacity];
        let pos = if capacity == 0 {
            0
        } else {
            rng.next_u32() as usize % capacity.min(4)
        };
        let mut cbuf = CBuffer {
            mbound: capacity as i32,
            pos: pos as i32,
            code: bytes.as_mut_ptr(),
        };
        let mask = 0x80u8 >> (trial % 8);
        let rack = rng.next_u32() as u8 & !((mask << 1).wrapping_sub(1));
        let mut cbit = CBitBuffer {
            buff: &mut cbuf,
            mask,
            rack,
        };
        let mut rbuf = rust::Buffer {
            code: bytes.clone(),
            pos,
        };
        let mut rbit = rust::BitBuffer {
            buff: &mut rbuf,
            mask,
            rack,
        };
        unsafe {
            bookCoder(
                &mut cbit as *mut _ as _,
                &mut cb as *mut _ as _,
                q.as_mut_ptr() as _,
                n as i32,
            );
        }
        rust::bookCoder(&mut rbit, book, &q);
        assert_eq!((cbit.mask, cbit.rack), (rbit.mask, rbit.rack));
        assert_eq!(cbuf.pos as usize, rbuf.pos);
        assert_eq!(bytes, rbuf.code);
    }
}

#[test]
fn book_decoder_matches_c() {
    let books = test_books();
    let mut rng = XorShift::new(52);
    for trial in 0..3000 {
        let nbit = trial % 15;
        let book = &books[nbit];
        let n = SIZES[trial % SIZES.len()];
        let mut q = symbols(&mut rng, n, nbit, trial);
        let mut codes: Vec<CEntCode> = book
            .code
            .iter()
            .map(|c| CEntCode {
                value: c.value,
                bits: c.bits,
            })
            .collect();
        let mut pref = book.pref2bit.clone();
        let mut cb = CBook {
            maxrun: book.maxrun,
            bitmaxrun: book.bitmaxrun,
            booksize: codes.len() as i32,
            code: codes.as_mut_ptr(),
            pref2bit: pref.as_mut_ptr(),
        };
        let mut bytes = vec![0; n * 4 + 32];
        let mut cbuf = CBuffer {
            mbound: bytes.len() as i32,
            pos: 0,
            code: bytes.as_mut_ptr(),
        };
        let mut cbit = CBitBuffer {
            buff: &mut cbuf,
            mask: 0x80,
            rack: 0,
        };
        unsafe {
            bookCoder(
                &mut cbit as *mut _ as _,
                &mut cb as *mut _ as _,
                q.as_mut_ptr() as _,
                n as i32,
            );
        }
        if cbit.mask != 0x80 {
            bytes[cbuf.pos as usize] = cbit.rack;
            cbuf.pos += 1;
        }
        let size = cbuf.pos as usize;
        let size = if trial % 5 == 0 { size / 2 } else { size };
        bytes.truncate(size);
        cbuf.mbound = size as i32;
        cbuf.pos = 0;
        cbit.mask = 0x80;
        cbit.rack = 0;
        let expected_n = n + usize::from(trial % 7 == 0);
        // The C routine checks count only after EOB. Supply guard elements.
        let mut cq = vec![-999; expected_n + 32];
        let mut rq = vec![-999; expected_n];
        let mut rb = rust::Buffer {
            code: bytes,
            pos: 0,
        };
        let mut rbit = rust::BitBuffer {
            buff: &mut rb,
            mask: 0x80,
            rack: 0,
        };
        let ret = unsafe {
            bookDecoder(
                &mut cbit as *mut _ as _,
                &mut cb as *mut _ as _,
                cq.as_mut_ptr() as _,
                expected_n as i32,
            )
        };
        assert_eq!(ret, rust::bookDecoder(&mut rbit, book, &mut rq));
        assert_eq!(&cq[..expected_n], rq);
        assert_eq!((cbit.mask, cbit.rack), (rbit.mask, rbit.rack));
        assert_eq!(cbuf.pos as usize, rb.pos);
        if ret == 1 {
            assert_eq!(rq, q);
        }
    }
}

#[test]
fn ent_coder_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(53);
    for trial in 0..2500 {
        let q = quant_blocks(&mut rng, trial);
        let mut cq = CQuants::new(&q);
        let capacity = if trial % 13 == 0 {
            trial % 17
        } else {
            q.iter().map(|q| q.numsmp * 4).sum::<usize>() + 128
        };
        let mut bytes = vec![0xa5; capacity];
        let pos = if capacity == 0 {
            0
        } else {
            trial % capacity.min(4)
        };
        let mut cb = CBuffer {
            mbound: capacity as i32,
            pos: pos as i32,
            code: bytes.as_mut_ptr(),
        };
        let mut rb = rust::Buffer {
            code: bytes.clone(),
            pos,
        };
        unsafe {
            entCoder(cq.ptr(), q.len() as i32, &mut cb as *mut _ as _);
        }
        rust::entCoder(&q, &mut rb);
        assert_eq!(cb.pos as usize, rb.pos);
        assert_eq!(bytes, rb.code);
        cq.assert_eq(&q);
    }
}

#[test]
fn ent_decoder_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(54);
    for trial in 0..2500 {
        let q = quant_blocks(&mut rng, trial);
        let mut source = CQuants::new(&q);
        let mut bytes = vec![0; q.iter().map(|q| q.numsmp * 4).sum::<usize>() + 128];
        let mut cb = CBuffer {
            mbound: bytes.len() as i32,
            pos: 0,
            code: bytes.as_mut_ptr(),
        };
        unsafe {
            entCoder(source.ptr(), q.len() as i32, &mut cb as *mut _ as _);
        }
        bytes.truncate(cb.pos as usize);
        if trial % 5 == 0 {
            bytes.truncate(bytes.len() / 2);
        }
        cb.mbound = bytes.len() as i32;
        cb.pos = 0;
        let mut rq = empty_quants(&q);
        let mut cq = CQuants::new(&rq);
        let mut rb = rust::Buffer {
            code: bytes,
            pos: 0,
        };
        let ret = unsafe { entDecoder(cq.ptr(), q.len() as i32, &mut cb as *mut _ as _) };
        assert_eq!(ret, rust::entDecoder(&mut rq, &mut rb));
        assert_eq!(cb.pos as usize, rb.pos);
        cq.assert_eq(&rq);
        if ret == 1 {
            for (q, r) in q.iter().zip(&rq) {
                assert_eq!(q.qx, r.qx);
            }
        }
    }
}

fn ccompressed(w: &mut rust::Compressed) -> CCompressed {
    CCompressed {
        numsmp: w.numsmp,
        nsize: w.code.len() as i32,
        step: w.overallstep,
        ave: w.overallave,
        code: w.code.as_mut_ptr(),
    }
}

#[test]
fn wpc1_encoder_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(55);
    for trial in 0..2500 {
        let q = quant_blocks(&mut rng, trial);
        let mut cq = CQuants::new(&q);
        let capacity = if trial % 13 == 0 {
            trial % 17
        } else {
            q.iter().map(|q| q.numsmp * 4).sum::<usize>() + 128
        };
        let mut r = rust::Compressed {
            numsmp: q.iter().map(|q| q.numsmp as i32).sum(),
            overallstep: -1.0,
            overallave: -1.0,
            code: vec![0xa5; capacity],
        };
        let mut storage = r.clone();
        let mut c = ccompressed(&mut storage);
        unsafe {
            wpc1Encoder(cq.ptr(), q.len() as i32, &mut c as *mut _ as _);
        }
        rust::wpc1Encoder(&q, &mut r);
        assert_eq!(c.nsize as usize, r.code.len());
        assert_eq!(c.step.to_bits(), r.overallstep.to_bits());
        assert_eq!(c.ave.to_bits(), r.overallave.to_bits());
        assert_eq!(&storage.code[..c.nsize as usize], r.code);
        cq.assert_eq(&q);
    }
}

#[test]
fn wpc1_decoder_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(56);
    for trial in 0..2500 {
        let q = quant_blocks(&mut rng, trial);
        let mut source = CQuants::new(&q);
        let mut w = rust::Compressed {
            numsmp: q.iter().map(|q| q.numsmp as i32).sum(),
            overallstep: 0.0,
            overallave: 0.0,
            code: vec![0; q.iter().map(|q| q.numsmp * 4).sum::<usize>() + 128],
        };
        let mut c = ccompressed(&mut w);
        unsafe {
            wpc1Encoder(source.ptr(), q.len() as i32, &mut c as *mut _ as _);
        }
        w.overallstep = c.step;
        w.overallave = c.ave;
        w.code.truncate(c.nsize as usize);
        if trial % 5 == 0 {
            w.code.truncate(w.code.len() / 2);
        }
        c.nsize = w.code.len() as i32;
        let mut rq = empty_quants(&q);
        // Test the caller's existing averages are preserved on the C zero-mean path.
        if trial % 11 == 0 || trial % 11 == 2 {
            for q in &mut rq {
                q.ave = -17.0;
            }
        }
        let mut cq = CQuants::new(&rq);
        let ret = unsafe { wpc1Decoder(cq.ptr(), q.len() as i32, &mut c as *mut _ as _) };
        assert_eq!(ret, rust::wpc1Decoder(&mut rq, &w));
        cq.assert_eq(&rq);
    }
}

const TRACE_SIZES: [usize; 22] = [
    1, 2, 3, 17, 31, 64, 127, 128, 129, 251, 255, 256, 257, 383, 384, 385, 509, 512, 513, 1024,
    4093, 8192,
];

#[test]
fn wpc1_compress_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(61);
    for trial in 0..3000 {
        let n = TRACE_SIZES[trial % TRACE_SIZES.len()];
        let mut x = samples(&mut rng, n, trial);
        let mut y = x.clone();
        let error = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0][trial % 6];
        let mut cn = -1;
        let mut rn = -1;
        unsafe {
            let p = wpc1Compress(x.as_mut_ptr(), n as i32, error, &mut cn);
            assert!(!p.is_null());
            let c = &*(p as *const CCompressed);
            let r = rust::wpc1Compress(&mut y, error, &mut rn);
            assert_eq!(cn, rn, "trial {trial}, n={n}");
            assert_eq!(c.numsmp, r.numsmp);
            assert_eq!(c.step.to_bits(), r.overallstep.to_bits());
            assert_eq!(c.ave.to_bits(), r.overallave.to_bits());
            assert_eq!(
                std::slice::from_raw_parts(c.code, cn as usize),
                r.code,
                "trial {trial}, n={n}"
            );
            bits_eq(&x, &y);
            wpc1Free(p);
            rust::wpc1Free(r);
        }
    }
}

// A deterministic C oracle for the original zero-mean path: the C top-level
// uncompress allocates quant structs without initializing ave, and its decoder
// leaves ave untouched when overallave<=FLT_MIN. Compose the exact same original
// C routines with ave initialized to zero. This also verifies dead traces without
// treating undefined malloc contents as reference samples.
fn c_uncompress_initialized(w: &rust::Compressed) -> (i32, Vec<f32>) {
    let (cc, _) = config(8, 4);
    let n = w.numsmp as usize;
    let left = n % 256;
    let regular = n / 256 + usize::from(left > 128);
    let mut q = vec![
        rust::Quant {
            numsmp: 256,
            qx: vec![0; 256],
            ..Default::default()
        };
        regular
    ];
    if left > 0 && left <= 128 {
        q.push(rust::Quant {
            numsmp: left,
            qx: vec![0; left],
            ..Default::default()
        });
    }
    let mut cq = CQuants::new(&q);
    let mut storage = w.clone();
    let mut cw = ccompressed(&mut storage);
    let mut f = vec![0.0; (regular * 256).max(n)];
    unsafe {
        let ret = wpc1Decoder(cq.ptr(), q.len() as i32, &mut cw as *mut _ as _);
        if ret != 0 {
            for ip in 0..regular {
                let block = f.as_mut_ptr().add(ip * 256);
                wpc1Dequant(block, 0, &cc as *const _ as _, cq.pointers[ip]);
                wpc1Trans(block, &cc as *const _ as _, 1);
            }
            if left > 0 && left <= 128 {
                wpc1Dequant(
                    f.as_mut_ptr().add(regular * 256),
                    left as i32,
                    &cc as *const _ as _,
                    cq.pointers[regular],
                );
            }
        }
        f.truncate(n);
        (ret, f)
    }
}

#[test]
fn wpc1_uncompress_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(62);
    for trial in 0..3000 {
        let n = TRACE_SIZES[trial % TRACE_SIZES.len()];
        let mut x = samples(&mut rng, n, trial);
        // Nonzero averages make the direct C API's quantization state defined.
        // Keep this separate from the all-zero cases below, which use C composition.
        let dc = [0.1f32, -1.0, 123.0][trial % 3];
        for v in &mut x {
            *v = *v * 1e-12 + dc;
        }
        let mut cn = 0;
        unsafe {
            let p = wpc1Compress(x.as_mut_ptr(), n as i32, 0.05, &mut cn);
            assert!(!p.is_null());
            let c = &*(p as *const CCompressed);
            assert!(c.ave > f32::MIN_POSITIVE);
            let mut w = rust::Compressed {
                numsmp: c.numsmp,
                overallstep: c.step,
                overallave: c.ave,
                code: std::slice::from_raw_parts(c.code, cn as usize).to_vec(),
            };
            let mut cx = vec![0.0; n];
            let mut rx = cx.clone();
            let cr = wpc1Uncompress(p, cx.as_mut_ptr());
            let rr = rust::wpc1Uncompress(&w, &mut rx);
            assert_eq!(cr, rr, "trial {trial}, n={n}");
            assert_eq!(cr, 1);
            bits_eq(&cx, &rx);
            wpc1Free(p);

            // Zero trace and a symmetric zero-mean trace, both evaluated on C's
            // initialized internal decoder/dequant/transform path.
            let mut zero = vec![0.0; n];
            let wz = rust::wpc1Compress(&mut zero, 0.05, &mut cn);
            let (ret, expected) = c_uncompress_initialized(&wz);
            assert_eq!(rust::wpc1Uncompress(&wz, &mut rx), ret);
            bits_eq(&expected, &rx);

            let mut symmetric = vec![0.0; n];
            for pair in symmetric.chunks_exact_mut(2) {
                let v = rng.uniform_f32(0.001, 10.0);
                pair[0] = v;
                pair[1] = -v;
            }
            let ws = rust::wpc1Compress(&mut symmetric, 0.05, &mut cn);
            let (ret, expected) = c_uncompress_initialized(&ws);
            assert_eq!(rust::wpc1Uncompress(&ws, &mut rx), ret);
            bits_eq(&expected, &rx);

            // Invalid sample count fails the original's initial consistency check.
            w.numsmp = -(trial as i32);
            let mut cw = ccompressed(&mut w);
            cx.fill(17.0);
            rx.fill(17.0);
            assert_eq!(
                wpc1Uncompress(&mut cw as *mut _ as _, cx.as_mut_ptr()),
                rust::wpc1Uncompress(&w, &mut rx)
            );
            bits_eq(&cx, &rx);
        }
    }
}

// Codebook fixtures transcribed from the original C builddict().
fn test_books() -> Vec<rust::CodeBook> {
    vec![
        rust::CodeBook {
            maxrun: 16,
            bitmaxrun: 4,
            code: vec![
                rust::EntCode { value: 0, bits: 1 },
                rust::EntCode { value: 2, bits: 2 },
                rust::EntCode { value: 3, bits: 2 },
            ],
            pref2bit: vec![0, 0, 1, 0],
        },
        rust::CodeBook {
            maxrun: 16,
            bitmaxrun: 4,
            code: vec![
                rust::EntCode { value: 0, bits: 2 },
                rust::EntCode { value: 1, bits: 2 },
                rust::EntCode { value: 2, bits: 2 },
                rust::EntCode { value: 3, bits: 2 },
            ],
            pref2bit: vec![0, 1, 2, 0],
        },
        rust::CodeBook {
            maxrun: 16,
            bitmaxrun: 4,
            code: vec![
                rust::EntCode { value: 0, bits: 2 },
                rust::EntCode { value: 1, bits: 2 },
                rust::EntCode { value: 2, bits: 2 },
                rust::EntCode { value: 6, bits: 3 },
                rust::EntCode { value: 7, bits: 3 },
            ],
            pref2bit: vec![0, 1, 2, 0, 0, 0, 3, 0],
        },
        rust::CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                rust::EntCode { value: 0, bits: 2 },
                rust::EntCode { value: 1, bits: 2 },
                rust::EntCode { value: 2, bits: 2 },
                rust::EntCode { value: 6, bits: 3 },
                rust::EntCode { value: 7, bits: 3 },
                rust::EntCode {
                    value: 127,
                    bits: 7,
                },
            ],
            pref2bit: vec![0, 1, 2, 0, 0, 0, 3, 4],
        },
        rust::CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                rust::EntCode { value: 0, bits: 2 },
                rust::EntCode { value: 1, bits: 2 },
                rust::EntCode { value: 4, bits: 3 },
                rust::EntCode { value: 5, bits: 3 },
                rust::EntCode { value: 6, bits: 3 },
                rust::EntCode { value: 7, bits: 3 },
                rust::EntCode {
                    value: 255,
                    bits: 8,
                },
            ],
            pref2bit: vec![0, 1, 0, 0, 2, 3, 4, 5],
        },
        rust::CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                rust::EntCode { value: 0, bits: 2 },
                rust::EntCode { value: 2, bits: 3 },
                rust::EntCode { value: 3, bits: 3 },
                rust::EntCode { value: 4, bits: 3 },
                rust::EntCode { value: 5, bits: 3 },
                rust::EntCode { value: 6, bits: 3 },
                rust::EntCode { value: 7, bits: 3 },
                rust::EntCode {
                    value: 511,
                    bits: 9,
                },
            ],
            pref2bit: vec![0, 0, 1, 2, 3, 4, 5, 6],
        },
        rust::CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                rust::EntCode { value: 0, bits: 3 },
                rust::EntCode { value: 1, bits: 3 },
                rust::EntCode { value: 2, bits: 3 },
                rust::EntCode { value: 3, bits: 3 },
                rust::EntCode { value: 4, bits: 3 },
                rust::EntCode { value: 5, bits: 3 },
                rust::EntCode { value: 6, bits: 3 },
                rust::EntCode { value: 7, bits: 3 },
                rust::EntCode {
                    value: 1023,
                    bits: 10,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 5, 6, 7],
        },
        rust::CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                rust::EntCode { value: 0, bits: 3 },
                rust::EntCode { value: 1, bits: 3 },
                rust::EntCode { value: 2, bits: 3 },
                rust::EntCode { value: 3, bits: 3 },
                rust::EntCode { value: 4, bits: 3 },
                rust::EntCode { value: 5, bits: 3 },
                rust::EntCode { value: 6, bits: 3 },
                rust::EntCode { value: 14, bits: 4 },
                rust::EntCode { value: 15, bits: 4 },
                rust::EntCode {
                    value: 4095,
                    bits: 12,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 7, 8],
        },
        rust::CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                rust::EntCode { value: 0, bits: 3 },
                rust::EntCode { value: 1, bits: 3 },
                rust::EntCode { value: 2, bits: 3 },
                rust::EntCode { value: 3, bits: 3 },
                rust::EntCode { value: 4, bits: 3 },
                rust::EntCode { value: 5, bits: 3 },
                rust::EntCode { value: 12, bits: 4 },
                rust::EntCode { value: 13, bits: 4 },
                rust::EntCode { value: 14, bits: 4 },
                rust::EntCode { value: 15, bits: 4 },
                rust::EntCode {
                    value: 8191,
                    bits: 13,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 6, 7, 8, 9],
        },
        rust::CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                rust::EntCode { value: 0, bits: 3 },
                rust::EntCode { value: 1, bits: 3 },
                rust::EntCode { value: 2, bits: 3 },
                rust::EntCode { value: 3, bits: 3 },
                rust::EntCode { value: 4, bits: 3 },
                rust::EntCode { value: 10, bits: 4 },
                rust::EntCode { value: 11, bits: 4 },
                rust::EntCode { value: 12, bits: 4 },
                rust::EntCode { value: 13, bits: 4 },
                rust::EntCode { value: 14, bits: 4 },
                rust::EntCode { value: 15, bits: 4 },
                rust::EntCode {
                    value: 16383,
                    bits: 14,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 5, 6, 7, 8, 9, 10],
        },
        rust::CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                rust::EntCode { value: 0, bits: 3 },
                rust::EntCode { value: 1, bits: 3 },
                rust::EntCode { value: 2, bits: 3 },
                rust::EntCode { value: 3, bits: 3 },
                rust::EntCode { value: 8, bits: 4 },
                rust::EntCode { value: 9, bits: 4 },
                rust::EntCode { value: 10, bits: 4 },
                rust::EntCode { value: 11, bits: 4 },
                rust::EntCode { value: 12, bits: 4 },
                rust::EntCode { value: 13, bits: 4 },
                rust::EntCode { value: 14, bits: 4 },
                rust::EntCode { value: 15, bits: 4 },
                rust::EntCode {
                    value: 32767,
                    bits: 15,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 0, 0, 0, 0, 4, 5, 6, 7, 8, 9, 10, 11],
        },
        rust::CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                rust::EntCode { value: 0, bits: 3 },
                rust::EntCode { value: 1, bits: 3 },
                rust::EntCode { value: 2, bits: 3 },
                rust::EntCode { value: 6, bits: 4 },
                rust::EntCode { value: 7, bits: 4 },
                rust::EntCode { value: 8, bits: 4 },
                rust::EntCode { value: 9, bits: 4 },
                rust::EntCode { value: 10, bits: 4 },
                rust::EntCode { value: 11, bits: 4 },
                rust::EntCode { value: 12, bits: 4 },
                rust::EntCode { value: 13, bits: 4 },
                rust::EntCode { value: 14, bits: 4 },
                rust::EntCode { value: 15, bits: 4 },
                rust::EntCode {
                    value: 65535,
                    bits: 16,
                },
            ],
            pref2bit: vec![0, 1, 2, 0, 0, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        },
        rust::CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                rust::EntCode { value: 0, bits: 3 },
                rust::EntCode { value: 1, bits: 3 },
                rust::EntCode { value: 4, bits: 4 },
                rust::EntCode { value: 5, bits: 4 },
                rust::EntCode { value: 6, bits: 4 },
                rust::EntCode { value: 7, bits: 4 },
                rust::EntCode { value: 8, bits: 4 },
                rust::EntCode { value: 9, bits: 4 },
                rust::EntCode { value: 10, bits: 4 },
                rust::EntCode { value: 11, bits: 4 },
                rust::EntCode { value: 12, bits: 4 },
                rust::EntCode { value: 13, bits: 4 },
                rust::EntCode { value: 14, bits: 4 },
                rust::EntCode { value: 15, bits: 4 },
                rust::EntCode {
                    value: 131071,
                    bits: 17,
                },
            ],
            pref2bit: vec![0, 1, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        },
        rust::CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                rust::EntCode { value: 0, bits: 3 },
                rust::EntCode { value: 2, bits: 4 },
                rust::EntCode { value: 3, bits: 4 },
                rust::EntCode { value: 4, bits: 4 },
                rust::EntCode { value: 5, bits: 4 },
                rust::EntCode { value: 6, bits: 4 },
                rust::EntCode { value: 7, bits: 4 },
                rust::EntCode { value: 8, bits: 4 },
                rust::EntCode { value: 9, bits: 4 },
                rust::EntCode { value: 10, bits: 4 },
                rust::EntCode { value: 11, bits: 4 },
                rust::EntCode { value: 12, bits: 4 },
                rust::EntCode { value: 13, bits: 4 },
                rust::EntCode { value: 14, bits: 4 },
                rust::EntCode { value: 15, bits: 4 },
                rust::EntCode {
                    value: 262143,
                    bits: 18,
                },
            ],
            pref2bit: vec![0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
        },
        rust::CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                rust::EntCode { value: 0, bits: 4 },
                rust::EntCode { value: 1, bits: 4 },
                rust::EntCode { value: 2, bits: 4 },
                rust::EntCode { value: 3, bits: 4 },
                rust::EntCode { value: 4, bits: 4 },
                rust::EntCode { value: 5, bits: 4 },
                rust::EntCode { value: 6, bits: 4 },
                rust::EntCode { value: 7, bits: 4 },
                rust::EntCode { value: 8, bits: 4 },
                rust::EntCode { value: 9, bits: 4 },
                rust::EntCode { value: 10, bits: 4 },
                rust::EntCode { value: 11, bits: 4 },
                rust::EntCode { value: 12, bits: 4 },
                rust::EntCode { value: 13, bits: 4 },
                rust::EntCode { value: 14, bits: 4 },
                rust::EntCode { value: 15, bits: 4 },
                rust::EntCode {
                    value: 524287,
                    bits: 19,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        },
    ]
}
