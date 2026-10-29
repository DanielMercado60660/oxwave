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

//! Port of comp/dwpt/1d/lib/wbuffalloc.c, comp/dwpt/1d/lib/wpc1.c, comp/dwpt/1d/lib/wpc1coding.c, comp/dwpt/1d/lib/wpc1quant.c, comp/dwpt/1d/lib/wpc1trans.c: bookCoder, bookDecoder, buffAlloc1, buffAlloc2, buffFree1, buffFree2, buffRealloc1, entCoder, entDecoder, wavePack1, wpc1Compress, wpc1Decoder, wpc1Dequant, wpc1Encoder, wpc1Free, wpc1Numsmp, wpc1Quant, wpc1Read, wpc1Size, wpc1Trans, wpc1Uncompress, wpc1Write.
//! Credits: Tong Chen (CWP), 1994–1995.
#![allow(non_snake_case)]

mod coding;
pub use coding::{
    bookCoder, bookDecoder, entCoder, entDecoder, wpc1Decoder, wpc1Encoder, BitBuffer, CodeBook,
    EntCode,
};

/// Wavelet configuration.
#[derive(Clone, Debug, Default)]
pub struct Config {
    pub powt: i32,
    pub tileszt: usize,
    pub staget: i32,
    pub nblock: i32,
    pub lblock: usize,
    pub ishift: i32,
    pub len: usize,
    pub filterh: [f32; 12],
    pub filterg: [f32; 12],
}

/// Quantization status.
#[derive(Clone, Debug, Default)]
pub struct Quant {
    pub step: f32,
    pub ave: f32,
    pub numsmp: usize,
    pub nbit: i32,
    pub qx: Vec<i32>,
}

/// Byte buffer and current position.
#[derive(Clone, Debug)]
pub struct Buffer {
    pub code: Vec<u8>,
    pub pos: usize,
}

/// Compressed trace and its quantization metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct Compressed {
    pub numsmp: i32,
    pub overallstep: f32,
    pub overallave: f32,
    pub code: Vec<u8>,
}

/// C prototype: `int wpc1Numsmp(void *wpc1);`
pub fn wpc1Numsmp(w: &Compressed) -> i32 {
    w.numsmp
}

/// C prototype: `int wpc1Size(void *wpc1);`
pub fn wpc1Size(w: &Compressed) -> i32 {
    w.code.len() as i32
}

/// C prototype: `void wpc1Free(void *wpc1);`
pub fn wpc1Free(w: Compressed) {
    drop(w);
}

/// C prototype: `void wpc1Write(void *wpc1, FILE *fp);`
pub fn wpc1Write(w: &Compressed, fp: &mut impl std::io::Write) -> std::io::Result<()> {
    fp.write_all(&w.numsmp.to_ne_bytes())?;
    fp.write_all(&(w.code.len() as i32).to_ne_bytes())?;
    fp.write_all(&w.overallstep.to_ne_bytes())?;
    fp.write_all(&w.overallave.to_ne_bytes())?;
    fp.write_all(&w.code)
}

/// C prototype: `void *wpc1Read(FILE *fp);`
pub fn wpc1Read(fp: &mut impl std::io::Read) -> std::io::Result<Compressed> {
    let mut word = [0u8; 4];
    fp.read_exact(&mut word)?;
    let numsmp = i32::from_ne_bytes(word);
    fp.read_exact(&mut word)?;
    let nsize = i32::from_ne_bytes(word);
    if nsize < 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "negative compressed size",
        ));
    }
    fp.read_exact(&mut word)?;
    let overallstep = f32::from_ne_bytes(word);
    fp.read_exact(&mut word)?;
    let overallave = f32::from_ne_bytes(word);
    let mut code = vec![0; nsize as usize];
    fp.read_exact(&mut code)?;
    Ok(Compressed {
        numsmp,
        overallstep,
        overallave,
        code,
    })
}

/// C prototype: `void *buffAlloc1(int n);`
pub fn buffAlloc1(n: usize) -> Buffer {
    Buffer {
        code: vec![0; n],
        pos: 0,
    }
}
/// C prototype: `void buffRealloc1(void *inbuff, int n);`
pub fn buffRealloc1(buff: &mut Buffer, n: usize) {
    buff.code.resize(n, 0);
    buff.pos = buff.pos.min(n);
}
/// C prototype: `void buffFree1(void *inbuff);`
pub fn buffFree1(buff: Buffer) {
    drop(buff);
}
/// C prototype: `void **buffAlloc2(int n1, int n2);`
pub fn buffAlloc2(n1: usize, n2: usize) -> Vec<Buffer> {
    (0..n1).map(|_| buffAlloc1(n2)).collect()
}
/// C prototype: `void buffFree2(void **inbuff, int n1);`
pub fn buffFree2(buff: Vec<Buffer>) {
    drop(buff);
}

/// C prototype: `void wpc1Quant(float *x, int n, float error, void *inconf, void *qstate);`
pub fn wpc1Quant(x: &mut [f32], n: usize, error: f32, config: &Config, quant: &mut Quant) {
    let nlen = if n == 0 { config.tileszt } else { n };
    let first = if n == 0 { nlen - config.lblock } else { 0 };
    let mut dev = 0.0f32;
    let mut fmax = 0.0f32;
    for &v in &x[..first] {
        let a = if v < 0.0 { -v } else { v };
        if fmax < a {
            fmax = a;
        }
        dev += a * a;
    }
    let mut average = 0.0f32;
    for &v in &x[first..nlen] {
        average += v;
    }
    average /= (nlen - first) as f32;
    for v in &mut x[first..nlen] {
        *v -= average;
        let a = if *v < 0.0 { -*v } else { *v };
        if fmax < a {
            fmax = a;
        }
        dev += a * a;
    }
    quant.numsmp = nlen;
    dev /= nlen as f32;
    dev = (dev as f64).sqrt() as f32;
    quant.ave = average;
    quant.qx.resize(nlen, 0);
    if dev <= f32::MIN_POSITIVE {
        quant.nbit = 0;
        quant.step = f32::MIN_POSITIVE;
        quant.qx.fill(0);
        return;
    }
    let loading = fmax / dev;
    let clip = ((error * dev) as f64 + 0.0 * (loading as f64)) as f32;
    let mut step_factor = (0.26 / clip as f64) as f32;
    let (nbit, level) = loop {
        let bound = fmax * step_factor;
        let level = (bound as f64 + 0.5) as i32;
        let mut nbit = 0;
        if level != 0 {
            while nbit < 31 && level >> nbit != 0 {
                nbit += 1;
            }
            nbit -= 1;
        }
        if nbit <= 14 {
            break (nbit, level);
        }
        eprintln!("Warning, cannot handle error so small, will do my best");
        step_factor = (step_factor as f64 * 0.5) as f32;
    };
    let levelfull = 1 << (nbit + 1);
    if level > 8 {
        if level == levelfull {
            step_factor = (step_factor as f64 * ((level as f64 - 2.0) / level as f64)) as f32;
        } else if level == levelfull - 1 {
            step_factor = (step_factor as f64 * ((level as f64 - 1.0) / level as f64)) as f32;
        }
    }
    quant.nbit = nbit;
    quant.step = step_factor;
    for i in 0..nlen {
        let a = x[i] * step_factor;
        quant.qx[i] = if a > 0.0 {
            (a as f64 + 0.5) as i32
        } else {
            (a as f64 - 0.5) as i32
        };
    }
}

/// C prototype: `void wpc1Dequant(float *x, int n, void *inconf, void *qstate);`
pub fn wpc1Dequant(x: &mut [f32], n: usize, config: &Config, quant: &Quant) {
    let nlen = if n == 0 { config.tileszt } else { n };
    if quant.step <= f32::MIN_POSITIVE {
        x[..nlen].fill(0.0);
    } else {
        let rstep = (1.0 / quant.step as f64) as f32;
        for i in 0..nlen {
            x[i] = quant.qx[i] as f32 * rstep;
        }
    }
    let first = if n == 0 { nlen - config.lblock } else { 0 };
    for v in &mut x[first..nlen] {
        *v += quant.ave;
    }
}

/// Wavelet filter used by the packet transform.
pub struct WaveFilter {
    pub ishift: i32,
    pub h: Vec<f32>,
    pub g: Vec<f32>,
}

/// C prototype: `void wavePack1(float *x, float *y, waveFILTER *filter, int twopow, int stage, int type);`
pub fn wavePack1(
    x: &mut [f32],
    y: &mut [f32],
    filter: &WaveFilter,
    twopow: i32,
    stage: i32,
    typ: i32,
) {
    let n = 1usize << twopow;
    let len = filter.h.len();
    assert_eq!(len, filter.g.len());
    let h = &filter.h;
    let g = &filter.g;
    if stage == 0 {
        if typ == 0 {
            y[..n].copy_from_slice(&x[..n]);
        } else {
            x[..n].copy_from_slice(&y[..n]);
        }
    }
    let mut tmp = vec![0.0f32; n];
    let position = |size: usize| -> usize {
        let mut pos = size as i32;
        while pos < filter.ishift.abs() {
            pos <<= 1;
        }
        (pos + filter.ishift) as usize
    };
    if typ == 0 {
        let pos = position(n);
        for j in 0..n / 2 {
            tmp[j] = 0.0;
            y[j] = 0.0;
            for i in (1..len).step_by(2) {
                let k = (i + 2 * j - 1 + pos) & (n - 1);
                tmp[j] += x[k] * h[i - 1] + x[k + 1] * h[i];
                y[j] += x[k] * g[i - 1] + x[k + 1] * g[i];
            }
        }
        for level in 1..stage {
            let nband = (1usize << level) >> 1;
            let nmid = 1usize << (twopow - level);
            let half = nmid >> 1;
            let pos = position(nmid);
            for band in 0..nband {
                let first = (band * nmid) << 1;
                let mut out = first + nmid;
                for j in 0..half {
                    let jj = out + j;
                    tmp[jj] = 0.0;
                    y[jj] = 0.0;
                    for i in (1..len).step_by(2) {
                        let k = ((i + 2 * j - 1 + pos) & (nmid - 1)) + first;
                        tmp[jj] += tmp[k] * h[i - 1] + tmp[k + 1] * h[i];
                        y[jj] += tmp[k] * g[i - 1] + tmp[k + 1] * g[i];
                    }
                }
                out += half;
                for j in 0..half {
                    let jj = out + j;
                    tmp[jj] = 0.0;
                    y[jj] = 0.0;
                    for i in (1..len).step_by(2) {
                        let k = ((i + 2 * j - 1 + pos) & (nmid - 1)) + first;
                        tmp[jj] += y[k] * g[i - 1] + y[k + 1] * g[i];
                        y[jj] += y[k] * h[i - 1] + y[k + 1] * h[i];
                    }
                }
                for j in 0..half {
                    tmp[first + j] = tmp[out + j];
                    y[first + j] = y[out + j];
                }
            }
        }
        let nmid = 1usize << (twopow - stage);
        for band in 0..((1usize << stage) >> 1) {
            let first = band * (nmid << 1);
            for i in 0..nmid {
                y[first + nmid + i] = tmp[first + i];
            }
        }
    } else {
        x[..n].copy_from_slice(&y[..n]);
        for level in (2..=stage).rev() {
            let nmid = 1usize << (twopow - level);
            let width = nmid << 1;
            let pos = position(width);
            for band in 0..((1usize << level) >> 1) {
                let first = band * width;
                if band % 2 == 0 {
                    for i in 0..nmid {
                        tmp[first + i] = x[first + nmid + i];
                    }
                    for i in 0..nmid {
                        tmp[first + nmid + i] = x[first + i];
                    }
                } else {
                    tmp[first..first + width].copy_from_slice(&x[first..first + width]);
                }
                x[first..first + width].fill(0.0);
                for j in 0..nmid {
                    for i in 0..len {
                        let k = ((2 * j + i + pos) & (width - 1)) + first;
                        x[k] += tmp[first + j] * g[i] + tmp[first + j + nmid] * h[i];
                    }
                }
            }
        }
        tmp[..n].copy_from_slice(&x[..n]);
        x[..n].fill(0.0);
        let pos = position(n);
        for j in 0..n / 2 {
            for i in 0..len {
                let k = (2 * j + i + pos) & (n - 1);
                x[k] += tmp[j] * g[i] + tmp[j + n / 2] * h[i];
            }
        }
    }
}

/// C prototype: `void wpc1Trans(float *x, void *inconf, int type);`
pub fn wpc1Trans(x: &mut [f32], config: &Config, typ: i32) {
    let filter = WaveFilter {
        ishift: config.ishift,
        h: config.filterh[..config.len].to_vec(),
        g: config.filterg[..config.len].to_vec(),
    };
    let mut z = vec![0.0; config.tileszt];
    if typ == 0 {
        wavePack1(x, &mut z, &filter, config.powt, config.staget, typ);
    } else {
        wavePack1(&mut z, x, &filter, config.powt, config.staget, typ);
    }
    x[..config.tileszt].copy_from_slice(&z);
}

fn init_config() -> Config {
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
    Config {
        powt: 8,
        tileszt: 256,
        staget: 4,
        nblock: 16,
        lblock: 16,
        ishift: -4,
        len: 12,
        filterh: h,
        filterg: g,
    }
}

/// C prototype: `void *wpc1Compress(float *f, int n1, float error, int *nbytes);`
/// Mutates unpadded input in place, as C does. `nbytes` excludes the 16 byte header.
pub fn wpc1Compress(f: &mut [f32], error: f32, nbytes: &mut i32) -> Compressed {
    assert!(!f.is_empty());
    let config = init_config();
    let n1 = f.len();
    let mut npartition = n1 / config.tileszt;
    let nleft = n1 - npartition * config.tileszt;
    let mut padded = Vec::new();
    let fpad = if nleft > config.tileszt / 2 {
        npartition += 1;
        padded.resize(npartition * config.tileszt, 0.0);
        padded[..n1].copy_from_slice(f);
        &mut padded[..]
    } else {
        &mut f[..]
    };
    let mut quant = Vec::new();
    for ip in 0..npartition {
        let block = &mut fpad[ip * config.tileszt..(ip + 1) * config.tileszt];
        wpc1Trans(block, &config, 0);
        let mut q = Quant::default();
        wpc1Quant(block, 0, error, &config, &mut q);
        quant.push(q);
    }
    if nleft != 0 && nleft <= config.tileszt / 2 {
        let mut q = Quant::default();
        wpc1Quant(
            &mut fpad[npartition * config.tileszt..],
            nleft,
            error,
            &config,
            &mut q,
        );
        quant.push(q);
    }
    let mut w = Compressed {
        numsmp: n1 as i32,
        overallstep: 0.0,
        overallave: 0.0,
        code: vec![0; n1 * 4],
    };
    wpc1Encoder(&quant, &mut w);
    *nbytes = w.code.len() as i32;
    w
}

/// C prototype: `int wpc1Uncompress(void *w1, float *f);`
/// Returns 1 for a consistent stream, 0 otherwise.
/// C leaves block averages uninitialized for zero-mean streams. This safe port
/// initializes them to zero; `wpc1Decoder` itself preserves caller-supplied values.
pub fn wpc1Uncompress(w: &Compressed, f: &mut [f32]) -> i32 {
    if w.numsmp <= 0 {
        return 0;
    }
    let n1 = w.numsmp as usize;
    assert!(f.len() >= n1);
    let config = init_config();
    let mut npartition = n1 / config.tileszt;
    let nleft = n1 - npartition * config.tileszt;
    let pad = nleft > config.tileszt / 2;
    if pad {
        npartition += 1;
    }
    let mut quant = Vec::new();
    for _ in 0..npartition {
        quant.push(Quant {
            numsmp: config.tileszt,
            qx: vec![0; config.tileszt],
            ..Default::default()
        });
    }
    if nleft != 0 && !pad {
        quant.push(Quant {
            numsmp: nleft,
            qx: vec![0; nleft],
            ..Default::default()
        });
    }
    let mut padded = vec![0.0; if pad { npartition * config.tileszt } else { 0 }];
    let fpad = if pad { &mut padded[..] } else { &mut f[..n1] };
    let retval = wpc1Decoder(&mut quant, w);
    if retval != 0 {
        for ip in 0..npartition {
            let block = &mut fpad[ip * config.tileszt..(ip + 1) * config.tileszt];
            wpc1Dequant(block, 0, &config, &quant[ip]);
            wpc1Trans(block, &config, 1);
        }
        if nleft != 0 && !pad {
            wpc1Dequant(
                &mut fpad[npartition * config.tileszt..],
                nleft,
                &config,
                &quant[npartition],
            );
        }
    }
    if pad {
        f[..n1].copy_from_slice(&padded[..n1]);
    }
    retval
}
