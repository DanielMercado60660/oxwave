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
//! SULFAF - low frequency array forming, port of su/main/filters/sulfaf.c (Revision 1.3).

use oxwave_core::prelude::*;
use oxwave_gather_io::*;
use oxwave_pfa_fft::*;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SULFAF -  Low frequency array forming\t\t\t\t\t",
    "\t\t\t  \t\t\t\t\t\t",
    "  sulfaf < stdin > stdout [optional parameters]\t\t\t",
    "\t\t\t  \t\t\t\t\t\t",
    " Optional Parameters:\t  \t\t\t\t\t\t",
    " key=ep\theader keyword describing the gathers\t\t\t",
    " f1=3\t\tlower frequency\tcutof\t\t\t\t\t",
    " f2=20\t\thigh frequency cutof\t\t\t\t\t",
    " fr=5\t\tfrequency ramp\t\t\t\t\t\t",
    "\t\t\t  \t\t\t\t\t\t",
    " vel=1000\tsurface wave velocity\t\t\t\t\t",
    " dx=10\t\ttrace spacing\t\t\t\t\t\t",
    " maxmix=tr.ntr\tdefault is the entire gather\t\t\t\t",
    " adb=1.0\tadd back ratio 1.0=pure filtered 0.0=origibal\t\t",
    "\t\t\t  \t\t\t\t\t\t",
    " tri=0\t\t1 Triangular weight in mixing window\t\t\t",
    "\t\t\t  \t\t\t\t\t\t",
    " Notes:\t\t  \t\t\t\t\t\t",
    " The traces transformed into the freqiency domain\t\t\t",
    " where a trace mix is performed in the specifed frequency range\t",
    " as Mix=vel/(freq*dx)\t\t\t\t\t\t\t",
    " This program uses \"get_gather\" and \"put_gather\" so requires that  ",
    " the  data be sorted into ensembles designated by \"key\", with the ntr",
    " field set to the number of traces in each respective ensemble.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t ",
    " susort ep offset < data.su > datasorted.su\t\t\t\t",
    " suputgthr dir=Data verbose=1 < datasorted.su\t\t\t  ",
    " sugetgthr dir=Data verbose=1 > dataupdated.su\t\t\t ",
    " sulfaf  < dataupdated.su > ccfiltdata.su\t\t\t\t",
    "\t\t\t  \t\t\t\t\t\t",
    " (Work in progress, editing required)                 \t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key = String::from("ep");
    let (mut f1, mut f2, mut dx, mut vel, mut fr, mut adb) = (3.0f32, 20.0f32, 10.0f32, 1000.0f32, 5.0f32, 1.0f32);
    let (mut maxmix, mut tri, mut verbose) = (-1, 0, 0);
    getparstring("key", &mut key);
    getparfloat("f1", &mut f1);
    getparfloat("f2", &mut f2);
    getparfloat("dx", &mut dx);
    getparfloat("vel", &mut vel);
    getparfloat("fr", &mut fr);
    getparint("maxmix", &mut maxmix);
    getparint("tri", &mut tri);
    getparfloat("adb", &mut adb);
    getparint("verbose", &mut verbose);
    let mut state = GatherState::default();
    let (mut type_, mut val, mut nt, mut ntr, mut dt, mut first) = ("", Value::default(), 0, 0, 0.0, 0);
    let mut rec = get_gather(&mut state, &key, &mut type_, &mut val, &mut nt, &mut ntr, &mut dt, &mut first);
    if ntr == 0 { err!("Can't get first record\n"); }
    let nfft = npfaro(nt, nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    let nf = (nfft / 2 + 1) as usize;
    let snfft = (1.0f64 / nfft as f64) as f32;
    let d1 = (1.0f64 / (nfft as f32 * dt) as f64) as f32;
    let mut ct = vec![Complex::default(); nf];
    let mut rt = vec![0.0f32; nfft as usize];
    let if1 = nint((f1 / d1) as f64);
    let if2 = nint((f2 / d1) as f64);
    let ifr = nint((fr / d1) as f64);
    // The C program leaves these array accesses unchecked. Reject only ranges
    // that would allocate an invalid size or access outside an array.
    let end = match if2.checked_add(ifr) {
        Some(end) if end >= 0 => end,
        _ => err!("Invalid frequency bin range"),
    };
    if (if1 < end && (if1 < 0 || end as usize > nf))
        || (if1 < if2 && (if1 < 0 || if2 > end || if2 as usize > nf))
        || (if2 < end && (if2 < 0 || end as usize > nf))
    {
        err!("Invalid frequency bin range");
    }
    let mut ng = 0;
    while let Some(mut records) = rec {
        if nt < 0 || nt > nfft {
            err!("Gather sample count exceeds FFT size");
        }
        if maxmix == -1 { maxmix = ntr; }
        ng += 1;
        let mut ffdr = vec![vec![0.0f32; nf]; ntr as usize];
        let mut ffdi = vec![vec![0.0f32; nf]; ntr as usize];
        let mut ffdrm = vec![vec![0.0f32; end as usize]; ntr as usize];
        let mut ffdim = vec![vec![0.0f32; end as usize]; ntr as usize];
        for itr in 0..ntr as usize {
            rt[..nt as usize].copy_from_slice(&records[itr].data[..nt as usize]);
            rt[nt as usize..].fill(0.0);
            pfarc(1, &rt, &mut ct);
            for iw in 0..nf { ffdr[itr][iw] = ct[iw].r; ffdi[itr][iw] = ct[iw].i; }
        }
        for iw in if1..end {
            let v = nint(((vel / iw as f32) * d1 * dx) as f64);
            let mut mix = if v < maxmix { v } else { maxmix };
            if mix % 2 == 0 { mix -= 1; }
            if verbose != 0 { warn!(" %f %d", iw as f32 * d1, mix); }
            for itr in 0..ntr {
                let ims = (itr - mix / 2).max(0);
                let ime = (ims + mix).min(ntr - 1);
                let (mut tmpr, mut tmpi, mut wh) = (0.0f32, 0.0f32, 1.0f32);
                let ws = mix / 2;
                let mut nmix = 0;
                let mut iww = -mix / 2;
                for itrm in ims..ime {
                    nmix += 1;
                    if tri != 0 { wh = (ws - iww.abs()) as f32; }
                    tmpr += ffdr[itrm as usize][iw as usize] * wh;
                    tmpi += ffdi[itrm as usize][iw as usize] * wh;
                    iww += 1;
                }
                ffdrm[itr as usize][iw as usize] = tmpr / nmix as f32;
                ffdim[itr as usize][iw as usize] = tmpi / nmix as f32;
            }
        }
        for iw in if1..if2 {
            for itr in 0..ntr as usize {
                let iw = iw as usize;
                ffdr[itr][iw] = ((ffdrm[itr][iw] * adb) as f64 + ffdr[itr][iw] as f64 * (1.0 - adb as f64)) as f32;
                ffdi[itr][iw] = ((ffdim[itr][iw] * adb) as f64 + ffdi[itr][iw] as f64 * (1.0 - adb as f64)) as f32;
            }
        }
        for iw in if2..end {
            let wh = (1.0 - ((iw-if2) as f32 / ifr as f32) as f64) as f32;
            for itr in 0..ntr as usize {
                let iw = iw as usize;
                ffdr[itr][iw] = (((wh * ffdrm[itr][iw]) as f64 + ffdr[itr][iw] as f64 * (1.0 - wh as f64)) * adb as f64 + ffdr[itr][iw] as f64 * (1.0 - adb as f64)) as f32;
                ffdi[itr][iw] = (((wh * ffdim[itr][iw]) as f64 + ffdi[itr][iw] as f64 * (1.0 - wh as f64)) * adb as f64 + ffdi[itr][iw] as f64 * (1.0 - adb as f64)) as f32;
            }
        }
        for itr in 0..ntr as usize {
            for iw in 0..nf { ct[iw].r = ffdr[itr][iw] * snfft; ct[iw].i = ffdi[itr][iw] * snfft; }
            pfacr(-1, &ct, &mut rt);
            records[itr].data[..nt as usize].copy_from_slice(&rt[..nt as usize]);
        }
        put_gather(&records, nt);
        rec = get_gather(&mut state, &key, &mut type_, &mut val, &mut nt, &mut ntr, &mut dt, &mut first);
    }
    warn!("Number of gathers %10d\n", ng);
    exit(CWP_Exit());
}
