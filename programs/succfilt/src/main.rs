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

//! SUCCFILT - FX domain correlation coefficient filter.
//! Port of su/main/filters/succfilt.c.
//! Credits: Potash Corporation: Balazs Nemeth, Saskatoon Canada, c. 2008.

use oxwave_core::prelude::*;
use oxwave_gather_io::*;
use oxwave_pfa_fft::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUCCFILT -  FX domain Correlation Coefficient FILTER\t\t\t",
    "                       \t\t\t\t\t\t",
    "   sucff < stdin > stdout [optional parameters]\t\t\t",
    "                       \t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " cch=1.0\t\tCorrelation coefficient high pass value\t\t",
    " ccl=0.3\t\tCorrelation coefficient low pass value\t\t",
    " key=ep\t\tensemble identifier\t\t\t\t",
    " padd=25\t\tFFT padding in percentage\t\t\t",
    "                       \t\t\t\t\t\t",
    " Notes:                       \t\t\t\t\t\t",
    " This program uses \"get_gather\" and \"put_gather\" so requires that\t",
    " the  data be sorted into ensembles designated by \"key\", with the ntr",
    " field set to the number of traces in each respective ensemble.  \t",
    "                       \t\t\t\t\t\t",
    " Example:                     \t\t\t\t\t\t",
    " susort ep offset < data.su > datasorted.su\t\t\t\t",
    " suputgthr dir=Data verbose=1 < datasorted.su\t\t\t\t",
    " sugetgthr dir=Data verbose=1 > dataupdated.su\t\t\t\t",
    " succfilt  < dataupdated.su > ccfiltdata.su\t\t\t\t",
    "                       \t\t\t\t\t\t",
    " (Work in progress, editing required)                 \t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key = String::new();
    if getparstring("key", &mut key) == 0 { key = "ep".to_string(); }
    let mut padd = 25.0f32;
    getparfloat("padd", &mut padd);
    padd = (1.0 + (padd as f64) / 100.0) as f32;
    let mut cch = 1.0f32;
    let mut ccl = 0.3f32;
    let mut verbose = 0;
    getparfloat("cch", &mut cch);
    getparfloat("ccl", &mut ccl);
    getparint("verbose", &mut verbose);
    let mut state = GatherState::default();
    let mut type_ = "";
    let mut val = Value::default();
    let mut nt = 0;
    let mut ntr = 0;
    let mut dt = 0.0;
    let mut first = 0;
    let mut rec_o = get_gather(&mut state, &key, &mut type_, &mut val, &mut nt, &mut ntr, &mut dt, &mut first);
    if ntr == 0 { err!("Can't get first record\n"); }
    let nfft = npfar((nt as f32 * padd) as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    // The C copy and zero-fill are out of bounds when padding makes nfft < nt.
    if nfft < nt { err!("Padded nt=%d--too small", nfft); }
    let nf = (nfft / 2 + 1) as usize;
    let snfft = (1.0 / nfft as f64) as f32;
    let mut rt = vec![0.0f32; nfft as usize];
    let mut ng = 0;
    loop {
        ng += 1;
        let rec = rec_o.as_mut().unwrap();
        let mut fd = vec![vec![Complex::default(); nf]; ntr as usize];
        let mut cc = vec![vec![0.0f32; nf]; ntr as usize];
        for itr in 0..ntr as usize {
            rt[..nt as usize].copy_from_slice(&rec[itr].data[..nt as usize]);
            rt[nt as usize..].fill(0.0);
            pfarc(1, &rt, &mut fd[itr]);
        }
        for itr in 0..ntr as usize - 1 {
            for ifr in 0..nf - 1 {
                let a = fd[itr][ifr];
                let b = fd[itr + 1][ifr];
                cc[itr][ifr] = ((a.i as f64).atan2(a.r as f64) - (b.i as f64).atan2(b.r as f64)).cos() as f32;
            }
        }
        for itr in 0..ntr as usize - 1 {
            for ifr in 0..nf - 1 {
                if cc[itr][ifr] > cch || cc[itr][ifr] < ccl {
                    fd[itr][ifr].r = 0.0;
                    fd[itr][ifr].i = 0.0;
                }
            }
        }
        for itr in 0..ntr as usize {
            pfacr(-1, &fd[itr], &mut rt);
            for it in 0..nt as usize { rec[itr].data[it] = rt[it] * snfft; }
        }
        put_gather(rec, nt);
        rec_o = get_gather(&mut state, &key, &mut type_, &mut val, &mut nt, &mut ntr, &mut dt, &mut first);
        if verbose != 0 { warn!(" %d %d\n", ng, ntr); }
        if ntr == 0 { break; }
    }
    warn!("Number of gathers %10d\n", ng);
    exit(CWP_Exit());
}
