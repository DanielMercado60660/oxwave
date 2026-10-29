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

//! SUFWMIX - FX domain multidimensional weighted mix.
//! Port of su/main/operations/sufwmix.c (Revision 1.4).
//! Credits: Potash Corporation, Balazs Nemeth, Saskatoon Saskatchewan CA,
//! given to CWP in 2008.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
use oxwave_gather_io::*;
use oxwave_pfa_fft::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUFWMIX -  FX domain multidimensional Weighted Mix\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tsufwmix < stdin > stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " key=key1,key2,..\tHeader words defining mixing dimension\t\t",
    " dx=d1,d2,..\t\tDistance units for each header word\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " keyg=ep\t\tHeader word indicating the start of gather\t",
    " vf=0\t\t\t=1 Do a frequency dependent mix\t\t\t",
    " vmin=5000\t\tVelocity of the reflection slope\t\t",
    "\t\t\tthan should not be attenuated\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Trace with the header word mark set to one will be\t\t\t",
    " the output trace \t\t\t\t\t\t\t",
    "  (a work in progress)\t\t\t\t\t\t\t",
];

fn n_distance(rec: &[Trace], index: &[i32], types: &[&str], dx: &[f32], imx: usize, itr: usize) -> f32 {
    let mut dist = 0.0f64;
    let mut val = Value::default();
    for i in 0..index.len() {
        gethval(&rec[imx], index[i], &mut val);
        let o = vtof(types[i], val);
        gethval(&rec[itr], index[i], &mut val);
        let x = vtof(types[i], val);
        let du = (o - x) / dx[i];
        dist += (du * du) as f64;
    }
    dist.sqrt() as f32
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut keyg = String::new();
    if getparstring("keyg", &mut keyg) == 0 { keyg = "ep".to_string(); }
    let mut vf = 1i32;
    getparint("vf", &mut vf);
    let mut vmin = 5000.0f32;
    getparfloat("vmin", &mut vmin);
    let mut padd = 25.0f32;
    getparfloat("padd", &mut padd);
    padd = (1.0 + padd as f64 / 100.0) as f32;
    let nd = countparval("key") as usize;
    let mut key = Vec::new();
    getparstringarray("key", &mut key);
    let mut types = Vec::new();
    let mut index = Vec::new();
    for k in key.iter().take(nd) { types.push(hdtype(k)); index.push(getindex(k)); }
    let mut dx = vec![0.0f32; nd];
    if getparfloat("dx", &mut dx) == 0 { err!("must specify %s=", "dx"); }
    let mut fac = 1.0f32;
    getparfloat("fac", &mut fac);
    if !(fac > 1.0) { fac = 1.0; }
    checkpars();
    let mut state = GatherState::default();
    let mut typeg = "";
    let mut valg = Value::default();
    let (mut nt, mut ntr, mut first, mut dt) = (0, 0, 0, 0.0f32);
    let mut rec = get_gather(&mut state, &keyg, &mut typeg, &mut valg, &mut nt, &mut ntr, &mut dt, &mut first);
    if ntr == 0 { err!("Can't get first record\n"); }
    let nfft = npfar((nt as f32 * padd) as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    let nf = (nfft / 2 + 1) as usize;
    let snfft = (1.0 / nfft as f64) as f32;
    let d1 = (1.0 / (nfft as f32 * dt) as f64) as f32;
    let mut rt = vec![0.0f32; nfft as usize];
    let mut tr = Trace::new();
    let mut ng = 0;
    while let Some(records) = rec {
        ng += 1;
        let mut fd = vec![vec![Complex::default(); nf]; ntr as usize];
        for itr in 0..ntr as usize {
            rt[..nt as usize].copy_from_slice(&records[itr].data[..nt as usize]);
            rt[nt as usize..].fill(0.0);
            pfarc(1, &rt, &mut fd[itr]);
        }
        let mut imx = 0;
        for itr in 0..ntr as usize {
            if records[itr].mark != 0 { imx = itr; break; }
        }
        let mut ctmix = fd[imx].clone();
        tr.hdr = records[imx].hdr;
        for itr in (0..imx).chain(imx + 1..ntr as usize) {
            let dist = n_distance(&records, &index, &types, &dx, imx, itr);
            let weight = (1.0 / ((dist / fac) as f64 + 1.0)) as f32;
            for ifr in 0..nf {
                let dlim = (vmin as f64 / (2.0 * ifr as f64 * d1 as f64)) as f32;
                if vf == 0 || dlim > dist {
                    ctmix[ifr].r += weight * fd[itr][ifr].r;
                    ctmix[ifr].i += weight * fd[itr][ifr].i;
                }
            }
        }
        // The C code discards crmul's return value, so normalization does nothing.
        pfacr(-1, &ctmix, &mut rt);
        for it in 0..nt as usize { tr.data[it] = rt[it] * snfft; }
        puttr(&mut tr);
        rec = get_gather(&mut state, &keyg, &mut typeg, &mut valg, &mut nt, &mut ntr, &mut dt, &mut first);
        fprintf!(stderr(), " %d %d\n", ng, ntr);
    }
    warn!("Number of gathers %10d\n", ng);
    exit(CWP_Exit());
}
