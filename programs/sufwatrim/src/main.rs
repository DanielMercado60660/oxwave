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

//! SUFWATRIM - FX domain alpha trim, port of su/main/filters/sufwatrim.c (Revision 1.2).
//! Credits: Potash Corporation of Saskatchewan, Balasz Nemeth; given to CWP in 2008.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::Complex;
use oxwave_gather_io::*;
use oxwave_pfa_fft::*;
use oxwave_sort::qkisort;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUFWATRIM - FX domain Alpha TRIM\t\t\t\t\t",
    "\t\t\t \t\t\t\t\t\t",
    "  sufwatrim  <stdin > stdout [optional parameters]\t\t\t",
    "\t\t\t \t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " key=key1,key2,..\tHeader words defining mixing dimesnion\t\t",
    " dx=d1,d2,..\t\tDistance units for each header word\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " keyg=ep\t\tHeader word indicating the start of gather\t",
    " vf=0\t\t\t=1 Do a frequency dependent mix\t\t\t",
    " vmin=5000\t\tVelocity of the reflection slope\t\t",
    "\t\t\tthan should not be attenuated\t\t\t",
    " Notes:\t\t \t\t\t\t\t\t",
    " Trace with the header word mark set to one will be \t\t\t",
    " the output trace \t\t\t\t\t\t\t",
];

fn alpha_trim_w(a: &[f32], w: &[f32], n: usize) -> f32 {
    let ist = (n as f64 * (15.0 / 100.0) + 0.5) as usize;
    let mut ind: Vec<i32> = (0..n as i32).collect();
    qkisort(a, &mut ind);
    let mut mean = 0.0f32;
    let mut ic = 0;
    for i in ist..n-ist {
        mean += a[ind[i] as usize] * w[ind[i] as usize];
        ic += 1;
    }
    mean / ic as f32
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut keyg = String::new();
    if getparstring("keyg", &mut keyg) == 0 { keyg = "ep".into(); }
    let mut vf = 1;
    getparint("vf", &mut vf);
    let mut vmin = 5000.0f32;
    getparfloat("vmin", &mut vmin);
    let mut padd = 25.0f32;
    getparfloat("padd", &mut padd);
    padd = (1.0 + padd as f64 / 100.0) as f32;
    let nd = countparval("key") as usize;
    let mut key = Vec::new();
    getparstringarray("key", &mut key);
    let mut index = Vec::new();
    let mut types = Vec::new();
    for k in &key { types.push(hdtype(k)); index.push(getindex(k)); }
    let mut dx = vec![0.0f32; nd];
    if getparfloat("dx", &mut dx) == 0 { err!("must specify %s=", "dx"); }
    let mut fac = 1.0f32;
    getparfloat("fac", &mut fac);
    if fac <= 1.0 { fac = 1.0; }
    let _ = (vf, vmin, fac);
    let mut state = GatherState::default();
    let mut typeg = "";
    let mut valg = Value::default();
    let (mut nt, mut ntr, mut first) = (0, 0, 0);
    let mut dt = 0.0;
    let mut rec = get_gather(&mut state, &keyg, &mut typeg, &mut valg, &mut nt, &mut ntr, &mut dt, &mut first);
    if ntr == 0 { err!("Can't get first record\n"); }
    let nfft = npfar((nt as f32 * padd) as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d--too big", nfft); }
    let nf = (nfft / 2 + 1) as usize;
    let snfft = (1.0f64 / nfft as f64) as f32;
    let mut rt = vec![0.0f32; nfft as usize];
    let mut tr = Trace::new();
    let mut ng = 0;
    while let Some(records) = rec {
        ng += 1;
        let n = ntr as usize;
        let mut fd = vec![vec![Complex::default(); nf]; n];
        for itr in 0..n {
            rt[..nt as usize].copy_from_slice(&records[itr].data[..nt as usize]);
            rt[nt as usize..].fill(0.0);
            pfarc(1, &rt, &mut fd[itr]);
        }
        let mut imx = 0;
        for itr in 0..n { if records[itr].mark != 0 { imx = itr; break; } }
        tr.hdr = records[imx].hdr;
        let mut wht = vec![1.0f32; n];
        for itr in 0..n {
            if itr == imx { continue; }
            let mut dist = 0.0f64;
            for i in 0..nd {
                let mut val = Value::default();
                gethval(&records[imx], index[i], &mut val);
                let o = vtof(types[i], val);
                gethval(&records[itr], index[i], &mut val);
                let x = vtof(types[i], val);
                let du = (o-x)/dx[i];
                dist += (du*du) as f64;
            }
            let _ = dist.sqrt();
            wht[itr] = 1.0;
        }
        let mut ctmix = vec![Complex::default(); nf];
        let mut itrm = vec![0.0f32; n];
        let mut rtrm = vec![0.0f32; n];
        for ifr in 0..nf {
            for itr in 0..n { itrm[itr] = fd[itr][ifr].i; rtrm[itr] = fd[itr][ifr].r; }
            ctmix[ifr].i = alpha_trim_w(&itrm, &wht, n);
            ctmix[ifr].r = alpha_trim_w(&rtrm, &wht, n);
        }
        pfacr(-1, &ctmix, &mut rt);
        for it in 0..nt as usize { tr.data[it] = rt[it] * snfft; }
        puttr(&mut tr);
        rec = get_gather(&mut state, &keyg, &mut typeg, &mut valg, &mut nt, &mut ntr, &mut dt, &mut first);
        fprintf!(stderr(), " %d %d\n", ng, ntr);
    }
    warn!("Number of gathers %10d\n", ng);
    exit(CWP_Exit());
}
