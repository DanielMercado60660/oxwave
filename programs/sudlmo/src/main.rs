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

//! SUDLMO - Dynamic Linear Move Out Correction for Surface Waves.
//! Port of su/main/stretching_moveout_resamp/sudlmo.c (Revision 1.1).
//! Credits: BRGM: Adnand Bitri, October 1998; based on sufrac.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cmplx, cmul, crmul, cwp_cexp, Complex};
use oxwave_index_search::xindex;
use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::{npfaro, pfarc, pfacr};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUDLMO -- Dynamic Linear Move Out Correction for Surface Waves  \t",
    "\t\t\t\t\t\t\t\t\t",
    " sudlmo <indata >outdata [optional parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "\tvnmo=400\tVelocity corresponding to  fnmo\t\t\t",
    "\tfnmo=0.0\tFrequency corresponding to  vnmo\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tdt=(from header)\ttime sample interval (in seconds)\t",
    "       cm = 0 \t\t\tfor offset in m   m                     ",
    "       \t\t       1 for offset in cm                     ",
    "       invert = 0             1 to perform invers DLMO        \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Algorithm:\t\t\t\t\t\t\t\t",
    "\tGdlmo(t,x) = Re[INVFTT{ ( (sign) iw*off/Cf*(Gdlmo(f,x)}]\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn interpovv(nt: usize, cdp: &[f32], ovv: &[Vec<f32>], cdpt: f32, ovvt: &mut [f32], indx: &mut i32) {
    if cdpt <= cdp[0] {
        for it in 0..nt { ovvt[it] = ovv[0][it]; }
    } else if cdpt >= cdp[cdp.len()-1] {
        for it in 0..nt { ovvt[it] = ovv[cdp.len()-1][it]; }
    } else {
        xindex(cdp, cdpt, indx);
        let i = *indx as usize;
        let a1 = (cdp[i+1]-cdpt)/(cdp[i+1]-cdp[i]);
        let a2 = (cdpt-cdp[i])/(cdp[i+1]-cdp[i]);
        for it in 0..nt { ovvt[it] = a1*ovv[i][it]+a2*ovv[i+1][it]; }
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if dt == 0.0 { err!("dt field is zero and not getparred"); }
    let mut invert = 0i32;
    let mut cm = 0i32;
    getparint("invert", &mut invert);
    getparint("cm", &mut cm);
    let nfft = npfaro(nt as i32, 2 * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d -- too big", nfft); }
    let nfft = nfft as usize;
    let nf = nfft/2+1;
    let onfft = (1.0f64 / nfft as f64) as f32;
    let domega = (2.0 * std::f64::consts::PI * onfft as f64 / dt as f64) as f32;
    let mut ncdp = countparval("cdp");
    if ncdp > 0 {
        if countparname("vnmo") != ncdp { err!("a vnmo array must be specified for each cdp"); }
        if countparname("fnmo") != ncdp { err!("a tnmo array must be specified for each cdp"); }
    } else {
        ncdp = 1;
        if countparname("vnmo") > 1 { err!("only one (or no) vnmo array must be specified"); }
        if countparname("fnmo") > 1 { err!("only one (or no) tnmo array must be specified"); }
    }
    let ncdp = ncdp as usize;
    let mut cdp = vec![0.0f32; ncdp];
    if getparfloat("cdp", &mut cdp) == 0 { cdp[0] = tr.cdp as f32; }
    // C passes nt rather than nf on later traces and accesses outside the
    // arrays. Extra storage avoids corruption; only nf entries are used.
    let mut ovv = vec![vec![0.0f32; nf.max(nt)]; ncdp];
    for icdp in 0..ncdp {
        let mut nvnmo = countnparval((icdp+1) as i32, "vnmo");
        let mut ntnmo = countnparval((icdp+1) as i32, "fnmo");
        if nvnmo != ntnmo && !(ncdp == 1 && nvnmo == 1 && ntnmo == 0) { err!("number of vnmo and tnmo values must be equal"); }
        if nvnmo == 0 { nvnmo = 1; }
        if ntnmo == 0 { ntnmo = nvnmo; }
        let mut vnmo = vec![0.0f32; nvnmo as usize];
        let mut fnmo = vec![0.0f32; nvnmo as usize];
        if getnparfloat((icdp+1) as i32, "vnmo", &mut vnmo) == 0 { vnmo[0] = 400.0; }
        if getnparfloat((icdp+1) as i32, "fnmo", &mut fnmo) == 0 { fnmo[0] = 0.0; }
        for it in 0..ntnmo as usize { fnmo[it] = (fnmo[it] as f64 * (2.0 * std::f64::consts::PI)) as f32; }
        for it in 1..ntnmo as usize { if fnmo[it] <= fnmo[it-1] { err!("tnmo values must increase monotonically"); } }
        let mut tn = 0.0f32;
        for it in 0..nf {
            let mut v = [0.0f32];
            intlin(&fnmo, &vnmo, vnmo[0], vnmo[nvnmo as usize-1], &[tn], &mut v);
            ovv[icdp][it] = (1.0f64 / v[0] as f64) as f32;
            tn += domega;
        }
    }
    for jcdp in 1..ncdp {
        let mut i = jcdp;
        while i > 0 && cdp[i-1] > cdp[i] { cdp.swap(i-1,i); ovv.swap(i-1,i); i -= 1; }
    }
    let mut ovvt = vec![0.0f32; nf.max(nt)];
    let mut indx = 0;
    interpovv(nf, &cdp, &ovv, tr.cdp as f32, &mut ovvt, &mut indx);
    let oldcdp = tr.cdp;
    let oldoffset = (tr.offset as i64) - 1;
    let mut rt = vec![0.0f32; nfft];
    let mut ct = vec![Complex::default(); nf];
    loop {
        let newsloth = tr.cdp != oldcdp && ncdp > 1;
        if newsloth { interpovv(nt, &cdp, &ovv, tr.cdp as f32, &mut ovvt, &mut indx); }
        if newsloth || tr.offset as i64 != oldoffset {
            let mut doffs = (tr.offset as f32).abs();
            if cm == 1 { doffs /= 100.0; }
            rt[..nt].copy_from_slice(&tr.data[..nt]);
            rt[nt..].fill(0.0);
            pfarc(1, &rt, &mut ct);
            for i in 0..nf {
                let omega = i as f32 * domega;
                let phase = if invert == 0 { (-1.0f64 * omega as f64 * ovvt[i] as f64 * doffs as f64) as f32 } else { (1.0f64 * omega as f64 * ovvt[i] as f64 * doffs as f64) as f32 };
                let filt = crmul(cwp_cexp(crmul(cmplx(0.0, 1.0), phase)), onfft);
                ct[i] = cmul(ct[i], filt);
            }
        }
        pfacr(-1, &ct, &mut rt);
        tr.data[..nt].copy_from_slice(&rt[..nt]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
