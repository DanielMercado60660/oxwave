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

//! Anisotropic Gaussian beam migration. Port of su/main/migration_inversion/sumiggbzoan.c (Revision 1.5).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: CWP: Tariq Alkhalifah, based on MIGGBZO by Dave Hale.
//! CWP: repackaged as an SU program by John Stockwell, April 2006.

#![allow(non_snake_case)]
use super::*;
pub(super) fn formBeams(
    bwh: f32,
    dxb: f32,
    fmin: f32,
    nt: i32,
    dt: f32,
    ft: f32,
    nx: i32,
    dx: f32,
    fx: f32,
    f: &[Vec<f32>],
    ntau: i32,
    _dtau: f32,
    ftau: f32,
    npx: i32,
    dpx: f32,
    fpx: f32,
    g: &mut [Vec<f32>],
) {
    let wmin = (2.0 * PI * fmin as f64) as f32;
    let pxmax = if dpx < 0.0 {
        fpx
    } else {
        fpx + (npx - 1) as f32 * dpx
    };
    let xmax = if dx < 0.0 {
        fx
    } else {
        fx + (nx - 1) as f32 * dx
    };
    let v = pxmax * xmax;
    let av = if v < 0.0 { -v } else { v };
    let ntpad = (av / dt) as i32;
    let ntfft = npfar((nt + ntpad).max(ntau));
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = 0.0f32;
    let fftscl = (1.0 / ntfft as f64) as f32;
    let mut fpad = vec![0.0f32; ntfft as usize];
    let mut gpad = fpad.clone();
    let mut cf = vec![vec![Complex::default(); nw as usize]; nx as usize];
    let mut cg = vec![vec![Complex::default(); nw as usize]; npx as usize];
    let mut cfx = vec![Complex::default(); nx as usize];
    let mut cgpx = vec![Complex::default(); npx as usize];
    for ix in 0..nx as usize {
        for it in 0..nt as usize {
            fpad[it] = f[ix][it];
        }
        for it in nt as usize..ntfft as usize {
            fpad[it] = 0.0;
        }
        pfarc(1, &fpad, &mut cf[ix]);
    }
    let mut w = fw;
    for iw in 0..nw as usize {
        let scale = (-0.5 * w as f64 / (wmin * bwh * bwh) as f64) as f32;
        let amp = ((fftscl * dpx * w) as f64 / (2.0 * PI)
            * (w as f64 / (PI * wmin as f64)).sqrt()
            * dxb as f64
            / bwh as f64) as f32;
        let phase = w * (ft - ftau);
        let a = (amp as f64 * (phase as f64).cos()) as f32;
        let b = (amp as f64 * (phase as f64).sin()) as f32;
        let mut x = fx;
        for ix in 0..nx as usize {
            let es = ((scale * x * x) as f64).exp() as f32;
            let as_ = a * es;
            let bs = b * es;
            let cfr = cf[ix][iw].r;
            let cfi = cf[ix][iw].i;
            cfx[ix].r = as_ * cfr - bs * cfi;
            cfx[ix].i = bs * cfr + as_ * cfi;
            x += dx;
        }
        xtop(w, nx, dx, fx, &cfx, npx, dpx, fpx, &mut cgpx);
        for ipx in 0..npx as usize {
            cg[ipx][iw].r = cgpx[ipx].r;
            cg[ipx][iw].i = cgpx[ipx].i;
        }
        w += dw;
    }
    for ipx in 0..npx as usize {
        pfacr(-1, &cg[ipx], &mut gpad);
        for itau in 0..ntau as usize {
            g[ipx][itau] = gpad[itau];
        }
    }
}

fn xtop(
    w: f32,
    nx: i32,
    dx: f32,
    fx: f32,
    g: &[Complex],
    np: i32,
    dp: f32,
    fp: f32,
    h: &mut [Complex],
) {
    let lwrap = 8usize;
    let nxfft = npfa((nx + lwrap as i32) * 2) as usize;
    let nk = nxfft;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let fk = (-PI / dx as f64) as f32;
    let ek = (PI / dx as f64) as f32;
    let fka = fk - lwrap as f32 * dk;
    let nka = lwrap + nk + lwrap;
    let mut gka = vec![Complex::default(); nka];
    let mut hp = vec![Complex::default(); np as usize];
    let mut kp = vec![0.0f32; np as usize];
    for ix in 0..nx as usize {
        gka[ix + lwrap].r = dx * g[ix].r;
        gka[ix + lwrap].i = dx * g[ix].i;
    }
    for ix in nx as usize..nxfft {
        gka[ix + lwrap].r = 0.0;
        gka[ix + lwrap].i = 0.0;
    }
    for ix in (1..nx as usize).step_by(2) {
        gka[ix + lwrap].r = -gka[ix + lwrap].r;
        gka[ix + lwrap].i = -gka[ix + lwrap].i;
    }
    pfacc(-1, &mut gka[lwrap..lwrap + nxfft]);
    for ik in 0..lwrap {
        gka[ik] = gka[ik + nk];
    }
    for ik in lwrap + nk..lwrap + nk + lwrap {
        gka[ik] = gka[ik - nk];
    }
    let xshift = (0.5 * (nx - 1) as f64 * dx as f64) as f32;
    let mut k = fka;
    for ik in 0..nka {
        let phase = k * xshift;
        let c = (phase as f64).cos() as f32;
        let s = (phase as f64).sin() as f32;
        let temp = gka[ik].r * c - gka[ik].i * s;
        gka[ik].i = gka[ik].r * s + gka[ik].i * c;
        gka[ik].r = temp;
        k += dk;
    }
    let mut p = fp;
    for ip in 0..np as usize {
        kp[ip] = w * p;
        if kp[ip] < fk && kp[ip] < ek {
            kp[ip] = (fk as f64 - 1000.0 * dk as f64) as f32;
        } else if kp[ip] > fk && kp[ip] > ek {
            kp[ip] = (ek as f64 + 1000.0 * dk as f64) as f32;
        }
        p += dp;
    }
    ints8c(
        dk,
        fka,
        &gka,
        cmplx(0.0, 0.0),
        cmplx(0.0, 0.0),
        &kp,
        &mut hp,
    );
    let xshift = (-fx as f64 - 0.5 * (nx - 1) as f64 * dx as f64) as f32;
    for ip in 0..np as usize {
        let phase = kp[ip] * xshift;
        let c = (phase as f64).cos() as f32;
        let s = (phase as f64).sin() as f32;
        h[ip].r = hp[ip].r * c - hp[ip].i * s;
        h[ip].i = hp[ip].r * s + hp[ip].i * c;
    }
}

pub(super) fn beamData(wmin: f32, nt: i32, dt: f32, ft: f32, f: &[f32]) -> BeamData {
    let ntpad = 25;
    let ntfft = npfaro(nt + ntpad, 2 * (nt + ntpad));
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = 0.0f32;
    let iwnyq = nw - 1;
    let ntrfft = npfao(4 * ntfft, 4 * ntfft + ntfft);
    let nwr = ntrfft;
    let dtr = dt * ntfft as f32 / ntrfft as f32;
    let ftr = ft;
    let ntr = (1.0f32 + (nt + ntpad - 1) as f32 * dt / dtr) as i32;
    let nti = 6;
    let dti = (-5.0 / (wmin * (nti - 1) as f32) as f64) as f32;
    let fti = 0.0f32;
    let mut cf = vec![vec![Complex::default(); ntr as usize]; nti as usize];
    let mut fa = vec![0.0f32; ntfft as usize];
    let mut ca = vec![Complex::default(); nw as usize];
    let mut cb = vec![Complex::default(); ntrfft as usize];
    for it in 0..nt as usize {
        fa[it] = f[it];
    }
    for it in nt as usize..ntfft as usize {
        fa[it] = 0.0;
    }
    pfarc(1, &fa, &mut ca);
    for iw in 1..iwnyq as usize {
        ca[iw].r = (ca[iw].r as f64 * 2.0) as f32;
        ca[iw].i = (ca[iw].i as f64 * 2.0) as f32;
    }
    let mut ti = fti;
    for iti in 0..nti as usize {
        let mut w = fw;
        for iw in 0..nw as usize {
            let scale = ((w * ti) as f64).exp() as f32;
            cb[iw].r = ca[iw].r * scale;
            cb[iw].i = ca[iw].i * scale;
            w += dw;
        }
        for iw in nw as usize..nwr as usize {
            cb[iw].r = 0.0;
            cb[iw].i = 0.0;
        }
        pfacc(-1, &mut cb);
        let scale = (1.0 / ntfft as f64) as f32;
        for itr in 0..ntr as usize {
            cf[iti][itr].r = scale * cb[itr].r;
            cf[iti][itr].i = scale * cb[itr].i;
        }
        ti += dti;
    }
    BeamData {
        ntr,
        dtr,
        ftr,
        nti,
        dti,
        fti,
        cf,
    }
}
