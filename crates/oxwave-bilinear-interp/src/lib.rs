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

//! Port of su/lib/bilinear.c. Credits: none listed in the C source.
use oxwave_linear_interp::{bhighi, linterpd};

/// Find upper grid indices and lower-point weights on nonempty sorted axes.
/// C: `void binterpfind(int kigi, int *mgi, int mgi_tot, int mgiextr, int kigc, int *mgc, int mgc_tot, int mgcextr, int *mgixo, int *mgcxo, double *wi, double *wc)`.
pub fn binterpfind(
    mut kigi: i32,
    mgi: &[i32],
    mgiextr: i32,
    mut kigc: i32,
    mgc: &[i32],
    mgcextr: i32,
    mgixo: &mut i32,
    mgcxo: &mut i32,
    wi: &mut f64,
    wc: &mut f64,
) {
    assert!(
        !mgi.is_empty()
            && !mgc.is_empty()
            && mgi.len() <= i32::MAX as usize
            && mgc.len() <= i32::MAX as usize
    );
    if mgi.len() == 1 && mgc.len() == 1 {
        *wi = 1.0;
        *wc = 0.0;
        *mgixo = 1;
        *mgcxo = 1;
        return;
    }
    let mgix = if kigi <= mgi[0] {
        if mgiextr == 0 || mgiextr == 3 {
            kigi = mgi[0];
        }
        1
    } else if kigi >= mgi[mgi.len() - 1] {
        if mgiextr == 0 || mgiextr == 2 {
            kigi = mgi[mgi.len() - 1];
        }
        mgi.len() - 1
    } else {
        bhighi(mgi, kigi) as usize
    };
    if mgc.len() == 1 {
        *wi = mgi[mgix].wrapping_sub(kigi) as f64 / mgi[mgix].wrapping_sub(mgi[mgix - 1]) as f64;
        *wc = 0.0;
        *mgixo = mgix as i32;
        *mgcxo = 1;
        return;
    }
    let mgcx = if kigc <= mgc[0] {
        if mgcextr == 0 || mgcextr == 3 {
            kigc = mgc[0];
        }
        1
    } else if kigc >= mgc[mgc.len() - 1] {
        if mgcextr == 0 || mgcextr == 2 {
            kigc = mgc[mgc.len() - 1];
        }
        mgc.len() - 1
    } else {
        bhighi(mgc, kigc) as usize
    };
    if mgi.len() == 1 {
        *wi = 0.0;
        *wc = mgc[mgcx].wrapping_sub(kigc) as f64 / mgc[mgcx].wrapping_sub(mgc[mgcx - 1]) as f64;
        *mgixo = 1;
        *mgcxo = mgcx as i32;
        return;
    }
    *wi = mgi[mgix].wrapping_sub(kigi) as f64 / mgi[mgix].wrapping_sub(mgi[mgix - 1]) as f64;
    *wc = mgc[mgcx].wrapping_sub(kigc) as f64 / mgc[mgcx].wrapping_sub(mgc[mgcx - 1]) as f64;
    *mgixo = mgix as i32;
    *mgcxo = mgcx as i32;
}

/// Apply grid weights to input arrays; output length controls the number processed.
/// C: `void binterpapply(double *lwitims, double *hiitims, int mgi_tot, double wi, double *lwctims, double *hictims, int mgc_tot, double wc, int lwinto, double *valsout)`.
pub fn binterpapply(
    lwitims: &[f64],
    hiitims: &[f64],
    mgi_tot: i32,
    wi: f64,
    lwctims: &[f64],
    hictims: &[f64],
    mgc_tot: i32,
    wc: f64,
    valsout: &mut [f64],
) {
    let len = valsout.len();
    if mgi_tot == 1 && mgc_tot == 1 {
        for n in 0..len {
            valsout[n] = lwitims[n];
        }
        return;
    }
    for n in 0..len {
        valsout[n] = 0.0;
    }
    if mgc_tot == 1 {
        if wi != 0.0 {
            for n in 0..len {
                valsout[n] += wi * lwitims[n];
            }
        }
        if wi != 1.0 {
            for n in 0..len {
                valsout[n] += (1.0 - wi) * hiitims[n];
            }
        }
        return;
    }
    if mgi_tot == 1 {
        if wc != 0.0 {
            for n in 0..len {
                valsout[n] += wc * lwctims[n];
            }
        }
        if wc != 1.0 {
            for n in 0..len {
                valsout[n] += (1.0 - wc) * hictims[n];
            }
        }
        return;
    }
    let aw = wc * wi;
    let bw = wc * (1.0 - wi);
    let cw = (1.0 - wc) * wi;
    let dw = (1.0 - wc) * (1.0 - wi);
    if aw != 0.0 {
        for n in 0..len {
            valsout[n] += aw * lwitims[n];
        }
    }
    if bw != 0.0 {
        for n in 0..len {
            valsout[n] += bw * hiitims[n];
        }
    }
    if cw != 0.0 {
        for n in 0..len {
            valsout[n] += cw * lwctims[n];
        }
    }
    if dw != 0.0 {
        for n in 0..len {
            valsout[n] += dw * hictims[n];
        }
    }
}

/// Interpolate each contributing offset table and apply the grid weights.
/// C: `void binterpvalue(double offset, int mgtextr, double *lwioffs, double *lwitims, int lwinto, double *hiioffs, double *hiitims, int hiinto, int mgi_tot, double wi, double *lwcoffs, double *lwctims, int lwcnto, double *hicoffs, double *hictims, int hicnto, int mgc_tot, double wc, double *timeout)`.
pub fn binterpvalue(
    offset: f64,
    mgtextr: i32,
    lwioffs: &[f64],
    lwitims: &[f64],
    hiioffs: &[f64],
    hiitims: &[f64],
    mgi_tot: i32,
    wi: f64,
    lwcoffs: &[f64],
    lwctims: &[f64],
    hicoffs: &[f64],
    hictims: &[f64],
    mgc_tot: i32,
    wc: f64,
    timeout: &mut f64,
) {
    if mgi_tot == 1 && mgc_tot == 1 {
        linterpd(offset, lwioffs, lwitims, mgtextr, timeout);
        return;
    }
    *timeout = 0.0;
    let mut time = 0.0;
    if mgc_tot == 1 {
        if wi != 0.0 {
            linterpd(offset, lwioffs, lwitims, mgtextr, &mut time);
            *timeout += wi * time;
        }
        if wi != 1.0 {
            linterpd(offset, hiioffs, hiitims, mgtextr, &mut time);
            *timeout += (1.0 - wi) * time;
        }
        return;
    }
    if mgi_tot == 1 {
        if wc != 0.0 {
            linterpd(offset, lwcoffs, lwctims, mgtextr, &mut time);
            *timeout += wc * time;
        }
        if wc != 1.0 {
            linterpd(offset, hicoffs, hictims, mgtextr, &mut time);
            *timeout += (1.0 - wc) * time;
        }
        return;
    }
    let aw = wc * wi;
    let bw = wc * (1.0 - wi);
    let cw = (1.0 - wc) * wi;
    let dw = (1.0 - wc) * (1.0 - wi);
    if aw != 0.0 {
        linterpd(offset, lwioffs, lwitims, mgtextr, &mut time);
        *timeout += aw * time;
    }
    if bw != 0.0 {
        linterpd(offset, hiioffs, hiitims, mgtextr, &mut time);
        *timeout += bw * time;
    }
    if cw != 0.0 {
        linterpd(offset, lwcoffs, lwctims, mgtextr, &mut time);
        *timeout += cw * time;
    }
    if dw != 0.0 {
        linterpd(offset, hicoffs, hictims, mgtextr, &mut time);
        *timeout += dw * time;
    }
}
