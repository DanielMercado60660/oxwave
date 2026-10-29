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

//! Ports of cwp/lib/intsinc8.c, inttable8.c, mksinc.c, ressinc8.c, and shfs8r.c.
//! Credits: Dave Hale, Colorado School of Mines, 06/02/89;
//! ressinc8.c: Dave Hale, Colorado School of Mines, 06/06/90.
//!
//! Array lengths replace the corresponding C counts. Interpolation writes
//! `xout.len()` values; resampling and shifting write `yout.len()` values.
//! Interpolation leaves any remaining output elements unchanged.
//! Input intervals must be nonzero, and coordinates must produce finite,
//! representable C `int` sample indices. Slice lengths must fit in C `int`.

use oxwave_cwp::complex::Complex;
use oxwave_cwp::cwp::nint;
use oxwave_sinc::dsinc;
use oxwave_toeplitz::stoepd;
use std::sync::OnceLock;

/// Compute least-squares optimal sinc interpolation coefficients.
///
/// C prototype: `void mksinc(float d, int lsinc, float sinc[]);`
///
/// `sinc.len()` is `lsinc`, normally even and at most 20; `d` is normally
/// in [0, 1]. Odd lengths and distances outside that range are evaluated by
/// the same C formula. Panics unless the length is in 1..=20.
pub fn mksinc(d: f32, sinc: &mut [f32]) {
    let lsinc = sinc.len();
    assert!((1..=20).contains(&lsinc));
    let mut s = [0.0f64; 20];
    let mut a = [0.0f64; 20];
    let mut c = [0.0f64; 20];
    let mut work = [0.0f64; 20];
    let fmax = 0.066 + 0.265 * (lsinc as f64).ln();
    let fmax = if fmax < 1.0 { fmax } else { 1.0 };
    for j in 0..lsinc {
        a[j] = dsinc(fmax * j as f64);
        // C rounds the integer-plus-float expression before multiplying
        // it by the double fmax.
        c[j] = dsinc(fmax * ((lsinc as i32 / 2 - j as i32 - 1) as f32 + d) as f64);
    }
    stoepd(
        &a[..lsinc],
        &c[..lsinc],
        &mut s[..lsinc],
        &mut work[..lsinc],
    );
    for j in 0..lsinc {
        sinc[j] = s[j] as f32;
    }
}

// The C callers initialize identical immutable tables. OnceLock shares
// those coefficients safely while preserving their computation order.
fn sinc_table() -> &'static [[f32; 8]; 513] {
    static TABLE: OnceLock<[[f32; 8]; 513]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [[0.0; 8]; 513];
        for jtable in 1..512 {
            let frac = jtable as f32 / 512.0f32;
            mksinc(frac, &mut table[jtable]);
        }
        for jtable in 0..8 {
            table[0][jtable] = 0.0;
            table[512][jtable] = 0.0;
        }
        table[0][3] = 1.0;
        table[512][4] = 1.0;
        table
    })
}

/// Interpolate uniformly sampled complex data with eight sinc coefficients.
///
/// C prototype: `void ints8c(int nxin, float dxin, float fxin, complex yin[],
/// complex yinl, complex yinr, int nxout, float xout[], complex yout[]);`
///
/// `yinl` and `yinr` extrapolate beyond the input. Panics if `yout` is short.
pub fn ints8c(
    dxin: f32,
    fxin: f32,
    yin: &[Complex],
    yinl: Complex,
    yinr: Complex,
    xout: &[f32],
    yout: &mut [Complex],
) {
    intt8c(sinc_table(), dxin, fxin, yin, yinl, yinr, xout, yout);
}

/// Interpolate uniformly sampled real data with eight sinc coefficients.
///
/// C prototype: `void ints8r(int nxin, float dxin, float fxin, float yin[],
/// float yinl, float yinr, int nxout, float xout[], float yout[]);`
///
/// `yinl` and `yinr` extrapolate beyond the input. Panics if `yout` is short.
pub fn ints8r(
    dxin: f32,
    fxin: f32,
    yin: &[f32],
    yinl: f32,
    yinr: f32,
    xout: &[f32],
    yout: &mut [f32],
) {
    intt8r(sinc_table(), dxin, fxin, yin, yinl, yinr, xout, yout);
}

/// Interpolate complex samples using the nearest row of an eight-point table.
///
/// C prototype: `void intt8c(int ntable, float table[][8], int nxin, float dxin,
/// float fxin, complex yin[], complex yinl, complex yinr, int nxout,
/// float xout[], complex yout[]);`
///
/// Rows represent fractions `i/(ntable-1)`; endpoint rows normally have unit
/// weights at positions 3 and 4. Panics if the table has fewer than two rows,
/// or if `yout` is short.
pub fn intt8c(
    table: &[[f32; 8]],
    dxin: f32,
    fxin: f32,
    yin: &[Complex],
    yinl: Complex,
    yinr: Complex,
    xout: &[f32],
    yout: &mut [Complex],
) {
    assert!(table.len() >= 2 && yout.len() >= xout.len());
    let ioutb = -3i32 - 8;
    let xoutf = fxin;
    let xouts = (1.0f64 / dxin as f64) as f32;
    let xoutb = (8.0f64 - (xoutf * xouts) as f64) as f32;
    let fntablem1 = (table.len() - 1) as f32;
    let nxinm8 = yin.len() as i32 - 8;
    for ixout in 0..xout.len() {
        let xoutn = xoutb + xout[ixout] * xouts;
        let ixoutn = xoutn as i32;
        let mut kyin = ioutb.wrapping_add(ixoutn);
        let frac = xoutn - ixoutn as f32;
        let ktable = if frac >= 0.0 {
            ((frac * fntablem1) as f64 + 0.5) as i32
        } else {
            ((frac as f64 + 1.0) * fntablem1 as f64 - 0.5) as i32
        };
        let ptable = &table[ktable as usize];
        if kyin >= 0 && kyin <= nxinm8 {
            let pyin = &yin[kyin as usize..];
            yout[ixout].r = pyin[0].r * ptable[0]
                + pyin[1].r * ptable[1]
                + pyin[2].r * ptable[2]
                + pyin[3].r * ptable[3]
                + pyin[4].r * ptable[4]
                + pyin[5].r * ptable[5]
                + pyin[6].r * ptable[6]
                + pyin[7].r * ptable[7];
            yout[ixout].i = pyin[0].i * ptable[0]
                + pyin[1].i * ptable[1]
                + pyin[2].i * ptable[2]
                + pyin[3].i * ptable[3]
                + pyin[4].i * ptable[4]
                + pyin[5].i * ptable[5]
                + pyin[6].i * ptable[6]
                + pyin[7].i * ptable[7];
        } else {
            let mut sumr = 0.0;
            let mut sumi = 0.0;
            for itable in 0..8 {
                let yini = if kyin < 0 {
                    yinl
                } else if kyin >= yin.len() as i32 {
                    yinr
                } else {
                    yin[kyin as usize]
                };
                sumr += yini.r * ptable[itable];
                sumi += yini.i * ptable[itable];
                kyin = kyin.wrapping_add(1);
            }
            yout[ixout].r = sumr;
            yout[ixout].i = sumi;
        }
    }
}

/// Interpolate real samples using the nearest row of an eight-point table.
///
/// C prototype: `void intt8r(int ntable, float table[][8], int nxin, float dxin,
/// float fxin, float yin[], float yinl, float yinr, int nxout,
/// float xout[], float yout[]);`
///
/// Rows represent fractions `i/(ntable-1)`; endpoint rows normally have unit
/// weights at positions 3 and 4. Panics if the table has fewer than two rows,
/// or if `yout` is short.
pub fn intt8r(
    table: &[[f32; 8]],
    dxin: f32,
    fxin: f32,
    yin: &[f32],
    yinl: f32,
    yinr: f32,
    xout: &[f32],
    yout: &mut [f32],
) {
    assert!(table.len() >= 2 && yout.len() >= xout.len());
    let ioutb = -3i32 - 8;
    let xoutf = fxin;
    let xouts = (1.0f64 / dxin as f64) as f32;
    let xoutb = (8.0f64 - (xoutf * xouts) as f64) as f32;
    let fntablem1 = (table.len() - 1) as f32;
    let nxinm8 = yin.len() as i32 - 8;
    for ixout in 0..xout.len() {
        let xoutn = xoutb + xout[ixout] * xouts;
        let ixoutn = xoutn as i32;
        let mut kyin = ioutb.wrapping_add(ixoutn);
        let frac = xoutn - ixoutn as f32;
        let ktable = if frac >= 0.0 {
            ((frac * fntablem1) as f64 + 0.5) as i32
        } else {
            ((frac as f64 + 1.0) * fntablem1 as f64 - 0.5) as i32
        };
        let ptable = &table[ktable as usize];
        if kyin >= 0 && kyin <= nxinm8 {
            let pyin = &yin[kyin as usize..];
            yout[ixout] = pyin[0] * ptable[0]
                + pyin[1] * ptable[1]
                + pyin[2] * ptable[2]
                + pyin[3] * ptable[3]
                + pyin[4] * ptable[4]
                + pyin[5] * ptable[5]
                + pyin[6] * ptable[6]
                + pyin[7] * ptable[7];
        } else {
            let mut sum = 0.0;
            for itable in 0..8 {
                let yini = if kyin < 0 {
                    yinl
                } else if kyin >= yin.len() as i32 {
                    yinr
                } else {
                    yin[kyin as usize]
                };
                sum += yini * ptable[itable];
                kyin = kyin.wrapping_add(1);
            }
            yout[ixout] = sum;
        }
    }
}

/// Resample complex data, first checking whether all samples lie on the grid.
///
/// C prototype: `void ress8c(int nxin, float dxin, float fxin, complex yin[],
/// complex yinl, complex yinr, int nxout, float dxout, float fxout,
/// complex yout[]);`
///
/// Distances within 0.01 input samples use direct copies if every output
/// qualifies. Otherwise the entire output is recomputed with `ints8c`.
pub fn ress8c(
    dxin: f32,
    fxin: f32,
    yin: &[Complex],
    yinl: Complex,
    yinr: Complex,
    dxout: f32,
    fxout: f32,
    yout: &mut [Complex],
) {
    let odxin = (1.0f64 / dxin as f64) as f32;
    let tiny = 0.01f32;
    let mut xout = fxout;
    let mut ixout = 0;
    while ixout < yout.len() {
        let xini = (xout - fxin) * odxin;
        let ixin = nint(xini as f64);
        let error = xini - ixin as f32;
        if (if error < 0.0 { -error } else { error }) > tiny {
            break;
        }
        yout[ixout] = if ixin < 0 {
            yinl
        } else if ixin >= yin.len() as i32 {
            yinr
        } else {
            yin[ixin as usize]
        };
        ixout += 1;
        xout += dxout;
    }
    if ixout == yout.len() {
        return;
    }
    let mut xoutt = vec![0.0; yout.len()];
    for ixout in 0..yout.len() {
        xoutt[ixout] = fxout + ixout as f32 * dxout;
    }
    ints8c(dxin, fxin, yin, yinl, yinr, &xoutt, yout);
}

/// Resample real data, first checking whether all samples lie on the grid.
///
/// C prototype: `void ress8r(int nxin, float dxin, float fxin, float yin[],
/// float yinl, float yinr, int nxout, float dxout, float fxout, float yout[]);`
///
/// Distances within 0.01 input samples use direct copies if every output
/// qualifies. Otherwise the entire output is recomputed with `ints8r`.
pub fn ress8r(
    dxin: f32,
    fxin: f32,
    yin: &[f32],
    yinl: f32,
    yinr: f32,
    dxout: f32,
    fxout: f32,
    yout: &mut [f32],
) {
    let odxin = (1.0f64 / dxin as f64) as f32;
    let tiny = 0.01f32;
    let mut xout = fxout;
    let mut ixout = 0;
    while ixout < yout.len() {
        let xini = (xout - fxin) * odxin;
        let ixin = nint(xini as f64);
        let error = xini - ixin as f32;
        if (if error < 0.0 { -error } else { error }) > tiny {
            break;
        }
        yout[ixout] = if ixin < 0 {
            yinl
        } else if ixin >= yin.len() as i32 {
            yinr
        } else {
            yin[ixin as usize]
        };
        ixout += 1;
        xout += dxout;
    }
    if ixout == yout.len() {
        return;
    }
    let mut xoutt = vec![0.0; yout.len()];
    for ixout in 0..yout.len() {
        xoutt[ixout] = fxout + ixout as f32 * dxout;
    }
    ints8r(dxin, fxin, yin, yinl, yinr, &xoutt, yout);
}

/// Shift real samples with the same interval using eight sinc coefficients.
///
/// C prototype: `void shfs8r(float dx, int nxin, float fxin, float yin[],
/// float yinl, float yinr, int nxout, float fxout, float yout[]);`
pub fn shfs8r(dx: f32, fxin: f32, yin: &[f32], yinl: f32, yinr: f32, fxout: f32, yout: &mut [f32]) {
    let table = sinc_table();
    let shift = (fxout - fxin) / dx;
    let mut ishift = shift as i32;
    if shift < 0.0 {
        ishift = ishift.wrapping_sub(1);
    }
    let dshift = shift - ishift as f32;
    let ktable = nint((dshift * 512.0f32) as f64);
    let nxout = yout.len() as i32;
    let nxin = yin.len() as i32;
    for ixout in 0..yout.len() {
        yout[ixout] = 0.0;
    }
    for itable in 0..8 {
        let tablei = table[ktable as usize][itable];
        let jshift = ishift
            .wrapping_add(itable as i32)
            .wrapping_sub(4)
            .wrapping_add(1);
        let ihi = nxout.min(jshift.wrapping_neg());
        let youti = yinl * tablei;
        for ixout in 0..ihi {
            yout[ixout as usize] += youti;
        }
        let ilo = 0.max(jshift.wrapping_neg());
        let ihi = nxout.min(nxin.wrapping_sub(jshift));
        for ixout in ilo..ihi {
            yout[ixout as usize] += yin[ixout.wrapping_add(jshift) as usize] * tablei;
        }
        let ilo = 0.max(nxin.wrapping_sub(jshift));
        let youti = yinr * tablei;
        for ixout in ilo..nxout {
            yout[ixout as usize] += youti;
        }
    }
}
