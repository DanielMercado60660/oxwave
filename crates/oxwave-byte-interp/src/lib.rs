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

//! Bilinear byte interpolation from cwp/lib/intl2b.c and cwp/lib/intlirr2b.c.
//! Author: Dave Hale, Colorado School of Mines, c. 1989-1991.
//! Author: Dave Hale, Colorado School of Mines, 07/02/89.
//! Modified: Dave Hale, Colorado School of Mines, 05/30/90: interpolate
//! unsigned rather than signed char (range 0 to 255).
//! Modified: Dave Hale, Colorado School of Mines, 06/01/91: avoid truncation
//! errors in the coefficient table that could exceed the input byte range.

use std::sync::OnceLock;

const ICMAX: i32 = 99;
const NTABLE: usize = (ICMAX + 1) as usize;
type Table = [[u8; 256]; NTABLE];

fn table() -> &'static Table {
    // Both C files build the same table. OnceLock preserves its exact contents
    // while making initialization safe for concurrent callers.
    static TABLE: OnceLock<Table> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [[0; 256]; NTABLE];
        for ic in 0..=ICMAX / 2 {
            let frac = ic as f32 / ICMAX as f32;
            for ib in 0..256 {
                table[ic as usize][ib] = (frac * ib as f32) as u8;
                table[(ICMAX - ic) as usize][ib] = ib as u8 - table[ic as usize][ib];
            }
        }
        table
    })
}

/// Bilinear interpolation of a contiguous row-major byte grid, with constant
/// extrapolation. Dimensions are retained to describe the flat slices.
///
/// C prototype:
/// ```c
/// void intl2b(int nxin, float dxin, float fxin,
///     int nyin, float dyin, float fyin, unsigned char *zin,
///     int nxout, float dxout, float fxout,
///     int nyout, float dyout, float fyout, unsigned char *zout);
/// ```
///
/// Input dimensions must be positive; output dimensions must be nonnegative.
/// The slices must contain at least the corresponding dimension products.
/// Negative input/output intervals and zero output intervals are supported.
/// Singleton input dimensions use constant interpolation: C's out-of-grid
/// reads with zero weight are omitted. No guard bytes are needed in Rust.
#[allow(clippy::too_many_arguments)]
pub fn intl2b(
    nxin: i32,
    dxin: f32,
    fxin: f32,
    nyin: i32,
    dyin: f32,
    fyin: f32,
    zin: &[u8],
    nxout: i32,
    dxout: f32,
    fxout: f32,
    nyout: i32,
    dyout: f32,
    fyout: f32,
    zout: &mut [u8],
) {
    check_dimensions(nxin, nyin, zin, nxout, nyout, zout);
    let table = table();
    let mut kzin = vec![0i32; nxout as usize];
    let mut kic = vec![0i32; nxout as usize];

    // Increment the float coordinate as C does, rather than multiplying ixout.
    let mut xout = fxout;
    for ixout in 0..nxout as usize {
        let rxin = (xout - fxin) / dxin;
        if rxin <= 0.0 {
            kzin[ixout] = 0;
            kic[ixout] = 0;
        } else if rxin >= (nxin - 1) as f32 {
            kzin[ixout] = nxin - 2;
            kic[ixout] = ICMAX * 256;
        } else {
            kzin[ixout] = rxin as i32;
            let frac = rxin - (rxin as i32) as f32;
            let ic = ((frac * ICMAX as f32) as f64 + 0.5) as i32;
            kic[ixout] = ic * 256;
        }
        xout += dxout;
    }
    interpolate_y(
        nxin, nyin, dyin, fyin, zin, nxout, nyout, dyout, fyout, &kzin, &kic, table, zout,
    );
}

/// Bilinear byte interpolation with irregular, strictly increasing input x
/// coordinates. `xin.len()` supplies C's `nxin`; other dimensions describe
/// the flat row-major byte slices. Extrapolation is constant.
///
/// C prototype:
/// ```c
/// void intlirr2b(int nxin, float *xin,
///     int nyin, float dyin, float fyin, unsigned char *zin,
///     int nxout, float dxout, float fxout,
///     int nyout, float dyout, float fyout, unsigned char *zout);
/// ```
///
/// `xin` must be nonempty and its length must fit a C int. Other dimension
/// and slice requirements, including singleton handling, are as in [`intl2b`].
#[allow(clippy::too_many_arguments)]
pub fn intlirr2b(
    xin: &[f32],
    nyin: i32,
    dyin: f32,
    fyin: f32,
    zin: &[u8],
    nxout: i32,
    dxout: f32,
    fxout: f32,
    nyout: i32,
    dyout: f32,
    fyout: f32,
    zout: &mut [u8],
) {
    assert!(xin.len() <= i32::MAX as usize);
    let nxin = xin.len() as i32;
    check_dimensions(nxin, nyin, zin, nxout, nyout, zout);
    assert!(xin.iter().all(|x| x.is_finite()));
    assert!(xin.windows(2).all(|w| w[0] < w[1]));
    let table = table();
    let mut kzin = vec![0i32; nxout as usize];
    let mut kic = vec![0i32; nxout as usize];
    let mut frac = 0.0f32;

    let mut xout = fxout;
    for ixout in 0..nxout as usize {
        // A NaN coordinate would make C's binary search fail to terminate.
        assert!(!xout.is_nan());
        if xout <= xin[0] {
            kzin[ixout] = 0;
            kic[ixout] = 0;
        } else if xout >= xin[(nxin - 1) as usize] {
            kzin[ixout] = nxin - 2;
            kic[ixout] = ICMAX * 256;
        } else {
            let mut found = false;
            let mut il = 0;
            let mut im = nxin / 2;
            let mut ir = nxin - 1;
            while !found {
                if xout >= xin[il as usize] && xout <= xin[im as usize] {
                    if im - il == 1 {
                        kzin[ixout] = il;
                        frac = (xout - xin[il as usize]) / (xin[im as usize] - xin[il as usize]);
                        found = true;
                    } else {
                        ir = im;
                        im = (il + im) / 2;
                    }
                } else if ir - im == 1 {
                    kzin[ixout] = im;
                    frac = (xout - xin[im as usize]) / (xin[ir as usize] - xin[im as usize]);
                    found = true;
                } else {
                    il = im;
                    im = (im + ir) / 2;
                }
                let ic = ((frac * ICMAX as f32) as f64 + 0.5) as i32;
                kic[ixout] = ic * 256;
            }
        }
        xout += dxout;
    }
    interpolate_y(
        nxin, nyin, dyin, fyin, zin, nxout, nyout, dyout, fyout, &kzin, &kic, table, zout,
    );
}

fn check_dimensions(nxin: i32, nyin: i32, zin: &[u8], nxout: i32, nyout: i32, zout: &[u8]) {
    assert!(nxin > 0 && nyin > 0 && nxout >= 0 && nyout >= 0);
    assert!(zin.len() >= nxin as usize * nyin as usize);
    assert!(zout.len() >= nxout as usize * nyout as usize);
}

// Identical y interpolation and temporary-row reuse in the two C files.
#[allow(clippy::too_many_arguments)]
fn interpolate_y(
    nxin: i32,
    nyin: i32,
    dyin: f32,
    fyin: f32,
    zin: &[u8],
    nxout: i32,
    nyout: i32,
    dyout: f32,
    fyout: f32,
    kzin: &[i32],
    kic: &[i32],
    table: &Table,
    zout: &mut [u8],
) {
    let mut temp1 = vec![0u8; nxout as usize];
    let mut temp2 = vec![0u8; nxout as usize];
    let mut iyinl = 1;
    let mut yout = fyout;
    for iyout in 0..nyout as usize {
        // Translate MAX/MIN as comparisons: f32::min/max differ for NaN.
        let ratio = (yout - fyin) / dyin;
        let upper = if ((nyin - 1) as f32) < ratio {
            (nyin - 1) as f32
        } else {
            ratio
        };
        let ryin = if 0.0 > upper { 0.0 } else { upper };
        let upper = if ((nyin - 2) as f32) < ryin {
            (nyin - 2) as f32
        } else {
            ryin
        };
        let iyin = (if 0.0 > upper { 0.0 } else { upper }) as i32;

        if iyin != iyinl || iyout == 0 {
            if nyin == 1 {
                // C also computes a nonexistent second row, subsequently
                // multiplied by table[0] in intl2by. Its bytes cannot affect
                // the result. Leave temp2 zero instead of reading that row.
                intl2bx(kzin, kic, nxin, table, zin, &mut temp1);
            } else if iyin == iyinl + 1 && iyout != 0 {
                std::mem::swap(&mut temp1, &mut temp2);
                intl2bx(
                    kzin,
                    kic,
                    nxin,
                    table,
                    &zin[((iyin + 1) as usize * nxin as usize)..],
                    &mut temp2,
                );
            } else if iyin == iyinl - 1 && iyout != 0 {
                std::mem::swap(&mut temp1, &mut temp2);
                intl2bx(
                    kzin,
                    kic,
                    nxin,
                    table,
                    &zin[(iyin as usize * nxin as usize)..],
                    &mut temp1,
                );
            } else {
                intl2bx(
                    kzin,
                    kic,
                    nxin,
                    table,
                    &zin[(iyin as usize * nxin as usize)..],
                    &mut temp1,
                );
                intl2bx(
                    kzin,
                    kic,
                    nxin,
                    table,
                    &zin[((iyin + 1) as usize * nxin as usize)..],
                    &mut temp2,
                );
            }
            iyinl = iyin;
        }
        let frac = ryin - iyin as f32;
        let ic = ((frac * ICMAX as f32) as f64 + 0.5) as i32;
        intl2by(
            ic,
            table,
            &temp1,
            &temp2,
            &mut zout[iyout * nxout as usize..],
        );
        yout += dyout;
    }
}

fn intl2bx(kzin: &[i32], kic: &[i32], nxin: i32, table: &Table, zin: &[u8], zout: &mut [u8]) {
    for ixout in 0..kzin.len() {
        if nxin == 1 {
            // kzin is either 0 (coefficient 0) or -1 (coefficient ICMAX).
            // The out-of-grid neighbor has zero weight in both cases.
            zout[ixout] = zin[0];
        } else {
            let jzin = kzin[ixout] as usize;
            let jic = (kic[ixout] / 256) as usize;
            zout[ixout] = (table[ICMAX as usize - jic][zin[jzin] as usize] as i32
                + table[jic][zin[jzin + 1] as usize] as i32) as u8;
        }
    }
}

fn intl2by(ic: i32, table: &Table, temp1: &[u8], temp2: &[u8], zout: &mut [u8]) {
    for ixout in 0..temp1.len() {
        zout[ixout] = (table[(ICMAX - ic) as usize][temp1[ixout] as usize] as i32
            + table[ic as usize][temp2[ixout] as usize] as i32) as u8;
    }
}
