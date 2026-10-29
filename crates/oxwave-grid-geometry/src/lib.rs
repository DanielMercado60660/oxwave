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

//! Grid geometry ported from `su/lib/gridxy.c`.
//! Credits: the C source names no individual author; `su/include/gridxy.h`
//! carries "Copyright (c) Colorado School of Mines, 2021. All rights reserved."
//!
//! Grid definitions use 18 doubles: 0 is the unrelated bin type, 1 is the
//! left/right multiplier, 2..10 hold corners A/B/C/D, 10..12 cell widths,
//! 12..14 cell counts, 14..16 first/last CDP, and 16..18 sine/cosine.
//! Supply corners A/B/C and widths to `gridset` before calling conversions.
//! Array pointers become slices and scalar outputs become mutable references.
//! Out-of-grid CDPs/indexes use -2147483645. Index-to-coordinate conversions
//! also accept outside indexes.

use oxwave_core::warn;

/// Prepare a grid, snapping its corners to cell centres.
/// Status: 0 success, 1/2 nonpositive B/C width, 3 B too close to A,
/// -1 usable grid with C reset to A (one row). Requires 18 elements.
/// C prototype: `void gridset(double *gvals, int *errwarn);`
pub fn gridset(gvals: &mut [f64], errwarn: &mut i32) {
    *errwarn = 0;
    if gvals[10] <= 0.0 {
        *errwarn = 1;
        return;
    }
    if gvals[11] <= 0.0 {
        *errwarn = 2;
        return;
    }
    let dab = ((gvals[2] - gvals[4]) * (gvals[2] - gvals[4])
        + (gvals[3] - gvals[5]) * (gvals[3] - gvals[5]))
        .sqrt();
    let mut nwb = (dab / gvals[10] + 0.5) as i32;
    let dabwb = nwb as f64 * gvals[10];
    if nwb < 1 {
        *errwarn = 3;
        return;
    }
    nwb = nwb.wrapping_add(1);
    gvals[4] = gvals[2] + dabwb / dab * (gvals[4] - gvals[2]);
    gvals[5] = gvals[3] + dabwb / dab * (gvals[5] - gvals[3]);
    let dac = ((gvals[2] - gvals[6]) * (gvals[2] - gvals[6])
        + (gvals[3] - gvals[7]) * (gvals[3] - gvals[7]))
        .sqrt();
    let mut nwc = (dac / gvals[11] + 0.5) as i32;
    if nwc < 1 {
        *errwarn = -1;
    }
    nwc = nwc.wrapping_add(1);
    gvals[12] = nwb as f64;
    gvals[13] = nwc as f64;
    gvals[14] = 1.0;
    // nwb*nwc is evaluated as a C int before being promoted to double.
    gvals[15] = gvals[14] + nwb.wrapping_mul(nwc) as f64 - 1.0;
    gvals[16] = (gvals[5] - gvals[3]) / dabwb;
    gvals[17] = (gvals[4] - gvals[2]) / dabwb;
    let dx = gvals[6] - gvals[2];
    let dy = gvals[7] - gvals[3];
    let ty = -dx * gvals[16] + dy * gvals[17];
    gvals[1] = 1.0;
    if ty < 0.0 {
        gvals[1] = -1.0;
    }
    let (mut rx, mut ry) = (0.0, 0.0);
    gridicrawxy(gvals, nwb, 1, &mut rx, &mut ry);
    gvals[4] = rx;
    gvals[5] = ry;
    gridicrawxy(gvals, 1, nwc, &mut rx, &mut ry);
    gvals[6] = rx;
    gvals[7] = ry;
    gridicrawxy(gvals, nwb, nwc, &mut rx, &mut ry);
    gvals[8] = rx;
    gvals[9] = ry;
}

/// Convert world coordinates to a CDP and grid indexes (retained outside grid).
/// C prototype: `void gridrawxycdpic(double *gvals, double dx, double dy, int *icdp, int *igi, int *igc);`
pub fn gridrawxycdpic(
    gvals: &[f64],
    mut dx: f64,
    mut dy: f64,
    icdp: &mut i32,
    igi: &mut i32,
    igc: &mut i32,
) {
    dx = dx - gvals[2];
    dy = dy - gvals[3];
    let rx = dx * gvals[17] + dy * gvals[16];
    let ry = (-dx * gvals[16] + dy * gvals[17]) * gvals[1];
    *igi = (rx / gvals[10] + 1.5).floor() as i32;
    *igc = (ry / gvals[11] + 1.5).floor() as i32;
    if *igi < 1 || *igi as f64 > gvals[12] || *igc < 1 || *igc as f64 > gvals[13] {
        *icdp = -2147483645;
    } else {
        *icdp = (gvals[14] + *igi as f64 - 1.0 + igc.wrapping_sub(1) as f64 * gvals[12]) as i32;
    }
}

/// Convert world coordinates to grid coordinates, including outside the grid.
/// C prototype: `void gridrawxygridxy(double *gvals, double dx, double dy, double *tx, double *ty);`
pub fn gridrawxygridxy(gvals: &[f64], mut dx: f64, mut dy: f64, tx: &mut f64, ty: &mut f64) {
    dx = dx - gvals[2];
    dy = dy - gvals[3];
    *tx = dx * gvals[17] + dy * gvals[16];
    *ty = (-dx * gvals[16] + dy * gvals[17]) * gvals[1];
}

/// Convert grid coordinates to world coordinates, including outside the grid.
/// C prototype: `void gridgridxyrawxy(double *gvals, double dx, double dy, double *tx, double *ty);`
pub fn gridgridxyrawxy(gvals: &[f64], dx: f64, mut dy: f64, tx: &mut f64, ty: &mut f64) {
    dy = dy * gvals[1];
    *tx = dx * gvals[17] - dy * gvals[16];
    *ty = dx * gvals[16] + dy * gvals[17];
    *tx = *tx + gvals[2];
    *ty = *ty + gvals[3];
}

/// Convert indexes to cell-centre grid coordinates, including outside indexes.
/// C prototype: `void gridicgridxy(double *gvals, int igi, int igc, double *dx, double *dy);`
pub fn gridicgridxy(gvals: &[f64], igi: i32, igc: i32, dx: &mut f64, dy: &mut f64) {
    *dx = igi.wrapping_sub(1) as f64 * gvals[10];
    *dy = igc.wrapping_sub(1) as f64 * gvals[11];
}

/// Convert indexes to cell-centre world coordinates, including outside indexes.
/// C prototype: `void gridicrawxy(double *gvals, int igi, int igc, double *dx, double *dy);`
pub fn gridicrawxy(gvals: &[f64], igi: i32, igc: i32, dx: &mut f64, dy: &mut f64) {
    let rx = igi.wrapping_sub(1) as f64 * gvals[10];
    let ry = igc.wrapping_sub(1) as f64 * gvals[11] * gvals[1];
    *dx = rx * gvals[17] - ry * gvals[16];
    *dy = rx * gvals[16] + ry * gvals[17];
    *dx = *dx + gvals[2];
    *dy = *dy + gvals[3];
}

/// Convert indexes to the CDP numbered along B; return -2147483645 outside.
/// C prototype: `void gridiccdp(double *gvals, int igi, int igc, int *icdp);`
pub fn gridiccdp(gvals: &[f64], igi: i32, igc: i32, icdp: &mut i32) {
    if igi < 1 || igi as f64 > gvals[12] || igc < 1 || igc as f64 > gvals[13] {
        *icdp = -2147483645;
    } else {
        *icdp = (gvals[14] + igi as f64 - 1.0 + igc.wrapping_sub(1) as f64 * gvals[12]) as i32;
    }
}

/// Convert indexes to the CDP numbered along C; return -2147483645 outside.
/// C prototype: `void gridiccdp90(double *gvals, int igi, int igc, int *jcdp);`
pub fn gridiccdp90(gvals: &[f64], igi: i32, igc: i32, jcdp: &mut i32) {
    if igi < 1 || igi as f64 > gvals[12] || igc < 1 || igc as f64 > gvals[13] {
        *jcdp = -2147483645;
    } else {
        *jcdp = (gvals[14] + igi.wrapping_sub(1) as f64 * gvals[13] + igc as f64 - 1.0) as i32;
    }
}

/// Convert a B-direction CDP to indexes; return -2147483645 for both outside.
/// C prototype: `void gridcdpic(double *gvals, int icdp, int *igi, int *igc);`
pub fn gridcdpic(gvals: &[f64], icdp: i32, igi: &mut i32, igc: &mut i32) {
    if (icdp as f64) < gvals[14] || icdp as f64 > gvals[15] {
        *igi = -2147483645;
        *igc = -2147483645;
        return;
    }
    let ncdp = (icdp as f64 - gvals[14]) as i32;
    let nwb = gvals[12] as i32;
    // ARM64 C division/remainder by zero produce 0/ncdp respectively.
    *igi = 1_i32.wrapping_add(ncdp.checked_rem(nwb).unwrap_or(ncdp));
    *igc = 1_i32.wrapping_add(ncdp.checked_div(nwb).unwrap_or(0));
}

/// Convert a C-direction CDP to indexes; return -2147483645 for both outside.
/// C prototype: `void gridcdpic90(double *gvals, int jcdp, int *igi, int *igc);`
pub fn gridcdpic90(gvals: &[f64], jcdp: i32, igi: &mut i32, igc: &mut i32) {
    if (jcdp as f64) < gvals[14] || jcdp as f64 > gvals[15] {
        *igi = -2147483645;
        *igc = -2147483645;
        return;
    }
    let ncdp = (jcdp as f64 - gvals[14]) as i32;
    let nwc = gvals[13] as i32;
    // ARM64 C division/remainder by zero produce 0/ncdp respectively.
    *igi = 1_i32.wrapping_add(ncdp.checked_div(nwc).unwrap_or(0));
    *igc = 1_i32.wrapping_add(ncdp.checked_rem(nwc).unwrap_or(ncdp));
}

/// Exercise conversions at all four corners, printing the C diagnostics.
/// icheck == 0 disables diagnostics; status is always zero.
/// C prototype: `void gridcheck(double *gvals, int icheck, int *errwarn);`
pub fn gridcheck(gvals: &[f64], icheck: i32, errwarn: &mut i32) {
    *errwarn = 0;
    if icheck == 0 {
        return;
    }

    let nwb = (gvals[12] + 0.1) as i32;
    let nwc = (gvals[13] + 0.1) as i32;
    let mut rx = 0.0;
    let mut ry = 0.0;
    let mut tx = 0.0;
    let mut ty = 0.0;
    let mut jcdp = 0;
    let mut jigi = 0;
    let mut jigc = 0;
    let mut kcdp = 0;
    let mut kigi = 0;
    let mut kigc = 0;

    gridicrawxy(gvals, 1, 1, &mut rx, &mut ry);
    warn!("gridicrawxy:     corner A raw  XYs= %.20f %.20f ", rx, ry);
    gridrawxygridxy(gvals, rx, ry, &mut tx, &mut ty);
    warn!("gridrawxygridxy: corner A grid XYs= %.20f %.20f ", tx, ty);
    gridicgridxy(gvals, 1, 1, &mut tx, &mut ty);
    warn!("gridicgridxy:    corner A grid XYs= %.20f %.20f ", tx, ty);
    gridgridxyrawxy(gvals, tx, ty, &mut rx, &mut ry);
    warn!("gridgridxyrawxy: corner A raw  XYs= %.20f %.20f ", rx, ry);
    gridrawxycdpic(gvals, rx, ry, &mut jcdp, &mut jigi, &mut jigc);
    gridcdpic(gvals, jcdp, &mut kigi, &mut kigc);
    gridiccdp(gvals, jigi, jigc, &mut kcdp);
    warn!(
        "gridrawxycdpic:          corner A cdp,igi,igc = %d %d %d ",
        jcdp, jigi, jigc
    );
    warn!(
        "gridcdpic and gridiccdp: corner A cdp,igi,igc = %d %d %d ",
        kcdp, kigi, kigc
    );

    gridicrawxy(gvals, nwb, 1, &mut rx, &mut ry);
    warn!("gridicrawxy:     corner B raw  XYs= %.20f %.20f ", rx, ry);
    gridrawxygridxy(gvals, rx, ry, &mut tx, &mut ty);
    warn!("gridrawxygridxy: corner B grid XYs= %.20f %.20f ", tx, ty);
    gridicgridxy(gvals, nwb, 1, &mut tx, &mut ty);
    warn!("gridicgridxy:    corner B grid XYs= %.20f %.20f ", tx, ty);
    gridgridxyrawxy(gvals, tx, ty, &mut rx, &mut ry);
    warn!("gridgridxyrawxy: corner B raw  XYs= %.20f %.20f ", rx, ry);
    gridrawxycdpic(gvals, rx, ry, &mut jcdp, &mut jigi, &mut jigc);
    gridcdpic(gvals, jcdp, &mut kigi, &mut kigc);
    gridiccdp(gvals, jigi, jigc, &mut kcdp);
    warn!(
        "gridrawxycdpic:          corner B cdp,igi,igc = %d %d %d ",
        jcdp, jigi, jigc
    );
    warn!(
        "gridcdpic and gridiccdp: corner B cdp,igi,igc = %d %d %d ",
        kcdp, kigi, kigc
    );

    gridicrawxy(gvals, 1, nwc, &mut rx, &mut ry);
    warn!("gridicrawxy:     corner C raw  XYs= %.20f %.20f ", rx, ry);
    gridrawxygridxy(gvals, rx, ry, &mut tx, &mut ty);
    warn!("gridrawxygridxy: corner C grid XYs= %.20f %.20f ", tx, ty);
    gridicgridxy(gvals, 1, nwc, &mut tx, &mut ty);
    warn!("gridicgridxy:    corner C grid XYs= %.20f %.20f ", tx, ty);
    gridgridxyrawxy(gvals, tx, ty, &mut rx, &mut ry);
    warn!("gridgridxyrawxy: corner C raw  XYs= %.20f %.20f ", rx, ry);
    gridrawxycdpic(gvals, rx, ry, &mut jcdp, &mut jigi, &mut jigc);
    gridcdpic(gvals, jcdp, &mut kigi, &mut kigc);
    gridiccdp(gvals, jigi, jigc, &mut kcdp);
    warn!(
        "gridrawxycdpic:          corner C cdp,igi,igc = %d %d %d ",
        jcdp, jigi, jigc
    );
    warn!(
        "gridcdpic and gridiccdp: corner C cdp,igi,igc = %d %d %d ",
        kcdp, kigi, kigc
    );

    gridicrawxy(gvals, nwb, nwc, &mut rx, &mut ry);
    warn!("gridicrawxy:     corner D raw  XYs= %.20f %.20f ", rx, ry);
    gridrawxygridxy(gvals, rx, ry, &mut tx, &mut ty);
    warn!("gridrawxygridxy: corner D grid XYs= %.20f %.20f ", tx, ty);
    gridicgridxy(gvals, nwb, nwc, &mut tx, &mut ty);
    warn!("gridicgridxy:    corner D grid XYs= %.20f %.20f ", tx, ty);
    gridgridxyrawxy(gvals, tx, ty, &mut rx, &mut ry);
    warn!("gridgridxyrawxy: corner D raw  XYs= %.20f %.20f ", rx, ry);
    gridrawxycdpic(gvals, rx, ry, &mut jcdp, &mut jigi, &mut jigc);
    gridcdpic(gvals, jcdp, &mut kigi, &mut kigc);
    gridiccdp(gvals, jigi, jigc, &mut kcdp);
    warn!(
        "gridrawxycdpic:          corner D cdp,igi,igc = %d %d %d ",
        jcdp, jigi, jigc
    );
    warn!(
        "gridcdpic and gridiccdp: corner D cdp,igi,igc = %d %d %d ",
        kcdp, kigi, kigc
    );
}
