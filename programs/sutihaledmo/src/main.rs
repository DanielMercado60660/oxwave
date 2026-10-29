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

//! SUTIHALEDMO - TI Hale Dip MoveOut (based on Hale's PhD thesis).
//! Port of su/main/dip_moveout/sutihaledmo.c (Revision 1.5).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: (Visitor to CSM from Mobil) John E. Anderson Spring 1994.
//! References: Anderson, J.E., and Tsvankin, I., 1994, Dip-moveout by
//! Fourier transform in anisotropic media, CWP-146.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{cadd, cmplx, cmul, crmul};
use oxwave_cwp::cwp::PI;
use oxwave_f_vnd::*;
use oxwave_pfa_fft::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUTIHALEDMO - TI Hale Dip MoveOut (based on Hale's PhD thesis)\t",
    "\t\t\t\t\t\t\t\t\t",
    "  sutihaledmo <infile >outfile [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " nxmax\t\t  maximum number of midpoints in common offset gather",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " option=1\t\t1 = traditional Hale DMO (from PhD thesis)\t",
    "\t\t\t2 = Bleistein's true amplitude DMO\t\t",
    "\t\t\t3 = Bleistein's cos*cos weighted DMO\t\t",
    "\t\t\t4 = Zhang's DMO\t\t\t\t\t",
    "\t\t\t5 = Tsvankin's anisotropic DMO\t\t\t",
    "\t\t\t6 = Tsvankin's VTI DMO weak anisotropy approximation",
    " dx=50.\t\t midpoint sampling interval between traces\t",
    "\t\t\tin a common offset gather.  (usually shot\t",
    "\t\t\tinterval in meters)\t\t\t\t",
    " v=1500.0\t\tvelocity (in meters/sec)\t\t\t",
    "\t\t\t(must enter a positive value for option=3)\t",
    "\t\t\t(for excluding evanescent energy)\t\t",
    " h=200.0\t\tsource-receiver half-offset (in meters)\t\t",
    " ntpad=0\t\tnumber of time samples to pad\t\t\t",
    " nxpad=h/dx\t\tnumber of midpoints to pad\t\t\t",
    " file=vnmo\t\tname of file with vnmo as a function of p\t",
    "\t\t\tused for option=5--otherwise not used\t\t",
    "\t\t\t(Generate this file by running program\t\t",
    "\t\t\tsutivel with appropriate list of Thomsen's\t",
    "\t\t\tparameters.)\t\t\t\t\t",
    " e=0.\t\t\tThompsen's epsilon\t\t\t\t",
    " d=0.\t\t\tThompsen's delta\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Note:\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " This module assumes a single common offset gather after NMO is\t",
    " to be input, DMO corrected, and output.  It is useful for computing\t",
    " theoretical DMO impulse responses.  The Hale algorithm is\t\t",
    " computationally intensive and not commonly used for bulk processing\t",
    " of all of the offsets on a 2-D line as there are cheaper alternative\t",
    " algorithms.  The Hale algorithm is commonly used in theoretical studies.",
    " Bulk processing for multiple common offset gathers is typically done\t",
    " using other modules.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Test run:   suspike | sutihaledmo nxmax=32 option=1 v=1500 | suxwigb & ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get info from first trace */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (0.000001 * tr.dt as f64) as f32;

    /* get optional parameters */
    let mut dx = 50.0f32;
    let mut v = 1500.0f32;
    let mut h = 200.0f32;
    let mut e = 0.0f32;
    let mut d = 0.0f32;
    getparfloat("dx", &mut dx);
    getparfloat("v", &mut v);
    getparfloat("h", &mut h);
    getparfloat("e", &mut e);
    getparfloat("d", &mut d);
    let mut nxpad = (h / dx) as i32;
    let mut ntpad = 0i32;
    let mut option = 1i32;
    let mut file = "vnmo".to_string();
    getparint("nxpad", &mut nxpad);
    getparint("ntpad", &mut ntpad);
    getparint("option", &mut option);
    getparstring("file", &mut file);
    if v <= 0.0 && (option == 3 || option == 6) {
        err!("must enter positive value for v when option=3");
    }
    let mut nxmax = 0i32;
    if getparint("nxmax", &mut nxmax) == 0 {
        err!("must enter nxmax");
    }
    checkpars();

    let hfp = etmpfile();
    let mut np = 0i32;
    let mut dp = 0.0f32;
    let mut vnmo = Vec::new();
    if option == 5 {
        let fp = efopen(&file, "r");
        np = fscanf_i32(&fp).unwrap_or_else(|| err!("invalid np in vnmo(p) file"));
        if np < 0 {
            err!("invalid np in vnmo(p) file");
        }
        vnmo = vec![0.0f32; np as usize];
        dp = fscanf_f32(&fp).unwrap_or_else(|| err!("missing dp in vnmo(p) file"));
        fprintf!(stderr(), "suhaledmo: np=%d dp=%e\n", np, dp);
        for ip in 0..np as usize {
            vnmo[ip] = fscanf_f32(&fp).unwrap_or_else(|| err!("missing velocity in vnmo(p) file"));
        }
        fclose(&fp);
    }

    let ntfft = npfao(
        (nt as i32).wrapping_add(ntpad),
        2i32.wrapping_mul((nt as i32).wrapping_add(ntpad)),
    ) as usize;
    let nxfft = npfar(nxmax.wrapping_add(nxpad)) as usize;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    fprintf!(
        stderr(),
        "suhaledmo: ntfft=%d nxfft=%d\n",
        ntfft as i32,
        nxfft as i32
    );
    let scale = (1.0 / nxfft as f64) as f32;
    let nmax = (nxfft / 2 + 1).max(nt);
    let mut crt = vec![Complex::default(); nmax];
    let mut ctemp = vec![Complex::default(); ntfft];
    let mut rt = vec![0.0f32; nxfft + 2];

    let fname = VNDtempname("suhaledmo");
    // SAFETY: VNDtempname returns an allocated, NUL-terminated ASCII name.
    // Like the C main, retain this allocation through program exit.
    let fname = unsafe { std::ffi::CStr::from_ptr(fname.cast()) }
        .to_str()
        .unwrap();
    let mut vnd = V2Dop(2, 100000, 4, fname, (nxfft + 2) as i64, nt as i64)
        .unwrap_or_else(|| err!("couldn't create VND panel"));
    let mut nx = 0usize;
    /* Main loop for saving input traces */
    loop {
        for it in 0..nt {
            tr.data[it] *= scale;
        }
        // More traces than real panel slots makes the C VND copy overrun
        // its buffer (its key check also accepts the one-past-end index).
        if nx >= vnd.n[0] as usize {
            err!("input trace exceeds VND panel dimensions");
        }
        V2Dw1(&mut vnd, nx as i64, tr.data_bytes_mut(), 1);
        efwrite(std::slice::from_ref(&tr.hdr), &hfp);
        nx += 1;
        if nx as i64 > nxmax as i64 {
            break;
        }
        if gettr(&mut tr) == 0 {
            break;
        }
    }

    /* Do forward FFT from x to k */
    for it in 0..nt {
        V2Dr0(
            &mut vnd,
            it as i64,
            oxwave_core::cio::as_bytes_mut(&mut rt),
            3,
        );
        for ix in nx..nxfft {
            rt[ix] = 0.0;
        }
        pfarc(1, &rt[..nxfft], &mut crt[..nxfft / 2 + 1]);
        V2Dw0(&mut vnd, it as i64, cbytes(&mut crt), 4);
    }

    /* change dimension 0 from real to complex; correct each wavenumber */
    VNDr2c(&mut vnd);
    for ik in 0..nxfft / 2 + 1 {
        V2Dr1(&mut vnd, ik as i64, cbytes(&mut crt), 5);
        let k = ik as f32 * dk;
        // Each active operator copies nt entries from its ntfft workspace.
        // C reads out of bounds when ntfft < nt. Reject only an operator
        // that reaches that copy; preserve the defined early-return paths.
        if ntfft < nt
            && (option == 3 || (matches!(option, 1 | 2 | 4 | 5 | 6) && !(h * k < 1.0e-20f32)))
        {
            VNDcl(vnd, 1);
            err!("ntfft=%d is smaller than nt=%d", ntfft as i32, nt as i32);
        }
        if option == 1 {
            haledmo(h, k, &mut crt, &mut ctemp, nt, ntfft, dt, v);
        }
        if option == 2 {
            bleisteindmo(h, k, &mut crt, &mut ctemp, nt, ntfft, dt, v);
        }
        if option == 3 {
            bleisteindmo2(h, k, &mut crt, &mut ctemp, nt, ntfft, dt, v);
        }
        if option == 4 {
            zhangdmo(h, k, &mut crt, &mut ctemp, nt, ntfft, dt, v);
        }
        if option == 5 {
            tsvankindmo(h, k, &mut crt, &mut ctemp, nt, ntfft, dt, np, dp, &vnmo);
        }
        if option == 6 {
            weakdmo(h, k, &mut crt, &mut ctemp, nt, ntfft, dt, e, d, v);
        }
        V2Dw1(&mut vnd, ik as i64, cbytes(&mut crt), 6);
    }

    /* Do inverse FFT from k to x */
    for it in 0..nt {
        V2Dr0(&mut vnd, it as i64, cbytes(&mut crt), 7);
        pfacr(-1, &crt[..nxfft / 2 + 1], &mut rt[..nxfft]);
        // C aliases the real and complex FFT buffers. Preserve its last
        // two real slots, which are untouched by the inverse transform.
        cbytes(&mut crt)[..nxfft * 4].copy_from_slice(oxwave_core::cio::as_bytes(&rt[..nxfft]));
        V2Dw0(&mut vnd, it as i64, cbytes(&mut crt), 8);
    }

    /* change dimension 0 back to real and output the DMO result */
    VNDc2r(&mut vnd);
    efseeko(&hfp, 0, SEEK_SET);
    for ix in 0..nx {
        V2Dr1(&mut vnd, ix as i64, tr.data_bytes_mut(), 7);
        efread(std::slice::from_mut(&mut tr.hdr), &hfp);
        puttr(&mut tr);
    }
    VNDcl(vnd, 1);
    if VNDtotalmem() != 0 {
        fprintf!(stderr(), "total VND memory at end = %ld\n", VNDtotalmem());
    }
    exit(CWP_Exit());
}

// C's (char *)crt view for the byte-oriented VND API.
fn cbytes(values: &mut [Complex]) -> &mut [u8] {
    // SAFETY: Complex is repr(C), two contiguous f32 fields with no padding;
    // every bit pattern is valid, and the borrow retains exclusive access.
    unsafe {
        std::slice::from_raw_parts_mut(values.as_mut_ptr().cast(), std::mem::size_of_val(values))
    }
}

// Hale fk DMO
#[allow(non_snake_case, clippy::too_many_arguments)]
fn haledmo(
    h: f32,
    k: f32,
    p: &mut [Complex],
    pwork: &mut [Complex],
    nt: usize,
    ntfft: usize,
    dt: f32,
    v: f32,
) {
    let hk = h * k;
    let eps = 1.0e-20f32;
    if hk < eps {
        return;
    }
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let scale = (1.0 / ntfft as f64) as f32;
    let iwtest = ntfft / 2 + 1;
    pwork[0] = cmplx(0.0, 0.0);
    for iw in 1..ntfft {
        let w = if iw < iwtest {
            iw as f32 * dw
        } else {
            (iw as i32 - ntfft as i32) as f32 * dw
        };
        let mut csum = cmplx(0.0, 0.0);
        if ((k * v) as f64) < 2.0 * (w as f64).abs() {
            let hkpw = hk / w;
            for it in 1..nt {
                let tn = it as f32 * dt;
                let hkpwt = hkpw / tn;
                let A = (1.0 + (hkpwt * hkpwt) as f64).sqrt() as f32;
                let phase = w * tn * A;
                let cfac = cmplx(
                    ((phase as f64).cos() / A as f64) as f32,
                    ((phase as f64).sin() / A as f64) as f32,
                );
                csum = cadd(csum, cmul(p[it], cfac));
            }
        }
        pwork[iw] = csum;
    }
    pfacc(-1, &mut pwork[..ntfft]);
    for it in 0..nt {
        p[it] = crmul(pwork[it], scale);
    }
}

// Bleistein true amplitude fk DMO
#[allow(non_snake_case, clippy::too_many_arguments)]
fn bleisteindmo(
    h: f32,
    k: f32,
    p: &mut [Complex],
    pwork: &mut [Complex],
    nt: usize,
    ntfft: usize,
    dt: f32,
    v: f32,
) {
    let hk = h * k;
    let eps = 1.0e-20f32;
    if hk < eps {
        return;
    }
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let scale = (1.0 / ntfft as f64) as f32;
    let iwtest = ntfft / 2 + 1;
    pwork[0] = cmplx(0.0, 0.0);
    for iw in 1..ntfft {
        let w = if iw < iwtest {
            iw as f32 * dw
        } else {
            (iw as i32 - ntfft as i32) as f32 * dw
        };
        let mut csum = cmplx(0.0, 0.0);
        if ((k * v) as f64) < 2.0 * (w as f64).abs() {
            let hkpw = hk / w;
            for it in 1..nt {
                let tn = it as f32 * dt;
                let hkpwt = hkpw / tn;
                let A = (1.0 + (hkpwt * hkpwt) as f64).sqrt() as f32;
                let B = (1.0 + (2.0f32 * hkpwt * hkpwt) as f64).sqrt() as f32;
                let phase = w * tn * A;
                let cfac = cmplx(
                    (B as f64 * (phase as f64).cos() / A as f64) as f32,
                    (B as f64 * (phase as f64).sin() / A as f64) as f32,
                );
                csum = cadd(csum, cmul(p[it], cfac));
            }
        }
        pwork[iw] = csum;
    }
    pfacc(-1, &mut pwork[..ntfft]);
    for it in 0..nt {
        p[it] = crmul(pwork[it], scale);
    }
}

// Bleistein cos*cos weighted amplitude fk DMO
#[allow(non_snake_case, clippy::too_many_arguments)]
fn bleisteindmo2(
    h: f32,
    k: f32,
    p: &mut [Complex],
    pwork: &mut [Complex],
    nt: usize,
    ntfft: usize,
    dt: f32,
    v: f32,
) {
    let hk = h * k;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let scale = (1.0 / ntfft as f64) as f32;
    let iwtest = ntfft / 2 + 1;
    pwork[0] = cmplx(0.0, 0.0);
    let hpv = (2.0 * h as f64 / v as f64) as f32;
    for iw in 1..ntfft {
        let w = if iw < iwtest {
            iw as f32 * dw
        } else {
            (iw as i32 - ntfft as i32) as f32 * dw
        };
        let mut csum = cmplx(0.0, 0.0);
        if ((k * v) as f64) < 2.0 * (w as f64).abs() {
            let hkpw = hk / w;
            for it in 1..nt {
                let tn = it as f32 * dt;
                let hkpwt = hkpw / tn;
                let hkpwt2 = hkpwt * hkpwt;
                let mut hpvt2 = hpv / tn;
                hpvt2 = hpvt2 * hpvt2;
                let A = (1.0 + hkpwt2 as f64).sqrt() as f32;
                let B = ((1.0 + (2.0f32 * hkpwt2) as f64) * (1.0f32 + hkpwt2) as f64
                    / (1.0f32 + hpvt2) as f64) as f32;
                let phase = w * tn * A;
                let cfac = cmplx(
                    (B as f64 * (phase as f64).cos() / A as f64) as f32,
                    (B as f64 * (phase as f64).sin() / A as f64) as f32,
                );
                csum = cadd(csum, cmul(p[it], cfac));
            }
        }
        pwork[iw] = csum;
    }
    pfacc(-1, &mut pwork[..ntfft]);
    for it in 0..nt {
        p[it] = crmul(pwork[it], scale);
    }
}

// Zhang fk DMO
#[allow(non_snake_case, clippy::too_many_arguments)]
fn zhangdmo(
    h: f32,
    k: f32,
    p: &mut [Complex],
    pwork: &mut [Complex],
    nt: usize,
    ntfft: usize,
    dt: f32,
    v: f32,
) {
    let hk = h * k;
    let eps = 1.0e-20f32;
    if hk < eps {
        return;
    }
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let scale = (1.0 / ntfft as f64) as f32;
    let iwtest = ntfft / 2 + 1;
    pwork[0] = cmplx(0.0, 0.0);
    for iw in 1..ntfft {
        let w = if iw < iwtest {
            iw as f32 * dw
        } else {
            (iw as i32 - ntfft as i32) as f32 * dw
        };
        let mut csum = cmplx(0.0, 0.0);
        if ((k * v) as f64) < 2.0 * (w as f64).abs() {
            let hkpw = hk / w;
            for it in 1..nt {
                let tn = it as f32 * dt;
                let hkpwt = hkpw / tn;
                let hkpwt2 = hkpwt * hkpwt;
                let A = (1.0 + hkpwt2 as f64).sqrt() as f32;
                let B = ((1.0 + (2.0f32 * hkpwt2) as f64) / (1.0f32 + hkpwt2) as f64) as f32;
                let phase = w * tn * A;
                let cfac = cmplx(
                    (B as f64 * (phase as f64).cos() / A as f64) as f32,
                    (B as f64 * (phase as f64).sin() / A as f64) as f32,
                );
                csum = cadd(csum, cmul(p[it], cfac));
            }
        }
        pwork[iw] = csum;
    }
    pfacc(-1, &mut pwork[..ntfft]);
    for it in 0..nt {
        p[it] = crmul(pwork[it], scale);
    }
}

// Tsvankin fk DMO, with tabulated NMO velocity as a function of ray parameter.
#[allow(non_snake_case, clippy::too_many_arguments)]
fn tsvankindmo(
    h: f32,
    k: f32,
    p: &mut [Complex],
    pwork: &mut [Complex],
    nt: usize,
    ntfft: usize,
    dt: f32,
    np: i32,
    dp: f32,
    vnmo: &[f32],
) {
    let hk = h * k;
    let eps = 1.0e-20f32;
    if hk < eps {
        return;
    }
    let npm = np - 2;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let scale = (1.0 / ntfft as f64) as f32;
    let iwtest = ntfft / 2 + 1;
    pwork[0] = cmplx(0.0, 0.0);
    for iw in 1..ntfft {
        let w = if iw < iwtest {
            iw as f32 * dw
        } else {
            (iw as i32 - ntfft as i32) as f32 * dw
        };
        let mut csum = cmplx(0.0, 0.0);
        let pp = (k as f64 / (2.0 * w as f64)).abs() as f32;
        let rp = pp / dp;
        let ip = rp as i32;
        if ip < npm {
            // A negative table index is undefined in C, not a Rust panic.
            if ip < 0 {
                err!("negative index in vnmo(p) file");
            }
            let wgt = rp - ip as f32;
            let v = (1.0 - wgt as f64) * vnmo[ip as usize] as f64
                + (wgt * vnmo[ip as usize + 1]) as f64;
            let v0 = vnmo[0] as f64;
            let fac = (4.0f32 * h * h) as f64 * (1.0 / (v0 * v0) - 1.0 / (v * v));
            let it1 = if fac < 0.0 {
                (2.0 + fac.abs().sqrt() / dt as f64) as i32
            } else {
                1
            };
            for it in it1..nt as i32 {
                let tn = it as f32 * dt;
                let tn2 = tn * tn;
                let A = (1.0 + fac / tn2 as f64).sqrt() as f32;
                let phase = w * tn * A;
                let cfac = cmplx(
                    ((phase as f64).cos() / A as f64) as f32,
                    ((phase as f64).sin() / A as f64) as f32,
                );
                csum = cadd(csum, cmul(p[it as usize], cfac));
            }
        }
        pwork[iw] = csum;
    }
    pfacc(-1, &mut pwork[..ntfft]);
    for it in 0..nt {
        p[it] = crmul(pwork[it], scale);
    }
}

// Tsvankin fk DMO weak anisotropy approximation.
#[allow(non_snake_case, clippy::too_many_arguments)]
fn weakdmo(
    h: f32,
    k: f32,
    p: &mut [Complex],
    pwork: &mut [Complex],
    nt: usize,
    ntfft: usize,
    dt: f32,
    e: f32,
    d: f32,
    v: f32,
) {
    let hk = h * k;
    let eps = 1.0e-20f32;
    if hk < eps {
        return;
    }
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let scale = (1.0 / ntfft as f64) as f32;
    let iwtest = ntfft / 2 + 1;
    pwork[0] = cmplx(0.0, 0.0);
    for iw in 1..ntfft {
        let w = if iw < iwtest {
            iw as f32 * dw
        } else {
            (iw as i32 - ntfft as i32) as f32 * dw
        };
        let mut csum = cmplx(0.0, 0.0);
        let pp = (k as f64 / (2.0 * w as f64)).abs() as f32;
        if pp < 1.0f32 / v {
            let y = pp * pp * v * v;
            let fac = (4.0f32
                * h
                * h
                * pp
                * pp
                * (1.0f32 + 2.0f32 * (e - d) * (4.0f32 * y * y - 9.0f32 * y + 6.0f32)))
                as f64;
            let it1 = if fac < 0.0 {
                (2.0 + fac.abs().sqrt() / dt as f64) as i32
            } else {
                1
            };
            for it in it1..nt as i32 {
                let tn = it as f32 * dt;
                let tn2 = tn * tn;
                let A = (1.0 + fac / tn2 as f64).sqrt() as f32;
                let phase = w * tn * A;
                let cfac = cmplx(
                    ((phase as f64).cos() / A as f64) as f32,
                    ((phase as f64).sin() / A as f64) as f32,
                );
                csum = cadd(csum, cmul(p[it as usize], cfac));
            }
        }
        pwork[iw] = csum;
    }
    pfacc(-1, &mut pwork[..ntfft]);
    for it in 0..nt {
        p[it] = crmul(pwork[it], scale);
    }
}
