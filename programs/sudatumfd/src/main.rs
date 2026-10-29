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

//! SUDATUMFD - 2D zero-offset finite difference acoustic wave-equation datuming.
//! Port of su/main/datuming/sudatumfd.c.
//! Author: Chris Robinson, 10/16/00, CWP, Colorado School of Mines.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_pfa_fft::{npfar, pfacr, pfarc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUDATUMFD - 2D zero-offset Finite Difference acoustic wave-equation\t",
    "\t\t DATUMing    \t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sudatumfd <stdin > stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t   \t",
    "\t\t\t\t\t\t\t\t\t",
    " nt=\t   number of time samples on each trace\t       \t\t\t",
    " nx=\t   number of receivers per shot gather\t\t\t\t",
    " nsx=\t  number of shot gathers\t\t\t\t    \t",
    " nz=\t   number of downward continuation depth steps\t\t\t",
    " dz=\t   depth sampling interval (in meters)\t\t\t\t",
    " mx=\t   number of horizontal samples in the velocity model\t\t",
    " mz=\t   number of vertical samples in the velocity model\t\t",
    " vfile1=       velocity file used for thin-lens term\t    \t\t",
    " vfile2=       velocity file used for diffraction term\t\t\t",
    " dx=           horizontal sampling interval (in meters)                ",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t   \t",
    "\t\t\t\t\t\t\t\t\t",
    " dt=.004       time sampling interval (in seconds)\t\t\t",
    " buff=5\tnumber of zero traces added to each side of each   \t",
    "\t     shot gather as a pad\t\t\t       \t\t",
    " tap_len=5     taper length (in number of traces)\t\t\t",
    " x_0=0.0       x coordinate of leftmost position in velocity model     ",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The algorithm is a 45-degree implicit-finite-difference scheme based  ",
    " on the one-way wave equation.  It works on poststack (zero-offset)    ",
    " data only.  The two velocity files, vfile1 and vfile2, are binary     ",
    " files containing floats with the format v[ix][iz].  There are two     ",
    " potentially different velocity files for the thin-lens and            ",
    " diffraction terms to allow for the use of a zero-velocity layer       ",
    " which allows for datuming from an irregular surface.                  ",
    "                                                                       ",
    " Source and receiver locations must be set in the header values in     ",
    " order for the datuming to work properly.  The leftmost position of    ",
    " of the velocity models given in vfile1 and vfile2 must also be given. ",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

// C alloc3complex data are contiguous, with frequency varying fastest.
type Volume = Vec<Complex>;

unsafe extern "C" {
    fn raise(sig: i32) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut nt, mut nx, mut nz, mut nsx) = (0i32, 0i32, 0i32, 0i32);
    let (mut v, mut dz, mut x_0) = (1000.0f32, 0.0f32, 0.0f32);
    let (mut buff, mut tap_len, mut mx, mut mz) = (0i32, 0i32, 0i32, 0i32);
    let (mut vfile1, mut vfile2) = (String::new(), String::new());
    if getparint("nt", &mut nt) == 0 {
        err!("must specify nt");
    }
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz");
    }
    if getparint("nsx", &mut nsx) == 0 {
        err!("must specify nsx");
    }
    if getparfloat("v", &mut v) == 0 {
        v = 1000.0;
    }
    if getparfloat("dz", &mut dz) == 0 {
        err!("must specify dz");
    }
    if getparfloat("x_0", &mut x_0) == 0 {
        err!("must specify x_0");
    }
    if getparint("buff", &mut buff) == 0 {
        err!("must specify buff");
    }
    if getparint("tap_len", &mut tap_len) == 0 {
        err!("must specify tap_len");
    }
    if getparint("mx", &mut mx) == 0 {
        err!("must specify mx");
    }
    if getparint("mz", &mut mz) == 0 {
        err!("must specify mz");
    }
    if getparstring("vfile1", &mut vfile1) == 0 {
        err!("must specify vfile1");
    }
    if getparstring("vfile2", &mut vfile2) == 0 {
        err!("must specify vfile2");
    }

    let vfp1 = match fopen(&vfile1, "r") {
        Some(fp) => fp,
        None => {
            // C passes the unchecked NULL FILE* to fread and dies with SIGSEGV.
            unsafe {
                signal(11, 0);
                raise(11);
            }
            exit(139);
        }
    };
    let mut vel1 = vec![0.0f32; (mz * mx) as usize];
    fread(&mut vel1, &vfp1);
    fclose(&vfp1);
    let vfp2 = match fopen(&vfile2, "r") {
        Some(fp) => fp,
        None => {
            // Preserve the original unchecked fopen failure, including no message.
            unsafe {
                signal(11, 0);
                raise(11);
            }
            exit(139);
        }
    };
    let mut vel2 = vec![0.0f32; (mz * mx) as usize];
    fread(&mut vel2, &vfp2);
    fclose(&vfp2);
    let nx_pad = nx + 2 * buff;
    let nt_alloc = nt as usize;
    let mut intime = vec![vec![vec![0.0f32; nt_alloc]; nsx as usize]; nx as usize];
    let mut gx = vec![vec![0i32; nsx as usize]; nx as usize];
    let mut sx = gx.clone();
    let mut cdp = gx.clone();
    let mut offset = gx.clone();
    let (mut ix, mut isx) = (0usize, 0usize);
    let (mut dt, mut dx) = (0.0f32, 0.0f32);
    // The header replaces nt, and dt/dx are queried for every input trace.
    while gettr(&mut tr) != 0 && ix < nx as usize {
        nt = tr.ns as i32;
        if getparfloat("dt", &mut dt) == 0 {
            if tr.dt != 0 {
                dt = (tr.dt as f64 / 1000000.0) as f32;
            } else {
                dt = 0.004;
                warn!("tr.dt not set, assuming dt=0.004");
            }
        }
        if getparfloat("dx", &mut dx) == 0 {
            if tr.d2 != 0.0 {
                dx = tr.d2;
            } else {
                dx = 1.0;
                warn!("tr.d2 not set, assuming d2=1.0");
            }
        }
        checkpars();
        // C writes beyond the allocated gather on excess traces. Reject that
        // undefined input before indexing either the headers or the samples.
        if isx >= nsx as usize {
            err!("input contains more than nx * nsx traces");
        }
        // The header replaces nt, but the original allocation does not grow.
        // C corrupts memory if ns exceeds that allocation; smaller ns is valid.
        if nt as usize > nt_alloc {
            err!("trace ns=%d exceeds nt=%d allocation", nt, nt_alloc as i32);
        }
        gx[ix][isx] = tr.gx;
        sx[ix][isx] = tr.sx;
        cdp[ix][isx] = tr.cdp;
        offset[ix][isx] = tr.offset;
        intime[ix][isx][..nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        ix += 1;
        if ix == nx as usize {
            ix = 0;
            isx += 1;
        }
    }

    // Temporal Fourier transform.
    let ntfft = npfar(2 * nt);
    let nw = ntfft / 2 + 1;
    let pi = std::f64::consts::PI as f32;
    let dw = (2.0 * pi as f64 / (((ntfft - 1) as f32 * dt) as f64)) as f32;
    let mut ifr = vec![Complex::default(); (nw * nsx * nx) as usize];
    for isx in 0..nsx as usize {
        for ix in 0..nx as usize {
            let mut storage1 = vec![0.0f32; nt as usize];
            for it in 0..nt as usize {
                storage1[it] = intime[ix][isx][it];
            }
            let mut p1 = vec![0.0f32; ntfft as usize];
            let mut cp1 = vec![Complex::default(); nw as usize];
            p1[..nt as usize].copy_from_slice(&storage1);
            for it in nt as usize..ntfft as usize {
                p1[it] = 0.0;
            }
            pfarc(1, &p1, &mut cp1);
            for iw in 0..nw as usize {
                ifr[((ix) as usize * nsx as usize + (isx) as usize) * nw as usize
                    + (iw) as usize] = cp1[iw];
            }
        }
    }
    drop(intime);
    let mut ifr1 = vec![Complex::default(); (nw * nsx * nx_pad) as usize];
    // The C array is uninitialized for nz<=0. Fresh allocations are zeroed
    // in the tested original; use deterministic zeroes rather than Rust UB.
    let mut ifr2 = ifr1.clone();
    padsub(&ifr, &mut ifr1, nsx, nx, nw, nx_pad, buff);
    drop(ifr);
    for j in 0..nz {
        taper(&ifr1, &mut ifr2, nsx, nx_pad, nw, dw, tap_len, buff, nx);
        tlens(
            &ifr2, &mut ifr1, nsx, nx_pad, nw, dw, dz, v, dx, &vel1, mz, j, &gx, buff, x_0,
        );
        diff(
            &ifr1, &mut ifr2, nsx, nx_pad, nw, dw, dz, v, dx, &vel2, mz, j, &gx, buff, x_0,
        );
        ifr1.clone_from(&ifr2);
    }
    drop(ifr1);
    let mut ifr3 = vec![Complex::default(); (nw * nsx * nx) as usize];
    unpadsub(&ifr2, &mut ifr3, nsx, nx, nw, nx_pad, buff);
    drop(ifr2);
    let ofr1 = ifr3.clone();
    drop(ifr3);

    // Temporal inverse Fourier transform.
    let mut outtime = vec![vec![vec![0.0f32; nt as usize]; nsx as usize]; nx as usize];
    for isx in 0..nsx as usize {
        for ix in 0..nx as usize {
            let mut p2 = vec![0.0f32; ntfft as usize];
            let mut cp2 = vec![Complex::default(); nw as usize];
            for iw in 0..nw as usize {
                cp2[iw] = ofr1
                    [((ix) as usize * nsx as usize + (isx) as usize) * nw as usize + (iw) as usize];
            }
            pfacr(-1, &cp2, &mut p2);
            for it in 0..nt as usize {
                outtime[ix][isx][it] = (p2[it] as f64 / (1.0 * ntfft as f64)) as f32;
            }
        }
    }
    // All other header fields come from the final input trace.
    for isx in 0..nsx as usize {
        for ix in 0..nx as usize {
            tr.ns = nt as u16;
            // The C conversion first converts to unsigned int, then stores 16 bits.
            tr.dt = (dt * 1000000.0f32) as u32 as u16;
            tr.gx = gx[ix][isx];
            tr.sx = sx[ix][isx];
            tr.cdp = cdp[ix][isx];
            tr.offset = offset[ix][isx];
            tr.data[..nt as usize].copy_from_slice(&outtime[ix][isx]);
            puttr(&mut tr);
        }
    }
    exit(CWP_Exit());
}

#[allow(clippy::too_many_arguments, unused_assignments)]
fn tlens(
    val_ifr_pad: &Volume,
    val_ifr1: &mut Volume,
    val_nsx: i32,
    val_nx_pad: i32,
    val_nw: i32,
    val_dw: f32,
    val_dz: f32,
    mut val_v: f32,
    val_dx: f32,
    val_vel: &[f32],
    mz: i32,
    val_zlev: i32,
    val_gx: &[Vec<i32>],
    val_buff: i32,
    val_x_0: f32,
) {
    let mut a = 0;
    for i in 0..val_nsx as usize {
        for j in 0..val_nx_pad {
            if j < val_buff {
                a = val_gx[0][i];
            }
            if j > val_nx_pad - 1 - val_buff {
                a = val_gx[(val_nx_pad - 1 - 2 * val_buff) as usize][i];
            }
            if j >= val_buff && j <= val_nx_pad - 1 - val_buff {
                a = val_gx[(j - val_buff) as usize][i];
            }
            val_v = (0.5
                * val_vel[(((a as f32 - val_x_0) / val_dx) as i32 * mz + val_zlev) as usize] as f64)
                as f32;
            let ig = ((-1.0 / val_v as f64) * val_dz as f64) as f32;
            for k in 0..val_nw as usize {
                let angle = (k as f32 * val_dw * ig) as f64;
                let mult = cmplx(angle.cos() as f32, angle.sin() as f32);
                val_ifr1[((j as usize) as usize * val_nsx as usize + (i) as usize)
                    * val_nw as usize
                    + (k) as usize] = cmul(
                    val_ifr_pad[((j as usize) as usize * val_nsx as usize + (i) as usize)
                        * val_nw as usize
                        + (k) as usize],
                    mult,
                );
            }
        }
    }
}

#[allow(clippy::too_many_arguments, unused_assignments)]
fn diff(
    val_ifr1: &Volume,
    val_ifr2: &mut Volume,
    val_nsx: i32,
    val_nx_pad: i32,
    val_nw: i32,
    val_dw: f32,
    val_dz: f32,
    mut val_v: f32,
    val_dx: f32,
    val_vel: &[f32],
    mz: i32,
    val_zlev: i32,
    val_gx: &[Vec<i32>],
    val_buff: i32,
    val_x_0: f32,
) {
    let sigma = 0.5f32;
    let tau = 0.25f32;
    let mut e = vec![Complex::default(); val_nx_pad as usize];
    let mut f = e.clone();
    let mut d = e.clone();
    let mut ohm = e.clone();
    let mut alpha = e.clone();
    let mut beta = e.clone();
    let mut gamma = e.clone();
    let mut a = 0;
    let last = (val_nx_pad - 1) as usize;
    for i in 0..val_nsx as usize {
        for k in 0..val_nw as usize {
            let mut w = k as f32 * val_dw;
            if w < 5.0 {
                w = 5.0;
            }
            for n in 0..val_nx_pad {
                if n < val_buff {
                    a = val_gx[0][i];
                }
                if n > val_nx_pad - 1 - val_buff {
                    a = val_gx[(val_nx_pad - 1 - 2 * val_buff) as usize][i];
                }
                if n >= val_buff && n <= val_nx_pad - 1 - val_buff {
                    a = val_gx[(n - val_buff) as usize][i];
                }
                val_v = (0.5
                    * val_vel[(((a as f32 - val_x_0) / val_dx) as i32 * mz + val_zlev) as usize]
                        as f64) as f32;
                let n = n as usize;
                alpha[n] = cmplx(
                    0.0,
                    (-1.0 * val_v as f64 * sigma as f64 * val_dz as f64
                        / ((val_dx * val_dx) as f64 * 2.0 * w as f64)) as f32,
                );
                beta[n] = cmul(
                    alpha[n],
                    cadd(
                        cmplx(1.0, 0.0),
                        cmplx(
                            0.0,
                            (-2.0 * val_v as f64 * tau as f64 / (w * sigma * val_dz) as f64) as f32,
                        ),
                    ),
                );
                gamma[n] = cmul(
                    alpha[n],
                    cadd(
                        cmplx(1.0, 0.0),
                        cmplx(
                            0.0,
                            (2.0 * val_v as f64 * tau as f64 / (w * sigma * val_dz) as f64) as f32,
                        ),
                    ),
                );
            }
            // Calculate d[n]'s. Keep d[0]/beta[0] at the far boundary, as in C.
            let num6 = cmul(
                cadd(cmplx(1.0, 0.0), cmul(cmplx(-2.0, 0.0), gamma[0])),
                val_ifr1[((0) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                    + (k) as usize],
            );
            let num7 = cmul(
                gamma[0],
                val_ifr1[((1) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                    + (k) as usize],
            );
            d[0] = cadd(num6, num7);
            ohm[0] = cmul(cmplx(-1.0, 0.0), cdiv(d[0], beta[0]));
            for n in 1..last {
                let num3 = cmul(
                    gamma[n],
                    val_ifr1[((n - 1) as usize * val_nsx as usize + (i) as usize)
                        * val_nw as usize
                        + (k) as usize],
                );
                let num4 = cmul(
                    cadd(cmplx(1.0, 0.0), cmul(cmplx(-2.0, 0.0), gamma[n])),
                    val_ifr1[((n) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                        + (k) as usize],
                );
                let num5 = cmul(
                    gamma[n],
                    val_ifr1[((n + 1) as usize * val_nsx as usize + (i) as usize)
                        * val_nw as usize
                        + (k) as usize],
                );
                d[n] = cadd(num3, cadd(num4, num5));
                ohm[n] = cmul(cmplx(-1.0, 0.0), cdiv(d[n], beta[n]));
            }
            let num8 = cmul(
                gamma[last],
                val_ifr1[((last - 1) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                    + (k) as usize],
            );
            let num9 = cmul(
                cadd(cmplx(1.0, 0.0), cmul(cmplx(-2.0, 0.0), gamma[last])),
                val_ifr1[((last) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                    + (k) as usize],
            );
            d[last] = cadd(num8, num9);
            ohm[last] = cmul(cmplx(-1.0, 0.0), cdiv(d[0], beta[0]));
            // Calculate e's and f's.
            e[1] = cadd(cdiv(cmplx(1.0, 0.0), beta[1]), cmplx(2.0, 0.0));
            f[1] = ohm[1];
            for m in 2..last {
                let num1 = cadd(cdiv(cmplx(1.0, 0.0), beta[m]), cmplx(2.0, 0.0));
                let num2 = cdiv(cmplx(-1.0, 0.0), e[m - 1]);
                e[m] = cadd(num1, num2);
                f[m] = cadd(ohm[m], cdiv(f[m - 1], e[m - 1]));
            }
            e[0] = cmplx(0.0, 0.0);
            e[last] = cmplx(0.0, 0.0);
            // Downward continue one frequency on one shot gather.
            val_ifr2[((last) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                + (k) as usize] = cmplx(0.0, 0.0);
            let num10 = cdiv(f[last - 1], e[last - 1]);
            val_ifr2[((last - 1) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                + (k) as usize] = cmul(cmplx(-1.0, 0.0), num10);
            for j in (1..val_nx_pad - 2).rev() {
                let j = j as usize;
                let num11 = cadd(
                    val_ifr2[((j + 1) as usize * val_nsx as usize + (i) as usize)
                        * val_nw as usize
                        + (k) as usize],
                    cmul(cmplx(-1.0, 0.0), f[j]),
                );
                val_ifr2[((j) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                    + (k) as usize] = cdiv(num11, e[j]);
            }
            val_ifr2[((0) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                + (k) as usize] = cmplx(0.0, 0.0);
        }
    }
}

fn padsub(
    val_ifr: &Volume,
    val_ifr_pad: &mut Volume,
    val_nsx: i32,
    val_nx: i32,
    val_nw: i32,
    val_nx_pad: i32,
    val_buff: i32,
) {
    for i in 0..val_nsx as usize {
        for j in 0..val_nx_pad {
            for k in 0..val_nw as usize {
                if j < val_nx + val_buff && j > -1 + val_buff {
                    val_ifr_pad[((j as usize) as usize * val_nsx as usize + (i) as usize)
                        * val_nw as usize
                        + (k) as usize] = val_ifr[(((j - val_buff) as usize) as usize
                        * val_nsx as usize
                        + (i) as usize)
                        * val_nw as usize
                        + (k) as usize];
                } else {
                    val_ifr_pad[((j as usize) as usize * val_nsx as usize + (i) as usize)
                        * val_nw as usize
                        + (k) as usize] = cmplx(0.0, 0.0);
                }
            }
        }
    }
}

fn unpadsub(
    val_ifr: &Volume,
    val_ifr_upad: &mut Volume,
    val_nsx: i32,
    val_nx: i32,
    val_nw: i32,
    _val_nx_pad: i32,
    val_buff: i32,
) {
    for i in 0..val_nsx as usize {
        for j in 0..val_nx as usize {
            for k in 0..val_nw as usize {
                val_ifr_upad[((j) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                    + (k) as usize] = val_ifr[(((j as i32 + val_buff) as usize) as usize
                    * val_nsx as usize
                    + (i) as usize)
                    * val_nw as usize
                    + (k) as usize];
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn taper(
    val_ifr_pad: &Volume,
    val_ifr_pt: &mut Volume,
    val_nsx: i32,
    val_nx_pad: i32,
    val_nw: i32,
    _val_dw: f32,
    val_tap_len: i32,
    val_buff: i32,
    val_nx: i32,
) {
    let pi = 3.1415927f32;
    // Negative tap_len leaves the C allocation untouched (zero in tested runs).
    let mut tap_func = vec![Complex::default(); val_nx_pad as usize];
    if val_tap_len == 0 {
        for m in 0..val_nx_pad as usize {
            tap_func[m] = cmplx(1.0, 0.0);
        }
    }
    if val_tap_len > 0 {
        for m in 0..val_nx_pad {
            if m < val_buff {
                tap_func[m as usize] = cmplx(0.0, 0.0);
            }
            if m >= val_buff && m < val_buff + val_tap_len - 1 {
                tap_func[m as usize] = cmplx(
                    (2.0 * pi as f64 * (m + 1 - val_buff) as f64 / (4.0 * val_tap_len as f64)).sin()
                        as f32,
                    0.0,
                );
            }
            if m >= val_buff + val_tap_len - 1 && m <= val_buff + val_nx - 1 - (val_tap_len - 1) {
                tap_func[m as usize] = cmplx(1.0, 0.0);
            }
            if m > val_buff + val_nx - 1 - (val_tap_len - 1) && m <= val_buff + val_nx - 1 {
                tap_func[m as usize] = cmplx(
                    (2.0 * pi as f64 * (m - (val_buff + val_nx - 1 - (val_tap_len - 1))) as f64
                        / (4.0 * val_tap_len as f64))
                        .cos() as f32,
                    0.0,
                );
            }
            if m > val_buff + val_nx - 1 {
                tap_func[m as usize] = cmplx(0.0, 0.0);
            }
        }
    }
    for i in 0..val_nsx as usize {
        for j in 0..val_nx_pad as usize {
            for k in 0..val_nw as usize {
                val_ifr_pt[((j) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                    + (k) as usize] = cmul(
                    tap_func[j],
                    val_ifr_pad[((j) as usize * val_nsx as usize + (i) as usize) * val_nw as usize
                        + (k) as usize],
                );
            }
        }
    }
}
