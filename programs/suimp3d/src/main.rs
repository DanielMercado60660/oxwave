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
//! SUIMP3D, port of synthetics_waveforms_testpatterns/suimp3d.c (revision 1.25).
//! Credits: CWP: Norm Bleistein, Jack K. Cohen; UHouston: Chris Liner 2010
//! (added direct arrival option).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{Complex, crmul};
use oxwave_cwp::cwp::nint;
use oxwave_pfa_fft::{npfaro, pfarc, pfacr};
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t",
    "SUIMP3D - generate inplane shot records for a point \t",
    "          scatterer embedded in three dimensions using\t",
    "          the Born integral equation\t\t\t",
    "\t\t\t\t\t\t\t",
    "suimp3d [optional parameters] >stdout \t\t\t",
    "\t\t\t\t\t\t\t",
    "Optional parameters\t\t\t\t\t",
    "\tnshot=1\t\tnumber of shots\t\t\t",
    "\tnrec=1\t\tnumber of receivers\t\t",
    "\tc=5000\t\tspeed\t\t\t\t",
    "\tdt=.004\t\tsampling rate\t\t\t",
    "\tnt=256\t\tnumber of samples\t\t",
    "\tx0=1000\t\tpoint scatterer location\t",
    "\ty0=0\t\tpoint scatterer location\t",
    "\tz0=1000\t\tpoint scatterer location\t",
    "   dir=0\t\tdo not include direct arrival\t",
    "\t            =1 include direct arrival\t",
    "\tsxmin=0\t\tfirst shot location\t\t",
    "\tsymin=0\t\tfirst shot location\t\t",
    "\tszmin=0\t\tfirst shot location\t\t",
    "\tgxmin=0\t\tfirst receiver location\t\t",
    "\tgymin=0\t\tfirst receiver location\t\t",
    "\tgzmin=0\t\tfirst receiver location\t\t",
    "\tdsx=100\t\tx-step in shot location\t\t",
    "\tdsy=0\t \ty-step in shot location\t\t",
    "\tdsz=0\t \tz-step in shot location\t\t",
    "\tdgx=100\t\tx-step in receiver location\t",
    "\tdgy=0\t\ty-step in receiver location\t",
    "\tdgz=0\t\tz-step in receiver location\t",
    "\t\t\t\t\t\t\t",
    " Example:                                              ",
    "       suimp3d nrec=32 | sufilter | supswigp | ...     ",
    "\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args(); initargs(&argv); requestdoc(0, SDOC);
    let (mut nshot, mut nrec, mut nt, mut dir) = (1, 1, 256, 0);
    let (mut c, mut dt, mut x0, mut y0, mut z0) = (5000.0f32, 0.004f32, 1000.0f32, 0.0f32, 1000.0f32);
    let (mut sxmin, mut symin, mut szmin) = (0.0f32, 0.0f32, 0.0f32);
    let (mut gxmin, mut gymin, mut gzmin) = (0.0f32, 0.0f32, 0.0f32);
    let (mut dsx, mut dsy, mut dsz) = (100.0f32, 0.0f32, 0.0f32);
    let (mut dgx, mut dgy, mut dgz) = (100.0f32, 0.0f32, 0.0f32);
    getparint("nshot", &mut nshot); getparint("nrec", &mut nrec); getparint("nt", &mut nt);
    getparfloat("c", &mut c); getparfloat("dt", &mut dt);
    getparfloat("x0", &mut x0); getparfloat("y0", &mut y0); getparfloat("z0", &mut z0); getparint("dir", &mut dir);
    getparfloat("sxmin", &mut sxmin); getparfloat("symin", &mut symin); getparfloat("szmin", &mut szmin);
    getparfloat("gxmin", &mut gxmin); getparfloat("gymin", &mut gymin); getparfloat("gzmin", &mut gzmin);
    getparfloat("dsx", &mut dsx); getparfloat("dsy", &mut dsy); getparfloat("dsz", &mut dsz);
    getparfloat("dgx", &mut dgx); getparfloat("dgy", &mut dgy); getparfloat("dgz", &mut dgz);
    let mut tr = Trace::new(); tr.ns = nt as u16; tr.dt = nint(dt as f64 * 1000000.0) as u16;
    // The C allocation converts a negative nt to an impossible size_t allocation.
    if nt < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut temp = vec![0.0f32; nt as usize];
    let nfft = npfaro(nt, 2 * nt);
    if nfft >= SU_NFLTS as i32 || nfft >= 720720 { err!("Padded nt=%d -- too big", nfft); }
    let nfby2p1 = nfft / 2 + 1;
    let mut rt = vec![0.0f32; nfft as usize]; let mut ct = vec![Complex::default(); nfby2p1 as usize];
    let k = (1.0 / (4.0 * c as f64 * c as f64 * dt as f64 * dt as f64 * dt as f64 * nfft as f64 * nfft as f64 * nfft as f64)) as f32;
    let mut tout = vec![0.0f32; nt as usize];
    for i in 0..nt { tout[i as usize] = i as f32 * dt; }
    let mut tracl = 0;
    for s in 0..nshot {
        let sx = sxmin + s as f32 * dsx; let sy = symin + s as f32 * dsy; let sz = szmin + s as f32 * dsz;
        let rs = (((sx-x0)*(sx-x0) + (sy-y0)*(sy-y0) + (sz-z0)*(sz-z0)) as f64).sqrt() as f32;
        for g in 0..nrec {
            tr.data[..nt as usize].fill(0.0);
            let gx = gxmin + g as f32 * dgx; let gy = gymin + g as f32 * dgy; let gz = gzmin + g as f32 * dgz;
            let rg = (((gx-x0)*(gx-x0) + (gy-y0)*(gy-y0) + (gz-z0)*(gz-z0)) as f64).sqrt() as f32;
            let d = rs + rg; let mut t = d/c; let mut spread = rs*rg; let mut amplitude = [k/spread];
            // The original yields zero for an infinite travel time with finite,
            // nonzero sampling; the sinc crate cannot index that coordinate.
            if t.is_infinite() && dt.is_finite() && dt != 0.0 {
                tr.data[..nt as usize].fill(0.0);
            } else {
                ints8r(dt, t, &amplitude, 0.0, 0.0, &tout, &mut tr.data[..nt as usize]);
            }
            if dir == 1 {
                let rd = (((gx-sx)*(gx-sx) + (gy-sy)*(gy-sy)) as f64).sqrt() as f32;
                t = rd/c; spread = rd; amplitude[0] = k/spread;
                if t.is_infinite() && dt.is_finite() && dt != 0.0 {
                    temp.fill(0.0);
                } else {
                    ints8r(dt, t, &amplitude, 0.0, 0.0, &tout, &mut temp);
                }
                for i in 0..nt as usize { tr.data[i] = tr.data[i] + temp[i]; }
            }
            rt[..nt as usize].copy_from_slice(&tr.data[..nt as usize]); rt[nt as usize..].fill(0.0);
            pfarc(1, &rt, &mut ct);
            for i in 0..nfby2p1 { ct[i as usize] = crmul(ct[i as usize], i.wrapping_mul(i) as f32); }
            pfacr(-1, &ct, &mut rt); tr.data[..nt as usize].copy_from_slice(&rt[..nt as usize]);
            tracl += 1; tr.tracl = tracl; tr.tracr = tracl; tr.fldr = 1+s; tr.tracf = 1+g;
            tr.sx = nint(sx as f64); tr.sy = nint(sy as f64); tr.selev = -nint(sz as f64);
            tr.gx = nint(gx as f64); tr.gy = nint(gy as f64); tr.gelev = -nint(gz as f64);
            tr.offset = (((sx-gx)*(sx-gx) + (sy-gy)*(sy-gy) + (sz-gz)*(sz-gz)) as f64).sqrt() as i32;
            if dgy == 0.0 && dgz == 0.0 { tr.offset = nint((if dsx > 0.0 { gx-sx } else { sx-gx }) as f64); }
            if dgx == 0.0 && dgz == 0.0 { tr.offset = nint((if dsy > 0.0 { gy-sy } else { sy-gy }) as f64); }
            if dgx == 0.0 && dgy == 0.0 { tr.offset = nint((if dsz > 0.0 { gz-sz } else { sz-gz }) as f64); }
            puttr(&mut tr);
        }
    }
    exit(CWP_Exit());
}
