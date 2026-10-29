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

//! SUSMGAUSS2 - Gaussian smoothing of a two-dimensional velocity array.
//! Port of su/main/filters/susmgauss2.c (Revision 1.2).
//! Credits: CWP: Carlos Pacheco, 2005.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{Complex, crmul};
use oxwave_pfa_fft::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSMGAUSS2 --- SMOOTH a uniformly sampled 2d array of velocities\t",
    "\t\tusing a Gaussian filter specified with correlation \t",
    "\tlengths a1 and a2.\t",
    "\t\t\t\t\t\t\t\t\t",
    " susmgauss2 < stdin [optional parameters ] > stdout\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " a1=0\t\t\tsmoothing parameter in the 1 direction\t\t",
    " a2=0\t\t\tsmoothing parameter in the 2 direction\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Larger a1 and a2 result in a smoother velocity. The velocities are\t",
    " first transformed to slowness and then a Gaussian filter is applied\t",
    " in the wavenumber domain.\t\t\t\t\t\t",
    "\t",
    " Input file must be in SU format. The output file is smoothed velocity ",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut a1 = 0.0f32;
    let mut a2 = 0.0f32;
    getparfloat("a1", &mut a1);
    getparfloat("a2", &mut a2);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    if tr.trid != 130 { warn!("tr.trid=%d", tr.trid); }
    let nx1 = tr.ns as usize;
    let mut headers = Vec::new();
    let mut v = Vec::new();
    loop {
        headers.push(tr.hdr);
        v.extend_from_slice(&tr.data[..nx1]);
        if gettr(&mut tr) == 0 { break; }
    }
    let nx2 = headers.len();
    let nx1fft = npfaro(nx1 as i32, 2 * nx1 as i32);
    let nx2fft = npfa(nx2 as i32);
    if nx1fft >= SU_NFLTS as i32 || nx1fft >= 720720 { err!("Padded nx1=%d--too big", nx1fft); }
    if nx2fft >= 720720 { err!("Padded nx2=%d--too big", nx2fft); }
    let n1 = nx1fft as usize;
    let n2 = nx2fft as usize;
    let nk1 = n1 / 2 + 1;
    let nk2 = n2 / 2 + 1;
    let mut rt = vec![0.0f32; n1 * n2];
    let mut ct = vec![Complex::default(); nk1 * n2];
    let mut kfilter = vec![0.0f32; n1 * n2];
    let mut kfilt1 = vec![0.0f32; nk1];
    let mut kfilt2 = vec![0.0f32; nk2];
    let pi = std::f64::consts::PI as f32;
    let dk1 = 2.0f32 * pi / nx1fft as f32;
    let dk2 = 2.0f32 * pi / nx2fft as f32;
    for ik1 in 0..nk1 {
        let c1 = a1 * ik1 as f32 * dk1 / 2.0f32;
        kfilt1[ik1] = (-(c1 as f64).powf(2.0)).exp() as f32;
    }
    for ik2 in 0..nk2 {
        let c2 = a2 * ik2 as f32 * dk2 / 2.0f32;
        kfilt2[ik2] = (-(c2 as f64).powf(2.0)).exp() as f32;
    }
    for ik2 in 0..nk2 {
        for ik1 in 0..nk1 { kfilter[ik2*n1+ik1] = kfilt2[ik2]*kfilt1[ik1]; }
    }
    for ik2 in nk2..n2 {
        for ik1 in 0..nk1 { kfilter[ik2*n1+ik1] = kfilt2[n2-ik2]*kfilt1[ik1]; }
    }
    let mut vmax = v[0];
    for ix2 in 0..nx2 {
        for ix1 in 0..nx1 {
            let x = v[ix2*nx1+ix1];
            vmax = if vmax > x { vmax } else { x };
        }
    }
    for ix2 in 0..nx2 { rt[ix2*n1..ix2*n1+nx1].copy_from_slice(&v[ix2*nx1..(ix2+1)*nx1]); }
    pfa2rc(-1, 1, nx1fft, nx2 as i32, &rt, &mut ct);
    pfa2cc(-1, 2, nk1 as i32, nx2fft, &mut ct);
    for ik2 in 0..n2 {
        for ik1 in 0..nk1 { ct[ik2*nk1+ik1] = crmul(ct[ik2*nk1+ik1], kfilter[ik2*n1+ik1]); }
    }
    pfa2cc(1, 2, nk1 as i32, nx2fft, &mut ct);
    pfa2cr(1, 1, nx1fft, nx2 as i32, &ct, &mut rt);
    let mut vfmax = rt[0];
    for ix2 in 0..nx2 {
        for ix1 in 0..nx1 {
            let x = rt[ix2*n1+ix1];
            vfmax = if vfmax > x { vfmax } else { x };
        }
    }
    for ix2 in 0..nx2 {
        tr.hdr = headers[ix2];
        for ix1 in 0..nx1 { tr.data[ix1] = rt[ix2*n1+ix1] * vmax / vfmax; }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
