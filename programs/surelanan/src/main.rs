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

//! SURELANAN - residual-moveout semblance analysis for anisotropic media.
//! Port of su/main/velocity_analysis/surelanan.c.
//! The C selfdoc credits Debashish Sarkar's anisotropic model building technique.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_f_cubicspline::{cakima, cmonot, csplin};
use oxwave_linear_interp::intlin;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SURELANAN - REsiduaL-moveout semblance ANalysis for ANisotropic media\t",
    "\t\t\t\t\t\t\t\t\t",
    " surelan refl= npicks=    [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " reflector file: reflec =\t\t\t\t\t\t",
    " number of points in the reflector file =\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nr1=51\t\tnumber of r1-parameter samples   \t\t",
    " dr1=0.01              r1-parameter sampling interval\t\t\t",
    " fr1=-0.25             first value of r1-parameter\t\t\t",
    " nr2=51\t\tnumber of r2-parameter samples   \t\t",
    " dr2=0.01              r2-parameter sampling interval\t\t\t",
    " fr2=-0.25             first value of r2-parameter\t\t\t",
    " dzratio=5             ratio of output to input depth sampling intervals",
    " nsmooth=dzratio*2+1   length of semblance num and den smoothing window",
    " verbose=0             =1 for diagnostic print on stderr\t\t",
    " method=linear\t\tfor linear interpolation of the interface       ",
    " \t\t\t=mono for monotonic cubic interpolation of interface",
    " \t\t\t=akima for Akima's cubic interpolation of interface ",
    " \t\t\t=spline for cubic spline interpolation of interface ",
    "\t\t\t\t\t\t\t\t\t",
    " Note: \t\t\t\t\t\t\t\t",
    " 1. This program is part of Debashish Sarkar's anisotropic model building",
    " technique. \t\t\t\t\t\t\t\t",
    " 2. Input migrated traces should be sorted by cdp - surelan outputs a \t",
    "    group of semblance traces every time cdp changes.  Therefore, the  ",
    "    output will be useful only if cdp gathers are input.  \t\t",
    " 3. The residual-moveout semblance for cdp gathers is based\t\t",
    "\ton z(h)*z(h) = z(0)*z(0) + r1*h^2 + r2*h^4/[h^2+z(0)^2] where z ",
    "\tdepth and h is the half-offset.   \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nz = tr.ns as usize;
    let dz = tr.d1;
    let fz = tr.f1;
    let mut cdp = tr.cdp as i64;
    let mut half_offset = (tr.offset / 2) as f32;
    let mut refl = String::new();
    let mut npicks = 0i32;
    getparstring("refl", &mut refl);
    getparint("npicks", &mut npicks);
    let (mut nr1, mut nr2) = (51i32, 51i32);
    let (mut dr1, mut dr2, mut fr1, mut fr2) = (0.01f32, 0.01f32, -0.25f32, -0.25f32);
    let mut dzratio = 5i32;
    let mut nsmooth = 0i32;
    let mut verbose = 0i32;
    let mut method = "linear".to_string();
    getparint("nr1", &mut nr1);
    getparint("nr2", &mut nr2);
    getparfloat("dr1", &mut dr1);
    getparfloat("dr2", &mut dr2);
    getparfloat("fr1", &mut fr1);
    getparfloat("fr2", &mut fr2);
    getparint("dzratio", &mut dzratio);
    if getparint("nsmooth", &mut nsmooth) == 0 {
        nsmooth = dzratio.wrapping_mul(2).wrapping_add(1);
    }
    getparint("verbose", &mut verbose);
    getparstring("method", &mut method);
    checkpars();
    // Reject empty pick arrays before interpolation, as in the UB-fixed C.
    if npicks <= 0 {
        err!("npicks must be positive");
    }
    let fp = fopen(&refl, "r").expect("reflector file could not be opened");
    // Allocate before reading the picks, as in the C program.
    let mut data = vec![0.0f32; nz];
    let size = (nr1 as usize)
        .wrapping_mul(nr2 as usize)
        .wrapping_mul(nsmooth as usize);
    let mut num = vec![0.0f32; size];
    let mut den = vec![0.0f32; size];
    let mut nnz = vec![0.0f32; size];
    let mut sem = vec![0.0f32; nr2 as usize];
    let mut zint = vec![0.0f32; npicks as usize];
    let mut xint = vec![0.0f32; npicks as usize];
    let xcdp = [cdp as f32];
    let mut zcdp = [0.0f32];

    // Input picked interface. The C fscanf has two conversions; a failed
    // first conversion prevents attempting the second one in that call.
    for i in 0..npicks as usize {
        if let Some(v) = fscanf_f32(&fp) {
            zint[i] = v;
            if let Some(v) = fscanf_f32(&fp) {
                xint[i] = v;
            }
        }
    }
    match method.as_bytes().first().copied() {
        Some(b'l') => intlin(
            &xint,
            &zint,
            zint[0],
            zint[npicks as usize - 1],
            &xcdp,
            &mut zcdp,
        ),
        Some(b'm') => {
            let mut zind = vec![[0.0f32; 4]; npicks as usize];
            cmonot(&xint, &zint, &mut zind);
            intcub(0, &xint, &zind, &xcdp, &mut zcdp);
        }
        Some(b'a') => {
            let mut zind = vec![[0.0f32; 4]; npicks as usize];
            cakima(&xint, &zint, &mut zind);
            intcub(0, &xint, &zind, &xcdp, &mut zcdp);
        }
        Some(b's') => {
            let mut zind = vec![[0.0f32; 4]; npicks as usize];
            csplin(&xint, &zint, &mut zind);
            intcub(0, &xint, &zind, &xcdp, &mut zcdp);
        }
        _ => err!("%s is an unknown interpolation method!\n", method),
    }
    let z = zcdp[0];
    warn!("%f \n", z);
    let mut gottrace = true;
    let mut cdpprev = tr.cdp as i64;
    loop {
        if gottrace {
            half_offset = (tr.offset / 2) as f32;
            cdp = tr.cdp as i64;
            data.copy_from_slice(&tr.data[..nz]);
        }
        if cdp != cdpprev || !gottrace {
            tr.offset = 0;
            tr.cdp = cdpprev as i32;
            tr.ns = nr2 as u16;
            tr.d1 = dr2;
            for ir1 in 0..nr1 as usize {
                for ir2 in 0..nr2 as usize {
                    let (mut nsum, mut dsum) = (0.0f32, 0.0f32);
                    for is in 0..nsmooth as usize {
                        let k = (ir1 * nr2 as usize + ir2) * nsmooth as usize + is;
                        nsum += num[k] * num[k];
                        dsum += nnz[k] * den[k];
                    }
                    sem[ir2] = if dsum != 0.0 { nsum / dsum } else { 0.0 };
                }
                tr.data[..nr2 as usize].copy_from_slice(&sem);
                puttr(&mut tr);
            }
            if verbose != 0 {
                warn!("tsemblance output for cdp=%d\n", cdpprev);
            }
            if !gottrace {
                break;
            }
            // The original retains the first CDP's reflector depth and
            // all accumulated sums when advancing to another gather.
            cdpprev = cdp;
        }
        let mut r1 = fr1;
        for ir1 in 0..nr1 as usize {
            let mut r2 = fr2;
            for ir2 in 0..nr2 as usize {
                let roffs2 = ((r1 * half_offset * half_offset) as f64
                    + (r2 as f64) * (half_offset as f64).powf(4.0)
                        / ((half_offset as f64).powf(2.0) + (z as f64).powf(2.0)))
                    as f32;
                let zmin = z - (nsmooth.wrapping_sub(1) as f32 * dz / 2.0f32);
                let mut zn = zmin;
                for iz in 0..nsmooth as usize {
                    let mut temp = zn * zn + roffs2;
                    let zi = if temp > fz * fz {
                        (((temp as f64).sqrt() - fz as f64) / dz as f64) as f32
                    } else {
                        0.0
                    };
                    let izi = zi as i32;
                    if izi >= 0 && izi < nz as i32 - 1 {
                        let frac = zi - izi as f32;
                        temp = ((1.0 - frac as f64) * data[izi as usize] as f64
                            + (frac * data[izi as usize + 1]) as f64)
                            as f32;
                        if temp != 0.0 {
                            let k = (ir1 * nr2 as usize + ir2) * nsmooth as usize + iz;
                            num[k] += temp;
                            den[k] += temp * temp;
                            nnz[k] += 1.0;
                        }
                    }
                    zn += dz;
                }
                r2 += dr2;
            }
            r1 += dr1;
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    exit(CWP_Exit());
}
