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
//! SUVEL2DF - compute two-velocity semblance at one time.
//! Port of su/main/velocity_analysis/suvel2df.c (Revision 1.4).
//! Credits: CWP: Tariq Alkhalifah, February 1997.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.

#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUVEL2DF - compute stacking VELocity semblance for a single time in   ",
    "\t\t\t    over Vnmo and eta in 2-D\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "    suvel2df <stdin >stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " tn\t\t\tzero-offset time of reflection\t\t\t",
    " offsetm\t\tMaximum offset considerd\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nv=50\t\t\tnumber of velocities\t\t\t\t",
    " dv=50.0\t\tvelocity sampling interval\t\t\t",
    " fv=1500.0\t\tfirst velocity\t\t\t\t\t",
    " nvh=50\t\tnumber of horizotal velocities\t\t\t",
    " dvh=50.0\t\thorizontal velocity sampling interval\t\t",
    " fvh=1500.0\t\tfirst horizontal velocity\t\t\t",
    " xod=1.5\t\tmaximum offset-to-depth ratio to resolve\t",
    " dtratio=5\t\tratio of output to input time sampling intervals",
    " nsmooth=dtratio*2+1\tlength of semblance num and den smoothing window",
    " verbose=0\t\t=1 for diagnostic print on stderr\t\t",
    " vavg=fv+0.5*(nv-1)*dv   average velocity used in the search\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Semblance is defined by the following quotient:\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t n-1\t\t \t\t\t\t\t",
    "\t\t[ sum q(t,j) ]^2\t\t\t\t\t",
    "\t\t j=0\t\t \t\t\t\t\t",
    "\ts(t) = ------------------\t\t\t\t\t",
    "\t\t n-1\t\t \t\t\t\t\t",
    "\t\tn sum [q(t,j)]^2\t\t\t\t\t",
    "\t\t j=0\t\t \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " where n is the number of non-zero samples after muting.\t\t",
    " Smoothing (nsmooth) is applied separately to the numerator and denominator",
    " before computing this semblance quotient.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Input traces should be sorted by cdp - suvel2df outputs a group of\t",
    " semblance traces every time cdp changes.  Therefore, the output will\t",
    " be useful only if cdp gathers are input.\t\t\t\t",
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
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut tnn = 0.0f32;
    let mut offsetm = 0.0f32;
    if getparfloat("tn", &mut tnn) == 0 {
        err!("must specify tn!");
    }
    if getparfloat("offsetm", &mut offsetm) == 0 {
        err!("must specify offsetm!");
    }
    let mut nv = 50i32;
    let mut dv = 50.0f32;
    let mut fv = 1500.0f32;
    let mut nvh = 50i32;
    let mut dvh = 50.0f32;
    let mut fvh = 1500.0f32;
    let mut xod = 1.5f32;
    let mut dtratio = 5i32;
    let mut nsmooth = 0i32;
    let mut verbose = 0i32;
    let mut vavg = 0.0f32;
    let mut offsetmax = 0.0f32;
    getparint("nv", &mut nv);
    getparfloat("dv", &mut dv);
    getparfloat("fv", &mut fv);
    getparint("nvh", &mut nvh);
    getparfloat("dvh", &mut dvh);
    getparfloat("fvh", &mut fvh);
    getparfloat("xod", &mut xod);
    getparint("dtratio", &mut dtratio);
    if getparint("nsmooth", &mut nsmooth) == 0 {
        nsmooth = dtratio.wrapping_mul(2).wrapping_add(1);
    }
    getparint("verbose", &mut verbose);
    if getparfloat("vavg", &mut vavg) == 0 {
        vavg = (fv as f64 + 0.5 * (nv - 1) as f64 * dv as f64) as f32;
    }
    if getparfloat("offsetmax", &mut offsetmax) == 0 {
        offsetmax = (0.5 * tnn as f64 * xod as f64 * vavg as f64) as f32;
    }
    checkpars();
    let size = (nv as usize)
        .wrapping_mul(nvh as usize)
        .wrapping_mul(nsmooth as usize);
    let mut num = vec![0.0f32; size];
    let mut den = vec![0.0f32; size];
    let mut nnz = vec![0.0f32; size];
    let mut data = vec![0.0f32; nt as usize];
    let mut sem = vec![0.0f32; nvh as usize];
    let mut gottrace = true;
    let mut cdp = tr.cdp;
    let mut cdpprev = cdp;
    let mut offset = tr.offset as f32;
    let mut max = 0.0f32;
    let mut vmax = fv;
    let mut vhmax = fvh;
    loop {
        if gottrace {
            offset = tr.offset as f32;
            cdp = tr.cdp;
            data.copy_from_slice(&tr.data[..nt as usize]);
        }
        if cdp != cdpprev || !gottrace {
            tr.offset = 0;
            tr.cdp = cdpprev;
            tr.ns = nvh as u16;
            tr.dt = dvh as u16;
            for iv in 0..nv {
                for ivh in 0..nvh {
                    let mut nsum = 0.0f32;
                    let mut dsum = 0.0f32;
                    for is in 0..nsmooth {
                        let idx = ((iv * nvh + ivh) * nsmooth + is) as usize;
                        nsum += num[idx] * num[idx];
                        dsum += nnz[idx] * den[idx];
                    }
                    sem[ivh as usize] = if dsum != 0.0 { nsum / dsum } else { 0.0 };
                    if sem[ivh as usize] > max {
                        max = sem[ivh as usize];
                        vhmax = fvh + ivh as f32 * dvh;
                        vmax = fv + iv as f32 * dv;
                    }
                }
                tr.data[..nvh as usize].copy_from_slice(&sem);
                puttr(&mut tr);
                for ivh in 0..nvh {
                    for it in 0..nsmooth {
                        let idx = ((iv * nvh + ivh) * nsmooth + it) as usize;
                        num[idx] = 0.0;
                        den[idx] = 0.0;
                        nnz[idx] = 0.0;
                    }
                }
            }
            warn!("max=%f vmax=%f vhmax=%f\n", max, vmax, vhmax);
            if verbose != 0 {
                warn!("semblance output for cdp=%ld\n", cdpprev as i64);
            }
            if !gottrace {
                break;
            }
            cdpprev = cdp;
        }
        if (if offset < 0.0 { -offset } else { offset }) > offsetmax {
            if gettr(&mut tr) == 0 {
                gottrace = false;
            }
            continue;
        }
        let offset2 = offset * offset;
        let mut v = fv;
        for iv in 0..nv {
            let mut vh = fvh;
            for ivh in 0..nvh {
                let v2 = v * v;
                let v4 = v2 * v2;
                let vh2 = vh * vh;
                let offovs = offset2 / v2;
                let fac1 = (vh2 - v2) * offset2 * offset2 / v2;
                let fac2 = vh2 * offset2;
                let mut ismin = (tnn / dt - (nsmooth / 2) as f32) as i32;
                let mut ismax = ismin + nsmooth;
                if ismin < 0 {
                    ismin = 0;
                }
                if ismax > nt - 1 {
                    ismax = nt - 1;
                }
                let mut tn = ft + ismin as f32 * dt;
                let mut is = 0;
                for _it in ismin..ismax {
                    let radicand = if offset == 0.0 {
                        tn * tn + offovs
                    } else {
                        tn * tn + offovs - fac1 / (tn * tn * v4 + fac2)
                    };
                    let ti = (((radicand as f64).sqrt() - ft as f64) / dt as f64) as f32;
                    let iti = ti as i32;
                    if iti < nt - 1 {
                        let frac = ti - iti as f32;
                        let temp = ((1.0 - frac as f64) * data[iti as usize] as f64
                            + (frac * data[(iti + 1) as usize]) as f64)
                            as f32;
                        if temp != 0.0 {
                            let idx = ((iv * nvh + ivh) * nsmooth + is) as usize;
                            num[idx] += temp;
                            den[idx] += temp * temp;
                            nnz[idx] += 1.0;
                        }
                    }
                    is += 1;
                    tn += dt;
                }
                vh += dvh;
            }
            v += dv;
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    exit(CWP_Exit());
}
