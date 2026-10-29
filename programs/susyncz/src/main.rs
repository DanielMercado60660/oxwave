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
//! SUSYNCZ - synthetic seismograms; port of synthetics_waveforms_testpatterns/susyncz.c (Revision 1.16).
//! Credits: Brian Sumner (Fortran); Stockwell and Cohen (C translation).

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUSYNCZ - SYNthetic seismograms for piecewise constant V(Z) function\t",
    "\t   True amplitude (primaries only) modeling for 2.5D\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  susyncz > outfile [parameters]\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " ninf=4        number of interfaces (not including upper surface)\t",
    " dip=5*i       dips of interfaces in degrees (i=1,2,3,4)\t\t",
    " zint=100*i    z-intercepts of interfaces at x=0 (i=1,2,3,4)\t\t",
    " v=1500+ 500*i velocities below surface & interfaces (i=0,1,2,3,4)\t",
    " rho=1,1,1,1,1 densities below surface & interfaces (i=0,1,2,3,4)\t",
    " nline=1\tnumber of (identical) lines\t\t\t\t",
    " ntr=32        number of traces\t\t\t\t\t",
    " dx=10         trace interval\t\t\t\t\t\t",
    " tdelay=0      delay in recording time after source initiation\t\t",
    " dt=0.004      time interval\t\t\t\t\t\t",
    " nt=128        number of time samples\t\t\t\t\t",
    " \t\t\t \t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The original purpose of this code was to create some nontrivial\t",
    " data for Brian Sumner's CZ suite.\t\t\t\t\t",
    " \t\t\t \t\t\t\t\t\t",
    " The program produces zero-offset data over dipping reflectors.\t",
    " \t\t\t \t\t\t\t\t\t",
    " In the original fortran code, some arrays had the index\t\t",
    " interval 1:ninf, as a natural way to index over the subsurface\t",
    " reflectors.  This indexing was preserved in this C translation.\t",
    " Consequently, some arrays in the code do not use the 0 \"slot\".\t",
    " \t\t\t \t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    "\tsusyncz | sufilter | sugain tpow=1 | display_program\t\t",
    " \t\t\t \t\t\t\t\t\t",
    " Trace header fields set: tracl, ns, dt, delrt, ntr, sx, gx\t\t",
];

use oxwave_core::prelude::*;
use oxwave_sinc_interp::ints8r;
use std::f64::consts::PI;

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut nline, mut ntr, mut nt, mut ninf) = (1, 32, 128, 4);
    let (mut dx, mut dt, mut tdelay) = (10.0f32, 0.004f32, 0.0f32);
    getparint("nline", &mut nline);
    getparint("ntr", &mut ntr);
    tr.ntr = ntr.wrapping_mul(nline);
    getparfloat("dx", &mut dx);
    getparint("nt", &mut nt);
    CHECK_NT("nt", nt);
    if nt < 0 {
        err!("invalid nt");
    }
    tr.ns = nt as u16;
    getparfloat("dt", &mut dt);
    tr.dt = (dt * 1000000.0f32) as u32 as u16;
    getparfloat("tdelay", &mut tdelay);
    tr.delrt = (tdelay * 1000.0f32) as i32 as i16;
    getparint("ninf", &mut ninf);
    if ninf < 0 || ninf == i32::MAX {
        err!("invalid ninf");
    }
    let n = (ninf + 1) as usize;
    let mut data = vec![0.0f32; nt as usize];
    let mut v = vec![0.0f32; n];
    let mut rho = vec![0.0f32; n];
    let mut zint = vec![0.0f32; n];
    let mut dip = vec![0.0f32; n];
    let mut meet = vec![0.0f32; n];
    let mut xl = vec![0.0f32; n];
    let mut xr = vec![0.0f32; n];
    let mut tc = vec![0.0f32; n];
    let mut theta = vec![vec![0.0f32; n]; n];
    let mut m = theta.clone();
    let mut b = theta.clone();
    let mut d = theta.clone();
    let mut w = theta.clone();
    let mut k = theta.clone();
    let mut dorow = vec![false; n];
    let mut tout = vec![0.0f32; nt as usize];
    for i in 0..nt as usize {
        tout[i] = i as f32 * dt;
    }
    if getparfloat("zint", &mut zint[1..]) == 0 {
        for i in 1..n {
            zint[i] = (100.0 * i as f64) as f32;
        }
    }
    if getparfloat("dip", &mut dip[1..]) == 0 {
        for i in 1..n {
            dip[i] = (5.0 * i as f64) as f32;
        }
    }
    if getparfloat("v", &mut v[..]) == 0 {
        for i in 0..n {
            v[i] = (1500.0 + 500.0 * i as f64) as f32;
        }
    }
    if getparfloat("rho", &mut rho[..]) == 0 {
        rho.fill(1.0);
    }
    for i in 1..n {
        if dip[i] < -90.0 || dip[i] > 90.0 {
            err!("dips1");
        }
        if (dip[i] - dip[i - 1]).abs() > 90.0 {
            err!("dips2");
        }
        dip[i] = (dip[i] as f64 * (PI / 180.0)) as f32;
        if zint[i] <= zint[i - 1] {
            err!("intercept");
        }
    }
    for i in 0..n - 1 {
        if v[i] <= 0.0 {
            err!("velocities");
        }
        if rho[i] <= 0.0 {
            err!("density");
        }
    }
    for j in 1..n - 1 {
        theta[j][j - 1] = -dip[j];
        for i in (1..j).rev() {
            let temp =
                (v[i - 1] as f64 * ((theta[j][i] + dip[i]) as f64).sin() / v[i] as f64) as f32;
            if temp.abs() > 1.0 {
                err!("specular %d", j as i32);
            }
            theta[j][i - 1] = ((temp as f64).asin() - dip[i] as f64) as f32;
        }
    }
    for j in 1..n {
        for i in 0..j {
            let s = (theta[j][i] as f64).sin() as f32;
            let temp = (1.0 - (dip[i + 1] as f64).tan() * s as f64) as f32;
            m[j][i] = ((1.0 - (dip[i] as f64).tan() * s as f64) / temp as f64) as f32;
            b[j][i] = (zint[i + 1] - zint[i]) * s / temp;
        }
    }
    let mut xmin = 0.0f32;
    let mut xmax = dx * ntr.wrapping_sub(1) as f32;
    for j in 1..n {
        xl[j] = xmin;
        xr[j] = xmax;
    }
    // The original second inner loop increments j rather than j2.
    let mut j = 0usize;
    while j < n - 1 {
        for j2 in j + 1..n {
            if dip[j2] != dip[j] {
                let intercept = ((zint[j2] - zint[j]) as f64
                    / ((dip[j] as f64).tan() - (dip[j2] as f64).tan()))
                    as f32;
                if intercept > xmin && intercept < xmax {
                    err!("intercept");
                }
            }
        }
        let j2 = j + 1;
        while j < n {
            xl[j2] = m[j2][j] * xl[j2] + b[j2][j];
            xr[j2] = m[j2][j] * xr[j2] + b[j2][j];
            if xl[j2] < xmin {
                xmin = xl[j2];
            }
            if xr[j2] > xmax {
                xmax = xr[j2];
            }
            j += 1;
        }
        j += 1;
    }
    for i in 0..n - 1 {
        if dip[i + 1] == dip[i] {
            for j in i + 1..n {
                d[j][i] = ((zint[i + 1] - zint[i]) as f64 * (dip[i] as f64).cos()
                    / ((theta[j][i] + dip[i]) as f64).cos()) as f32;
            }
        } else {
            dorow[i] = true;
            meet[i] = ((zint[i + 1] - zint[i]) as f64
                / ((dip[i] as f64).tan() - (dip[i + 1] as f64).tan())) as f32;
            for j in i + 1..n {
                k[j][i] = (((dip[i] - dip[i + 1]) as f64).sin()
                    / (((theta[j][i] + dip[i + 1]) as f64).cos() * (dip[i] as f64).cos()))
                    as f32;
            }
        }
    }
    tabtrcoefs(n, &rho, &v, &theta, &dip, &mut tc);
    for j in 1..n {
        w[j][0] = 1.0;
        for i in 1..j {
            let c1 = ((theta[j][i - 1] + dip[i]) as f64).cos() as f32;
            let c2 = ((theta[j][i] + dip[i]) as f64).cos() as f32;
            w[j][i] = w[j][i - 1] * c1 * c1 / (c2 * c2);
        }
    }
    for _ in 0..nline {
        let mut x = 0.0f32;
        for itr in 0..ntr {
            tr.data[..nt as usize].fill(0.0);
            data.fill(0.0);
            for j in 1..n {
                let mut ex = x;
                let mut dist = if dorow[0] {
                    (meet[0] - ex) * k[j][0]
                } else {
                    d[j][0]
                };
                let mut t0 = dist / v[0];
                let mut tmp = dist * v[0];
                let mut p1 = tmp;
                let mut p2 = tmp;
                for i in 1..j {
                    ex = m[j][i - 1] * ex + b[j][i - 1];
                    dist = if dorow[i] {
                        (meet[i] - ex) * k[j][i]
                    } else {
                        d[j][i]
                    };
                    t0 += dist / v[i];
                    tmp = dist * v[i];
                    p1 += tmp;
                    p2 += tmp * w[j][i];
                }
                let amp = [(v[0] as f64 * tc[j] as f64
                    / (dt as f64 * 4.0 * PI * 2.0 * ((p1 * p2) as f64).sqrt()))
                    as f32];
                let time = (2.0 * t0 as f64 - tdelay as f64) as f32;
                let q = time / dt;
                let itime = (q as f64 + if q >= 0.0 { 0.5 } else { -0.5 }) as i32;
                if itime >= 0 && itime < nt {
                    data[itime as usize] = amp[0];
                }
                ints8r(dt, time, &amp, 0.0, 0.0, &tout, &mut data);
                for it in 0..nt as usize {
                    tr.data[it] += data[it];
                }
            }
            tr.tracl = itr + 1;
            tr.gx = x as i32;
            tr.sx = x as i32;
            puttr(&mut tr);
            x += dx;
        }
    }
    exit(CWP_Exit());
}
fn tabtrcoefs(n: usize, rho: &[f32], v: &[f32], theta: &[Vec<f32>], dip: &[f32], tc: &mut [f32]) {
    let mut temp = vec![0.0f32; n];
    let mut tmp = vec![vec![0.0f32; n]; n];
    for i in 0..n {
        temp[i] = rho[i] * v[i];
    }
    for j in 1..n {
        for i in 1..j {
            let t1 = (temp[i] as f64 * ((theta[j][i - 1] + dip[i]) as f64).cos()) as f32;
            let t2 = (temp[i - 1] as f64 * ((theta[j][i] + dip[i]) as f64).cos()) as f32;
            let r = (t1 - t2) / (t1 + t2);
            tmp[j][i] = (1.0 - (r * r) as f64) as f32;
        }
    }
    for j in 1..n {
        let t1 = temp[j];
        let t2 = temp[j - 1];
        let mut p = (t1 - t2) / (t1 + t2);
        for i in 1..j {
            p *= tmp[j][i];
        }
        tc[j] = p;
    }
}
