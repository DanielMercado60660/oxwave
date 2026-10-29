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

//! SUNMO_a - NMO for arbitrary velocity and anisotropy.
//! Port of stretching_moveout_resamp/sunmo_a.c, revision 1.5.
//! Credits: SEP: Shuki, Chuck Sword; CWP: Shuki, Jack, Dave Hale, Bjoern Rommel; Carlos E. Theodoro, Sang-yong Suh, John Stockwell.

use oxwave_core::prelude::*;
use oxwave_index_search::xindex;
use oxwave_inverse_interp::yxtoxy;
use oxwave_linear_interp::intlin;
use oxwave_sinc_interp::ints8r;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUNMO_a - NMO for an arbitrary velocity function of time and CDP with\t     ",
    "\t\texperimental Anisotropy options\t\t\t\t     ",
    "  sunmo <stdin >stdout [optional parameters]\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:\t\t\t\t\t\t\t     ",
    " tnmo=0,...\t\tNMO times corresponding to velocities in vnmo\t     ",
    " vnmo=1500,..\t\tNMO velocities corresponding to times in tnmo\t     ",
    " anis1=0\t\ttwo anisotropy coefficients making up quartic term   ",
    " anis2=0\t\tin traveltime curve, corresponding to times in tnmo  ",
    " cdp=\t\t\tCDPs for which vnmo & tnmo are specified (see Notes) ",
    " smute=1.5\t\tsamples with NMO stretch exceeding smute are zeroed  ",
    " lmute=25\t\tlength (in samples) of linear ramp for stretch mute  ",
    " sscale=1\t\t=1 to divide output samples by NMO stretch factor    ",
    " invert=0\t\t=1 to perform (approximate) inverse NMO\t\t     ",
    " upward=0\t\t=1 to scan upward to find first sample to kill\t",
    "\t\t\t\t\t\t\t\t\t     ",
    " Notes:\t\t\t\t\t\t\t\t     ",
    " For constant-velocity NMO, specify only one vnmo=constant and omit tnmo.   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " NMO interpolation error is less than 1% for frequencies less than 60% of   ",
    " the Nyquist frequency.\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Exact inverse NMO is impossible, particularly for early times at large     ",
    " offsets and for frequencies near Nyquist with large interpolation errors.  ",
    " \t\t\t\t\t\t\t\t     \t     ",
    " The \"offset\" header field must be set.\t\t\t\t     ",
    " Use suazimuth to set offset header field when sx,sy,gx,gy are all\t     ",
    " nonzero. \t\t\t\t\t\t\t   \t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " For NMO with a velocity function of time only, specify the arrays\t     ",
    "\t   vnmo=v1,v2,... tnmo=t1,t2,...\t\t\t\t     ",
    " where v1 is the velocity at time t1, v2 is the velocity at time t2, ...    ",
    " The times specified in the tnmo array must be monotonically increasing.    ",
    " Linear interpolation and constant extrapolation of the specified velocities",
    " is used to compute the velocities at times not specified.\t\t     ",
    " The same holds for the anisotropy coefficients as a function of time only. ",
    "\t\t\t\t\t\t\t\t\t     ",
    " For NMO with a velocity function of time and CDP, specify the array\t     ",
    "\t   cdp=cdp1,cdp2,...\t\t\t\t\t\t     ",
    " and, for each CDP specified, specify the vnmo and tnmo arrays as described ",
    " above. The first (vnmo,tnmo) pair corresponds to the first cdp, and so on. ",
    " Linear interpolation and constant extrapolation of 1/velocity^2 is used    ",
    " to compute velocities at CDPs not specified.\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Anisotropy option:\t\t\t\t\t\t\t     ",
    " Caveat, this is an experimental option,\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The anisotropy coefficients anis1, anis2 permit non-hyperbolicity due\t     ",
    " to layering, mode conversion, or anisotropy. Default is isotropic NMO.     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The same holds for the anisotropy coefficients as a function of time and   ",
    " CDP.\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Moveout is defined by\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   1\t\t anis1\t\t\t\t\t\t\t     ",
    "  --- x^2 + ------------- x^4.\t\t\t\t\t\t     ",
    "  v^2\t     1 + anis2 x^2\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Note: In general, the user should set the cdp parameter.  The default is   ",
    "\tto use tr.cdp from the first trace and assume only one cdp.\t  ",
    " Caveat:\t\t\t\t\t\t\t\t     ",
    " Nmo cannot handle negative moveout as in triplication caused by\t     ",
    " anisotropy. But negative moveout happens necessarily for negative anis1 at ",
    " sufficiently large offsets. Then the error-negative moveout- is printed.   ",
    " Check anis1. An error (anis2 too small) is also printed if the\t     ",
    " denominator of the quartic term becomes negative. Check anis2. These errors",
    " are prompted even if they occur in traces which would not survive the\t     ",
    " NMO-stretch threshold. Chop off enough far-offset traces (e.g. with suwind)",
    " if anis1, anis2 are fine for near-offset traces.\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
];
fn interp(cdp: &[f32], a: &[Vec<Vec<f32>>], x: f32, o: &mut [Vec<f32>], idx: &mut i32) {
    let n = cdp.len();
    if x <= cdp[0] {
        for k in 0..3 {
            o[k].copy_from_slice(&a[k][0]);
        }
    } else if x >= cdp[n - 1] {
        for k in 0..3 {
            o[k].copy_from_slice(&a[k][n - 1]);
        }
    } else {
        xindex(cdp, x, idx);
        let j = *idx as usize;
        let p = (cdp[j + 1] - x) / (cdp[j + 1] - cdp[j]);
        let q = (x - cdp[j]) / (cdp[j + 1] - cdp[j]);
        for k in 0..3 {
            for i in 0..o[k].len() {
                o[k][i] = p * a[k][j][i] + q * a[k][j + 1][i];
            }
        }
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    // A zero-length trace cannot supply the samples required by moveout or
    // inverse interpolation. The original accesses invalid memory in this case.
    if nt == 0 {
        err!("trace has no samples");
    }
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut nc = countparval("cdp");
    if nc > 0 {
        if countparname("vnmo") != nc {
            err!("a vnmo array must be specified for each cdp");
        }
        if countparname("tnmo") != nc {
            err!("a tnmo array must be specified for each cdp");
        }
        if countparname("anis1") != nc && countparname("anis1") != 0 {
            err!("an anis1 array must be specified for each cdp, or omitted at all");
        }
        if countparname("anis2") != nc && countparname("anis2") != 0 {
            err!("an anis2 array must be specified for each cdp, or omitted at all");
        }
    } else {
        nc = 1;
        if countparname("vnmo") > 1 {
            err!("only one (or no) vnmo array must be specified");
        }
        if countparname("tnmo") > 1 {
            err!("only one (or no) tnmo array must be specified");
        }
        if countparname("anis1") > 1 {
            err!("only one (or no) anis1 array must be specified");
        }
        if countparname("anis2") > 1 {
            err!("only one (or no) anis2 array must be specified");
        }
    }
    let nc = nc as usize;
    let mut cdp = vec![0.0; nc];
    if getparfloat("cdp", &mut cdp) == 0 {
        cdp[0] = tr.cdp as f32;
    }
    let mut a = vec![vec![vec![0.0; nt]; nc]; 3];
    for j in 0..nc {
        let n = (j + 1) as i32;
        let mut nv = countnparval(n, "vnmo");
        let mut nn = countnparval(n, "tnmo");
        let mut n1 = countnparval(n, "anis1");
        let mut n2 = countnparval(n, "anis2");
        if nv != nn && !(nc == 1 && nv == 1 && nn == 0) {
            err!("number of vnmo and tnmo values must be equal");
        }
        if n1 != nv && n1 != 0 {
            err!("number of vnmo and anis1 values must be equal");
        }
        if n2 != nv && n2 != 0 {
            err!("number of vnmo and anis2 values must be equal");
        }
        if nv == 0 {
            nv = 1;
        }
        if nn == 0 {
            nn = nv;
        }
        if n1 == 0 {
            n1 = nv;
        }
        if n2 == 0 {
            n2 = nv;
        }
        let mut v = vec![0.0; nv as usize];
        let mut t = v.clone();
        let mut b = v.clone();
        let mut c = v.clone();
        if getnparfloat(n, "vnmo", &mut v) == 0 {
            v[0] = 1500.0;
        }
        if getnparfloat(n, "tnmo", &mut t) == 0 {
            t[0] = 0.0;
        }
        getnparfloat(n, "anis1", &mut b);
        getnparfloat(n, "anis2", &mut c);
        for i in 1..nn as usize {
            if t[i] <= t[i - 1] {
                err!("tnmo values must increase monotonically");
            }
        }
        for k in 0..3 {
            let vals = if k == 0 {
                &v
            } else if k == 1 {
                &b
            } else {
                &c
            };
            let end = if k == 0 {
                nv
            } else if k == 1 {
                n1
            } else {
                n2
            } as usize;
            let mut tn = ft;
            for it in 0..nt {
                let mut z = [0.0];
                intlin(
                    &t[..nn as usize],
                    vals,
                    vals[0],
                    vals[end - 1],
                    &[tn],
                    &mut z,
                );
                a[k][j][it] = if k == 0 {
                    (1.0f64 / (z[0] * z[0]) as f64) as f32
                } else {
                    z[0]
                };
                tn += dt;
            }
        }
    }
    for j in 1..nc {
        let mut i = j;
        while i > 0 && cdp[i - 1] > cdp[i] {
            cdp.swap(i - 1, i);
            for x in &mut a {
                x.swap(i - 1, i);
            }
            i -= 1;
        }
    }
    let mut smute = 0.0;
    if getparfloat("smute", &mut smute) == 0 {
        smute = 1.5;
    }
    if smute <= 0.0 {
        err!("smute must be greater than 0.0");
    }
    let mut lmute = 0;
    if getparint("lmute", &mut lmute) == 0 {
        lmute = 25;
    }
    let mut sscale = 0;
    if getparint("sscale", &mut sscale) == 0 {
        sscale = 1;
    }
    let mut invert = 0;
    getparint("invert", &mut invert);
    let mut upward = 0;
    getparint("upward", &mut upward);
    checkpars();
    let mut o = vec![vec![0.0; nt]; 3];
    let mut idx = 0;
    interp(&cdp, &a, tr.cdp as f32, &mut o, &mut idx);
    let mut tt = vec![0.0f32; nt.max(2)];
    let mut an = tt.clone();
    let mut ti = tt.clone();
    let mut at = tt.clone();
    let mut q = vec![0.0; nt];
    let mut oldcdp = tr.cdp;
    let mut oldoff = tr.offset.wrapping_sub(1);
    let mut mute = 0;
    loop {
        let changed = tr.cdp != oldcdp && nc > 1;
        if changed {
            interp(&cdp, &a, tr.cdp as f32, &mut o, &mut idx);
        }
        if changed || tr.offset != oldoff {
            let offset = tr.offset as f32;
            let temp = offset * offset / (dt * dt);
            let mut tn = ft / dt;
            for i in 0..nt {
                let tsq = ((temp * o[0][i]) as f64
                    + (o[1][i] * temp * temp) as f64 / (1.0 + (o[2][i] * temp) as f64))
                    as f32;
                if tsq < 0.0 {
                    err!("negative moveout; check anis1, anis2, or suwind far-offset traces");
                }
                if 1.0 + (o[2][i] * temp) as f64 <= 0.0 {
                    err!("anis2 negative and too small; check anis2, or suwind far-offset traces");
                }
                tt[i] = ((tn * tn + tsq) as f64).sqrt() as f32;
                tn = (tn as f64 + 1.0) as f32;
            }
            an[0] = tt[1] - tt[0];
            for i in 1..nt {
                an[i] = tt[i] - tt[i - 1];
            }
            let os = (1.0f64 / smute as f64) as f32;
            let mut i = if upward == 0 { 0 } else { nt - 1 };
            if upward == 0 {
                while i < nt - 1 && an[i] < os {
                    i += 1;
                }
            } else {
                while i > 0 && an[i] >= os {
                    i -= 1;
                }
            }
            mute = i;
            if invert != 0 {
                let start = ((ft / dt) as f64 + mute as f64) as f32;
                yxtoxy(
                    1.0,
                    start,
                    &tt[mute..nt],
                    1.0,
                    start,
                    ((ft / dt) as f64 - nt as f64) as f32,
                    ((ft / dt) as f64 + nt as f64) as f32,
                    &mut ti[mute..nt],
                );
                mute = ((1.0 + tt[mute] as f64 - (ft / dt) as f64) as i32)
                    .min(nt as i32 - 2)
                    .max(0) as usize;
                if sscale != 0 {
                    for i in mute + 1..nt {
                        at[i] = ti[i] - ti[i - 1];
                    }
                    at[mute] = at[mute + 1];
                }
            }
        }
        if invert == 0 {
            ints8r(
                1.0,
                ft / dt,
                &tr.data[..nt],
                0.0,
                0.0,
                &tt[mute..nt],
                &mut q[mute..nt],
            );
            q[..mute].fill(0.0);
            for i in mute..nt {
                if (i as i32) < mute as i32 + lmute {
                    q[i] *= (i - mute + 1) as f32 / lmute as f32;
                }
            }
            if sscale != 0 {
                for i in mute..nt {
                    q[i] *= an[i];
                }
            }
        } else {
            ints8r(
                1.0,
                ft / dt,
                &tr.data[..nt],
                0.0,
                0.0,
                &ti[mute..nt],
                &mut q[mute..nt],
            );
            q[..mute].fill(0.0);
            if sscale != 0 {
                for i in mute..nt {
                    q[i] *= at[i];
                }
            }
        }
        tr.data[..nt].copy_from_slice(&q);
        puttr(&mut tr);
        oldoff = tr.offset;
        oldcdp = tr.cdp;
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
