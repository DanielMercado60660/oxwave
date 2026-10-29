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

//! SUINTERPFOWLER - interpolate images from velocity panels.
//! Port of su/main/interp_extrap/suinterpfowler.c; John E. Anderson, Spring 1994.

use oxwave_core::prelude::*;
use oxwave_linear_interp::intlin;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUINTERPFOWLER - interpolate output image from constant velocity panels",
    "\t   built by SUTIFOWLER or CVS\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " These parameters should be specified the same as in SUTIFOWLER:\t",
    " vmin=1500.\t\tminimum velocity\t\t\t\t",
    " vmax=2500.\t\tmaximum velocity\t\t\t\t",
    " nv=21\t\t\tnumber of velocity panels\t\t\t",
    " etamin=0.10\t\tminimum eta value\t\t\t\t",
    " etamax=0.25\t\tmaximum eta value\t\t\t\t",
    " neta=11\t\tnumber of eta values\t\t\t\t",
    " ncdps=1130\t\tnumber of cdp points\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " If these parameters are specified so that nvstack>5, then the input \t",
    " data are assumed to come from CVS and the SUTIFOWLER parameters are ignored.",
    " nvstack=0\t\tnumber of constant velocity stack panels output by CVS",
    " vminstack=1450\tminimum velocity specified for CVS\t\t",
    " vscale=1.0\t\tscale factor for velocity functions\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " These parameters specify the desired output (time,velocity,eta) model\t",
    " at each cdp location. The sequential cdp numbers should be specified in",
    " increasing order from 0 to 'ncdps-1' at from 1 to 'ncdps' control point",
    " locations. (Time values are in seconds.)\t\t\t\t",
    " cdp=0\t\t\tcdp number for (t,v,eta) triplets (specify more than",
    " \t\t\t\tonce if needed)\t\t\t\t",
    " t=0.\t\t\tarray of times for (t,v,eta) triplets (specify more",
    "\t\t\t\tthan once if needed)\t\t\t",
    " v=1500.\t\tarray of velocities for (t,v,eta) triplets (specify",
    "\t\t\t\tmore than once if needed)\t\t",
    " eta=0.\t\tarray of etas for (t,v,eta) triplets (specify more",
    "\t\t\t\tthan once if needed)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: This is a simple research code based on linear interpolation.\t",
    " There are no protections against aliasing built into the code beyond\t",
    " suggesting that this program have a knowledgable user. A final version",
    " should do a better job taking care of endpoint conditions.",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut nv, mut neta, mut ncdps, mut nvstack) = (21, 11, 1130, 0);
    let (mut vmin, mut vmax, mut etamin, mut etamax, mut vminstack, mut vscale) =
        (1500.0f32, 2500.0f32, 0.10f32, 0.25f32, 0.25f32, 1.0f32);
    getparint("nv", &mut nv);
    getparint("neta", &mut neta);
    getparint("ncdps", &mut ncdps);
    getparfloat("vmin", &mut vmin);
    getparfloat("vmax", &mut vmax);
    getparfloat("etamin", &mut etamin);
    getparfloat("etamax", &mut etamax);
    getparint("nvstack", &mut nvstack);
    getparfloat("vminstack", &mut vminstack);
    getparfloat("vscale", &mut vscale);
    checkpars();
    let mut dv = (vmax - vmin) / (nv - 1).max(1) as f32;
    let mut deta = (etamax - etamin) / (neta - 1).max(1) as f32;
    if deta == 0.0 {
        deta = 1.0;
        neta = 1;
    }
    if nvstack > 5 {
        etamin = 0.0;
        deta = 1.0;
        vmin = 0.0;
        nv = nvstack;
        neta = 1;
        dv = (1.0f64 / ((vminstack * vminstack * (nvstack - 5) as f32) as f64)) as f32;
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (0.000001f64 * tr.dt as f64) as f32;
    fprintf!(
        stderr(),
        "supaint: sample rate = %f s, number of samples = %d\n",
        dt,
        nt as i32
    );
    let mut p = vec![0.0f32; (nv * neta) as usize * nt];
    let nc = countparname("cdp") as usize;
    if countparname("v") != nc as i32 {
        err!("v array must be specified for each cdp");
    }
    if countparname("t") != nc as i32 {
        err!("t array must be specified for each cdp");
    }
    if nvstack < 6 && countparname("eta") != nc as i32 {
        err!("eta array must be specified for each cdp");
    }
    let mut vel = vec![0.0f32; nc * nt];
    let mut eta = vec![0.0f32; nc * nt];
    let mut time = vec![0.0f32; nt];
    let mut c = vec![0.0f32; nc];
    for it in 0..nt {
        time[it] = (it as f32) * dt;
    }
    for ic in 0..nc {
        getnparfloat((ic + 1) as i32, "cdp", &mut c[ic]);
        let nvel = countnparval((ic + 1) as i32, "v") as usize;
        if nvstack > 5 {
            fprintf!(
                stderr(),
                "cdp=%f, number of (time,velocity) values = %d\n",
                c[ic],
                nvel as i32
            );
            if countnparval((ic + 1) as i32, "t") != nvel as i32 {
                err!("For each cdp control point, t and v must have the same number of values");
            }
        } else {
            fprintf!(
                stderr(),
                "cdp=%f, number of triplet values = %d\n",
                c[ic],
                nvel as i32
            );
            if countnparval((ic + 1) as i32, "t") != nvel as i32
                || countnparval((ic + 1) as i32, "eta") != nvel as i32
            {
                err!("For each cdp control point, t,v,and eta must have the same number of values");
            }
        }
        // The C code dereferences empty arrays in this case.
        if nvel == 0 {
            err!("control point must have at least one time and velocity value");
        }
        let mut v = vec![0.0f32; nvel];
        let mut t = vec![0.0f32; nvel];
        getnparfloat((ic + 1) as i32, "t", &mut t);
        getnparfloat((ic + 1) as i32, "v", &mut v);
        for x in &mut v {
            *x *= vscale;
        }
        if nvstack > 5 {
            for j in 0..nvel {
                fprintf!(stderr(), "\ttime=%f velocity=%f\n", t[j], v[j]);
            }
        } else {
            let mut e = vec![0.0f32; nvel];
            getnparfloat((ic + 1) as i32, "eta", &mut e);
            for j in 0..nvel {
                fprintf!(stderr(), "\ttime=%f velocity=%f eta=%f\n", t[j], v[j], e[j]);
            }
            intlin(
                &t,
                &e,
                e[0],
                e[nvel - 1],
                &time,
                &mut eta[ic * nt..(ic + 1) * nt],
            );
        }
        intlin(
            &t,
            &v,
            v[0],
            v[nvel - 1],
            &time,
            &mut vel[ic * nt..(ic + 1) * nt],
        );
    }
    if nvstack > 5 {
        for j in 0..nt * nc {
            vel[j] = (1.0f64 / (vel[j] * vel[j]) as f64) as f32;
            eta[j] = 0.0;
        }
    }
    for icmp in 0..ncdps {
        for j in 0..(nv * neta) as usize {
            if !(icmp == 0 && j == 0) && gettr(&mut tr) == 0 {
                err!("can't get input trace");
            }
            for it in 0..nt {
                p[it + j * nt] = tr.data[it];
            }
        }
        // The original dereferences c[0] when there are no control points.
        if nc == 0 {
            err!("at least one cdp control point must be specified");
        }
        let (mut ic1, mut ic2, mut wc) = (0usize, 0usize, 0.0f32);
        if icmp as f32 <= c[0] {
        } else if icmp as f32 >= c[nc - 1] {
            ic1 = nc - 1;
            ic2 = nc - 1;
        } else {
            for ic in 1..nc {
                ic1 = ic - 1;
                ic2 = ic;
                wc = (c[ic2] - icmp as f32) / (c[ic2] - c[ic1]);
                if icmp as f32 > c[ic1] && icmp as f32 <= c[ic2] {
                    break;
                }
            }
        }
        for it in 0..nt {
            let ee = ((wc * eta[it + ic1 * nt]) as f64
                + (1.0 - wc as f64) * eta[it + ic2 * nt] as f64
                - etamin as f64)
                / deta as f64;
            let vv = ((wc * vel[it + ic1 * nt]) as f64
                + (1.0 - wc as f64) * vel[it + ic2 * nt] as f64
                - vmin as f64)
                / dv as f64;
            let ee = ee as f32;
            let vv = vv as f32;
            let ie = ee as i32;
            let we = ee - ie as f32;
            let iv = vv as i32;
            let wv = vv - iv as f32;
            tr.data[it] = 0.0;
            if neta == 1 {
                if iv >= 0 && iv < nv - 1 {
                    tr.data[it] =
                        (1.0 - wv) * p[it + iv as usize * nt] + wv * p[it + (iv + 1) as usize * nt];
                }
            } else if iv >= 0 && iv < nv - 1 && ie >= 0 && ie < neta - 1 {
                tr.data[it] = (1.0 - we) * (1.0 - wv) * p[it + (ie * nv + iv) as usize * nt]
                    + (1.0 - we) * wv * p[it + (ie * nv + iv + 1) as usize * nt]
                    + we * (1.0 - wv) * p[it + ((ie + 1) * nv + iv) as usize * nt]
                    + we * wv * p[it + ((ie + 1) * nv + iv + 1) as usize * nt];
            }
        }
        tr.f2 = 0.0;
        tr.d2 = 1.0;
        tr.f1 = 0.0;
        tr.d1 = dt;
        tr.offset = 0;
        tr.cdp = icmp;
        tr.igc = 0;
        tr.igi = 0;
        tr.sx = (0.5f64 * tr.sx.wrapping_add(tr.gx) as f64) as i32;
        tr.gx = tr.sx;
        tr.cdpt = 0;
        puttr(&mut tr);
        if icmp == 20 * (icmp / 20) {
            fprintf!(stderr(), "supaint: have output cmp %d\n", icmp);
        }
    }
    exit(CWP_Exit());
}
