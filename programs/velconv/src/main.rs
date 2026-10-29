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

//! VELCONV - velocity conversion. Port of par/main/velocity_model_building/velconv.c (Revision 1.17).
//! Author: Dave Hale, Colorado School of Mines, 07/07/89.

use oxwave_core::prelude::*;
use oxwave_inverse_interp::yxtoxy;
use oxwave_linear_interp::intlin;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " VELCONV - VELocity CONVersion\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " velconv <infile >outfile intype= outtype= [optional parameters]",
    " \t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t",
    " intype=                input data type (see valid types below)",
    " outtype=               output data type (see valid types below)",
    " \t\t\t\t\t\t\t\t",
    " Valid types for input and output data are:\t\t\t",
    " vintt          interval velocity as a function of time\t",
    " vrmst          RMS velocity as a function of time\t\t",
    " vintz          velocity as a function of depth\t\t",
    " zt             depth as a function of time\t\t\t",
    " tz             time as a function of depth\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " These are active only if otherfile= is specified \t\t",
    " othert \tother as a function of time\t\t\t",
    " otherz \tother as a function of depth\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " nt=all                 number of time samples\t\t\t",
    " dt=1.0                 time sampling interval\t\t\t",
    " ft=0.0                 first time\t\t\t\t",
    " nz=all                 number of depth samples\t\t",
    " dz=1.0                 depth sampling interval\t\t",
    " fz=0.0                 first depth\t\t\t\t",
    " nx=all                 number of traces\t\t\t",
    " verbose=0\t\t silent, =1 chatty \t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Other input: \t\t\t\t\t\t\t",
    " otherfile=\tinput file containg other input\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " If type=othert, then vint(t) is assumed. If type=otherz then     ",
    "  vint(z) is assumed. \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Example:  \"intype=vintz outtype=vrmst\" converts an interval velocity",
    "           function of depth to an RMS velocity function of time.",
    " \t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t",
    " The variables t, dt, and ft are used only for input and output ",
    " functions of time; you need specify these only for vintt, vrmst,",
    " orzt.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Likewise, nz, dz, and fz are used only for input and output\t",
    " functions of depth.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " The otherfile= option is an input file with the same number of",
    " samples as the input interval velocity as a function of depth",
    " or time. The otherfile may contain impedence, density, or other",
    " model that is desired to be converted from depth to time or time to depth.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " The input and output data formats are C-style binary floats.\t",
];

fn zttz(dt: f32, ft: f32, zt: &[f32], vft: f32, vlt: f32, dz: f32, fz: f32, tz: &mut [f32]) {
    let nt = zt.len();
    let nz = tz.len();
    let lt = ft + (nt as i32 - 1) as f32 * dt;
    let lz = fz + (nz as i32 - 1) as f32 * dz;
    yxtoxy(dt, ft, zt, dz, fz, 0.0, 0.0, tz);
    let (mut i, mut z) = (0, fz);
    while i < nz && z <= zt[0] {
        tz[i] = (2.0 * z as f64 / vft as f64) as f32;
        i += 1;
        z += dz;
    }
    let (mut i, mut z) = (nz as isize - 1, lz);
    while i >= 0 && z >= zt[nt - 1] {
        tz[i as usize] = (lt as f64 + 2.0 * (z - zt[nt - 1]) as f64 / vlt as f64) as f32;
        i -= 1;
        z -= dz;
    }
}
fn tzzt(dz: f32, fz: f32, tz: &[f32], vfz: f32, vlz: f32, dt: f32, ft: f32, zt: &mut [f32]) {
    let nt = zt.len();
    let nz = tz.len();
    let lt = ft + (nt as i32 - 1) as f32 * dt;
    let lz = fz + (nz as i32 - 1) as f32 * dz;
    yxtoxy(dz, fz, tz, dt, ft, 0.0, 0.0, zt);
    let (mut i, mut t) = (0, ft);
    while i < nt && t <= tz[0] {
        zt[i] = (0.5 * t as f64 * vfz as f64) as f32;
        i += 1;
        t += dt;
    }
    let (mut i, mut t) = (nt as isize - 1, lt);
    while i >= 0 && t >= tz[nz - 1] {
        zt[i as usize] = (lz as f64 + 0.5 * (t - tz[nz - 1]) as f64 * vlz as f64) as f32;
        i -= 1;
        t -= dt;
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(2, SDOC);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let (mut intype, mut outtype) = (String::new(), String::new());
    if getparstring("intype", &mut intype) == 0 {
        err!("Must specify intype!\n");
    }
    if getparstring("outtype", &mut outtype) == 0 {
        err!("Must specify outtype!\n");
    }
    let infp = stdin();
    let outfp = stdout();
    let (mut nt, mut nz) = (0, 0);
    let (mut dt, mut ft, mut dz, mut fz) = (1.0f32, 0.0f32, 1.0f32, 0.0f32);
    if getparint("nt", &mut nt) == 0 {
        if getparint("nz", &mut nz) != 0 {
            nt = nz;
        } else {
            efseek(&infp, 0, 2);
            nt = (eftell(&infp) / 4) as i32;
        }
    }
    getparfloat("dt", &mut dt);
    getparfloat("ft", &mut ft);
    if getparint("nz", &mut nz) == 0 {
        if getparint("nt", &mut nt) != 0 {
            nz = nt;
        } else {
            efseek(&infp, 0, 2);
            nz = (eftell(&infp) / 4) as i32;
        }
    }
    getparfloat("dz", &mut dz);
    getparfloat("fz", &mut fz);
    let mut otherfile = String::new();
    getparstring("otherfile", &mut otherfile);
    let other = !otherfile.is_empty();
    if other && intype != "otherz" && intype != "othert" {
        err!("Must have intype=otherz or othert if otherfile= set!");
    }
    if other && outtype != "otherz" && outtype != "othert" {
        err!("Must have outtype=othert or otherz if otherfile= set!");
    }
    let otherfp = if other {
        Some(efopen(&otherfile, "r"))
    } else {
        None
    };
    let nin = match intype.as_str() {
        "vintt" | "vrmst" | "zt" => nt,
        "vintz" | "otherz" | "tz" => nz,
        "othert" if other => {
            if verbose != 0 {
                warn!("assuming vintt input velocity file!");
            }
            nt
        }
        _ => err!("invalid intype=%s!\n", intype),
    } as usize;
    let nout = match outtype.as_str() {
        "othert" | "vintt" | "vrmst" | "zt" => nt,
        "otherz" | "vintz" | "tz" => nz,
        _ => err!("invalid outtype=%s!\n", outtype),
    } as usize;
    let mut nx = -1;
    getparint("nx", &mut nx);
    checkpars();
    let nt = nt as usize;
    let nz = nz as usize;
    let mut tz = vec![0.0f32; nz];
    let mut zt = vec![0.0f32; nt];
    let mut din = vec![0.0f32; nin];
    let mut dout = vec![0.0f32; nout];
    let mut otherin = vec![0.0f32; nin];
    efseek(&infp, 0, 0);
    if let Some(ref fp) = otherfp {
        efseek(fp, 0, 0);
    }
    let mut ix = 0;
    while ix < nx || nx < 0 {
        dout.fill(0.0);
        if efread(&mut din, &infp) != nin {
            break;
        }
        if let Some(ref fp) = otherfp {
            if efread(&mut otherin, fp) != nin {
                break;
            }
        }
        match intype.as_str() {
            "vintt" | "othert" => {
                if zt[0] < 0.0 {
                    err!("Negative value detected! Only positive velocities allowed");
                }
                zt[0] = (0.5 * ft as f64 * din[0] as f64) as f32;
                for i in 1..nt {
                    if zt[i] < 0.0 {
                        err!("Negative value detected! Only positive velocities allowed");
                    }
                    zt[i] = (zt[i - 1] as f64 + 0.5 * dt as f64 * din[i - 1] as f64) as f32;
                }
                zttz(dt, ft, &zt, din[0], din[nt - 1], dz, fz, &mut tz);
            }
            "vrmst" => {
                if din[0] < 0.0 {
                    err!("Negative value detected! Only positive velocities allowed");
                }
                let mut v = din[0];
                let tiny = (0.00001 * v as f64 * v as f64) as f32;
                zt[0] = (0.5 * ft as f64 * v as f64) as f32;
                let mut t = ft + dt;
                for i in 1..nt {
                    if din[i] < 0.0 {
                        err!("Negative value detected! Only positive velocities allowed");
                    }
                    let s = (t * din[i] * din[i] - (t - dt) * din[i - 1] * din[i - 1]) / dt;
                    v = (if s > tiny { s } else { tiny } as f64).sqrt() as f32;
                    zt[i] = (zt[i - 1] as f64 + 0.5 * dt as f64 * v as f64) as f32;
                    t += dt;
                }
                zttz(dt, ft, &zt, din[0], v, dz, fz, &mut tz);
            }
            "zt" => {
                zt.copy_from_slice(&din);
                let a = (2.0 * (zt[1] - zt[0]) as f64 / dt as f64) as f32;
                let b = (2.0 * (zt[nt - 1] - zt[nt - 2]) as f64 / dt as f64) as f32;
                zttz(dt, ft, &zt, a, b, dz, fz, &mut tz);
            }
            "vintz" | "otherz" => {
                if din[0] < 0.0 {
                    err!("Negative value detected! Only positive velocities allowed");
                }
                tz[0] = (2.0 * fz as f64 / din[0] as f64) as f32;
                for i in 1..nz {
                    if din[i] < 0.0 {
                        err!("Negative value detected! Only positive velocities allowed");
                    }
                    tz[i] = (tz[i - 1] as f64 + 2.0 * dz as f64 / din[i - 1] as f64) as f32;
                }
                tzzt(dz, fz, &tz, din[0], din[nz - 1], dt, ft, &mut zt);
            }
            "tz" => {
                tz.copy_from_slice(&din);
                let a = (2.0 * dz as f64 / (tz[1] - tz[0]) as f64) as f32;
                let b = (2.0 * dz as f64 / (tz[nz - 1] - tz[nz - 2]) as f64) as f32;
                tzzt(dz, fz, &tz, a, b, dt, ft, &mut zt);
            }
            _ => unreachable!(),
        }
        match outtype.as_str() {
            "vintt" => {
                for i in 0..nt - 1 {
                    dout[i] = (2.0 * (zt[i + 1] - zt[i]) as f64 / dt as f64) as f32;
                }
                dout[nt - 1] = dout[nt - 2];
            }
            "vrmst" => {
                let mut v = (2.0 * (zt[1] - zt[0]) as f64 / dt as f64) as f32;
                let mut sum = ft * v * v;
                dout[0] = v;
                let mut t = ft + dt;
                for i in 1..nt {
                    v = (2.0 * (zt[i] - zt[i - 1]) as f64 / dt as f64) as f32;
                    sum += dt * v * v;
                    dout[i] = ((sum / t) as f64).sqrt() as f32;
                    t += dt;
                }
            }
            "zt" => dout.copy_from_slice(&zt),
            "vintz" => {
                for i in 0..nz - 1 {
                    dout[i] = (2.0 * dz as f64 / (tz[i + 1] - tz[i]) as f64) as f32;
                }
                dout[nz - 1] = dout[nz - 2];
            }
            "tz" => dout.copy_from_slice(&tz),
            "othert" | "otherz" if other => {
                let x = if outtype == "othert" { &tz } else { &zt };
                intlin(
                    &x[..nin],
                    &otherin,
                    otherin[0],
                    otherin[nin - 1],
                    &x[..nout],
                    &mut dout,
                );
            }
            _ => err!("invalid outtype=%s!\n", outtype),
        }
        efwrite(&dout, &outfp);
        ix += 1;
    }
    exit(CWP_Exit());
}
