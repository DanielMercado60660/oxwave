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

//! DZDV: par/main/ray_theory/dzdv.c, revision 1.6; CWP: Zhenyue Liu, 1995.
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
/// Self-documentation generated from dzdv.c.
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " DZDV - determine depth derivative with respect to the velocity\t",
    "  parameter, dz/dv,  by ratios of migrated data with the primary \t",
    "  amplitude and those with the extra amplitude\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " dzdv <infile afile=afile dfile=dfile>outfile [parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " infile=\tinput of common image gathers with primary amplitude\t",
    " afile=\tinput of common image gathers with extra amplitude\t",
    " dfile=\toutput of imaged depths in common image gathers \t",
    " nx= \t        number of migrated traces \t\t\t\t",
    " nz=\t        number of points in migrated traces \t\t\t",
    " dx=\t\thorizontal spacing of migrated trace \t\t\t",
    " dz=\t        vertical spacing of output trace \t\t\t",
    " fx=\t        x-coordinate of first migrated trace \t\t\t",
    " fz=\t        z-coordinate of first point in migrated trace \t\t",
    " off0=         first offset in common image gathers \t\t\t",
    " noff=\t        number of offsets in common image gathers  \t\t",
    " doff=\t        offset increment in common image gathers  \t\t",
    " cip=x1,z1,r1,..., cip=xn,zn,rn         description of input CIGS\t",
    "\tx\tx-value of a common image point\t\t\t\t",
    "\tz\tz-value of a common image point\tat zero offset\t\t",
    "\tr\tr-parameter in a common image gather\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " outfile\toutput of dz/dv at the imaged points\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nxw, nzw=0\t\twindow widths along x- and z-directions in \t",
    "\t\t\twhich points are contributed in solving dz/dv. \t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This program is used as part of the velocity analysis technique developed",
    " by Zhenyue Liu, CWP:1995.\t\t\t\t\t\t",
];

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_linear_interp::intlin;

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut nx, mut nz, mut noff) = (0, 0, 0);
    let (mut dx, mut dz, mut fx, mut fz, mut off0, mut doff) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx!\n");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!\n");
    }
    if getparfloat("dx", &mut dx) == 0 {
        err!("must specify dx!\n");
    }
    if getparfloat("dz", &mut dz) == 0 {
        err!("must specify dz!\n");
    }
    getparfloat("fx", &mut fx);
    getparfloat("fz", &mut fz);
    if getparfloat("off0", &mut off0) == 0 {
        err!("must specify off0!\n");
    }
    if getparfloat("doff", &mut doff) == 0 {
        err!("must specify doff!\n");
    }
    if getparint("noff", &mut noff) == 0 {
        err!("must specify noff!\n");
    }
    let mut afile = String::from("afile");
    let mut dfile = String::from("dfile");
    getparstring("afile", &mut afile);
    getparstring("dfile", &mut dfile);
    let (mut nxw, mut nzw) = (0, 0);
    getparint("nxw", &mut nxw);
    if nxw > nx - 1 {
        nxw = nx - 1;
    }
    getparint("nzw", &mut nzw);
    if nzw > nz - 1 {
        nzw = nz - 1;
    }
    let ncip = countparname("cip");
    if ncip < 1 {
        err!("Number of CIGS must be greater 0!\n");
    }
    checkpars();
    let mut xcip = vec![0.0; ncip as usize];
    let mut zcip = xcip.clone();
    let mut rcip = xcip.clone();
    let count = (nx * nz * noff) as usize;
    let mut g0 = vec![0.0f32; count];
    let mut g = g0.clone();
    let mut x = vec![0.0; nx as usize];
    let mut z0 = x.clone();
    let mut r = x.clone();
    let mut rat = vec![0.0f32; noff as usize];
    let mut z = rat.clone();
    if fread(&mut g0, &stdin()) != count {
        err!("cannot read %d values from file %s\n", count, "");
    }
    let afp = fopen(&afile, "r").unwrap_or_else(|| crash());
    if fread(&mut g, &afp) != count {
        err!("cannot read %d values from file %s\n", count, "");
    }
    let mut temp = vec![0.0f32; 3];
    for icdp in 0..ncip as usize {
        getnparfloat(icdp as i32 + 1, "cip", &mut temp);
        xcip[icdp] = temp[0];
        zcip[icdp] = temp[1];
        rcip[icdp] = temp[2];
    }
    let dfp = fopen(&dfile, "w").unwrap_or_else(|| crash());
    for ix in 0..nx as usize {
        x[ix] = fx + ix as f32 * dx;
    }
    intlin(&xcip, &zcip, zcip[0], zcip[ncip as usize - 1], &x, &mut z0);
    intlin(&xcip, &rcip, rcip[0], rcip[ncip as usize - 1], &x, &mut r);
    for ix in 0..nx {
        let ixl = if ix < nxw { 0 } else { ix - nxw };
        let ixh = if ix > nx - 1 - nxw { nx - 1 } else { ix + nxw };
        let z2 = z0[ix as usize] * z0[ix as usize];
        for ioff in 0..noff {
            let offs = off0 + ioff as f32 * doff;
            let mut temp1 = z2 + r[ix as usize] * offs * offs;
            if temp1 < 0.0 {
                temp1 = 0.0;
            }
            z[ioff as usize] = (temp1 as f64).sqrt() as f32;
            let iz = nint(((z[ioff as usize] - fz) / dz) as f64);
            let mut izl = if iz < nzw { 0 } else { iz - nzw };
            if izl > nz - 1 {
                izl = nz - 1;
            }
            let mut izh = if iz > nz - 1 - nzw { nz - 1 } else { iz + nzw };
            if izh < 0 {
                izh = 0;
            }
            let (mut fden, mut fnum) = (0.0f32, 0.0f32);
            for jx in ixl..=ixh {
                for iz in izl..=izh {
                    let i = ((jx * noff + ioff) * nz + iz) as usize;
                    fden += g0[i] * g0[i];
                    fnum += g[i] * g0[i];
                }
            }
            rat[ioff as usize] = if fden == 0.0 { 0.0 } else { fnum / fden };
        }
        fwrite(&rat, &stdout());
        fwrite(&z, &dfp);
    }
    exit(CWP_Exit());
}

// The original passes a null FILE pointer to libc when fopen fails.
fn crash() -> ! {
    unsafe extern "C" {
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        raise(11);
    }
    std::process::abort()
}
