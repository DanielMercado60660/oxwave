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

//! WKBJ - Compute WKBJ ray theoretic parameters, via finite differencing.
//! Port of par/main/ray_theory/wkbj.c (Revision 1.10).
//! Credits (from the C source): CWP: Zhenyue Liu, Dave Hale, pre 1992.

use oxwave_core::prelude::*;
use oxwave_f_upweik::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " WKBJ - Compute WKBJ ray theoretic parameters, via finite differencing\t",
    "\t\t\t\t\t\t\t\t\t",
    " wkbj <vfile >tfile nx= nz= xs= zs= [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " <vfile\tfile containing velocities v[nx][nz]\t\t\t",
    " nx=\t\tnumber of x samples (2nd dimension)\t\t\t",
    " nz=\t\tnumber of z samples (1st dimension)\t\t\t",
    " xs=\t\tx coordinate of source\t\t\t\t\t",
    " zs=\t\tz coordinate of source\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dx=1.0\t\t x sampling interval\t\t\t\t",
    " fx=0.0\t\t first x sample\t\t\t\t\t",
    " dz=1.0\t\t z sampling interval\t\t\t\t",
    " fz=0.0\t\t first z sample\t\t\t\t\t",
    " sfile=sfile\tfile containing sigmas sg[nx][nz]\t\t\t",
    " bfile=bfile\tfile containing incident angles bet[nx][nz]\t\t",
    " afile=afile\tfile containing propagation angles a[nx][nz]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Traveltimes, propagation angles, sigmas, and incident angles in WKBJ\t",
    " by finite differences  in polar coordinates. Traveltimes are calculated",
    " by upwind scheme; sigmas and incident angles by a Crank-Nicolson scheme.",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let (mut nx, mut nz) = (0i32, 0i32);
    let (mut fx, mut fz, mut dx, mut dz, mut xs, mut zs) =
        (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let vfp = stdin();
    let tfp = stdout();
    let (mut bfile, mut sfile, mut afile) = (String::new(), String::new(), String::new());

    /* hook up getpar to handle the parameters */
    initargs(&argv);
    requestdoc(0, SDOC);

    /* get required parameters */
    if getparint("nx", &mut nx) == 0 {
        err!("must specify nx!\n");
    }
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz!\n");
    }
    if getparfloat("xs", &mut xs) == 0 {
        err!("must specify xs!\n");
    }
    if getparfloat("zs", &mut zs) == 0 {
        err!("must specify zs!\n");
    }

    /* get optional parameters */
    if getparfloat("dx", &mut dx) == 0 {
        dx = 1.0;
    }
    if getparfloat("fx", &mut fx) == 0 {
        fx = 0.0;
    }
    if getparfloat("dz", &mut dz) == 0 {
        dz = 1.0;
    }
    if getparfloat("fz", &mut fz) == 0 {
        fz = 0.0;
    }
    if getparstring("sfile", &mut sfile) == 0 {
        sfile = "sfile".to_string();
    }
    if getparstring("bfile", &mut bfile) == 0 {
        bfile = "bfile".to_string();
    }
    if getparstring("afile", &mut afile) == 0 {
        afile = "afile".to_string();
    }
    checkpars();

    let sfp = fopen(&sfile, "w").unwrap_or_else(|| err!("cannot open sfile=%s", sfile));
    let bfp = fopen(&bfile, "w").unwrap_or_else(|| err!("cannot open bfile=%s", bfile));
    let afp = fopen(&afile, "w").unwrap_or_else(|| err!("cannot open afile=%s", afile));

    /* ensure source is in grid */
    let ex = fx + nx.wrapping_sub(1) as f32 * dx;
    let ez = fz + nz.wrapping_sub(1) as f32 * dz;
    if fx > xs || ex < xs || fz > zs || ez < zs {
        err!("source lies outside of specified (x,z) grid\n");
    }

    /* allocate space: the C grids are contiguous, indexed [x][z]. */
    let n = nx as usize * nz as usize;
    let mut v = vec![0.0f32; n];
    let mut t = vec![0.0f32; n];
    let mut sg = vec![0.0f32; n];
    let mut a = vec![0.0f32; n];
    let mut bet = vec![0.0f32; n];

    /* read velocities; like C, leave the fread count unchecked. */
    fread(&mut v, &vfp);

    /* compute times, angles, sigma, and betas */
    eiktam(
        xs,
        zs,
        dz,
        fz,
        dx,
        fx,
        &v.chunks(nz as usize).collect::<Vec<_>>(),
        &mut t.chunks_mut(nz as usize).collect::<Vec<_>>(),
        &mut a.chunks_mut(nz as usize).collect::<Vec<_>>(),
        &mut sg.chunks_mut(nz as usize).collect::<Vec<_>>(),
        &mut bet.chunks_mut(nz as usize).collect::<Vec<_>>(),
    );

    /* write first-arrival times */
    fwrite(&t, &tfp);
    /* write sigma */
    fwrite(&sg, &sfp);
    /* write angle */
    fwrite(&a, &afp);
    /* write beta */
    fwrite(&bet, &bfp);

    /* close files */
    fclose(&sfp);
    fclose(&afp);
    fclose(&bfp);

    exit(CWP_Exit());
}
