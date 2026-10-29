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

//! MAKEVEL - make a velocity function. Port of par/main/velocity_model_building/makevel.c (Revision 1.12).
//! Author: Dave Hale.

use oxwave_core::prelude::*;
use oxwave_f_frannor::frannor;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " MAKEVEL - MAKE a VELocity function v(x,y,z)\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " makevel > outfile nx= nz= [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " nx=                    number of x samples (3rd dimension)\t\t",
    " nz=                    number of z samples (1st dimension)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " ny=1                   number of y samples (2nd dimension)\t\t",
    " dx=1.0                 x sampling interval\t\t\t\t",
    " fx=0.0                 first x sample\t\t\t\t\t",
    " dy=1.0                 y sampling interval\t\t\t\t",
    " fy=0.0                 first y sample\t\t\t\t\t",
    " dz=1.0                 z sampling interval\t\t\t\t",
    " fz=0.0                 first z sample\t\t\t\t\t",
    " v000=2.0               velocity at (x=0,y=0,z=0)\t\t\t",
    " dvdx=0.0               velocity gradient with respect to x\t\t",
    " dvdy=0.0               velocity gradient with respect to y\t\t",
    " dvdz=0.0               velocity gradient with respect to z\t\t",
    " vlens=0.0              velocity perturbation in parabolic lens\t",
    " tlens=0.0              thickness of parabolic lens\t\t\t",
    " dlens=0.0              diameter of parabolic lens\t\t\t",
    " xlens=                 x coordinate of center of parabolic lens\t",
    " ylens=                 y coordinate of center of parabolic lens\t",
    " zlens=                 z coordinate of center of parabolic lens\t",
    " lambda=1.0             make lambda larger to sharpen edge of lens \t",
    " vran=0.0\t\tstandard deviation of random perturbation\t",
    " vzfile=                file containing v(z) profile\t\t\t",
    " vzran=0.0              standard deviation of random perturbation to v(z)",
    " vzc=0.0                v(z) chirp amplitude\t\t\t\t",
    " z1c=fz                 z at which to begin chirp\t\t\t",
    " z2c=fz+(nz-1)*dz       z at which to end chirp\t\t\t",
    " l1c=dz                 wavelength at beginning of chirp\t\t",
    " l2c=dz                 wavelength at end of chirp\t\t\t",
    " exc=1.0                exponent of chirp\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut nx = 0; let mut nz = 0;
    if getparint("nx", &mut nx) == 0 { err!("must specify nx!\n"); }
    if getparint("nz", &mut nz) == 0 { err!("must specify nz!\n"); }
    let mut ny = 1; getparint("ny", &mut ny);
    let mut dx = 1.0f32; let mut dy = 1.0f32; let mut dz = 1.0f32;
    let mut fx = 0.0f32; let mut fy = 0.0f32; let mut fz = 0.0f32;
    getparfloat("dx", &mut dx); getparfloat("dy", &mut dy); getparfloat("dz", &mut dz);
    getparfloat("fx", &mut fx); getparfloat("fy", &mut fy); getparfloat("fz", &mut fz);
    let mut v000 = 2.0f32; let mut dvdx = 0.0f32; let mut dvdy = 0.0f32; let mut dvdz = 0.0f32;
    getparfloat("v000", &mut v000); getparfloat("dvdx", &mut dvdx);
    getparfloat("dvdy", &mut dvdy); getparfloat("dvdz", &mut dvdz);
    let mut xlens = fx; let mut ylens = fy; let mut zlens = fz;
    getparfloat("xlens", &mut xlens); getparfloat("ylens", &mut ylens); getparfloat("zlens", &mut zlens);
    let mut vlens = 0.0f32; let mut dlens = 1.0f32; let mut tlens = 1.0f32; let mut lambda = 1.0f32;
    getparfloat("vlens", &mut vlens); getparfloat("dlens", &mut dlens);
    getparfloat("tlens", &mut tlens); getparfloat("lambda", &mut lambda);
    let mut vran = 0.0f32; getparfloat("vran", &mut vran);
    let mut vzfile = String::new(); getparstring("vzfile", &mut vzfile);
    let mut vzran = 0.0f32; let mut vzc = 0.0f32;
    getparfloat("vzran", &mut vzran); getparfloat("vzc", &mut vzc);
    let mut z1c = fz; let mut z2c = fz + (nz.wrapping_sub(1) as f32) * dz;
    let mut l1c = dz; let mut l2c = dz; let mut exc = 1.0f32;
    getparfloat("z1c", &mut z1c); getparfloat("z2c", &mut z2c);
    getparfloat("l1c", &mut l1c); getparfloat("l2c", &mut l2c); getparfloat("exc", &mut exc);
    checkpars();

    let bc = (std::f64::consts::PI / ((z2c-z1c) as f64) * (1.0 / l2c as f64 - 1.0 / l1c as f64)) as f32;
    let ac = (2.0 * std::f64::consts::PI / l1c as f64 - 2.0 * bc as f64 * z1c as f64) as f32;
    // Negative sizes make alloc1float fail in C; the empty loops do not use it.
    let nalloc = if nz > 0 { nz as usize } else { 0 };
    let mut v = vec![0.0f32; nalloc];
    let mut vz = vec![0.0f32; nalloc];
    if !vzfile.is_empty() {
        let fp = match fopen(&vzfile, "r") { Some(fp) => fp, None => err!("error opening vzfile=%s", vzfile) };
        if fread(&mut vz, &fp) != nz as usize { err!("error reading vzfile=%s", vzfile); }
        fclose(&fp);
    }
    for iz in 0..nz { vz[iz as usize] += vzran * frannor(); }
    let mut x = fx;
    for _ix in 0..nx {
        let mut y = fy;
        for _iy in 0..ny {
            let mut z = fz;
            for iz in 0..nz {
                let i = iz as usize;
                v[i] = vz[i];
                v[i] += v000 + x*dvdx + y*dvdy + z*dvdz;
                let xn = (2.0 * (x-xlens) as f64 / dlens as f64) as f32;
                let ynn = (2.0 * (y-ylens) as f64 / dlens as f64) as f32;
                let zn = (2.0 * (z-zlens) as f64 / tlens as f64) as f32;
                v[i] = (v[i] as f64 + vlens as f64 * ((-lambda * (xn*xn + ynn*ynn + zn*zn)) as f64).exp()) as f32;
                if z > z1c && z < z2c {
                    let vtemp = (((ac + bc*z)*z) as f64).sin() as f32;
                    if vtemp < 0.0 {
                        v[i] = (v[i] as f64 - vzc as f64 * (-(vtemp as f64)).powf(exc as f64)) as f32;
                    } else {
                        v[i] = (v[i] as f64 + vzc as f64 * (vtemp as f64).powf(exc as f64)) as f32;
                    }
                }
                v[i] += vran * frannor();
                z += dz;
            }
            fwrite(&v, &stdout());
            y += dy;
        }
        x += dx;
    }
    exit(CWP_Exit());
}
