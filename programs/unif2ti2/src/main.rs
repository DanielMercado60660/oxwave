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
// Copyright (c) Colorado School of Mines, 2010. All rights reserved.
//! UNIF2TI2, port of par/main/velocity_model_building/unif2ti2.c.
//! Modified by Pengfei Cai (CWP), Dec 2011, and Xiaoxiang Wang (CWP), Aug 2010.
//! Based on unif2aniso by John Stockwell, 2005.
use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_f_cubicspline::{cakima, cmonot, csplin};
use oxwave_f_refaniso::{thom2stiffTI, Stiff2D};
use oxwave_linear_interp::intlin;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " UNIF2TI2 - generate a 2-D UNIFormly sampled profile of stiffness \t",
    " \tcoefficients of a layered medium (only for P waves). \t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  unif2ti2 < infile [Parameters]\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " none \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " ninf=5\tnumber of interfaces\t\t\t\t\t",
    " nx=100\tnumber of x samples (2nd dimension)\t\t\t",
    " nz=100\tnumber of z samples (1st dimension)\t\t\t",
    " dx=10\t\tx sampling interval\t\t\t\t\t",
    " dz=10\t\tz sampling interval\t\t\t\t\t",
    " ns=0          number of samples in vertical direction for smoothing   ",
    "\t\tparameters across boundary\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " npmax=201\tmaximum number of points on interfaces\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " fx=0.0\tfirst x sample\t\t\t\t\t\t",
    " fz=0.0\tfirst z sample\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " x0=0.0,0.0,..., \tdistance x at which vp00 and vs00 are specified\t",
    " z0=0.0,0.0,..., \tdepth z at which vp00 and vs00 are specified\t",
    " \t\t\t\t\t\t\t\t\t",
    " vp00=1500,2000,...,\tP-velocity at each x0,z0 (m/sec)\t\t",
    " eps00=0,0,0...,\tThomsen parameter epsilon at each x0,z0\t\t",
    " delta00=0,0,0...,\tThomsen\tparameter delta at each x0,z0\t\t",
    " rho00=1000,1100,...,\tdensity at each x0,z0 (kg/m^3)\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dvpdx=0.0,0.0,...,\tx-derivative of P-velocity (dvp/dx)\t\t",
    " dvpdz=0.0,0.0,...,\tz-derivative of P-velocity (dvs/dz)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dedx=0.0,0.0,...,\tx-derivative of epsilon (de/dx)\t\t\t",
    " dedz=0.0,0.0,...,\tz-derivative of epsilon with depth z (de/dz)\t",
    " \t\t\t\t\t\t\t\t\t",
    " dddx=0.0,0.0,...,\tx-derivative of delta (dd/dx)\t\t\t",
    " dddz=0.0,0.0,...,\tz-derivative of delta (dd/dz)\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " drdx=0.0,0.0,...,\tx-derivative of density (d rho/dx)\t\t",
    " drdz=0.0,0.0,...,\tz-derivative of density (d rho/dz)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " nufile= \t\tbinary file containning tilt value at each grid point",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " ...output filenames \t\t\t\t\t\t\t",
    " c11_file=c11_file\toutput filename for c11 values\t \t\t",
    " c13_file=c13_file\toutput filename for c13 values\t \t\t",
    " c15_file=c15_file\toutput filename for c15 values\t \t\t",
    " c33_file=c33_file\toutput filename for c33 values\t \t\t",
    " c35_file=c35_file\toutput filename for c35 values\t \t\t",
    " c55_file=c55_file\toutput filename for c55 values\t \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " method=linear\t\tfor linear interpolation of interface\t\t",
    " \t\t\t=mono for monotonic cubic interpolation of interface",
    "\t\t\t=akima for Akima's cubic interpolation of interface",
    "\t\t\t=spline for cubic spline interpolation of interface",
    " \t\t\t\t\t\t\t\t\t",
    " tfile=\t\t=testfilename  if set, a sample input dataset is",
    " \t\t\t output to \"testfilename\".\t\t\t",
    " \t\t\t \t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The input file is an ASCII file containing x z values representing a\t",
    " piecewise continuous velocity model with a flat surface on top.\t",
    " \t\t\t \t\t\t\t\t\t",
    " The surface and each successive boundary between media is represented ",
    " by a list of selected x z pairs written column form. The first and\t",
    " last x values must be the same for all boundaries. Use the entry\t",
    " 1.0  -99999  to separate the entries for successive boundaries. No\t",
    " boundary may cross another. Note that the choice of the method of\t",
    " interpolation may cause boundaries to cross that do not appear to\t",
    " cross in the input data file.\t\t\t\t\t\t",
    " \t\t\t \t\t\t\t\t\t",
    " The number of interfaces is specified by the parameter \"ninf\". This ",
    " number does not include the top surface of the model. The input data\t",
    " format is the same as a CSHOT model file with all comments removed.\t",
    " \t\t\t \t\t\t\t\t\t",
    " The algorithm works by transforming the P-wavespeed, Thomsen parameters ",
    " epsilon and delta, and the tilt of the symmetry axis into density-normalized",
    " stiffness coefficients. \t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " At this stage, the tilt-field file can be prepared using the \t\t",
    " Matlab M-file nu_mod.m based on 2D interpolation between interfaces.  ",
    " The binary file contains nu values at each grid point.\t\t",
    " The interfaces are obtained by interpolation on the picked ones \t",
    " stored in the infile, and the symmetry axis at each point of interface",
    " is assumed to be parallel to the normal direction.\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Common ranges of Thomsen parameters are\t\t\t\t",
    "  epsilon:  0.0 -> 0.5\t\t\t\t\t\t\t",
    "  delta:   -0.2 -> 0.4\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " If the tilt-field file is not given, the model will be assumed to \t",
    " have VTI symmetry. \t\t\t\t \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example using test input file generating feature:\t\t\t",
    " unif2aniso tfile=testfilename  produces a 5 interface demonstration model",
    " unif2aniso < testfilename \t\t\t\t\t\t",
    " ximage < c11_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c13_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c15_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c33_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c35_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c55_file n1=100 n2=100\t\t\t\t\t",
    " ximage < rho_file n1=100 n2=100\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];
const SAMPLE: &[&str] = &[
    "   0         0",
    "1000         0",
    "   1    -99999",
    "   0       200",
    "1000       200",
    "   1    -99999",
    "   0       400",
    " 100.      400",
    " 400       560",
    " 550.      560",
    "1000.      360",
    "   1    -99999",
    "   0       800",
    " 100.      800",
    " 400.      750.",
    " 500.      780",
    " 550.      780",
    "1000.      500",
    "   1.   -99999",
    "   0       950",
    " 100       950",
    " 400       900.",
    " 500       800.",
    " 550       800.",
    " 750       920",
    "1000       930",
    "   1    -99999",
    "   0      1000",
    "1000      1000",
    "   1    -99999",
];
fn layer(name: &str, n: usize, base: f32, step: f32) -> Vec<f32> {
    let mut v = vec![0.; n];
    if getparfloat(name, &mut v) == 0 {
        for i in 0..n {
            v[i] = (base as f64 + step as f64 * i as f64) as f32;
        }
    }
    v
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tfile = String::new();
    getparstring("tfile", &mut tfile);
    if !tfile.is_empty() {
        let fp = efopen(&tfile, "w");
        for s in SAMPLE {
            fprintf!(&fp, "%s\n", *s);
        }
        efclose(&fp);
        exit(CWP_Exit());
    }
    let (mut ninf, mut nx, mut nz, mut npmax, mut ns) = (5, 100, 100, 201, 0);
    let (mut dx, mut dz, mut fx, mut fz) = (10f32, 10f32, 0f32, 0f32);
    getparint("ninf", &mut ninf);
    getparint("nx", &mut nx);
    getparint("nz", &mut nz);
    getparfloat("dx", &mut dx);
    getparfloat("dz", &mut dz);
    getparfloat("fx", &mut fx);
    getparfloat("fz", &mut fz);
    getparint("npmax", &mut npmax);
    getparint("ns", &mut ns);
    let mut nufile = String::new();
    let vti = getparstring("nufile", &mut nufile) == 0;
    if vti {
        printf!("Since user does not give tilt-field file, the medium is assumed to be VTI\n");
    }
    let mut method = "spline".to_string();
    getparstring("method", &mut method);
    let names = [
        "c11_file",
        "c13_file",
        "c33_file",
        "c15_file",
        "c35_file",
        "c55_file",
        "vp_file",
        "vs_file",
        "epsilon_file",
        "delta_file",
        "rho_file",
    ];
    let mut files = Vec::new();
    for name in names {
        let mut path = name.to_string();
        getparstring(name, &mut path);
        files.push(efopen(&path, "w"));
    }
    let (nx, nz, ninf, npmax) = (nx as usize, nz as usize, ninf as usize, npmax as usize);
    let xs: Vec<f32> = (0..nx).map(|i| fx + i as f32 * dx).collect();
    let mut zs = vec![vec![0f32; nx]; ninf + 1];
    let mut xint = vec![0f32; npmax];
    let mut zint = vec![0f32; npmax];
    for i in 0..=ninf {
        let mut j = 0;
        loop {
            if j >= npmax {
                err!(
                    "The point number on the %ith interface exceeds %i!",
                    i as i32,
                    npmax as i32
                );
            }
            if let Some(x) = fscanf_f32(&stdin()) {
                xint[j] = x;
            }
            if let Some(z) = fscanf_f32(&stdin()) {
                zint[j] = z;
            }
            if zint[j] <= -9999. {
                break;
            }
            j += 1;
        }
        let npt = j;
        if method.starts_with('l') || nx == 1 {
            intlin(
                &xint[..npt],
                &zint[..npt],
                zint[0],
                zint[npt - 1],
                &xs,
                &mut zs[i],
            );
        } else {
            let mut yd = vec![[0f32; 4]; npt];
            if method.starts_with('m') {
                cmonot(&xint[..npt], &zint[..npt], &mut yd);
            } else if method.starts_with('a') {
                cakima(&xint[..npt], &zint[..npt], &mut yd);
            } else if method.starts_with('s') {
                csplin(&xint[..npt], &zint[..npt], &mut yd);
            } else {
                err!("%s is an unknown interpolation method!\n", method);
            }
            intcub(0, &xint[..npt], &yd, &xs, &mut zs[i]);
        }
        for ix in 0..nx {
            if i > 0 && zs[i][ix] < zs[i - 1][ix] {
                err!(
                    "the %ith interface is above the %ith one at position \n\t\t(%g,%g)! \n",
                    i as i32,
                    (i - 1) as i32,
                    xs[ix],
                    zs[i][ix]
                );
            }
        }
    }
    let n = ninf + 1;
    let x0 = layer("x0", n, 0., 0.);
    let z0 = layer("z0", n, 0., 0.);
    let vp00 = layer("vp00", n, 1500., 500.);
    let vs00 = layer("vs00", n, 1500., 500.);
    let eps00 = layer("eps00", n, 0., 0.);
    let delta00 = layer("delta00", n, 0., 0.);
    let rho00 = layer("rho00", n, 1000., 10.);
    let drdx = layer("drdx", n, 0., 0.);
    let drdz = layer("drdz", n, 0., 0.);
    let dvpdx = layer("dvpdx", n, 0., 0.);
    let dvpdz = layer("dvpdz", n, 0., 0.);
    let dvsdx = layer("dvsdx", n, 0., 0.);
    let dvsdz = layer("dvsdz", n, 0., 0.);
    let dedx = layer("dedx", n, 0., 0.);
    let dedz = layer("dedz", n, 0., 0.);
    let dddx = layer("dddx", n, 0., 0.);
    let dddz = layer("dddz", n, 0., 0.);
    let mut nu = vec![0f32; nx * nz];
    if !vti {
        let bytes = match std::fs::read(&nufile) {
            Ok(b) => b,
            Err(_) => err!("error opening nufile=%s!\n", nufile),
        };
        if bytes.len() < nx * nz * 4 {
            err!("error reading nufile=%s!\n", nufile);
        }
        for k in 0..nx * nz {
            nu[k] = f32::from_ne_bytes(bytes[4 * k..4 * k + 4].try_into().unwrap());
        }
    }
    let mut grids = vec![vec![0f32; nx * nz]; 11];
    for ix in 0..nx {
        let mut z = fz;
        for iz in 0..nz {
            let mut chosen = None;
            if z >= zs[ninf][ix] {
                chosen = Some(ninf);
            } else {
                for j in 0..ninf {
                    if z >= zs[j][ix] && z < zs[j + 1][ix] {
                        chosen = Some(j);
                    }
                }
            }
            if let Some(j) = chosen {
                let k = ix * nz + iz;
                let x = xs[ix] - x0[j];
                let zz = z - z0[j];
                grids[6][k] = vp00[j] + x * dvpdx[j] + zz * dvpdz[j];
                grids[7][k] = vs00[j] + x * dvsdx[j] + zz * dvsdz[j];
                grids[8][k] = eps00[j] + x * dedx[j] + zz * dedz[j];
                grids[9][k] = delta00[j] + x * dddx[j] + zz * dddz[j];
                grids[10][k] = rho00[j] + x * drdx[j] + zz * drdz[j];
            }
            z += dz;
        }
    }
    if ns != 0 {
        for ix in 0..nx {
            for j in 1..=ninf {
                let q = (zs[j][ix] - fz) / dz;
                let (nl, nh) = if ns % 2 == 0 {
                    (
                        (q - (ns / 2) as f32 + 1.) as i32,
                        (q + (ns / 2) as f32 + 1.) as i32,
                    )
                } else {
                    (
                        (q - (ns / 2) as f32) as i32,
                        (q + (ns / 2) as f32 + 1.) as i32,
                    )
                };
                for iz in nl + 1..nh {
                    for g in 6..11 {
                        let a = ix * nz + nl as usize;
                        let b = ix * nz + nh as usize;
                        let k = ix * nz + iz as usize;
                        grids[g][k] = grids[g][a]
                            + (iz - nl) as f32 * (grids[g][b] - grids[g][a]) / ns as f32;
                    }
                }
            }
        }
    }
    let mut spar = Stiff2D::default();
    for k in 0..nx * nz {
        let vp = grids[6][k];
        let vs = grids[7][k];
        let e = grids[8][k];
        let d = grids[9][k];
        let rho = grids[10][k];
        if vti {
            let c33 = vp * vp * rho;
            let c55 = (0.25f64 * c33 as f64) as f32;
            let c11 = c33 * (1. + 2. * e);
            let c = c33 * (2. * c55 * (1. + d) - c33 * (1. + 2. * d));
            let b = 2. * c55;
            let c13 = ((-b as f64 + ((b * b - 4. * c) as f64).sqrt()) / 2.) as f32;
            grids[0][k] = c11;
            grids[1][k] = c13;
            grids[2][k] = c33;
            grids[5][k] = c55;
        } else {
            thom2stiffTI(
                vp as f64,
                vs as f64,
                e as f64,
                d as f64,
                e as f64,
                nu[k] as f64,
                &mut spar,
                1,
            );
            grids[0][k] = spar.a1111 as f32 * rho;
            grids[1][k] = spar.a1133 as f32 * rho;
            grids[2][k] = spar.a3333 as f32 * rho;
            grids[3][k] = spar.a1113 as f32 * rho;
            grids[4][k] = spar.a3313 as f32 * rho;
            grids[5][k] = spar.a1313 as f32 * rho;
        }
    }
    for i in 0..11 {
        efwrite(&grids[i], &files[i]);
    }
    exit(CWP_Exit());
}
