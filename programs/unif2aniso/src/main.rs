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

//! UNIF2ANISO, port of par/main/velocity_model_building/unif2aniso.c (Revision 1.16).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: John Stockwell; based on unif2 by Zhenyue Liu; contributions by
//! Zhu Zhaolin, Tongji University, and Cristiano Calonaci, CINECA.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_f_cubicspline::{cakima, cmonot, csplin};
use oxwave_f_refaniso::{rottens2D, Stiff2D};
use oxwave_linear_interp::intlin;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t ",
    " UNIF2ANISO - generate a 2-D UNIFormly sampled profile of elastic\t",
    "\tconstants from a layered model.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  unif2aniso < infile [Parameters]\t\t\t\t\t",
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
    " vs00=866,1155...,\tS-velocity at each x0,z0 (m/sec)\t\t",
    " rho00=1000,1100,...,\tdensity at each x0,z0 (kg/m^3)\t\t\t",
    " q00=110,120,130,..,\t\tattenuation Q, at each x0,z0 (kg/m^3)\t",
    " \t\t\t\t\t\t\t\t\t",
    " eps00=0,0,0...,\tThomsen or Sayers epsilon\t\t\t",
    " delta00=0,0,0...,\tThomsen or Sayers delta\t\t\t\t",
    " gamma00=0,0,0...,\tThomsen or Sayers gamma\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dqdx=0.0,0.0,...,\tx-derivative of Q (d q/dx)\t\t\t",
    " dqdz=0.0,0.0,...,\tz-derivative of Q (d q/dz)\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " drdx=0.0,0.0,...,\tx-derivative of density (d rho/dx)\t\t",
    " drdz=0.0,0.0,...,\tz-derivative of density (d rho/dz)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dvpdx=0.0,0.0,...,\tx-derivative of P-velocity (dvp/dx)\t\t",
    " dvpdz=0.0,0.0,...,\tz-derivative of P-velocity (dvs/dz)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dvsdx=0.0,0.0,...,\tx-derivative of S-velocity (dvs/dx)\t\t",
    " dvsdz=0.0,0.0,...,\tz-derivative of S-velocity (dvs/dz)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dedx=0.0,0.0,...,\tx-derivative of epsilon (de/dx)\t\t\t",
    " dedz=0.0,0.0,...,\tz-derivative of epsilon with depth z (de/dz)\t",
    " \t\t\t\t\t\t\t\t\t",
    " dddx=0.0,0.0,...,\tx-derivative of delta (dd/dx)\t\t\t",
    " dddz=0.0,0.0,...,\tz-derivative of delta (dd/dz)\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " dgdz=0.0,0.0,...,\tx-derivative of gamma (dg/dz)\t\t\t",
    " dgdx=0.0,0.0,...,\tz-derivative of gamma (dg/dx)\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " phi00=0,0,..., \trotation angle(s) in each layer\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " ...output filenames \t\t\t\t\t\t\t",
    " c11_file=c11_file\toutput filename for c11 values\t \t\t",
    " c13_file=c13_file\toutput filename for c13 values\t \t\t",
    " c15_file=c15_file\toutput filename for c15 values\t \t\t",
    " c33_file=c33_file\toutput filename for c33 values\t \t\t",
    " c35_file=c35_file\toutput filename for c35 values\t \t\t",
    " c44_file=c44_file\toutput filename for c44 values\t \t\t",
    " c55_file=c55_file\toutput filename for c55 values\t \t\t",
    " c66_file=c66_file\toutput filename for c66 values\t \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " rho_file=rho_file\toutput filename for density values \t\t",
    " q_file=q_file\t\toutput filename for Q values\t \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " paramtype=1   =1 Thomsen parameters, =0 Sayers parameters(see below)\t",
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
    " The algorithm works by transforming the P-wavespeed , S-wavespeed,\t",
    " density and the Thomsen or Sayers parameters epsilon, delta, and gamma",
    " into elastic stiffness coefficients. Furthermore, the\tuser can specify",
    " rotations, phi, to the elasticity tensor in each layer.\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Common ranges of Thomsen parameters are\t\t\t\t",
    "  epsilon:  0.0 -> 0.5\t\t\t\t\t\t\t",
    "  delta:   -0.2 -> 0.4\t\t\t\t\t\t\t",
    "  gamma:\t0.0 -> 0.4\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " If only P-wave, S-wave velocities and density is given as input,\t",
    " the model is, by definition,  isotropic.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " If files containing Thomsen/Sayers parameters are given, the model\t",
    " will be assumed to have VTI symmetry.\t\t \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example using test input file generating feature:\t\t\t",
    " unif2aniso tfile=testfilename  produces a 5 interface demonstration model",
    " unif2aniso < testfilename \t\t\t\t\t\t",
    " ximage < c11_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c13_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c15_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c33_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c35_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c44_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c55_file n1=100 n2=100\t\t\t\t\t",
    " ximage < c66_file n1=100 n2=100\t\t\t\t\t",
    " ximage < rho_file n1=100 n2=100\t\t\t\t\t",
    " ximage < q_file   n1=100 n2=100\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

const SAMPLE: &[&str] = &[
    "   0	 0",
    "1000	 0",
    "   1	-99999",
    "   0	200",
    "1000	200",
    "   1	-99999",
    "   0	400",
    " 100.	400",
    " 400	560",
    " 550.	560",
    "1000.	360",
    "   1	-99999",
    "   0	800",
    " 100.	800",
    " 400.	750.",
    " 500.	780",
    " 550.	780",
    "1000.	500",
    "   1.   -99999",
    "   0	950",
    " 100	950",
    " 400	900.",
    " 500	800.",
    " 550	800.",
    " 750	920",
    "1000	930",
    "   1	-99999",
    "   0	1000",
    "1000	1000",
    "   1	-99999",
];
fn floats(name: &str, defaults: Vec<f32>) -> Vec<f32> {
    let mut v = defaults;
    getparfloat(name, &mut v);
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
    let (mut ninf, mut nx, mut nz, mut npmax) = (5, 100, 100, 201);
    let (mut dx, mut dz, mut fx, mut fz) = (10.0f32, 10.0f32, 0.0f32, 0.0f32);
    let mut method = "linear".to_string();
    let mut paramtype = 1;
    getparint("ninf", &mut ninf);
    getparint("nx", &mut nx);
    getparint("nz", &mut nz);
    getparfloat("dx", &mut dx);
    getparfloat("dz", &mut dz);
    getparfloat("fx", &mut fx);
    getparfloat("fz", &mut fz);
    getparint("npmax", &mut npmax);
    getparstring("method", &mut method);
    getparint("paramtype", &mut paramtype);
    getparint("verbose", &mut paramtype);
    let names = [
        "c11_file", "c13_file", "c33_file", "c15_file", "c35_file", "c55_file", "c44_file",
        "c66_file", "rho_file", "q_file",
    ];
    let mut files = Vec::new();
    for name in names {
        let mut path = name.to_string();
        getparstring(name, &mut path);
        files.push(efopen(&path, "w"));
    }
    let (nx, nz, ni, np) = (nx as usize, nz as usize, ninf as usize, npmax as usize);
    let mut xs = vec![0.0; nx];
    for i in 0..nx {
        xs[i] = fx + (i as f32) * dx;
    }
    let mut zs = vec![vec![0.0; nx]; ni + 1];
    let mut xi = vec![0.0; np];
    let mut zi = vec![0.0; np];
    for i in 0..=ni {
        let mut j = 0;
        loop {
            if j >= np {
                err!(
                    "The point number on the %ith interface exceeds %i!",
                    i as i32,
                    np as i32
                );
            }
            if let Some(x) = fscanf_f32(&stdin()) {
                xi[j] = x;
                if let Some(z) = fscanf_f32(&stdin()) {
                    zi[j] = z;
                }
            }
            if zi[j] <= -9999.0 {
                break;
            }
            j += 1;
        }
        let n = j;
        if method.starts_with('l') || nx == 1 {
            intlin(&xi[..n], &zi[..n], zi[0], zi[n - 1], &xs, &mut zs[i]);
        } else {
            let mut yd = vec![[0.0; 4]; n];
            match method.chars().next() {
                Some('m') => cmonot(&xi[..n], &zi[..n], &mut yd),
                Some('a') => cakima(&xi[..n], &zi[..n], &mut yd),
                Some('s') => csplin(&xi[..n], &zi[..n], &mut yd),
                _ => err!("%s is an unknown interpolation method!\n", method.as_str()),
            }
            intcub(0, &xi[..n], &yd, &xs, &mut zs[i]);
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
    let len = ni + 1;
    let x0 = floats("x0", vec![0.0; len]);
    let z0 = floats("z0", vec![0.0; len]);
    let vp = floats("vp00", (0..len).map(|i| (1500 + 500 * i) as f32).collect());
    let vs = floats(
        "vs00",
        (0..len)
            .map(|i| (3.0f64.sqrt() * (1500 + 500 * i) as f64) as f32)
            .collect(),
    );
    let rho = floats("rho00", (0..len).map(|i| (1000 + 10 * i) as f32).collect());
    let q = floats("q00", (0..len).map(|i| (100 + 10 * i) as f32).collect());
    let eps = floats("eps00", vec![0.0; len]);
    let delta = floats("delta00", vec![0.0; len]);
    let gamma = floats("gamma00", vec![0.0; len]);
    let dvpdx = floats("dvpdx", vec![0.0; len]);
    let dvpdz = floats("dvpdz", vec![0.0; len]);
    let dvsdx = floats("dvsdx", vec![0.0; len]);
    let dvsdz = floats("dvsdz", vec![0.0; len]);
    let dqdx = floats("dqdx", vec![0.0; len]);
    let dqdz = floats("dqdz", vec![0.0; len]);
    let drdx = floats("drdx", vec![0.0; len]);
    let drdz = floats("drdz", vec![0.0; len]);
    let dedx = floats("dedx", vec![0.0; len]);
    let dedz = floats("dedz", vec![0.0; len]);
    let dddx = floats("dddx", vec![0.0; len]);
    let dddz = floats("dddz", vec![0.0; len]);
    let dgdx = floats("dgdx", vec![0.0; len]);
    let dgdz = floats("dgdz", vec![0.0; len]);
    let mut phi = vec![0.0f64; len];
    getpardouble("phi00", &mut phi);
    let mut out = vec![vec![0.0f32; nx * nz]; 10];
    for ix in 0..nx {
        let mut j = 1usize;
        let mut z = fz;
        for iz in 0..nz {
            let k = if z >= zs[ni][ix] {
                ni
            } else if z <= zs[1][ix] {
                0
            } else {
                while z > zs[j][ix] {
                    j += 1;
                }
                j -= 1;
                j
            };
            let x = xs[ix] - x0[k];
            let zz = z - z0[k];
            let v = vp[k] + x * dvpdx[k] + zz * dvpdz[k];
            let s = vs[k] + x * dvsdx[k] + zz * dvsdz[k];
            let r = rho[k] + x * drdx[k] + zz * drdz[k];
            let qq = q[k] + x * dqdx[k] + zz * dqdz[k];
            let e = eps[k] + x * dedx[k] + zz * dedz[k];
            let d = delta[k] + x * dddx[k] + zz * dddz[k];
            let g = gamma[k] + x * dgdx[k] + zz * dgdz[k];
            let c33 = v * v * r;
            let c44 = s * s * r;
            let c11 = c33 * (1.0 + 2.0 * e);
            let c66 = c44 * (1.0 + 2.0 * g);
            let c13 = if paramtype == 1 {
                let c = c33 * (2.0 * c44 * (1.0 + d) - c33 * (1.0 + 2.0 * d));
                let b = 2.0 * c44;
                ((-(b as f64) + ((b * b - 4.0 * c) as f64).sqrt()) / 2.0) as f32
            } else {
                (d + 1.0) * c33 - 2.0 * c44
            };
            let mut a = Stiff2D {
                a1111: c11 as f64,
                a1133: c13 as f64,
                a3333: c33 as f64,
                a2323: c44 as f64,
                a1212: c66 as f64,
                a1313: c44 as f64,
                ..Default::default()
            };
            rottens2D(&mut a, phi[k]);
            let vals = [
                a.a1111, a.a1133, a.a3333, a.a1113, a.a3313, a.a1313, a.a2323, a.a1212,
            ];
            let pos = ix * nz + iz;
            for h in 0..8 {
                out[h][pos] = if vals[h].abs() < 1.0e-3 {
                    0.0
                } else {
                    vals[h] as f32
                };
            }
            out[8][pos] = r;
            out[9][pos] = qq;
            z += dz;
        }
    }
    for i in 0..10 {
        efwrite(&out[i], &files[i]);
    }
    for fp in &files {
        efclose(fp);
    }
    exit(CWP_Exit());
}
