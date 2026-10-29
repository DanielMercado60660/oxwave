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
//! VEL2STIFF: par/main/velocity_model_building/vel2stiff.c, version 1.0.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: CWP: Sverre Brandsberg-Dahl 1999. Extended by CWP: Stig-Kyrre
//! Foss 2001 for Sayers parameters; C. M. Sayers, Geophysics 1995, pp. 1933-1935.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " VEL2STIFF - Transforms VELocities, densities, and Thomsen or Sayers   ",
    "\t\tparameters to elastic STIFFnesses \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " vel2stiff  [Required parameters] [Optional Parameters] > stdout\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " vpfile=\tfile with P-wave velocities\t\t\t\t",
    " vsfile=\tfile with S-wave velocities\t\t\t\t",
    " rhofile=\tfile with densities\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " epsfile=\tfile with Thomsen/Sayers epsilon\t\t\t",
    " deltafile=\tfile with Thomsen/Sayers delta\t\t\t \t",
    " gammafile=\tfile with Thomsen/Sayers gamma\t\t\t \t",
    " phi_file=\tangle of axis of symmetry from vertical (radians)\t",
    "\t\t\t\t\t\t\t\t\t",
    " c11_file=c11_file     output filename for c11 values                  ",
    " c13_file=c13_file     output filename for c13 values                  ",
    " c15_file=c15_file     output filename for c15 values                  ",
    " c33_file=c33_file     output filename for c33 values                  ",
    " c35_file=c35_file     output filename for c35 values                  ",
    " c44_file=c44_file     output filename for c44 values                  ",
    " c55_file=c55_file     output filename for c55 values                  ",
    " c66_file=c66_file     output filename for c66 values                  ",
    " \t\t\t\t\t\t\t\t\t",
    " paramtype=1  (1) Thomsen parameters, (0) Sayers parameters(see below) ",
    " \t\t\t\t\t\t\t\t\t",
    " nx=101\tnumber of x samples 2nd (slow) dimension\t\t",
    " nz=101\tnumber of z samples 1st (fast) dimension\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " Transforms velocities, density and Thomsen/Sayers parameters\t\t",
    " epsilon, delta, and gamma into elastic stiffness coefficients.\t",
    " \t\t\t\t\t\t\t\t\t",
    " If only P-wave, S-wave velocities and density is given as input,\t",
    " the model is assumed to be isotropic.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " If files containing Thomsen/Sayers parameters are given, the model\t",
    " will be assumed to have VTI symmetry.\t\t \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " All input files  vpfile, vsfile, rhofile etc. are assumed to consist  ",
    " only of C style binary floating point numbers representing the        ",
    " corresponding  material values of vp, vs, rho etc. Similarly, the output",
    " files consist of the coresponding stiffnesses as C style binary floats. ",
    " If the output files are to be used as input for a modeling program,   ",
    " such as suea2df, then further, the contents are assumed be arrays of  ",
    " floating point numbers of the form of   Array[n2][n1], where the fast ",
    " dimension, dimension 1, represents depth.                             ",
];

use oxwave_core::prelude::*;
use oxwave_f_refaniso::{rottens2D, Stiff2D};

fn required_grid(name: &str, n: usize, label: &str) -> Vec<f32> {
    let fp = efopen(name, "r");
    let mut data = vec![0.0; n];
    if efread(&mut data, &fp) != n {
        err!("error reading %s=%s!\n", label, name);
    }
    fclose(&fp);
    data
}
fn optional_grid(name: &str, n: usize, label: &str, verbose: i32) -> Vec<f32> {
    let mut data = vec![0.0; n];
    if !name.is_empty() {
        let fp = match fopen(name, "r") {
            Some(fp) => fp,
            None => err!("cannot open %sfile=%s", label, name),
        };
        if fread(&mut data, &fp) != n {
            let read_label = if label == "eps" { "epsilon" } else { label };
            err!("error reading %sfile=%s", read_label, name);
        }
        fclose(&fp);
    } else if verbose != 0 {
        let warning_label = if label == "eps" { "epsilon" } else { label };
        warn!(
            "No input file for %s, setting %s to zero!!!",
            warning_label, warning_label
        );
    }
    data
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut vpfile, mut vsfile, mut rhofile) = (String::new(), String::new(), String::new());
    mustgetparstring("vpfile", &mut vpfile);
    mustgetparstring("vsfile", &mut vsfile);
    mustgetparstring("rhofile", &mut rhofile);
    let (mut nx, mut nz, mut verbose, mut paramtype) = (101, 101, 1, 1);
    getparint("nx", &mut nx);
    getparint("nz", &mut nz);
    getparint("verbose", &mut verbose);
    getparint("paramtype", &mut paramtype);
    let (mut epsfile, mut deltafile, mut gammafile, mut phifile) =
        (String::new(), String::new(), String::new(), String::new());
    getparstring("epsfile", &mut epsfile);
    getparstring("deltafile", &mut deltafile);
    getparstring("gammafile", &mut gammafile);
    getparstring("phifile", &mut phifile);
    let n = nx.wrapping_mul(nz) as usize;
    let vp = required_grid(&vpfile, n, "vpfile");
    let vs = required_grid(&vsfile, n, "vsfile");
    let rho = required_grid(&rhofile, n, "rhofile");
    let eps = optional_grid(&epsfile, n, "eps", verbose);
    let delta = optional_grid(&deltafile, n, "delta", verbose);
    let gamma = optional_grid(&gammafile, n, "gamma", verbose);
    // C reads n floats into double storage, leaving the remaining bytes uninitialized.
    let mut phi = vec![0.0f64; n];
    if !phifile.is_empty() {
        let fp = match fopen(&phifile, "r") {
            Some(fp) => fp,
            None => err!("cannot open phifile=%s", phifile),
        };
        let mut bytes = vec![0u8; n * 4];
        if fread(&mut bytes, &fp) != n * 4 {
            err!("error reading phifile=%s", phifile);
        }
        fclose(&fp);
        for (i, chunk) in bytes.chunks(8).enumerate() {
            let mut b = [0u8; 8];
            b[..chunk.len()].copy_from_slice(chunk);
            phi[i] = f64::from_ne_bytes(b);
        }
    } else if verbose != 0 {
        warn!("No input file for phi, setting phi to zero!!!");
    }
    let names = [
        "c11_file", "c13_file", "c33_file", "c44_file", "c66_file", "c15_file", "c35_file",
        "c55_file",
    ];
    let mut filenames = Vec::new();
    let mut files = Vec::new();
    for name in names {
        let mut filename = String::new();
        if getparstring(name, &mut filename) == 0 {
            filename = name.to_string();
        }
        files.push(efopen(&filename, "w"));
        filenames.push(filename);
    }
    checkpars();
    let mut out = vec![vec![0.0f32; n]; 8];
    for ix in 0..nx {
        for iz in 0..nz {
            let i = (ix * nz + iz) as usize;
            let c33 = vp[i] * vp[i] * rho[i];
            let c44 = vs[i] * vs[i] * rho[i];
            let c11 = c33 * (1.0f32 + 2.0f32 * eps[i]);
            let c66 = c44 * (1.0f32 + 2.0f32 * gamma[i]);
            let c13 = if paramtype == 1 {
                let c =
                    c33 * (2.0f32 * c44 * (1.0f32 + delta[i]) - c33 * (1.0f32 + 2.0f32 * delta[i]));
                let b = 2.0f32 * c44;
                ((-(b as f64) + ((b * b - 4.0f32 * c) as f64).sqrt()) / 2.0) as f32
            } else {
                (delta[i] + 1.0f32) * c33 - 2.0f32 * c44
            };
            out[0][i] = c11;
            out[1][i] = c13;
            out[2][i] = c33;
            out[3][i] = c44;
            out[4][i] = c66;
            out[5][i] = 0.0;
            out[6][i] = 0.0;
            out[7][i] = c44;
        }
        for iz in 0..nz {
            let i = (ix * nz + iz) as usize;
            let mut spar = Stiff2D {
                a1111: out[0][i] as f64,
                a1133: out[1][i] as f64,
                a3333: out[2][i] as f64,
                a2323: out[3][i] as f64,
                a1212: out[4][i] as f64,
                a1113: out[5][i] as f64,
                a3313: out[6][i] as f64,
                a1313: out[7][i] as f64,
                a1223: 0.0,
            };
            rottens2D(&mut spar, phi[i]);
            out[0][i] = spar.a1111 as f32;
            out[1][i] = spar.a1133 as f32;
            out[2][i] = spar.a3333 as f32;
            out[3][i] = spar.a2323 as f32;
            out[4][i] = spar.a1212 as f32;
            out[5][i] = spar.a1113 as f32;
            out[6][i] = spar.a3313 as f32;
            out[7][i] = spar.a1313 as f32;
        }
    }
    for i in 0..8 {
        efwrite(&out[i], &files[i]);
    }
    if verbose != 0 {
        for i in 0..5 {
            warn!(
                "Output file for %s parameter plane : %s ",
                &names[i][..3],
                &filenames[i]
            );
        }
    }
    exit(CWP_Exit());
}
