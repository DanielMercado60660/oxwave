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
//! STIFF2VEL - transforms 2D elastic stiffnesses to velocities and anisotropy.
//! Port of par/main/velocity_model_building/stiff2vel.c (version 1.0).
//! Credits: Aramco: Chris Liner, 9/25/2005 (based on vel2stiff.c).

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " STIFF2VEL - Transforms 2D elastic stiffnesses to (vp,vs,epsilon,delta) ",
    " \t\t\t\t\t\t\t\t\t",
    " stiff2vel  nx=  nz=  [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " nx=\t\t \t number of x samples (2nd dimension)\t\t",
    " nz=\t\t \t number of z samples (1st dimension)\t\t",
    " rho_file='rho.bin'     input file containing rho(x,z)\t\t\t",
    " c11_file='c11.bin'     input file containing c11(x,z)\t\t\t",
    " c13_file='c13.bin'     input file containing c13(x,z)\t\t\t",
    " c33_file='c33.bin'     input file containing c33(x,z)\t\t\t",
    " c44_file='c44.bin'     input file containing c44(x,z)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " vp_file='vp.bin'\t output file containing P-wave velocities\t",
    " vs_file='vs.bin'\t output file containing S-wave velocities\t",
    " rho_file='rho.bin'\t output file containing densities\t\t",
    " eps_file='eps.bin'\t output file containing Thomsen epsilon\t       \t",
    " delta_file='delta.bin' output file containing Thomsen delta\t       \t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " 1. All quantities in MKS units\t\t\t\t\t",
    " 2. Isotropy implied by c11(x,z)=c33(x,z)=0 \t\t\t",
    " 3. Vertical symmetry axis is assumed.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut nx = 0;
    let mut nz = 0;
    mustgetparint("nx", &mut nx);
    mustgetparint("nz", &mut nz);
    let mut rho_file = String::new();
    let mut c11_file = String::new();
    let mut c13_file = String::new();
    let mut c33_file = String::new();
    let mut c44_file = String::new();
    let mut vp_file = String::new();
    let mut vs_file = String::new();
    let mut eps_file = String::new();
    let mut delta_file = String::new();
    let mut verbose = 1;
    if getparstring("rho_file", &mut rho_file) == 0 { rho_file = "rho.bin".into(); }
    if getparstring("c11_file", &mut c11_file) == 0 { c11_file = "c11.bin".into(); }
    if getparstring("c13_file", &mut c13_file) == 0 { c13_file = "c13.bin".into(); }
    if getparstring("c33_file", &mut c33_file) == 0 { c33_file = "c33.bin".into(); }
    if getparstring("c44_file", &mut c44_file) == 0 { c44_file = "c44.bin".into(); }
    if getparstring("vp_file", &mut vp_file) == 0 { vp_file = "vp.bin".into(); }
    if getparstring("vs_file", &mut vs_file) == 0 { vs_file = "vs.bin".into(); }
    if getparstring("eps_file", &mut eps_file) == 0 { eps_file = "eps.bin".into(); }
    if getparstring("delta_file", &mut delta_file) == 0 { delta_file = "delta.bin".into(); }
    getparint("verbose", &mut verbose);
    checkpars();
    let n = nx.wrapping_mul(nz) as usize;
    let mut rho = vec![0.0f32; n];
    let mut c11 = vec![0.0f32; n];
    let mut c13 = vec![0.0f32; n];
    let mut c33 = vec![0.0f32; n];
    let mut c44 = vec![0.0f32; n];
    let mut vp = vec![0.0f32; n];
    let mut vs = vec![0.0f32; n];
    let mut eps = vec![0.0f32; n];
    let mut delta = vec![0.0f32; n];
    let rhofp = efopen(&rho_file, "r");
    if efread(&mut rho, &rhofp) != n { err!("error reading rho_file=%s!\n", rho_file); }
    let c11fp = efopen(&c11_file, "r");
    if efread(&mut c11, &c11fp) != n { err!("error reading c11_file=%s!\n", c11_file); }
    let c13fp = efopen(&c13_file, "r");
    if efread(&mut c13, &c13fp) != n { err!("error reading c13_file=%s!\n", c13_file); }
    let c33fp = efopen(&c33_file, "r");
    if efread(&mut c33, &c33fp) != n { err!("error reading c33_file=%s!\n", c33_file); }
    let c44fp = efopen(&c44_file, "r");
    if efread(&mut c44, &c44fp) != n { err!("error reading c44_file=%s!\n", c44_file); }
    fclose(&rhofp); fclose(&c11fp); fclose(&c13fp); fclose(&c33fp); fclose(&c44fp);
    let vpfp = fopen(&vp_file, "w").unwrap();
    let vsfp = fopen(&vs_file, "w").unwrap();
    let epsfp = fopen(&eps_file, "w").unwrap();
    let deltafp = fopen(&delta_file, "w").unwrap();
    for ix in 0..nx {
        for iz in 0..nz {
            let i = (ix * nz + iz) as usize;
            vp[i] = ((c33[i] / rho[i]) as f64).sqrt() as f32;
            vs[i] = ((c44[i] / rho[i]) as f64).sqrt() as f32;
            eps[i] = (c11[i] - c33[i]) / (2.0f32 * c33[i]);
            let mut tmp = (c13[i] + c44[i]) * (c13[i] + c44[i]);
            tmp = tmp - (c33[i] - c44[i]) * (c33[i] - c44[i]);
            delta[i] = tmp / (2.0f32 * c33[i] * (c33[i] - c44[i]));
        }
    }
    efwrite(&vp, &vpfp);
    efwrite(&vs, &vsfp);
    efwrite(&eps, &epsfp);
    efwrite(&delta, &deltafp);
    if verbose != 0 {
        warn!("Output file for vp : %s ", vp_file);
        warn!("Output file for vs : %s ", vs_file);
        warn!("Output file for epsilon : %s ", eps_file);
        warn!("Output file for delta : %s ", delta_file);
    }
    exit(CWP_Exit());
}
