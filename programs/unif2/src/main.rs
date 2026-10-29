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

//! UNIF2 - uniformly sampled velocity profile from a layered model.
//! Port of par/main/velocity_model_building/unif2.c (Revision 1.12).
//! Credits: CWP: Zhenyue Liu, 1994; John Stockwell, 1994.
use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_f_cubicspline::{cakima, cmonot, csplin};
use oxwave_linear_interp::intlin;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t ",
    " UNIF2 - generate a 2-D UNIFormly sampled velocity profile from a layered",
    "  \t model. In each layer, velocity is a linear function of position.",
    " \t\t\t\t\t\t\t\t\t",
    "  unif2 < infile > outfile [parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
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
    " x0=0.0,0.0,..., \tdistance x at which v00 is specified\t\t",
    " z0=0.0,0.0,..., \tdepth z at which v00 is specified\t\t",
    " v00=1500,2000,2500...,\tvelocity at each x0,z0 (m/sec)\t\t",
    " dvdx=0.0,0.0,...,\tderivative of velocity with distance x (dv/dx)\t",
    " dvdz=0.0,0.0,...,\tderivative of velocity with depth z (dv/dz)\t",
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
    " piecewise continuous velocity model with a flat surface on top. The surface",
    " and each successive boundary between media are represented by a list of",
    " selected x z pairs written column form. The first and last x values must",
    " be the same for all boundaries. Use the entry   1.0  -99999  to separate",
    " entries for successive boundaries. No boundary may cross another. Note",
    " that the choice of the method of interpolation may cause boundaries \t",
    " to cross that do not appear to cross in the input data file.\t\t",
    " The number of interfaces is specified by the parameter \"ninf\". This ",
    " number does not include the top surface of the model. The input data\t",
    " format is the same as a CSHOT model file with all comments removed.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example using test input file generating feature:\t\t\t",
    " unif2 tfile=testfilename    produces a 5 interface demonstration model",
    " unif2 < testfilename | psimage n1=100 n2=100 d1=10 d2=10 | ...\t",
    "\t\t\t\t\t\t\t\t\t",
];
const SAMPLE: &[&str] = &[
"   0         0", "1000         0", "   1    -99999", "   0       200", "1000       200", "   1    -99999", "   0       400", " 100.      400", " 400       560", " 550.      560", "1000.      360", "   1    -99999", "   0       800", " 100.      800", " 400.      750.", " 500.      780", " 550.      780", "1000.      500", "   1.   -99999", "   0       950", " 100       950", " 400       900.", " 500       800.", " 550       800.", " 750       920", "1000       930", "   1    -99999", "   0      1000", "1000      1000", "   1    -99999"];
fn main() {
 let argv = args(); initargs(&argv); requestdoc(1, SDOC);
 let mut tfile = String::new(); getparstring("tfile", &mut tfile);
 if !tfile.is_empty() { let fp = efopen(&tfile, "w"); for s in SAMPLE { fprintf!(fp, "%s\n", *s); } efclose(&fp); exit(0); }
 let (mut ninf, mut nx, mut nz, mut npmax) = (5, 100, 100, 201);
 let (mut dx, mut dz, mut fx, mut fz) = (10.0f32, 10.0f32, 0.0f32, 0.0f32);
 getparint("ninf", &mut ninf); getparint("nx", &mut nx); getparint("nz", &mut nz);
 getparfloat("dx", &mut dx); getparfloat("dz", &mut dz); getparfloat("fx", &mut fx); getparfloat("fz", &mut fz); getparint("npmax", &mut npmax);
 let mut method = String::from("linear"); getparstring("method", &mut method);
 let (ni, nx, nz, np) = (ninf as usize, nx as usize, nz as usize, npmax as usize);
 let mut v = vec![0.0f32; nz]; let mut x0 = vec![0.0f32; ni+1]; let mut z0 = x0.clone(); let mut v00 = x0.clone(); let mut dvdx = x0.clone(); let mut dvdz = x0.clone();
 let mut xint = vec![0.0f32; np]; let mut zint = xint.clone(); let mut xs = vec![0.0f32; nx]; let mut zs = vec![vec![0.0f32; nx]; ni+1];
 for ix in 0..nx { xs[ix] = fx + ix as f32 * dx; }
 let input = stdin();
 for i in 0..=ni {
  let mut j = 0usize;
  loop { if j >= np { err!("The point number on the %ith interface exceeds %i!", i as i32, npmax); }
   if let Some(x) = fscanf_f32(&input) { xint[j] = x; } if let Some(z) = fscanf_f32(&input) { zint[j] = z; }
   if zint[j] <= -9999.0 { break; } j += 1;
  }
  let npt = j;
  if method.starts_with('l') || nx == 1 { intlin(&xint[..npt], &zint[..npt], zint[0], zint[npt-1], &xs, &mut zs[i]); }
  else if method.starts_with('m') || method.starts_with('a') || method.starts_with('s') {
   let mut zind = vec![[0.0f32; 4]; npt];
   if method.starts_with('m') { cmonot(&xint[..npt], &zint[..npt], &mut zind); } else if method.starts_with('a') { cakima(&xint[..npt], &zint[..npt], &mut zind); } else { csplin(&xint[..npt], &zint[..npt], &mut zind); }
   intcub(0, &xint[..npt], &zind, &xs, &mut zs[i]);
  } else { err!("%s is an unknown interpolation method!\n", method); }
  for ix in 0..nx { if i > 0 && zs[i][ix] < zs[i-1][ix] { err!("the %ith interface is above the %ith one at position \n\t\t(%g,%g)! \n", i as i32, (i-1) as i32, xs[ix], zs[i][ix]); } }
 }
 getparfloat("x0", &mut x0); getparfloat("z0", &mut z0);
 if getparfloat("v00", &mut v00) == 0 { for i in 0..=ni { v00[i] = (1500.0f64 + (500*i) as f64) as f32; } }
 getparfloat("dvdx", &mut dvdx); getparfloat("dvdz", &mut dvdz); checkpars();
 for ix in 0..nx { let mut j = 1usize; let mut z = fz;
  for iz in 0..nz { let layer = if z >= zs[ni][ix] { ni } else if z <= zs[1][ix] { 0 } else { while z > zs[j][ix] { j += 1; } j -= 1; j };
   v[iz] = v00[layer] + (xs[ix]-x0[layer])*dvdx[layer] + (z-z0[layer])*dvdz[layer]; z += dz;
  } efwrite(&v, &stdout());
 }
 exit(CWP_Exit());
}
