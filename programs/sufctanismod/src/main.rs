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

//! SUFCTANISMOD - Flux-corrected transport for anisotropic elastic waves.
//! Port of su/main/synthetics_waveforms_testpatterns/sufctanismod.c (Revision 1.8).
//! Author: Tong Fei, CWP, Colorado School of Mines, Dec 1993.
//! Additional features: Stig-Kyrre Foss, CWP, Oct 2001.
//! Modifications: Chris Liner, U Houston, Mar 2010.

#![allow(unused_mut, unused_variables, unused_assignments)]

use oxwave_core::prelude::*;
mod solver;
use solver::*;

// alloc2float provides one contiguous allocation with [outer][inner] access.
#[derive(Clone, Default)]
struct Grid {
    values: Vec<f32>,
    n1: usize,
}

impl Grid {
    fn new(n1: i32, n2: i32) -> Self {
        Self {
            values: vec![0.0f32; n1 as usize * n2 as usize],
            n1: n1 as usize,
        }
    }

    // C's row pointers address a contiguous allocation. An inner coordinate
    // may cross a row boundary while still referring to an allocated element.
    fn c_index(&self, ix: i32, iz: i32) -> Option<usize> {
        if self.n1 == 0 || ix < 0 || ix as usize >= self.values.len() / self.n1 {
            return None;
        }
        let index = ix as i64 * self.n1 as i64 + iz as i64;
        if index < 0 || index >= self.values.len() as i64 {
            None
        } else {
            Some(index as usize)
        }
    }

    fn receiver(&self, ix: i32, depth: i32) -> f32 {
        let Some(index) = self.c_index(ix, depth) else {
            // The C read is undefined only once it leaves the allocation.
            err!(
                "Receiver depth=%i is outside the allocated grid. Exit",
                depth
            );
        };
        self.values[index]
    }
}

impl std::ops::Index<usize> for Grid {
    type Output = [f32];
    fn index(&self, i: usize) -> &[f32] {
        &self.values[i * self.n1..(i + 1) * self.n1]
    }
}

impl std::ops::IndexMut<usize> for Grid {
    fn index_mut(&mut self, i: usize) -> &mut [f32] {
        &mut self.values[i * self.n1..(i + 1) * self.n1]
    }
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUFCTANISMOD - Flux-Corrected Transport correction applied to the 2D",
    "\t  elastic wave equation for finite difference modeling in \t",
    "\t  anisotropic media\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sufctanismod > outfile [optional parameters]\t\t\t\t",
    "\t\toutfile is the final wavefield snapshot x-component\t",
    "\t\tx-component of wavefield snapshot is in snapshotx.data\t",
    "\t\ty-component of wavefield snapshot is in snapshoty.data\t",
    "\t\tz-component of wavefield snapshot is in snapshotz.data\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Output Files:\t\t\t\t\t\t",
    " reflxfile=\treflection seismogram file name for x-component\t\t",
    "\t\tno output produced if no name specified\t \t\t",
    " reflyfile=\treflection seismogram file name for y-component\t\t",
    "\t\tno output produced if no name specified\t \t\t",
    " reflzfile=\treflection seismogram file name for z-component\t\t",
    "\t\tno output produced if no name specified\t \t\t",
    " vspxfile=\tVSP seismogram file name for x-component\t\t",
    "\t\tno output produced if no name specified\t \t\t",
    " vspyfile=\tVSP seismogram file name for y-component\t\t",
    "\t\tno output produced if no name specified\t \t\t",
    " vspzfile=\tVSP seismogram file name for z-component\t\t",
    "\t\tno output produced if no name specified\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suhead=1      To get SU-header output seismograms (else suhead=0)\t",
    "\t\t\t\t\t\t\t\t\t",
    " New parameter:\t\t\t\t\t\t\t",
    "\t               \t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " mt=1          number of time steps per output snapshot  \t\t",
    " dofct=1 \t1 do the FCT correction\t\t\t\t\t",
    "\t\t0 do not do the FCT correction \t\t\t\t",
    " FCT Related parameters:\t\t\t\t\t\t",
    " eta0=0.03\tdiffusion coefficient\t\t\t\t\t",
    "\t\ttypical values ranging from 0.008 to 0.06\t\t",
    "\t\tabout 0.03 for the second-order method \t\t\t",
    "\t\tabout 0.012 for the fourth-order method \t\t",
    " eta=0.04\tanti-diffusion coefficient \t\t\t\t",
    "\t\ttypical values ranging from 0.008 to 0.06\t\t",
    "\t\tabout 0.04 for the second-order method  \t\t",
    "\t\tabout 0.015 for the fourth-order method \t\t",
    " fctxbeg=0 \tx coordinate to begin applying the FCT correction\t",
    " fctzbeg=0 \tz coordinate to begin applying the FCT correction\t",
    " fctxend=nx \tx coordinate to stop applying the FCT correction\t",
    " fctzend=nz \tz coordinate to stop applying the FCT correction\t",
    " \t\t\t\t\t\t\t\t\t",
    " deta0dx=0.0\tgradient of eta0 in x-direction  d(eta0)/dx\t\t",
    " deta0dz=0.0\tgradient of eta0 in z-direction  d(eta0)/dz\t\t",
    " detadx=0.0\tgradient of eta in x-direction \t d(eta)/dx\t\t",
    " detadz=0.0\tgradient of eta in z-direction \t d(eta)/dz\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " General Parameters:\t\t\t\t\t\t\t",
    " order=2\t2 second-order finite-difference \t\t\t",
    "\t\t4 fourth-order finite-difference \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " nt=200        number of time steps \t\t\t \t\t",
    " dt=0.004\ttime step  \t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " nx=100 \tnumber of grid points in x-direction \t\t\t",
    " nz=100 \tnumber of grid points in z-direction \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " dx=0.02\tspatial step in x-direction \t\t\t\t",
    " dz=0.02\tspatial step in z-direction \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sx=nx/2\tsource x-coordinate (in gridpoints)\t\t\t",
    " sz=nz/2\tsource z-coordinate (in gridpoints)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " fpeak=20\tpeak frequency of the wavelet \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " receiverdepth=sz  depth of horizontal receivers (in gridpoints)      ",
    " vspnx=sx\t\t\tx grid loc of vsp\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " verbose=0     silent operation\t\t\t\t\t\t\t",
    "\t\t\t\t=1 for diagnostic messages, =2 for more\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " wavelet=1\t1 AKB wavelet\t\t\t\t\t\t",
    " \t\t2 Ricker wavelet \t\t\t\t\t",
    "\t\t3 impulse \t\t\t\t\t\t",
    "\t\t4 unity \t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " isurf=2\t1 absorbing surface condition \t\t\t\t",
    "\t\t2 free surface condition \t\t\t\t",
    "\t\t3 zero surface condition \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " source=1\t1 point source \t\t\t\t\t\t",
    " \t\t2 sources are located on a given refelector \t        ",
    "\t\t\t(two horizontal and one dipping reflectors) \t",
    " \t\t3 sources are located on a given dipping refelector     ",
    "\t\t\t\t\t\t\t\t\t",
    " sfile= \tthe name of input source file, if no name specified then",
    "\t\tuse default source location. (source=1 or 2) \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Density and Elastic Parameters:\t\t\t\t\t",
    " dfile= \tthe name of input density file,                         ",
    "               if no name specified then                             ",
    "\t\tassume a linear density profile with ...\t\t",
    " rho00=2.0\tdensity at (0, 0) \t\t\t\t\t",
    " drhodx=0.0\tdensity gradient in x-direction  d(rho)/dx\t\t",
    " drhodz=0.0\tdensity gradient in z-direction  d(rho)/dz\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " afile= \tname of input elastic param.  (c11) aa file, if no name ",
    "\t\tspecified then, assume a linear profile with ...\t",
    " aa00=2.0\telastic parameter at (0, 0) \t\t\t\t",
    " daadx=0.0\tparameter gradient in x-direction  d(aa)/dx\t\t",
    " daadz=0.0\tparameter gradient in z-direction  d(aa)/dz\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " cfile= \tname of input elastic param. (c33)  cc file, if no name ",
    "\t\tspecified then, assume a linear profile with ...\t",
    " cc00=2.0\telastic parameter at (0, 0) \t\t\t\t",
    " dccdx=0.0\tparameter gradient in x-direction  d(cc)/dx\t\t",
    " dccdz=0.0\tparameter gradient in z-direction  d(cc)/dz\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " ffile= \tname of input elastic param.  (c13) ff file, if no name ",
    "\t\tspecified then, assume a linear profile with ...\t",
    " ff00=2.0\telastic parameter at (0, 0) \t\t\t\t",
    " dffdx=0.0\tparameter gradient in x-direction  d(ff)/dx\t\t",
    " dffdz=0.0\tparameter gradient in z-direction  d(ff)/dz\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " lfile= \tname of input elastic param.  (c44) ll file, if no name ",
    "\t\tspecified then, assume a linear profile with ...\t",
    " ll00=2.0\telastic parameter at (0, 0) \t\t\t\t",
    " dlldx=0.0\tparameter gradient in x-direction  d(ll)/dx\t\t",
    " dlldz=0.0\tparameter gradient in z-direction  d(ll)/dz\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " nfile= \tname of input elastic param. (c66)  nn file, if no name ",
    "\t\tspecified then, assume a linear profile with ...\t",
    " nn00=2.0\telastic parameter at (0, 0) \t\t\t\t",
    " dnndx=0.0\tparameter gradient in x-direction  d(nn)/dx\t\t",
    " dnndz=0.0\tparameter gradient in z-direction  d(nn)/dz\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optimizations:\t\t\t\t\t\t\t",
    " The moving boundary option permits the user to restrict the computations",
    " of the wavefield to be confined to a specific range of spatial coordinates.",
    " The boundary of this restricted area moves with the wavefield\t\t",
    " movebc=0\t0 do not use moving boundary optimization\t\t",
    "\t\t1 use moving boundaries\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut nt: i32 = 250;
    getparint("nt", &mut nt);
    let mut nx: i32 = 300;
    getparint("nx", &mut nx);
    let mut nz: i32 = 250;
    getparint("nz", &mut nz);
    let mut sx: i32 = nx / 2;
    getparint("sx", &mut sx);
    if sx > nx {
        err!("Source sx=%i is greater than nx=%i. Exit", sx, nx);
    }
    let mut sz: i32 = nz / 2;
    getparint("sz", &mut sz);
    // C passes an int to %f; on this build its bits print as 0.000000.
    if sz > nz {
        err!("Source sz=%f is greater than nz=%i. Exit", 0.0, nx);
    }
    let mut receiverdepth: i32 = sz;
    getparint("receiverdepth", &mut receiverdepth);
    let mut vspnx: i32 = sx;
    getparint("vspnx", &mut vspnx);
    let mut impulse: i32 = 0;
    getparint("impulse", &mut impulse);
    let mut source: i32 = 0;
    getparint("source", &mut source);
    let mut isurf: i32 = 1;
    getparint("isurf", &mut isurf);
    let mut dofct: i32 = 0;
    getparint("dofct", &mut dofct);
    let mut fctxbeg: i32 = 0;
    getparint("fctxbeg", &mut fctxbeg);
    let mut fctzbeg: i32 = 0;
    getparint("fctzbeg", &mut fctzbeg);
    let mut fctxend: i32 = nx;
    getparint("fctxend", &mut fctxend);
    let mut fctzend: i32 = nz;
    getparint("fctzend", &mut fctzend);
    let mut indexux: i32 = 0;
    getparint("indexux", &mut indexux);
    let mut indexuy: i32 = 0;
    getparint("indexuy", &mut indexuy);
    let mut indexuz: i32 = 0;
    getparint("indexuz", &mut indexuz);
    let mut wavelet: i32 = 1;
    getparint("wavelet", &mut wavelet);
    let mut movebc: i32 = 0;
    getparint("movebc", &mut movebc);
    let mut order: i32 = 2;
    getparint("order", &mut order);
    let mut suhead: i32 = 1;
    getparint("suhead", &mut suhead);
    let mut mt: i32 = 1;
    getparint("mt", &mut mt);
    let mut verbose: i32 = 0;
    getparint("verbose", &mut verbose);
    let mut dx: f32 = 0.02;
    getparfloat("dx", &mut dx);
    let mut dz: f32 = 0.02;
    getparfloat("dz", &mut dz);
    let mut dt: f32 = 0.002;
    getparfloat("dt", &mut dt);
    let mut fpeak: f32 = 20.0;
    getparfloat("fpeak", &mut fpeak);
    let mut aa00: f32 = 2.0;
    getparfloat("aa00", &mut aa00);
    let mut cc00: f32 = 2.0;
    getparfloat("cc00", &mut cc00);
    let mut ff00: f32 = 2.0;
    getparfloat("ff00", &mut ff00);
    let mut ll00: f32 = 2.0;
    getparfloat("ll00", &mut ll00);
    let mut nn00: f32 = 2.0;
    getparfloat("nn00", &mut nn00);
    let mut daadx: f32 = 0.0;
    getparfloat("daadx", &mut daadx);
    let mut daadz: f32 = 0.0;
    getparfloat("daadz", &mut daadz);
    let mut dccdx: f32 = 0.0;
    getparfloat("dccdx", &mut dccdx);
    let mut dccdz: f32 = 0.0;
    getparfloat("dccdz", &mut dccdz);
    let mut dffdx: f32 = 0.0;
    getparfloat("dffdx", &mut dffdx);
    let mut dffdz: f32 = 0.0;
    getparfloat("dffdz", &mut dffdz);
    let mut dlldx: f32 = 0.0;
    getparfloat("dlldx", &mut dlldx);
    let mut dlldz: f32 = 0.0;
    getparfloat("dlldz", &mut dlldz);
    let mut dnndx: f32 = 0.0;
    getparfloat("dnndx", &mut dnndx);
    let mut dnndz: f32 = 0.0;
    getparfloat("dnndz", &mut dnndz);
    let mut drhodx: f32 = 0.0;
    getparfloat("drhodx", &mut drhodx);
    let mut drhodz: f32 = 0.0;
    getparfloat("drhodz", &mut drhodz);
    let mut rho00: f32 = 1.0;
    getparfloat("rho00", &mut rho00);
    let mut eta0: f32 = 0.03;
    getparfloat("eta0", &mut eta0);
    let mut eta: f32 = 0.04;
    getparfloat("eta", &mut eta);
    let mut deta0dx: f32 = 0.0;
    getparfloat("deta0dx", &mut deta0dx);
    let mut deta0dz: f32 = 0.0;
    getparfloat("deta0dz", &mut deta0dz);
    let mut detadx: f32 = 0.0;
    getparfloat("detadx", &mut detadx);
    let mut detadz: f32 = 0.0;
    getparfloat("detadz", &mut detadz);
    let mut sfile = String::new();
    getparstring("sfile", &mut sfile);
    let mut dfile = String::new();
    getparstring("dfile", &mut dfile);
    let mut afile = String::new();
    getparstring("afile", &mut afile);
    let mut cfile = String::new();
    getparstring("cfile", &mut cfile);
    let mut ffile = String::new();
    getparstring("ffile", &mut ffile);
    let mut lfile = String::new();
    getparstring("lfile", &mut lfile);
    let mut nfile = String::new();
    getparstring("nfile", &mut nfile);
    getparstring("nfile", &mut nfile);
    let mut reflxfile = String::new();
    getparstring("reflxfile", &mut reflxfile);
    let mut reflyfile = String::new();
    getparstring("reflyfile", &mut reflyfile);
    let mut reflzfile = String::new();
    getparstring("reflzfile", &mut reflzfile);
    let mut vspxfile = String::new();
    getparstring("vspxfile", &mut vspxfile);
    let mut vspyfile = String::new();
    getparstring("vspyfile", &mut vspyfile);
    let mut vspzfile = String::new();
    getparstring("vspzfile", &mut vspzfile);

    // Arrays are [x][z]; seismograms are [receiver][time].
    let mut u = Grid::new(nz, nx);
    let mut v = u.clone();
    let mut w = u.clone();
    let mut e11 = u.clone();
    let mut e33 = u.clone();
    let mut e23 = u.clone();
    let mut e12 = u.clone();
    let mut e13 = u.clone();
    let mut xzsource = u.clone();
    let mut aa = u.clone();
    let mut cc = u.clone();
    let mut ff = u.clone();
    let mut ll = u.clone();
    let mut nn = u.clone();
    let mut rho = u.clone();
    // The spike helper writes element 1 even when nt=1.
    let mut fx = vec![0.0f32; (nt as usize).max(2)];
    let mut fy = fx.clone();
    let mut fz = fx.clone();
    let mut t = vec![0.0f32; nt as usize];

    let mut refl_x = if !reflxfile.is_empty() {
        Grid::new(nt, nx)
    } else {
        Grid::default()
    };
    let mut refl_y = if !reflyfile.is_empty() {
        Grid::new(nt, nx)
    } else {
        Grid::default()
    };
    let mut refl_z = if !reflzfile.is_empty() {
        Grid::new(nt, nx)
    } else {
        Grid::default()
    };
    let mut vsp_x = if !vspxfile.is_empty() {
        Grid::new(nt, nz)
    } else {
        Grid::default()
    };
    let mut vsp_y = if !vspyfile.is_empty() {
        Grid::new(nt, nz)
    } else {
        Grid::default()
    };
    let mut vsp_z = if !vspzfile.is_empty() {
        Grid::new(nt, nz)
    } else {
        Grid::default()
    };

    warn!(
        "Memory usage approximately %f MB",
        (15 * nz * nx * 4 + 4 * nt * 4 + 6 * nt * nx * 4) as f32 / (1024 * 1024) as f32
    );
    for it in 0..nt as usize {
        t[it] = it as f32 * dt;
    }
    // The component flags select identical time functions independently.
    for (index, force) in [(indexux, &mut fx), (indexuy, &mut fy), (indexuz, &mut fz)] {
        if index != 0 {
            if wavelet == 1 {
                tforce_akb(nt, force, dt, fpeak);
            }
            if wavelet == 2 {
                tforce_ricker(nt, force, dt, fpeak);
            }
            if wavelet == 3 {
                tforce_spike(nt, force, dt, fpeak);
            }
            if wavelet == 4 {
                tforce_unit(nt, force, dt, fpeak);
            }
        }
    }
    warn!("Source function set");

    let mut vmax = 0.0f32;
    read_parameter(nx, nz, dx, dz, rho00, drhodx, drhodz, &dfile, &mut rho);
    read_parameter(nx, nz, dx, dz, aa00, daadx, daadz, &afile, &mut aa);
    read_parameter(nx, nz, dx, dz, cc00, dccdx, dccdz, &cfile, &mut cc);
    read_parameter(nx, nz, dx, dz, ff00, dffdx, dffdz, &ffile, &mut ff);
    read_parameter(nx, nz, dx, dz, ll00, dlldx, dlldz, &lfile, &mut ll);
    read_parameter(nx, nz, dx, dz, nn00, dnndx, dnndz, &nfile, &mut nn);
    warn!("Parameter grids read");
    for ix in 0..nx as usize {
        for iz in 0..nz as usize {
            rho[ix][iz] = (1.0 / rho[ix][iz] as f64) as f32;
            if cc[ix][iz] * rho[ix][iz] > vmax {
                vmax = cc[ix][iz] * rho[ix][iz];
            }
        }
    }
    vmax = (vmax as f64).sqrt() as f32;
    warn!("vmax = %f", vmax);
    if dz <= dx {
        warn!(
            "stability check:  vmax*dt/(sqrt(2)*min(dx,dz)) = %f (should be < 1)",
            (vmax * dt) as f64 / (2.0f64.sqrt() * dz as f64)
        );
    } else {
        warn!(
            "stability check:  vmax*dt/(sqrt(2)*min(dx,dz)) = %f (should be < 1)",
            (vmax * dt) as f64 / (2.0f64.sqrt() * dx as f64)
        );
    }
    if !sfile.is_empty() {
        read_parameter(nx, nz, dx, dz, aa00, daadx, daadz, &sfile, &mut xzsource);
    } else {
        locate_source(nx, nz, sx, sz, &mut xzsource, source);
    }

    // These are unchecked fopen calls in C, not efopen. A directory yields
    // NULL without a diagnostic. Retain None until a snapshot is written:
    // with nt=0 the loop is skipped, the final grid is still written to
    // stdout, and the program succeeds even if all three opens failed.
    let outfpx = fopen("snapshotx.data", "w");
    let outfpy = fopen("snapshoty.data", "w");
    let outfpz = fopen("snapshotz.data", "w");
    warn!("Snapshot files open");
    let (mut mbx1, mut mbz1, mut mbx2, mut mbz2) = (0, 0, 0, 0);
    for it in 0..nt {
        if verbose == 2 {
            fprintf!(stderr(), "it= %d\n", it);
        }
        moving_bc(
            it, nx, nz, sx, sz, dx, dz, impulse, movebc, &t, vmax, &mut mbx1, &mut mbz1, &mut mbx2,
            &mut mbz2,
        );
        moving_fctbc(
            mbx1,
            mbz1,
            mbx2,
            mbz2,
            0,
            0,
            nx,
            nz,
            &mut fctxbeg,
            &mut fctzbeg,
            &mut fctxend,
            &mut fctzend,
        );
        // As in C, order is parsed but anis_solver2 is always used.
        anis_solver2(
            it, &mut u, &mut v, &mut w, &mut e11, &mut e33, &mut e23, &mut e12, &mut e13, &aa, &cc,
            &ff, &ll, &nn, &rho, &xzsource, &fx, &fy, &fz, dx, dz, dt, nx, nz, dofct, isurf, eta0,
            eta, deta0dx, deta0dz, detadx, detadz, mbx1, mbx2, mbz1, mbz2,
        );
        for ix in 0..nx as usize {
            if !reflxfile.is_empty() {
                refl_x[ix][it as usize] = u.receiver(ix as i32, receiverdepth);
            }
        }
        for ix in 0..nx as usize {
            if !reflyfile.is_empty() {
                refl_y[ix][it as usize] = v.receiver(ix as i32, receiverdepth);
            }
        }
        for ix in 0..nx as usize {
            if !reflzfile.is_empty() {
                if w.receiver(ix as i32, receiverdepth) > 1.0 {
                    warn!("problems at it=%d x=%d \n", it, ix);
                }
                refl_z[ix][it as usize] = w.receiver(ix as i32, receiverdepth);
            }
        }
        if vspnx >= 0 && vspnx < nx {
            for iz in 0..nz as usize {
                if !vspxfile.is_empty() {
                    vsp_x[iz][it as usize] = u[vspnx as usize][iz];
                }
            }
            for iz in 0..nz as usize {
                if !vspyfile.is_empty() {
                    vsp_y[iz][it as usize] = v[vspnx as usize][iz];
                }
            }
            for iz in 0..nz as usize {
                if !vspzfile.is_empty() {
                    vsp_z[iz][it as usize] = w[vspnx as usize][iz];
                }
            }
        }
        // arm64 integer division by zero produces zero, hence remainder it.
        if it.checked_rem(mt).unwrap_or(it) == 0 {
            // Passing a failed open to fwrite is undefined in C. Report
            // that path through the normal SU error API instead.
            write_grid(
                &u,
                outfpx
                    .as_ref()
                    .unwrap_or_else(|| err!("Cannot open file=%s", "snapshotx.data")),
            );
            write_grid(
                &v,
                outfpy
                    .as_ref()
                    .unwrap_or_else(|| err!("Cannot open file=%s", "snapshoty.data")),
            );
            write_grid(
                &w,
                outfpz
                    .as_ref()
                    .unwrap_or_else(|| err!("Cannot open file=%s", "snapshotz.data")),
            );
            if verbose >= 1 {
                warn!("snap at it=%i", it);
            }
        }
    }
    write_grid(&u, &stdout());
    warn!("Source function set");
    output_seismogram(&reflxfile, nx, sx, sz, nt, dx, dt, suhead, &refl_x, true);
    output_seismogram(&reflyfile, nx, sx, sz, nt, dx, dt, suhead, &refl_y, true);
    output_seismogram(&reflzfile, nx, sx, sz, nt, dx, dt, suhead, &refl_z, true);
    output_seismogram(&vspxfile, nz, sx, sz, nt, dx, dt, suhead, &vsp_x, true);
    // C leaves these two streams open until exit. Keep their buffers alive:
    // closing early changes the final data when output paths alias snapshots.
    let _outvspy = output_seismogram(&vspyfile, nz, sx, sz, nt, dx, dt, suhead, &vsp_y, false);
    let _outvspz = output_seismogram(&vspzfile, nz, sx, sz, nt, dx, dt, suhead, &vsp_z, false);
    exit(CWP_Exit());
}

fn write_grid(data: &Grid, fp: &CFile) {
    fwrite(&data.values, fp);
}

fn read_parameter(
    nx: i32,
    nz: i32,
    dx: f32,
    dz: f32,
    p00: f32,
    dpdx: f32,
    dpdz: f32,
    file: &str,
    pp: &mut Grid,
) {
    if !file.is_empty() {
        let Some(infp) = fopen(file, "r") else {
            fprintf!(stderr(), "Cannot open file=%s\n", file);
            exit(1);
        };
        // C reads the contiguous allocation in one call and ignores short reads.
        fread(&mut pp.values, &infp);
    } else {
        for ix in 0..nx {
            for iz in 0..nz {
                pp[ix as usize][iz as usize] =
                    p00 + dpdx * (ix - 1) as f32 * dx + dpdz * (iz - 1) as f32 * dz;
            }
        }
    }
}

fn output_seismogram(
    file: &str,
    nx: i32,
    sx: i32,
    sz: i32,
    nt: i32,
    dx: f32,
    dt: f32,
    suhead: i32,
    data: &Grid,
    close: bool,
) -> Option<CFile> {
    if file.is_empty() {
        return None;
    }
    let Some(fp) = fopen(file, "w") else {
        fprintf!(stderr(), "Cannot open file=%s\n", file);
        exit(1);
    };
    if suhead == 1 {
        su_output(nx, sx, sz, nt, dx, dt, &fp, data);
    } else {
        write_grid(data, &fp);
    }
    if close {
        fclose(&fp);
        None
    } else {
        Some(fp)
    }
}

fn su_output(
    nx: i32,
    sx: i32,
    sdepth: i32,
    ns: i32,
    dx: f32,
    dt: f32,
    outputfile: &CFile,
    data: &Grid,
) {
    let mut sudata = Trace::new();
    sudata.trid = 1;
    sudata.sx = sx;
    sudata.sdepth = sdepth;
    sudata.ns = ns as u16;
    // The target C conversion truncates to int before narrowing to short.
    sudata.dt = (1000000.0f32 * dt) as i32 as i16 as u16;
    sudata.d2 = dx;
    let (mut tracl, mut tracr) = (0, 0);
    for ix in 0..nx {
        tracr += 1;
        tracl += 1;
        sudata.offset = (ix as f32 * dx - sx as f32) as i32;
        sudata.tracl = tracl;
        sudata.tracr = tracr;
        for it in 0..ns {
            sudata.data[it as usize] = data[ix as usize][it as usize];
        }
        fputtr(outputfile, &mut sudata);
    }
}
