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

//! TETRAMOD - TETRAhedron MODel builder.
//! Port of tetra/main/tetramod.c (Revision 1.5).
//! Credits: CWP: Zhaobo Meng, 1996.
//! Reference: Zhaobo Meng and Norman Bleistein, Wavefront Construction (WF)
//! Ray Tracing in Tetrahedral Models, CWP report 251, 1997.

#![allow(unused_parens, unused_assignments, unused_mut)]
use oxwave_core::prelude::*;
use oxwave_lib_tetra::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t \t\t\t",
    " TETRAMOD - TETRAhedron MODel builder. In each layer, velocity gradient",
    "       is constant or a 2-D grid; horizons could be a uniform grid and/or",
    "       added by a 2-D grid specified.   \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tetramod [parameters] > tetrafile\t\t \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "                                                                       ",
    " nxhz=\t\tnumber of samples (2nd dimension) for horizons\t\t",
    " nyhz=\t\tnumber of samples (1st dimension) for horizons\t\t",
    " hzfile=\toutput xhz,yhz,zhz,v0hz,v1hz for viewer3                ",
    "                                                                       ",
    " Optional parameters:\t\t\t\t\t\t\t",
    "                                                                       ",
    " xmin=0       \tx of the lower left point in the model\t        \t",
    " ymin=0        y of the lower left point in the model                  ",
    " xmax=2        x of the upper right point in the model                 ",
    " ymax=2        y of the upper right point in the model                 ",
    " zmax=(max z)  max z in the model                                      ",
    " blt=1.0       bottom layer thickness                                  ",
    " nhz=1\t\tnumber of layers in the model (except the model base)\t",
    " ficth=-1      ficticious horizons (no velocity interpolation based on them)",
    " The following four numbers define the four corners for 3-D model;\t",
    " z00=0,0.6,1.2,... \tz at (xmin,ymin) on each horizon\t\t",
    " z01=0,0.6,1.2,... \tz at (xmin,ymax) on each horizon\t\t",
    " z10=0,0.6,1.2,... \tz at (xmax,ymin) on each horizon                ",
    " z11=0,0.6,1.2,... \tz at (xmax,ymax) on each horizon                ",
    " v00=1.0,2.0,3.0,...\tv at (xmin,ymin) on each horizon\t\t",
    " v01=1.0,2.0,3.0,...\tv at (xmin,ymax) on each horizon\t\t",
    " v10=1.0,2.0,3.0,...\tv at (xmax,xmin) on each horizon                ",
    " v11=1.0,2.0,3.0,...\tv at (xmax,xmax) on each horizon                ",
    " dvdz00=0,0,0,...   \tdvdz at (xmin,ymin) on each horizon\t\t",
    " dvdz01=0,0,0,...   \tdvdz at (xmin,ymax) on each horizon\t\t",
    " dvdz10=0,0,0,...   \tdvdz at (xmax,xmin) on each horizon             ",
    " dvdz11=0,0,0,...   \tdvdz at (xmax,xmax) on each horizon             ",
    " x0file=   \tx grid for horizon 0 \t\t\t\t\t",
    " y0file=   \ty grid for horizon 0 \t\t\t\t\t",
    " z0file=\tz grid for horizon 0 added to 4-z interpolation   \t",
    " v0file=\tv grid for horizon 0 added to 4-v interpolation   \t",
    " dvdz0file=\tdvdz grid for horizon 00 added to 4-dvdz interpolation\t",
    "\t\t\t\t\t\t\t\t\t",
    " x1file=\tx grid for horizon 1\t\t\t\t\t",
    " y1file=\ty grid for horizon 1\t\t\t\t\t",
    " z1file=\tz grid for horizon 1 added to 4-z interpolation  \t",
    " v1file=\tv grid for horizon 1 added to 4-v interpolation  \t",
    " dvdz1file=\tdvdz grid for horizon 1 added to 4-dvdz interpolation\t",
    " ...\t\tfor horizon  #2, #3 #4, etc... \t\t\t\t",
    " verbose=0\t=1 print some useful iinformation\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Remarks:\t\t\t\t\t\t\t\t",
    " TETRAMOD defines its own grammar to describe a 3-D model (including a ",
    " tetra file for ray tracer SUTETRARAY and a horizon file for VIEWER3   ",
    "\t\t\t\t\t\t\t\t\t",
    " Disclaimer:                                                           ",
    " This is a research code that will take considerable work to get into  ",
    " the form of a a production level 3D migration code. The code is       ",
    " offered as is, along with tetramod and sukdmig3d, to provide a starting",
    " point for researchers who wish to write their own 3D migration codes.",
    "                                                                       ",
];

#[derive(Clone, Copy, Default)]
struct Point {
    x: [f32; 3],
    s: f32,
    n: [f32; 3],
}
#[derive(Clone, Copy, Default)]
struct Facet {
    ip: [i32; 3],
    cn: [f32; 3],
    itetra: [i32; 2],
    area: f32,
    ct: [f32; 10],
}
#[derive(Clone, Copy, Default)]
struct Tetra {
    ip: [i32; 4],
    ifacet: [i32; 4],
    ireg: i32,
    gs: [f32; 3],
    v: f32,
}

fn main() {
    let argv = args();
    let mut nhz: i32 = 0;
    let mut nxhz: i32 = 0;
    let mut nyhz: i32 = 0;
    let mut dxhz: f32 = 0.0;
    let mut dyhz: f32 = 0.0;
    let mut xhz: Vec<Grid2> = Vec::new();
    let mut yhz: Vec<Grid2> = Vec::new();
    let mut zhz: Vec<Grid2> = Vec::new();
    let mut v0hz: Vec<Grid2> = Vec::new();
    let mut v0hztem: Vec<Grid2> = Vec::new();
    let mut v1hz: Vec<Grid2> = Vec::new();
    let mut dvdzhz: Vec<Grid2> = Vec::new();
    let mut dvdzhztem: Vec<Grid2> = Vec::new();
    let mut n1hz: Vec<Grid2> = Vec::new();
    let mut n2hz: Vec<Grid2> = Vec::new();
    let mut n3hz: Vec<Grid2> = Vec::new();
    let mut n1ll: Grid2 = Grid2::default();
    let mut n1ur: Grid2 = Grid2::default();
    let mut n2ll: Grid2 = Grid2::default();
    let mut n2ur: Grid2 = Grid2::default();
    let mut n3ll: Grid2 = Grid2::default();
    let mut n3ur: Grid2 = Grid2::default();
    let mut npoint: i32 = 0;
    let mut point: Vec<Point> = Vec::new();
    let mut nfacet: i32 = 0;
    let mut facet: Vec<Facet> = Vec::new();
    let mut ntetra: i32 = 0;
    let mut tetra: Vec<Tetra> = Vec::new();
    let mut z00: Vec<f32> = Vec::new();
    let mut z01: Vec<f32> = Vec::new();
    let mut z10: Vec<f32> = Vec::new();
    let mut z11: Vec<f32> = Vec::new();
    let mut v00: Vec<f32> = Vec::new();
    let mut v01: Vec<f32> = Vec::new();
    let mut v10: Vec<f32> = Vec::new();
    let mut v11: Vec<f32> = Vec::new();
    let mut dvdz00: Vec<f32> = Vec::new();
    let mut dvdz01: Vec<f32> = Vec::new();
    let mut dvdz10: Vec<f32> = Vec::new();
    let mut dvdz11: Vec<f32> = Vec::new();
    let mut xfile: String = String::new();
    let mut yfile: String = String::new();
    let mut zfile: String = String::new();
    let mut vfile: String = String::new();
    let mut dvdzfile: String = String::new();
    let mut hzfile: String = String::new();
    let hzfp: CFile;
    let mut zflat: Vec<Grid2> = Vec::new();
    let mut vtemp: f32 = 0.0;
    let mut ficth: Vec<i32> = Vec::new();
    let mut ifict: i32 = 0;
    let mut nficth: i32 = 0;
    let mut xmin: f32 = 0.0;
    let mut ymin: f32 = 0.0;
    let mut xmax: f32 = 0.0;
    let mut ymax: f32 = 0.0;
    let mut zmax: f32 = ((-999999.9f64) as f32);
    let mut zmin: f32 = 0.0f32;
    let mut a: [[f32; 9]; 9] = [[0.0; 9]; 9];
    let mut b: [f32; 9] = [0.0; 9];
    let mut ixhz: i32 = 0;
    let mut iyhz: i32 = 0;
    let mut ihz: i32 = 0;
    let mut iihz: Vec<i32> = Vec::new();
    let mut iiihz: i32 = 0;
    let mut verbose: i32 = 0;
    let mut fz0: f32 = 0.0;
    let mut fz1: f32 = 0.0;
    let mut fv0: f32 = 0.0;
    let mut fv1: f32 = 0.0;
    let mut fdvdz0: f32 = 0.0;
    let mut fdvdz1: f32 = 0.0;
    let mut beta: f32 = 0.0;
    let mut alpha: f32 = 0.0;
    let mut dis0: f32 = 0.0f32;
    let mut dis1: f32 = 0.0f32;
    let mut dis2: f32 = 0.0f32;
    let mut i000: i32 = 0;
    let mut i001: i32 = 0;
    let mut i010: i32 = 0;
    let mut i100: i32 = 0;
    let mut i110: i32 = 0;
    let mut i101: i32 = 0;
    let mut i011: i32 = 0;
    let mut i111: i32 = 0;
    let mut if1: i32 = 0;
    let mut if2: i32 = 0;
    let mut if3: i32 = 0;
    let mut if4: i32 = 0;
    let mut if5: i32 = 0;
    let mut if6: i32 = 0;
    let mut if7: i32 = 0;
    let mut if8: i32 = 0;
    let mut if9: i32 = 0;
    let mut if10: i32 = 0;
    let mut if11: i32 = 0;
    let mut if12: i32 = 0;
    let mut it1: i32 = 0;
    let mut it2: i32 = 0;
    let mut it3: i32 = 0;
    let mut it4: i32 = 0;
    let mut it5: i32 = 0;
    let mut it6: i32 = 0;
    let mut ntrap: i32 = 0;
    let mut nhztrap: i32 = 0;
    let mut nhztrap1: i32 = 0;
    let mut nhzfacetx: i32 = 0;
    let mut nfacetx: i32 = 0;
    let mut nhzfacety: i32 = 0;
    let mut nfacety: i32 = 0;
    let mut iwh: i32 = 0;
    let mut ip: i32 = 0;
    let mut nhzpoint: i32 = 0;
    let mut ihz0: i32 = 0;
    let mut blt: f32 = 0.0;
    let mut ip0: i32 = 0;
    let mut ip1: i32 = 0;
    let mut ip2: i32 = 0;
    let mut ip3: i32 = 0;
    let mut itetra: i32 = 0;
    let mut jumphz: i32 = 0;
    let mut if1p: i32 = 0;
    let mut if2p: i32 = 0;
    let mut ifs: i32 = 0;
    // Hook up getpar.
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get parameters.
    if !(getparint("nhz", &mut nhz) != 0) {
        nhz = 1;
    }
    if !(getparint("nxhz", &mut nxhz) != 0) {
        err!("Must specify nxhz");
    }
    if !(getparint("nyhz", &mut nyhz) != 0) {
        err!("Must specify nyhz");
    }
    if !(getparstring("hzfile", &mut hzfile) != 0) {
        err!("Must specify hzfile");
    }
    if !(getparint("verbose", &mut verbose) != 0) {
        verbose = 0;
    }
    nficth = countparval("ficth");
    ficth = vec![0; (nficth as u64) as usize];
    getparint("ficth", &mut ficth);
    if ((nficth == 1) && (ficth[0] == (-1))) {
        nficth = 0;
    }
    if (verbose != 0) {
        fprintf!(stderr(), "nficth=%d\n", nficth);
        ifict = 0;
        while (ifict < nficth) {
            fprintf!(stderr(), "ficth[%d]=%d", ifict, ficth[ifict as usize]);
            ifict += 1;
        }
    }
    if !(getparfloat("blt", &mut blt) != 0) {
        blt = 1.0f32;
    }
    if (verbose != 0) {
        warn!("hzfile=%s\n", hzfile);
    }
    if !(getparfloat("xmin", &mut xmin) != 0) {
        xmin = 0.0f32;
    }
    if !(getparfloat("xmax", &mut xmax) != 0) {
        xmax = 2.0f32;
    }
    if !(getparfloat("ymin", &mut ymin) != 0) {
        ymin = 0.0f32;
    }
    if !(getparfloat("ymax", &mut ymax) != 0) {
        ymax = 2.0f32;
    }
    if !(getparfloat("zmax", &mut zmax) != 0) {
        zmax = 0.0f32;
    }
    hzfp = match fopen(&hzfile, "w") {
        Some(fp) => fp,
        None => err!("Can not open hzfile"),
    };
    if (verbose != 0) {
        warn!("nhz=%d", nhz);
        warn!(
            "xmin=%g\nxmax=%g\nymin=%g\nymax=%g\n",
            (xmin as f64),
            (xmax as f64),
            (ymin as f64),
            (ymax as f64)
        );
    }
    // Negative nhz makes the C ealloc1int fail after the unchecked allocations.
    if nhz < 0 {
        unsafe extern "C" {
            fn __error() -> *mut i32;
        }
        unsafe {
            *__error() = 12;
        }
        syserr!("%s: malloc failed", "ealloc.c");
    }
    z00 = vec![0.0f32; (nhz as u64) as usize];
    z01 = vec![0.0f32; (nhz as u64) as usize];
    z10 = vec![0.0f32; (nhz as u64) as usize];
    z11 = vec![0.0f32; (nhz as u64) as usize];
    v00 = vec![0.0f32; (nhz as u64) as usize];
    v01 = vec![0.0f32; (nhz as u64) as usize];
    v10 = vec![0.0f32; (nhz as u64) as usize];
    v11 = vec![0.0f32; (nhz as u64) as usize];
    dvdz00 = vec![0.0f32; (nhz as u64) as usize];
    dvdz01 = vec![0.0f32; (nhz as u64) as usize];
    dvdz10 = vec![0.0f32; (nhz as u64) as usize];
    dvdz11 = vec![0.0f32; (nhz as u64) as usize];
    iihz = vec![0; (nhz as u64) as usize];
    // Default depths: ihz*0.6.
    if !(getparfloat("z00", &mut z00) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            z00[ihz as usize] = (((ihz as f64) * 0.6f64) as f32);
            ihz += 1;
        }
    }
    if !(getparfloat("z01", &mut z01) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            z01[ihz as usize] = (((ihz as f64) * 0.6f64) as f32);
            ihz += 1;
        }
    }
    if !(getparfloat("z10", &mut z10) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            z10[ihz as usize] = (((ihz as f64) * 0.6f64) as f32);
            ihz += 1;
        }
    }
    if !(getparfloat("z11", &mut z11) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            z11[ihz as usize] = (((ihz as f64) * 0.6f64) as f32);
            ihz += 1;
        }
    }
    if (verbose != 0) {
        ihz = 0;
        while (ihz < nhz) {
            warn!(
                "layer %d z: %f %f %f %f",
                ihz,
                (z00[ihz as usize] as f64),
                (z01[ihz as usize] as f64),
                (z10[ihz as usize] as f64),
                (z11[ihz as usize] as f64)
            );
            ihz += 1;
        }
    }
    // Default velocities: ihz+1 (km/s).
    if !(getparfloat("v00", &mut v00) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            v00[ihz as usize] = ((ihz + 1) as f32);
            ihz += 1;
        }
    }
    if !(getparfloat("v01", &mut v01) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            v01[ihz as usize] = ((ihz + 1) as f32);
            ihz += 1;
        }
    }
    if !(getparfloat("v10", &mut v10) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            v10[ihz as usize] = ((ihz + 1) as f32);
            ihz += 1;
        }
    }
    if !(getparfloat("v11", &mut v11) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            v11[ihz as usize] = ((ihz + 1) as f32);
            ihz += 1;
        }
    }
    if (verbose != 0) {
        ihz = 0;
        while (ihz < nhz) {
            warn!(
                "layer %d v: %f %f %f %f",
                ihz,
                (v00[ihz as usize] as f64),
                (v01[ihz as usize] as f64),
                (v10[ihz as usize] as f64),
                (v11[ihz as usize] as f64)
            );
            ihz += 1;
        }
    }
    if !(getparfloat("dvdz00", &mut dvdz00) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            dvdz00[ihz as usize] = 0.0f32;
            ihz += 1;
        }
    }
    if !(getparfloat("dvdz01", &mut dvdz01) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            dvdz01[ihz as usize] = 0.0f32;
            ihz += 1;
        }
    }
    if !(getparfloat("dvdz10", &mut dvdz10) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            dvdz10[ihz as usize] = 0.0f32;
            ihz += 1;
        }
    }
    if !(getparfloat("dvdz11", &mut dvdz11) != 0) {
        ihz = 0;
        while (ihz < nhz) {
            dvdz11[ihz as usize] = 0.0f32;
            ihz += 1;
        }
    }
    if (verbose != 0) {
        ihz = 0;
        while (ihz < nhz) {
            warn!(
                "layer %d dvdz: %f %f %f %f",
                ihz,
                (dvdz00[ihz as usize] as f64),
                (dvdz01[ihz as usize] as f64),
                (dvdz10[ihz as usize] as f64),
                (dvdz11[ihz as usize] as f64)
            );
            ihz += 1;
        }
    }
    // Horizon grids need not be uniform. Validate their dimensions.
    warn!("nxhz=%d", nxhz);
    warn!("nyhz=%d", nyhz);
    if (nxhz <= 1) {
        err!("nxhz must >1");
    }
    if (nyhz <= 2) {
        err!("nyhz must >2");
    }
    if (((nyhz / 2) * 2) == nyhz) {
        err!("nyhz must be an odd number\n");
    }
    n1hz = vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    n2hz = vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    n3hz = vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    n1ll = Grid2::new(((nxhz - 1) as u64) as usize, ((nyhz - 1) as u64) as usize);
    n1ur = Grid2::new(((nxhz - 1) as u64) as usize, ((nyhz - 1) as u64) as usize);
    n2ll = Grid2::new(((nxhz - 1) as u64) as usize, ((nyhz - 1) as u64) as usize);
    n2ur = Grid2::new(((nxhz - 1) as u64) as usize, ((nyhz - 1) as u64) as usize);
    n3ll = Grid2::new(((nxhz - 1) as u64) as usize, ((nyhz - 1) as u64) as usize);
    n3ur = Grid2::new(((nxhz - 1) as u64) as usize, ((nyhz - 1) as u64) as usize);
    xhz = vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    yhz = vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    zhz = vec![
        Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize);
        ((nhz + 1) as u64) as usize
    ];
    v0hz = vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    v0hztem =
        vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    v1hz = vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    dvdzhz =
        vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    dvdzhztem =
        vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    zflat = vec![Grid2::new((nxhz as u64) as usize, (nyhz as u64) as usize); (nhz as u64) as usize];
    // These spacings are used only when xhz/yhz grids are not supplied.
    dxhz = ((xmax - xmin) / ((nxhz - 1) as f32));
    dyhz = ((ymax - ymin) / ((nyhz - 1) as f32));
    // Read the optional coordinate and additive depth/velocity/gradient grids.
    // C's getparstring repoints each filename into its argument buffer; the
    // next horizon's sprintf then writes there and can corrupt later argument
    // names when filenames are short. Keep owned strings and read the named
    // files: reproducing that undefined aliasing is deliberately not attempted.
    ihz = 0;
    while (ihz < nhz) {
        {
            xfile = sprintf!("x%dfile", ihz);
            yfile = sprintf!("y%dfile", ihz);
            zfile = sprintf!("z%dfile", ihz);
            vfile = sprintf!("v%dfile", ihz);
            dvdzfile = sprintf!("dvdz%dfile", ihz);
            if (verbose == 1) {
                warn!("%s %s %s %s %s", xfile, yfile, zfile, vfile, dvdzfile);
            }
            if (getparstring(&xfile.clone(), &mut xfile) != 0) {
                if ((read_grid(&mut xhz[ihz as usize], &xfile) as u64) != ((nxhz * nyhz) as u64)) {
                    err!("Can not read in %s\n", xfile);
                }
                warn!("read in %s", xfile);
                iyhz = 0;
                while (iyhz < nyhz) {
                    {
                        if (((xhz[ihz as usize][iyhz as usize][0] - xmin) as f64).abs() > 1.0E-4f64)
                        {
                            err!("xhz[%d][%d][0] != xmin, not allowed", ihz, iyhz);
                        } else {
                            xhz[ihz as usize][iyhz as usize][0] = xmin;
                        }
                        if (((xhz[ihz as usize][iyhz as usize][(nxhz - 1) as usize] - xmax) as f64)
                            .abs()
                            > 1.0E-4f64)
                        {
                            err!(
                                "xhz[%d][%d][%d]=%e != xmax=%e,not allowed",
                                ihz,
                                iyhz,
                                (nxhz - 1),
                                (xhz[ihz as usize][iyhz as usize][(nxhz - 1) as usize] as f64),
                                (xmax as f64)
                            );
                        } else {
                            xhz[ihz as usize][iyhz as usize][(nxhz - 1) as usize] = xmax;
                        }
                    }
                    iyhz += 1;
                }
            } else {
                warn!("Using uniform Cartesian grids for x-horizon%d", ihz);
                ixhz = 0;
                while (ixhz < nxhz) {
                    iyhz = 0;
                    while (iyhz < nyhz) {
                        xhz[ihz as usize][iyhz as usize][ixhz as usize] =
                            (xmin + (dxhz * (ixhz as f32)));
                        iyhz += 1;
                    }
                    ixhz += 1;
                }
            }
            if (getparstring(&yfile.clone(), &mut yfile) != 0) {
                if ((read_grid(&mut yhz[ihz as usize], &yfile) as u64) != ((nyhz * nxhz) as u64)) {
                    err!("Can not read in %s\n", yfile);
                }
                warn!("read in %s", yfile);
                ixhz = 0;
                while (ixhz < nxhz) {
                    {
                        if (((yhz[ihz as usize][0][ixhz as usize] - ymin) as f64).abs() > 1.0E-4f64)
                        {
                            err!(
                                "yhz[%d][0][%d]=%e != ymin=%e,not allowed",
                                ihz,
                                ixhz,
                                (yhz[ihz as usize][0][ixhz as usize] as f64),
                                (ymin as f64)
                            );
                        } else {
                            yhz[ihz as usize][0][ixhz as usize] = ymin;
                        }
                        if (((yhz[ihz as usize][(nyhz - 1) as usize][ixhz as usize] - ymax) as f64)
                            .abs()
                            > 1.0E-4f64)
                        {
                            err!(
                                "yhz[%d][%d][%d]=%e!=ymax=%e,not allowed",
                                ihz,
                                (nyhz - 1),
                                ixhz,
                                (yhz[ihz as usize][(nyhz - 1) as usize][ixhz as usize] as f64),
                                (ymax as f64)
                            );
                        } else {
                            yhz[ihz as usize][(nyhz - 1) as usize][ixhz as usize] = ymax;
                        }
                    }
                    ixhz += 1;
                }
            } else {
                warn!("Using uniform Cartesian grids for y-horizon%d", ihz);
                ixhz = 0;
                while (ixhz < nxhz) {
                    iyhz = 0;
                    while (iyhz < nyhz) {
                        yhz[ihz as usize][iyhz as usize][ixhz as usize] =
                            (ymin + ((iyhz as f32) * dyhz));
                        iyhz += 1;
                    }
                    ixhz += 1;
                }
            }
            if (getparstring(&zfile.clone(), &mut zfile) != 0) {
                if ((read_grid(&mut zhz[ihz as usize], &zfile) as u64) != ((nyhz * nxhz) as u64)) {
                    err!("Can not read in %s\n", zfile);
                }
                warn!("read in %s", zfile);
            } else {
                warn!("Using uniform cartesian grids for z-horizon%d", ihz);
                ixhz = 0;
                while (ixhz < nxhz) {
                    iyhz = 0;
                    while (iyhz < nyhz) {
                        zhz[ihz as usize][iyhz as usize][ixhz as usize] = (0 as f32);
                        iyhz += 1;
                    }
                    ixhz += 1;
                }
            }
            if (getparstring(&vfile.clone(), &mut vfile) != 0) {
                if ((read_grid(&mut v0hz[ihz as usize], &vfile) as u64) != ((nxhz * nyhz) as u64)) {
                    err!("Can not read in %s\n", vfile);
                }
                warn!("read in %s", vfile);
            } else {
                warn!("Using default for v-horizon %d", ihz);
                ixhz = 0;
                while (ixhz < nxhz) {
                    iyhz = 0;
                    while (iyhz < nyhz) {
                        v0hz[ihz as usize][iyhz as usize][ixhz as usize] = (0 as f32);
                        iyhz += 1;
                    }
                    ixhz += 1;
                }
            }
            if (getparstring(&dvdzfile.clone(), &mut dvdzfile) != 0) {
                if ((read_grid(&mut dvdzhz[ihz as usize], &dvdzfile) as u64)
                    != ((nyhz * nxhz) as u64))
                {
                    err!("Can not read in %s\n", dvdzfile);
                }
                warn!("read in %s", dvdzfile);
            } else {
                warn!("Using default for dvdz-horizon %d", ihz);
                ixhz = 0;
                while (ixhz < nxhz) {
                    iyhz = 0;
                    while (iyhz < nyhz) {
                        dvdzhz[ihz as usize][iyhz as usize][ixhz as usize] = (0 as f32);
                        iyhz += 1;
                    }
                    ixhz += 1;
                }
            }
            checkpars();
            ixhz = 0;
            while (ixhz < nxhz) {
                {
                    iyhz = 0;
                    while (iyhz < nyhz) {
                        {
                            alpha = ((xhz[ihz as usize][iyhz as usize][ixhz as usize] - xmin)
                                / (xmax - xmin));
                            fz0 = ((alpha * z01[ihz as usize])
                                + (((1 as f32) - alpha) * z00[ihz as usize]));
                            fz1 = ((alpha * z11[ihz as usize])
                                + (((1 as f32) - alpha) * z10[ihz as usize]));
                            beta = ((yhz[ihz as usize][iyhz as usize][ixhz as usize] - ymin)
                                / (ymax - ymin));
                            zflat[ihz as usize][iyhz as usize][ixhz as usize] =
                                ((beta * fz1) + (((1 as f32) - beta) * fz0));
                            zhz[ihz as usize][iyhz as usize][ixhz as usize] = zhz[ihz as usize]
                                [iyhz as usize][ixhz as usize]
                                + zflat[ihz as usize][iyhz as usize][ixhz as usize];
                        }
                        iyhz += 1;
                    }
                }
                ixhz += 1;
            }
        }
        ihz += 1;
    }
    // Set the model base beneath the maximum horizon depth.
    if nhz == 0 {
        c_segv();
    }
    iyhz = 0;
    while (iyhz < nyhz) {
        ixhz = 0;
        while (ixhz < nxhz) {
            zmax = (if ((zmax) > (zhz[(nhz - 1) as usize][iyhz as usize][ixhz as usize])) {
                (zmax)
            } else {
                (zhz[(nhz - 1) as usize][iyhz as usize][ixhz as usize])
            });
            ixhz += 1;
        }
        iyhz += 1;
    }
    zmax = zmax + blt;
    iyhz = 0;
    while (iyhz < nyhz) {
        ixhz = 0;
        while (ixhz < nxhz) {
            zhz[nhz as usize][iyhz as usize][ixhz as usize] = zmax;
            ixhz += 1;
        }
        iyhz += 1;
    }
    // iihz is the previous non-fictitious horizon.
    ihz = 0;
    while (ihz < nhz) {
        {
            if ((nficth != 0) && (ihz != 0)) {
                iiihz = ihz;
                while (iiihz >= 0) {
                    {
                        ifict = 0;
                        while (ifict < nficth) {
                            if (iiihz == ficth[ifict as usize]) {
                                break;
                            }
                            ifict += 1;
                        }
                        if (ifict == nficth) {
                            break;
                        }
                    }
                    iiihz -= 1;
                }
                iihz[ihz as usize] = iiihz;
                if (iiihz == nhz) {
                    iihz[ihz as usize] = ihz;
                }
            } else {
                iihz[ihz as usize] = ihz;
            }
            fprintf!(stderr(), "iihz[%d]=%d\n", ihz, iihz[ihz as usize]);
        }
        ihz += 1;
    }
    // Keep copies of the additive velocity and gradient grids.
    ihz = 0;
    while (ihz < nhz) {
        ixhz = 0;
        while (ixhz < nxhz) {
            iyhz = 0;
            while (iyhz < nyhz) {
                {
                    v0hztem[ihz as usize][iyhz as usize][ixhz as usize] =
                        v0hz[ihz as usize][iyhz as usize][ixhz as usize];
                    dvdzhztem[ihz as usize][iyhz as usize][ixhz as usize] =
                        dvdzhz[ihz as usize][iyhz as usize][ixhz as usize];
                }
                iyhz += 1;
            }
            ixhz += 1;
        }
        ihz += 1;
    }
    // Interpolate velocity and gradients, then construct horizon normals.
    ihz = 0;
    while (ihz < nhz) {
        {
            ixhz = 0;
            while (ixhz < nxhz) {
                {
                    iyhz = 0;
                    while (iyhz < nyhz) {
                        {
                            iiihz = iihz[ihz as usize];
                            if (ihz < (nhz - 1)) {
                                if (zhz[ihz as usize][iyhz as usize][ixhz as usize]
                                    > zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize])
                                {
                                    {
                                        let mut temp: f32 =
                                            (zhz[ihz as usize][iyhz as usize][ixhz as usize]);
                                        (zhz[ihz as usize][iyhz as usize][ixhz as usize]) =
                                            (zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize]);
                                        (zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize]) =
                                            temp;
                                    }

                                    iiihz = iihz[(ihz + 1) as usize];
                                }
                            }
                            // ficth can leave no non-fictitious predecessor,
                            // including after the depth swap above. C then reads
                            // index -1; match its usual SIGSEGV without a Rust
                            // panic. Occasional C models built from heap garbage
                            // on this path cannot be reproduced meaningfully.
                            if iiihz < 0 {
                                c_segv();
                            }
                            alpha = ((xhz[ihz as usize][iyhz as usize][ixhz as usize] - xmin)
                                / (xmax - xmin));
                            fv0 = ((alpha * v01[iiihz as usize])
                                + (((1 as f32) - alpha) * v00[iiihz as usize]));
                            fv1 = ((alpha * v11[iiihz as usize])
                                + (((1 as f32) - alpha) * v10[iiihz as usize]));
                            fdvdz0 = ((alpha * dvdz01[iiihz as usize])
                                + (((1 as f32) - alpha) * dvdz00[iiihz as usize]));
                            fdvdz1 = ((alpha * dvdz11[iiihz as usize])
                                + (((1 as f32) - alpha) * dvdz10[iiihz as usize]));
                            beta = ((yhz[ihz as usize][iyhz as usize][ixhz as usize] - ymin)
                                / (ymax - ymin));
                            vtemp = ((v0hztem[iiihz as usize][iyhz as usize][ixhz as usize]
                                + (beta * fv1))
                                + (((1 as f32) - beta) * fv0));
                            dvdzhz[ihz as usize][iyhz as usize][ixhz as usize] = ((dvdzhztem
                                [iiihz as usize][iyhz as usize][ixhz as usize]
                                + (beta * fdvdz1))
                                + (((1 as f32) - beta) * fdvdz0));
                            v0hz[ihz as usize][iyhz as usize][ixhz as usize] = (vtemp
                                + (dvdzhz[ihz as usize][iyhz as usize][ixhz as usize]
                                    * (zhz[ihz as usize][iyhz as usize][ixhz as usize]
                                        - zflat[iiihz as usize][iyhz as usize][ixhz as usize])));
                            v1hz[ihz as usize][iyhz as usize][ixhz as usize] = (vtemp
                                + (dvdzhz[ihz as usize][iyhz as usize][ixhz as usize]
                                    * (zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize]
                                        - zflat[iiihz as usize][iyhz as usize][ixhz as usize])));
                        }
                        iyhz += 1;
                    }
                }
                ixhz += 1;
            }
            ixhz = 0;
            while (ixhz < (nxhz - 1)) {
                {
                    iyhz = 0;
                    while (iyhz < (nyhz - 1)) {
                        {
                            let nn = tm_normal(
                                xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                            );
                            n1ll[iyhz as usize][ixhz as usize] = nn[0];
                            n2ll[iyhz as usize][ixhz as usize] = nn[1];
                            n3ll[iyhz as usize][ixhz as usize] = nn[2];
                            let nn = tm_normal(
                                xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                xhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                yhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                zhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                            );
                            n1ur[iyhz as usize][ixhz as usize] = nn[0];
                            n2ur[iyhz as usize][ixhz as usize] = nn[1];
                            n3ur[iyhz as usize][ixhz as usize] = nn[2];
                        }
                        iyhz += 1;
                    }
                }
                ixhz += 1;
            }
            n1hz[ihz as usize][0][0] = n1ll[0][0];
            n2hz[ihz as usize][0][0] = n2ll[0][0];
            n3hz[ihz as usize][0][0] = n3ll[0][0];
            n1hz[ihz as usize][0][(nxhz - 1) as usize] = n1ur[0][(nxhz - 2) as usize];
            n2hz[ihz as usize][0][(nxhz - 1) as usize] = n2ur[0][(nxhz - 2) as usize];
            n3hz[ihz as usize][0][(nxhz - 1) as usize] = n3ur[0][(nxhz - 2) as usize];
            n1hz[ihz as usize][(nyhz - 1) as usize][0] = n1ll[(nyhz - 2) as usize][0];
            n2hz[ihz as usize][(nyhz - 1) as usize][0] = n2ll[(nyhz - 2) as usize][0];
            n3hz[ihz as usize][(nyhz - 1) as usize][0] = n3ll[(nyhz - 2) as usize][0];
            n1hz[ihz as usize][(nyhz - 1) as usize][(nxhz - 1) as usize] =
                n1ur[(nyhz - 2) as usize][(nxhz - 2) as usize];
            n2hz[ihz as usize][(nyhz - 1) as usize][(nxhz - 1) as usize] =
                n2ur[(nyhz - 2) as usize][(nxhz - 2) as usize];
            n3hz[ihz as usize][(nyhz - 1) as usize][(nxhz - 1) as usize] =
                n3ur[(nyhz - 2) as usize][(nxhz - 2) as usize];
            ixhz = 1;
            while (ixhz < (nxhz - 1)) {
                {
                    n1hz[ihz as usize][0][ixhz as usize] = ((n1ll[0][(ixhz - 1) as usize]
                        + n1ur[0][(ixhz - 1) as usize])
                        + n1ll[0][ixhz as usize]);
                    n2hz[ihz as usize][0][ixhz as usize] = ((n2ll[0][(ixhz - 1) as usize]
                        + n2ur[0][(ixhz - 1) as usize])
                        + n2ll[0][ixhz as usize]);
                    n3hz[ihz as usize][0][ixhz as usize] = ((n3ll[0][(ixhz - 1) as usize]
                        + n3ur[0][(ixhz - 1) as usize])
                        + n1ll[0][ixhz as usize]);
                    n1hz[ihz as usize][(nyhz - 1) as usize][ixhz as usize] = ((n1ur
                        [(nyhz - 2) as usize][(ixhz - 1) as usize]
                        + n1ll[(nyhz - 2) as usize][ixhz as usize])
                        + n1ur[(nyhz - 2) as usize][ixhz as usize]);
                    n2hz[ihz as usize][(nyhz - 1) as usize][ixhz as usize] = ((n2ur
                        [(nyhz - 2) as usize][(ixhz - 1) as usize]
                        + n2ll[(nyhz - 2) as usize][ixhz as usize])
                        + n2ur[(nyhz - 2) as usize][ixhz as usize]);
                    n3hz[ihz as usize][(nyhz - 1) as usize][ixhz as usize] = ((n3ur
                        [(nyhz - 2) as usize][(ixhz - 1) as usize]
                        + n3ll[(nyhz - 2) as usize][ixhz as usize])
                        + n3ur[(nyhz - 2) as usize][ixhz as usize]);
                }
                ixhz += 1;
            }
            iyhz = 1;
            while (iyhz < (nyhz - 1)) {
                {
                    n1hz[ihz as usize][iyhz as usize][0] = ((n1ll[(iyhz - 1) as usize][0]
                        + n1ur[(iyhz - 1) as usize][0])
                        + n1ll[iyhz as usize][0]);
                    n2hz[ihz as usize][iyhz as usize][0] = ((n2ll[(iyhz - 1) as usize][0]
                        + n2ur[(iyhz - 1) as usize][0])
                        + n2ll[iyhz as usize][0]);
                    n3hz[ihz as usize][iyhz as usize][0] = ((n3ll[(iyhz - 1) as usize][0]
                        + n3ur[(iyhz - 1) as usize][0])
                        + n3ll[iyhz as usize][0]);
                    n1hz[ihz as usize][iyhz as usize][(nxhz - 1) as usize] = ((n1ur
                        [(iyhz - 1) as usize][(nxhz - 2) as usize]
                        + n1ll[iyhz as usize][(nxhz - 2) as usize])
                        + n1ur[iyhz as usize][(nxhz - 2) as usize]);
                    n2hz[ihz as usize][iyhz as usize][(nxhz - 1) as usize] = ((n2ur
                        [(iyhz - 1) as usize][(nxhz - 2) as usize]
                        + n2ll[iyhz as usize][(nxhz - 2) as usize])
                        + n2ur[iyhz as usize][(nxhz - 2) as usize]);
                    n3hz[ihz as usize][iyhz as usize][(nxhz - 1) as usize] = ((n3ur
                        [(iyhz - 1) as usize][(nxhz - 2) as usize]
                        + n3ll[iyhz as usize][(nxhz - 2) as usize])
                        + n3ur[iyhz as usize][(nxhz - 2) as usize]);
                }
                iyhz += 1;
            }
            iyhz = 1;
            while (iyhz < (nyhz - 1)) {
                {
                    ixhz = 1;
                    while (ixhz < (nxhz - 1)) {
                        {
                            n1hz[ihz as usize][iyhz as usize][ixhz as usize] = (((((n1ll
                                [iyhz as usize][(ixhz - 1) as usize]
                                + n1ur[iyhz as usize][(ixhz - 1) as usize])
                                + n1ll[iyhz as usize][ixhz as usize])
                                + n1ur[(iyhz - 1) as usize][ixhz as usize])
                                + n1ll[(iyhz - 1) as usize][ixhz as usize])
                                + n1ur[(iyhz - 1) as usize][(ixhz - 1) as usize]);
                            n2hz[ihz as usize][iyhz as usize][ixhz as usize] = (((((n2ll
                                [iyhz as usize][(ixhz - 1) as usize]
                                + n2ur[iyhz as usize][(ixhz - 1) as usize])
                                + n2ll[iyhz as usize][ixhz as usize])
                                + n2ur[(iyhz - 1) as usize][ixhz as usize])
                                + n2ll[(iyhz - 1) as usize][ixhz as usize])
                                + n2ur[(iyhz - 1) as usize][(ixhz - 1) as usize]);
                            n3hz[ihz as usize][iyhz as usize][ixhz as usize] = (((((n3ll
                                [iyhz as usize][(ixhz - 1) as usize]
                                + n3ur[iyhz as usize][(ixhz - 1) as usize])
                                + n3ll[iyhz as usize][ixhz as usize])
                                + n3ur[(iyhz - 1) as usize][ixhz as usize])
                                + n3ll[(iyhz - 1) as usize][ixhz as usize])
                                + n3ur[(iyhz - 1) as usize][(ixhz - 1) as usize]);
                        }
                        ixhz += 1;
                    }
                }
                iyhz += 1;
            }
            iyhz = 0;
            while (iyhz < nyhz) {
                ixhz = 0;
                while (ixhz < nxhz) {
                    let nn = tm_normalize([
                        n1hz[ihz as usize][iyhz as usize][ixhz as usize],
                        n2hz[ihz as usize][iyhz as usize][ixhz as usize],
                        n3hz[ihz as usize][iyhz as usize][ixhz as usize],
                    ]);
                    n1hz[ihz as usize][iyhz as usize][ixhz as usize] = nn[0];
                    n2hz[ihz as usize][iyhz as usize][ixhz as usize] = nn[1];
                    n3hz[ihz as usize][iyhz as usize][ixhz as usize] = nn[2];
                    ixhz += 1;
                }
                iyhz += 1;
            }
        }
        ihz += 1;
    }
    // Write the native binary model header to hzfile.
    (fwrite(&[nhz], &hzfp) as u64);
    (fwrite(&[nxhz], &hzfp) as u64);
    (fwrite(&[nyhz], &hzfp) as u64);
    (fwrite(&[xmin], &hzfp) as u64);
    (fwrite(&[xmax], &hzfp) as u64);
    (fwrite(&[ymin], &hzfp) as u64);
    (fwrite(&[ymax], &hzfp) as u64);
    (fwrite(&[zmin], &hzfp) as u64);
    (fwrite(&[zmax], &hzfp) as u64);
    // Write horizon coordinates and velocities to hzfile.
    ihz = 0;
    while (ihz < nhz) {
        {
            if ((write_grid(&xhz[ihz as usize], &hzfp) as u64) != ((nxhz * nyhz) as u64)) {
                err!("Can not write xhz to hzfp");
            }
            if ((write_grid(&yhz[ihz as usize], &hzfp) as u64) != ((nxhz * nyhz) as u64)) {
                err!("Can not write yhz to hzfp");
            }
            if ((write_grid(&zhz[ihz as usize], &hzfp) as u64) != ((nxhz * nyhz) as u64)) {
                err!("Can not write zhz to hzfp");
            }
            if ((write_grid(&v0hz[ihz as usize], &hzfp) as u64) != ((nxhz * nyhz) as u64)) {
                err!("Can not write v0hz to hzfp");
            }
            if ((write_grid(&v1hz[ihz as usize], &hzfp) as u64) != ((nxhz * nyhz) as u64)) {
                err!("Can not write v1hz to hzfp");
            }
        }
        ihz += 1;
    }
    // Count trapezoids, facets, and points for the connectivity tables.
    ntrap = ((nxhz - 1) * (nyhz - 1));
    nhztrap = (nhz * ntrap);
    nhztrap1 = ((nhz + 1) * ntrap);
    nfacetx = ((nyhz - 1) * nxhz);
    nhzfacetx = (nhz * nfacetx);
    nfacety = (nyhz * (nxhz - 1));
    nhzfacety = (nhz * nfacety);
    nhzpoint = ((nhz * nxhz) * nyhz);
    if (verbose != 0) {
        warn!("ntrap=%d nhztrap=%d nhztrap1=%d", ntrap, nhztrap, nhztrap1);
        warn!("nfacetx=%d nhzfacetx=%d", nfacetx, nhzfacetx);
        warn!("nfacety=%d nfacety=%d", nfacety, nhzfacety);
    }
    npoint = (2 * nhzpoint);
    point = vec![Point::default(); npoint as usize];
    if (verbose != 0) {
        fprintf!(stderr(), "npoint=%d\n", npoint);
    }
    nfacet = ((((2 * nhztrap1) + (6 * nhztrap)) + (2 * nhzfacetx)) + (2 * nhzfacety));
    facet = vec![Facet::default(); nfacet as usize];
    if (verbose != 0) {
        fprintf!(stderr(), "nfacet=%d\n", nfacet);
    }
    ntetra = (6 * nhztrap);
    tetra = vec![Tetra::default(); ntetra as usize];
    if (verbose != 0) {
        fprintf!(stderr(), "ntetra=%d\n", ntetra);
    }
    // Assign upper and lower control points for each layer.
    ihz = 0;
    while (ihz < nhz) {
        {
            fprintf!(stderr(), "processing: ihz=%d\n", ihz);
            ihz0 = (if ((nhz - 1) < (ihz + 1)) {
                (nhz - 1)
            } else {
                (ihz + 1)
            });
            iyhz = 0;
            while (iyhz < nyhz) {
                {
                    ixhz = 0;
                    while (ixhz < nxhz) {
                        {
                            ip = ((((ihz * nxhz) * nyhz) + (iyhz * nxhz)) + ixhz);
                            point[ip as usize].x[0] =
                                xhz[ihz as usize][iyhz as usize][ixhz as usize];
                            point[ip as usize].x[1] =
                                yhz[ihz as usize][iyhz as usize][ixhz as usize];
                            point[ip as usize].x[2] =
                                zhz[ihz as usize][iyhz as usize][ixhz as usize];
                            point[ip as usize].s = (((1f64
                                / (v0hz[ihz as usize][iyhz as usize][ixhz as usize] as f64))
                                / (v0hz[ihz as usize][iyhz as usize][ixhz as usize] as f64))
                                as f32);
                            point[ip as usize].n[0] =
                                n1hz[ihz as usize][iyhz as usize][ixhz as usize];
                            point[ip as usize].n[1] =
                                n2hz[ihz as usize][iyhz as usize][ixhz as usize];
                            point[ip as usize].n[2] =
                                n3hz[ihz as usize][iyhz as usize][ixhz as usize];
                            ip = ip + nhzpoint;
                            point[ip as usize].x[0] =
                                xhz[ihz0 as usize][iyhz as usize][ixhz as usize];
                            point[ip as usize].x[1] =
                                yhz[ihz0 as usize][iyhz as usize][ixhz as usize];
                            point[ip as usize].x[2] =
                                zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize];
                            point[ip as usize].s = (((1f64
                                / (v1hz[ihz as usize][iyhz as usize][ixhz as usize] as f64))
                                / (v1hz[ihz as usize][iyhz as usize][ixhz as usize] as f64))
                                as f32);
                            if (ihz == (nhz - 1)) {
                                point[ip as usize].n[0] = 0.0f32;
                                point[ip as usize].n[1] = 0.0f32;
                                point[ip as usize].n[2] = 1.0f32;
                            } else {
                                point[ip as usize].n[0] =
                                    n1hz[(ihz + 1) as usize][iyhz as usize][ixhz as usize];
                                point[ip as usize].n[1] =
                                    n2hz[(ihz + 1) as usize][iyhz as usize][ixhz as usize];
                                point[ip as usize].n[2] =
                                    n3hz[(ihz + 1) as usize][iyhz as usize][ixhz as usize];
                            }
                        }
                        ixhz += 1;
                    }
                }
                iyhz += 1;
            }
        }
        ihz += 1;
    }
    // For each trapezoid, use the original even/odd-row tetrahedron split.
    ihz = 0;
    while (ihz < nhz) {
        {
            ihz0 = (if ((nhz - 1) < (ihz + 1)) {
                (nhz - 1)
            } else {
                (ihz + 1)
            });
            iyhz = 0;
            while (iyhz < (nyhz - 1)) {
                {
                    ixhz = 0;
                    while (ixhz < (nxhz - 1)) {
                        {
                            if1 = (((ihz * ntrap) + (iyhz * (nxhz - 1))) + ixhz);
                            iwh = nhztrap1;
                            if2 = (iwh + if1);
                            iwh = (nhztrap1 * 2);
                            if3 = (iwh + if1);
                            iwh = ((nhztrap1 * 2) + nhztrap);
                            if4 = (iwh + if1);
                            iwh = ((nhztrap1 * 2) + (nhztrap * 2));
                            if5 = (iwh + if1);
                            iwh = ((nhztrap1 * 2) + (nhztrap * 3));
                            if6 = (iwh + if1);
                            iwh = ((nhztrap1 * 2) + (nhztrap * 4));
                            if7 = (iwh + if1);
                            iwh = ((nhztrap1 * 2) + (nhztrap * 5));
                            if8 = (iwh + if1);
                            iwh = ((nhztrap1 * 2) + (nhztrap * 6));
                            if9 = (((iwh + (ihz * nfacetx)) + (iyhz * nxhz)) + ixhz);
                            if10 = (nhzfacetx + if9);
                            iwh = (((nhztrap1 * 2) + (nhztrap * 6)) + (nhzfacetx * 2));
                            if11 = (((iwh + (ihz * nfacety)) + (iyhz * (nxhz - 1))) + ixhz);
                            if12 = (nhzfacety + if11);
                            it1 = (((ihz * ntrap) + (iyhz * (nxhz - 1))) + ixhz);
                            it2 = (nhztrap + it1);
                            it3 = ((nhztrap * 2) + it1);
                            it4 = ((nhztrap * 3) + it1);
                            it5 = ((nhztrap * 4) + it1);
                            it6 = ((nhztrap * 5) + it1);
                            i000 = ((((ihz * nxhz) * nyhz) + (iyhz * nxhz)) + ixhz);
                            i010 = (i000 + nxhz);
                            i100 = (i000 + 1);
                            i110 = ((i000 + nxhz) + 1);
                            i001 = (nhzpoint + i000);
                            i011 = (nhzpoint + i010);
                            i101 = (nhzpoint + i100);
                            i111 = (nhzpoint + i110);
                            tetra[it1 as usize].ireg = ihz;
                            tetra[it2 as usize].ireg = ihz;
                            tetra[it3 as usize].ireg = ihz;
                            tetra[it4 as usize].ireg = ihz;
                            tetra[it5 as usize].ireg = ihz;
                            tetra[it6 as usize].ireg = ihz;
                            if ((iyhz % 2) == 0) {
                                tetra[it1 as usize].ip[0] = i000;
                                tetra[it1 as usize].ip[1] = i010;
                                tetra[it1 as usize].ip[2] = i100;
                                tetra[it1 as usize].ip[3] = i101;
                                tetra[it2 as usize].ip[0] = i000;
                                tetra[it2 as usize].ip[1] = i010;
                                tetra[it2 as usize].ip[2] = i001;
                                tetra[it2 as usize].ip[3] = i101;
                                tetra[it3 as usize].ip[0] = i001;
                                tetra[it3 as usize].ip[1] = i011;
                                tetra[it3 as usize].ip[2] = i101;
                                tetra[it3 as usize].ip[3] = i010;
                                tetra[it4 as usize].ip[0] = i010;
                                tetra[it4 as usize].ip[1] = i110;
                                tetra[it4 as usize].ip[2] = i100;
                                tetra[it4 as usize].ip[3] = i101;
                                tetra[it5 as usize].ip[0] = i010;
                                tetra[it5 as usize].ip[1] = i011;
                                tetra[it5 as usize].ip[2] = i101;
                                tetra[it5 as usize].ip[3] = i110;
                                tetra[it6 as usize].ip[0] = i011;
                                tetra[it6 as usize].ip[1] = i101;
                                tetra[it6 as usize].ip[2] = i111;
                                tetra[it6 as usize].ip[3] = i110;
                                tetra[it1 as usize].ifacet[0] = if1;
                                tetra[it1 as usize].ifacet[1] = if3;
                                tetra[it1 as usize].ifacet[2] = if11;
                                tetra[it1 as usize].ifacet[3] = if5;
                                tetra[it2 as usize].ifacet[0] = if5;
                                tetra[it2 as usize].ifacet[1] = if9;
                                tetra[it2 as usize].ifacet[2] = if12;
                                tetra[it2 as usize].ifacet[3] = if7;
                                tetra[it3 as usize].ifacet[0] = if7;
                                tetra[it3 as usize].ifacet[1] = if10;
                                tetra[it3 as usize].ifacet[2] = if4;
                                tetra[it3 as usize].ifacet[3] = (if1 + ntrap);
                                tetra[it4 as usize].ifacet[0] = if2;
                                tetra[it4 as usize].ifacet[1] = if3;
                                tetra[it4 as usize].ifacet[2] = (if9 + 1);
                                tetra[it4 as usize].ifacet[3] = if6;
                                tetra[it5 as usize].ifacet[0] = if6;
                                tetra[it5 as usize].ifacet[1] = ((if11 + nxhz) - 1);
                                tetra[it5 as usize].ifacet[2] = if4;
                                tetra[it5 as usize].ifacet[3] = if8;
                                tetra[it6 as usize].ifacet[0] = if8;
                                tetra[it6 as usize].ifacet[1] = (if10 + 1);
                                tetra[it6 as usize].ifacet[2] = ((if12 + nxhz) - 1);
                                tetra[it6 as usize].ifacet[3] = (if2 + ntrap);
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                );
                                facet[if1 as usize].cn[0] = nn[0];
                                facet[if1 as usize].cn[1] = nn[1];
                                facet[if1 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                );
                                facet[if2 as usize].cn[0] = nn[0];
                                facet[if2 as usize].cn[1] = nn[1];
                                facet[if2 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                );
                                facet[if3 as usize].cn[0] = nn[0];
                                facet[if3 as usize].cn[1] = nn[1];
                                facet[if3 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                );
                                facet[if4 as usize].cn[0] = nn[0];
                                facet[if4 as usize].cn[1] = nn[1];
                                facet[if4 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                );
                                facet[if5 as usize].cn[0] = nn[0];
                                facet[if5 as usize].cn[1] = nn[1];
                                facet[if5 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                );
                                facet[if6 as usize].cn[0] = nn[0];
                                facet[if6 as usize].cn[1] = nn[1];
                                facet[if6 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                );
                                facet[if7 as usize].cn[0] = nn[0];
                                facet[if7 as usize].cn[1] = nn[1];
                                facet[if7 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                );
                                facet[if8 as usize].cn[0] = nn[0];
                                facet[if8 as usize].cn[1] = nn[1];
                                facet[if8 as usize].cn[2] = nn[2];
                                facet[if1 as usize].ip[0] = i000;
                                facet[if1 as usize].ip[1] = i100;
                                facet[if1 as usize].ip[2] = i010;
                                facet[if2 as usize].ip[0] = i010;
                                facet[if2 as usize].ip[1] = i110;
                                facet[if2 as usize].ip[2] = i100;
                                facet[if3 as usize].ip[0] = i010;
                                facet[if3 as usize].ip[1] = i100;
                                facet[if3 as usize].ip[2] = i101;
                                facet[if4 as usize].ip[0] = i010;
                                facet[if4 as usize].ip[1] = i011;
                                facet[if4 as usize].ip[2] = i101;
                                facet[if5 as usize].ip[0] = i000;
                                facet[if5 as usize].ip[1] = i010;
                                facet[if5 as usize].ip[2] = i101;
                                facet[if6 as usize].ip[0] = i010;
                                facet[if6 as usize].ip[1] = i110;
                                facet[if6 as usize].ip[2] = i101;
                                facet[if7 as usize].ip[0] = i001;
                                facet[if7 as usize].ip[1] = i010;
                                facet[if7 as usize].ip[2] = i101;
                                facet[if8 as usize].ip[0] = i011;
                                facet[if8 as usize].ip[1] = i110;
                                facet[if8 as usize].ip[2] = i101;
                            } else {
                                tetra[it1 as usize].ip[0] = i000;
                                tetra[it1 as usize].ip[1] = i010;
                                tetra[it1 as usize].ip[2] = i110;
                                tetra[it1 as usize].ip[3] = i111;
                                tetra[it2 as usize].ip[0] = i000;
                                tetra[it2 as usize].ip[1] = i010;
                                tetra[it2 as usize].ip[2] = i011;
                                tetra[it2 as usize].ip[3] = i111;
                                tetra[it3 as usize].ip[0] = i001;
                                tetra[it3 as usize].ip[1] = i011;
                                tetra[it3 as usize].ip[2] = i111;
                                tetra[it3 as usize].ip[3] = i000;
                                tetra[it4 as usize].ip[0] = i000;
                                tetra[it4 as usize].ip[1] = i100;
                                tetra[it4 as usize].ip[2] = i110;
                                tetra[it4 as usize].ip[3] = i111;
                                tetra[it5 as usize].ip[0] = i000;
                                tetra[it5 as usize].ip[1] = i001;
                                tetra[it5 as usize].ip[2] = i111;
                                tetra[it5 as usize].ip[3] = i100;
                                tetra[it6 as usize].ip[0] = i100;
                                tetra[it6 as usize].ip[1] = i001;
                                tetra[it6 as usize].ip[2] = i111;
                                tetra[it6 as usize].ip[3] = i101;
                                tetra[it1 as usize].ifacet[0] = if1;
                                tetra[it1 as usize].ifacet[1] = if3;
                                tetra[it1 as usize].ifacet[2] = ((if11 + nxhz) - 1);
                                tetra[it1 as usize].ifacet[3] = if5;
                                tetra[it2 as usize].ifacet[0] = if5;
                                tetra[it2 as usize].ifacet[1] = if9;
                                tetra[it2 as usize].ifacet[2] = ((if12 + nxhz) - 1);
                                tetra[it2 as usize].ifacet[3] = if7;
                                tetra[it3 as usize].ifacet[0] = if7;
                                tetra[it3 as usize].ifacet[1] = if10;
                                tetra[it3 as usize].ifacet[2] = if4;
                                tetra[it3 as usize].ifacet[3] = (if1 + ntrap);
                                tetra[it4 as usize].ifacet[0] = if2;
                                tetra[it4 as usize].ifacet[1] = if3;
                                tetra[it4 as usize].ifacet[2] = (if9 + 1);
                                tetra[it4 as usize].ifacet[3] = if6;
                                tetra[it5 as usize].ifacet[0] = if6;
                                tetra[it5 as usize].ifacet[1] = if4;
                                tetra[it5 as usize].ifacet[2] = if11;
                                tetra[it5 as usize].ifacet[3] = if8;
                                tetra[it6 as usize].ifacet[0] = if8;
                                tetra[it6 as usize].ifacet[1] = (if10 + 1);
                                tetra[it6 as usize].ifacet[2] = if12;
                                tetra[it6 as usize].ifacet[3] = (if2 + ntrap);
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                );
                                facet[if1 as usize].cn[0] = nn[0];
                                facet[if1 as usize].cn[1] = nn[1];
                                facet[if1 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                );
                                facet[if2 as usize].cn[0] = nn[0];
                                facet[if2 as usize].cn[1] = nn[1];
                                facet[if2 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize],
                                );
                                facet[if3 as usize].cn[0] = nn[0];
                                facet[if3 as usize].cn[1] = nn[1];
                                facet[if3 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize],
                                );
                                facet[if4 as usize].cn[0] = nn[0];
                                facet[if4 as usize].cn[1] = nn[1];
                                facet[if4 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize],
                                );
                                facet[if5 as usize].cn[0] = nn[0];
                                facet[if5 as usize].cn[1] = nn[1];
                                facet[if5 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize],
                                );
                                facet[if6 as usize].cn[0] = nn[0];
                                facet[if6 as usize].cn[1] = nn[1];
                                facet[if6 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize],
                                );
                                facet[if7 as usize].cn[0] = nn[0];
                                facet[if7 as usize].cn[1] = nn[1];
                                facet[if7 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize],
                                );
                                facet[if8 as usize].cn[0] = nn[0];
                                facet[if8 as usize].cn[1] = nn[1];
                                facet[if8 as usize].cn[2] = nn[2];
                                facet[if1 as usize].ip[0] = i000;
                                facet[if1 as usize].ip[1] = i010;
                                facet[if1 as usize].ip[2] = i110;
                                facet[if2 as usize].ip[0] = i000;
                                facet[if2 as usize].ip[1] = i100;
                                facet[if2 as usize].ip[2] = i110;
                                facet[if3 as usize].ip[0] = i000;
                                facet[if3 as usize].ip[1] = i110;
                                facet[if3 as usize].ip[2] = i111;
                                facet[if4 as usize].ip[0] = i000;
                                facet[if4 as usize].ip[1] = i001;
                                facet[if4 as usize].ip[2] = i111;
                                facet[if5 as usize].ip[0] = i000;
                                facet[if5 as usize].ip[1] = i010;
                                facet[if5 as usize].ip[2] = i111;
                                facet[if6 as usize].ip[0] = i000;
                                facet[if6 as usize].ip[1] = i100;
                                facet[if6 as usize].ip[2] = i111;
                                facet[if7 as usize].ip[0] = i000;
                                facet[if7 as usize].ip[1] = i011;
                                facet[if7 as usize].ip[2] = i111;
                                facet[if8 as usize].ip[0] = i100;
                                facet[if8 as usize].ip[1] = i001;
                                facet[if8 as usize].ip[2] = i111;
                            }
                            if (ihz == 0) {
                                facet[if1 as usize].itetra[0] = (-1);
                                facet[if2 as usize].itetra[0] = (-1);
                            } else {
                                facet[if1 as usize].itetra[0] =
                                    ((((nhztrap * 2) + ((ihz - 1) * ntrap)) + (iyhz * (nxhz - 1)))
                                        + ixhz);
                                facet[if2 as usize].itetra[0] =
                                    (facet[if1 as usize].itetra[0] + (3 * nhztrap));
                            }
                            facet[if1 as usize].itetra[1] = it1;
                            facet[if2 as usize].itetra[1] = it4;
                            facet[if3 as usize].itetra[0] = it1;
                            facet[if3 as usize].itetra[1] = it4;
                            facet[if4 as usize].itetra[0] = it3;
                            facet[if4 as usize].itetra[1] = it5;
                            facet[if5 as usize].itetra[0] = it1;
                            facet[if5 as usize].itetra[1] = it2;
                            facet[if6 as usize].itetra[0] = it4;
                            facet[if6 as usize].itetra[1] = it5;
                            facet[if7 as usize].itetra[0] = it2;
                            facet[if7 as usize].itetra[1] = it3;
                            facet[if8 as usize].itetra[0] = it5;
                            facet[if8 as usize].itetra[1] = it6;
                        }
                        ixhz += 1;
                    }
                }
                iyhz += 1;
            }
        }
        ihz += 1;
    }
    // Construct the facets facing x and y.
    ihz = 0;
    while (ihz < nhz) {
        {
            iyhz = 0;
            while (iyhz < (nyhz - 1)) {
                {
                    ixhz = 0;
                    while (ixhz < nxhz) {
                        {
                            iwh = ((nhztrap1 * 2) + (nhztrap * 6));
                            if9 = (((iwh + (ihz * nfacetx)) + (iyhz * nxhz)) + ixhz);
                            if10 = (nhzfacetx + if9);
                            it1 = (((ihz * ntrap) + (iyhz * (nxhz - 1))) + ixhz);
                            it2 = (nhztrap + it1);
                            it3 = ((nhztrap * 2) + it1);
                            it4 = ((nhztrap * 3) + it1);
                            it5 = ((nhztrap * 4) + it1);
                            it6 = ((nhztrap * 5) + it1);
                            i000 = ((((ihz * nxhz) * nyhz) + (iyhz * nxhz)) + ixhz);
                            i010 = (i000 + nxhz);
                            i100 = (i000 + 1);
                            i110 = ((i000 + nxhz) + 1);
                            i001 = (nhzpoint + i000);
                            i011 = (nhzpoint + i010);
                            i101 = (nhzpoint + i100);
                            i111 = (nhzpoint + i110);
                            if ((iyhz % 2) == 0) {
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                );
                                facet[if9 as usize].cn[0] = nn[0];
                                facet[if9 as usize].cn[1] = nn[1];
                                facet[if9 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                );
                                facet[if10 as usize].cn[0] = nn[0];
                                facet[if10 as usize].cn[1] = nn[1];
                                facet[if10 as usize].cn[2] = nn[2];
                                facet[if9 as usize].ip[0] = i000;
                                facet[if9 as usize].ip[1] = i010;
                                facet[if9 as usize].ip[2] = i001;
                                facet[if10 as usize].ip[0] = i010;
                                facet[if10 as usize].ip[1] = i001;
                                facet[if10 as usize].ip[2] = i011;
                            } else {
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize],
                                );
                                facet[if9 as usize].cn[0] = nn[0];
                                facet[if9 as usize].cn[1] = nn[1];
                                facet[if9 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][(iyhz + 1) as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][(iyhz + 1) as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                );
                                facet[if10 as usize].cn[0] = nn[0];
                                facet[if10 as usize].cn[1] = nn[1];
                                facet[if10 as usize].cn[2] = nn[2];
                                facet[if9 as usize].ip[0] = i000;
                                facet[if9 as usize].ip[1] = i010;
                                facet[if9 as usize].ip[2] = i011;
                                facet[if10 as usize].ip[0] = i000;
                                facet[if10 as usize].ip[1] = i001;
                                facet[if10 as usize].ip[2] = i011;
                            }
                            if (ixhz == 0) {
                                facet[if9 as usize].itetra[0] = (-1);
                                facet[if10 as usize].itetra[0] = (-1);
                            } else {
                                facet[if9 as usize].itetra[0] = (it4 - 1);
                                facet[if10 as usize].itetra[0] = (it6 - 1);
                            }
                            if (ixhz == (nxhz - 1)) {
                                facet[if9 as usize].itetra[1] = (-1);
                                facet[if10 as usize].itetra[1] = (-1);
                            } else {
                                facet[if9 as usize].itetra[1] = it2;
                                facet[if10 as usize].itetra[1] = it3;
                            }
                        }
                        ixhz += 1;
                    }
                }
                iyhz += 1;
            }
            iyhz = 0;
            while (iyhz < nyhz) {
                {
                    ixhz = 0;
                    while (ixhz < (nxhz - 1)) {
                        {
                            iwh = (((nhztrap1 * 2) + (nhztrap * 6)) + (nhzfacetx * 2));
                            if11 = (((iwh + (ihz * nfacety)) + (iyhz * (nxhz - 1))) + ixhz);
                            if12 = (nhzfacety + if11);
                            it1 = (((ihz * ntrap) + (iyhz * (nxhz - 1))) + ixhz);
                            it2 = (nhztrap + it1);
                            it3 = ((nhztrap * 2) + it1);
                            it4 = ((nhztrap * 3) + it1);
                            it5 = ((nhztrap * 4) + it1);
                            it6 = ((nhztrap * 5) + it1);
                            i000 = ((((ihz * nxhz) * nyhz) + (iyhz * nxhz)) + ixhz);
                            i010 = (i000 + nxhz);
                            i100 = (i000 + 1);
                            i110 = ((i000 + nxhz) + 1);
                            i001 = (nhzpoint + i000);
                            i011 = (nhzpoint + i010);
                            i101 = (nhzpoint + i100);
                            i111 = (nhzpoint + i110);
                            if ((iyhz % 2) == 0) {
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                );
                                facet[if11 as usize].cn[0] = nn[0];
                                facet[if11 as usize].cn[1] = nn[1];
                                facet[if11 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                );
                                facet[if12 as usize].cn[0] = nn[0];
                                facet[if12 as usize].cn[1] = nn[1];
                                facet[if12 as usize].cn[2] = nn[2];
                                facet[if11 as usize].ip[0] = i000;
                                facet[if11 as usize].ip[1] = i100;
                                facet[if11 as usize].ip[2] = i101;
                                facet[if12 as usize].ip[0] = i000;
                                facet[if12 as usize].ip[1] = i001;
                                facet[if12 as usize].ip[2] = i101;
                                if (iyhz == 0) {
                                    facet[if11 as usize].itetra[0] = (-1);
                                    facet[if12 as usize].itetra[0] = (-1);
                                } else {
                                    facet[if11 as usize].itetra[0] = (it1 - (nxhz - 1));
                                    facet[if12 as usize].itetra[0] = (it2 - (nxhz - 1));
                                }
                                if (iyhz == (nyhz - 1)) {
                                    facet[if11 as usize].itetra[1] = (-1);
                                    facet[if12 as usize].itetra[1] = (-1);
                                } else {
                                    facet[if11 as usize].itetra[1] = it1;
                                    facet[if12 as usize].itetra[1] = it2;
                                }
                            } else {
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[ihz as usize][iyhz as usize][ixhz as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                );
                                facet[if11 as usize].cn[0] = nn[0];
                                facet[if11 as usize].cn[1] = nn[1];
                                facet[if11 as usize].cn[2] = nn[2];
                                let nn = tm_normal(
                                    xhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    xhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    yhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz0 as usize][iyhz as usize][(ixhz + 1) as usize],
                                    yhz[ihz0 as usize][iyhz as usize][ixhz as usize],
                                    zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][(ixhz + 1) as usize],
                                    zhz[(ihz + 1) as usize][iyhz as usize][ixhz as usize],
                                );
                                facet[if12 as usize].cn[0] = nn[0];
                                facet[if12 as usize].cn[1] = nn[1];
                                facet[if12 as usize].cn[2] = nn[2];
                                facet[if11 as usize].ip[0] = i000;
                                facet[if11 as usize].ip[1] = i100;
                                facet[if11 as usize].ip[2] = i001;
                                facet[if12 as usize].ip[0] = i100;
                                facet[if12 as usize].ip[1] = i001;
                                facet[if12 as usize].ip[2] = i101;
                                facet[if11 as usize].itetra[0] = (it5 - (nxhz - 1));
                                facet[if12 as usize].itetra[0] = (it6 - (nxhz - 1));
                                if (iyhz == (nyhz - 1)) {
                                    facet[if11 as usize].itetra[1] = (-1);
                                    facet[if12 as usize].itetra[1] = (-1);
                                } else {
                                    facet[if11 as usize].itetra[1] = it5;
                                    facet[if12 as usize].itetra[1] = it6;
                                }
                            }
                        }
                        ixhz += 1;
                    }
                }
                iyhz += 1;
            }
        }
        ihz += 1;
    }
    // Construct the bottom facets.
    iyhz = 0;
    while (iyhz < (nyhz - 1)) {
        {
            ixhz = 0;
            while (ixhz < (nxhz - 1)) {
                {
                    if1 = (((nhz * ntrap) + (iyhz * (nxhz - 1))) + ixhz);
                    if2 = (nhztrap1 + if1);
                    i000 = (((nhzpoint + (((nhz - 1) * nxhz) * nyhz)) + (iyhz * nxhz)) + ixhz);
                    i010 = (i000 + nxhz);
                    i100 = (i000 + 1);
                    i110 = ((i000 + nxhz) + 1);
                    if ((iyhz % 2) == 0) {
                        facet[if1 as usize].ip[0] = i000;
                        facet[if1 as usize].ip[1] = i010;
                        facet[if1 as usize].ip[2] = i100;
                        facet[if2 as usize].ip[0] = i010;
                        facet[if2 as usize].ip[1] = i100;
                        facet[if2 as usize].ip[2] = i110;
                    } else {
                        facet[if1 as usize].ip[0] = i000;
                        facet[if1 as usize].ip[1] = i010;
                        facet[if1 as usize].ip[2] = i110;
                        facet[if2 as usize].ip[0] = i000;
                        facet[if2 as usize].ip[1] = i100;
                        facet[if2 as usize].ip[2] = i110;
                    }
                    facet[if1 as usize].cn[0] = 0.0f32;
                    facet[if1 as usize].cn[1] = 0.0f32;
                    facet[if1 as usize].cn[2] = ((-1f64) as f32);
                    facet[if2 as usize].cn[0] = 0.0f32;
                    facet[if2 as usize].cn[1] = 0.0f32;
                    facet[if2 as usize].cn[2] = ((-1f64) as f32);
                    it3 = ((((nhztrap * 2) + ((nhz - 1) * ntrap)) + (iyhz * (nxhz - 1))) + ixhz);
                    it6 = ((nhztrap * 3) + it3);
                    facet[if1 as usize].itetra[0] = it3;
                    facet[if1 as usize].itetra[1] = (-2);
                    facet[if2 as usize].itetra[0] = it6;
                    facet[if2 as usize].itetra[1] = (-2);
                }
                ixhz += 1;
            }
        }
        iyhz += 1;
    }
    // Calculate sloth gradients from the four control points.
    itetra = 0;
    while (itetra < ntetra) {
        {
            ip0 = tetra[itetra as usize].ip[0];
            ip1 = tetra[itetra as usize].ip[1];
            ip2 = tetra[itetra as usize].ip[2];
            ip3 = tetra[itetra as usize].ip[3];
            tetra[itetra as usize].gs = gradient(
                point[ip0 as usize].x,
                point[ip0 as usize].s,
                point[ip1 as usize].x,
                point[ip1 as usize].s,
                point[ip2 as usize].x,
                point[ip2 as usize].s,
                point[ip3 as usize].x,
                point[ip3 as usize].s,
            );
        }
        itetra += 1;
    }
    fprintf!(stdout(), "%f=xmin\n", (xmin as f64));
    fprintf!(stdout(), "%f=xmax\n", (xmax as f64));
    fprintf!(stdout(), "%f=ymin\n", (ymin as f64));
    fprintf!(stdout(), "%f=ymax\n", (ymax as f64));
    fprintf!(stdout(), "%f=zmin\n", (zmin as f64));
    fprintf!(stdout(), "%f=zmax\n", (zmax as f64));
    write_point(&mut point, npoint);
    // Check whether layers are smeared out; reconnect coincident facets.
    ixhz = 0;
    while (ixhz < (nxhz - 1)) {
        {
            iyhz = 0;
            while (iyhz < (nyhz - 1)) {
                {
                    ihz = 0;
                    while (ihz < nhz) {
                        {
                            iiihz = (ihz + 1);
                            while (iiihz <= nhz) {
                                {
                                    dis0 = (((zhz[ihz as usize][iyhz as usize][ixhz as usize]
                                        - zhz[iiihz as usize][iyhz as usize][ixhz as usize])
                                        as f64)
                                        .abs() as f32);
                                    dis1 = (((zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize]
                                        - zhz[iiihz as usize][(iyhz + 1) as usize][ixhz as usize])
                                        as f64)
                                        .abs() as f32);
                                    dis2 = (((zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize]
                                        - zhz[iiihz as usize][iyhz as usize][(ixhz + 1) as usize])
                                        as f64)
                                        .abs() as f32);
                                    if ((((dis0 as f64) > 1.0E-4f64)
                                        || ((dis1 as f64) > 1.0E-4f64))
                                        || ((dis2 as f64) > 1.0E-4f64))
                                    {
                                        break;
                                    }
                                }
                                iiihz += 1;
                            }
                            jumphz = ((iiihz - 1) - ihz);
                            if (jumphz == 0) {
                                ihz += 1;
                                continue;
                            }
                            if ((((dis0 as f64) < 1.0E-4f64) && ((dis1 as f64) < 1.0E-4f64))
                                && ((dis2 as f64) < 1.0E-4f64))
                            {
                                if1 = (((ihz * ntrap) + (iyhz * (nxhz - 1))) + ixhz);
                                if1p = (if1 + (((iiihz - 1) - ihz) * ntrap));
                                facet[if1 as usize].itetra[1] = facet[if1p as usize].itetra[1];
                                facet[if1p as usize].itetra[0] = facet[if1 as usize].itetra[0];
                            }
                        }
                        ihz += 1;
                    }
                }
                iyhz = iyhz + 2;
            }
        }
        ixhz += 1;
    }
    ixhz = 0;
    while (ixhz < (nxhz - 1)) {
        {
            iyhz = 0;
            while (iyhz < (nyhz - 1)) {
                {
                    ihz = 0;
                    while (ihz < nhz) {
                        {
                            iiihz = (ihz + 1);
                            while (iiihz <= nhz) {
                                {
                                    dis0 = (((zhz[ihz as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize]
                                        - zhz[iiihz as usize][(iyhz + 1) as usize]
                                            [(ixhz + 1) as usize])
                                        as f64)
                                        .abs() as f32);
                                    dis1 = (((zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize]
                                        - zhz[iiihz as usize][(iyhz + 1) as usize][ixhz as usize])
                                        as f64)
                                        .abs() as f32);
                                    dis2 = (((zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize]
                                        - zhz[iiihz as usize][iyhz as usize][(ixhz + 1) as usize])
                                        as f64)
                                        .abs() as f32);
                                    if ((((dis0 as f64) > 1.0E-4f64)
                                        || ((dis1 as f64) > 1.0E-4f64))
                                        || ((dis2 as f64) > 1.0E-4f64))
                                    {
                                        break;
                                    }
                                }
                                iiihz += 1;
                            }
                            jumphz = ((iiihz - 1) - ihz);
                            if (jumphz == 0) {
                                ihz += 1;
                                continue;
                            }
                            if ((((dis0 as f64) < 1.0E-4f64) && ((dis1 as f64) < 1.0E-4f64))
                                && ((dis2 as f64) < 1.0E-4f64))
                            {
                                if2 = (((nhztrap1 + (ihz * ntrap)) + (iyhz * (nxhz - 1))) + ixhz);
                                if2p = (if2 + (jumphz * ntrap));
                                facet[if2 as usize].itetra[1] = facet[if2p as usize].itetra[1];
                                facet[if2p as usize].itetra[0] = facet[if2 as usize].itetra[0];
                            }
                        }
                        ihz += 1;
                    }
                }
                iyhz = iyhz + 2;
            }
        }
        ixhz += 1;
    }
    ixhz = 0;
    while (ixhz < (nxhz - 1)) {
        {
            iyhz = 1;
            while (iyhz < (nyhz - 1)) {
                {
                    ihz = 0;
                    while (ihz < nhz) {
                        {
                            iiihz = (ihz + 1);
                            while (iiihz <= nhz) {
                                {
                                    dis0 = (((zhz[ihz as usize][iyhz as usize][ixhz as usize]
                                        - zhz[iiihz as usize][iyhz as usize][ixhz as usize])
                                        as f64)
                                        .abs() as f32);
                                    dis1 = (((zhz[ihz as usize][(iyhz + 1) as usize][ixhz as usize]
                                        - zhz[iiihz as usize][(iyhz + 1) as usize][ixhz as usize])
                                        as f64)
                                        .abs() as f32);
                                    dis2 = (((zhz[ihz as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize]
                                        - zhz[iiihz as usize][(iyhz + 1) as usize]
                                            [(ixhz + 1) as usize])
                                        as f64)
                                        .abs() as f32);
                                    if ((((dis0 as f64) > 1.0E-4f64)
                                        || ((dis1 as f64) > 1.0E-4f64))
                                        || ((dis2 as f64) > 1.0E-4f64))
                                    {
                                        break;
                                    }
                                }
                                iiihz += 1;
                            }
                            jumphz = ((iiihz - 1) - ihz);
                            if (jumphz == 0) {
                                ihz += 1;
                                continue;
                            }
                            if ((((dis0 as f64) < 1.0E-4f64) && ((dis1 as f64) < 1.0E-4f64))
                                && ((dis2 as f64) < 1.0E-4f64))
                            {
                                if1 = (((ihz * ntrap) + (iyhz * (nxhz - 1))) + ixhz);
                                if1p = (if1 + (((iiihz - 1) - ihz) * ntrap));
                                facet[if1 as usize].itetra[1] = facet[if1p as usize].itetra[1];
                                facet[if1p as usize].itetra[0] = facet[if1 as usize].itetra[0];
                            }
                        }
                        ihz += 1;
                    }
                }
                iyhz = iyhz + 2;
            }
        }
        ixhz += 1;
    }
    ixhz = 0;
    while (ixhz < (nxhz - 1)) {
        {
            iyhz = 1;
            while (iyhz < (nyhz - 1)) {
                {
                    ihz = 0;
                    while (ihz < nhz) {
                        {
                            iiihz = (ihz + 1);
                            while (iiihz <= nhz) {
                                {
                                    dis0 = (((zhz[ihz as usize][(iyhz + 1) as usize]
                                        [(ixhz + 1) as usize]
                                        - zhz[iiihz as usize][(iyhz + 1) as usize]
                                            [(ixhz + 1) as usize])
                                        as f64)
                                        .abs() as f32);
                                    dis1 = (((zhz[ihz as usize][iyhz as usize][ixhz as usize]
                                        - zhz[iiihz as usize][iyhz as usize][ixhz as usize])
                                        as f64)
                                        .abs() as f32);
                                    dis2 = (((zhz[ihz as usize][iyhz as usize][(ixhz + 1) as usize]
                                        - zhz[iiihz as usize][iyhz as usize][(ixhz + 1) as usize])
                                        as f64)
                                        .abs() as f32);
                                    if ((((dis0 as f64) > 1.0E-4f64)
                                        || ((dis1 as f64) > 1.0E-4f64))
                                        || ((dis2 as f64) > 1.0E-4f64))
                                    {
                                        break;
                                    }
                                }
                                iiihz += 1;
                            }
                            jumphz = ((iiihz - 1) - ihz);
                            if (jumphz == 0) {
                                ihz += 1;
                                continue;
                            }
                            if ((((dis0 as f64) < 1.0E-4f64) && ((dis1 as f64) < 1.0E-4f64))
                                && ((dis2 as f64) < 1.0E-4f64))
                            {
                                if2 = (((nhztrap1 + (ihz * ntrap)) + (iyhz * (nxhz - 1))) + ixhz);
                                if2p = (if2 + (((iiihz - 1) - ihz) * ntrap));
                                facet[if2 as usize].itetra[1] = facet[if2p as usize].itetra[1];
                                facet[if2p as usize].itetra[0] = facet[if2 as usize].itetra[0];
                            }
                        }
                        ihz += 1;
                    }
                }
                iyhz = iyhz + 2;
            }
        }
        ixhz += 1;
    }
    // Calculate facet areas with the tetra library.
    ifs = 0;
    while (ifs < nfacet) {
        facet[ifs as usize].area = area3d(
            &point[(facet[ifs as usize].ip[0]) as usize].x,
            &point[(facet[ifs as usize].ip[1]) as usize].x,
            &point[(facet[ifs as usize].ip[2]) as usize].x,
        );
        ifs += 1;
    }
    // Orient facet normals consistently.
    ifs = 0;
    while (ifs < nfacet) {
        {
            if ((facet[ifs as usize].cn[2] as f64) > 0f64) {
                facet[ifs as usize].cn[0] = (-facet[ifs as usize].cn[0]);
                facet[ifs as usize].cn[1] = (-facet[ifs as usize].cn[1]);
                facet[ifs as usize].cn[2] = (-facet[ifs as usize].cn[2]);
            }
            if ((facet[ifs as usize].cn[2] as f64) == 0f64) {
                if ((facet[ifs as usize].cn[0] as f64) < 0f64) {
                    facet[ifs as usize].cn[0] = (-facet[ifs as usize].cn[0]);
                    facet[ifs as usize].cn[1] = (-facet[ifs as usize].cn[1]);
                    facet[ifs as usize].cn[2] = (-facet[ifs as usize].cn[2]);
                } else {
                    if ((facet[ifs as usize].cn[0] as f64) == 0f64) {
                        if ((facet[ifs as usize].cn[1] as f64) < 0f64) {
                            facet[ifs as usize].cn[0] = (-facet[ifs as usize].cn[0]);
                            facet[ifs as usize].cn[1] = (-facet[ifs as usize].cn[1]);
                            facet[ifs as usize].cn[2] = (-facet[ifs as usize].cn[2]);
                        }
                    }
                }
            }
        }
        ifs += 1;
    }
    // Initialize planar coefficients, then fit curved tiles between regions.
    ifs = 0;
    while (ifs < nfacet) {
        {
            facet[ifs as usize].ct[0] = 0.0f32;
            facet[ifs as usize].ct[1] = 0.0f32;
            facet[ifs as usize].ct[2] = 0.0f32;
            facet[ifs as usize].ct[3] = 0.0f32;
            facet[ifs as usize].ct[4] = 0.0f32;
            facet[ifs as usize].ct[5] = 0.0f32;
            facet[ifs as usize].ct[6] = facet[ifs as usize].cn[0];
            facet[ifs as usize].ct[7] = facet[ifs as usize].cn[1];
            facet[ifs as usize].ct[8] = facet[ifs as usize].cn[2];
            facet[ifs as usize].ct[9] = ((((-point[(facet[ifs as usize].ip[0]) as usize].x[0])
                * facet[ifs as usize].cn[0])
                - (point[(facet[ifs as usize].ip[0]) as usize].x[1] * facet[ifs as usize].cn[1]))
                - (point[(facet[ifs as usize].ip[0]) as usize].x[2] * facet[ifs as usize].cn[2]));
            if ((region(&tetra, &point, &facet, facet[ifs as usize].itetra[0])
                == region(&tetra, &point, &facet, facet[ifs as usize].itetra[1]))
                || ((facet[ifs as usize].cn[2] as f64).abs() < 0.02f64))
            {
                ifs += 1;
                continue;
            }
            a[0][0] = (((2f64 * (point[(facet[ifs as usize].ip[0]) as usize].x[0] as f64))
                * (point[(facet[ifs as usize].ip[0]) as usize].n[2] as f64))
                as f32);
            a[0][1] = (point[(facet[ifs as usize].ip[0]) as usize].x[1]
                * point[(facet[ifs as usize].ip[0]) as usize].n[2]);
            a[0][2] = (point[(facet[ifs as usize].ip[0]) as usize].x[0]
                * (point[(facet[ifs as usize].ip[0]) as usize].n[2]
                    - point[(facet[ifs as usize].ip[0]) as usize].n[0]));
            a[0][3] = ((-point[(facet[ifs as usize].ip[0]) as usize].x[1])
                * point[(facet[ifs as usize].ip[0]) as usize].n[0]);
            a[0][4] = 0.0f32;
            a[0][5] = ((((-2f64) * (point[(facet[ifs as usize].ip[0]) as usize].x[2] as f64))
                * (point[(facet[ifs as usize].ip[0]) as usize].n[0] as f64))
                as f32);
            a[0][6] = point[(facet[ifs as usize].ip[0]) as usize].n[2];
            a[0][7] = 0.0f32;
            a[0][8] = 0.0f32;
            a[1][0] = 0.0f32;
            a[1][1] = (point[(facet[ifs as usize].ip[0]) as usize].x[0]
                * point[(facet[ifs as usize].ip[0]) as usize].n[2]);
            a[1][2] = ((-point[(facet[ifs as usize].ip[0]) as usize].x[0])
                * point[(facet[ifs as usize].ip[0]) as usize].n[1]);
            a[1][3] = ((point[(facet[ifs as usize].ip[0]) as usize].x[2]
                * point[(facet[ifs as usize].ip[0]) as usize].n[2])
                - (point[(facet[ifs as usize].ip[0]) as usize].x[1]
                    * point[(facet[ifs as usize].ip[0]) as usize].n[1]));
            a[1][4] = (((2f64 * (point[(facet[ifs as usize].ip[0]) as usize].x[1] as f64))
                * (point[(facet[ifs as usize].ip[0]) as usize].n[2] as f64))
                as f32);
            a[1][5] = ((((-2f64) * (point[(facet[ifs as usize].ip[0]) as usize].x[2] as f64))
                * (point[(facet[ifs as usize].ip[0]) as usize].n[1] as f64))
                as f32);
            a[1][6] = 0.0f32;
            a[1][7] = point[(facet[ifs as usize].ip[0]) as usize].n[2];
            a[1][8] = 0.0f32;
            a[2][0] = (((2f64 * (point[(facet[ifs as usize].ip[1]) as usize].x[0] as f64))
                * (point[(facet[ifs as usize].ip[1]) as usize].n[2] as f64))
                as f32);
            a[2][1] = (point[(facet[ifs as usize].ip[1]) as usize].x[1]
                * point[(facet[ifs as usize].ip[1]) as usize].n[2]);
            a[2][2] = (point[(facet[ifs as usize].ip[1]) as usize].x[0]
                * (point[(facet[ifs as usize].ip[1]) as usize].n[2]
                    - point[(facet[ifs as usize].ip[1]) as usize].n[0]));
            a[2][3] = ((-point[(facet[ifs as usize].ip[1]) as usize].x[1])
                * point[(facet[ifs as usize].ip[1]) as usize].n[0]);
            a[2][4] = 0.0f32;
            a[2][5] = ((((-2f64) * (point[(facet[ifs as usize].ip[1]) as usize].x[2] as f64))
                * (point[(facet[ifs as usize].ip[1]) as usize].n[0] as f64))
                as f32);
            a[2][6] = point[(facet[ifs as usize].ip[1]) as usize].n[2];
            // Preserve the repeated a[2][6] assignment in C; a[2][7] is untouched.
            a[2][6] = 0.0f32;
            a[2][8] = 0.0f32;
            a[3][0] = 0.0f32;
            a[3][1] = (point[(facet[ifs as usize].ip[1]) as usize].x[0]
                * point[(facet[ifs as usize].ip[1]) as usize].n[2]);
            a[3][2] = ((-point[(facet[ifs as usize].ip[1]) as usize].x[0])
                * point[(facet[ifs as usize].ip[1]) as usize].n[1]);
            a[3][3] = ((point[(facet[ifs as usize].ip[1]) as usize].x[2]
                * point[(facet[ifs as usize].ip[1]) as usize].n[2])
                - (point[(facet[ifs as usize].ip[1]) as usize].x[1]
                    * point[(facet[ifs as usize].ip[1]) as usize].n[1]));
            a[3][4] = (((2f64 * (point[(facet[ifs as usize].ip[1]) as usize].x[1] as f64))
                * (point[(facet[ifs as usize].ip[1]) as usize].n[2] as f64))
                as f32);
            a[3][5] = ((((-2f64) * (point[(facet[ifs as usize].ip[1]) as usize].x[2] as f64))
                * (point[(facet[ifs as usize].ip[1]) as usize].n[1] as f64))
                as f32);
            a[3][6] = 0.0f32;
            a[3][7] = point[(facet[ifs as usize].ip[1]) as usize].n[2];
            a[3][8] = 0.0f32;
            a[4][0] = (((2f64 * (point[(facet[ifs as usize].ip[2]) as usize].x[0] as f64))
                * (point[(facet[ifs as usize].ip[2]) as usize].n[2] as f64))
                as f32);
            a[4][1] = (point[(facet[ifs as usize].ip[2]) as usize].x[1]
                * point[(facet[ifs as usize].ip[2]) as usize].n[2]);
            a[4][2] = (point[(facet[ifs as usize].ip[2]) as usize].x[0]
                * (point[(facet[ifs as usize].ip[2]) as usize].n[2]
                    - point[(facet[ifs as usize].ip[2]) as usize].n[0]));
            a[4][3] = ((-point[(facet[ifs as usize].ip[2]) as usize].x[1])
                * point[(facet[ifs as usize].ip[2]) as usize].n[0]);
            a[4][4] = 0.0f32;
            a[4][5] = ((((-2f64) * (point[(facet[ifs as usize].ip[2]) as usize].x[2] as f64))
                * (point[(facet[ifs as usize].ip[2]) as usize].n[0] as f64))
                as f32);
            a[4][6] = point[(facet[ifs as usize].ip[2]) as usize].n[2];
            a[4][7] = 0.0f32;
            a[4][8] = 0.0f32;
            a[5][0] = 0.0f32;
            a[5][1] = (point[(facet[ifs as usize].ip[2]) as usize].x[0]
                * point[(facet[ifs as usize].ip[2]) as usize].n[2]);
            a[5][2] = ((-point[(facet[ifs as usize].ip[2]) as usize].x[2])
                * point[(facet[ifs as usize].ip[2]) as usize].n[1]);
            a[5][3] = ((point[(facet[ifs as usize].ip[2]) as usize].x[2]
                * point[(facet[ifs as usize].ip[2]) as usize].n[2])
                - (point[(facet[ifs as usize].ip[2]) as usize].x[1]
                    * point[(facet[ifs as usize].ip[2]) as usize].n[1]));
            a[5][4] = (((2f64 * (point[(facet[ifs as usize].ip[2]) as usize].x[1] as f64))
                * (point[(facet[ifs as usize].ip[2]) as usize].n[2] as f64))
                as f32);
            a[5][5] = ((((-2f64) * (point[(facet[ifs as usize].ip[2]) as usize].x[2] as f64))
                * (point[(facet[ifs as usize].ip[2]) as usize].n[1] as f64))
                as f32);
            a[5][6] = 0.0f32;
            a[5][7] = point[(facet[ifs as usize].ip[2]) as usize].n[2];
            a[5][8] = 0.0f32;
            a[6][0] = (point[(facet[ifs as usize].ip[0]) as usize].x[0]
                * point[(facet[ifs as usize].ip[0]) as usize].x[0]);
            a[6][1] = (point[(facet[ifs as usize].ip[0]) as usize].x[0]
                * point[(facet[ifs as usize].ip[0]) as usize].x[1]);
            a[6][2] = (point[(facet[ifs as usize].ip[0]) as usize].x[0]
                * point[(facet[ifs as usize].ip[0]) as usize].x[2]);
            a[6][3] = (point[(facet[ifs as usize].ip[0]) as usize].x[1]
                * point[(facet[ifs as usize].ip[0]) as usize].x[2]);
            a[6][4] = (point[(facet[ifs as usize].ip[0]) as usize].x[1]
                * point[(facet[ifs as usize].ip[0]) as usize].x[1]);
            a[6][5] = (point[(facet[ifs as usize].ip[0]) as usize].x[2]
                * point[(facet[ifs as usize].ip[0]) as usize].x[2]);
            a[6][6] = point[(facet[ifs as usize].ip[0]) as usize].x[0];
            a[6][7] = point[(facet[ifs as usize].ip[0]) as usize].x[1];
            a[6][8] = 1.0f32;
            a[7][0] = (point[(facet[ifs as usize].ip[1]) as usize].x[0]
                * point[(facet[ifs as usize].ip[1]) as usize].x[0]);
            a[7][1] = (point[(facet[ifs as usize].ip[1]) as usize].x[0]
                * point[(facet[ifs as usize].ip[1]) as usize].x[1]);
            a[7][2] = (point[(facet[ifs as usize].ip[1]) as usize].x[0]
                * point[(facet[ifs as usize].ip[1]) as usize].x[2]);
            a[7][3] = (point[(facet[ifs as usize].ip[1]) as usize].x[1]
                * point[(facet[ifs as usize].ip[1]) as usize].x[2]);
            a[7][4] = (point[(facet[ifs as usize].ip[1]) as usize].x[1]
                * point[(facet[ifs as usize].ip[1]) as usize].x[1]);
            a[7][5] = (point[(facet[ifs as usize].ip[1]) as usize].x[2]
                * point[(facet[ifs as usize].ip[1]) as usize].x[2]);
            a[7][6] = point[(facet[ifs as usize].ip[1]) as usize].x[0];
            a[7][7] = point[(facet[ifs as usize].ip[1]) as usize].x[1];
            a[7][8] = 1.0f32;
            a[8][0] = (point[(facet[ifs as usize].ip[2]) as usize].x[0]
                * point[(facet[ifs as usize].ip[2]) as usize].x[0]);
            a[8][1] = (point[(facet[ifs as usize].ip[2]) as usize].x[0]
                * point[(facet[ifs as usize].ip[2]) as usize].x[1]);
            a[8][2] = (point[(facet[ifs as usize].ip[2]) as usize].x[0]
                * point[(facet[ifs as usize].ip[2]) as usize].x[2]);
            a[8][3] = (point[(facet[ifs as usize].ip[2]) as usize].x[1]
                * point[(facet[ifs as usize].ip[2]) as usize].x[2]);
            a[8][4] = (point[(facet[ifs as usize].ip[2]) as usize].x[1]
                * point[(facet[ifs as usize].ip[2]) as usize].x[1]);
            a[8][5] = (point[(facet[ifs as usize].ip[2]) as usize].x[2]
                * point[(facet[ifs as usize].ip[2]) as usize].x[2]);
            a[8][6] = point[(facet[ifs as usize].ip[2]) as usize].x[0];
            a[8][7] = point[(facet[ifs as usize].ip[2]) as usize].x[1];
            a[8][8] = 1.0f32;
            b[0] = point[(facet[ifs as usize].ip[0]) as usize].n[0];
            b[1] = point[(facet[ifs as usize].ip[0]) as usize].n[1];
            b[2] = point[(facet[ifs as usize].ip[1]) as usize].n[0];
            b[3] = point[(facet[ifs as usize].ip[1]) as usize].n[1];
            b[4] = point[(facet[ifs as usize].ip[2]) as usize].n[0];
            b[5] = point[(facet[ifs as usize].ip[2]) as usize].n[1];
            b[6] = point[(facet[ifs as usize].ip[0]) as usize].x[2];
            b[7] = point[(facet[ifs as usize].ip[1]) as usize].x[2];
            b[8] = point[(facet[ifs as usize].ip[2]) as usize].x[2];
            tetracgr(&a, &mut facet[ifs as usize].ct, &b);
        }
        ifs += 1;
    }
    write_facet(&facet, nfacet);
    // Calculate tetrahedron volumes; mark vanishing tetrahedra.
    itetra = 0;
    while (itetra < ntetra) {
        {
            tetra[itetra as usize].v = tetra_volume(
                &point[(tetra[itetra as usize].ip[0]) as usize].x,
                &point[(tetra[itetra as usize].ip[1]) as usize].x,
                &point[(tetra[itetra as usize].ip[2]) as usize].x,
                &point[(tetra[itetra as usize].ip[3]) as usize].x,
            );
            if ((tetra[itetra as usize].v as f64) < 1.0E-4f64) {
                tetra[itetra as usize].v = ((-9999.9f64) as f32);
            }
        }
        itetra += 1;
    }
    write_tetra(&tetra, ntetra);

    fclose(&hzfp);
    exit(CWP_Exit());
}

fn write_point(point: &mut [Point], npoint: i32) {
    let mut ip: i32 = 0;
    ip = 0;
    while (ip < npoint) {
        {
            if ((point[ip as usize].n[2] as f64) > 0f64) {
                point[ip as usize].n[0] = (-point[ip as usize].n[0]);
                point[ip as usize].n[1] = (-point[ip as usize].n[1]);
                point[ip as usize].n[2] = (-point[ip as usize].n[2]);
            }
            if ((point[ip as usize].n[2] as f64) == 0f64) {
                if ((point[ip as usize].n[0] as f64) < 0f64) {
                    point[ip as usize].n[0] = (-point[ip as usize].n[0]);
                    point[ip as usize].n[1] = (-point[ip as usize].n[1]);
                    point[ip as usize].n[2] = (-point[ip as usize].n[2]);
                } else {
                    if ((point[ip as usize].n[0] as f64) == 0f64) {
                        if ((point[ip as usize].n[1] as f64) < 0f64) {
                            point[ip as usize].n[0] = (-point[ip as usize].n[0]);
                            point[ip as usize].n[1] = (-point[ip as usize].n[1]);
                            point[ip as usize].n[2] = (-point[ip as usize].n[2]);
                        }
                    }
                }
            }
        }
        ip += 1;
    }
    fprintf!(stdout(), "%d=npoint : x[3] : s : n[3]\n", npoint);
    ip = 0;
    while (ip < npoint) {
        {
            fprintf!(
                stdout(),
                "%d : %f %f %f : %f : %f %f %f\n",
                ip,
                (point[ip as usize].x[0] as f64),
                (point[ip as usize].x[1] as f64),
                (point[ip as usize].x[2] as f64),
                (point[ip as usize].s as f64),
                (point[ip as usize].n[0] as f64),
                (point[ip as usize].n[1] as f64),
                (point[ip as usize].n[2] as f64)
            );
        }
        ip += 1;
    }
    fprintf!(stderr(), "%d point written\n", npoint);
}

fn write_facet(facet: &[Facet], nfacet: i32) {
    let mut ifs: i32 = 0;
    fprintf!(
        stdout(),
        "%d=nfacet : area : ct[10] : ip[3] : cn[3] : itetra[2]\n",
        nfacet
    );
    ifs = 0;
    while (ifs < nfacet) {
        {
            fprintf!(
                stdout(),
                "%d : %f : %f %f %f %f %f %f %f %f %f %f : %d %d %d : %f %f %f : %d %d\n",
                ifs,
                (facet[ifs as usize].area as f64),
                (facet[ifs as usize].ct[0] as f64),
                (facet[ifs as usize].ct[1] as f64),
                (facet[ifs as usize].ct[2] as f64),
                (facet[ifs as usize].ct[3] as f64),
                (facet[ifs as usize].ct[4] as f64),
                (facet[ifs as usize].ct[5] as f64),
                (facet[ifs as usize].ct[6] as f64),
                (facet[ifs as usize].ct[7] as f64),
                (facet[ifs as usize].ct[8] as f64),
                (facet[ifs as usize].ct[9] as f64),
                facet[ifs as usize].ip[0],
                facet[ifs as usize].ip[1],
                facet[ifs as usize].ip[2],
                (facet[ifs as usize].cn[0] as f64),
                (facet[ifs as usize].cn[1] as f64),
                (facet[ifs as usize].cn[2] as f64),
                facet[ifs as usize].itetra[0],
                facet[ifs as usize].itetra[1]
            );
        }
        ifs += 1;
    }
    fprintf!(stderr(), "%d facet written\n", nfacet);
}

fn write_tetra(tetra: &[Tetra], ntetra: i32) {
    let mut it: i32 = 0;
    fprintf!(
        stdout(),
        "%d=ntetra : v : ip[4] : ifacet[4] : ireg : gs[3]\n",
        ntetra
    );
    it = 0;
    while (it < ntetra) {
        {
            fprintf!(
                stdout(),
                "%d : %f : %d %d %d %d : %d %d %d %d : %d : %f %f %f\n",
                it,
                (tetra[it as usize].v as f64),
                tetra[it as usize].ip[0],
                tetra[it as usize].ip[1],
                tetra[it as usize].ip[2],
                tetra[it as usize].ip[3],
                tetra[it as usize].ifacet[0],
                tetra[it as usize].ifacet[1],
                tetra[it as usize].ifacet[2],
                tetra[it as usize].ifacet[3],
                tetra[it as usize].ireg,
                (tetra[it as usize].gs[0] as f64),
                (tetra[it as usize].gs[1] as f64),
                (tetra[it as usize].gs[2] as f64)
            );
        }
        it += 1;
    }
    fprintf!(stderr(), "%d tetra written\n", ntetra);
}

// Exterior facet indices are -1/-2, but C nevertheless reads their ireg.
// In the reference allocator, blocks up to 8192 bytes follow the preceding
// small allocation (facet, or point when facet uses a separate large block).
// Reproduce the initialized words and zero padding that those reads reach,
// without performing an out-of-bounds Rust access.
fn region(tetra: &[Tetra], point: &[Point], facet: &[Facet], i: i32) -> i32 {
    if i >= 0 {
        return tetra[i as usize].ireg;
    }
    if tetra.len() * 52 > 8192 {
        return 0;
    }
    let facet_bytes = facet.len() * 76;
    let point_bytes = point.len() * 28;
    let bytes = if facet_bytes <= 8192 {
        facet_bytes
    } else {
        point_bytes
    };
    if bytes > 8192 {
        return 0;
    }
    let quantum = (bytes.next_power_of_two() / 8).max(16);
    let allocated = bytes.div_ceil(quantum) * quantum;
    let offset = allocated as i32 + i * 52 + 32;
    if offset < 0 || offset as usize >= bytes {
        return 0;
    }
    let word = offset as usize / 4;
    if facet_bytes <= 8192 {
        let f = facet[word / 19];
        match word % 19 {
            0..=2 => f.ip[word % 19],
            3..=5 => f.cn[word % 19 - 3].to_bits() as i32,
            6..=7 => f.itetra[word % 19 - 6],
            8 => f.area.to_bits() as i32,
            _ => f.ct[word % 19 - 9].to_bits() as i32,
        }
    } else {
        let p = point[word / 7];
        match word % 7 {
            0..=2 => p.x[word % 7].to_bits() as i32,
            3 => p.s.to_bits() as i32,
            _ => p.n[word % 7 - 4].to_bits() as i32,
        }
    }
}

#[derive(Clone, Default)]
struct Grid2 {
    data: Vec<f32>,
    nx: usize,
}

impl Grid2 {
    fn new(nx: usize, ny: usize) -> Self {
        Self {
            data: vec![0.0; nx * ny],
            nx,
        }
    }
}
impl std::ops::Index<usize> for Grid2 {
    type Output = [f32];
    fn index(&self, iy: usize) -> &[f32] {
        &self.data[iy * self.nx..(iy + 1) * self.nx]
    }
}
impl std::ops::IndexMut<usize> for Grid2 {
    fn index_mut(&mut self, iy: usize) -> &mut [f32] {
        &mut self.data[iy * self.nx..(iy + 1) * self.nx]
    }
}

fn read_grid(grid: &mut Grid2, name: &str) -> usize {
    let fp = match fopen(name, "r") {
        Some(fp) => fp,
        None => c_segv(), // C calls fread on the null fopen result.
    };
    let n = fread(&mut grid.data, &fp);
    // As in C, input streams remain open until process exit.
    n
}

fn write_grid(grid: &Grid2, fp: &CFile) -> usize {
    fwrite(&grid.data, fp)
}

fn c_segv() -> ! {
    unsafe extern "C" {
        fn raise(sig: i32) -> i32;
        fn signal(sig: i32, handler: usize) -> usize;
    }
    // Rust installs a stack-overflow handler. The C binary has SIG_DFL.
    unsafe {
        signal(11, 0);
        raise(11);
    }
    exit(139)
}

// tm_normal - normal to the plane defined by the three vertices.
fn tm_normal(
    x0: f32,
    x1: f32,
    x2: f32,
    y0: f32,
    y1: f32,
    y2: f32,
    z0: f32,
    z1: f32,
    z2: f32,
) -> [f32; 3] {
    tm_normalize([
        (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0),
        (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0),
        (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0),
    ])
}

fn tm_normalize(mut n: [f32; 3]) -> [f32; 3] {
    let mut norm = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    norm = (norm as f64).sqrt() as f32;
    if norm != 0.0 {
        n[0] = n[0] / norm;
        n[1] = n[1] / norm;
        n[2] = n[2] / norm;
    } else {
        n = [0.0, 0.0, -1.0];
    }
    n
}

fn gradient(
    x0: [f32; 3],
    s0: f32,
    x1: [f32; 3],
    s1: f32,
    x2: [f32; 3],
    s2: f32,
    x3: [f32; 3],
    s3: f32,
) -> [f32; 3] {
    let mut a = [[0.0f32; 3]; 3];
    a[0][0] = x1[0] - x0[0];
    a[0][1] = x1[1] - x0[1];
    a[0][2] = x1[2] - x0[2];
    a[1][0] = x2[0] - x0[0];
    a[1][1] = x2[1] - x0[1];
    a[1][2] = x2[2] - x0[2];
    a[2][0] = x3[0] - x0[0];
    a[2][1] = x3[1] - x0[1];
    a[2][2] = x3[2] - x0[2];
    let mut gs = [s1 - s0, s2 - s0, s3 - s0];
    gaussj(&mut a, &mut gs);
    gs
}

// Gaussian Jordan solver, including the original pivot order.
fn gaussj(a: &mut [[f32; 3]; 3], b: &mut [f32; 3]) -> i32 {
    let mut icol = 0;
    let mut irow = 0;
    let mut indxc = [0usize; 3];
    let mut indxr = [0usize; 3];
    let mut ipiv = [0; 3];
    for j in 0..3 {
        ipiv[j] = 0;
    }
    for i in 0..3 {
        let mut big = 0.0f32;
        for j in 0..3 {
            if ipiv[j] != 1 {
                for k in 0..3 {
                    if ipiv[k] == 0 {
                        if (a[j][k] as f64).abs() >= big as f64 {
                            big = (a[j][k] as f64).abs() as f32;
                            irow = j;
                            icol = k;
                        }
                    } else if ipiv[k] > 1 {
                        for l in 0..3 {
                            b[l] = 0.0;
                        }
                        return -1;
                    }
                }
            }
        }
        ipiv[icol] += 1;
        if irow != icol {
            for l in 0..3 {
                let temp = a[irow][l];
                a[irow][l] = a[icol][l];
                a[icol][l] = temp;
            }
            let temp = b[irow];
            b[irow] = b[icol];
            b[icol] = temp;
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if a[icol][icol] == 0.0 {
            for l in 0..3 {
                b[l] = 0.0;
            }
            return -2;
        }
        let pivinv = (1.0f64 / a[icol][icol] as f64) as f32;
        a[icol][icol] = 1.0;
        for l in 0..3 {
            a[icol][l] *= pivinv;
        }
        b[icol] *= pivinv;
        for ll in 0..3 {
            if ll != icol {
                let dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for l in 0..3 {
                    a[ll][l] -= a[icol][l] * dum;
                }
                b[ll] -= b[icol] * dum;
            }
        }
    }
    for l in (0..3).rev() {
        if indxr[l] != indxc[l] {
            for k in 0..3 {
                let temp = a[k][indxr[l]];
                a[k][indxr[l]] = a[k][indxc[l]];
                a[k][indxc[l]] = temp;
            }
        }
    }
    1
}

// Conjugate gradient solver; np=9, maxiterations=5 in tetramod.c.
fn tetracgr(a: &[[f32; 9]; 9], ct: &mut [f32; 10], rhs: &[f32; 9]) {
    let mut x = [0.0f32; 9];
    let mut y = [0.0f32; 9];
    let mut s = [0.0f32; 9];
    let mut p = [0.0f32; 9];
    let mut r = [0.0f32; 9];
    let mut q = [0.0f32; 9];
    for ix in 0..6 {
        x[ix] = 0.0;
    }
    x[6] = ct[6];
    x[7] = ct[7];
    x[8] = ct[9];
    ax(a, &x, &mut y);
    for iy in 0..9 {
        s[iy] = rhs[iy] - y[iy];
    }
    ax(a, &s, &mut p);
    for iy in 0..9 {
        r[iy] = p[iy];
    }
    let mut normold = norm2(&r, 9);
    let mut iter = 0;
    while normold as f64 > 0.0001 {
        iter += 1;
        if iter == 5 {
            for ix in 0..6 {
                x[ix] = 0.0;
            }
            break;
        }
        ax(a, &p, &mut q);
        let alpha = normold / norm2(&q, 9);
        for ix in 0..9 {
            x[ix] += alpha * p[ix];
        }
        ax(a, &x, &mut y);
        let mut error = 0.0f32;
        for ix in 0..9 {
            let work = y[ix] - rhs[ix];
            error += work * work;
        }
        if (error as f64) < 0.0001 {
            break;
        }
        for ix in 0..9 {
            s[ix] += -alpha * q[ix];
        }
        ax(a, &s, &mut r);
        let normnew = norm2(&r, 9);
        let beta = normnew / normold;
        for ix in 0..9 {
            p[ix] = r[ix] + beta * p[ix];
        }
        normold = normnew;
    }
    let normnew = ((x[6] * x[6] + x[7] * x[7]) as f64 + 1.0).sqrt() as f32;
    for ix in 0..8 {
        ct[ix] = x[ix] / normnew;
    }
    ct[8] = (-1.0f64 / normnew as f64) as f32;
    ct[9] = x[8] / normnew;
}

fn ax(a: &[[f32; 9]; 9], x: &[f32; 9], y: &mut [f32; 9]) {
    for iy in 0..9 {
        y[iy] = 0.0;
        for ix in 0..9 {
            y[iy] += a[iy][ix] * x[ix];
        }
    }
}

fn norm2(xy: &[f32], n: usize) -> f32 {
    let mut sum = 0.0f32;
    for i in 0..n {
        sum += xy[i] * xy[i];
    }
    sum
}
