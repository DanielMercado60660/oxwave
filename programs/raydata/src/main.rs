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

//! RAYDATA - simple program to display ray data from elaray.
//! Port of Trielas/main/raydata.c.
//!
//! Credits (from the C source):
//! AUTHOR: Andreas Rueger, Colorado School of Mines, 02/15/94.

use oxwave_core::prelude::*;
use std::f64::consts::PI;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " RAYDATA - simple program to display ray data from elaray\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "raydata <rayends   [optional parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Optional Parameters: \t\t\t\t\t\t\t",
    "kend=      Index of interest\t\t\t\t\t\t",
    "t=0    =1 output of x_t requested \t\t\t\t\t",
    "px=0   =1 output of x_px requested \t\t\t\t\t",
    "pz=0   =1 output of x_pz requested \t\t\t\t\t",
    "vgx=0  =1 output of x_vgx requested \t\t\t\t\t",
    "vgz=0  =1 output of x_vgz requested \t\t\t\t\t",
    "pxvx=0 =1 output of px_vgx requested \t\t\t\t\t",
    "pol=0  =1 output of x_polangle requested \t\t\t\t",
    "ascci=0    binary output\t\t\t\t\t\t",
    "\t=1  ascci output\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Simple program to display some of the computed raydata.                ",
    "Output is written into files <*.data>                                  ",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

/// Native RayEnd layout from Trielas/include/elastic.h (not an XDR trace).
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct RayEnd {
    sigma: f32,
    x: f32,
    z: f32,
    px: f32,
    pz: f32,
    t: f32,
    q1: f32,
    p1: f32,
    q2: f32,
    p2: f32,
    kmah: i32,
    nref: i32,
    vgx: f32,
    vgz: f32,
    ab: f32,
    kend: i32,
    ampli: f32,
    ampliphase: f32,
    dangle: f32,
    mode: i32,
    g11: f32,
    g13: f32,
    g33: f32,
    nameref: i32,
    num: i32,
    rhob: f32,
    rhoe: f32,
}

// SAFETY: repr(C), only 32-bit floats/integers, no padding, all bit patterns valid.
unsafe impl oxwave_core::cio::Pod for RayEnd {}
const _: () = assert!(std::mem::size_of::<RayEnd>() == 108);

fn main() {
    let argv = args();
    let mut tout = 0;
    let mut ascci = 0;
    let mut pxout = 0;
    let mut pzout = 0;
    let mut vgxout = 0;
    let mut vgzout = 0;
    let mut pxvxout = 0;
    let mut polout = 0;
    let mut tdata = Vec::<f32>::new();
    let mut pxdata = Vec::<f32>::new();
    let mut pzdata = Vec::<f32>::new();
    let mut vgxdata = Vec::<f32>::new();
    let mut vgzdata = Vec::<f32>::new();
    let mut poldata = Vec::<f32>::new();
    let mut g11data = Vec::<f32>::new();
    let mut g33data = Vec::<f32>::new();
    let mut g13data = Vec::<f32>::new();
    let mut kend = 0;

    /* hook up getpar to handle the parameters */
    initargs(&argv);
    requestdoc(0, SDOC);

    /* get parameters */
    if getparint("kend", &mut kend) == 0 {
        kend = i32::MAX;
    }
    if getparint("t", &mut tout) == 0 {
        tout = 0;
    }
    if getparint("px", &mut pxout) == 0 {
        pxout = 0;
    }
    if getparint("pz", &mut pzout) == 0 {
        pzout = 0;
    }
    if getparint("vgx", &mut vgxout) == 0 {
        vgxout = 0;
    }
    if getparint("vgz", &mut vgzout) == 0 {
        vgzout = 0;
    }
    if getparint("pxvx", &mut pxvxout) == 0 {
        pxvxout = 0;
    }
    if getparint("pol", &mut polout) == 0 {
        polout = 0;
    }
    if getparint("ascci", &mut ascci) == 0 {
        ascci = 0;
    }
    checkpars();

    /* output file control */
    let txfp = if tout > 0 {
        fopen("x_t.data", "w")
    } else {
        None
    };
    let pxfp = if pxout != 0 {
        fopen("x_px.data", "w")
    } else {
        None
    };
    let pzfp = if pzout != 0 {
        fopen("x_pz.data", "w")
    } else {
        None
    };
    let vgxfp = if vgxout != 0 {
        fopen("x_vgx.data", "w")
    } else {
        None
    };
    let vgzfp = if vgzout != 0 {
        fopen("x_vgz.data", "w")
    } else {
        None
    };
    let pxvxfp = if pxvxout != 0 {
        fopen("px_vgx.data", "w")
    } else {
        None
    };
    let polfp = if polout != 0 {
        fopen("x_pol.data", "w")
    } else {
        None
    };
    let outparfp = if ascci == 0 {
        Some(efopen("outpar", "w"))
    } else {
        None
    };

    /* read rayends */
    let mut nre: i32 = 0;
    let mut nri: i32 = 0;
    let mut nrealloc: i32 = 301;
    let mut re = vec![RayEnd::default(); nrealloc as usize];
    while fread(&mut re[nre as usize..nre as usize + 1], &stdin()) == 1 {
        nre += 1;
        if nre == nrealloc {
            nrealloc += 100;
            re.resize(nrealloc as usize, RayEnd::default());
        }
    }
    if kend == i32::MAX {
        nri = nre;
    } else {
        /* how many rayends are of interest */
        for ir in 0..nre {
            if re[ir as usize].kend == kend {
                nri += 1;
            }
        }
    }

    /* allocate space for data files */
    if tout > 0 {
        tdata = vec![0.0; nri as usize];
    }
    if pxout != 0 || pxvxout != 0 {
        pxdata = vec![0.0; nri as usize];
    }
    if pzout != 0 {
        pzdata = vec![0.0; nri as usize];
    }
    if vgxout != 0 || pxvxout != 0 {
        vgxdata = vec![0.0; nri as usize];
    }
    if vgzout != 0 {
        vgzdata = vec![0.0; nri as usize];
    }
    if polout != 0 {
        // C allocates nri but computes nre polarizations, overrunning all four
        // arrays when kend filters records. Pad the unused tail to avoid UB;
        // only the first nri values are written. C also leaves poldata
        // uninitialized if none of its four conditions holds (zeroed here).
        poldata = vec![0.0; nre as usize];
        g11data = vec![0.0; nre as usize];
        g13data = vec![0.0; nre as usize];
        g33data = vec![0.0; nre as usize];
    }
    let mut xdata = vec![0.0f32; nri as usize];
    let mut iri: i32 = 0;

    /* read in data into files */
    for ir in 0..nre {
        let ir = ir as usize;
        if re[ir].kend == kend || kend == i32::MAX {
            let i = iri as usize;
            if tout > 0 {
                tdata[i] = re[ir].t;
            }
            if pxout != 0 || pxvxout != 0 {
                pxdata[i] = re[ir].px;
            }
            if pzout != 0 {
                pzdata[i] = re[ir].pz;
            }
            if vgxout != 0 || pxvxout != 0 {
                vgxdata[i] = re[ir].vgx;
            }
            if vgzout != 0 {
                vgzdata[i] = re[ir].vgz;
            }
            if polout != 0 {
                g11data[i] = re[ir].g11;
                g33data[i] = re[ir].g33;
                g13data[i] = re[ir].g13;
            }
            xdata[i] = re[ir].x;
            iri += 1;
        }
    }

    /* compute polarization */
    if polout != 0 {
        for ir in 0..nre {
            let ir = ir as usize;
            // C divides floats before promoting the result for sqrt and atan.
            if g13data[ir] > 0.0 {
                poldata[ir] = ((g11data[ir] / g33data[ir]) as f64).sqrt().atan() as f32;
            } else if g13data[ir] < 0.0 {
                poldata[ir] =
                    (PI / 2.0 + ((g33data[ir] / g11data[ir]) as f64).sqrt().atan()) as f32;
            } else if g13data[ir] == 0.0 && g33data[ir] > 0.0 {
                poldata[ir] = 0.0;
            } else if g13data[ir] == 0.0 && g33data[ir] < 0.0 {
                poldata[ir] = (PI / 2.0) as f32;
            }
        }
    }

    /* ASCII output for x_t */
    if ascci == 1 && tout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fprintf!(txfp.as_ref().unwrap(), "%f\t%f\n", xdata[ir], tdata[ir]);
        }
    /* Binary output for x_t */
    } else if ascci == 0 && tout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fwrite(&[tdata[ir]], txfp.as_ref().unwrap());
            fwrite(&[xdata[ir]], txfp.as_ref().unwrap());
        }
    }

    /* ASCII output for x_px */
    if ascci == 1 && pxout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fprintf!(pxfp.as_ref().unwrap(), "%f\t%f\n", xdata[ir], pxdata[ir]);
        }
    /* Binary output for x_px */
    } else if ascci == 0 && pxout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fwrite(&[pxdata[ir]], pxfp.as_ref().unwrap());
            fwrite(&[xdata[ir]], pxfp.as_ref().unwrap());
        }
    }

    /* ASCII output for x_pz */
    if ascci == 1 && pzout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fprintf!(pzfp.as_ref().unwrap(), "%f\t%f\n", xdata[ir], pzdata[ir]);
        }
    /* Binary output for x_pz */
    } else if ascci == 0 && pzout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fwrite(&[pzdata[ir]], pzfp.as_ref().unwrap());
            fwrite(&[xdata[ir]], pzfp.as_ref().unwrap());
        }
    }

    /* ASCII output for x_vgx */
    if ascci == 1 && vgxout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fprintf!(vgxfp.as_ref().unwrap(), "%f\t%f\n", xdata[ir], vgxdata[ir]);
        }
    /* Binary output for x_vgx */
    } else if ascci == 0 && vgxout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fwrite(&[vgxdata[ir]], vgxfp.as_ref().unwrap());
            fwrite(&[xdata[ir]], vgxfp.as_ref().unwrap());
        }
    }

    /* ASCII output for x_vgz */
    if ascci == 1 && vgzout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fprintf!(vgzfp.as_ref().unwrap(), "%f\t%f\n", xdata[ir], vgzdata[ir]);
        }
    /* Binary output for x_vgz */
    } else if ascci == 0 && vgzout == 1 {
        for ir in 0..nri {
            let ir = ir as usize;
            fwrite(&[vgzdata[ir]], vgzfp.as_ref().unwrap());
            fwrite(&[xdata[ir]], vgzfp.as_ref().unwrap());
        }
    }

    /* ASCII output for px_vgx */
    if ascci == 1 && pxvxout != 0 {
        for ir in 0..nri {
            let ir = ir as usize;
            fprintf!(
                pxvxfp.as_ref().unwrap(),
                "%f\t%f\n",
                pxdata[ir],
                vgxdata[ir]
            );
        }
    /* Binary output for px_vgx */
    } else if ascci == 0 && pxvxout != 0 {
        for ir in 0..nri {
            let ir = ir as usize;
            fwrite(&[pxdata[ir]], pxvxfp.as_ref().unwrap());
            fwrite(&[vgxdata[ir]], pxvxfp.as_ref().unwrap());
        }
    }

    /* ASCII output for x_pol */
    if ascci == 1 && polout != 0 {
        for ir in 0..nri {
            let ir = ir as usize;
            fprintf!(polfp.as_ref().unwrap(), "%f\t%f\n", xdata[ir], poldata[ir]);
        }
    /* Binary output for x_pol */
    } else if ascci == 0 && polout != 0 {
        for ir in 0..nri {
            let ir = ir as usize;
            fwrite(&[xdata[ir]], polfp.as_ref().unwrap());
            fwrite(&[poldata[ir]], polfp.as_ref().unwrap());
        }
    }

    if ascci == 0 {
        fprintf!(outparfp.as_ref().unwrap(), "%i\n", nri);
    }
    // The C program returns EXIT_FAILURE even after writing every output file.
    exit(EXIT_FAILURE);
}
