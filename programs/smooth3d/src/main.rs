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

//! SMOOTH3D - 3D grid velocity smoothing by damped least squares.
//! Port of par/main/smoothing/smooth3d.c (Revision 1.10).
//! Author: CWP: Zhenyue Liu, March 1995.
//! Reference: Liu, Z., 1994, A velocity smoothing technique based on damped
//! least squares, Project Reveiew, May 10, 1994, Consortium Project on
//! Seismic Inverse Methods for Complex Stuctures.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SMOOTH3D - 3D grid velocity SMOOTHing by the damped least squares\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tsmooth3d <infile >outfile [parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1=\t\tnumber of samples along 1st dimension\t\t\t",
    " n2=\t\tnumber of samples along 2nd dimension\t\t\t",
    " n3=\t\tnumber of samples along 3rd dimension\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Smoothing parameters (0 = no smoothing)\t\t\t\t",
    " r1=0.0\toperator length in 1st dimension\t\t\t",
    " r2=0.0\toperator length in 2nd dimension\t\t\t",
    " r3=0.0\toperator length in 3rd dimension\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Sample intervals:\t\t\t\t\t\t\t",
    " d1=1.0\t1st dimension\t\t\t\t\t\t",
    " d2=1.0\t2nd dimension\t\t\t\t\t\t",
    " d3=1.0\t3rd dimension\t\t\t\t\t\t",
    " iter=2\tnumber of iteration used\t\t\t\t",
    " time=0\twhich dimension the time axis is (0 = no time axis)\t",
    " depth=1\twhich dimension the depth axis is (ignored when time>0)\t",
    " mu=1\t\tthe relative weight at maximum depth (or time)\t\t",
    " verbose=0\t=1 for printing minimum wavelengths\t\t\t",
    " slowness=0\t=1 smoothing on slowness; =0 smoothing on velocity\t",
    " vminc=0\tvelocity values below it are cliped before smoothing\t",
    " vmaxc=99999\tvelocity values above it are cliped before smoothing\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " 1. The larger the smoothing parameters, the smoother the output velocity.",
    "    These parameters are lengths of smoothing operators in each dimension.",
    " 2. iter controls the orders of derivatives to be smoothed in the output",
    "    velocity. e.g., iter=2 means derivatives until 2nd order smoothed.\t",
    " 3. mu is the multipler of smoothing parameters at the bottom compared to",
    "    those at the surface.\t\t\t\t\t\t",
    " 4. Minimum wavelengths of each dimension and the total may be printed\t",
    "   for the resulting output velocity is. To compute these parameters for",
    "   the input velocity, use r1=r2=r3=0.\t\t\t\t\t",
    " 5. Smoothing directly on slowness works better to preserve traveltime.",
    "   So the program optionally converts the input velocity into slowness\t",
    "   and smooths the slowness, then converts back into velocity.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let invp = stdin();
    let outvp = stdout();
    initargs(&argv);
    requestdoc(0, SDOC);

    /* Required parameters. */
    let (mut n1, mut n2, mut n3) = (0i32, 0i32, 0i32);
    if getparint("n1", &mut n1) == 0 {
        n1 = 0;
    }
    if n1 <= 0 {
        err!("sample number of 1st dimension invalid");
    }
    if getparint("n2", &mut n2) == 0 {
        n2 = 0;
    }
    if n2 <= 0 {
        err!("sample number of 2nd dimension invalid");
    }
    if getparint("n3", &mut n3) == 0 {
        n3 = 0;
    }
    if n3 <= 0 {
        err!("sample number of 3rd dimension invalid");
    }

    /* Sample intervals and smoothing parameters. */
    let (mut d1, mut d2, mut d3) = (1.0f32, 1.0f32, 1.0f32);
    if getparfloat("d1", &mut d1) == 0 {
        d1 = 1.0;
    }
    if getparfloat("d2", &mut d2) == 0 {
        d2 = 1.0;
    }
    if getparfloat("d3", &mut d3) == 0 {
        d3 = 1.0;
    }
    let (mut r1, mut r2, mut r3) = (0.0f32, 0.0f32, 0.0f32);
    if getparfloat("r1", &mut r1) == 0 || n1 < 4 {
        r1 = 0.0;
    }
    if getparfloat("r2", &mut r2) == 0 || n2 < 4 {
        r2 = 0.0;
    }
    if getparfloat("r3", &mut r3) == 0 || n3 < 4 {
        r3 = 0.0;
    }
    r1 = if d1 > 0.0 { r1 / d1 } else { 0.0 };
    r2 = if d2 > 0.0 { r2 / d2 } else { 0.0 };
    r3 = if d3 > 0.0 { r3 / d3 } else { 0.0 };
    r1 = (0.5 * r1 as f64 * r1 as f64) as f32;
    r2 = (0.5 * r2 as f64 * r2 as f64) as f32;
    r3 = (0.5 * r3 as f64 * r3 as f64) as f32;

    let mut iter = 2i32;
    if getparint("iter", &mut iter) == 0 {
        iter = 2;
    }
    if iter <= 0 || iter > 3 {
        err!("\t iter must be between 1 and 3!\n");
    }
    let (mut time, mut depth) = (0i32, 1i32);
    if getparint("time", &mut time) == 0 {
        time = 0;
    }
    if getparint("depth", &mut depth) == 0 {
        depth = 1;
    }
    if time != 0 {
        depth = time;
    }
    let mut mu = 1.0f32;
    if getparfloat("mu", &mut mu) == 0 {
        mu = 1.0;
    }
    if mu < 1.0 {
        err!("mu must not be less than 1 \n");
    }
    let mut slowness = 0i32;
    if getparint("slowness", &mut slowness) == 0 {
        slowness = 0;
    }
    let (mut vminc, mut vmaxc) = (0.0f32, 99999.0f32);
    if getparfloat("vminc", &mut vminc) == 0 {
        vminc = 0.0;
    }
    if getparfloat("vmaxc", &mut vmaxc) == 0 {
        vmaxc = 99999.0;
    }

    /* Contiguous storage: C reads and writes vel[0][0] as a single block. */
    let mut vel = vec![0.0f32; n1 as usize * n2 as usize * n3 as usize];
    let mut wl = [0.0f32; 4];
    efread(&mut vel, &invp);
    vsm3d(
        &mut vel,
        n3 as usize,
        n2 as usize,
        n1 as usize,
        iter,
        depth,
        r3,
        r2,
        r1,
        mu,
        slowness,
        vminc,
        vmaxc,
    );
    efwrite(&vel, &outvp);

    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    if verbose != 0 {
        wavel(
            n1 as usize,
            n2 as usize,
            n3 as usize,
            d1,
            d2,
            d3,
            time,
            &mut wl,
            &vel,
        );
        fprintf!(stderr(), "minimum wavelengths of smoothed velocity:\n");
        fprintf!(stderr(), "\tlambda1 = %g,\n", wl[1]);
        fprintf!(stderr(), "\tlambda2 = %g,\n", wl[2]);
        fprintf!(stderr(), "\tlambda3 = %g,\n", wl[3]);
        fprintf!(stderr(), "\tlambda = %g,\n", wl[0]);
    }
    checkpars();
    efclose(&invp);
    efclose(&outvp);
    exit(CWP_Exit());
}

/// Compute minimum wavelengths for individual dimensions and total.
#[allow(clippy::too_many_arguments)]
fn wavel(
    n1: usize,
    n2: usize,
    n3: usize,
    d1: f32,
    d2: f32,
    d3: f32,
    time: i32,
    wl: &mut [f32; 4],
    v: &[f32],
) {
    let (mut k, mut k1, mut k2, mut k3) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let mut scal = 0.0f32;
    let od1 = (1.0 / d1 as f64) as f32;
    let od2 = (1.0 / d2 as f64) as f32;
    let od3 = (1.0 / d3 as f64) as f32;
    let (o2d1, o2d2, o2d3) = (od1 * od1, od2 * od2, od3 * od3);
    for i3 in 0..n3 {
        for i2 in 0..n2 {
            for i1 in 0..n1 {
                let p = (i3 * n2 + i2) * n1 + i1;
                let vel = v[p];
                let ovel = (1.0 / vel as f64) as f32;
                let (mut d2v1, mut ddv1);
                if i1 == 0 || i1 == n1 - 1 {
                    d2v1 = 0.0;
                    ddv1 = 0.0;
                } else {
                    let dv1 = (vel - v[p - 1]) * od1 * ovel;
                    d2v1 = dv1 * dv1;
                    ddv1 = ((v[p - 1] as f64 - 2.0 * vel as f64 + v[p + 1] as f64)
                        * o2d1 as f64
                        * ovel as f64) as f32;
                    if ddv1 < 0.0 {
                        ddv1 = -ddv1;
                    }
                }
                let (mut d2v2, mut ddv2);
                if i2 == 0 || i2 == n2 - 1 {
                    d2v2 = 0.0;
                    ddv2 = 0.0;
                } else {
                    let dv2 = (vel - v[p - n1]) * od2 * ovel;
                    d2v2 = dv2 * dv2;
                    ddv2 = ((v[p - n1] as f64 - 2.0 * vel as f64 + v[p + n1] as f64)
                        * o2d2 as f64
                        * ovel as f64) as f32;
                    if ddv2 < 0.0 {
                        ddv2 = -ddv2;
                    }
                }
                let (mut d2v3, mut ddv3);
                if i3 == 0 || i3 == n3 - 1 {
                    d2v3 = 0.0;
                    ddv3 = 0.0;
                } else {
                    let dv3 = (vel - v[p - n1 * n2]) * od3 * ovel;
                    d2v3 = dv3 * dv3;
                    ddv3 = ((v[p - n1 * n2] as f64 - 2.0 * vel as f64 + v[p + n1 * n2] as f64)
                        * o2d3 as f64
                        * ovel as f64) as f32;
                    if ddv3 < 0.0 {
                        ddv3 = -ddv3;
                    }
                }
                if time != 0 {
                    scal = (2000.0 * ovel as f64) as f32;
                    scal = scal * scal;
                }
                if time == 1 {
                    d2v1 *= scal;
                    ddv1 *= scal;
                } else if time == 2 {
                    d2v2 *= scal;
                    ddv2 *= scal;
                } else if time == 3 {
                    d2v3 *= scal;
                    ddv3 *= scal;
                }
                if k1 < d2v1 {
                    k1 = d2v1;
                }
                if k1 < ddv1 {
                    k1 = ddv1;
                }
                if k2 < d2v2 {
                    k2 = d2v2;
                }
                if k2 < ddv2 {
                    k2 = ddv2;
                }
                if k3 < d2v3 {
                    k3 = d2v3;
                }
                if k3 < ddv3 {
                    k3 = ddv3;
                }
                let dv = d2v1 + d2v2 + d2v3;
                let ddv = ddv1 + ddv2 + ddv3;
                if k < dv {
                    k = dv;
                }
                if k < ddv {
                    k = ddv;
                }
            }
        }
    }
    wl[0] = if k != 0.0 {
        (2.0 * PI / (k as f64).sqrt()) as f32
    } else {
        -1.0
    };
    wl[1] = if k1 != 0.0 {
        (2.0 * PI / (k1 as f64).sqrt()) as f32
    } else {
        -1.0
    };
    wl[2] = if k2 != 0.0 {
        (2.0 * PI / (k2 as f64).sqrt()) as f32
    } else {
        -1.0
    };
    wl[3] = if k3 != 0.0 {
        (2.0 * PI / (k3 as f64).sqrt()) as f32
    } else {
        -1.0
    };
}

/// Smooth 3d velocity. The three passes deliberately run in order 2, 3, 1.
#[allow(clippy::too_many_arguments)]
fn vsm3d(
    v: &mut [f32],
    n3: usize,
    n2: usize,
    n1: usize,
    iter: i32,
    depth: i32,
    mut r3: f32,
    mut r2: f32,
    mut r1: f32,
    mut mu: f32,
    sl: i32,
    vmin: f32,
    vmax: f32,
) {
    let mut ww = 1.0f32;
    /* Compute the weight function. */
    let mut w = vec![0.0f32; n1 + n2 + n3 - 2];
    if depth == 1 {
        mu = (((mu * mu) as f64 - 1.0) / (n1 as i32).wrapping_mul(n1 as i32) as f64) as f32;
        for i1 in 0..n1 {
            w[i1] =
                (1.0 / (1.0f32 + (i1 as i32).wrapping_mul(i1 as i32) as f32 * mu) as f64) as f32;
        }
    }
    if depth == 2 {
        mu = (((mu * mu) as f64 - 1.0) / (n2 as i32).wrapping_mul(n2 as i32) as f64) as f32;
        for i2 in 0..n2 {
            w[i2] =
                (1.0 / (1.0f32 + (i2 as i32).wrapping_mul(i2 as i32) as f32 * mu) as f64) as f32;
        }
    }
    if depth == 3 {
        mu = (((mu * mu) as f64 - 1.0) / (n3 as i32).wrapping_mul(n3 as i32) as f64) as f32;
        for i3 in 0..n3 {
            w[i3] =
                (1.0 / (1.0f32 + (i3 as i32).wrapping_mul(i3 as i32) as f32 * mu) as f64) as f32;
        }
    }
    /* Scale smoothing parameters according to the iteration number. */
    if iter == 1 {
        r1 = (r1 as f64 / (3.39 * 3.39)) as f32;
        r2 = (r2 as f64 / (3.39 * 3.39)) as f32;
        r3 = (r3 as f64 / (3.39 * 3.39)) as f32;
    } else if iter == 2 {
        r1 = (r1 as f64 / (5.19 * 5.19)) as f32;
        r2 = (r2 as f64 / (5.19 * 5.19)) as f32;
        r3 = (r3 as f64 / (5.19 * 5.19)) as f32;
    } else {
        r1 = (r1 as f64 / (6.60 * 6.60)) as f32;
        r2 = (r2 as f64 / (6.60 * 6.60)) as f32;
        r3 = (r3 as f64 / (6.60 * 6.60)) as f32;
    }
    /* Clip velocity. */
    for i3 in 0..n3 {
        for i2 in 0..n2 {
            for i1 in 0..n1 {
                let p = (i3 * n2 + i2) * n1 + i1;
                if v[p] > vmax {
                    v[p] = vmax;
                }
                if v[p] < vmin {
                    v[p] = vmin;
                }
            }
        }
    }
    if sl != 0 {
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                for i1 in 0..n1 {
                    let p = (i3 * n2 + i2) * n1 + i1;
                    v[p] = (1.0 / v[p] as f64) as f32;
                }
            }
        }
    }

    if r2 > 0.0 {
        /* Smoothing velocity in the second direction. */
        let mut d = vec![vec![0.0f32; n1]; n2];
        let mut e = vec![vec![0.0f32; n1]; n2];
        let mut f = vec![vec![0.0f32; n1]; n2];
        for i3 in 0..n3 {
            if depth == 3 {
                ww = w[i3];
            }
            for i2 in 0..n2 - 1 {
                if depth == 2 {
                    ww = w[i2 + 1];
                }
                for i1 in 0..n1 {
                    if depth == 1 {
                        ww = w[i1];
                    }
                    d[i2][i1] = (ww as f64 + r2 as f64 * 2.0) as f32;
                    e[i2][i1] = -r2;
                    f[i2][i1] = ww * v[(i3 * n2 + i2 + 1) * n1 + i1];
                }
            }
            for i1 in 0..n1 {
                d[n2 - 2][i1] -= r2;
                f[0][i1] += r2 * v[i3 * n2 * n1 + i1];
            }
            tripd2(&mut d, &mut e, &mut f, n2 - 1, n1);
            for _i in 1..iter {
                for i2 in 0..n2 - 1 {
                    if depth == 2 {
                        ww = w[i2 + 1];
                    }
                    for i1 in 0..n1 {
                        if depth == 1 {
                            ww = w[i1];
                        }
                        d[i2][i1] = (ww as f64 + r2 as f64 * 2.0) as f32;
                        e[i2][i1] = -r2;
                        f[i2][i1] *= ww;
                    }
                }
                for i1 in 0..n1 {
                    d[n2 - 2][i1] -= r2;
                    f[0][i1] += r2 * v[i3 * n2 * n1 + i1];
                }
                tripd2(&mut d, &mut e, &mut f, n2 - 1, n1);
            }
            for i2 in 0..n2 - 1 {
                for i1 in 0..n1 {
                    v[(i3 * n2 + i2 + 1) * n1 + i1] = f[i2][i1];
                }
            }
        }
    }

    if r3 > 0.0 {
        /* Smooth velocity in the third direction. */
        let mut d = vec![vec![0.0f32; n1]; n3];
        let mut e = vec![vec![0.0f32; n1]; n3];
        let mut f = vec![vec![0.0f32; n1]; n3];
        for i2 in 0..n2 {
            if depth == 2 {
                ww = w[i2];
            }
            for i3 in 0..n3 - 1 {
                if depth == 3 {
                    ww = w[i3 + 1];
                }
                for i1 in 0..n1 {
                    if depth == 1 {
                        ww = w[i1];
                    }
                    d[i3][i1] = (ww as f64 + 2.0 * r3 as f64) as f32;
                    e[i3][i1] = -r3;
                    f[i3][i1] = ww * v[((i3 + 1) * n2 + i2) * n1 + i1];
                }
            }
            for i1 in 0..n1 {
                d[n3 - 2][i1] -= r3;
                f[0][i1] += r3 * v[i2 * n1 + i1];
            }
            tripd2(&mut d, &mut e, &mut f, n3 - 1, n1);
            for _i in 1..iter {
                for i3 in 0..n3 - 1 {
                    if depth == 3 {
                        ww = w[i3 + 1];
                    }
                    for i1 in 0..n1 {
                        if depth == 1 {
                            ww = w[i1];
                        }
                        d[i3][i1] = (ww as f64 + 2.0 * r3 as f64) as f32;
                        e[i3][i1] = -r3;
                        f[i3][i1] *= ww;
                    }
                }
                for i1 in 0..n1 {
                    d[n3 - 2][i1] -= r3;
                    f[0][i1] += r3 * v[i2 * n1 + i1];
                }
                tripd2(&mut d, &mut e, &mut f, n3 - 1, n1);
            }
            for i3 in 0..n3 - 1 {
                for i1 in 0..n1 {
                    v[((i3 + 1) * n2 + i2) * n1 + i1] = f[i3][i1];
                }
            }
        }
    }

    if r1 > 0.0 {
        /* Smooth velocity in the first direction. */
        let mut d = vec![vec![0.0f32; 1]; n1];
        let mut e = vec![vec![0.0f32; 1]; n1];
        let mut f = vec![vec![0.0f32; 1]; n1];
        for i3 in 0..n3 {
            if depth == 3 {
                ww = w[i3];
            }
            for i2 in 0..n2 {
                if depth == 2 {
                    ww = w[i2];
                }
                for i1 in 0..n1 - 1 {
                    if depth == 1 {
                        ww = w[i1 + 1];
                    }
                    d[i1][0] = (ww as f64 + r1 as f64 * 2.0) as f32;
                    e[i1][0] = -r1;
                    f[i1][0] = ww * v[(i3 * n2 + i2) * n1 + i1 + 1];
                }
                d[n1 - 2][0] -= r1;
                f[0][0] += r1 * v[(i3 * n2 + i2) * n1];
                tripd2(&mut d, &mut e, &mut f, n1 - 1, 1);
                for _i in 1..iter {
                    for i1 in 0..n1 - 1 {
                        if depth == 1 {
                            ww = w[i1 + 1];
                        }
                        d[i1][0] = (ww as f64 + r1 as f64 * 2.0) as f32;
                        e[i1][0] = -r1;
                        f[i1][0] *= ww;
                    }
                    d[n1 - 2][0] -= r1;
                    f[0][0] += r1 * v[(i3 * n2 + i2) * n1];
                    tripd2(&mut d, &mut e, &mut f, n1 - 1, 1);
                }
                for i1 in 0..n1 - 1 {
                    v[(i3 * n2 + i2) * n1 + i1 + 1] = f[i1][0];
                }
            }
        }
    }
    if sl != 0 {
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                for i1 in 0..n1 {
                    let p = (i3 * n2 + i2) * n1 + i1;
                    v[p] = (1.0 / v[p] as f64) as f32;
                }
            }
        }
    }
}

/// Solve m symmetric, positive definite tridiagonal n-by-n systems.
/// This helper is defined in smooth3d.c itself, not in a CWP library.
fn tripd2(d: &mut [Vec<f32>], e: &mut [Vec<f32>], b: &mut [Vec<f32>], n: usize, m: usize) {
    /* Decomposition. */
    for k in 1..n {
        for i in 0..m {
            let temp = e[k - 1][i];
            e[k - 1][i] = temp / d[k - 1][i];
            d[k][i] -= temp * e[k - 1][i];
        }
    }
    /* Substitution. */
    for k in 1..n {
        for i in 0..m {
            b[k][i] -= e[k - 1][i] * b[k - 1][i];
        }
    }
    for i in 0..m {
        b[n - 1][i] /= d[n - 1][i];
    }
    for k in (1..n).rev() {
        for i in 0..m {
            b[k - 1][i] = b[k - 1][i] / d[k - 1][i] - e[k - 1][i] * b[k][i];
        }
    }
}
