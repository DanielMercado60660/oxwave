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

//! ADDRVL3D - Add a random velocity layer to a gridded velocity model.
//! Port of par/main/velocity_model_building/addrvl3d.c (Revision 1.1).
//! Copyright (c) Colorado School of Mines, 2005. All rights reserved.
//! Author: Saudi Aramco: Chris Liner Jan/Feb 2005.
//! Based on smooth3d (CWP: Zhenyue Liu March 1995).

use oxwave_core::prelude::*;
use oxwave_f_franuni::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " ADDRVL3D - Add a random velocity layer (RVL) to a gridded             ",
    "            v(x,y,z) velocity model                                    ",
    "\t\t\t\t\t\t\t\t\t",
    "\taddrvl3d <infile n1= n2= >outfile [parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1=\t\tnumber of samples along 1st dimension\t\t\t",
    " n2=\t\tnumber of samples along 2nd dimension\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " n3=1          number of samples along 3rd dimension\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " mode=1             add single layer populated with random vels\t",
    "                    =2 add nrvl layers of random thickness and vel     ",
    " seed=from_clock    random number seed (integer)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "         ********** mode=1 params **********\t\t\t\t",
    " ---->New layer geometry info\t\t\t\t\t\t",
    " i1beg=1       1st dimension beginning sample \t\t\t\t",
    " i1end=n1/5    1st dimension ending sample \t\t\t\t",
    " i2beg=1       2nd dimension beginning sample \t\t\t\t",
    " i2end=n2      2nd dimension ending sample \t\t\t\t",
    " i3beg=1       3rd dimension beginning sample \t\t\t\t",
    " i3end=n3      3rd dimension ending sample \t\t\t\t",
    " ---->New layer velocity info\t\t\t\t\t\t",
    " vlsd=v/3     range (std dev) of random velocity in layer, \t\t",
    "               where v=v(0,0,i1) and i1=(i1beg+i1end)/2 \t \t",
    " add=1         add random vel to original vel (v_orig) at that point \t",
    "               =0 replace vel at that point with (v_orig+v_rand) \t",
    " how=0         random vels can be higher or lower than v_orig\t\t",
    "               =1 random vels are always lower than v_orig\t\t",
    "               =2 random vels are always higher than v_orig\t\t",
    " cvel=2000     layer filled with constant velocity cvel \t\t",
    "               (overides vlsd,add,how params)\t\t\t",
    " ---->Smoothing parameters (0 = no smoothing)\t\t\t\t",
    " r1=0.0\t1st dimension operator length in samples\t\t",
    " r2=0.0\t2nd dimension operator length in samples\t\t",
    " r3=0.0\t3rd dimension operator length in samples\t\t",
    " slowness=0\t=1 smoothing on slowness; =0 smoothing on velocity\t",
    "\t\t\t\t\t\t\t\t\t",
    "         ********** mode=2 params **********\t\t\t\t",
    " nrvl=n1/10    number of const velocity layers to add     \t\t",
    " pdv=10.       percentage velocity deviation (max) from input model\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " *** mode=1 **********\t\t\t\t\t\t\t",
    " 1. Smoothing radii usually fall in the range of [0,20].\t\t",
    " 2. Smoothing radii can be used to set aspect ratio of random velocity ",
    "    anomalies in the new layer.  For example (r1=5,r2=0,r3=0) will     ",
    "    result in vertical vel streaks that mimick vertical fracturing.    ",
    " 3. Smoothing on slowness works better to preserve traveltimes relative",
    "    to the unsmoothed case.\t\t\t\t\t\t",
    " 4. Default case is a random velocity (+/-30%) near surface layer whose",
    "    thickness is 20% of the total 2D model thickness.\t\t\t",
    " *** mode=2 **********\t\t\t                          \t",
    " 5. Each layer vel is a random perturbation on input model at that level.",
    " *** both modes ******\t\t\t\t",
    " 6. The depth dimension is assumed to be along axis 1.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    " 1. 2D RVL with no smoothing\t\t\t\t\t\t",
    "   makevel nz=250 nx=200 | addrvl3d n1=250 n2=200 | ximage n1=250      ",
    " 2. 3D RVL with no smoothing\t\t\t\t\t\t",
    "   makevel nz=250 nx=200 ny=220 |\t\t\t\t\t",
    "   addrvl3d n1=250 n2=200 n3=220 | \t\t\t\t\t",
    "   xmovie n1=250 n2=200\t\t\t\t\t    \t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut n1 = 0i32;
    let mut n2 = 0i32;
    getparint("n1", &mut n1);
    if n1 <= 0 {
        err!("sample number of 1st dimension invalid");
    }
    getparint("n2", &mut n2);
    if n2 <= 0 {
        err!("sample number of 2nd dimension invalid");
    }
    let mut n3 = 1i32;
    if getparint("n3", &mut n3) == 0 {
        n3 = 1;
    }
    let mut mode = 1i32;
    if getparint("mode", &mut mode) == 0 {
        mode = 1;
    }
    if mode < 1 || mode > 2 {
        err!("mode must be 1 or 2");
    }

    let mut seed = 0u32;
    if getparuint("seed", &mut seed) == 0 {
        seed = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|t| t.as_secs() as u32)
            .unwrap_or(u32::MAX);
        if seed == u32::MAX {
            err!("time() failed to set seed");
        }
    }
    sranuni(seed as i32);
    if mode == 2 {
        warn!("Random seed for layers = %i", seed as i32);
    }

    let mut i1beg = 0i32;
    let mut i1end = n1 / 5;
    let mut i2beg = 0i32;
    let mut i2end = n2;
    let mut i3beg = 0i32;
    let mut i3end = n3;
    if getparint("i1beg", &mut i1beg) == 0 {
        i1beg = 0;
    }
    if getparint("i1end", &mut i1end) == 0 {
        i1end = n1 / 5;
    }
    if getparint("i2beg", &mut i2beg) == 0 {
        i2beg = 0;
    }
    if getparint("i2end", &mut i2end) == 0 {
        i2end = n2;
    }
    if getparint("i3beg", &mut i3beg) == 0 {
        i3beg = 0;
    }
    if getparint("i3end", &mut i3end) == 0 {
        i3end = n3;
    }
    let mut cvel = 0.0f32;
    if getparfloat("cvel", &mut cvel) == 0 {
        cvel = 0.0;
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
    let mut nrvl = n1 / 10;
    if getparint("nrvl", &mut nrvl) == 0 {
        nrvl = n1 / 10;
    }
    if nrvl > n1 - 1 {
        err!("nrvl cannot be greater than n1");
    }
    let mut pdv = 10.0f32;
    if getparfloat("pdv", &mut pdv) == 0 {
        pdv = 10.0;
    }
    r1 = (0.5 * r1 as f64 * r1 as f64) as f32;
    r2 = (0.5 * r2 as f64 * r2 as f64) as f32;
    r3 = (0.5 * r3 as f64 * r3 as f64) as f32;
    let iter = 2;
    let depth = 1;
    let mu = 1.0;
    let mut slowness = 0i32;
    if getparint("slowness", &mut slowness) == 0 {
        slowness = 0;
    }
    let vminc = 0.0;
    let vmaxc = 99999.0;

    // alloc3float's data are contiguous. Unread malloc storage and the never
    // initialized v2 buffer below are represented by zeroed storage, as in the
    // original's fresh allocations. Allocator-dependent garbage is undefined.
    let count = n1.wrapping_mul(n2).wrapping_mul(n3) as usize;
    let mut vel = vec![0.0f32; count];
    let invp = stdin();
    let outvp = stdout();
    efseeko(&invp, 0, SEEK_SET);
    fread(&mut vel, &invp);

    let i1 = i1beg.wrapping_add(i1end) / 2;
    let vorig = vel[i1 as usize];
    let mut vlsd = 0.0f32;
    if getparfloat("vlsd", &mut vlsd) == 0 {
        vlsd = (0.3 * vorig as f64) as f32;
    }
    let mut add = 1i32;
    if getparint("add", &mut add) == 0 {
        add = 1;
    }
    if add < 0 || add > 1 {
        err!("add must be 0 or 1");
    }
    let mut how = 0i32;
    if getparint("how", &mut how) == 0 {
        how = 0;
    }
    if how < 0 || how > 2 {
        err!("how must be 0,1, or 2");
    }

    if mode == 1 {
        let nn1 = i1end - i1beg;
        let nn2 = i2end - i2beg;
        let nn3 = i3end - i3beg;
        // Reversed bounds leave both cell loops empty. C's alloc3float can
        // return NULL for these negative extents; that buffer is never used.
        // Keep an absent layer empty instead of casting its size to usize.
        let mut v2 = if nn1 > 0 && nn2 > 0 && nn3 > 0 {
            vec![0.0f32; (nn1 * nn2 * nn3) as usize]
        } else {
            Vec::new()
        };
        // The original writes the perturbations into vel, rather than v2,
        // and starts at the origin even when the requested layer is offset.
        for i3 in 0..nn3 {
            for i2 in 0..nn2 {
                for i1 in 0..nn1 {
                    let tmp = if how == 0 {
                        (vlsd as f64 * (0.5 - franuni() as f64)) as f32
                    } else if how == 1 {
                        -vlsd * franuni()
                    } else {
                        vlsd * franuni()
                    };
                    let k = index(i3, i2, i1, n2, n1);
                    if cvel > 0.0 {
                        vel[k] = cvel;
                    } else {
                        vel[k] += tmp;
                    }
                }
            }
        }
        // Smoothing an absent layer cannot affect any velocity sample.
        if nn1 > 0 && nn2 > 0 && nn3 > 0 {
            vsm3d(
                &mut v2, nn3, nn2, nn1, iter, depth, r3, r2, r1, mu, slowness, vminc, vmaxc,
            );
        }
        for i3 in i3beg..i3end {
            let j3 = i3 - i3beg;
            for i2 in i2beg..i2end {
                let j2 = i2 - i2beg;
                for i1 in i1beg..i1end {
                    let j1 = i1 - i1beg;
                    let k = index(i3, i2, i1, n2, n1);
                    let q = index(j3, j2, j1, nn2, nn1);
                    if add == 1 {
                        if how == 1 {
                            vel[k] -= v2[q];
                        } else {
                            vel[k] += v2[q];
                        }
                    } else if add == 0 {
                        if how == 1 {
                            vel[k] = vorig - v2[q];
                        } else {
                            vel[k] = vorig + v2[q];
                        }
                    }
                    if cvel > 0.0 {
                        vel[k] = cvel;
                    }
                }
            }
        }
    }
    if mode == 2 {
        let mut i1 = 0;
        for ilay in 0..nrvl {
            let i3 = if n3 > 2 { n3 / 2 } else { 0 };
            let i2 = n2 / 2;
            // Random thicknesses can exceed n1. C keeps the original strides,
            // so those samples alias the following columns. Retain the same
            // addressing and provide backing storage for writes past the last
            // column; only the original grid is written to stdout.
            let anchor = index(i3, i2, i1, n2, n1);
            if anchor >= vel.len() {
                vel.resize(anchor + 1, 0.0);
            }
            let tmp = ((0.5 - franuni() as f64) * pdv as f64 * vel[anchor] as f64 / 100.0) as f32;
            let ih = (franuni() * n1 as f32 / (nrvl - ilay) as f32 + 1.0) as i32;
            let end = index(n3 - 1, n2 - 1, i1 + ih, n2, n1);
            if end > vel.len() {
                vel.resize(end, 0.0);
            }
            for _j in 0..ih {
                for i3 in 0..n3 {
                    for i2 in 0..n2 {
                        vel[index(i3, i2, i1, n2, n1)] += tmp;
                    }
                }
                i1 += 1;
            }
        }
    }
    fwrite(&vel[..count], &outvp);
    fclose(&invp);
    fclose(&outvp);
    exit(CWP_Exit());
}

fn index(i3: i32, i2: i32, i1: i32, n2: i32, n1: i32) -> usize {
    ((i3 * n2 + i2) * n1 + i1) as usize
}

// Static program helpers from addrvl3d.c, not library routines.
#[allow(clippy::too_many_arguments)]
fn vsm3d(
    v: &mut [f32],
    n3: i32,
    n2: i32,
    n1: i32,
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
    // The C smoother indexes d[n-2] and calls tripd2 with n-1 rows.
    // Reject singleton active dimensions as in the UB verification fix.
    if (r1 > 0.0 && n1 < 2) || (r2 > 0.0 && n2 < 2) || (r3 > 0.0 && n3 < 2) {
        err!("smoothing requires at least two samples in each active dimension");
    }
    let mut ww = 1.0f32;
    let mut w = vec![0.0f32; (n1 + n2 + n3 - 2) as usize];
    if depth == 1 {
        mu = (((mu * mu) as f64 - 1.0) / (n1 * n1) as f64) as f32;
        for i1 in 0..n1 {
            w[i1 as usize] = (1.0 / (1.0f32 + (i1 * i1) as f32 * mu) as f64) as f32;
        }
    }
    if depth == 2 {
        mu = (((mu * mu) as f64 - 1.0) / (n2 * n2) as f64) as f32;
        for i2 in 0..n2 {
            w[i2 as usize] = (1.0 / (1.0f32 + (i2 * i2) as f32 * mu) as f64) as f32;
        }
    }
    if depth == 3 {
        mu = (((mu * mu) as f64 - 1.0) / (n3 * n3) as f64) as f32;
        for i3 in 0..n3 {
            w[i3 as usize] = (1.0 / (1.0f32 + (i3 * i3) as f32 * mu) as f64) as f32;
        }
    }
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
    for i3 in 0..n3 {
        for i2 in 0..n2 {
            for i1 in 0..n1 {
                let k = index(i3, i2, i1, n2, n1);
                if v[k] > vmax {
                    v[k] = vmax;
                }
                if v[k] < vmin {
                    v[k] = vmin;
                }
            }
        }
    }
    if sl != 0 {
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                for i1 in 0..n1 {
                    let k = index(i3, i2, i1, n2, n1);
                    v[k] = (1.0 / v[k] as f64) as f32;
                }
            }
        }
    }
    if r2 > 0.0 {
        let mut d = vec![vec![0.0f32; n1 as usize]; n2 as usize];
        let mut e = d.clone();
        let mut f = d.clone();
        for i3 in 0..n3 {
            if depth == 3 {
                ww = w[i3 as usize];
            }
            for i2 in 0..n2 - 1 {
                if depth == 2 {
                    ww = w[(i2 + 1) as usize];
                }
                for i1 in 0..n1 {
                    if depth == 1 {
                        ww = w[i1 as usize];
                    }
                    d[i2 as usize][i1 as usize] = (ww as f64 + r2 as f64 * 2.0) as f32;
                    e[i2 as usize][i1 as usize] = -r2;
                    f[i2 as usize][i1 as usize] = ww * v[index(i3, i2 + 1, i1, n2, n1)];
                }
            }
            for i1 in 0..n1 {
                d[(n2 - 2) as usize][i1 as usize] -= r2;
                f[0][i1 as usize] += r2 * v[index(i3, 0, i1, n2, n1)];
            }
            tripd2(&mut d, &mut e, &mut f, n2 - 1, n1);
            for _i in 1..iter {
                for i2 in 0..n2 - 1 {
                    if depth == 2 {
                        ww = w[(i2 + 1) as usize];
                    }
                    for i1 in 0..n1 {
                        if depth == 1 {
                            ww = w[i1 as usize];
                        }
                        d[i2 as usize][i1 as usize] = (ww as f64 + r2 as f64 * 2.0) as f32;
                        e[i2 as usize][i1 as usize] = -r2;
                        f[i2 as usize][i1 as usize] *= ww;
                    }
                }
                for i1 in 0..n1 {
                    d[(n2 - 2) as usize][i1 as usize] -= r2;
                    f[0][i1 as usize] += r2 * v[index(i3, 0, i1, n2, n1)];
                }
                tripd2(&mut d, &mut e, &mut f, n2 - 1, n1);
            }
            for i2 in 0..n2 - 1 {
                for i1 in 0..n1 {
                    v[index(i3, i2 + 1, i1, n2, n1)] = f[i2 as usize][i1 as usize];
                }
            }
        }
    }
    if r3 > 0.0 {
        let mut d = vec![vec![0.0f32; n1 as usize]; n3 as usize];
        let mut e = d.clone();
        let mut f = d.clone();
        for i2 in 0..n2 {
            if depth == 2 {
                ww = w[i2 as usize];
            }
            for i3 in 0..n3 - 1 {
                if depth == 3 {
                    ww = w[(i3 + 1) as usize];
                }
                for i1 in 0..n1 {
                    if depth == 1 {
                        ww = w[i1 as usize];
                    }
                    d[i3 as usize][i1 as usize] = (ww as f64 + 2.0 * r3 as f64) as f32;
                    e[i3 as usize][i1 as usize] = -r3;
                    f[i3 as usize][i1 as usize] = ww * v[index(i3 + 1, i2, i1, n2, n1)];
                }
            }
            for i1 in 0..n1 {
                d[(n3 - 2) as usize][i1 as usize] -= r3;
                f[0][i1 as usize] += r3 * v[index(0, i2, i1, n2, n1)];
            }
            tripd2(&mut d, &mut e, &mut f, n3 - 1, n1);
            for _i in 1..iter {
                for i3 in 0..n3 - 1 {
                    if depth == 3 {
                        ww = w[(i3 + 1) as usize];
                    }
                    for i1 in 0..n1 {
                        if depth == 1 {
                            ww = w[i1 as usize];
                        }
                        d[i3 as usize][i1 as usize] = (ww as f64 + 2.0 * r3 as f64) as f32;
                        e[i3 as usize][i1 as usize] = -r3;
                        f[i3 as usize][i1 as usize] *= ww;
                    }
                }
                for i1 in 0..n1 {
                    d[(n3 - 2) as usize][i1 as usize] -= r3;
                    f[0][i1 as usize] += r3 * v[index(0, i2, i1, n2, n1)];
                }
                tripd2(&mut d, &mut e, &mut f, n3 - 1, n1);
            }
            for i3 in 0..n3 - 1 {
                for i1 in 0..n1 {
                    v[index(i3 + 1, i2, i1, n2, n1)] = f[i3 as usize][i1 as usize];
                }
            }
        }
    }
    if r1 > 0.0 {
        let mut d = vec![vec![0.0f32; 1]; n1 as usize];
        let mut e = d.clone();
        let mut f = d.clone();
        for i3 in 0..n3 {
            if depth == 3 {
                ww = w[i3 as usize];
            }
            for i2 in 0..n2 {
                if depth == 2 {
                    ww = w[i2 as usize];
                }
                for i1 in 0..n1 - 1 {
                    if depth == 1 {
                        ww = w[(i1 + 1) as usize];
                    }
                    d[i1 as usize][0] = (ww as f64 + r1 as f64 * 2.0) as f32;
                    e[i1 as usize][0] = -r1;
                    f[i1 as usize][0] = ww * v[index(i3, i2, i1 + 1, n2, n1)];
                }
                d[(n1 - 2) as usize][0] -= r1;
                f[0][0] += r1 * v[index(i3, i2, 0, n2, n1)];
                tripd2(&mut d, &mut e, &mut f, n1 - 1, 1);
                for _i in 1..iter {
                    for i1 in 0..n1 - 1 {
                        if depth == 1 {
                            ww = w[(i1 + 1) as usize];
                        }
                        d[i1 as usize][0] = (ww as f64 + r1 as f64 * 2.0) as f32;
                        e[i1 as usize][0] = -r1;
                        f[i1 as usize][0] *= ww;
                    }
                    d[(n1 - 2) as usize][0] -= r1;
                    f[0][0] += r1 * v[index(i3, i2, 0, n2, n1)];
                    tripd2(&mut d, &mut e, &mut f, n1 - 1, 1);
                }
                for i1 in 0..n1 - 1 {
                    v[index(i3, i2, i1 + 1, n2, n1)] = f[i1 as usize][0];
                }
            }
        }
    }
    if sl != 0 {
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                for i1 in 0..n1 {
                    let k = index(i3, i2, i1, n2, n1);
                    v[k] = (1.0 / v[k] as f64) as f32;
                }
            }
        }
    }
}

fn tripd2(d: &mut [Vec<f32>], e: &mut [Vec<f32>], b: &mut [Vec<f32>], n: i32, m: i32) {
    for k in 1..n as usize {
        for i in 0..m as usize {
            let temp = e[k - 1][i];
            e[k - 1][i] = temp / d[k - 1][i];
            d[k][i] -= temp * e[k - 1][i];
        }
    }
    for k in 1..n as usize {
        for i in 0..m as usize {
            b[k][i] -= e[k - 1][i] * b[k - 1][i];
        }
    }
    for i in 0..m as usize {
        b[(n - 1) as usize][i] /= d[(n - 1) as usize][i];
    }
    for k in (1..n as usize).rev() {
        for i in 0..m as usize {
            b[k - 1][i] = b[k - 1][i] / d[k - 1][i] - e[k - 1][i] * b[k][i];
        }
    }
}
