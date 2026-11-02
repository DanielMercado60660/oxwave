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

//! RANDVEL3D - Add a random velocity layer to a gridded velocity model.
//! Port of par/main/velocity_model_building/randvel3d.c (Revision 1.4).
//! Original source: Copyright (c) Colorado School of Mines, 2011.
//! All rights reserved.
//! Author: U Houston: Chris Liner c. 2008.
//! Based on smooth3d (CWP: Zhenyue Liu, March 1995).

use oxwave_core::prelude::*;
use oxwave_f_franuni::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " RANDVEL3D - Add a random velocity layer (RVL) to a gridded             ",
    "            v(x,y,z) velocity model                                    ",
    "\t\t\t\t\t\t\t\t\t",
    "\trandvel3d <infile n1= n2= >outfile [parameters]\t\t\t",
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
    "   makevel nz=250 nx=200 | randvel3d n1=250 n2=200 | ximage n1=250      ",
    " 2. 3D RVL with no smoothing\t\t\t\t\t\t",
    "   makevel nz=250 nx=200 ny=220 |\t\t\t\t\t",
    "   randvel3d n1=250 n2=200 n3=220 | \t\t\t\t\t",
    "   xmovie n1=250 n2=200\t\t\t\t\t    \t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let invp = stdin();
    let outvp = stdout();
    initargs(&argv);
    requestdoc(0, SDOC);

    // Required dimensions.
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
    getparint("n3", &mut n3);
    // Reject an empty input grid before alloc3float's zero-sized pointer
    // tables and the unconditional vorig sample read (the repaired C does so).
    if n3 <= 0 {
        err!("sample number of 3rd dimension invalid");
    }
    let mut mode = 1i32;
    getparint("mode", &mut mode);
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
    getparint("i1beg", &mut i1beg);
    getparint("i1end", &mut i1end);
    getparint("i2beg", &mut i2beg);
    getparint("i2end", &mut i2end);
    getparint("i3beg", &mut i3beg);
    getparint("i3end", &mut i3end);

    let mut cvel = 0.0f32;
    getparfloat("cvel", &mut cvel);
    let mut r1 = 0.0f32;
    let mut r2 = 0.0f32;
    let mut r3 = 0.0f32;
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
    getparint("nrvl", &mut nrvl);
    if nrvl > n1 - 1 {
        err!("nrvl cannot be greater than n1");
    }
    let mut pdv = 10.0f32;
    getparfloat("pdv", &mut pdv);

    r1 = (0.5 * r1 as f64 * r1 as f64) as f32;
    r2 = (0.5 * r2 as f64 * r2 as f64) as f32;
    r3 = (0.5 * r3 as f64 * r3 as f64) as f32;
    let iter = 2;
    let depth = 1;
    let mu = 1.0f32;
    let mut slowness = 0i32;
    getparint("slowness", &mut slowness);
    let vminc = 0.0f32;
    let vmaxc = 99999.0f32;

    // alloc3float is contiguous. Fresh SU allocations are normally zero-filled
    // by the OS; C leaves any unread samples and the entire v2 uninitialized.
    let nvel = n1.wrapping_mul(n2).wrapping_mul(n3) as usize;
    let mut vel = vec![0.0f32; nvel];
    let ix = |i3: i32, i2: i32, i1: i32| ((i3 * n2 + i2) * n1 + i1) as usize;
    efseeko(&invp, 0, SEEK_SET);
    fread(&mut vel, &invp);

    let i1 = (i1beg + i1end) / 2;
    let vorig = vel[ix(0, 0, i1)];
    let mut vlsd = (0.3 * vorig as f64) as f32;
    getparfloat("vlsd", &mut vlsd);
    let mut add = 1i32;
    getparint("add", &mut add);
    if add < 0 || add > 1 {
        err!("add must be 0 or 1");
    }
    let mut how = 0i32;
    getparint("how", &mut how);
    if how < 0 || how > 2 {
        err!("how must be 0,1, or 2");
    }
    checkpars();

    if mode == 1 {
        let nn1 = i1end - i1beg;
        let nn2 = i2end - i2beg;
        let nn3 = i3end - i3beg;
        // A reversed bound makes the corresponding signed C loops empty.
        // alloc3float then fails its enormous size_t allocation; when the
        // signed loops skip the buffer, its NULL result is never dereferenced.
        let mut v2 = if nn1 < 0 || nn2 < 0 || nn3 < 0 {
            Vec::new()
        } else {
            vec![0.0f32; (nn1 * nn2 * nn3) as usize]
        };
        for i3 in 0..nn3 {
            for i2 in 0..nn2 {
                for i1 in 0..nn1 {
                    let tmp = if how == 0 {
                        (vlsd as f64 * (0.5 - franuni() as f64)) as f32
                    } else {
                        vlsd * franuni()
                    };
                    // Preserve the C source's use of vel here, rather than v2.
                    if cvel > 0.0 {
                        vel[ix(i3, i2, i1)] = cvel;
                    } else {
                        vel[ix(i3, i2, i1)] += tmp;
                    }
                }
            }
        }
        vsm3d(
            &mut v2, nn3, nn2, nn1, iter, depth, r3, r2, r1, mu, slowness, vminc, vmaxc,
        );
        for i3 in i3beg..i3end {
            let j3 = i3 - i3beg;
            for i2 in i2beg..i2end {
                let j2 = i2 - i2beg;
                for i1 in i1beg..i1end {
                    let j1 = i1 - i1beg;
                    let k = ((j3 * nn2 + j2) * nn1 + j1) as usize;
                    if add == 1 {
                        if how == 1 {
                            vel[ix(i3, i2, i1)] -= v2[k];
                        } else {
                            vel[ix(i3, i2, i1)] += v2[k];
                        }
                    } else if add == 0 {
                        if how == 1 {
                            vel[ix(i3, i2, i1)] = vorig - v2[k];
                        } else {
                            vel[ix(i3, i2, i1)] = vorig + v2[k];
                        }
                    }
                    if cvel > 0.0 {
                        vel[ix(i3, i2, i1)] = cvel;
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
            // The C layer thicknesses can carry i1 past n1. Retain the
            // resulting contiguous row aliasing, with zero-backed padding
            // beyond the allocation rather than a Rust bounds panic. Values
            // in that padding are undefined in C, like unread input samples.
            let k = ix(i3, i2, i1);
            if k >= vel.len() {
                vel.resize(k + 1, 0.0);
            }
            let tmp = ((0.5 - franuni() as f64) * pdv as f64 * vel[k] as f64 / 100.0) as f32;
            let ih = (franuni() * n1 as f32 / (nrvl - ilay) as f32 + 1.0f32) as i32;
            for _j in 0..ih {
                for i3 in 0..n3 {
                    for i2 in 0..n2 {
                        let k = ix(i3, i2, i1);
                        if k >= vel.len() {
                            vel.resize(k + 1, 0.0);
                        }
                        vel[k] += tmp;
                    }
                }
                i1 += 1;
            }
        }
    }
    fwrite(&vel[..nvel], &outvp);
    fclose(&invp);
    fclose(&outvp);
    exit(CWP_Exit());
}

/// Smooth 3-D velocity (the program's own vsm3d helper).
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
    // Keep the dimensions signed: C's loops do not run for reversed bounds.
    let ix = |i3: i32, i2: i32, i1: i32| ((i3 * n2 + i2) * n1 + i1) as usize;
    let mut ww = 1.0f32;
    // C's n1+n2+n3-2 underallocates when two layer axes are empty. Reserve
    // enough for each weight loop, as in the repaired C. For the empty layer
    // no velocity samples use these weights, so the original output is kept.
    let nw = (n1 + n2 + n3 - 2).max(n1).max(n2).max(n3);
    let mut w = if nw < 0 {
        Vec::new()
    } else {
        vec![0.0f32; nw as usize]
    };
    // Compute the weight function.
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
    // Scale smoothing parameters according to the iteration number.
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
    // Clip velocity.
    for i3 in 0..n3 {
        for i2 in 0..n2 {
            for i1 in 0..n1 {
                let k = ix(i3, i2, i1);
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
                    let k = ix(i3, i2, i1);
                    v[k] = (1.0 / v[k] as f64) as f32;
                }
            }
        }
    }

    if r2 > 0.0 {
        // Smoothing velocity in the second direction.
        let mut d = vec![0.0f32; (n1 * n2) as usize];
        let mut e = vec![0.0f32; (n1 * n2) as usize];
        let mut f = vec![0.0f32; (n1 * n2) as usize];
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
                    let k = i2 * n1 + i1;
                    d[k as usize] = (ww as f64 + r2 as f64 * 2.0) as f32;
                    e[k as usize] = -r2;
                    f[k as usize] = ww * v[ix(i3, i2 + 1, i1)];
                }
            }
            for i1 in 0..n1 {
                d[((n2 - 2) * n1 + i1) as usize] -= r2;
                f[i1 as usize] += r2 * v[ix(i3, 0, i1)];
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
                        let k = i2 * n1 + i1;
                        d[k as usize] = (ww as f64 + r2 as f64 * 2.0) as f32;
                        e[k as usize] = -r2;
                        f[k as usize] *= ww;
                    }
                }
                for i1 in 0..n1 {
                    d[((n2 - 2) * n1 + i1) as usize] -= r2;
                    f[i1 as usize] += r2 * v[ix(i3, 0, i1)];
                }
                tripd2(&mut d, &mut e, &mut f, n2 - 1, n1);
            }
            for i2 in 0..n2 - 1 {
                for i1 in 0..n1 {
                    v[ix(i3, i2 + 1, i1)] = f[(i2 * n1 + i1) as usize];
                }
            }
        }
    }
    if r3 > 0.0 {
        // Smoothing velocity in the third direction.
        let mut d = vec![0.0f32; (n1 * n3) as usize];
        let mut e = vec![0.0f32; (n1 * n3) as usize];
        let mut f = vec![0.0f32; (n1 * n3) as usize];
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
                    let k = i3 * n1 + i1;
                    d[k as usize] = (ww as f64 + 2.0 * r3 as f64) as f32;
                    e[k as usize] = -r3;
                    f[k as usize] = ww * v[ix(i3 + 1, i2, i1)];
                }
            }
            for i1 in 0..n1 {
                d[((n3 - 2) * n1 + i1) as usize] -= r3;
                f[i1 as usize] += r3 * v[ix(0, i2, i1)];
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
                        let k = i3 * n1 + i1;
                        d[k as usize] = (ww as f64 + 2.0 * r3 as f64) as f32;
                        e[k as usize] = -r3;
                        f[k as usize] *= ww;
                    }
                }
                for i1 in 0..n1 {
                    d[((n3 - 2) * n1 + i1) as usize] -= r3;
                    f[i1 as usize] += r3 * v[ix(0, i2, i1)];
                }
                tripd2(&mut d, &mut e, &mut f, n3 - 1, n1);
            }
            for i3 in 0..n3 - 1 {
                for i1 in 0..n1 {
                    v[ix(i3 + 1, i2, i1)] = f[(i3 * n1 + i1) as usize];
                }
            }
        }
    }
    if r1 > 0.0 {
        // Smoothing velocity in the first direction.
        let mut d = vec![0.0f32; n1 as usize];
        let mut e = vec![0.0f32; n1 as usize];
        let mut f = vec![0.0f32; n1 as usize];
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
                    d[i1 as usize] = (ww as f64 + r1 as f64 * 2.0) as f32;
                    e[i1 as usize] = -r1;
                    f[i1 as usize] = ww * v[ix(i3, i2, i1 + 1)];
                }
                d[(n1 - 2) as usize] -= r1;
                f[0] += r1 * v[ix(i3, i2, 0)];
                tripd2(&mut d, &mut e, &mut f, n1 - 1, 1);
                for _i in 1..iter {
                    for i1 in 0..n1 - 1 {
                        if depth == 1 {
                            ww = w[(i1 + 1) as usize];
                        }
                        d[i1 as usize] = (ww as f64 + r1 as f64 * 2.0) as f32;
                        e[i1 as usize] = -r1;
                        f[i1 as usize] *= ww;
                    }
                    d[(n1 - 2) as usize] -= r1;
                    f[0] += r1 * v[ix(i3, i2, 0)];
                    tripd2(&mut d, &mut e, &mut f, n1 - 1, 1);
                }
                for i1 in 0..n1 - 1 {
                    v[ix(i3, i2, i1 + 1)] = f[i1 as usize];
                }
            }
        }
    }
    if sl != 0 {
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                for i1 in 0..n1 {
                    let k = ix(i3, i2, i1);
                    v[k] = (1.0 / v[k] as f64) as f32;
                }
            }
        }
    }
}

/// Overwrite b with the solutions of m symmetric tridiagonal systems.
fn tripd2(d: &mut [f32], e: &mut [f32], b: &mut [f32], n: i32, m: i32) {
    // Decomposition.
    for k in 1..n {
        for i in 0..m {
            let temp = e[((k - 1) * m + i) as usize];
            e[((k - 1) * m + i) as usize] = temp / d[((k - 1) * m + i) as usize];
            d[(k * m + i) as usize] -= temp * e[((k - 1) * m + i) as usize];
        }
    }
    // Substitution.
    for k in 1..n {
        for i in 0..m {
            b[(k * m + i) as usize] -=
                e[((k - 1) * m + i) as usize] * b[((k - 1) * m + i) as usize];
        }
    }
    for i in 0..m {
        b[((n - 1) * m + i) as usize] /= d[((n - 1) * m + i) as usize];
    }
    for k in (1..n).rev() {
        for i in 0..m {
            b[((k - 1) * m + i) as usize] = b[((k - 1) * m + i) as usize]
                / d[((k - 1) * m + i) as usize]
                - e[((k - 1) * m + i) as usize] * b[(k * m + i) as usize];
        }
    }
}
