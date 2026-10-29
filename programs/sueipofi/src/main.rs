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

//! SUEIPOFI - Eigenimage (SVD) based polarization filter for three-component data.
//! Port of su/main/multicomponent/sueipofi.c (Revision 1.6).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: Nils Maercklin, GeoForschungsZentrum (GFZ) Potsdam, Germany, 2001.
//! E-mail: nils@gfz-potsdam.de
//!
//! References (from the C source):
//! Franco, R. de, and Musacchio, G., 2000: Polarization Filter with Singular
//! Value Decomposition, submitted to Geophysics, published electronically.
//! Jurkevics, A., 1988: Polarization analysis of three-component array data,
//! Bulletin of the Seismological Society of America, vol. 78, no. 5.
//! Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P.,
//! 1996: Numerical Recipes in C - The Art of Scientific Computing,
//! Cambridge University Press, Cambridge.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::cwp::nint;
use oxwave_f_cubicspline::csplin;
use oxwave_linear_interp::intlin;
use oxwave_svd::{compute_svd, svd_sort};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUEIPOFI - EIgenimage (SVD) based POlarization FIlter for             ",
    "            three-component data                                       ",
    "                                                                       ",
    " sueipofi <stdin >stdout [optional parameters]                         ",
    "                                                                       ",
    " Required parameters:                                                  ",
    "    none                                                               ",
    "                                                                       ",
    " Optional parameters:                                                  ",
    "    dt=(from header)  time sampling intervall in seconds               ",
    "    wl=0.1            SVD time window length in seconds                ",
    "    pwr=1.0           exponent of filter weights                       ",
    "    interp=cubic      interpolation between initially calculated       ",
    "                      weights, choose \"cubic\" or \"linear\"          ",
    "    verbose=0         1 = echo additional information                  ",
    "                                                                       ",
    "    file=polar        base name for additional output file(s) of       ",
    "                      filter weights (see flags below)                 ",
    "    rl1=0             1 = rectilinearity along first principal axis    ",
    "    rl2=0             1 = rectilinearity along second principal axis   ",
    "    pln=0             1 = planarity                                    ",
    "                                                                       ",
    "                                                                       ",
    " Notes:                                                                ",
    "    Three adjacent traces are considered as one three-component        ",
    "    dataset.                                                           ",
    "                                                                       ",
    "    The filter is the sum of the first two eigenimages of the singular ",
    "    value decomposition (SVD) of the signal matrix (time window).      ",
    "    Weighting functions depending on linearity and planarity of the    ",
    "    signal are applied, additionally. To avoid edge effects, these are ",
    "    interpolated linearily or via cubic splines between initially      ",
    "    calculated values of non-overlapping time windows.                 ",
    "    The algorithm is based on the assumption that the particle motion  ",
    "    trajectory is essentially 2D (elliptical polarization).            ",
    "                                                                       ",
    " Caveat:                                                               ",
    "    Cubic spline interpolation may result in filter weights exceeding  ",
    "    the set of values of initial weights. Weights outside the valid    ",
    "    interval [0.0, 1.0] are clipped.                                   ",
    "                                                                       ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get info from first trace before reading parameters, as in C.
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let mut wl = 0.1f32;
    let mut dt = (tr.dt as f64 / 1000000.0) as f32;
    let mut pwr = 1.0f32;
    let mut interp = "cubic".to_string();
    let mut file = "polar".to_string();
    let mut rl1 = 0i32;
    let mut rl2 = 0i32;
    let mut pln = 0i32;
    let mut verbose = 0i32;
    getparfloat("wl", &mut wl);
    getparfloat("dt", &mut dt);
    getparfloat("pwr", &mut pwr);
    getparstring("interp", &mut interp);
    getparstring("file", &mut file);
    getparint("rl1", &mut rl1);
    getparint("rl2", &mut rl2);
    getparint("pln", &mut pln);
    getparint("verbose", &mut verbose);
    checkpars();

    let iscubic = if interp == "cubic" {
        true
    } else if interp == "linear" {
        false
    } else {
        err!("unknown interp=%s", interp);
    };
    if dt == 0.0 {
        dt = 0.004;
        warn!("dt not set, assuming dt=0.004");
    }
    let iwl = nint((wl / dt) as f64);
    if iwl < 3 {
        err!("wl=%g must be >= %g sec. (3 samples)", wl, 3.0 * dt as f64);
    }
    if iwl > nt as i32 {
        err!("wl=%g too long for trace", wl);
    }
    if pwr < 0.0 {
        err!("pwr=%g must be positive", pwr);
    }
    let iwl = iwl as usize;
    let nw = nt / iwl;
    if verbose != 0 {
        warn!("%d windows of %d samples per trace\n", nw, iwl);
    }

    let headerfp = etmpfile();
    let rl1fp = if rl1 != 0 {
        Some(efopen(&format!("{file}.rl1"), "w"))
    } else {
        None
    };
    let rl2fp = if rl2 != 0 {
        Some(efopen(&format!("{file}.rl2"), "w"))
    } else {
        None
    };
    let plnfp = if pln != 0 {
        Some(efopen(&format!("{file}.pln"), "w"))
    } else {
        None
    };

    let mut data = vec![vec![0.0f32; nt]; 3];
    // svd_sort's supplied API explicitly requires an extra column and value
    // for the original routine's inclusive upper bound. The allocator's
    // observed zero padding is represented explicitly, without unsafe reads.
    let mut udata = vec![vec![0.0f32; 4]; iwl];
    let mut vdata = vec![vec![0.0f32; 4]; 3];
    let mut sval = vec![0.0f32; 4];
    let mut eigen = vec![vec![vec![0.0f32; nt]; 3]; 3];
    let mut polar = vec![vec![0.0f32; nt]; 3];
    // ealloc2float stores rows contiguously. do_intlin reads iny[nin]:
    // for the first two attributes this is the next attribute's first value;
    // the final extra float models the observed zero allocator padding.
    let mut initpa = vec![0.0f32; nw * 3 + 1];
    let mut icomp = 0;
    let mut nstat = 0;
    loop {
        efwrite(tr.hdr.as_bytes(), &headerfp);
        data[icomp].copy_from_slice(&tr.data[..nt]);
        icomp += 1;
        if icomp == 3 {
            erewind(&headerfp);
            icomp = 0;
            nstat += 1;
            if verbose != 0 {
                fprintf!(stderr(), "%s: processing station %d \r", argv[0], nstat);
            }
            let mut jt = 0;
            for iw in 0..nw {
                for i in 0..3 {
                    for it in 0..iwl {
                        udata[it][i] = data[i][it + jt];
                    }
                }
                compute_svd(&mut udata, &mut sval[..3], &mut vdata);
                sort_contiguous(&mut udata, &mut sval, &mut vdata);
                for j in 0..3 {
                    for i in 0..3 {
                        for it in jt..jt + iwl {
                            eigen[j][i][it] = udata[it - jt][j] * vdata[i][j] * sval[j];
                        }
                    }
                }
                for j in 0..2 {
                    initpa[j * nw + iw] = calc_svdrl(&sval, j, pwr);
                }
                initpa[2 * nw + iw] = calc_planarity(&sval, pwr);
                jt += iwl;
            }
            for i in 0..3 {
                if iscubic {
                    do_intcub(&initpa[i * nw..i * nw + nw], &mut polar[i], iwl, iwl / 2);
                    for it in 0..nt {
                        if polar[i][it] > 1.0 {
                            polar[i][it] = 1.0;
                        }
                        if polar[i][it] < 0.0 {
                            polar[i][it] = 0.0;
                        }
                    }
                } else {
                    do_intlin(
                        &initpa[i * nw..i * nw + nw + 1],
                        &mut polar[i],
                        nw,
                        iwl,
                        iwl / 2,
                    );
                }
            }
            for i in 0..3 {
                for it in 0..nw * iwl {
                    data[i][it] = eigen[0][i][it] * polar[0][it] + eigen[1][i][it] * polar[1][it];
                    data[i][it] *= polar[2][it];
                }
            }
            fputdata3c(&stdout(), &headerfp, &data, &mut tr);
            if let Some(fp) = &rl1fp {
                fputdata(fp, &headerfp, &polar[0], &mut tr);
            }
            if let Some(fp) = &rl2fp {
                fputdata(fp, &headerfp, &polar[1], &mut tr);
            }
            if let Some(fp) = &plnfp {
                fputdata(fp, &headerfp, &polar[2], &mut tr);
            }
        }
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    if verbose != 0 {
        fprintf!(stderr(), "\n");
        if icomp != 0 {
            warn!("last %d trace(s) skipped", icomp);
        }
    }
    if let Some(fp) = &rl1fp {
        efclose(fp);
    }
    if let Some(fp) = &rl2fp {
        efclose(fp);
    }
    if let Some(fp) = &plnfp {
        efclose(fp);
    }
    exit(CWP_Exit());
}

// Rectilinearity from descending singular values (Jurkevics, 1988;
// Franco and Musacchio, 2000).
fn calc_svdrl(w: &[f32], axis: usize, pwr: f32) -> f32 {
    let rl = if w[axis] * w[axis] != 0.0 {
        (1.0 - ((w[2] * w[2]) / (w[axis] * w[axis])) as f64) as f32
    } else {
        0.0f32
    };
    (rl as f64).powf(pwr as f64) as f32
}

fn calc_planarity(w: &[f32], pwr: f32) -> f32 {
    let p = if w[0] + w[1] != 0.0 {
        (1.0 - (2.0 * w[2] as f64 * w[2] as f64) / (w[0] * w[0] + w[1] * w[1]) as f64) as f32
    } else {
        0.0f32
    };
    (p as f64).powf(pwr as f64) as f32
}

// The C allocation has three contiguous columns, so the library sort's
// extra column aliases the next row's first element. Ask the supplied sort
// for the column permutation, then apply its exchanges to that layout.
// Replaying both exchanges also handles zero windows, where column 3 can
// participate in both selection steps. The last extra element represents
// the observed zero allocator padding and persists between windows.
fn sort_contiguous(u: &mut [Vec<f32>], w: &mut [f32], v: &mut [Vec<f32>]) {
    let mut tags = vec![vec![0.0f32, 1.0, 2.0, 3.0]];
    let mut probe_v = vec![vec![0.0f32; 4]; 3];
    svd_sort(&mut tags, w, &mut probe_v, 3);
    let mut order = [0usize, 1, 2, 3];
    for i in 0..2 {
        let k = (i..4).find(|&k| order[k] == tags[0][i] as usize).unwrap();
        order.swap(i, k);
        if k != i {
            exchange_contiguous(v, i, k);
            exchange_contiguous(u, i, k);
        }
    }
}

fn exchange_contiguous(matrix: &mut [Vec<f32>], i: usize, k: usize) {
    for j in 0..matrix.len() {
        if k == 3 && j + 1 < matrix.len() {
            let saved = matrix[j][i];
            matrix[j][i] = matrix[j + 1][0];
            matrix[j + 1][0] = saved;
        } else {
            matrix[j].swap(i, k);
        }
    }
}

fn do_intcub(iny: &[f32], outy: &mut [f32], din: usize, inx0: usize) {
    let nin = iny.len();
    let mut yd = vec![[0.0f32; 4]; nin];
    let mut inx = vec![0.0f32; nin];
    let mut outx = vec![0.0f32; outy.len()];

    // Initially calculated attributes are assigned to window centers.
    for ix in 0..nin {
        inx[ix] = (inx0 + ix * din) as f32;
    }
    for ix in 0..outx.len() {
        outx[ix] = ix as f32;
    }

    csplin(&inx, iny, &mut yd);
    intcub(0, &inx, &yd, &outx, outy);
}

fn do_intlin(iny: &[f32], outy: &mut [f32], nin: usize, din: usize, inx0: usize) {
    let mut inx = vec![0.0f32; nin];
    let mut outx = vec![0.0f32; outy.len()];
    for ix in 0..nin {
        inx[ix] = (inx0 + ix * din) as f32;
    }
    for ix in 0..outx.len() {
        outx[ix] = ix as f32;
    }
    intlin(&inx, &iny[..nin], iny[0], iny[nin], &outx, outy);
}

fn fputdata(fileptr: &CFile, headerptr: &CFile, outdata: &[f32], tr: &mut Trace) {
    efread(tr.hdr.as_bytes_mut(), headerptr);
    erewind(headerptr);
    tr.data[..outdata.len()].copy_from_slice(outdata);
    fputtr(fileptr, tr);
}

fn fputdata3c(fileptr: &CFile, headerptr: &CFile, outdata3c: &[Vec<f32>], tr: &mut Trace) {
    for i in 0..3 {
        efread(tr.hdr.as_bytes_mut(), headerptr);
        tr.data[..outdata3c[i].len()].copy_from_slice(&outdata3c[i]);
        fputtr(fileptr, tr);
    }
    erewind(headerptr);
}
