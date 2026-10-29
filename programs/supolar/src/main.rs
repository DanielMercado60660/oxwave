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

//! SUPOLAR - POLarization analysis of three-component data.
//! Port of su/main/multicomponent/supolar.c (Revision 1.9).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: Nils Maercklin, GeoForschungsZentrum (GFZ) Potsdam, Germany,
//! 1998-2001. E-mail: nils@gfz-potsdam.de.
//! Covariance modification: Mega G. Baitoff, 6 May 2015.
//! References: Jurkevics (1988), Kanasewich (1981, 1990), Meyer (1988),
//! Press, Teukolsky, Vetterling and Flannery, Numerical Recipes in C (1996),
//! Samson (1973), and Sheriff, Encyclopedic Dictionary of Exploration
//! Geophysics, 3rd ed., Society of Exploration Geophysicists (1991).

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{nint, PI};
use oxwave_symmetric_eigen::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUPOLAR - POLarization analysis of three-component data               ",
    "                                                                       ",
    " supolar <stdin [optional parameters]                                  ",
    "                                                                       ",
    " Required parameters:                                                  ",
    "    none                                                               ",
    "                                                                       ",
    " Optional parameters:                                                  ",
    "    dt=(from header)  time sampling intervall in seconds               ",
    "    wl=0.1            correlation window length in seconds             ",
    "    win=boxcar        correlation window shape, choose \"boxcar\",     ",
    "                      \"hanning\", \"bartlett\", or \"welsh\"          ",
    "    file=polar        base of output file name(s)                      ",
    "    rl=1              1 = rectilinearity evaluating 2 eigenvalues,     ",
    "                      2, 3 = rectilinearity evaluating 3 eigenvalues   ",
    "    rlq=1.0           contrast parameter for rectilinearity            ",
    "    dir=1             1 = 3 components of direction of polarization    ",
    "                          (the only three-component output file)       ",
    "    tau=0             1 = global polarization parameter                ",
    "    ellip=0           1 = principal, subprincipal, and transverse      ",
    "                          ellipticities e21, e31, and e32              ",
    "    pln=0             1 = planarity measure                            ",
    "    f1=0              1 = flatness or oblateness coefficient           ",
    "    l1=0              1 = linearity coefficient                        ",
    "    amp=0             1 = amplitude parameters: instantaneous,         ",
    "                          quadratic, and eigenresultant ir, qr, and er ",
    "    theta=0           1, 2, 3 = incidence angle of principal axis      ",
    "    phi=0             1, 2, 3 = horizontal azimuth of principal axis   ",
    "    angle=rad         unit of angles theta and phi, choose \"rad\",    ",
    "                      \"deg\", or \"gon\"                              ",
    "    all=0             1, 2, 3 = set all output flags to that value     ",
    "    verbose=0         1 = echo additional information                  ",
    "                                                                       ",
    "                                                                       ",
    " Notes:                                                                ",
    "    Three adjacent traces are considered as one three-component        ",
    "    dataset.                                                           ",
    "    Correct calculation of angles theta and phi requires the first of  ",
    "    these traces to be the vertical component, followed by the two     ",
    "    horizontal components (e.g. Z, N, E, or Z, inline, crossline).     ",
    "    Significant signal energy on Z is necessary to resolve the 180 deg ",
    "    ambiguity of phi (options phi=2,3 only).                           ",
    "                                                                       ",
    "    Each calculated polarization attribute is written into its own     ",
    "    SU file. These files get the same base name (set with \"file=\")   ",
    "    and the parameter flag as an extension (e.g. polar.rl).            ",
    "                                                                       ",
    "    In case of a tapered correlation window, the window length wl may  ",
    "    have to be increased compared to the boxcar case, because of their ",
    "    smaller effective widths (Bartlett, Hanning: 1/2, Welsh: 1/3).     ",
    "                                                                       ",
    " Range of values:                                                      ",
    "    parameter     option  interval                                     ",
    "    rl            1, 2    0.0 ... 1.0   (1.0: linear polarization)     ",
    "    rl            3      -1.0 ... 1.0                                  ",
    "    tau, l1       1       0.0 ... 1.0   (1.0: linear polarization)     ",
    "    pln, f1       1       0.0 ... 1.0   (1.0: planar polarization)     ",
    "    e21, e31, e32 1       0.0 ... 1.0   (0.0: linear polarization)     ",
    "    theta         1      -pi/2... pi/2  rad                            ",
    "    theta         2, 3    0.0 ... pi/2  rad                            ",
    "    phi           1      -pi/2... pi/2  rad                            ",
    "    phi           2      -pi  ... pi    rad   (see notes above)        ",
    "    phi           3       0.0 ... 2 pi  rad   (see notes above)        ",
    "                                                                       ",
    "                                                                       ",
];

const WBOXCAR: i32 = 0;
const WBARTLETT: i32 = 1;
const WHANNING: i32 = 2;
const WWELSH: i32 = 3;

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;

    let mut file = String::new();
    if getparstring("file", &mut file) == 0 {
        file = "polar".into();
    }
    let mut angle = String::new();
    if getparstring("angle", &mut angle) == 0 {
        angle = "rad".into();
    }
    let mut win = String::new();
    if getparstring("win", &mut win) == 0 {
        win = "boxcar".into();
    }
    let mut wl: f32 = 0.0;
    if getparfloat("wl", &mut wl) == 0 {
        wl = 0.1;
    }
    let mut dt: f32 = 0.0;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    let mut rlq: f32 = 0.0;
    if getparfloat("rlq", &mut rlq) == 0 {
        rlq = 1.0;
    }
    let mut verbose = 0;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    let mut all = 0;
    if getparint("all", &mut all) == 0 {
        all = 0;
    }
    let mut rl = 0;
    if getparint("rl", &mut rl) == 0 {
        rl = if all != 0 { all } else { 1 };
    }
    let mut dir = 0;
    if getparint("dir", &mut dir) == 0 {
        dir = 1;
    }
    let mut theta = 0;
    if getparint("theta", &mut theta) == 0 {
        theta = all;
    }
    let mut phi = 0;
    if getparint("phi", &mut phi) == 0 {
        phi = all;
    }
    let mut tau = 0;
    if getparint("tau", &mut tau) == 0 {
        tau = i32::from(all != 0);
    }
    let mut ellip = 0;
    if getparint("ellip", &mut ellip) == 0 {
        ellip = i32::from(all != 0);
    }
    let mut pln = 0;
    if getparint("pln", &mut pln) == 0 {
        pln = i32::from(all != 0);
    }
    let mut f1 = 0;
    if getparint("f1", &mut f1) == 0 {
        f1 = i32::from(all != 0);
    }
    let mut l1 = 0;
    if getparint("l1", &mut l1) == 0 {
        l1 = i32::from(all != 0);
    }
    let mut amp = 0;
    if getparint("amp", &mut amp) == 0 {
        amp = i32::from(all != 0);
    }
    checkpars();

    let iwin = match win.as_str() {
        "boxcar" => WBOXCAR,
        "bartlett" => WBARTLETT,
        "hanning" => WHANNING,
        "welsh" => WWELSH,
        _ => err!("unknown win=%s", win),
    };
    let fangle: f32 = match angle.as_str() {
        "rad" => 1.0,
        "deg" => (180.0 / PI) as f32,
        "gon" => (200.0 / PI) as f32,
        _ => err!("unknown angle=%s", angle),
    };
    if dt == 0.0 {
        dt = 0.004;
        warn!("dt not set, assuming dt=0.004");
    }
    let iwl = nint((wl / dt) as f64);
    if iwl < 1 {
        err!("wl=%g must be positive", wl);
    }
    if iwl > nt as i32 {
        err!("wl=%g too long for trace", wl);
    }
    if file.is_empty() {
        err!("file= not set and default overridden");
    }
    if verbose != 0 && (theta != 0 || phi != 0) {
        warn!("computing angles in %s", angle);
    }
    if verbose != 0 {
        warn!("%s window length = %d samples\n", win, iwl);
    }
    let iwl = iwl as usize;
    let headerfp = etmpfile();

    let rlfp = if rl != 0 {
        Some(efopen(&format!("{}.rl", file), "w"))
    } else {
        None
    };
    let thetafp = if theta != 0 {
        Some(efopen(&format!("{}.theta", file), "w"))
    } else {
        None
    };
    let phifp = if phi != 0 {
        Some(efopen(&format!("{}.phi", file), "w"))
    } else {
        None
    };
    let taufp = if tau != 0 {
        Some(efopen(&format!("{}.tau", file), "w"))
    } else {
        None
    };
    let e21fp = if ellip != 0 {
        Some(efopen(&format!("{}.e21", file), "w"))
    } else {
        None
    };
    let e31fp = if ellip != 0 {
        Some(efopen(&format!("{}.e31", file), "w"))
    } else {
        None
    };
    let e32fp = if ellip != 0 {
        Some(efopen(&format!("{}.e32", file), "w"))
    } else {
        None
    };
    let plnfp = if pln != 0 {
        Some(efopen(&format!("{}.pln", file), "w"))
    } else {
        None
    };
    let f1fp = if f1 != 0 {
        Some(efopen(&format!("{}.f1", file), "w"))
    } else {
        None
    };
    let l1fp = if l1 != 0 {
        Some(efopen(&format!("{}.l1", file), "w"))
    } else {
        None
    };
    let dirfp = if dir != 0 {
        Some(efopen(&format!("{}.dir", file), "w"))
    } else {
        None
    };
    let erfp = if amp != 0 {
        Some(efopen(&format!("{}.er", file), "w"))
    } else {
        None
    };
    let irfp = if amp != 0 {
        Some(efopen(&format!("{}.ir", file), "w"))
    } else {
        None
    };
    let qrfp = if amp != 0 {
        Some(efopen(&format!("{}.qr", file), "w"))
    } else {
        None
    };

    // alloc2float stores the three components contiguously. ampparams
    // intentionally keeps the C loop bounds, including cross-component
    // reads. Its final-component overread has undefined values in C; use
    // explicit zero padding here instead of guessing the surrounding heap.
    // Only those undefined qr samples can depend on this padding.
    let mut data3c = vec![0.0f32; 3 * nt + iwl];
    let mut a = vec![vec![0.0f32; 3]; 3];
    let mut v = vec![vec![0.0f32; 3]; 3];
    let mut d = vec![0.0f32; 3];
    let mut w = vec![0.0f32; iwl];
    calc_window(&mut w, iwl, iwin);
    let mut data_rl = if rl != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_theta = if theta != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_phi = if phi != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_tau = if tau != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_e21 = if ellip != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_e31 = if ellip != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_e32 = if ellip != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_pln = if pln != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_f1 = if f1 != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_l1 = if l1 != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_er = if amp != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_ir = if amp != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data_qr = if amp != 0 {
        vec![0.0f32; nt]
    } else {
        Vec::new()
    };
    let mut data3c_dir = if dir != 0 {
        vec![vec![0.0f32; nt]; 3]
    } else {
        Vec::new()
    };
    let mut icomp = 0;
    let mut nstat = 0;
    loop {
        efwrite(tr.hdr.as_bytes(), &headerfp);
        icomp += 1;
        data3c[(icomp - 1) * nt..icomp * nt].copy_from_slice(&tr.data[..nt]);
        if icomp == 3 {
            erewind(&headerfp);
            icomp = 0;
            nstat += 1;
            if verbose != 0 {
                fprintf!(stderr(), "%s: analyzing station %d \r", argv[0], nstat);
            }
            for it in iwl / 2..nt - iwl / 2 {
                for i in 0..3 {
                    for j in i..3 {
                        let cov =
                            covar(&data3c[i * nt..], &data3c[j * nt..], it - iwl / 2, iwl, &w);
                        a[i][j] = cov;
                        a[j][i] = cov;
                    }
                }
                eig_jacobi(&mut a, &mut d, &mut v);
                sort_eigenvalues(&mut d, &mut v);
                if rl != 0 {
                    data_rl[it] = calc_rl(&d, rlq, rl);
                }
                if theta != 0 {
                    data_theta[it] = calc_theta(&v, theta) * fangle;
                }
                if phi != 0 {
                    data_phi[it] = calc_phi(&v, phi) * fangle;
                }
                if tau != 0 {
                    data_tau[it] = calc_tau(&d);
                }
                if ellip != 0 {
                    data_e21[it] = calc_ellip(&d, 2, 1);
                    data_e31[it] = calc_ellip(&d, 3, 1);
                    data_e32[it] = calc_ellip(&d, 3, 2);
                }
                if pln != 0 {
                    data_pln[it] = calc_plan(&d);
                }
                if f1 != 0 {
                    data_f1[it] = calc_f1(&d);
                }
                if l1 != 0 {
                    data_l1[it] = calc_l1(&d);
                }
                if amp != 0 {
                    data_er[it] = calc_er(&d);
                }
                if dir != 0 {
                    calc_dir(&mut data3c_dir, &v, it);
                }
            }
            if amp != 0 {
                ampparams(&data3c, &mut data_ir, &mut data_qr, nt, iwl);
            }
            if let Some(fp) = &rlfp {
                fputdata(fp, &headerfp, &data_rl, nt, &mut tr);
            }
            if let Some(fp) = &thetafp {
                fputdata(fp, &headerfp, &data_theta, nt, &mut tr);
            }
            if let Some(fp) = &phifp {
                fputdata(fp, &headerfp, &data_phi, nt, &mut tr);
            }
            if let Some(fp) = &taufp {
                fputdata(fp, &headerfp, &data_tau, nt, &mut tr);
            }
            if let Some(fp) = &e21fp {
                fputdata(fp, &headerfp, &data_e21, nt, &mut tr);
            }
            if let Some(fp) = &e31fp {
                fputdata(fp, &headerfp, &data_e31, nt, &mut tr);
            }
            if let Some(fp) = &e32fp {
                fputdata(fp, &headerfp, &data_e32, nt, &mut tr);
            }
            if let Some(fp) = &plnfp {
                fputdata(fp, &headerfp, &data_pln, nt, &mut tr);
            }
            if let Some(fp) = &f1fp {
                fputdata(fp, &headerfp, &data_f1, nt, &mut tr);
            }
            if let Some(fp) = &l1fp {
                fputdata(fp, &headerfp, &data_l1, nt, &mut tr);
            }
            if let Some(fp) = &erfp {
                fputdata(fp, &headerfp, &data_er, nt, &mut tr);
            }
            if let Some(fp) = &irfp {
                fputdata(fp, &headerfp, &data_ir, nt, &mut tr);
            }
            if let Some(fp) = &qrfp {
                fputdata(fp, &headerfp, &data_qr, nt, &mut tr);
            }
            if let Some(fp) = &dirfp {
                fputdata3c(fp, &headerfp, &data3c_dir, nt, &mut tr);
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
    efclose(&headerfp);
    for fp in [
        &rlfp, &thetafp, &phifp, &taufp, &e21fp, &e31fp, &e32fp, &plnfp, &f1fp, &l1fp, &erfp,
        &irfp, &qrfp, &dirfp,
    ]
    .into_iter()
    .flatten()
    {
        efclose(fp);
    }
    exit(CWP_Exit());
}

// Window weights, after Sheriff (1991), p. 335; tapered windows suggested
// by Jurkevics (1988).
fn calc_window(w: &mut [f32], iwl: usize, iwin: i32) {
    let m = iwl as f32 / 2.0f32;
    let mut j = 0;
    for i in -(iwl as i32) / 2..iwl as i32 / 2 {
        w[j] = match iwin {
            WBOXCAR => 1.0,
            WBARTLETT => (1.0 - (i as f32 as f64).abs() / m as f64) as f32,
            WHANNING => (0.5 + 0.5 * (PI * (i as f32 as f64).abs() / m as f64).cos()) as f32,
            WWELSH => ((i as f32 as f64).abs() / m as f64 - 1.0).abs() as f32,
            _ => unreachable!(),
        };
        if iwin == WWELSH {
            w[j] *= w[j];
        }
        j += 1;
    }
    if j < iwl {
        w[iwl - 1] = w[0];
    }
}

fn covar(data1: &[f32], data2: &[f32], istart: usize, iwl: usize, w: &[f32]) -> f32 {
    let mut cov = 0.0f32;
    let mut mean1 = 0.0f32;
    let mut mean2 = 0.0f32;
    for (j, i) in (istart..istart + iwl).enumerate() {
        mean1 += data1[i] * w[j];
        mean2 += data2[i] * w[j];
    }
    mean1 /= iwl as f32;
    mean2 /= iwl as f32;
    for (j, i) in (istart..istart + iwl).enumerate() {
        // Modified by Mega G. Baitoff, 6 May 2015.
        cov += (data1[i] * w[j] - mean1) * (data2[i] * w[j] - mean2);
    }
    cov /= iwl as f32;
    cov
}

fn calc_rl(d: &[f32], rlq: f32, opt: i32) -> f32 {
    if d[0] != 0.0 {
        match opt {
            2 => {
                (1.0 - (0.5 * (((d[1] / d[0]) as f64).abs() + ((d[2] / d[0]) as f64).abs()))
                    .powf(rlq as f64)) as f32
            }
            3 => {
                (1.0 - (((d[1] / d[0]) as f64).abs() + ((d[2] / d[0]) as f64).abs())
                    .powf(rlq as f64)) as f32
            }
            _ => (1.0 - ((d[1] / d[0]) as f64).abs().powf(rlq as f64)) as f32,
        }
    } else {
        0.0
    }
}

fn calc_theta(v: &[Vec<f32>], opt: i32) -> f32 {
    match opt {
        2 | 3 => (v[0][0] as f64).abs().acos() as f32,
        _ => {
            if v[0][0] != 0.0 {
                let horiz = ((v[1][0] * v[1][0] + v[2][0] * v[2][0]) as f64).sqrt() as f32;
                ((horiz / v[0][0]) as f64).atan() as f32
            } else {
                0.0
            }
        }
    }
}

fn calc_phi(v: &[Vec<f32>], opt: i32) -> f32 {
    let vsign = if v[0][0] < 0.0 { -1.0f64 } else { 1.0f64 };
    match opt {
        2 | 3 => {
            let mut phi = if v[1][0] != 0.0 {
                (v[2][0] as f64 * vsign).atan2(v[1][0] as f64 * vsign) as f32
            } else if v[2][0] > 0.0 {
                (0.5 * PI * vsign) as f32
            } else {
                (-0.5 * PI * vsign) as f32
            };
            if phi < 0.0 && opt == 3 {
                phi = (phi as f64 + 2.0 * PI) as f32;
            }
            phi
        }
        _ => {
            // Preserve the C test of v[1][2].
            if v[0][1] != 0.0 {
                ((v[2][0] / v[1][0]) as f64).atan() as f32
            } else if v[2][0] > 0.0 {
                (0.5 * PI) as f32
            } else {
                (-0.5 * PI) as f32
            }
        }
    }
}

fn calc_tau(d: &[f32]) -> f32 {
    if d[0] != 0.0 {
        let x1 = (1.0 - ((d[1] / d[0]) as f64).abs()).powf(2.0) as f32;
        let x2 = (1.0 - ((d[2] / d[0]) as f64).abs()).powf(2.0) as f32;
        let x3 = (((d[1] / d[0]) as f64).abs() - ((d[2] / d[0]) as f64).abs()).powf(2.0) as f32;
        let x4 =
            (1.0 + ((d[1] / d[0]) as f64).abs() + ((d[2] / d[0]) as f64).abs()).powf(2.0) as f32;
        (((x1 + x2 + x3) / (2.0f32 * x4)) as f64).sqrt() as f32
    } else {
        0.0
    }
}

fn calc_ellip(d: &[f32], d1: usize, d2: usize) -> f32 {
    if d[d2 - 1] != 0.0 {
        ((d[d1 - 1] / d[d2 - 1]) as f64).abs().sqrt() as f32
    } else {
        0.0
    }
}

fn calc_plan(d: &[f32]) -> f32 {
    if d[0] + d[1] != 0.0 {
        (1.0 - 2.0 * d[2] as f64 / (d[0] + d[1]) as f64) as f32
    } else {
        0.0
    }
}

fn calc_f1(d: &[f32]) -> f32 {
    let x1 = (3.0 * calc_ellip(d, 3, 1) as f64) as f32;
    let x2 = (1.0 + calc_ellip(d, 2, 1) as f64 + calc_ellip(d, 3, 1) as f64) as f32;
    (1.0 - (x1 / x2) as f64) as f32
}

fn calc_l1(d: &[f32]) -> f32 {
    let x1 = (3.0 * (calc_ellip(d, 2, 1) + calc_ellip(d, 3, 1)) as f64) as f32;
    let x2 = (2.0 * (1.0f32 + calc_ellip(d, 2, 1) + calc_ellip(d, 3, 1)) as f64) as f32;
    (1.0 - (x1 / x2) as f64) as f32
}

fn calc_dir(data3c_dir: &mut [Vec<f32>], v: &[Vec<f32>], it: usize) {
    for i in 0..3 {
        data3c_dir[i][it] = v[i][0];
    }
}

fn calc_er(d: &[f32]) -> f32 {
    (d[0] as f64).abs().sqrt() as f32
}

fn ampparams(indata: &[f32], data_ir: &mut [f32], data_qr: &mut [f32], nt: usize, iwl: usize) {
    for it in 0..nt {
        data_ir[it] = ((indata[it] * indata[it]
            + indata[nt + it] * indata[nt + it]
            + indata[2 * nt + it] * indata[2 * nt + it]) as f64)
            .sqrt() as f32;
        if it >= iwl / 2 && it < nt - iwl / 2 {
            let mut sqrsum = 0.0f32;
            // The C upper bound is it+iwl, not the symmetric window end.
            // For it+iwl > nt the third component enters our explicit
            // padding; C instead reads unrelated allocations (undefined).
            for i in it - iwl / 2..it + iwl {
                sqrsum += indata[i] * indata[i]
                    + indata[nt + i] * indata[nt + i]
                    + indata[2 * nt + i] * indata[2 * nt + i];
            }
            data_qr[it] = sqrsum / iwl as f32;
        }
    }
}

fn fputdata(fileptr: &CFile, headerptr: &CFile, outdata: &[f32], nt: usize, tr: &mut Trace) {
    efread(tr.hdr.as_bytes_mut(), headerptr);
    erewind(headerptr);
    tr.data[..nt].copy_from_slice(&outdata[..nt]);
    fputtr(fileptr, tr);
}

fn fputdata3c(
    fileptr: &CFile,
    headerptr: &CFile,
    outdata3c: &[Vec<f32>],
    nt: usize,
    tr: &mut Trace,
) {
    for i in 0..3 {
        efread(tr.hdr.as_bytes_mut(), headerptr);
        tr.data[..nt].copy_from_slice(&outdata3c[i][..nt]);
        fputtr(fileptr, tr);
    }
    erewind(headerptr);
}
