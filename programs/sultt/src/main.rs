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

//! SULTT - trace by trace, sample by sample rotation of shear wave data.
//! Port of su/main/multicomponent/sultt.c (Revision 1.4).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//!
//! Credits (from the C source):
//! CWP/RCP: Rodrigo Felicio Fuck
//! Code based on algorithms presented in Li & Crampin (1993) and
//! Li & MacBeth (1997).
//! Li, X.Y., and Crampin, S., 1993, Linear-transform techniques for
//! processing shear-wave anisotropy in four-component seismic data,
//! Geophysics, 58, 240-256.
//! Li, X.Y., and MacBeth, C., 1997, Data-Matrix asymmetry and polarization
//! changes from multicomponent surface seismics.

#![allow(non_snake_case)]
#![forbid(unsafe_code)]
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;
use oxwave_f_convolution::convolve_cwp;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SULTT - trace by trace, sample by sample, rotation of shear wave data ",
    "\t  volumes using the Linear Transform Technique of Li & Crampin  ",
    "\t  (1993)\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sultt inS11=file1 inS22=file2 inS12=file3 inS21=file4 [optional       ",
    " parameters]                                                           ",
    " \t\t\t\t\t\t\t\t\t",
    " optional parameters:\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " mode\t\tdetermines what the linear transform will compute \t",
    "\t\t\tmode=1, computes asymmetry indexes\t\t",
    " \t\t\tmode=2, computes Polarization and main       \t",
    "\t\t\t\treflectivity series.\t\t\t",
    "\t\t\tmode=3, same as above, but using eigenvalues    ",
    " \t\t\t\t\t\t\t\t\t",
    " \t\tmode=3 is more robust estimation for Polarization angle ",
    "\t\tthan mode=2. In both cases the reflectivity series is   ",
    "\t\tcomputed in the same way. mode=1 outputs two other data ",
    "\t\tvolumes, each containing the asymmetry parameters theta,",
    "\t\tand gamma. The other two modes only output an extra da- ",
    "\t\tta volume, the instant polarization alpha\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " outSij\tdefines the names of the output seismic data files,     ",
    "\t\ti and j equal either 1 or 2\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " alpha, gamma\tname for optional output volumes: instanteneous polarity",
    " theta\t\t, alpha; theta, for angle misalignment between source   ",
    "\t\tand receiver; gamma, the medium asymmetric response     ",
    "\t\tcoefficient\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " wl\t\trunning window acting on traces (in samples)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " ntraces\tnumber of traces to be average for each location\t",
    " \t\t\t\t\t\t\t\t\t",
    " CAVEAT\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Naming convention for off-diagonal volumes:\t\t\t\t",
    " S12 - Inline source, Xline receiver\t\t\t\t\t",
    " S21 - Xline source, Inline receiver\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " the running will always have an odd number of samples, despite the    ",
    " input length.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

// The original globals are owned here and passed explicitly to the helpers.
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut ntraces = 5i32;
    if getparint("ntraces", &mut ntraces) == 0 {
        ntraces = 5;
    }

    let mut inS11 = String::new();
    let mut inS12 = String::new();
    let mut inS21 = String::new();
    let mut inS22 = String::new();
    if getparstring("inS11", &mut inS11) == 0 {
        err!("no file specified for S11 volume. Aborting...");
    }
    if getparstring("inS12", &mut inS12) == 0 {
        err!("no file specified for S12 volume. Aborting...");
    }
    if getparstring("inS21", &mut inS21) == 0 {
        err!("no file specified for S21 volume. Aborting...");
    }
    if getparstring("inS22", &mut inS22) == 0 {
        err!("no file specified for S22 volume. Aborting...");
    }

    // The component output statements are commented out in the C source.
    let mut outS11 = String::new();
    let mut outS12 = String::new();
    let mut outS21 = String::new();
    let mut outS22 = String::new();
    let mut alpha = String::new();
    let mut theta = String::new();
    let mut gamma = String::new();
    if getparstring("outS11", &mut outS11) == 0 {
        outS11 = "S11_out".into();
    }
    if getparstring("outS12", &mut outS12) == 0 {
        outS12 = "S12_out".into();
    }
    if getparstring("outS21", &mut outS21) == 0 {
        outS21 = "S21_out".into();
    }
    if getparstring("outS22", &mut outS22) == 0 {
        outS22 = "S22_out".into();
    }
    if getparstring("alpha", &mut alpha) == 0 {
        alpha = "alpha".into();
    }
    if getparstring("theta", &mut theta) == 0 {
        theta = "theta".into();
    }
    if getparstring("gamma", &mut gamma) == 0 {
        gamma = "gamma".into();
    }
    let _ = (&outS11, &outS12, &outS21, &outS22);

    let mut in_S11_tr = Trace::new();
    let mut in_S12_tr = Trace::new();
    let mut in_S21_tr = Trace::new();
    let mut in_S22_tr = Trace::new();
    let (iS11fp, nS11, dtS11) = open_vol(&inS11, &mut in_S11_tr);
    let (iS12, nS12, _dtS12) = open_vol(&inS12, &mut in_S12_tr);
    let (iS21, nS21, _dtS21) = open_vol(&inS21, &mut in_S21_tr);
    let (iS22, nS22, _dtS22) = open_vol(&inS22, &mut in_S22_tr);

    if nS11 != nS22 || nS11 != nS12 || nS11 != nS21 {
        err!("trace length for S11 is not equal to other traces \n");
    }
    if nS22 != nS12 || nS22 != nS21 {
        err!("trace length for S22 is not equal to other traces\n");
    }
    if nS12 != nS21 {
        err!("trace length for off diagonals is not the same\n");
    }

    // Negative sizes become huge size_t values and fail ealloc2float in C.
    if ntraces < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    if ntraces == 0 {
        // C dereferences Sij[0] after allocating zero rows. Preserve the
        // shell-visible SIGSEGV status without calling C. Direct process
        // observers can distinguish exit(139) from signal 11; reproducing
        // signal termination requires a signal helper in oxwave-core.
        exit(139);
    }
    let ntraces = ntraces as usize;
    // C uses each ealloc2float panel as one contiguous block when clearing it.
    // Fallible reservations preserve ealloc2float's diagnostic and exit code
    // for positive dimensions whose panels cannot fit in the address space.
    let mut S11 = Vec::<f32>::new();
    if S11.try_reserve_exact(nS11 * ntraces).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let mut S12 = Vec::<f32>::new();
    if S12.try_reserve_exact(nS12 * ntraces).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let mut S21 = Vec::<f32>::new();
    if S21.try_reserve_exact(nS21 * ntraces).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let mut S22 = Vec::<f32>::new();
    if S22.try_reserve_exact(nS22 * ntraces).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    S11.resize(nS11 * ntraces, 0.0);
    S12.resize(nS12 * ntraces, 0.0);
    S21.resize(nS21 * ntraces, 0.0);
    S22.resize(nS22 * ntraces, 0.0);
    S11[..nS11].copy_from_slice(&in_S11_tr.data[..nS11]);
    S12[..nS12].copy_from_slice(&in_S12_tr.data[..nS11]);
    S21[..nS21].copy_from_slice(&in_S21_tr.data[..nS11]);
    S22[..nS22].copy_from_slice(&in_S22_tr.data[..nS11]);

    let mut Total = 0i32;
    let mut wl = 3i32;
    let mut mode = 2i32;
    if getparint("Total", &mut Total) == 0 {
        Total = 0;
    }
    if getparint("wl", &mut wl) == 0 {
        wl = 3;
    }
    if getparint("mode", &mut mode) == 0 {
        mode = 2;
    }
    checkpars();

    let (mut oalpha, mut otheta, mut ogamma) = (None, None, None);
    if mode == 1 {
        otheta = Some(check_output(&theta));
        ogamma = Some(check_output(&gamma));
    } else {
        oalpha = Some(check_output(&alpha));
    }

    let mut mix = vec![0.0f32; ntraces];
    for i in 0..ntraces {
        mix[i] = 1.0f32 / ntraces as f32;
    }
    let mut mixS11 = vec![0.0f32; nS11];
    let mut mixS12 = vec![0.0f32; nS12];
    let mut mixS21 = vec![0.0f32; nS21];
    let mut mixS22 = vec![0.0f32; nS22];
    let mut Alpha = vec![0.0f32; nS11];
    let mut Theta = vec![0.0f32; nS11];
    let mut Gamma = vec![0.0f32; nS11];
    let mut alpha_tr = Trace::new();
    let mut theta_tr = Trace::new();
    let mut gamma_tr = Trace::new();
    let mut counter = 0i32;

    loop {
        if counter != 0 {
            S22[..nS22].copy_from_slice(&in_S22_tr.data[..nS22]);
        }
        Alpha.fill(0.0);
        Theta.fill(0.0);
        Gamma.fill(0.0);
        mixS11.fill(0.0);
        mixS12.fill(0.0);
        mixS21.fill(0.0);
        mixS22.fill(0.0);
        for i in 0..nS11 {
            for j in 0..ntraces {
                mixS11[i] += S11[j * nS11 + i] * mix[j];
                mixS12[i] += S12[j * nS12 + i] * mix[j];
                mixS21[i] += S21[j * nS21 + i] * mix[j];
                mixS22[i] += S22[j * nS22 + i] * mix[j];
            }
        }

        LTT(
            &mixS11, &mixS12, &mixS21, &mixS22, nS11, dtS11, &mut Alpha, &mut Theta, &mut Gamma,
            wl, mode,
        );

        if mode == 1 {
            gamma_tr.data[..nS11].copy_from_slice(&Gamma);
            theta_tr.data[..nS11].copy_from_slice(&Theta);
            gamma_tr.hdr = in_S11_tr.hdr.clone();
            theta_tr.hdr = in_S11_tr.hdr.clone();
            fputtr(otheta.as_ref().unwrap(), &mut theta_tr);
            fputtr(ogamma.as_ref().unwrap(), &mut gamma_tr);
        } else {
            alpha_tr.data[..nS11].copy_from_slice(&Alpha);
            alpha_tr.hdr = in_S11_tr.hdr.clone();
            fputtr(oalpha.as_ref().unwrap(), &mut alpha_tr);
        }

        for i in (1..ntraces).rev() {
            for j in 0..nS11 {
                S11[i * nS11 + j] = S11[(i - 1) * nS11 + j];
                S12[i * nS12 + j] = S12[(i - 1) * nS12 + j];
                S21[i * nS21 + j] = S21[(i - 1) * nS21 + j];
                S22[i * nS22 + j] = S22[(i - 1) * nS22 + j];
            }
        }
        if fgettr(&iS11fp, &mut in_S11_tr) == 0 {
            break;
        }
        if fgettr(&iS12, &mut in_S12_tr) == 0 {
            break;
        }
        if fgettr(&iS21, &mut in_S21_tr) == 0 {
            break;
        }
        S11[..nS11].copy_from_slice(&in_S11_tr.data[..nS11]);
        S12[..nS12].copy_from_slice(&in_S12_tr.data[..nS12]);
        S21[..nS21].copy_from_slice(&in_S21_tr.data[..nS21]);
        counter = counter.wrapping_add(1);
        if Total != 0 {
            let f_counter = counter as f32;
            printf!("%2.f %% done\r", 100.0f32 * (f_counter / Total as f32));
        }
        if fgettr(&iS22, &mut in_S22_tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}

#[allow(clippy::too_many_arguments)]
fn LTT(
    S11: &[f32],
    S12: &[f32],
    S21: &[f32],
    S22: &[f32],
    ns: usize,
    _dt: f32,
    Alpha: &mut [f32],
    Theta: &mut [f32],
    Gamma: &mut [f32],
    wls: i32,
    flag: i32,
) {
    let mut sumalphas = 0.0f32;
    if wls < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let mut w = vec![0.0f32; wls as usize];
    for i in 0..wls as usize {
        w[i] = 1.0;
    }

    let (mut ZETA, mut CHI, mut ZETA2, mut CHI2, mut ZETACHI) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut XI, mut ETA, mut XI2, mut ETA2, mut XIETA) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    if flag == 1 {
        ZETA = vec![0.0f32; ns];
        CHI = vec![0.0f32; ns];
        ZETA2 = vec![0.0f32; ns];
        CHI2 = vec![0.0f32; ns];
        ZETACHI = vec![0.0f32; ns];
    } else if flag > 1 {
        XI = vec![0.0f32; ns];
        ETA = vec![0.0f32; ns];
    }
    if flag == 3 {
        XI2 = vec![0.0f32; ns];
        ETA2 = vec![0.0f32; ns];
        XIETA = vec![0.0f32; ns];
    }
    for i in 0..ns {
        if flag == 1 {
            ZETA[i] = S11[i] + S22[i];
            // Opposite Sij convention to Li & Crampin (1993).
            CHI[i] = S21[i] - S12[i];
            ZETA2[i] = ZETA[i] * ZETA[i];
            CHI2[i] = CHI[i] * CHI[i];
            ZETACHI[i] = ZETA[i] * CHI[i];
        }
        if flag > 1 {
            XI[i] = S11[i] - S22[i];
            ETA[i] = S12[i] + S21[i];
        }
        if flag == 3 {
            XI2[i] = XI[i] * XI[i];
            ETA2[i] = ETA[i] * ETA[i];
            XIETA[i] = XI[i] * ETA[i];
        }
    }
    if flag == 1 {
        sliding_window(ns, &mut ZETA, wls, &w);
        sliding_window(ns, &mut CHI, wls, &w);
        sliding_window(ns, &mut ZETA2, wls, &w);
        sliding_window(ns, &mut CHI2, wls, &w);
        sliding_window(ns, &mut ZETACHI, wls, &w);
    } else if flag > 1 {
        sliding_window(ns, &mut XI, wls, &w);
        sliding_window(ns, &mut ETA, wls, &w);
    }
    if flag == 3 {
        sliding_window(ns, &mut XI2, wls, &w);
        sliding_window(ns, &mut ETA2, wls, &w);
        sliding_window(ns, &mut XIETA, wls, &w);
    }
    for i in 0..ns {
        if flag == 1 {
            let (diffalphas, ratio) = Asym_analysis(ZETA2[i], ZETACHI[i], CHI2[i]);
            Gamma[i] = ratio;
            Theta[i] = diffalphas;
        } else if flag == 2 {
            sumalphas = (((ETA[i] / XI[i]) as f64).atan() * 180.0 / PI) as f32;
        } else if flag == 3 {
            sumalphas = eig_2x2(XI2[i], XIETA[i], ETA2[i]);
        }
        if flag > 1 {
            Alpha[i] = (0.5 * sumalphas as f64) as f32;
        }
    }
}

#[allow(dead_code)]
fn AngleAv(angle1: f32, angle2: f32) -> f32 {
    let x1 = (angle1 as f64 * PI / 90.0).cos() as f32;
    let y1 = (angle1 as f64 * PI / 90.0).sin() as f32;
    let x2 = (angle2 as f64 * PI / 90.0).cos() as f32;
    let y2 = (angle2 as f64 * PI / 90.0).sin() as f32;
    let xav = (x1 + x2) / 2.0f32;
    let yav = (y1 + y2) / 2.0f32;
    (((yav / xav) as f64).atan() * 90.0 / PI) as f32
}

#[allow(dead_code)]
fn TimeSeries(x: f32, y: f32) -> f32 {
    if x != 0.0 {
        let signx = (x as f64 / (x as f64).abs()) as f32;
        (signx as f64 * ((x * x + y * y) as f64).sqrt()) as f32
    } else {
        0.0
    }
}

fn Asym_analysis(a11: f32, a12: f32, a22: f32) -> (f32, f32) {
    let threshold = 1.0e-10f32;
    let mut D = (((a11 - a22) * (a11 - a22) - 4.0f32 * a12 * a12) as f64).sqrt() as f32;
    let major = (0.5 * (a11 + a22 + D) as f64) as f32;
    let minor = (0.5 * (a11 + a22 - D) as f64) as f32;
    if major <= threshold {
        D = 0.0;
    } else {
        D = minor / major;
    }
    let angle = (0.5 * ((2.0f32 * a12) as f64).atan2((a11 - a22) as f64)) as f32;
    // Multiplication by the integer 180 is performed in float before / PI.
    let phi = ((angle * 180.0f32) as f64 / PI) as f32;
    (phi, D)
}

fn eig_2x2(a11: f32, a12: f32, a22: f32) -> f32 {
    let angle = (0.5 * ((2.0f32 * a12) as f64).atan2((a11 - a22) as f64)) as f32;
    ((angle * 180.0f32) as f64 / PI) as f32
}

fn sliding_window(ns: usize, vector: &mut [f32], nsfilter: i32, filter: &[f32]) {
    let ntemp = (ns as i32).wrapping_add(nsfilter).wrapping_sub(1);
    let mut temp = vec![0.0f32; ntemp as usize];
    let bound = nsfilter / 2;
    convolve_cwp(0, &vector[..ns], -bound, filter, 0, &mut temp);
    // Preserve the original shifted window and even-length end behavior.
    for i in bound..(ntemp - bound) {
        vector[(i - bound) as usize] = temp[i as usize];
    }
}

fn open_vol(vol_name: &str, tr: &mut Trace) -> (CFile, usize, f32) {
    let fp = match fopen(vol_name, "r") {
        Some(fp) => fp,
        None => err!("cannot open volume %s\n", vol_name),
    };
    let mut tr2 = Trace::new();
    if fgettr(&fp, &mut tr2) == 0 {
        err!("cannot read first trace in %s volume\n", vol_name);
    }
    let nt = tr2.ns as usize;
    let dt = (tr2.dt as f64 / 1000000.0) as f32;
    // C checks the non-null dt pointer rather than *dt; dt=0 is accepted.
    *tr = tr2;
    (fp, nt, dt)
}

fn check_output(filename: &str) -> CFile {
    match fopen(filename, "w") {
        Some(fp) => fp,
        None => err!("cannot open output volume %s file \n", filename),
    }
}
