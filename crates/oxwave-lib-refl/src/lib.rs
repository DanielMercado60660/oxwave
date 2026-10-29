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

//! Port of Refl/lib/auxil.c, Refl/lib/compute.c, Refl/lib/convolve.c, Refl/lib/psvreflect.c, Refl/lib/shreflect.c, Refl/lib/synthetic.c: apply_earth_flattening, apply_filters, compute_Hanning_window, compute_al_at, compute_alpha_betha, compute_block_size, compute_gl_gt_gam, compute_moment_tensor, compute_p_aux_arrays, compute_prs, compute_pwin, compute_reflectivities, compute_sigmas, compute_slowness, compute_synthetics, compute_w_aux_arrays, construct_tx_trace, convolve_wavelet, parameter_interpolation, psv_reflectivities, random_velocity_layers, red_vel_factor, sh_reflectivities, source_receiver_type.
//! Credits: Colorado School of Mines; the source files list no individual authors.

#![allow(non_snake_case)]

mod bounds;
mod exact;
mod interfaces;
mod kernels;
pub use interfaces::*;

use oxwave_cwp::complex::{cdiv, cmplx, cmul, crmul, csub, Complex};

/// C prototype: `void compute_block_size(int np, int *block_size, int *nblock, int *left);`
pub fn compute_block_size(np: i32, block_size: &mut i32, nblock: &mut i32, left: &mut i32) {
    if np <= 100 {
        *block_size = np;
        *nblock = 1;
        *left = 0;
    } else {
        *block_size = 100;
        *left = np % 100;
        *nblock = np / 100 + 1;
    }
}

/// C prototype: `void source_receiver_type(int nor, int nlayers, int lsource, int *lobs, float *cl, float *ct, int *acoustic, int *flag);`
pub fn source_receiver_type(
    lsource: i32,
    lobs: &[i32],
    cl: &[f32],
    ct: &[f32],
    acoustic: &mut [i32],
    flag: &mut [i32],
) {
    let mut k = 0;
    for iz in 0..lobs.len() {
        let jkl = (lsource - 1) as usize;
        let lmn = (lobs[iz] - 1) as usize;
        acoustic[iz] = 0;
        flag[iz] = 0;
        if cl[cl.len() - 1] == cl[lmn] && ct[cl.len() - 1] == ct[lmn] {
            flag[iz] = 1;
        }
        if ct[jkl] == 0.0 {
            if ct[lmn] == 0.0 {
                acoustic[iz] = 1;
            } else {
                if lmn > jkl {
                    k += 1;
                }
                acoustic[iz] = if k == 1 { 2 } else { 3 };
            }
        } else if ct[lmn] == 0.0 {
            acoustic[iz] = 1;
        }
    }
    if lobs.len() <= 1 && lsource > lobs[0] && ct[0] == 0.0 {
        acoustic[0] = 4;
    }
}

/// C prototype: `void apply_earth_flattening(int nlayers, float z0, float *cl, float *ct, float *rho, float *t);`
pub fn apply_earth_flattening(
    z0: f32,
    cl: &mut [f32],
    ct: &mut [f32],
    rho: &mut [f32],
    t: &mut [f32],
) {
    let n = cl.len();
    let mut zz = vec![0.0; n];
    let mut z = z0;
    zz[0] = z;
    for il in 0..n - 1 {
        zz[il + 1] = z + t[il];
        z = zz[il + 1];
    }
    for il in 0..n {
        let y0 = (6371.0 / (6371.0 - zz[il] as f64)) as f32;
        zz[il] = (6371.0 * (y0 as f64).ln()) as f32;
        let w0 = (6371.0 / (6371.0 - zz[il] as f64)) as f32;
        cl[il] *= w0;
        ct[il] *= w0;
        rho[il] /= w0;
    }
    for il in 1..n {
        t[il - 1] = zz[il] - zz[il - 1];
    }
    t[n - 1] = 10000.0;
}

/// C prototype: `void compute_Hanning_window(int iwin, int iw, int if1, int if2, int nf, complex *win);`
pub fn compute_Hanning_window(iwin: i32, iw: i32, if1: i32, if2: i32, nf: i32, win: &mut Complex) {
    let pi = std::f64::consts::PI;
    let r = if iwin == 1 {
        if iw < if1 - 1 {
            (0.5 * (1.0 + (pi * (if1 - iw + 1) as f64 / if1 as f64).cos())) as f32
        } else if iw <= if2 - 1 {
            1.0
        } else if iw < nf - 1 {
            (0.5 * (1.0 + (pi * (iw - if2 + 1) as f64 / (nf - if2) as f64).cos())) as f32
        } else {
            0.0
        }
    } else if iw < nf - 1 {
        (0.5 * (1.0 + (pi * (iw + 1) as f64 / nf as f64).cos())) as f32
    } else {
        0.0
    };
    *win = Complex { r, i: 0.0 };
}

/// C prototype: `void compute_prs(int nor, int *lobs, complex *al, complex *at, float *rho, complex *prs);`
pub fn compute_prs(lobs: &[i32], al: &[Complex], at: &[Complex], rho: &[f32], prs: &mut [Complex]) {
    for ijk in 0..lobs.len() {
        let j = (lobs[ijk] - 1) as usize;
        prs[ijk] = if at[j].r == 0.0 {
            crmul(cdiv(cmplx(1.0, 0.0), cmul(al[j], al[j])), rho[j])
        } else {
            let ctemp = cmul(cdiv(al[j], at[j]), cdiv(al[j], at[j]));
            crmul(
                csub(cmplx(1.0, 0.0), crmul(ctemp, (4.0_f64 / 3.0) as f32)),
                rho[j],
            )
        };
    }
}

/// C prototype: `void compute_alpha_betha(int lsource, int block_size, complex ats, float rhos, complex *gl, complex *gt, complex *alpha, complex *betha);`
pub fn compute_alpha_betha(
    lsource: i32,
    ats: Complex,
    rhos: f32,
    gl: &[Complex],
    gt: &[Complex],
    alpha: &mut [Complex],
    betha: &mut [Complex],
) {
    let block_size = alpha.len();
    let ik1 = (lsource - 1) as usize * block_size;
    for ip in 0..block_size {
        let j = ik1 + ip;
        alpha[ip] = cdiv(cmplx(1.0, 0.0), crmul(crmul(gl[j], rhos), 2.0));
        betha[ip] = if ats.r == 0.0 && ats.i == 0.0 {
            cmplx(0.0, 0.0)
        } else {
            cdiv(cmplx(1.0, 0.0), crmul(crmul(gt[j], rhos), 2.0))
        };
    }
}

/// C prototype: `void red_vel_factor(float x, float red_vel, float tlag, complex wpie, complex wsq, complex *rvfac);`
pub fn red_vel_factor(
    x: f32,
    red_vel: f32,
    tlag: f32,
    wpie: Complex,
    wsq: Complex,
    rvfac: &mut Complex,
) {
    let lag = if red_vel == 0.0 {
        tlag
    } else {
        tlag - x / red_vel
    };
    *rvfac = cmul(
        oxwave_f_complexf::cwp_cexp1(crmul(cmul(cmplx(0.0, 1.0), wpie), lag)),
        wsq,
    );
}

/// C prototype: `void compute_moment_tensor(int wtype, float phi, float lambda, float delta, float phis, float m0, float *m1, float *m2, float *m3);`
pub fn compute_moment_tensor(
    wtype: i32,
    phi: f32,
    mut lambda: f32,
    mut delta: f32,
    mut phis: f32,
    m0: f32,
    m1: &mut f32,
    m2: &mut f32,
    m3: &mut f32,
) {
    phis -= phi;
    delta = (delta as f64 * (std::f64::consts::PI / 180.0)) as f32;
    lambda = (lambda as f64 * (std::f64::consts::PI / 180.0)) as f32;
    phis = (phis as f64 * (std::f64::consts::PI / 180.0)) as f32;
    let (p, l, d) = (phis as f64, lambda as f64, delta as f64);
    if wtype == 1 {
        let a = p.sin() as f32;
        let b = a * a;
        let c = (d.sin() * l.cos() * (2.0 * p).sin()) as f32;
        let d1 = ((2.0 * d).sin() * l.sin()) as f32;
        let e = (d.cos() * l.cos() * p.cos()) as f32;
        let f = ((2.0 * d).cos() * l.sin() * a as f64) as f32;
        *m1 = -m0 * (c + d1 * b);
        *m2 = -m0 * (e + f);
        *m3 = m0 * d1;
    } else if wtype == 2 {
        *m1 = (m0 as f64
            * (d.sin() * l.cos() * (2.0 * p).cos()
                + 0.5 * (2.0 * d).sin() * l.sin() * (2.0 * p).sin())) as f32;
        *m2 = (-m0 as f64 * (d.cos() * l.cos() * p.sin() - (2.0 * d).cos() * l.sin() * p.cos()))
            as f32;
        *m3 = 0.0;
    }
}

/// C prototype: `void compute_pwin(int wtype, int block_size, float bp, float dp, float w, float decay, float pw1, float pw2, float pw3, float pw4, complex *pwin, complex *p, complex *pp);`
pub fn compute_pwin(
    wtype: i32,
    bp: f32,
    dp: f32,
    w: f32,
    decay: f32,
    pw1: f32,
    pw2: f32,
    pw3: f32,
    pw4: f32,
    pwin: &mut [Complex],
    p: &mut [Complex],
    pp: &mut [Complex],
) {
    let mut win = 0.0f32;
    for ip in 0..pwin.len() {
        let rp = bp + ip as f32 * dp;
        if wtype == 1 {
            p[ip] = cmplx(rp, 0.0);
        } else if wtype == 2 {
            p[ip] = cmplx(rp, -decay * rp / w);
        }
        pp[ip] = cmul(p[ip], p[ip]);
        if rp >= pw1 && rp < pw2 {
            win = (0.5
                * (1.0 + (std::f64::consts::PI * (((rp - pw1) / (pw2 - pw1)) as f64 - 1.0)).cos()))
                as f32;
        }
        if rp >= pw2 && rp <= pw3 {
            win = 1.0;
        }
        if rp > pw3 && rp < pw4 {
            win = (0.5 * (1.0 + (std::f64::consts::PI * ((rp - pw3) / (pw4 - pw3)) as f64).cos()))
                as f32;
        }
        if rp < pw1 || rp >= pw4 {
            win = 0.0;
        }
        pwin[ip] = crmul(oxwave_cwp::complex::cwp_csqrt(p[ip]), win);
    }
}

/// C prototype: `void compute_gl_gt_gam(int wtype, int nlayers, int block_size, complex *al, complex *at, float *rho, complex *pp, complex *gl, complex *gt, complex *gam);`
pub fn compute_gl_gt_gam(
    wtype: i32,
    al: &[Complex],
    at: &[Complex],
    rho: &[f32],
    pp: &[Complex],
    gl: &mut [Complex],
    gt: &mut [Complex],
    gam: &mut [Complex],
) {
    use oxwave_cwp::complex::{cneg, cwp_csqrt};
    for il in 0..al.len() {
        let p1 = cmul(al[il], al[il]);
        let p2 = cmul(at[il], at[il]);
        for ip in 0..pp.len() {
            let j = il * pp.len() + ip;
            gl[j] = cwp_csqrt(csub(p1, pp[ip]));
            if p2.r != 0.0 || p2.i != 0.0 {
                gt[j] = cwp_csqrt(csub(p2, pp[ip]));
                if gl[j].i < 0.0 {
                    gl[j].i *= -1.0;
                }
                if gt[j].i < 0.0 {
                    gt[j].i *= -1.0;
                }
                if wtype == 2 {
                    if gl[j].i < 0.0 {
                        gl[j].r *= -1.0;
                    }
                    if gt[j].i < 0.0 {
                        gt[j].r *= -1.0;
                    }
                }
                gam[j] = crmul(csub(cmplx(1.0, 0.0), crmul(cdiv(pp[ip], p2), 2.0)), rho[il]);
            } else {
                gt[j] = cmplx(0.0, 0.0);
                if gl[0].i < 0.0 {
                    gl[j] = cneg(gl[j]);
                }
                gam[j] = cmplx(rho[il], 0.0);
            }
        }
    }
}

/// C prototype: `void parameter_interpolation(int nlayers, int *intlayers, int *nintlayers, float *intlayth, float *cl, float *ql, float *ct, float *qt, float *rho, float *t);`
pub fn parameter_interpolation(
    intlayers: &[i32],
    nintlayers: &[i32],
    intlayth: &[f32],
    cl: &mut [f32],
    ql: &mut [f32],
    ct: &mut [f32],
    qt: &mut [f32],
    rho: &mut [f32],
    t: &mut [f32],
) {
    let (mut i, mut ii) = (0usize, 0usize);
    let mut crust = 0.0f32;
    while i < cl.len() {
        if intlayers[ii] == i as i32 {
            let n = nintlayers[ii] as usize;
            let mut thickness = intlayth[ii];
            if thickness < 0.0 {
                thickness = -thickness - crust;
            }
            crust += thickness;
            if ct[i + n] == -1.0 {
                ct[i + n] = (cl[i + n] as f64 / 3.0f64.sqrt()) as f32;
            }
            if qt[i + n] == -1.0 {
                qt[i + n] = (ql[i + n] as f64 / 2.0) as f32;
            }
            if rho[i + n - 1] == -1.0 {
                rho[i + n] = ((cl[i + n] as f64 + 1.5) / 3.0) as f32;
            }
            let dcl = (cl[i + n] - cl[i - 1]) / (n + 1) as f32;
            let dql = (ql[i + n] - ql[i - 1]) / (n + 1) as f32;
            let dct = (ct[i + n] - ct[i - 1]) / (n + 1) as f32;
            let dqt = (qt[i + n] - qt[i - 1]) / (n + 1) as f32;
            let drho = (rho[i + n] - rho[i - 1]) / (n + 1) as f32;
            for j in 0..n {
                cl[j + i - 1] = cl[i - 1] + (j + 1) as f32 * dcl;
                ql[j + i - 1] = ql[i - 1] + (j + 1) as f32 * dql;
                ct[j + i - 1] = ct[i - 1] + (j + 1) as f32 * dct;
                qt[j + i - 1] = qt[i - 1] + (j + 1) as f32 * dqt;
                rho[j + i - 1] = rho[i - 1] + (j + 1) as f32 * drho;
                t[j + i - 1] = thickness / n as f32;
            }
            i += n;
            ii += 1;
        } else {
            if ct[i] == -1.0 {
                ct[i] = (cl[i] as f64 / 3.0f64.sqrt()) as f32;
            }
            if qt[i] == -1.0 {
                qt[i] = (ql[i] as f64 / 2.0) as f32;
            }
            if rho[i] == -1.0 {
                rho[i] = ((cl[i] as f64 + 1.5) / 3.0) as f32;
            }
            if t[i] < 0.0 {
                t[i] = -t[i] - crust;
            }
            if i < cl.len() - 1 {
                crust += t[i];
            } else {
                t[i] = 10000.0;
            }
            i += 1;
        }
    }
}

/// C prototype: `void compute_slowness(int wtype, int lsource, float w, float p2w, int *np, float *fp, float fs, float e4, float dk, float decay, float xmax, float *dp, float *pw1, float *pw2, float *pw3, float *pw4, complex *cdp, float *cl, float *ct, float *t);`
pub fn compute_slowness(
    wtype: i32,
    lsource: i32,
    w: f32,
    p2w: f32,
    np: &mut i32,
    fp: &mut f32,
    fs: f32,
    e4: f32,
    dk: f32,
    decay: f32,
    xmax: f32,
    dp: &mut f32,
    pw1: &mut f32,
    pw2: &mut f32,
    pw3: &mut f32,
    pw4: &mut f32,
    cdp: &mut Complex,
    cl: &[f32],
    ct: &[f32],
    t: &[f32],
) {
    let fq = (w as f64 / (std::f64::consts::PI * 2.0)) as f32;
    let index = if wtype == 1 {
        (lsource - 1) as usize
    } else {
        0
    };
    let velocities = if ct[index] <= 0.0005 { cl } else { ct };
    let s = (lsource - 1) as usize;
    let f = |i: usize| {
        ((1.0 / velocities[i] as f64).powf(2.0)
            + (50.0 / (w as f64 * (t[i] as f64 + 0.00001))).powf(2.0))
        .sqrt() as f32
    };
    let (fp1, fp2) = (f(s), f(s - 1));
    *fp = if fp1 > fp2 { fp1 } else { fp2 };
    let dkt = (dk as f64 * std::f64::consts::PI / xmax as f64) as f32;
    *dp = ((dkt / w) as f64 * (1.0 + (fq / fs) as f64).powf(e4 as f64 / 2.0)) as f32;
    if *dp > 0.002 {
        *dp = 0.002;
    }
    if p2w > 0.0 && *fp >= p2w {
        *fp = p2w;
    }
    *np = (*fp / *dp) as i32;
    *fp = *np as f32 * *dp;
    *cdp = cmplx(*dp, -(*dp) * (decay / w));
    if *pw2 <= 0.0 && *pw3 <= 0.0 {
        let numw = (*np as f64 * 0.15 - 1.0) as i32;
        *pw1 = 0.0;
        *pw2 = *dp * numw as f32;
        *pw3 = *fp - *dp * numw as f32;
        *pw4 = *fp;
    }
}
