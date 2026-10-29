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

//! SUMIGPREPSPI - 2-D prestack common-shot phase-shift-plus-interpolation migration.
//! Port of su/main/migration_inversion/sumigprepspi.c (Vision 1.00).
//! Credits: CWP, Baoniu Han, bhan@dix.mines.edu, April 19th, 1998.
//! Modified: Chris Stolk, 11 Dec 2005, removed erroneous time delay.
//! Modified: CWP, John Stockwell, 26 Sept 2006, replaced goto loops and
//! fixed honoring sx, gx, and scalco.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::{min, sgn, PI};
use oxwave_pfa_fft::{npfa, npfar, pfa2cc, pfa2rc, pfarc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGPREPSPI --- The 2-D PREstack commom-shot Phase-Shift-Plus \t",
    "\t\t\tinterpolation depth MIGration.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   sumigprepspi <indata >outfile [parameters] \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t   \t",
    "\t\t\t\t\t\t\t\t\t",
    " nxo=     number of total horizontal output samples\t\t\t",
    " nxshot=  number of shot gathers to be migrated\t\t   \t",
    " nz=      number of depth samples\t\t\t\t \t",
    " dx=      horizontal sampling interval\t\t\t  \t",
    " dz=      depth sampling interval\t\t\t\t \t",
    " vfile=   velocity profile, it must be binary format.                 ",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t   \t",
    " fmax=25    the peak frequency of Ricker wavelet used as source wavelet",
    " f1=5,f2=10,f3=40,f4=50     frequencies to build a Hamming window     ",
    " lpad=9999,rpad=9999        number of zero traces padded on both\t",
    "                            sides of depth section to determine the\t",
    "                            migration aperture, the default values    ",
    "                            are using the full aperture.              ",
    " nflag=0    normalization of cross-correlation:                       ",
    "            0: none, 1: by source wave field                          ",
    " verbose=0  silent, =1 additional runtime information\t                ",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The input velocity file 'vfile' consists of C-style binary floats.\t",
    " The structure of this file is vfile[iz][ix]. Note that this means that",
    " the x-direction is the fastest direction instead of z-direction! Such a",
    " structure is more convenient for the downward continuation type\t",
    " migration algorithm than using z as fastest dimension as in other SU  ",
    " programs.\t\t\t\t\t\t   \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Because most of the tools in the SU package (such as  unif2, unisam2, ",
    " and makevel) produce output with the structure vfile[ix][iz], you will",
    " need to transpose the velocity files created by these programs. You may",
    " use the SU program 'transp' in SU to transpose such files into the  ",
    " required vfile[iz][ix] structure.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  \t\t",
    " denotes a v(z,x) array, the opposite of what Matlab and Fortran\t",
    " programmers may expect.)\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Also, sx must be monotonically increasing throughout the dataset, and ",
    " and gx must be monotonically increasing within a shot. You may resort ",
    " your data with 'susort', accordingly.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The scalco header field is honored so this field must be set correctly.",
    " See selfdocs of 'susort', 'suchw'. Also:   sukeyword scalco\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let (mut nz, mut dz, mut vfile, mut nxo, mut nxshot) = (0i32, 0f32, String::new(), 0i32, 0i32);
    mustgetparint("nz", &mut nz);
    mustgetparfloat("dz", &mut dz);
    mustgetparstring("vfile", &mut vfile);
    mustgetparint("nxo", &mut nxo);
    mustgetparint("nxshot", &mut nxshot);
    let (mut fmax, mut f1, mut f2, mut f3, mut f4) = (25f32, 10f32, 20f32, 40f32, 50f32);
    getparfloat("fmax", &mut fmax);
    getparfloat("f1", &mut f1);
    getparfloat("f2", &mut f2);
    getparfloat("f3", &mut f3);
    getparfloat("f4", &mut f4);
    let (mut lpad, mut rpad, mut dip, mut nflag, mut verbose) = (9999, 9999, 65, 0, 0);
    getparint("lpad", &mut lpad);
    getparint("rpad", &mut rpad);
    getparint("dip", &mut dip);
    getparint("nflag", &mut nflag);
    getparint("verbose", &mut verbose);

    let mut cresult = vec![0f32; (nxo * nz) as usize];
    let mut vp = vec![0f32; (nxo * nz) as usize];
    let vfp = efopen(&vfile, "r");
    efread(&mut vp, &vfp);
    efclose(&vfp);

    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let (mut sx, mut gx) = get_sx_gx(&tr);
    let min_sx_gx = min(sx, gx);
    sx -= min_sx_gx;
    gx -= min_sx_gx;
    let (mut dt, mut dx) = (0.004f32, 0f32);
    if getparfloat("dt", &mut dt) == 0 {
        if tr.dt != 0 {
            dt = (tr.dt as f64 / 1000000.0) as f32;
        } else {
            dt = 0.004;
            warn!("tr.dt not set, assuming dt=0.004");
        }
    }
    if getparfloat("dx", &mut dx) == 0 {
        if tr.d2 != 0.0 {
            dx = tr.d2;
        } else {
            dx = 1.0;
            warn!("tr.d2 not set, assuming d2=1.0");
        }
    }
    checkpars();

    let (mut gxmin, mut gxmax) = (0f32, 0f32);
    let (mut oldigx, mut oldisx, mut isx) = (0i32, 0i32, 0i32);
    loop {
        // do ... while(nxshot)
        let ntfft = npfar(nt);
        let mut nw = ntfft / 2 + 1;
        let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
        let fw = (2.0 * PI * f1 as f64) as f32;
        let nf1 = ((fw / dw) as f64 + 0.5) as i32;
        let fw = (2.0 * PI * f2 as f64) as f32;
        let nf2 = ((fw / dw) as f64 + 0.5) as i32;
        let fw = (2.0 * PI * f3 as f64) as f32;
        let nf3 = ((fw / dw) as f64 + 0.5) as i32;
        let fw = (2.0 * PI * f4 as f64) as f32;
        let nf4 = ((fw / dw) as f64 + 0.5) as i32;
        let truenw = nf4 - nf1 + 1;
        let fw = (0.0 + (nf1 as f32 * dw) as f64) as f32;
        if verbose == 2 {
            warn!(
                "nf1=%d nf2=%d nf3=%d nf4=%d nw=%d",
                nf1, nf2, nf3, nf4, truenw
            );
        }
        let mut wl = vec![0f32; ntfft as usize];
        let mut wlsp = vec![Complex::default(); nw as usize];
        let wtmp = ricker(fmax, dt);
        let ntw = wtmp.len() as i32;
        for it in 0..ntw {
            wl[((it - ntw / 2 + ntfft) % ntfft) as usize] = wtmp[it as usize];
        }
        pfarc(-1, &wl, &mut wlsp);

        let mut p = vec![0f32; (ntfft * nxo) as usize];
        let mut cq = vec![Complex::default(); (nw * nxo) as usize];
        let mut nx = 0;
        let mut igx = 0i32;
        let oldsx = sx;
        let oldgx = gx;
        loop {
            // copy before updating igx, exactly as in the C source
            let off = (igx * ntfft) as usize;
            p[off..off + nt as usize].copy_from_slice(&tr.data[..nt as usize]);
            (sx, gx) = get_sx_gx(&tr);
            sx -= min_sx_gx;
            gx -= min_sx_gx;
            igx = (gx / dx) as i32;
            if igx == oldigx {
                warn!("repeated igx!!! check dx or scalco value!!!");
            }
            oldigx = igx;
            if gxmin > gx {
                gxmin = gx;
            }
            if gxmax < gx {
                gxmax = gx;
            }
            if verbose == 2 {
                warn!(
                    " inside loop:  min_sx_gx %f isx %d igx %d gx %f sx %f",
                    min_sx_gx, isx, igx, gx, sx
                );
            }
            if !(oldsx <= sx) {
                err!("sx field must be monotonically increasing!");
            }
            if !(oldgx <= gx) {
                err!("gx field must be monotonically increasing!");
            }
            nx += 1;
            if gettr(&mut tr) == 0 || sx != oldsx {
                break;
            }
        }
        isx = (oldsx / dx) as i32;
        if isx == oldisx {
            warn!("repeated isx!!! check dx or scalco value!!!");
        }
        oldisx = isx;
        let ixshot = isx;
        if verbose == 2 {
            warn!(
                "sx %f, gx %f , gxmin %f  gxmax %f nx %d",
                sx, gx, gxmin, gxmax, nx
            );
            warn!("isx %d igx %d ixshot %d", isx, igx, ixshot);
        }
        pfa2rc(1, 1, ntfft, nxo, &p, &mut cq);
        let ix1 = (oldsx / dx) as i32;
        let mut ix2 = (gxmin / dx) as i32;
        let mut ix3 = (gxmax / dx) as i32;
        if ix1 >= ix3 {
            ix3 = ix1;
        }
        if ix1 <= ix2 {
            ix2 = ix1;
        }
        ix2 -= lpad;
        ix3 += rpad;
        if ix2 < 0 {
            ix2 = 0;
        }
        if ix3 > nxo - 1 {
            ix3 = nxo - 1;
        }
        nx = ix3 - ix2 + 1;
        let fullnw = nw;
        nw = truenw;
        let nxfft = npfa(nx);
        let nk = nxfft;
        let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
        let fk = (-PI / dx as f64) as f32;
        let mut v = vec![0f32; (nz * nx) as usize];
        for iz in 0..nz {
            for ix in 0..nx {
                v[(iz * nx + ix) as usize] = vp[(iz * nxo + ix + ix2) as usize];
            }
        }
        let mut cp = vec![Complex::default(); (nx * nw) as usize];
        let mut cp1 = vec![Complex::default(); (nx * nw) as usize];
        for ix in 0..nx {
            for iw in 0..nw {
                let q = cq[((ix + ix2) * fullnw + iw + nf1) as usize];
                cp[(iw * nx + ix) as usize] = if iw >= 0 && iw <= nf2 - nf1 {
                    let mut tmpp = (PI / (nf2 - nf1) as f64) as f32;
                    let tmppp = ((tmpp * (iw - nf1) as f32) as f64 - PI) as f32;
                    tmpp = (0.54 + 0.46 * (tmppp as f64).cos()) as f32;
                    crmul(q, tmpp)
                } else if iw >= nf3 - nf1 && iw <= nf4 - nf1 {
                    let mut tmpp = (PI / (nf4 - nf3) as f64) as f32;
                    let tmppp = tmpp * (iw - nf3) as f32;
                    tmpp = (0.54 + 0.46 * (tmppp as f64).cos()) as f32;
                    crmul(q, tmpp)
                } else {
                    q
                };
                cp1[(iw * nx + ix) as usize] = cmplx(0.0, 0.0);
            }
        }
        for iw in 0..nw {
            cp1[(iw * nx + ixshot - ix2) as usize] = wlsp[(iw + nf1) as usize];
        }
        if verbose == 2 {
            warn!(
                "ixshot %d ix %d ix1 %d ix2 %d ix3 %d",
                ixshot, nx, ix1, ix2, ix3
            );
            warn!("oldsx %f ", oldsx);
        }
        drop(p);
        drop(cq);
        drop(wl);
        drop(wlsp);
        let mut cq = vec![Complex::default(); (nk * nw) as usize];
        let mut cq1 = vec![Complex::default(); (nk * nw) as usize];
        for iz in 0..nz {
            let mut eps = 0f32;
            if nflag != 0 {
                for ix in 0..nx {
                    for iw in 0..nw {
                        let rtmp = rcabs(cp1[(iw * nx + ix) as usize]);
                        eps = if eps > rtmp { eps } else { rtmp };
                    }
                }
                eps = (eps as f64 * 10.0) as f32;
                eps = eps * eps;
            }
            for ix in 0..nx {
                let mut _w = fw;
                for iw in 0..nw {
                    let (tmp, tmpd);
                    let j = (iw * nx + ix) as usize;
                    if ((ix + ix2 - ixshot).abs() as f32 * dx) < 10f32 * iz as f32 * dz {
                        tmp = cmul(cp[j], cp1[j]);
                        tmpd = if nflag != 0 {
                            cmul(cp1[j], conjg(cp1[j]))
                        } else {
                            Complex::default()
                        };
                    } else {
                        tmp = cmplx(0.0, 0.0);
                        tmpd = cmplx(1.0, 0.0);
                    }
                    if nflag != 0 {
                        cresult[((ix + ix2) * nz + iz) as usize] +=
                            tmp.r / (tmpd.r + eps) / ntfft as f32;
                    } else {
                        cresult[((ix + ix2) * nz + iz) as usize] += tmp.r / ntfft as f32;
                    }
                    _w += dw;
                }
            }
            for ik in 0..nx {
                let mut w = fw;
                for iw in 0..nw {
                    let j = (iw * nx + ik) as usize;
                    cp[j] = cmul(
                        cp[j],
                        cwp_cexp(cmplx(0.0, -w * dz / v[(iz * nx + ik) as usize])),
                    );
                    cp1[j] = cmul(
                        cp1[j],
                        cwp_cexp(cmplx(0.0, -w * dz / v[(iz * nx + ik) as usize])),
                    );
                    cq[(iw * nk + ik) as usize] = if ik % 2 != 0 { cneg(cp[j]) } else { cp[j] };
                    cq1[(iw * nk + ik) as usize] = if ik % 2 != 0 { cneg(cp1[j]) } else { cp1[j] };
                    w += dw;
                }
            }
            for ik in nx..nk {
                for iw in 0..nw {
                    cq[(iw * nk + ik) as usize] = cmplx(0.0, 0.0);
                    cq1[(iw * nk + ik) as usize] = cmplx(0.0, 0.0);
                }
            }
            pfa2cc(-1, 1, nk, nw, &mut cq);
            pfa2cc(-1, 1, nk, nw, &mut cq1);
            let (mut vmax, mut vmin) = (v[(iz * nx) as usize], v[(iz * nx) as usize]);
            for ix in 0..nx {
                if v[(iz * nx + ix) as usize] >= vmax {
                    vmax = v[(iz * nx + ix) as usize];
                }
                if v[(iz * nx + ix) as usize] <= vmin {
                    vmin = v[(iz * nx + ix) as usize];
                }
            }
            let v1 = vmin;
            let v2 = vmax;
            if (((v2 - v1) / v1) as f64) < 0.01 {
                let mut k = fk;
                for ik in 0..nk {
                    let mut w = fw;
                    for iw in 0..nw {
                        if w == 0.0 {
                            w = (1.0e-10 / dt as f64) as f32;
                        }
                        let kz1 = 1.0 - ((v1 * k / w) as f64).powf(2.0);
                        let cshift1 = if kz1 >= 0.0 {
                            let phase1 = -(w as f64) * kz1.sqrt() * dz as f64 / v1 as f64
                                + (w * dz / v1) as f64;
                            cmplx(phase1.cos() as f32, phase1.sin() as f32)
                        } else {
                            let phase1 = -(w as f64) * (-kz1).sqrt() * dz as f64 / v1 as f64;
                            cwp_cexp(cmplx(phase1 as f32, w * dz / v1))
                        };
                        let j = (iw * nk + ik) as usize;
                        cq[j] = cmul(cq[j], cshift1);
                        cq1[j] = cmul(cq1[j], cshift1);
                        w += dw;
                    }
                    k += dk;
                }
                pfa2cc(1, 1, nk, nw, &mut cq);
                pfa2cc(1, 1, nk, nw, &mut cq1);
                for ix in 0..nx {
                    for iw in 0..nw {
                        let q = (iw * nk + ix) as usize;
                        let j = (iw * nx + ix) as usize;
                        cq[q] = crmul(cq[q], (1.0 / nxfft as f64) as f32);
                        cp[j] = if ix % 2 != 0 { cneg(cq[q]) } else { cq[q] };
                        cq1[q] = crmul(cq1[q], (1.0 / nxfft as f64) as f32);
                        cp1[j] = if ix % 2 != 0 { cneg(cq1[q]) } else { cq1[q] };
                    }
                }
            } else {
                // Adaptive reference velocities from the C program's own histogram.
                const L: usize = 10;
                let mut c = [0f32; 101];
                let mut p = [0f32; 100];
                let mut y = [0f32; 101];
                for ik in 0..=L {
                    c[ik] = (vmin as f64 + ik as f64 * 1.0 * (vmax - vmin) as f64 / L as f32 as f64)
                        as f32;
                }
                for ix in 0..nx {
                    for ik in 0..L {
                        let vel = v[(iz * nx + ix) as usize];
                        if (vel >= c[ik] && vel < c[ik + 1]) || (ik == L - 1 && vel == vmax) {
                            p[ik] = (p[ik] as f64 + 1.0 / nx as f64) as f32;
                            break;
                        }
                    }
                }
                let mut sz = 0f32;
                for ik in 0..L {
                    if p[ik] != 0.0 {
                        sz = (sz as f64 - p[ik] as f64 * (p[ik] as f64).ln()) as f32;
                    }
                }
                let bz = ((sz as f64).exp() + 0.5) as usize;
                y[0] = 0.0;
                y[L] = 1.0;
                for ik in 1..L {
                    y[ik] = 0.0;
                    for ix in 0..ik {
                        y[ik] = y[ik] + p[ix];
                    }
                }
                let mut reference = vec![0f32; bz + 1];
                reference[0] = vmin;
                for ix in 1..=bz {
                    for ik in 0..L {
                        if ix as f64 * 1.0 / bz as f64 > y[ik] as f64
                            && ix as f64 * 1.0 / bz as f64 <= y[ik + 1] as f64
                        {
                            reference[ix] = (c[ik] as f64
                                + (ix as f64 * 1.0 / bz as f64 - y[ik] as f64)
                                    * (c[ik + 1] - c[ik]) as f64
                                    / (y[ik + 1] - y[ik]) as f64)
                                as f32;
                            break;
                        }
                    }
                }
                reference[bz] = (reference[bz] as f64 * 1.005) as f32;
                let mut cq2 = vec![vec![Complex::default(); (nk * nw) as usize]; bz + 1];
                let mut cq3 = vec![vec![Complex::default(); (nk * nw) as usize]; bz + 1];
                for ix in 0..=bz {
                    let mut w = fw;
                    for iw in 0..nw {
                        let mut k = fk;
                        for ik in 0..nk {
                            if w == 0.0 {
                                w = (1.0e-10 / dt as f64) as f32;
                            }
                            let kz1 = 1.0 - ((reference[ix] * k / w) as f64).powf(2.0);
                            let cshift1 = if kz1 >= 0.0 {
                                let phase1 = -(w as f64) * kz1.sqrt() * dz as f64
                                    / reference[ix] as f64
                                    + (w * dz / reference[ix]) as f64;
                                cwp_cexp(cmplx(0.0, phase1 as f32))
                            } else {
                                let phase1 =
                                    -(w as f64) * (-kz1).sqrt() * dz as f64 / reference[ix] as f64;
                                cwp_cexp(cmplx(phase1 as f32, w * dz / reference[ix]))
                            };
                            let j = (iw * nk + ik) as usize;
                            cq2[ix][j] = cmul(cq[j], cshift1);
                            cq3[ix][j] = cmul(cq1[j], cshift1);
                            k += dk;
                        }
                        w += dw;
                    }
                    pfa2cc(1, 1, nk, nw, &mut cq2[ix]);
                    pfa2cc(1, 1, nk, nw, &mut cq3[ix]);
                    for ik in 0..nx {
                        for iw in 0..nw {
                            // The original computes a taper, then unconditionally sets g=1.
                            let g = 1f32;
                            let j = (iw * nk + ik) as usize;
                            cq2[ix][j] = crmul(cq2[ix][j], (g as f64 * 1.0 / nxfft as f64) as f32);
                            if ik % 2 != 0 {
                                cq2[ix][j] = cneg(cq2[ix][j]);
                            }
                            cq3[ix][j] = crmul(cq3[ix][j], (g as f64 * 1.0 / nxfft as f64) as f32);
                            if ik % 2 != 0 {
                                cq3[ix][j] = cneg(cq3[ix][j]);
                            }
                        }
                    }
                }
                for ix in 0..nx {
                    for ik in 0..bz {
                        let vel = v[(iz * nx + ix) as usize];
                        if vel >= reference[ik] && vel < reference[ik + 1] {
                            let v1 = reference[ik];
                            let v2 = reference[ik + 1];
                            for iw in 0..nw {
                                let j = (iw * nk + ix) as usize;
                                let (a1, a2) = (cq2[ik][j].r as f64, cq2[ik + 1][j].r as f64);
                                let (theta1, theta2) =
                                    (cq2[ik][j].i as f64, cq2[ik + 1][j].i as f64);
                                let a = a1 * (v2 - vel) as f64 / (v2 - v1) as f64
                                    + a2 * (vel - v1) as f64 / (v2 - v1) as f64;
                                let theta = theta1 * (v2 - vel) as f64 / (v2 - v1) as f64
                                    + theta2 * (vel - v1) as f64 / (v2 - v1) as f64;
                                cp[(iw * nx + ix) as usize] = cmplx(a as f32, theta as f32);
                                let (a1, a2) = (cq3[ik][j].r as f64, cq3[ik + 1][j].r as f64);
                                let (theta1, theta2) =
                                    (cq3[ik][j].i as f64, cq3[ik + 1][j].i as f64);
                                let a = a1 * (v2 - vel) as f64 / (v2 - v1) as f64
                                    + a2 * (vel - v1) as f64 / (v2 - v1) as f64;
                                let theta = theta1 * (v2 - vel) as f64 / (v2 - v1) as f64
                                    + theta2 * (vel - v1) as f64 / (v2 - v1) as f64;
                                cp1[(iw * nx + ix) as usize] = cmplx(a as f32, theta as f32);
                            }
                            break;
                        }
                    }
                }
            }
        }
        nxshot -= 1;
        if verbose == 1 {
            warn!("nxshot %d", nxshot);
        }
        if nxshot == 0 {
            break;
        }
    }
    for ix in 0..nxo {
        tro.ns = nz as u16;
        tro.d1 = dz;
        tro.d2 = dx;
        tro.offset = 0;
        tro.trid = 130;
        tro.cdp = ix;
        tro.tracl = ix;
        tro.tracr = ix;
        tro.data[..nz as usize]
            .copy_from_slice(&cresult[(ix * nz) as usize..((ix + 1) * nz) as usize]);
        puttr(&mut tro);
    }
    exit(CWP_Exit());
}

fn ricker(mut freq: f32, mut dt: f32) -> Vec<f32> {
    if freq == 0.0 {
        freq = 30.0;
    }
    if dt == 0.0 {
        dt = 0.004;
    }
    let bpar = (6f64.sqrt() / (PI * freq as f64)) as f32;
    let n = (1.35 * bpar as f64 / dt as f64).ceil() as i32;
    let mut amp = vec![0f32; (2 * n + 1) as usize];
    amp[n as usize] = 1.0;
    for i in 1..=n {
        let t = dt * i as f32;
        let u = (2.0 * 6f64.sqrt() * t as f64 / bpar as f64) as f32;
        let value = (0.5 * (2.0 - (u * u) as f64) * ((-u * u) as f64 / 4.0).exp()) as f32;
        amp[(n + i) as usize] = value;
        amp[(n - i) as usize] = value;
    }
    amp
}

fn get_sx_gx(tr: &Trace) -> (f32, f32) {
    let (mut sx, mut gx, sy, gy);
    if tr.scalco != 0 {
        if tr.scalco > 0 {
            sx = tr.sx as f32 * tr.scalco as f32;
            gx = tr.gx as f32 * tr.scalco as f32;
            sy = tr.sy as f32 * tr.scalco as f32;
            gy = tr.gy as f32 * tr.scalco as f32;
        } else {
            let scale = (tr.scalco as i32).abs() as f32;
            sx = tr.sx as f32 / scale;
            gx = tr.gx as f32 / scale;
            sy = tr.sy as f32 / scale;
            gy = tr.gy as f32 / scale;
        }
    } else {
        sx = tr.sx as f32;
        gx = tr.gx as f32;
        sy = tr.sy as f32;
        gy = tr.gy as f32;
    }
    sx = (sgn(sx - sy) * ((sx * sx + sy * sy) as f64).sqrt()) as f32;
    gx = (sgn(gx - gy) * ((gx * gx + gy * gy) as f64).sqrt()) as f32;
    (sx, gx)
}
