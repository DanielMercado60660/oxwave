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

//! SUMIGPRESP - 2-D prestack common-shot split-step Fourier migration.
//! Port of su/main/migration_inversion/sumigpresp.c (Vision 1.00).
//! Credits: CWP, Baoniu Han, bhan@dix.mines.edu, April 19th, 1998.
//! Modified: Chris Stolk, 11 Dec 2005, changed data input to remove erroneous
//! time delay; CWP, John Stockwell, 26 Sept 2006, replaced goto loops and
//! fixed honoring sx, gx, and scalco.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::{min, sgn, PI};

use oxwave_pfa_fft::{npfa, npfar, pfa2cc, pfa2rc, pfarc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGPRESP - The 2-D prestack common-shot split-step Fourier\t\t",
    "\t\tmigration \t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   sumigpresp <indata >outfile [parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " nxo=\t   number of total horizontal output samples\t\t\t",
    " nxshot=\tnumber of shot gathers to be migrated\t\t\t",
    " nz=\t    number of depth samples\t\t\t\t\t",
    " dx=\t    horizontal sampling interval\t\t\t\t",
    " dz=\t    depth sampling interval\t\t\t\t\t",
    " vfile=\t velocity profile, it must be binary format.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t   \t",
    " fmax=25\tThe peak frequency of Ricker wavelet used as source wavelet",
    " f1=5,f2=10,f3=40,f4=50\t frequencies to build a Hamming window\t",
    " lpad=9999,rpad=9999\t    number of zero traces padded on both    \t",
    "\t\t\t\tsides of depth section to determine the ",
    "\t\t\t\tmigration aperature, the default values ",
    "\t\t\t\tare using the full aperature.\t\t",
    " verbose=0             silent, =1 additional runtime information       ",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The input velocity file consists of C-style binary floats.\t    \t",
    " The structure of this file is vfile[iz][ix]. Note that this means that",
    " the x-direction is the fastest direction instead of z-direction! Such a",
    " structure is more convenient for the downward continuation type\t",
    " migration algorithm than using z as fastest dimension as in other SU  ",
    " programs.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Because most of the tools in the SU package (such as  unif2, unisam2, ",
    " and makevel) produce output with the structure vfile[ix][iz], you will",
    " need to transpose the velocity files created by these programs. You may",
    " use the SU program 'transp' in SU to transpose such files into the  ",
    " required vfile[iz][ix] structure.\t\t\t\t\t",
    " (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  \t\t",
    " denotes a v(z,x) array, the opposite of what Matlab and Fortran\t",
    " programmers may expect.)\t\t\t\t\t\t",
    "                                                                       ",
    " Also, sx must be monotonically increasing throughout the dataset, and ",
    " and gx must be monotonically increasing within a shot. You may resort ",
    " your data with 'susort', accordingly.                               ",
    "                                                                       ",
    " The scalco header field is honored so this field must be set correctly.",
    " See selfdocs of 'susort', 'suchw'. Also:   sukeyword scalco       ",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let (mut nz, mut nxo, mut nxshot) = (0i32, 0i32, 0i32);
    let (mut dz, mut dx, mut dt) = (0.0f32, 0.0f32, 0.004f32);
    let mut vfile = String::new();
    let (mut fmax, mut f1, mut f2, mut f3, mut f4) = (25.0f32, 10.0f32, 20.0f32, 40.0f32, 50.0f32);
    let (mut lpad, mut rpad, mut dip, mut verbose) = (9999i32, 9999i32, 65i32, 0i32);
    let mut oldisx = 0i32;
    let mut isx = 0i32;

    initargs(&argv);
    requestdoc(1, SDOC);

    mustgetparint("nz", &mut nz);
    mustgetparfloat("dz", &mut dz);
    mustgetparstring("vfile", &mut vfile);
    mustgetparint("nxo", &mut nxo);
    mustgetparint("nxshot", &mut nxshot);
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = 25.0;
    }
    if getparfloat("f1", &mut f1) == 0 {
        f1 = 10.0;
    }
    if getparfloat("f2", &mut f2) == 0 {
        f2 = 20.0;
    }
    if getparfloat("f3", &mut f3) == 0 {
        f3 = 40.0;
    }
    if getparfloat("f4", &mut f4) == 0 {
        f4 = 50.0;
    }
    if getparint("lpad", &mut lpad) == 0 {
        lpad = 9999;
    }
    if getparint("rpad", &mut rpad) == 0 {
        rpad = 9999;
    }
    if getparint("dip", &mut dip) == 0 {
        dip = 65;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    let _ = dip; // Parsed but otherwise unused in the original.

    // cresult[ix][iz] and vp[iz][ix], contiguous as in alloc2float.
    let mut cresult = vec![0.0f32; (nz * nxo) as usize];
    let mut vp = vec![0.0f32; (nz * nxo) as usize];
    let vfp = efopen(&vfile, "r");
    efread(&mut vp, &vfp);
    efclose(&vfp);
    cresult.fill(0.0);

    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let (mut sx, mut gx) = get_sx_gx(&tr);
    let min_sx_gx = min(sx, gx);
    let (mut gxmin, mut gxmax) = (gx, gx);
    if matches!(
        filestat(STDIN),
        FileType::Pipe | FileType::Fifo | FileType::Socket
    ) {
        // libc rewind's failed seek preserves read-ahead on a pipe and
        // clears error/EOF. Core erewind currently discards that read-ahead.
        clearerr(&stdin());
    } else {
        erewind(&stdin());
    }

    if getparfloat("dt", &mut dt) == 0 {
        if tr.dt != 0 {
            dt = (tr.dt as f64 / 1000000.0) as f32;
        } else {
            dt = 0.004;
            if verbose == 1 {
                warn!("tr.dt not set, assuming dt=0.004");
            }
        }
    }
    if getparfloat("dx", &mut dx) == 0 {
        if tr.d2 != 0.0 {
            dx = tr.d2;
        } else {
            dx = 1.0;
            if verbose == 1 {
                warn!("tr.d2 not set, assuming d2=1.0");
            }
        }
    }
    checkpars();

    loop {
        let ntfft = npfar(nt);
        let mut nw = ntfft / 2 + 1;
        let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
        let mut fw = (2.0 * PI * f1 as f64) as f32;
        let nf1 = ((fw / dw) as f64 + 0.5) as i32;
        fw = (2.0 * PI * f2 as f64) as f32;
        let nf2 = ((fw / dw) as f64 + 0.5) as i32;
        fw = (2.0 * PI * f3 as f64) as f32;
        let nf3 = ((fw / dw) as f64 + 0.5) as i32;
        fw = (2.0 * PI * f4 as f64) as f32;
        let nf4 = ((fw / dw) as f64 + 0.5) as i32;
        let truenw = nf4 - nf1 + 1;
        fw = (0.0 + (nf1 as f32 * dw) as f64) as f32;
        if verbose == 2 {
            warn!(
                "nf1=%d nf2=%d nf3=%d nf4=%d nw=%d",
                nf1, nf2, nf3, nf4, truenw
            );
        }

        let mut wl = vec![0.0f32; ntfft as usize];
        let mut wlsp = vec![Complex::default(); nw as usize];
        let wtmp = ricker(fmax, dt);
        wl.fill(0.0);
        // The original still uses the uncentered wavelet: its proposed
        // time-zero correction is inside a C comment.
        // C sumigpresp.c:255 overruns wl for short traces. Bounding this
        // copy leaves the nt=4 regression's C output byte-identical.
        for it in 0..wtmp.len().min(wl.len()) {
            wl[it] = wtmp[it];
        }
        drop(wtmp);
        pfarc(-1, &wl, &mut wlsp);

        let mut p = vec![0.0f32; (ntfft * nxo) as usize];
        let mut cq = vec![Complex::default(); (nw * nxo) as usize];
        p.fill(0.0);
        let mut nx = 0i32;
        let mut igx = if gx < 0.0 {
            (gx / dx + nxo as f32) as i32
        } else {
            (gx / dx) as i32
        };
        let mut oldigx = igx;
        let oldsx = sx;
        let oldgx = gx;
        // The C initializes oldgx once per shot, not once per receiver.
        while gettr(&mut tr) != 0 {
            (sx, gx) = get_sx_gx(&tr);
            if verbose == 2 {
                warn!("%d nx=%d", igx, nx);
            }
            igx = if gx < 0.0 {
                (gx / dx + nxo as f32) as i32
            } else {
                (gx / dx) as i32
            };
            if igx == oldigx {
                warn!("repeated igx!!! check dx or scalco value!!!");
            }
            oldigx = igx;
            // Deliberately compare the raw header to scaled oldsx, as C does.
            if tr.sx as f32 != oldsx {
                efseeko(&stdin(), (-240 - nt * 4) as i64, SEEK_CUR);
                break;
            }
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
            let start = (igx * ntfft) as usize;
            p[start..start + nt as usize].copy_from_slice(&tr.data[..nt as usize]);
            nx += 1;
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
        let il = ix2;
        let ir = ix3;
        ix2 -= lpad;
        ix3 += rpad;
        if ix2 < 0 {
            ix2 = 0;
        }
        if ix3 > nxo - 1 {
            ix3 = nxo - 1;
        }
        nx = ix3 - ix2 + 1;
        let input_nw = nw;
        nw = truenw;
        let nxfft = npfa(nx);
        let nk = nxfft;
        let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
        let fk = (-PI / dx as f64) as f32;

        let mut v = vec![vec![0.0f32; nx as usize]; nz as usize];
        for iz in 0..nz {
            for ix in 0..nx {
                v[iz as usize][ix as usize] = vp[(iz * nxo + ix + ix2) as usize];
            }
        }
        let mut cp = vec![vec![Complex::default(); nx as usize]; nw as usize];
        let mut cp1 = vec![vec![Complex::default(); nx as usize]; nw as usize];
        // Transpose cq[ix][iw] and apply the original Hamming window.
        for ix in 0..nx {
            for iw in 0..nw {
                // C sumigpresp.c:405-413 can read beyond the input spectrum.
                // Zero extension matches the minimal C bounds fix, whose
                // f4=200 regression output is identical to the original.
                let frequency = iw + nf1;
                let value = if frequency >= 0 && frequency < input_nw {
                    cq[((ix + ix2) * input_nw + frequency) as usize]
                } else {
                    cmplx(0.0, 0.0)
                };
                cp[iw as usize][ix as usize] = if iw >= 0 && iw <= nf2 - nf1 {
                    let mut tmpp = (PI / (nf2 - nf1) as f64) as f32;
                    let tmppp = ((tmpp * (iw - nf1) as f32) as f64 - PI) as f32;
                    tmpp = (0.54 + 0.46 * (tmppp as f64).cos()) as f32;
                    crmul(value, tmpp)
                } else if iw >= nf3 - nf1 && iw <= nf4 - nf1 {
                    let mut tmpp = (PI / (nf4 - nf3) as f64) as f32;
                    let tmppp = tmpp * (iw - nf3) as f32;
                    tmpp = (0.54 + 0.46 * (tmppp as f64).cos()) as f32;
                    crmul(value, tmpp)
                } else {
                    value
                };
                cp1[iw as usize][ix as usize] = cmplx(0.0, 0.0);
            }
        }
        for iw in 0..nw {
            // Apply the same bounds fix to the C wavelet read at line 420.
            let frequency = iw + nf1;
            cp1[iw as usize][(ixshot - ix2) as usize] =
                if frequency >= 0 && frequency < input_nw {
                    wlsp[frequency as usize]
                } else {
                    cmplx(0.0, 0.0)
                };
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
        let mut cq = vec![Complex::default(); (nxfft * nw) as usize];
        let mut cq1 = vec![Complex::default(); (nxfft * nw) as usize];

        for iz in 0..nz {
            // Imaging condition, before continuing either wavefield.
            for ix in 0..nx {
                let mut w = fw;
                for iw in 0..nw {
                    let ratio = 10.0f32;
                    let tmp = if (ix + ix2 - ixshot).abs() as f32 * dx < ratio * iz as f32 * dz {
                        cmul(cp[iw as usize][ix as usize], cp1[iw as usize][ix as usize])
                    } else {
                        cmplx(0.0, 0.0)
                    };
                    cresult[((ix + ix2) * nz + iz) as usize] += tmp.r / ntfft as f32;
                    w += dw;
                }
                let _ = w;
            }
            let mut vmin = 0.0f32;
            for ix in il - ix2..=ir - ix2 {
                vmin = (vmin as f64
                    + 1.0 / v[iz as usize][ix as usize] as f64 / (ir - il + 1) as f64)
                    as f32;
            }
            vmin = (1.0 / vmin as f64) as f32;

            for ik in 0..nx {
                for iw in 0..nw {
                    cq[(iw * nxfft + ik) as usize] = if ik % 2 != 0 {
                        cneg(cp[iw as usize][ik as usize])
                    } else {
                        cp[iw as usize][ik as usize]
                    };
                    cq1[(iw * nxfft + ik) as usize] = if ik % 2 != 0 {
                        cneg(cp1[iw as usize][ik as usize])
                    } else {
                        cp1[iw as usize][ik as usize]
                    };
                }
            }
            for ik in nx..nk {
                for iw in 0..nw {
                    cq[(iw * nxfft + ik) as usize] = cmplx(0.0, 0.0);
                    cq1[(iw * nxfft + ik) as usize] = cmplx(0.0, 0.0);
                }
            }
            pfa2cc(-1, 1, nk, nw, &mut cq);
            pfa2cc(-1, 1, nk, nw, &mut cq1);
            let v1 = vmin;
            let mut k = fk;
            for ik in 0..nk {
                let mut w = fw;
                for iw in 0..nw {
                    if w == 0.0 {
                        w = (1.0e-10 / dt as f64) as f32;
                    }
                    let kz1 = 1.0 - ((v1 * k / w) as f64).powf(2.0);
                    let index = (iw * nxfft + ik) as usize;
                    if kz1 > 0.15 {
                        let phase1 = -(w as f64) * kz1.sqrt() * dz as f64 / v1 as f64;
                        let cshift1 = cmplx(phase1.cos() as f32, phase1.sin() as f32);
                        cq[index] = cmul(cq[index], cshift1);
                        cq1[index] = cmul(cq1[index], cshift1);
                    } else {
                        cq[index] = cmplx(0.0, 0.0);
                        cq1[index] = cmplx(0.0, 0.0);
                    }
                    w += dw;
                }
                k += dk;
            }
            pfa2cc(1, 1, nk, nw, &mut cq);
            pfa2cc(1, 1, nk, nw, &mut cq1);
            for ix in 0..nx {
                let mut w = fw;
                for iw in 0..nw {
                    let a = 0.015f32;
                    let mut g = 1.0f32;
                    let i = 10;
                    if ix <= i {
                        g = ((-a * (i - ix) as f32 * (i - ix) as f32) as f64).exp() as f32;
                    }
                    if ix >= nx - i {
                        g = ((-a * (-nx + i + ix) as f32 * (-nx + i + ix) as f32) as f64).exp()
                            as f32;
                    }
                    let _ = g; // The original computes but never applies it.
                    let index = (iw * nxfft + ix) as usize;
                    cq[index] = crmul(cq[index], (1.0 / nxfft as f64) as f32);
                    if ix % 2 != 0 {
                        cq[index] = cneg(cq[index]);
                    }
                    let kz2 = (1.0 / v1 as f64 - 1.0 / v[iz as usize][ix as usize] as f64)
                        * w as f64
                        * dz as f64;
                    let cshift2 = cmplx(kz2.cos() as f32, kz2.sin() as f32);
                    cp[iw as usize][ix as usize] = cmul(cq[index], cshift2);
                    cq1[index] = crmul(cq1[index], (1.0 / nxfft as f64) as f32);
                    if ix % 2 != 0 {
                        cq1[index] = cneg(cq1[index]);
                    }
                    cp1[iw as usize][ix as usize] = cmul(cq1[index], cshift2);
                    w += dw;
                }
            }
        }
        drop(cp);
        drop(cp1);
        drop(cq);
        drop(cq1);
        drop(v);
        nxshot -= 1;
        if verbose == 1 {
            warn!("nsxhot %d", nxshot);
        }
        if nxshot == 0 {
            break;
        }
    }
    for ix in 0..nxo {
        tr.ns = nz as u16;
        tr.d1 = dz;
        tr.d2 = dx;
        tr.offset = 0;
        tr.cdp = ix;
        tr.tracl = ix;
        tr.data[..nz as usize]
            .copy_from_slice(&cresult[(ix * nz) as usize..((ix + 1) * nz) as usize]);
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

/// The Ricker helper belongs to sumigpresp.c, rather than libcwp.
fn ricker(mut freq: f32, mut dt: f32) -> Vec<f32> {
    if freq == 0.0 {
        freq = 30.0;
    }
    if dt == 0.0 {
        dt = 0.004;
    }
    let bpar = (6.0f64.sqrt() / (PI * freq as f64)) as f32;
    let n = (1.35 * bpar as f64 / dt as f64).ceil() as i32;
    let mut amp = vec![0.0f32; (2 * n + 1) as usize];
    amp[n as usize] = 1.0;
    for i in 1..=n {
        let t = dt * i as f32;
        let u = (2.0 * 6.0f64.sqrt() * t as f64 / bpar as f64) as f32;
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
            let scale = -(tr.scalco as i32);
            sx = tr.sx as f32 / scale as f32;
            gx = tr.gx as f32 / scale as f32;
            sy = tr.sy as f32 / scale as f32;
            gy = tr.gy as f32 / scale as f32;
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
