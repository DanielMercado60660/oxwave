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

//! SUMIGPREFD - 2-D prestack common-shot finite-difference depth migration.
//! Port of su/main/migration_inversion/sumigprefd.c (Vision 1.00, 2015/08/07).
//! Credits: CWP, Baoniu Han, bhan@dix.mines.edu, April 19th, 1998
//! Modified: Chris Stolk, 11 Dec 2005, - changed data input
//! to remove erroneous time delay.
//! Modified: CWP, John Stockwell 26 Sept 2006 - replaced Han's
//! "goto-loop" in two places with "do { }while loops".
//! Fixed it so that sx, gx, and scalco are honored.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::{nint, sgn, PI};
use oxwave_pfa_fft::{npfar, pfa2rc, pfarc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGPREFD --- The 2-D prestack common-shot 45-90 degree\t\t",
    "\t\t\tfinite-difference depth migration. \t\t",
    "\t\t\t\t\t\t\t\t\t",
    "    sumigprefd <indata >outfile [parameters] \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " nxo=\t\tnumber of total horizontal output samples\t\t",
    " nxshot=\tnumber of shot gathers to be migrated\t\t\t",
    " nz=\t\tnumber of depth samples\t\t\t\t\t",
    " dx=\t\thorizontal sampling interval\t\t\t\t",
    " dz=\t\tdepth sampling interval\t\t\t\t \t",
    " vfile=\tvelocity profile, it must be binary format (see Notes)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dip=79\tthe maximum dip to migrate, possible values are:\t",
    "\t\t45,65,79,80,87,89,90 degrees\t\t\t\t",
    "\t\tThe computation cost is 45=65=79<80<87<89<90\t\t",
    " fmax=25\tpeak frequency of Ricker wavelet used as source wavelet\t",
    " f1=5,f2=10,f3=40,f4=50\t frequencies to build a Hamming window\t",
    "\t\t\t\t\t\t\t\t\t",
    " lpad=9999,rpad=9999\tnumber of zero traces padded on both\t\t",
    "\t\t\tsides of depth section to determine the\t\t",
    "\t\t\tmigration aperature, the default values \t",
    "\t\t\tare using the full aperature.\t\t\t",
    " verbose=0\t\tsilent, =1 additional runtime information\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The input velocity file 'vfile' consists of C-style binary floats.  ",
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
    "\t\t\t\t\t\t\t\t\t",
    " Also, sx must be monotonically increasing throughout the dataset, and ",
    " and gx must be monotonically increasing within a shot. You may resort\t",
    " your data with 'susort', accordingly.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The scalco header field is honored so this field must be set correctly.",
    " See selfdocs of 'susort', 'suchw'. Also:   sukeyword scalco\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut nz, mut nxo, mut nxshot) = (0i32, 0i32, 0i32);
    let mut dz = 0.0f32;
    let mut vfile = String::new();
    mustgetparint("nz", &mut nz);
    mustgetparint("nxo", &mut nxo);
    mustgetparfloat("dz", &mut dz);
    mustgetparstring("vfile", &mut vfile);
    mustgetparint("nxshot", &mut nxshot);
    let (mut fmax, mut f1, mut f2, mut f3, mut f4) = (25.0f32, 10.0f32, 20.0f32, 40.0f32, 50.0f32);
    getparfloat("fmax", &mut fmax);
    getparfloat("f1", &mut f1);
    getparfloat("f2", &mut f2);
    getparfloat("f3", &mut f3);
    getparfloat("f4", &mut f4);
    let (mut lpad, mut rpad, mut dip, mut verbose) = (9999i32, 9999i32, 79i32, 0i32);
    getparint("lpad", &mut lpad);
    getparint("rpad", &mut rpad);
    getparint("dip", &mut dip);
    getparint("verbose", &mut verbose);
    // These sizes otherwise produce invalid C allocations or overrun tr.data.
    if nz < 0 || nz as usize > SU_NFLTS || nxo <= 0 {
        err!("invalid migration dimensions: nz=%d nxo=%d", nz, nxo);
    }
    let nvelocity = nxo.checked_mul(nz).unwrap_or_else(|| {
        err!("velocity dimensions exceed integer range");
    });
    let mut cresult = vec![vec![0.0f32; nz as usize]; nxo as usize];
    let mut vp = vec![0.0f32; nvelocity as usize];
    let vfp = efopen(&vfile, "r");
    efread(&mut vp, &vfp);
    efclose(&vfp);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let (mut sx, mut gx) = get_sx_gx(&tr);
    let min_sx_gx = if sx < gx { sx } else { gx };
    sx -= min_sx_gx;
    gx -= min_sx_gx;
    let (mut dt, mut dx) = (0.004f32, 0.0f32);
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
    let (mut gxmin, mut gxmax) = (0.0f32, 0.0f32);
    let mut isx = 0i32;
    let mut oldisx = 0i32;
    loop {
        let ntfft = npfar(nt as i32) as usize;
        let nw = ntfft / 2 + 1;
        let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
        let mut fw = (2.0 * PI * f1 as f64) as f32;
        let nf1 = ((fw / dw) as f64 + 0.5) as i32;
        fw = (2.0 * PI * f2 as f64) as f32;
        let nf2 = ((fw / dw) as f64 + 0.5) as i32;
        fw = (2.0 * PI * f3 as f64) as f32;
        let nf3 = ((fw / dw) as f64 + 0.5) as i32;
        fw = (2.0 * PI * f4 as f64) as f32;
        let nf4 = ((fw / dw) as f64 + 0.5) as i32;
        // Only nf1..=nf4 index cq and wlsp. Check the rounded bins, so
        // sub-bin negative corners and the Nyquist bin retain C behavior.
        if nf1 < 0 || nf4 < nf1 || nf4 as usize >= nw {
            err!(
                "frequency window must select FFT bins between 0 and %d",
                nw - 1
            );
        }
        let truenw = nf4 - nf1 + 1;
        fw = (0.0 + (nf1 as f32 * dw) as f64) as f32;
        if verbose == 2 {
            warn!(
                "nf1=%d nf2=%d nf3=%d nf4=%d nw=%d",
                nf1, nf2, nf3, nf4, truenw
            );
        }
        let mut wl = vec![0.0f32; ntfft];
        let mut wlsp = vec![Complex::default(); nw];
        let wtmp = ricker(fmax, dt, ntfft);
        let ntw = wtmp.len() as i32;
        // Stolk's centered source wavelet. ricker checks that adding ntfft
        // makes every centered index nonnegative, as required by C's remainder.
        for it in 0..ntw {
            wl[((it - ntw / 2 + ntfft as i32) % ntfft as i32) as usize] = wtmp[it as usize];
        }
        pfarc(-1, &wl, &mut wlsp);
        let mut p = vec![0.0f32; ntfft * nxo as usize];
        let mut cq = vec![Complex::default(); nw * nxo as usize];
        let mut nx = 0i32;
        let mut igx = 0i32;
        let mut oldigx = 0i32;
        let oldsx = sx;
        let oldgx = gx;
        // Preserve the C's delayed receiver placement and shot-header updates.
        loop {
            if igx < 0 || igx >= nxo {
                err!("receiver index %d outside output section", igx);
            }
            let start = igx as usize * ntfft;
            p[start..start + nt].copy_from_slice(&tr.data[..nt]);
            (sx, gx) = get_sx_gx(&tr);
            sx -= min_sx_gx;
            gx -= min_sx_gx;
            igx = nint((gx / dx) as f64);
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
        isx = nint((oldsx / dx) as f64);
        let ixshot = isx;
        if isx == oldisx {
            warn!("repeated isx!!! check dx or scalco value!!!");
        }
        oldisx = isx;
        if verbose == 2 {
            warn!(
                "sx %f, gx %f , gxmin %f  gxmax %f nx %d",
                sx, gx, gxmin, gxmax, nx
            );
            warn!("isx %d igx %d ixshot %d", isx, igx, ixshot);
        }
        pfa2rc(1, 1, ntfft as i32, nxo, &p, &mut cq);
        let ix1 = nint((oldsx / dx) as f64);
        let mut ix2 = nint((gxmin / dx) as f64);
        let mut ix3 = nint((gxmax / dx) as f64);
        if ix1 >= ix3 {
            ix3 = ix1;
        }
        if ix1 <= ix2 {
            ix2 = ix1;
        }
        ix2 = ix2.checked_sub(lpad).unwrap_or_else(|| {
            err!("migration aperture index exceeds integer range");
        });
        ix3 = ix3.checked_add(rpad).unwrap_or_else(|| {
            err!("migration aperture index exceeds integer range");
        });
        if ix2 < 0 {
            ix2 = 0;
        }
        if ix3 > nxo - 1 {
            ix3 = nxo - 1;
        }
        if ix2 > ix3 || ixshot < ix2 || ixshot > ix3 {
            err!("migration aperture must contain the source within the output section");
        }
        nx = ix3 - ix2 + 1;
        // fdmig reads data[2] and data[nx-3] at each frequency. C accesses
        // outside its allocations for narrower apertures; report an error.
        if nz > 0 && nx < 3 {
            err!("migration aperture must contain at least 3 positions");
        }
        let nw = truenw as usize;
        let nx = nx as usize;
        let mut v = vec![vec![0.0f32; nx]; nz as usize];
        for iz in 0..nz as usize {
            for ix in 0..nx {
                v[iz][ix] = vp[iz * nxo as usize + ix + ix2 as usize];
            }
        }
        let mut cp = vec![vec![Complex::default(); nx]; nw];
        let mut cp1 = vec![vec![Complex::default(); nx]; nw];
        for ix in 0..nx {
            for iw in 0..nw {
                let i = iw as i32;
                let mut tmpp;
                let tmppp;
                if i >= 0 && i <= nf2 - nf1 {
                    tmpp = (PI / (nf2 - nf1) as f64) as f32;
                    tmppp = ((tmpp * (i - nf1) as f32) as f64 - PI) as f32;
                    tmpp = (0.54 + 0.46 * (tmppp as f64).cos()) as f32;
                    cp[iw][ix] = crmul(
                        cq[(ix + ix2 as usize) * (ntfft / 2 + 1) + (i + nf1) as usize],
                        tmpp,
                    );
                } else if i >= nf3 - nf1 && i <= nf4 - nf1 {
                    tmpp = (PI / (nf4 - nf3) as f64) as f32;
                    tmppp = tmpp * (i - nf3) as f32;
                    tmpp = (0.54 + 0.46 * (tmppp as f64).cos()) as f32;
                    cp[iw][ix] = crmul(
                        cq[(ix + ix2 as usize) * (ntfft / 2 + 1) + (i + nf1) as usize],
                        tmpp,
                    );
                } else {
                    cp[iw][ix] = cq[(ix + ix2 as usize) * (ntfft / 2 + 1) + (i + nf1) as usize];
                }
                cp1[iw][ix] = cmplx(0.0, 0.0);
            }
        }
        for iw in 0..nw {
            cp1[iw][(ixshot - ix2) as usize] = wlsp[iw + nf1 as usize];
        }
        if verbose == 2 {
            warn!(
                "ixshot %d ix %d ix1 %d ix2 %d ix3 %d",
                ixshot, nx as i32, ix1, ix2, ix3
            );
            warn!("oldsx %f ", oldsx);
        }
        for iz in 0..nz as usize {
            // Imaging condition: ordered frequency accumulation in single precision.
            for ix in 0..nx {
                let mut _w = fw;
                for iw in 0..nw {
                    let ratio = 10.0f32;
                    let tmp =
                        if (ix as i32 + ix2 - ixshot).abs() as f32 * dx < ratio * iz as f32 * dz {
                            cmul(cp[iw][ix], cp1[iw][ix])
                        } else {
                            cmplx(0.0, 0.0)
                        };
                    cresult[ix + ix2 as usize][iz] += tmp.r / ntfft as f32;
                    _w += dw;
                }
            }
            let mut v1 = 0.0f32;
            for ix in 0..nx {
                v1 += v[iz][ix] / nx as f32;
            }
            for ix in 0..nx {
                let mut w = fw;
                for iw in 0..nw {
                    let kz2 = -(1.0 / v1 as f64) * w as f64 * dz as f64;
                    let cshift2 = cmplx(kz2.cos() as f32, kz2.sin() as f32);
                    cp[iw][ix] = cmul(cp[iw][ix], cshift2);
                    cp1[iw][ix] = cmul(cp1[iw][ix], cshift2);
                    w += dw;
                }
            }
            fdmig(&mut cp, nx, nw, &v[iz], fw, dw, dz, dx, dt, dip);
            fdmig(&mut cp1, nx, nw, &v[iz], fw, dw, dz, dx, dt, dip);
            for ix in 0..nx {
                let mut w = fw;
                for iw in 0..nw {
                    let kz2 = -(1.0 / v[iz][ix] as f64 - 1.0 / v1 as f64) * w as f64 * dz as f64;
                    let cshift2 = cmplx(kz2.cos() as f32, kz2.sin() as f32);
                    cp[iw][ix] = cmul(cp[iw][ix], cshift2);
                    cp1[iw][ix] = cmul(cp1[iw][ix], cshift2);
                    w += dw;
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
    for ix in 0..nxo as usize {
        tr.ns = nz as u16;
        tr.d1 = dz;
        tr.d2 = dx;
        tr.offset = 0;
        tr.cdp = ix as i32;
        tr.tracl = ix as i32;
        tr.data[..nz as usize].copy_from_slice(&cresult[ix]);
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

fn ricker(mut freq: f32, mut dt: f32, ntfft: usize) -> Vec<f32> {
    if freq == 0.0 {
        freq = 30.0;
    }
    if dt == 0.0 {
        dt = 0.004;
    }
    let bpar = (6.0f64.sqrt() / (PI * freq as f64)) as f32;
    let n = (1.35 * bpar as f64 / dt as f64).ceil() as i32;
    // A half-length above ntfft makes the first centered C remainder
    // negative. Validate before allocating, also rejecting negative lengths.
    if n < 0 || n as usize > ntfft {
        err!("Ricker wavelet length is invalid for FFT length %d", ntfft);
    }
    let mut amp = vec![0.0f32; (2 * n + 1) as usize];
    amp[n as usize] = 1.0;
    for i in 1..=n {
        let t = dt * i as f32;
        let u = (2.0 * 6.0f64.sqrt() * t as f64 / bpar as f64) as f32;
        let value = (0.5 * (2.0 - (u * u) as f64) * (((-u * u) as f64) / 4.0).exp()) as f32;
        amp[(n + i) as usize] = value;
        amp[(n - i) as usize] = value;
    }
    amp
}

#[allow(clippy::too_many_arguments)]
fn fdmig(
    cp: &mut [Vec<Complex>],
    nx: usize,
    nw: usize,
    v: &[f32],
    fw: f32,
    dw: f32,
    dz: f32,
    dx: f32,
    dt: f32,
    dip: i32,
) {
    if !matches!(dip, 45 | 65 | 79 | 80 | 87 | 89 | 90) {
        err!("unsupported migration dip %d", dip);
    }
    let mut s1 = vec![0.0f32; nx];
    let mut s2 = vec![0.0f32; nx];
    let mut data = vec![Complex::default(); nx];
    let mut d = vec![Complex::default(); nx];
    let mut a = vec![Complex::default(); nx];
    let mut b = vec![Complex::default(); nx];
    let mut c = vec![Complex::default(); nx];
    // Only the seven documented dips have defined C coefficients.
    let mut coefa = [0.0f32; 5];
    let mut coefb = [0.0f32; 5];
    let mut step = 1;
    if dip == 45 {
        coefa[0] = 0.5;
        coefb[0] = 0.25;
        step = 1;
    }
    if dip == 65 {
        coefa[0] = 0.478242060;
        coefb[0] = 0.376369527;
        step = 1;
    }
    if dip == 79 {
        coefa[0] = 0.4575;
        coefb[0] = 0.4575;
        step = 1;
    }
    if dip == 80 {
        coefa[1] = 0.040315157;
        coefb[1] = 0.873981642;
        coefa[0] = 0.457289566;
        coefb[0] = 0.222691983;
        step = 2;
    }
    if dip == 87 {
        coefa[2] = 0.00421042;
        coefb[2] = 0.972926132;
        coefa[1] = 0.081312882;
        coefb[1] = 0.744418059;
        coefa[0] = 0.414236605;
        coefb[0] = 0.150843924;
        step = 3;
    }
    if dip == 89 {
        coefa[3] = 0.000523275;
        coefb[3] = 0.994065088;
        coefa[2] = 0.014853510;
        coefb[2] = 0.919432661;
        coefa[1] = 0.117592008;
        coefb[1] = 0.614520676;
        coefa[0] = 0.367013245;
        coefb[0] = 0.105756624;
        step = 4;
    }
    if dip == 90 {
        coefa[4] = 0.000153427;
        coefb[4] = 0.997370236;
        coefa[3] = 0.004172967;
        coefb[3] = 0.964827992;
        coefa[2] = 0.033860918;
        coefb[2] = 0.824918565;
        coefa[1] = 0.143798076;
        coefb[1] = 0.483340757;
        coefa[0] = 0.318013812;
        coefb[0] = 0.073588213;
        step = 5;
    }
    let v1 = v[0];
    let vn = v[nx - 1];
    let trick = 0.1f32;
    loop {
        step -= 1;
        let mut w = fw;
        for iw in 0..nw {
            if (w as f64).abs() <= 1.0e-10 {
                w = (1.0e-10 / dt as f64) as f32;
            }
            for ix in 0..nx {
                s1[ix] = (v[ix] * v[ix]) * coefb[step] / (dx * dx * w * w) + trick;
                s2[ix] = ((-v[ix] * dz * coefa[step] / (w * dx * dx)) as f64 * 0.5) as f32;
            }
            for ix in 0..nx {
                data[ix] = cp[iw][ix];
            }
            let cp2 = data[1];
            let cp3 = data[2];
            let cpnm1 = data[nx - 2];
            let cpnm2 = data[nx - 3];
            let mut a1 = crmul(cmul(cp2, conjg(cp3)), 2.0);
            let b1 = cadd(cmul(cp2, conjg(cp2)), cmul(cp3, conjg(cp3)));
            if b1.r == 0.0 && b1.i == 0.0 {
                a1 = cwp_cexp(cmplx(0.0, ((-w * dx) as f64 * 0.5 / v1 as f64) as f32));
            } else {
                a1 = cdiv(a1, b1);
            }
            if a1.i > 0.0 {
                a1 = cwp_cexp(cmplx(0.0, ((-w * dx) as f64 * 0.5 / v1 as f64) as f32));
            }
            let mut a2 = crmul(cmul(cpnm1, conjg(cpnm2)), 2.0);
            let b2 = cadd(cmul(cpnm1, conjg(cpnm1)), cmul(cpnm2, conjg(cpnm2)));
            if b2.r == 0.0 && b2.i == 0.0 {
                a2 = cwp_cexp(cmplx(0.0, ((-w * dx) as f64 * 0.5 / vn as f64) as f32));
            } else {
                a2 = cdiv(a2, b2);
            }
            if a2.i > 0.0 {
                a2 = cwp_cexp(cmplx(0.0, ((-w * dx) as f64 * 0.5 / vn as f64) as f32));
            }
            for ix in 0..nx {
                a[ix] = cmplx(s1[ix], s2[ix]);
                b[ix] = cmplx(
                    (1.0 - 2.0 * s1[ix] as f64) as f32,
                    (-2.0 * s2[ix] as f64) as f32,
                );
            }
            for ix in 1..nx - 1 {
                d[ix] = cadd(
                    cadd(cmul(data[ix + 1], a[ix + 1]), cmul(data[ix - 1], a[ix - 1])),
                    cmul(data[ix], b[ix]),
                );
            }
            d[0] = cadd(
                cmul(cadd(b[0], cmul(a[0], a1)), data[0]),
                cmul(data[1], a[1]),
            );
            d[nx - 1] = cadd(
                cmul(cadd(b[nx - 1], cmul(a[nx - 1], a2)), data[nx - 1]),
                cmul(data[nx - 2], a[nx - 2]),
            );
            for ix in 0..nx {
                data[ix] = cmplx(s1[ix], -s2[ix]);
                b[ix] = cmplx(
                    (1.0 - 2.0 * s1[ix] as f64) as f32,
                    (2.0 * s2[ix] as f64) as f32,
                );
            }
            let endl = cadd(b[0], cmul(data[0], a1));
            let endr = cadd(b[nx - 1], cmul(data[nx - 1], a2));
            for ix in 1..nx - 1 {
                a[ix] = data[ix + 1];
                c[ix] = data[ix - 1];
            }
            a[0] = data[1];
            c[nx - 1] = data[nx - 2];
            retris(&mut data, &a, &c, &b, endl, endr, nx, &d);
            for ix in 0..nx {
                cp[iw][ix] = data[ix];
            }
            w += dw;
        }
        if step == 0 {
            break;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn retris(
    data: &mut [Complex],
    a: &[Complex],
    c: &[Complex],
    b: &[Complex],
    endl: Complex,
    endr: Complex,
    nx: usize,
    d: &[Complex],
) {
    let mut e = vec![Complex::default(); nx];
    let mut f = vec![Complex::default(); nx];
    e[0] = cdiv(cneg(a[0]), endl);
    f[0] = cdiv(d[0], endl);
    for ix in 1..nx - 1 {
        let den = cadd(b[ix], cmul(c[ix], e[ix - 1]));
        e[ix] = cdiv(cneg(a[ix]), den);
        f[ix] = cdiv(csub(d[ix], cmul(f[ix - 1], c[ix])), den);
    }
    data[nx - 1] = cdiv(
        csub(d[nx - 1], cmul(f[nx - 2], c[nx - 2])),
        cadd(endr, cmul(c[nx - 2], e[nx - 2])),
    );
    for ix in (0..=nx - 2).rev() {
        data[ix] = cadd(cmul(data[ix + 1], e[ix]), f[ix]);
    }
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
            sx = tr.sx as f32 / (tr.scalco as i32).abs() as f32;
            gx = tr.gx as f32 / (tr.scalco as i32).abs() as f32;
            sy = tr.sy as f32 / (tr.scalco as i32).abs() as f32;
            gy = tr.gy as f32 / (tr.scalco as i32).abs() as f32;
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
