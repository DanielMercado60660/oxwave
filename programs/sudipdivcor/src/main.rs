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

//! SUDIPDIVCOR - Dip-dependent Divergence (spreading) correction.
//! Port of su/main/amplitudes/sudipdivcor.c (Revision 1.4).
//! Credits: CWP: Francesca Fazzari, May 1992.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::*;
use oxwave_cwp::cwp::{abs, min, nint, PI};
use oxwave_f_cubicspline::*;
use oxwave_linear_interp::*;
use oxwave_pfa_fft::*;

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get information from the first header */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;

    /* get parameters */
    let mut dx = 0.0f32;
    if getparfloat("dxcdp", &mut dx) == 0 {
        err!("dxcdp required");
    }
    let mut np = 0i32;
    let mut trans = 0i32;
    let mut norm = 0i32;
    let mut conv = 0i32;
    let mut verbose = 0i32;
    if getparint("np", &mut np) == 0 {
        np = 50;
    }
    if getparint("trans", &mut trans) == 0 {
        trans = 0;
    }
    if getparint("norm", &mut norm) == 0 {
        norm = 1;
    }
    if getparint("conv", &mut conv) == 0 {
        conv = 0;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* get velocity function */
    let mut vt = vec![0.0f32; nt];
    let mut tt = vec![0.0f32; nt];
    for it in 0..nt {
        tt[it] = it as f32 * dt;
    }
    let mut vfile = String::new();
    if getparstring("vfile", &mut vfile) == 0 {
        let mut ntau = countparval("tmig");
        if ntau == 0 {
            ntau = 1;
        }
        let mut tmig = vec![0.0f32; ntau as usize];
        if getparfloat("tmig", &mut tmig) == 0 {
            tmig[0] = 0.0;
        }
        let mut nv = countparval("vmig");
        if nv == 0 {
            nv = 1;
        }
        if nv != ntau {
            err!("number of tmig and vmig must be equal");
        }
        let mut vmig = vec![0.0f32; nv as usize];
        if getparfloat("vmig", &mut vmig) == 0 {
            vmig[0] = 1500.0;
        }
        for itau in 1..ntau as usize {
            if tmig[itau] <= tmig[itau - 1] {
                err!("tmig must increase monotonically");
            }
        }
        let mut t = 0.0f32;
        for it in 0..nt {
            intlin(
                &tmig,
                &vmig,
                vmig[0],
                vmig[ntau as usize - 1],
                &[t],
                &mut vt[it..it + 1],
            );
            t += dt;
        }
        if ntau as usize != nt {
            let mut vsind = vec![[0.0f32; 4]; ntau as usize];
            cmonot(&tmig, &vmig, &mut vsind);
            intcub(0, &tmig, &vsind, &tt, &mut vt);
        }
    } else {
        // The C source calls fread on an unchecked fopen result.
        let fp = fopen(&vfile, "r").unwrap_or_else(|| {
            unsafe extern "C" {
                fn raise(sig: i32) -> i32;
            }
            unsafe {
                raise(11);
            }
            exit(139);
        });
        if fread(&mut vt, &fp) != nt {
            err!("Not %d velocities in file %s", nt, vfile);
        }
    }
    checkpars();

    /* copy traces and headers to temporary files */
    let tfp = etmpfile();
    let hfp = etmpfile();
    let mut nx = 0usize;
    loop {
        nx += 1;
        fwrite(std::slice::from_ref(&tr.hdr), &hfp);
        fwrite(&tr.data[..nt], &tfp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    fseeko(&hfp, 0, SEEK_SET);
    fseeko(&tfp, 0, SEEK_SET);
    if verbose != 0 {
        fprintf!(stderr(), "\t%d traces input\n", nx);
    }

    /* determine wavenumber and frequency sampling */
    let nxfft = npfar(nx as i32);
    let ntfft = npfa(nt as i32);
    let nk = nxfft / 2 + 1;
    dx = (dx as f64 * 0.001) as f32;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;

    // Separate buffers permit safe slices for the same C transforms.
    let mut ptk = vec![Complex::default(); nt * nk as usize];
    let mut ptx = vec![0.0f32; nt * nxfft as usize];
    let mut divcor = vec![vec![0.0f32; nt.max(2)]; np as usize];
    // C reads divcor[0][1] even when nt=1; reserve the allocator padding.
    divcortable(nt, np, dt, &tt, &mut vt, &mut divcor, trans, norm);

    /* apply conventional correction if required */
    if conv == 1 {
        for ix in 0..nx {
            efread(&mut ptx[ix * nt..(ix + 1) * nt], &tfp);
            for it in 0..nt {
                ptx[ix * nt + it] *= divcor[0][it];
            }
        }
    } else {
        let kscl = (1.0 / nxfft as f64) as f32;
        for ix in 0..nx {
            efread(&mut ptx[ix * nt..(ix + 1) * nt], &tfp);
            for it in 0..nt {
                ptx[ix * nt + it] *= kscl;
            }
        }
        for ix in nx..nxfft as usize {
            for it in 0..nt {
                ptx[ix * nt + it] = 0.0;
            }
        }
        pfa2rc(-1, 2, nt as i32, nxfft, &ptx, &mut ptk);
        if verbose != 0 {
            fprintf!(stderr(), "\tFourier transform done\n");
        }
        let nkmax = min(nk, nint(PI / dt as f64 / vt[0] as f64 / dk as f64));
        let dpx = (1.0 / (np - 1) as f64 / vt[0] as f64) as f32;
        fprintf!(
            stderr(),
            "nkmax %d nk %d dk %f dpx %f \n",
            nkmax,
            nk,
            dk,
            dpx
        );
        for it in 0..nt {
            ptk[it].r *= divcor[0][it];
            ptk[it].i *= divcor[0][it];
        }
        let mut k = dk;
        for ik in 1..nkmax {
            if verbose != 0 && ik % (if nkmax / 10 > 0 { nkmax / 10 } else { 1 }) == 0 {
                fprintf!(stderr(), "\t%d of %d wavenumbers done\n", ik, nkmax);
            }
            dipfilt(
                k,
                dpx,
                dt,
                np,
                ntfft as usize,
                nt,
                &divcor,
                &mut ptk[ik as usize * nt..(ik as usize + 1) * nt],
            );
            k += dk;
        }
        pfa2cr(1, 2, nt as i32, nxfft, &ptk, &mut ptx);
        if verbose != 0 {
            fprintf!(stderr(), "\tinverse Fourier transform done\n");
        }
    }
    for ix in 0..nx {
        efread(std::slice::from_mut(&mut tr.hdr), &hfp);
        tr.data[..nt].copy_from_slice(&ptx[ix * nt..(ix + 1) * nt]);
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

/// Jacubowicz filter to apply dip-dependent divergence correction.
#[allow(clippy::too_many_arguments)]
fn dipfilt(
    k: f32,
    dpx: f32,
    dt: f32,
    np: i32,
    nw: usize,
    nt: usize,
    div: &[Vec<f32>],
    p: &mut [Complex],
) {
    let mut kq = vec![Complex::default(); nw];
    let mut qq = vec![Complex::default(); nw];
    let dw = (2.0 * PI / (nw as f32 * dt) as f64) as f32;
    let wny = (PI / dt as f64) as f32;
    let pscl = 0.5f32;
    let fftscl = (1.0 / nw as f64) as f32;
    let pmin = k / wny;
    for iw in 0..nw {
        qq[iw].r = 0.0;
        qq[iw].i = 0.0;
    }
    for ip in (0..np).rev() {
        let ph = dpx * (ip + 1) as f32;
        let pm = dpx * ip as f32;
        let pl = dpx * (ip - 1) as f32;
        let mut iwl = (k / ph / dw + 1.0f32) as i32;
        let mut iwh = (k / pl / dw) as i32;
        if (pl as f64) < 1.01 * pmin as f64 {
            iwh = if nw % 2 != 0 {
                (nw as i32 - 1) / 2
            } else {
                nw as i32 / 2 - 1
            };
        }
        if (pm as f64) < 1.01 * pmin as f64 {
            iwh = if nw % 2 != 0 {
                nw as i32 / 2
            } else {
                nw as i32 / 2 - 1
            };
        }
        if iwh >= iwl {
            for it in 0..nt {
                kq[it].r = p[it].r * div[ip as usize][it];
                kq[it].i = p[it].i * div[ip as usize][it];
            }
            for it in nt..nw {
                kq[it].r = 0.0;
                kq[it].i = 0.0;
            }
            pfacc(1, &mut kq);
            for iw in iwl..=iwh {
                qq[iw as usize].r += kq[iw as usize].r * pscl;
                qq[iw as usize].i += kq[iw as usize].i * pscl;
            }
            iwl = nw as i32 - iwl;
            iwh = nw as i32 - iwh;
            for iw in iwh..=iwl {
                qq[iw as usize].r += kq[iw as usize].r * pscl;
                qq[iw as usize].i += kq[iw as usize].i * pscl;
            }
        }
        if (pm as f64) < 1.01 * pmin as f64 {
            break;
        }
    }
    pfacc(-1, &mut qq);
    for it in 0..nt {
        p[it].r = qq[it].r * fftscl;
        p[it].i = qq[it].i * fftscl;
    }
}

/// Divergence correction for arbitrary v(z) as a function A(px,t).
#[allow(clippy::too_many_arguments)]
fn divcortable(
    nt: usize,
    np: i32,
    dt: f32,
    tt: &[f32],
    vt: &mut [f32],
    divcor: &mut [Vec<f32>],
    trans: i32,
    norm: i32,
) {
    let mut vind = vec![[0.0f32; 4]; nt];
    /* use v/2 to compensate for v(two-way t) */
    for it in 0..nt {
        vt[it] = (vt[it] as f64 * 0.0005) as f32;
    }
    cmonot(tt, vt, &mut vind);
    for ip in 0..np as usize {
        for it in 0..nt {
            divcor[ip][it] = 0.0;
        }
    }
    let v0 = vind[0][0];
    let dpx = (1.0 / v0 as f64 / (np - 1) as f64) as f32;
    for ip in 0..np {
        let px = ip as f32 * dpx;
        let mut vel = v0;
        let mut vel1 = vind[0][1];
        let mut vel2 = vind[0][2];
        let mut tau = 0.0f32;
        let px2 = px * px;
        let mut gamma =
            (px2 as f64 * (((vel1 / vel) as f64).powf(2.0) - (vel2 / vel) as f64)) as f32;
        let mut p = (1.0 / v0 as f64) as f32;
        let mut q = 0.0f32;
        let mut sigma = 0.0f32;
        let (mut pz, mut flag) = if ip == np - 1 {
            (0.0f32, 0)
        } else {
            (((1.0f32 / vel / vel - px2) as f64).sqrt() as f32, 1)
        };
        for it in 1..nt {
            tau += vel * pz * dt;
            if tau >= 0.0 {
                let velold = vel;
                let vel1old = vel1;
                intcub(0, tt, &vind, &[tau], std::slice::from_mut(&mut vel));
                intcub(1, tt, &vind, &[tau], std::slice::from_mut(&mut vel1));
                intcub(2, tt, &vind, &[tau], std::slice::from_mut(&mut vel2));
                if flag == 1 {
                    let pz1 = pz - vel1old / velold / velold * dt;
                    pz = ((1.0f32 / vel / vel - px2) as f64).sqrt() as f32;
                    if pz1 < 0.0 || px * vel > 1.0 {
                        flag = 0;
                        pz = pz1;
                    }
                } else {
                    pz += -vel1old / velold / velold * dt;
                }
                sigma += velold * velold * dt;
                let gammaold = gamma;
                gamma =
                    (px2 as f64 * (((vel1 / vel) as f64).powf(2.0) - (vel2 / vel) as f64)) as f32;
                let alpha = (1.0 - gamma as f64 * ((dt / 2.0f32 * vel) as f64).powf(2.0)) as f32;
                let qold = q;
                let pold = p;
                q = ((pold * (dt / 2.0f32 / alpha) * (vel * vel + velold * velold)) as f64
                    + (qold / alpha) as f64
                        * (1.0 + ((dt * vel / 2.0f32) as f64).powf(2.0) * gammaold as f64))
                    as f32;
                p = ((qold * (dt / 2.0f32 / alpha) * (gamma + gammaold)) as f64
                    + (pold / alpha) as f64
                        * (1.0 + ((dt * velold / 2.0f32) as f64).powf(2.0) * gamma as f64))
                    as f32;
                divcor[ip as usize][it] =
                    ((2.0f32 * sigma * abs(q) / if trans == 1 { vel } else { v0 }) as f64).sqrt()
                        as f32;
            } else {
                vel = v0;
                sigma += vel * vel * dt;
                let gammaold = gamma;
                gamma = 0.0;
                let _alpha = 1.0f32;
                let qold = q;
                let pold = p;
                q = ((pold * (dt * vel * vel)) as f64
                    + qold as f64
                        * (1.0 + ((dt * vel / 2.0f32) as f64).powf(2.0) * gammaold as f64))
                    as f32;
                p = qold * dt / 2.0f32 * gammaold + pold;
                divcor[ip as usize][it] =
                    ((2.0f32 * sigma * abs(q) / if trans == 1 { vel } else { v0 }) as f64).sqrt()
                        as f32;
            }
        }
    }
    if norm == 1 {
        let denom = (1.0 / divcor[0][1] as f64) as f32;
        for ip in (0..np as usize).rev() {
            for it in 1..nt {
                divcor[ip][it] *= denom;
            }
        }
    }
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUDIPDIVCOR - Dip-dependent Divergence (spreading) correction\t",
    " \t\t\t\t\t\t\t\t",
    "\tsudipdivcor <stdin >stdout  [optional parms]\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t",
    "\tdxcdp\tdistance between sucessive cdps\tin meters\t",
    " \t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    "\tnp=50\t\tnumber of slopes\t\t\t",
    "\ttmig=0.0\ttimes corresponding to rms velocities in vmig",
    "\tvmig=1500.0\trms velocities corresponding to times in tmig",
    "\tvfile=binary\t(non-ascii) file containing velocities vmig(t) ",
    "\tconv=0\t\t=1 to apply the conventional divergence correction",
    "\ttrans=0\t\t=1 to include transmission factors \t",
    "\tverbose=0\t=1 for diagnostic print\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The tmig, vmig arrays specify an rms velocity function of time.",
    " Linear interpolation and constant extrapolation is used to determine",
    " rms velocities at times not specified.  Values specified in tmig",
    " must increase monotonically.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Alternatively, rms velocities may be stored in a binary file",
    " containing one velocity for every time sample.  If vfile is\t",
    " specified, then the tmig and vmig arrays are ignored.\t\t",
    " The time of the first sample is assumed to be constant, and is",
    " taken as the value of the first trace header field delrt. \t",
    " \t\t\t\t\t\t\t\t",
    " Whereas the conventional divergence correction (sudivcor) is\t",
    " valid only for horizontal reflectors, which have zero reflection",
    " slope, the dip-dependent divergence correction is valid for any",
    " reflector dip or reflection slope.  Only the conventional\t",
    " correction will be applied to the data if conv=1 is specified. ",
    " Note that the conventional correction over-amplifies\t\t",
    " reflections from dipping beds\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " The transmission factor should be applied when the divergence ",
    " corrected data is to be migrated with a reverse time migration ",
    " based on the constant density wave equation.\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Trace header fields accessed:  ns, dt, delrt\t\t\t",
    " \t\t\t\t\t\t\t\t",
];
